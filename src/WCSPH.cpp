/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/***			        + XSPH Correction (Also described in Monaghan)                    ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)   (Currently inactive)           ***/
/*** Smoothing Kernel: Wendland's C2 ***/
/*** Integrator: Newmark-Beta ****/
/*** Variable Timestep Criteria: CFL + Monaghan, J.J. (1989) conditions ***/

#include <chrono>

#include "Var.h"
#include "IO.h"
#include "Neighbours.h"
#include "Kernel.h"
#include "Init.h"
#include "Add.h"
#include "Resid.h"
#include "Crossing.h"
#include "Newmark_Beta.h"

using namespace std::chrono;
using namespace nanoflann;


void First_Step(SIM& svar, const FLUID& fvar, const AERO& avar, const outl& outlist, State& pnp1, State& airP)
{
	const uint start = svar.bndPts;
	const uint end = svar.totPts;

	#if DEBUG 
		dbout << "Starting first step. ";
		dbout << "  Start index: " << start << "  End index: " << end << endl;
	#endif

	std::vector<std::vector<Part>> neighb;
	neighb.reserve(end);
	for(uint ii = 0; ii < start; ++ii)
		neighb.emplace_back();

	/*Check if a particle is running low on neighbours, and add ficticious particles*/
	vector<vector<Part>> air;
	#pragma omp parallel shared(svar, pnp1, outlist)
	{
		std::vector<std::vector<Part>> localN;
		vector<vector<Part>> localA;
		#pragma omp for schedule(static) nowait 
		for (uint ii = start; ii < end; ++ii)
		{
			std::vector<Part> temp;
			if(svar.ghost == 1 && pnp1[ii].b == 2 && outlist[ii].size() < avar.nfull &&
				outlist[ii].size() > 0.4*avar.nfull)
				temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,ii,pnp1,outlist);

			localA.emplace_back(temp);

			for(auto j:outlist[ii])
				temp.emplace_back(Part(pnp1[j]));

			localN.emplace_back(temp);
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),localN.begin(),localN.end());
    		air.insert(air.end(),localA.begin(),localA.end());
    	}
	}
	airP.clear();

	for(uint ii = 0; ii < air.size(); ++ii)
		for(uint jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	#pragma omp parallel for shared(outlist)
	for(uint ii = start; ii < end; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}

	/*Previous State for error calc*/
	vector<StateVecD> xih(svar.totPts);
	#pragma omp parallel for shared(pnp1)
	for (uint  ii=0; ii < end; ++ii)
		xih[ii] = pnp1[ii].xi;
	
	vector<StateVecD> res(svar.totPts,StateVecD::Zero());
	vector<StateVecD> Af(end,StateVecD::Zero());
	vector<ldouble> Rrho(svar.totPts,0.0);
	Forces(svar,fvar,avar,pnp1,neighb,outlist,res,Rrho,Af); /*Guess force at time n+1*/

	/*Find maximum safe timestep*/
	vector<StateVecD>::iterator maxfi = std::max_element(res.begin(),res.end(),
		[](StateVecD p1, StateVecD p2){return p1.norm() < p2.norm();});
	ldouble maxf = maxfi->norm();
	ldouble dtf = sqrt(fvar.H/maxf);
	ldouble dtcv = fvar.H/(fvar.Cs+svar.maxmu);
	const ldouble dt = 0.3*std::min(dtf,dtcv);

	#pragma omp parallel for shared(res, Rrho)
	for(uint ii = 0; ii < end; ++ii)
	{
		pnp1[ii].f = res[ii];
		pnp1[ii].Rrho = Rrho[ii];
		pnp1[ii].Af = Af[ii]; 
		xih[ii] = pnp1[ii].xi + dt*pnp1[ii].v;
	}
#if DEBUG 
	ldouble errsum = 0.0;

	for (uint ii = start; ii < end; ++ii)
	{
		StateVecD r = xih[ii]-pnp1[ii].xi;
		errsum += r.squaredNorm();
	}

	ldouble error=log10(sqrt(errsum/(double(svar.totPts))));

	
		dbout << "Exiting first step. Error: " << error << endl;
#endif
}

void Find_MinMax(SIM& svar, const State& pnp1)
{
	/*Find the max and min positions*/
		auto xC = std::minmax_element(pnp1.begin(),pnp1.end(),
					[](Particle p1, Particle p2){return p1.xi(0)< p2.xi(0);});
		auto yC = std::minmax_element(pnp1.begin(),pnp1.end(),
					[](Particle p1, Particle p2){return p1.xi(1)< p2.xi(1);});
		#if SIMDIM == 3
			auto zC = std::minmax_element(pnp1.begin(),pnp1.end(),
					[](Particle p1, Particle p2){return p1.xi(2)< p2.xi(2);});

			StateVecD minC = StateVecD(xC.first->xi(0),yC.first->xi(1),zC.first->xi(2));
			StateVecD maxC = StateVecD(xC.second->xi(0),yC.second->xi(1),zC.second->xi(2));
		#else 
			StateVecD minC = StateVecD(xC.first->xi(0),yC.first->xi(1));
			StateVecD maxC = StateVecD(xC.second->xi(0),yC.second->xi(1));
		#endif

		/*Check if they are more than the previous step*/
		for (uint dim = 0; dim < SIMDIM; ++dim)
		{
			if(minC(dim) < svar.minC(dim))
				svar.minC(dim) = minC(dim);

			if(maxC(dim) > svar.maxC(dim))
				svar.maxC(dim) = maxC(dim);
		}	
}


int main(int argc, char *argv[])
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
	Eigen::initParallel();
	// omp_set_num_threads(NTHREADS);

    double duration;
    double error = 0;
    cout.width(13);
    cout << std::scientific << std::left;
	
	write_header();
    

    // cout << MEPSILON << endl;
    // for(double ii = 1; ii < 200; ii++)
    // {
    // 	if(double(1.0+pow(2.0,-ii)) == double(1.0))
    // 	{
    // 		cout << ii << "  " << pow(2.0,-ii) << endl;
    // 		break;
    // 	}
    // }

    /******* Define the global simulation parameters ******/
	SIM svar;
	FLUID fvar;
	AERO avar;
	outl outlist;
	MESH cells;

	GetInput(argc,argv,svar,fvar,avar);
	
	if(MakeOutputDir(argc,argv,svar))
	{
		cout << "Couldn't make output directory. Please check permissions." << endl;
		exit(-1);
	}

	if(svar.Bcase == 6)
	{
		#if SIMDIM == 3
		Read_TAUMESH_FACE(svar,cells,fvar);
		// Read_TAUMESH(svar,cells,fvar);
		#else
		Read_TAUMESH(svar,cells,fvar);
		#endif
	}	

	// Check if cells have been initialsed before making a tree off it
	if(cells.cCentre.size() == 0)
		cells.cCentre.emplace_back(StateVecD::Zero());

	Vec_Tree CELL_INDEX(SIMDIM,cells.cCentre,20);
	CELL_INDEX.index->buildIndex();
	if(svar.Bcase == 6)
	{
		cout << "Building cell neighbours..." << endl;
		FindCellNeighbours(CELL_INDEX, cells.cCentre, cells.cNeighb);
	}


	/*Make a guess of how many there will be...*/
	int partCount = ParticleCount(svar);
    ///****** Initialise the particles memory *********/
	State pn;	    /*Particles at n   */
	State pnp1; 	/*Particles at n+1 */
	State airP;

	cout << "Final particle count:  " << partCount << endl;
	svar.finPts = partCount;
	pn.reserve(partCount);
  	pnp1.reserve(partCount);
	
	
	// if (svar.Bcase == 6){
	// 	Write_Mesh_Data(svar,cells);
	// }	

	InitSPH(svar,fvar,avar,pn,pnp1);

	///********* Tree algorithm stuff ************/
	Sim_Tree NP1_INDEX(SIMDIM,pnp1,20);
	NP1_INDEX.index->buildIndex();
	FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);


	///*** Perform an iteration to populate the vectors *****/
	// std::vector<std::vector<uint>>::iterator nfull = 
	// 	std::max_element(outlist.begin(),outlist.end(),
	// 	[](std::vector<uint> p1, std::vector<uint> p2){return p1.size()< p2.size();});

	

	#if SIMDIM == 3
		// avar.nfull = (2.0/3.0) * double(nfull->size());
		avar.nfull = 1.713333e+02;
		svar.nfull = 257;
	#endif
	#if SIMDIM == 2
		// avar.nfull = (2.0/3.0) * double(nfull->size());
		// avar.nfull = 32.67;
		avar.nfull = 37;
		svar.nfull = 48;
	#endif

	// cout << avar.nfull << endl;
	First_Step(svar,fvar,avar,outlist,pnp1,airP);
	///*************** Open simulation files ***************/
	std::ofstream f1,f2,f3,fb,fg;
	// NcFile* fn;
	// h5_file_t fh5;

	string framef = svar.outfolder;
	framef.append("frame.info");
	f2.open(framef, std::ios::out);

	f1 << std::scientific << std::setprecision(6);
	f3 << std::scientific << std::setw(10);
	uint ghost_strand = 0;
	if(svar.boutform == 0)
		ghost_strand = 2;
	else
		ghost_strand = 3;


	if(svar.outtype == 0 )
	{
		if (svar.outform > 4)
		{	
			cout << "Output type not within design. Outputting fluid data..." << endl;
			svar.outform = 1;
		}

		
		/*Write sim particles*/
		
		Init_Binary_PLT(svar,"Fuel.szplt","Simulation Particles");
		Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1);


		if (svar.Bcase != 0 && svar.Bcase !=5)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"Boundary.szplt","Boundary Particles");
			Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2); 
			if(svar.boutform == 0)
				TECEND142();			
		}	

		if (svar.ghost == 1 && svar.gout == 1)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"Ghost.szplt","Ghost Particles");		
		}	
	}
	else if (svar.outtype == 1)
	{

		if (svar.Bcase != 0 && svar.Bcase != 5)
		{	/*If the boundary exists, write it.*/
			string bfile = svar.outfolder;
			bfile.append("Boundary.plt");
			fb.open(bfile, std::ios::out);
			if(fb.is_open())
			{
				if(svar.boutform == 0)
				{
					Write_Boundary_ASCII(fb,svar,pnp1);
					fb.close();
				}
				else
				{
					State empty;
					Write_ASCII_header(fb,svar);
					Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
				}
			}
			else
			{
				cerr << "Error opening boundary file." << endl;
				exit(-1);
			}
		}

		/* Write first timestep */
		string mainfile = svar.outfolder;
		mainfile.append("Fuel.plt");
		f1.open(mainfile, std::ios::out);
		if(f1.is_open())
		{
			Write_ASCII_header(f1,svar);
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,"Fuel");
			// if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
			// 	Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,airP);

		}
		else
		{
			cerr << "Failed to open fuel.plt. Stopping." << endl;
			exit(-1);
		}

		if(svar.ghost == 1 && svar.gout == 1)
		{
			string ghostfile = svar.outfolder;
			ghostfile.append("Ghost.plt");
			fg.open(ghostfile,std::ios::out);
			if(fg.is_open())
			{
				Write_ASCII_header(fg,svar);
				Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
			}
		}
	}
	// else if (svar.outtype == 2)
	// {
	// 	string mainfile = svar.outfolder;
	// 	mainfile.append("/Fuel.h5part");
	// 	fn = new NcFile(mainfile, NcFile::replace);
	// 	Write_CDF_File(*fn,svar,pnp1);
	// 	if (svar.Bcase != 0 && svar.Bcase !=5)
	// 		Write_Boundary_CDF(svar, pnp1);
	// }
	// else if (svar.outtype == 3)
	// {
	// 	fh5 = H5OpenFile("testfile.h5part", H5_O_WRONLY, H5_PROP_DEFAULT);
	// 	H5SetStepNameFormat(fh5,"Step",6);
	// 	Write_H5_File(fh5,svar,pnp1);
	// }
	else
	{
		cerr << "Output type ambiguous. Please select 0, 1 or 2 for output data type." << endl;
		exit(-1);
	}

	#if SIMDIM == 3
		if(svar.Bcase == 4)
			svar.vortex.write_VLM_Panels(svar.outfolder);		
	#endif
	/*Timing calculation + error sum output*/
	t2 = high_resolution_clock::now();
	duration = duration_cast<microseconds>(t2-t1).count()/1e6;
	cout << "Frame: " << 0 << "  Sim Time: " << svar.t << "  Compute Time: "
	<< duration <<"  Error: " << error << endl;
	
	f2 << "Frame: " << 0 << endl;
	f2 << "Total Points: " << svar.totPts << " Boundary Points: " << svar.bndPts 
		<< " Fluid Points: " << svar.simPts << endl;
	f2 << "Sim Time:  " << std::scientific << std::setprecision(6) << svar.t 
	<< " Comp Time: " << std::fixed << duration << " Error: " << 0 << " Sub-iterations: " 
    << 0 << endl;
	
	Find_MinMax(svar,pnp1);	
	
	f2 << "Minimum Coords:" << endl << std::scientific << std::setprecision(8) 
		<< svar.minC(0) << " " << svar.minC(1) << " ";
	#if SIMDIM ==3
	f2 << svar.minC(2) << " ";
	#endif
	f2 << endl << "Maximum Coords:" << endl  << svar.maxC(0) << " " << svar.maxC(1);
	#if SIMDIM ==3
	f2 <<  " " << svar.maxC(2); 
	#endif
	f2 << endl;

	///************************* MAIN LOOP ********************/
	svar.frame = 0;
#ifdef DEBUG
	const ldouble a = 1 - svar.gamma;
	const ldouble b = svar.gamma;
	const ldouble c = 0.5*(1-2*svar.beta);
	const ldouble d = svar.beta;
	const ldouble B = fvar.B;
	const ldouble gam = fvar.gam;
	dbout << "Newmark Beta integration parameters" << endl;
	dbout << "a: " << a << "  b: " << b << endl;
	dbout << "c: " << c << "  d: " << d << endl;
	dbout << "B: " << B << "  gam: " << gam << endl << endl; 
#endif

	for (uint frame = 1; frame<= svar.Nframe; ++frame)
	{
		int stepits=0;
		double stept=0.0;
		while (stept<svar.framet)
		{
		    error = Newmark_Beta(NP1_INDEX,CELL_INDEX,svar,fvar,avar,cells,pn,pnp1,airP,outlist);
		    stept+=svar.dt;
		    ++stepits;
		    //cout << svar.t << "  " << svar.dt << endl;
		}
		++svar.frame;

		Find_MinMax(svar,pnp1);		

		t2= high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;

		/*Write each frame info to file (Useful to debug for example)*/
		f2 << endl;
		f2 << "Frame: " << frame << endl;

		f2 << "Total Points: " << svar.totPts << " Boundary Points: " << svar.bndPts 
			<< " Fluid Points: " << svar.simPts << endl;
		f2 << "Sim Time:  " << std::scientific << std::setprecision(6) << svar.t 
		<< " Comp Time: " << std::fixed << duration << " Error: " << error << " Sub-iterations: " 
	    << stepits << endl;
		
		
		f2 << "Minimum Coords:" << endl << std::scientific << std::setprecision(8) 
			<< svar.minC(0) << " " << svar.minC(1) << " ";
		#if SIMDIM ==3
		f2 << svar.minC(2) << " ";
		#endif
		f2 << endl << "Maximum Coords:" << endl  << svar.maxC(0) << " " << svar.maxC(1);
		#if SIMDIM ==3
		f2 <<  " " << svar.maxC(2);
		#endif
		f2 << endl;
		

		if(svar.outframe !=0)
		{
			if (frame % svar.outframe == 0 )
			{	/*Output to console every 20 or so steps*/
			  	cout << "Frame: " << frame << "  Sim Time: " << svar.t-svar.dt << "  Compute Time: "
			  	<< duration <<"  Error: " << error << endl;
			}
		}


		// DensityReinit(fvar, pnp1, outlist);
		if (svar.outtype == 0)
		{
			if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
			{	/*Write boundary particles*/
				Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2); 
			}
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1); /*Write sim particles*/
			if(svar.ghost == 1 && svar.gout == 1 && airP.size() != 0)
				Write_Binary_Timestep(svar,airP,0,airP.size(),"Ghost",ghost_strand);
		} 
		else if (svar.outtype == 1)
		{
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,"Fuel");
			if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
			{
				State empty;
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
			}

			if(svar.ghost == 1 && svar.gout == 1)
				Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
		}
		// else if (svar.outtype == 2)
		// {
		// 	Write_CDF_File(*fn,svar,pnp1);
		// }
		// else if (svar.outtype == 3)
		// {
		// 	Write_H5_File(fh5,svar,pnp1);
		// }
	}

	/*Wrap up simulation files and close them*/
	if (f1.is_open())
		f1.close();
	if (f2.is_open())
		f2.close();
	if (f3.is_open())
		f3.close();
	if (fb.is_open())
		fb.close();
	if (fg.is_open())
		fg.close();
	
	if(svar.outtype == 0)
	{
		if(TECEND142())
			exit(-1);
	}
	// else if (svar.outtype == 3)
	// {
	// 	if(H5CloseFile(fh5))
	// 		exit(-1);
	// }
	
	cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.t << " seconds" << endl;

    /*Perform post processing if desired.*/
    if(svar.afterSim == 1)
    {
    	fvar.H = svar.postRadius;
	  	fvar.HSQ = fvar.H*fvar.H; 
		fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/

    	TECMESH grid;
    	grid.DoPostProcessing(svar,fvar);
    }
	return 0;
}
