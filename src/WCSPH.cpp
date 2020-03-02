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


void First_Step(SIM& svar, const FLUID& fvar, const AERO& avar, const MESH& cells, const outl& outlist, State& pnp1, State& airP)
{
	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;

	#if DEBUG 
		dbout << "Starting first step. ";
		dbout << "  Start index: " << start << "  End index: " << end << endl;
	#endif

	std::vector<std::vector<Part>> neighb;
	neighb.reserve(end);
	for(size_t ii = 0; ii < start; ++ii)
		neighb.emplace_back();


	/*Check if a particle is running low on neighbours, and add ficticious particles*/
	vector<vector<Part>> air;
	#pragma omp parallel shared(svar, pnp1, outlist)
	{
		std::vector<std::vector<Part>> localN;
		vector<vector<Part>> localA;
		#pragma omp for schedule(static) nowait 
		for (size_t ii = start; ii < end; ++ii)
		{
			std::vector<Part> temp;
			if(svar.ghost == 1 && pnp1[ii].b == FREE && outlist[ii].size() < avar.nfull &&
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

	for(size_t ii = 0; ii < air.size(); ++ii)
		for(size_t jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	#pragma omp parallel for shared(outlist)
	for(size_t ii = start; ii < end; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}

	/*Previous State for error calc*/
	vector<StateVecD> xih(svar.totPts);
	#pragma omp parallel for shared(pnp1)
	for (size_t  ii=0; ii < end; ++ii)
		xih[ii] = pnp1[ii].xi;
	
	// vector<StateVecD> vPert(end,StateVecD::Zero());
	vector<StateVecD> res(end,StateVecD::Zero());
	vector<StateVecD> Af(end,StateVecD::Zero());
	vector<real> Rrho(svar.totPts,0.0);
	Forces(svar,fvar,avar,cells,pnp1,neighb,outlist/*,vPert*/,res,Rrho,Af); /*Guess force at time n+1*/

	/*Find maximum safe timestep*/
	vector<StateVecD>::iterator maxfi = std::max_element(res.begin(),res.end(),
		[](StateVecD p1, StateVecD p2){return p1.norm() < p2.norm();});
	real maxf = maxfi->norm();
	real dtf = sqrt(fvar.H/maxf);
	real dtcv = fvar.H/(fvar.Cs+svar.maxmu);
	const real dt = 0.3*std::min(dtf,dtcv);

	#pragma omp parallel for shared(res, Rrho)
	for(size_t ii = 0; ii < end; ++ii)
	{
		pnp1[ii].f = res[ii];
		pnp1[ii].Rrho = Rrho[ii];
		pnp1[ii].Af = Af[ii]; 
		xih[ii] = pnp1[ii].xi + dt*pnp1[ii].v;
	}
#if DEBUG 
	real errsum = 0.0;

	for (size_t ii = start; ii < end; ++ii)
	{
		StateVecD r = xih[ii]-pnp1[ii].xi;
		errsum += r.squaredNorm();
	}

	real error=log10(sqrt(errsum/(real(svar.totPts))));

	
		dbout << "Exiting first step. Error: " << error << endl;
#endif
}

std::pair<StateVecD,StateVecD> Find_MinMax(SIM& svar, const State& pnp1)
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

	return std::pair<StateVecD,StateVecD>(minC,maxC);
}

uint Check_Pipe(Vec_Tree& CELL_INDEX, const MESH& cells, const StateVecD& testp)
{
    StateVecD rayp;
#if SIMDIM == 3    
    rayp = testp;
    rayp(0) += 1e+10;	    
#endif
	const size_t num_results = 20;
    vector<size_t> ret_indexes(num_results);
    vector<real> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<real> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    
    CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));
	
	uint inside_flag=0;
    for(auto index:ret_indexes)
    {   
        inside_flag = 0;
        uint line_flag = 0;

        for (uint const& findex:cells.cFaces[index] ) 
        {
            const vector<size_t>& face = cells.faces[findex];

#if SIMDIM == 3            
            if(Crossings3D(cells.verts,face,testp,rayp))
#else
            if(Crossings2D(cells.verts,face,testp))
#endif                
            {   
                inside_flag=!inside_flag;
                if ( line_flag ) break; //Convex assumption

                //  note that one edge has been hit by the ray's line 
                line_flag = TRUE;
            }
        }

        if(inside_flag == TRUE)
        {
        	// cout << "Found cell." << endl;
        	break;
        }
	}

	if(inside_flag == 0)
	{
		Write_Containment(ret_indexes, cells, testp);
	}

	return inside_flag;
}


int main(int argc, char *argv[])
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
	Eigen::initParallel();
	// omp_set_num_threads(NTHREADS);

    real duration;
    real error = 0;
    cout.width(13);
    cout << std::scientific << std::left << std::setprecision(2);
	
	write_header();
    
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
		Read_TAUMESH_EDGE(svar,cells,fvar);
		#endif
	}	

	cout << std::setprecision(5);
	cout << "Adjusted Start Coordinates: " << endl;
	cout << svar.Start(0) << "  " << svar.Start(1);
#if SIMDIM == 3
	cout << "  " << svar.Start(2); 
#endif
	cout << endl << endl;
	/*Make a guess of how many there will be...*/
	// int partCount = ParticleCount(svar);
    ///****** Initialise the particles memory *********/
	State pn;	    /*Particles at n   */
	State pnp1; 	/*Particles at n+1 */
	State airP;

	cout << "Final particle count:  " << svar.nmax << endl;
	svar.finPts = svar.nmax;
	pn.reserve(svar.nmax);
  	pnp1.reserve(svar.nmax);

	
	
	// if (svar.Bcase == 6){
	// 	Write_Mesh_Data(svar,cells);
	// }	

  	if(svar.restart == 1)
  	{
  		Restart(svar,pn,cells);
  		pnp1 = pn;
  	}
  	else
  		InitSPH(svar,fvar,avar,pn,pnp1);

	// Check if cells have been initialsed before making a tree off it
	if(cells.cCentre.size() == 0)
		cells.cCentre.emplace_back(StateVecD::Zero());
	if (cells.bVerts.size() == 0)
		cells.bVerts.emplace_back(StateVecD::Zero());

	// Vec_Tree CELL_INDEX(SIMDIM,cells.cCentre,20);
		
	cout << "Starting counts: " << endl;
	cout << "Boundary: " << svar.bndPts << "  Sim: " << svar.simPts << endl;
	 

	///********* Tree algorithm stuff ************/
	KDTREE TREE(pnp1,cells);
	// Sim_Tree NP1_INDEX(SIMDIM,pnp1,20);
	TREE.CELL.index->buildIndex();
	TREE.NP1.index->buildIndex();

	// if(svar.Bcase == 6)
	// {
	// 	cout << "Building cell neighbours..." << endl;
	// 	FindCellNeighbours(TREE.CELL, cells.cCentre, cells.cNeighb);
	// }

	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	#if SIMDIM == 3
		// avar.nfull = (2.0/3.0) * real(nfull->size());
		avar.nfull = 1.713333e+02;
		svar.nfull = 257;
	#endif
	#if SIMDIM == 2
		// avar.nfull = (2.0/3.0) * real(nfull->size());
		// avar.nfull = 32.67;
		avar.nfull = 37;
		svar.nfull = 48;
	#endif

	///*** Perform an iteration to populate the vectors *****/
	First_Step(svar,fvar,avar,cells,outlist,pnp1,airP);

	///*************** Open simulation files ***************/
	std::fstream f1,f2,f3,fb,fg;
	// NcFile* fn;
	// h5_file_t fh5;

	string framef = svar.outfolder;
	framef.append("frame.info");
	if(svar.restart == 1)
		f2.open(framef, std::ios::out | std::ios::app);
	else
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
		/*Write sim particles*/
		
		Init_Binary_PLT(svar,"Fuel.szplt","Simulation Particles");
		// if(svar.restart == 0)
		// {
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1);
		// }

		if (svar.Bcase != 0 && svar.Bcase !=5 && svar.restart == 0)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"Boundary.szplt","Boundary Particles");

			if(svar.boutform == 0)
			{   //Don't write a strand, so zone is static. 
				Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",0); 
				TECEND142();			
			}
			else
				Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2); 
		}	

		if (svar.ghost == 1 && svar.gout == 1)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"Ghost.szplt","Ghost Particles");		
		}	
	}
	else if (svar.outtype == 1)
	{

		if (svar.Bcase != 0 && svar.Bcase != 5 && svar.restart == 0)
		{	/*If the boundary exists, write it.*/
			string bfile = svar.outfolder;
			bfile.append("Boundary.plt");
			fb.open(bfile, std::ios::out);
			if(fb.is_open())
			{
				Write_ASCII_header(fb,svar);
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
				if(svar.boutform == 0)
					fb.close();
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
		cerr << "Output type ambiguous. Please select 0 or 1 for output data type." << endl;
		exit(-1);
	}

	if(svar.Bcase == 6)
	{
				// Check if the pipe is inside the mesh
		real holeD = svar.Jet(0)+8*svar.dx; /*Diameter of hole (or width)*/
		real stepb = (svar.Pstep*svar.Bstep);
		real r = 0.5*holeD;
		#if SIMDIM == 3
    	real dtheta = atan((stepb)/(r));
		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		{
			StateVecD xi(r*sin(theta), 0.0, r*cos(theta));
			/*Apply Rotation...*/
			xi = svar.Rotate*xi;
			xi += svar.Start;
		    if(!Check_Pipe(TREE.CELL, cells, xi))
		    {

		    	cout << "Some of the pipe is outside of the simulation mesh." << endl;
		    	cout << "Fuel will be excessively close to begin." << endl;
		    	exit(-1);
		    }

    	}
    	#else
    	for(real x = -r; x <= r; x+=stepb)
    	{

    		StateVecD xi(x,0.0);
    		xi = svar.Rotate*xi;
			xi += svar.Start;

			// cout << "Checking point: " << xi(0) << "  " << xi(1) << endl;
		    if(!Check_Pipe(TREE.CELL, cells, xi))
		    {
		    	cout << "Some of the pipe is outside of the simulation mesh." << endl;
		    	cout << "Fuel will be excessively close to begin." << endl;
		    	exit(-1);
		    }
    	}
    	#endif
	}

	#if SIMDIM == 3
		if(svar.Bcase == 4)
			svar.vortex.write_VLM_Panels(svar.outfolder);		
	#endif
	/*Timing calculation + error sum output*/
	t2 = high_resolution_clock::now();
	duration = duration_cast<microseconds>(t2-t1).count()/1e6;
	
	if(svar.restart != 1)
		svar.frame = 0;
	else
	{
		cout << "Restarting simulation..." << endl;
	}

	cout << "Frame: " << svar.frame << "  Sim Time: " << svar.t << "  Compute Time: "
	<< duration <<"  Error: " << error << endl;

	if(svar.restart != 1)
	{
		f2 << "Frame: " << svar.frame << endl;
		f2 << "Total Points: " << svar.totPts << " Boundary Points: " << svar.bndPts 
			<< " Fluid Points: " << svar.simPts << endl;
		f2 << "Sim Time:  " << std::scientific << std::setprecision(6) << svar.t 
		<< " Comp Time: " << std::fixed << duration << " Error: " << 0 << " Sub-iterations: " 
	    << 0 << endl;
		f2 << "Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;
		

		// f2 << "Minimum Coords:" << endl << std::scientific << std::setprecision(8) 
		// 	<< svar.minC(0) << " " << svar.minC(1) << " ";
		// #if SIMDIM ==3
		// f2 << svar.minC(2) << " ";
		// #endif
		// f2 << endl << "Maximum Coords:" << endl  << svar.maxC(0) << " " << svar.maxC(1);
		// #if SIMDIM ==3
		// f2 <<  " " << svar.maxC(2); 
		// #endif
		// f2 << endl;
	}
	///************************* MAIN LOOP ********************/
	
#ifdef DEBUG
	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	const real B = fvar.B;
	const real gam = fvar.gam;
	dbout << "Newmark Beta integration parameters" << endl;
	dbout << "a: " << a << "  b: " << b << endl;
	dbout << "c: " << c << "  d: " << d << endl;
	dbout << "B: " << B << "  gam: " << gam << endl << endl; 
#endif

	for (uint frame = svar.frame+1; frame<= svar.Nframe; ++frame)
	{
		int stepits=0;
		real stept=0.0;
		while (stept<svar.framet)
		{
		    error = Newmark_Beta(TREE,svar,fvar,avar,cells,pn,pnp1,airP,outlist);
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
	    f2 << "Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;
		
		// f2 << "Minimum Coords:" << endl << std::scientific << std::setprecision(8) 
		// 	<< svar.minC(0) << " " << svar.minC(1);
		// #if SIMDIM ==3
		// f2 << " " << svar.minC(2);
		// #endif
		// f2 << endl << "Maximum Coords:" << endl  << svar.maxC(0) << " " << svar.maxC(1);
		// #if SIMDIM ==3
		// f2 <<  " " << svar.maxC(2);
		// #endif
		// f2 << endl;
		

		if(svar.outframe !=0)
		{
			if (frame % svar.outframe == 0 )
			{	/*Output to console every 20 or so steps*/
			  	cout << "Frame: " << frame << "  Sim Time: " << svar.t-svar.dt << "  Compute Time: "
			  	<< duration <<"  Error: " << error << endl;
			  	cout << "Boundary particles:  " << svar.bndPts << " Sim particles: " << svar.totPts-svar.bndPts
			  	<< " Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;
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

	return 0;
}
