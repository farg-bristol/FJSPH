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
    	for(int i=0; i<NTHREADS; i++)
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
	
	vector<StateVecD> res(svar.totPts,StateVecD::Zero());
	vector<ldouble> Rrho(svar.totPts,0.0);
	Forces(svar,fvar,avar,pnp1,neighb,outlist,res,Rrho); /*Guess force at time n+1*/

	#pragma omp parallel for shared(res, Rrho)
	for(uint ii = 0; ii < end; ++ii)
	{
		pnp1[ii].f = res[ii];
		pnp1[ii].Rrho = Rrho[ii]; 
	}

	#if DEBUG 
		dbout << "Exiting first step." << endl;
	#endif
}


int main(int argc, char *argv[])
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
	Eigen::initParallel();
	omp_set_num_threads(NTHREADS);

    double duration;
    double error = 0;
    cout.width(13);
    cout << std::scientific << std::left;
    write_header();

    /******* Define the global simulation parameters ******/
	SIM svar;
	FLUID fvar;
	AERO avar;
	outl outlist;
	MESH cells;

	GetInput(argc,argv,svar,fvar,avar);
	
	if(svar.Bcase == 6)
	{
		#if SIMDIM == 3
		Read_TAUMESH(svar,cells,fvar);
		#else
		string meshfile = svar.infolder;
		meshfile.append(svar.meshfile);
		Read_TAUPLT(meshfile,cells,fvar);
		#endif
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
	
	MakeOutputDir(argc,argv,svar);
	

	// if (svar.Bcase == 6){
	// 	Write_Mesh_Data(svar,cells);
	// }	

	InitSPH(svar,fvar,avar,pn,pnp1);

	// Check if cells have been initialsed before making a tree off it
	if(cells.cCentre.size() == 0)
		cells.cCentre.emplace_back(StateVecD::Zero());


	///********* Tree algorithm stuff ************/
	Vec_Tree CELL_INDEX(SIMDIM,cells.cCentre,20);
	Sim_Tree NP1_INDEX(SIMDIM,pnp1,20);
	CELL_INDEX.index->buildIndex();
	NP1_INDEX.index->buildIndex();

	FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);
	if(svar.Bcase == 6)
	{
		cout << "Building cell neighbours..." << endl;
		FindCellNeighbours(CELL_INDEX, cells.cCentre, cells.cNeighb);
	}
		// cells.Build_cPolys();
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
		avar.nfull = 32.67;
		svar.nfull = 48;
	#endif

	
	// cout << avar.nfull << endl;
	First_Step(svar,fvar,avar,outlist,pnp1,airP);

	///*************** Open simulation files ***************/
	std::ofstream f1,f2,f3,fb;


	if (svar.frameout == 1)
	{
		string framef = svar.outfolder;
		framef.append("/frame.info");
		f2.open(framef, std::ios::out);
	}

	if(svar.frameout ==2)
		f3.open("Crossflow.txt", std::ios::out);

	f1 << std::scientific << std::setprecision(6);
	f2 << std::scientific << std::setw(10);
	f3 << std::scientific << std::setw(10);

	if(svar.outtype == 0 )
	{
		if(svar.outform < 3)
		{	
			if (svar.Bcase != 0 || svar.Bcase !=5)
				Write_Boundary_Binary(svar,pnp1);
			
			Init_Binary_PLT(svar);
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,2); /*Write sim particles*/
		}
		else if (svar.outform > 2)
		{
			cout << "Output type not within design. Outputting basic data..." << endl;
			svar.outform = 0;
			if (svar.Bcase != 0 || svar.Bcase !=5)
				Write_Boundary_Binary(svar,pnp1);

			Init_Binary_PLT(svar);
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,2); /*Write sim particles*/
		}
	}
	else if (svar.outtype == 1)
	{

		if (svar.Bcase != 0 && svar.Bcase != 5)
		{	/*If the boundary exists, write it.*/
			string bfile = svar.outfolder;
			bfile.append("/Boundary.plt");
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
					Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,empty);
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
		mainfile.append("/Fuel.plt");
		f1.open(mainfile, std::ios::out);
		if(f1.is_open())
		{
			Write_ASCII_header(f1,svar);
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,airP);
			if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,airP);

		}
		else
		{
			cerr << "Failed to open fuel.plt. Stopping." << endl;
			exit(-1);
		}
	}
	else
	{
		cerr << "Output type ambiguous. Please select 0 or 1 for output data type." << endl;
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
	f2 << "Frame:  B-Points: S-Points:  Sim Time:       Comp Time:     Error:       Its:" << endl;
	f2 << 0 << "        " << svar.bndPts << "  " << svar.simPts << "    " << svar.t << "    " << duration
		<< "    " << error << "  " << 0 << endl;

	///************************* MAIN LOOP ********************/
	svar.frame = 0;
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

		t2= high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;

		/*Write each frame info to file (Useful to debug for example)*/
		if (svar.frameout == 1)
		{
			f2 << frame << "        " << svar.bndPts << "  " << svar.simPts << "    " << svar.t << "    " << duration
				<< "    " << error << "  " << stepits << endl;
			
		}

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
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,2); /*Write sim particles*/
		} 
		else if (svar.outtype == 1)
		{
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,airP);
			if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,airP);

		}
	}
	if (f1.is_open())
		f1.close();
	if (f2.is_open())
		f2.close();
	if (f3.is_open())
		f3.close();
	if (fb.is_open())
		fb.close();
			
	if(svar.outtype == 0)
	{
		if(TECEND142())
			exit(-1);
	}
	
	cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.t << " seconds" << endl;
	return 0;
}
