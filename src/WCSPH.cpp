/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/***			        + XSPH Correction (Also described in Monaghan)                    ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)                                  ***/
/*** Density Reinitialisation: Colagrossi, A. and Landrini, M. (2003): MLS                ***/
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

// using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace nanoflann;

///**************** Integration loop **************///
ldouble Newmark_Beta(Sim_Tree& NP1_INDEX, SIM& svar, const FLUID& fvar, const AERO& avar, 
	const MESH& cells, State& pn, State& pnp1/*, State& airP*/, outl& outlist)
{
	// cout << "Entered Newmark_Beta" << endl;
	const uint start = svar.bndPts;
	const uint end = svar.totPts;
	const uint piston = svar.psnPts;
		
	double errsum = 1.0;
	double logbase = 0.0;
	unsigned int k = 0;
	double error1 = 0.0;
	double error2 = 0.0;

	/*Find maximum safe timestep*/
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});
	ldouble maxf = maxfi->f.norm();
	ldouble dtf = sqrt(fvar.H/maxf);
	ldouble dtcv = fvar.H/(fvar.Cs+svar.maxmu);


	// if (std::isinf(maxf))
	// {
	// 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
	// 	exit(-1);
	// }

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	svar.dt = 0.3*min(dtf,dtcv)/* dtf*/;
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

	if (svar.dt > svar.framet)
		svar.dt = svar.framet;


	/*Check if the particle has moved to a new cell*/
	if (svar.Bcase == 6)
	{
		FindCell(start,end,avar.nfull,pnp1,cells,outlist);
	}
	// cout << svar.Start(0) << "  " << svar.Start(1) << endl;
	// cout << svar.Transp << endl;
	/*Check if a particle is running low on neighbours, and add ficticious particles*/
	std::vector<std::vector<Part>> air(start);
	if(svar.ghost == 1 )
	{
		#pragma omp parallel shared(pnp1, outlist)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (uint ii = start; ii < end; ++ii)
			{
				StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
				std::vector<Part> temp;
				if(vec(1) > 0.0 && pnp1[ii].b == 2 &&
				outlist[ii].size() > 0.4*avar.nfull && outlist[ii].size() < avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,ii,pnp1,outlist);
				}
				local.emplace_back(temp);
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i<NTHREADS; i++)
	    	{
	    		#pragma omp ordered
	    		air.insert(air.end(),local.begin(),local.end());
	    	}
		}
	}
	

	#pragma omp parallel for shared(outlist)
	for(uint ii = svar.bndPts; ii < svar.totPts; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}

	vector<StateVecD> xih(svar.totPts);
	const ldouble a = 1 - svar.gamma;
	const ldouble b = svar.gamma;
	const ldouble c = 0.5*(1-2*svar.beta);
	const ldouble d = svar.beta;
	const ldouble B = fvar.B;
	const ldouble gam = fvar.gam;

	// int RestartCount = 0;
	while (log10(sqrt(errsum/(double(svar.totPts)))) - logbase > -7.0)
	{
		/****** UPDATE TREE ***********/
		NP1_INDEX.index->buildIndex();
		FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);
		
		// airP.clear();
		// cout << "Creating neighb list" << endl;
		std::vector<std::vector<Part>> neighb;
		neighb.reserve(end);
		// for(uint ii = 0; ii < start; ++ii)
		// {
		// 	neighb.emplace_back();
		// }

		#pragma omp parallel shared(pnp1, outlist, air)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (uint ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				if(air[ii].size()>0)
				{
					temp.insert(temp.end(), air[ii].begin(), air[ii].end());
				}
				local.push_back(temp);
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i<NTHREADS; i++)
	    	{
	    		#pragma omp ordered
	    		neighb.insert(neighb.end(),local.begin(),local.end());
	    	}
		}
		

		// for(uint ii = start; ii < end; ++ii)
		// {
		// 	std::vector<Part> temp;
		// 	for(auto j:outlist[ii])
		// 			temp.emplace_back(Part(pnp1[j])); 

		// 	neighb.emplace_back(temp);
		// }
		
		// cout << "K: " << k << endl;
		// cout << "Calculating forces" << endl;
 		Forces(svar,fvar,avar,pnp1,neighb,outlist); /*Guess force at time n+1*/

		// #pragma omp parallel
	
		/*Previous State for error calc*/
		#pragma omp parallel for shared(pnp1)
		for (uint  ii=0; ii < svar.totPts; ++ii)
			xih[ii] = pnp1[ii].xi;


		const ldouble dt = svar.dt;
		const ldouble dt2 = dt*dt;

		/*Update the state at time n+1*/
		#pragma omp parallel shared(pn)
		{
			// #pragma omp for
			// for (uint ii = 0; ii < piston; ++ii)
			// {	/****** PISTON PARTICLES *************/
			// 	pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+pnp1[ii].Rrho);
			// 	pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
			// 	pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
			// }


			#pragma omp for 
			for (uint ii=piston; ii < start; ++ii)
			{	/****** BOUNDARY PARTICLES ***********/
				pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+pnp1[ii].Rrho);
				pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
			}
			

			#pragma omp for
			for (uint ii=start; ii < end; ++ii)
			{	/****** FLUID PARTICLES **************/
				
				/*Check if the particle is clear of the starting area*/
				if(pnp1[ii].b == 1)
				{   /* boundary particle value of 1 means its a                    */
					/* fluid particle that isn't clear of the starting area        */
					/* 2 means it is clear, and can receive a force aerodynamically*/
					/* 3 is an air particle. Stored only for one timestep at a time*/
					StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
					if(vec(1) > svar.clear)
					{	/*Tag it as clear if it's higher than the plane of the exit*/
						pnp1[ii].b=2;
						if(svar.Bcase == 6)
						{
							/*Search through the cells to find which cell it's in*/
							if (outlist[ii].size() < avar.nfull)
							{
								#if SIMDIM == 3
									uint found = 0;
									uint jj = 0;
									StateVecD testp = pnp1[ii].xi;
									// cout << cells.cFaces.size() << endl;
									for(std::vector<std::vector<StateVecD>> const& cell:cells.cFaces)
									{
										if(Crossings3D(cell,testp))
										{
											pnp1[ii].cellID = jj;
											pnp1[ii].cellV = cells.cVel[jj];
											pnp1[ii].cellP = cells.cellP[jj];
											pnp1[ii].cellRho = cells.cellRho[jj];
											found = 1;
											// cout << "Cell found!" << endl;
											break;
											
										}
										++jj;
									}

									if(found != 1)
									{
										cout << "Containing cell not found. Something is wrong." << endl;
										exit(-1);
									}	
								#else
									uint found = 0;
									uint jj = 0;
									StateVecD testp = pnp1[ii].xi;
									/*Do a cell containment*/
									for(std::vector<StateVecD> const& cell:cells.cVerts)
									{
										if(Crossings2D(cell,testp))
										{
											pnp1[ii].cellID = jj;
											pnp1[ii].cellV = cells.cVel[jj];
											pnp1[ii].cellP = cells.cellP[jj];
											pnp1[ii].cellRho = cells.cellRho[jj];
											found = 1;
											// cout << "Found the containing cell!" << endl;
											// cout << jj << "  " << cells.cVel[jj][0] << endl;
											break;
											
										}
										++jj;
									}

									if(found != 1)
									{
										cout << "Containing cell not found. Something is wrong." << endl;
										exit(-1);
									}
								#endif
							}
						}
					}
					else
					{	
						/*For the particles marked 1, perform a prescribed motion*/
						pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
						pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
						pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					}
				}

				if(pnp1[ii].b==2)
				{	/*For the particles marked 2, intergrate as normal*/
					pnp1[ii].v =  pn[ii].v +dt*(a*pn[ii].f+b*pnp1[ii].f);
					pnp1[ii].xi = pn[ii].xi+dt*pn[ii].v+dt2*(c*pn[ii].f+d*pnp1[ii].f);
					pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
					pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
				}
			} /*End fluid particles*/
		
		/****** FIND ERROR ***********/
		errsum = 0.0;
		#pragma omp parallel for reduction(+:errsum)
		for (uint ii = start; ii < end; ++ii)
		{
			StateVecD r = pnp1[ii].xi-xih[ii];
			errsum += r.squaredNorm();
		}
		
		}/*End pragma omp parallel*/

		if(k == 0)
			logbase=log10(sqrt(errsum/(double(svar.totPts))));


		error1 = log10(sqrt(errsum/(double(svar.totPts)))) - logbase;
		// cout << RestartCount << "  " << k << "  " << error1  << "  " << svar.dt << endl;
		// cout << k << "  " << error1 << "  " << svar.dt << endl;

		if (error1-error2 > 0.0 || std::isnan(error1))
		{	/*If simulation starts diverging, then reduce the timestep and try again.*/
			// cout << "Unstable timestep. Reducing timestep..." << endl;
			pnp1 = pn;
			svar.dt = 0.5*svar.dt;
			k = 0;
			error1 = 0.0;
			// RestartCount++;
		}	/*Check if we've exceeded the maximum iteration count*/
		else if (k > svar.subits)
			break;
		else
		{	/*Otherwise, roll forwards*/
			++k;
		}

		error2 = error1;
		
	} /*End of subits*/

	// RestartCount = 0;
	/*Add time to global*/
	svar.t+=svar.dt;

	// cout << "Timestep Params: " << maxf << " " << fvar.Cs + maxmu << " " << dtf << " " << dtcv << endl;
	// cout << "New Time: " << svar.t << endl;
	

	/*Check if more particles need to be created*/
	if(svar.Bcase >= 2 && svar.Bcase !=5)
	{
		switch(svar.Bclosed)
		{
			case 0:
			{
				int refresh = 1;
				for (uint ii = svar.totPts-svar.nrefresh; ii<svar.totPts; ++ii)
				{	/*Check that the starting area is clear first...*/
					StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
					if(svar.Bcase == 2)
					{
						if(vec[1]< svar.dx-svar.Jet(1)*3)
							refresh = 0;
					}
					else
					{
						if(vec[1]< svar.dx-svar.Jet(1))
							refresh = 0;
					}
				}

				if(refresh == 1)
				{	/*...If it is, then check if we've exceeded the max points...*/
					if (svar.addcount < svar.nmax)
					{	/*...If we havent, then add points. */
						// cout << "adding points..." << endl;
						StateVecD vec = svar.Transp*(pnp1[svar.totPts-1].xi-svar.Start);
						AddPoints(vec[1]-svar.dx, svar, fvar, avar, pn, pnp1);
					}
					else
					{	/*...If we have, then check we're sufficiently
						clear to close the boundary*/
						cout << "End of adding particle rounds." << endl;
						svar.Bclosed = 1;
						
					}
					NP1_INDEX.index->buildIndex();
					FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);
				}
				break;
			}
			case 1:
				break;
		}
	}

	/****** UPDATE TIME N ***********/
	if(svar.totPts != pnp1.size())
	{
		cout << "Size mismatch. Total points not equal to array size. Stopping" << endl;
		exit(-1);
	}

	pn = pnp1;

	return log10(sqrt(errsum/(double(svar.totPts))))-logbase;
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
		string meshfile = svar.infolder;
		meshfile.append(svar.meshfile);
		Read_TAUMESH(meshfile,cells,fvar);
	}

	/*Make a guess of how many there will be...*/
	int partCount = ParticleCount(svar);
    ///****** Initialise the particles memory *********/
	State pn;	    /*Particles at n   */
	State pnp1; 	/*Particles at n+1 */

	cout << "Final particle count:  " << partCount << endl;
	svar.finPts = partCount;
	pn.reserve(partCount);
  	pnp1.reserve(partCount);
	
	MakeOutputDir(argc,argv,svar);

	if (svar.Bcase == 6)
		Write_Mesh_Data(svar,cells);

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
		avar.nfull = 32.67;
		svar.nfull = 48;
	#endif

	
	// cout << avar.nfull << endl;

	/*Check if a particle is running low on neighbours, and add ficticious particles*/
	const uint start = svar.bndPts;
	const uint end = svar.totPts;

	std::vector<std::vector<Part>> neighb;
	neighb.reserve(end);
	for(uint ii = 0; ii < start; ++ii)
		neighb.emplace_back();

	#pragma omp parallel shared(svar, pnp1, outlist)
	{
		std::vector<std::vector<Part>> local;
		#pragma omp for schedule(static) nowait 
		for (uint ii = start; ii < end; ++ii)
		{
			std::vector<Part> temp;
			if(svar.ghost == 1 && pnp1[ii].b == 2 && outlist[ii].size() < avar.nfull)
				temp = PoissonSample::generatePoissonPoints(svar,fvar,ii,pnp1,outlist);

			for(auto j:outlist[ii])
				temp.emplace_back(Part(pnp1[j]));

			local.emplace_back(temp);
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<NTHREADS; i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),local.begin(),local.end());
    	}
	}
	

	// for(uint ii = start; ii < end; ++ii)
	// {
	// 	std::vector<Part> temp;
	// 	for(auto j:outlist[ii])
	// 			temp.emplace_back(Part(pnp1[j])); 

	// 	neighb.emplace_back(temp);
	// }

	#pragma omp parallel for shared(outlist)
	for(uint ii = svar.bndPts; ii < svar.totPts; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}


	Forces(svar,fvar,avar,pnp1,neighb,outlist);

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

	f1 << std::scientific << setprecision(6);
	f2 << std::scientific<< setw(10);
	f3 << std::scientific << setw(10);

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
					Write_ASCII_header(fb,svar);
					Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts);
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
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts);
			if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts);

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
		    error = Newmark_Beta(NP1_INDEX,svar,fvar,avar,cells,pn,pnp1,outlist);
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
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts);
			if(svar.Bcase != 0 && svar.Bcase != 5 && svar.boutform == 1)
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts);

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
