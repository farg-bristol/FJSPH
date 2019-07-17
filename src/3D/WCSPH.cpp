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
#include "Kernel.h"
#include "Init.h"
#include "Cross.h"
#include "Neighbours.h"
#include "Resid.h"

// using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace nanoflann;


///**************** Integration loop **************///
ldouble Newmark_Beta(Sim_Tree &NP1_INDEX, SIM &svar, FLUID &fvar, CROSS &cvar,
	 State &pn, State &pnp1, outl &outlist)
{

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


	if (std::isinf(maxf))
	{
		std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
		exit(-1);
	}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	svar.dt = 0.3*min(dtf,dtcv)/* dtf*/;
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

	if (svar.dt > svar.framet)
		svar.dt = svar.framet;

	vector<StateVecD> xih(svar.totPts);
	// int RestartCount = 0;
	while (log10(sqrt(errsum/(double(svar.totPts)))) - logbase > -7.0)
	{
		if(svar.totPts != pnp1.size())
		{
			cout << "Size mismatch. Total points not equal to array size. Stopping" << endl;
			exit(-1);
		}
		/****** UPDATE TREE ***********/
		// cout << "Timestep: " << svar.dt << " Sim Time: " << svar.t << " Error: " << error1 << endl;
		NP1_INDEX.index->buildIndex();
		FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);

		// cout << "K: " << k << endl;
		Forces(svar,fvar,cvar,pnp1,outlist); /*Guess force at time n+1*/

		// #pragma omp parallel
	
		/*Previous State for error calc*/
		// #pragma omp parallel for
		for (uint  i=0; i< svar.totPts; ++i)
			xih[i] = pnp1[i].xi;

		ldouble a = 1 - svar.gamma;
		ldouble b = svar.gamma;
		ldouble c = 1-2*svar.beta;
		ldouble d = svar.beta;
		ldouble dt = svar.dt;
		ldouble dt2 = dt*dt;
		ldouble B = fvar.B;
		ldouble gam = fvar.gam;

		/*Update the state at time n+1*/
		// #pragma omp parallel shared(a,b,c,d,dt,dt2,B,gam) num_threads(nthreads)
		{
			// #pragma omp for 
			for (uint i=0; i <= svar.totPts; ++i)
			{	/****** BOUNDARY PARTICLES ***********/
				#pragma omp critical
				{
				pnp1[i].rho = pn[i].rho+0.5*dt*(pn[i].Rrho+pnp1[i].Rrho);
				pnp1[i].p = B*(pow(pnp1[i].rho/fvar.rho0,gam)-1);
				}
			}
			

			// #pragma omp for
			for (uint i=svar.bndPts; i < svar.totPts ; ++i)
			{	/****** FLUID PARTICLES ***********/
				#pragma omp critical
				{
				pnp1[i].v =  pn[i].v +dt*(a*pn[i].f+b*pnp1[i].f);
				pnp1[i].xi = pn[i].xi+dt*pn[i].v+dt2*(0.5*c*pn[i].f+d*pnp1[i].f);
				}
				// pnp1[i].rho = pn[i].rho+dt*(a*pn[i].Rrho+b*pnp1[i].Rrho);
				// pnp1[i].p = fvar.B*(pow(pnp1[i].rho/fvar.rho0,fvar.gam)-1);
			}
		}

			/****** FIND ERROR ***********/
			errsum = 0.0;
			// #pragma omp parallel for reduction(+:errsum)
			for (uint i=svar.bndPts; i < svar.totPts; ++i)
			{
				StateVecD r = pnp1[i].xi-xih[i];
				errsum += r.squaredNorm();
			}
		

		if(k == 0)
			logbase=log10(sqrt(errsum/(double(svar.totPts))));


		error1 = log10(sqrt(errsum/(double(svar.totPts)))) - logbase;
		// cout << RestartCount << "  " << k << "  " << error1  << "  " << svar.dt << endl;

		if (error1-error2 > 0.0)
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
	}

	// RestartCount = 0;
	/*Add time to global*/
	svar.t+=svar.dt;

// cout << "Timestep Params: " << maxf << " " << fvar.Cs + maxmu << " " << dtf << " " << dtcv << endl;
// cout << "New Timestep: " << svar.dt << endl;

	/*Check if more particles need to be created*/
	if(svar.Bcase >= 3 && svar.Bcase !=5)
	{
		switch(svar.Bclosed)
		{
			case 0:
			{
				int refresh = 1;
				for (uint i = svar.totPts-svar.nrefresh; i<svar.totPts; ++i)
				{	/*Check that the starting area is clear first...*/
					StateVecD vec = svar.Rotate.transpose()*(pn[i].xi-svar.Start);
					if(vec[1]< svar.dx-svar.Jet(1))
						refresh = 0;
				}

				if(refresh == 1)
				{	/*...If it is, then check if we've exceeded the max points...*/
					if (svar.addcount < svar.nmax)
					{	/*...If we havent, then add points. */
						// cout << "adding points..." << endl;
						StateVecD vec = svar.Rotate.transpose()*(pn[svar.totPts-1].xi-svar.Start);
						AddPoints(vec[1]-svar.dx, svar, fvar, cvar, pn, pnp1);
					}
					else
					{	/*...If we have, then check we're sufficiently
						clear to close the boundary*/
						cout << "End of adding particle rounds. Stopping..." << endl;
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
	pn = pnp1;

	return log10(sqrt(errsum/(double(svar.totPts))))-logbase;
}

int main(int argc, char *argv[])
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
	Eigen::initParallel();
	omp_set_num_threads(nthreads);

    double duration;
    double error = 0;
    cout << setw(5);

    write_header();

    /******* Define the global simulation parameters ******/
	SIM svar;
	FLUID fvar;
	CROSS cvar;
	outl outlist;

	GetInput(argc,argv,svar,fvar,cvar);

	/*Make a guess of how many there will be...*/
	int partCount = ParticleCount(svar);
    ///****** Initialise the particles memory *********/
	State pn;	    /*Particles at n   */
	State pnp1; 	/*Particles at n+1 */

	cout << "Final particle count:  " << partCount << endl;
	pn.reserve(partCount);
  	pnp1.reserve(partCount);

  	std::ofstream f1;
	MakeOutputDir(argc,argv,svar,f1);
	
	InitSPH(svar,fvar,cvar, pn, pnp1);

	///********* Tree algorithm stuff ************/
	Sim_Tree NP1_INDEX(simDim,pnp1,20);
	NP1_INDEX.index->buildIndex();
	FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);

	std::vector<std::vector<uint>>::iterator nfull = 
		std::max_element(outlist.begin(),outlist.end(),
		[](std::vector<uint> p1, std::vector<uint> p2){return p1.size()< p2.size();});

	fvar.avar.nfull = double(nfull->size());
	// cout << fvar.avar.nfull << endl;
	///*** Perform an iteration to populate the vectors *****/
	Forces(svar,fvar,cvar,pnp1, outlist);


	///*************** Open simulation files ***************/
	std::ofstream f2,f3;
	if (svar.frameout == 1)
		f2.open("frame.info", std::ios::out);

	if(svar.outform ==2)
		f3.open("Crossflow.txt", std::ios::out);


	if (f1.is_open())
	{
		f1 << std::scientific << setprecision(5);
		f2 << std::scientific<< setw(10);
		f3 << std::scientific << setw(10);


		/* Write file header defining variable names */
		write_file_header(f1,svar,pnp1);
		if(svar.Bcase == 4)
			svar.vortex.write_VLM_Panels(svar.outfolder);

		/*Timing calculation + error sum output*/
		t2 = high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;
		cout << "Frame: " << 0 << "  Sim Time: " << svar.t << "  Compute Time: "
		<< duration <<"  Error: " << error << endl;
		f2 << "Frame:  Points:   Sim Time:       Comp Time:     Error:       Its:" << endl;
		f2 << 0 << "        " << svar.totPts << "    " << svar.t << "    " << duration
			<< "    " << error << "  " << 0 << endl;

		///************************* MAIN LOOP ********************/

		for (uint frame = 1; frame<= svar.Nframe; ++frame)
		{
			int stepits=0;
			double stept=0.0;
			while (stept<svar.framet)
			{
			    error = Newmark_Beta(NP1_INDEX,svar,fvar,cvar,pn,pnp1,outlist);
			    stept+=svar.dt;
			    ++stepits;
			    //cout << svar.t << "  " << svar.dt << endl;
			}


			t2= high_resolution_clock::now();
			duration = duration_cast<microseconds>(t2-t1).count()/1e6;

			/*Write each frame info to file (Useful to debug for example)*/
			if (svar.frameout == 1)
			{
				f2 << frame << "         " << svar.totPts << "  " << svar.t << "  " << duration
				<< "  " << error << "  " << stepits << endl;
			}

			if(svar.outframe !=0)
			{
				if (frame % svar.outframe == 0 )
				{	/*Output to console every 20 or so steps*/
				  	cout << "Frame: " << frame << "  Sim Time: " << svar.t-svar.dt << "  Compute Time: "
				  	<< duration <<"  Error: " << error << endl;
				}
			}


			DensityReinit(fvar, pnp1, outlist);
			switch (svar.outform)
			{
				case 1:
					write_fluid_data(f1, svar, pnp1);
					break;
				case 2:
					write_research_data(f1, svar, pnp1);
					break;
				case 3:
					write_basic_data(f1, svar, pnp1);
					break;
			}

		}
		f1.close();
		if (f2.is_open())
			f2.close();
		if(f3.is_open())
			f3.close();
	}
	else
	{
		cerr << "Error opening frame output file." << endl;
		exit(-1);
	}

	cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.t << " seconds" << endl;
	return 0;
}
