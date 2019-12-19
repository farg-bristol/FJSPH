/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEWMARK_BETA_H
#define NEWMARK_BETA_H


#include <chrono>
#include "Var.h"
#include "IO.h"
#include "Neighbours.h"
#include "Resid.h"

///**************** Integration loop **************///
ldouble Newmark_Beta(Sim_Tree& NP1_INDEX, Vec_Tree& CELL_INDEX, SIM& svar, const FLUID& fvar, const AERO& avar, 
	MESH& cells, State& pn, State& pnp1, State& airP, outl& outlist)
{
	// cout << "Entered Newmark_Beta" << endl;
	const uint start = svar.bndPts;
	const uint end = svar.totPts;
	const uint piston = svar.psnPts;

	uint   k = 0;	
	double errsum = 1.0;
	double logbase = 0.0;
	double error1 = 0.0;
	double error2 = 0.0;

	/*Find maximum safe timestep*/
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});
	ldouble maxf = maxfi->f.norm();
	ldouble dtf = sqrt(fvar.H/maxf);
	ldouble dtcv = fvar.H/(fvar.Cs+svar.maxmu);

	/*Only use if -fno-finite-math is on*/
	// if (std::isinf(maxf))
	// {
	// 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
	// 	exit(-1);
	// }

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	svar.dt = 0.3*std::min(dtf,dtcv);
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/


	if (svar.dt > svar.framet)
	{
		svar.dt = svar.framet;
	}
	else
	{
		uint frac = ceil(svar.framet/svar.dt);
		svar.dt = svar.framet/frac;
	}


	/*Check if the particle has moved to a new cell*/
	if (svar.Bcase == 6)
	{
		FindCell(start,end,avar.nfull,CELL_INDEX,cells,outlist,pnp1);
	}
	// cout << svar.Start(0) << "  " << svar.Start(1) << endl;
	// cout << svar.Transp << endl;
	/*Check if a particle is running low on neighbours, and add ficticious particles*/
	std::vector<std::vector<Part>> air(start);
	if(svar.ghost == 1)
	{
		#pragma omp parallel shared(pnp1, outlist)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (uint ii = start; ii < end; ++ii)
			{
				std::vector<Part> temp;
				if(pnp1[ii].b == 2 && outlist[ii].size() > 0.4*avar.nfull 
					&& outlist[ii].size() < avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,ii,pnp1,outlist);
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
	airP.clear();

	for(uint ii = 0; ii < air.size(); ++ii)
		for(uint jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	#pragma omp parallel for shared(outlist)
	for(uint ii = start; ii < end; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}

	vector<StateVecD> xih(svar.totPts);
	const static ldouble a = 1 - svar.gamma;
	const static ldouble b = svar.gamma;
	const static ldouble c = 0.5*(1-2*svar.beta);
	const static ldouble d = svar.beta;
	const static ldouble B = fvar.B;
	const static ldouble gam = fvar.gam;

	while (log10(sqrt(errsum/(double(svar.totPts)))) - logbase > -7.0)
	{
		/****** UPDATE TREE ***********/
		NP1_INDEX.index->buildIndex();
		FindNeighbours(NP1_INDEX, fvar, pnp1, outlist);
		
		// airP.clear();
		std::vector<std::vector<Part>> neighb;
		neighb.reserve(end);
		

		#pragma omp parallel shared(pnp1, outlist, air)
		{
			std::vector<std::vector<Part>> local;
			if(svar.ghost == 1 )
			{
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
			}
			else
			{
				#pragma omp for schedule(static) nowait
				for (uint ii = 0; ii < end; ++ii)
				{
					std::vector<Part> temp;
					for(auto jj:outlist[ii])
						temp.push_back(Part(pnp1[jj])); 

					local.push_back(temp);
				}
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i<NTHREADS; i++)
	    	{
	    		#pragma omp ordered
	    		neighb.insert(neighb.end(),local.begin(),local.end());
	    	}
		}
				
		// cout << "K: " << k << endl;
		// cout << "Calculating forces" << endl;
		vector<StateVecD> res(end,StateVecD::Zero());
		vector<ldouble> Rrho(end,0.0);
 		Forces(svar,fvar,avar,pnp1,neighb,outlist,res,Rrho); /*Guess force at time n+1*/

		// #pragma omp parallel
	
		/*Previous State for error calc*/
		#pragma omp parallel for shared(pnp1)
		for (uint  ii=0; ii < end; ++ii)
			xih[ii] = pnp1[ii].xi;


		const ldouble dt = svar.dt;
		const ldouble dt2 = dt*dt;

		/*Update the state at time n+1*/
		#pragma omp parallel shared(pn,res,Rrho,svar,fvar)
		{
			#pragma omp for
			for (uint ii = 0; ii < piston; ++ii)
			{	/****** PISTON PARTICLES *************/
				pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+Rrho[ii]);
				pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
				pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
				pnp1[ii].f = res[ii];
				pnp1[ii].Rrho = Rrho[ii];
			}


			#pragma omp for 
			for (uint ii=piston; ii < start; ++ii)
			{	/****** BOUNDARY PARTICLES ***********/
				// pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+pnp1[ii].Rrho);
				// pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
				pnp1[ii].p = fvar.pPress;
				pnp1[ii].rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
				pnp1[ii].f = res[ii];
				pnp1[ii].Rrho = Rrho[ii];
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
								FirstCell(start,end,ii, CELL_INDEX, cells, outlist, pnp1);
							}
						}
					}
					else
					{	
						/*For the particles marked 1, perform a prescribed motion*/
						pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
						pnp1[ii].f = res[ii];
						pnp1[ii].Rrho = Rrho[ii];
						// pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
						// pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					}
				}

				if(pnp1[ii].b==2)
				{	/*For the particles marked 2, intergrate as normal*/
					pnp1[ii].v =  pn[ii].v +dt*(a*pn[ii].f+b*res[ii]);
					pnp1[ii].xi = pn[ii].xi+dt*pn[ii].v+dt2*(c*pn[ii].f+d*res[ii]);
					pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*Rrho[ii]);
					pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					pnp1[ii].f = res[ii];
					pnp1[ii].Rrho = Rrho[ii];
				}
			} /*End fluid particles*/
		
		/****** FIND ERROR ***********/
		errsum = 0.0;
		#pragma omp parallel for reduction(+:errsum) shared(pnp1,xih)
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

		if (error1-error2 > 0.0 /*|| std::isnan(error1)*/)
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
	// cout << "Out of the parallel loop." << endl;
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

#endif