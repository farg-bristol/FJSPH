/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEWMARK_BETA_H
#define NEWMARK_BETA_H

#include <chrono>
#include <unordered_set>
#include "Var.h"
#include "Resid.h"
#include "Crossing.h"

int Check_Error(KDTREE& TREE, SIM& svar, FLUID const& fvar, size_t const& start, size_t const& end, 
		real& error1, real& error2, real& logbase, 
		vector<size_t>& cellsused, outl& outlist, vector<StateVecD> const& xih,
		State& pn, State& pnp1, uint& k, uint& nUnstab)
{
	/****** FIND ERROR ***********/
	real errsum = 0.0;
	#pragma omp parallel for reduction(+:errsum) schedule(static) shared(pnp1,xih)
	for (size_t ii = start; ii < end; ++ii)
	{
		StateVecD r = pnp1[ii].xi-xih[ii-start];
		errsum += r.squaredNorm();
	}

	if(k == 0)
		logbase=log10(sqrt(errsum/(real(svar.totPts))));

	error1 = log10(sqrt(errsum/(real(svar.totPts)))) - logbase;
	// cout << RestartCount << "  " << k << "  " << error1  << "  " << svar.dt << endl;
	// cout << k << "  " << error1 << "  " << error1 - error2 << "  " << svar.dt << endl;

	if (k > svar.subits)
	{
	// if (error1-error2 > 0.0 /*|| std::isnan(error1)*/)
	// {	/*If simulation starts diverging, then reduce the timestep and try again.*/
		// nUnstab++;

		// if(nUnstab > 4)
		if(error1 > 0.0)
		{
			pnp1 = pn;
			
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

			if(svar.Asource == 2)
			{
				cellsused.clear();
				#pragma omp parallel shared(pnp1)
				{
					std::vector<size_t> local;
					#pragma omp for schedule(static) nowait
					for (size_t ii = start; ii < end; ++ii)
					{	
						if(pnp1[ii].b == PartState.FREE_)
							local.emplace_back(pnp1[ii].cellID);
					}

					#pragma omp for schedule(static) ordered
					for(int i=0; i<omp_get_num_threads(); i++)
					{
						#pragma omp ordered
						cellsused.insert(cellsused.end(),local.begin(),local.end());
					}
				}

				// Sort the vector and remove unique values.
				std::unordered_set<size_t> s;
				for(size_t const& ii:cellsused)
					s.insert(ii);

				cellsused.assign(s.begin(),s.end());
				std::sort(cellsused.begin(),cellsused.end());
			}

			svar.dt = 0.5*svar.dt;
			cout << "Unstable timestep. New dt: " << svar.dt << endl;
			k = 0;
			error1 = 0.0;
			// RestartCount++;
			return -1;
		}

	// 	return 0;
	// }	/*Check if we've exceeded the maximum iteration count*/
	
		return 1;
	}

	/*Otherwise, roll forwards*/
	nUnstab = 0;
	return 0;
	
}

void Do_NB_Iter(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, 
	real const& a, real const& b, real const& c, real const& d, real const& B, real const& gam,
	MESH& cells, vector<size_t> const& cellsused, outl const& outlist, DELTAP const& dp,
	State& pn, State& pnp1, StateVecD& Force, StateVecD& dropVel)
{
	/*Update the state at time n+1*/
	
	// cout << "Calculating forces" << endl;
	vector<StateVecD> res;
	vector<StateVecD> Af;
	vector<real> Rrho;
	vector<real> curve;

	Force = StateVecD::Zero();
	svar.AForce = StateVecD::Zero();

	Get_Boundary_Pressure(fvar,start,outlist,pnp1);

	Forces(svar,fvar,avar,cells,pnp1,outlist,dp,res,Rrho,Af,Force); /*Guess force at time n+1*/

	#pragma omp parallel shared(pn,svar,fvar,res,Af,Rrho,curve) /*reduction(+:Force,dropVel)*/
	{
		const real dt = svar.dt;
		const real dt2 = dt*dt;

		if(svar.Asource == 2)
		{
			#pragma omp for schedule(static) /*nowait*/
			for(size_t const& ii : cellsused)
			{
				cells.fNum[ii] = 0;
				cells.fMass[ii] = 0.0;
				cells.vFn[ii] = StateVecD::Zero();
				cells.vFnp1[ii] = StateVecD::Zero();
			}
		}

		// #pragma omp for schedule(static) nowait
		// for (size_t ii=0; ii < start; ++ii)
		// {	/****** BOUNDARY PARTICLES ***********/
		// 	pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+Rrho[ii]);
		// 	pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
		// 	// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
		// 	// pnp1[ii].p = fvar.pPress;
		// 	// pnp1[ii].rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
		// 	// pnp1[ii].f = res[ii];
		// 	pnp1[ii].Rrho = Rrho[ii];
		// }


		#pragma omp for nowait 
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* BUFFER = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			

			if(pnp1[ii].b > PartState.BUFFER_)
			{	
				StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.sim_start);

				if(pnp1[ii].b == PartState.PIPE_)
				{	/*Do a check to see if it needs to be given an aero force*/	
					if(vec(1) > 0.0)
					{
						pnp1[ii].b = PartState.FREE_;
						if(dp.lam_ng[ii] < 0.75 && (svar.Asource == 1 || svar.Asource == 2))
						{
							/*retrieve the cell it's in*/
							uint to_del = 0;
							FirstCell(svar,TREE.CELL, cells, pnp1[ii], to_del);

							if(to_del)
							{
								#pragma omp critical
								{
								// cout << "Particle has crossed an outer boundary!" << endl;
								// cout << "Particle will be deleted." << endl;
								pnp1.erase(pnp1.begin()+ii);
								pn.erase(pn.begin()+ii);
								}
								#pragma omp atomic
									svar.totPts--;
								#pragma omp atomic
									end--;
							}
						}
					}	
				}

				/*For any other particles, intergrate as normal*/
				#ifdef NOALE
				pnp1[ii].xi = pn[ii].xi+dt*pn[ii].v+dt2*(c*pn[ii].f+d*res[ii]);
				#else 
				pnp1[ii].xi = pn[ii].xi+dt*(pn[ii].v + pnp1[ii].vPert)+dt2*(c*pn[ii].f+d*res[ii]);
				#endif
				pnp1[ii].v =  (pn[ii].v/* + pnp1[ii].vPert*/) +dt*(a*pn[ii].f+b*res[ii]);
				pnp1[ii].f = res[ii];
				pnp1[ii].Af = Af[ii];
				pnp1[ii].Rrho = Rrho[ii];

				pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*Rrho[ii]);
				pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
				// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);

				Force += res[ii]*pnp1[ii].m;
				dropVel += (pnp1[ii].v + pnp1[ii].vPert);
				
				// pnp1[ii].theta = dp.lam[ii];
				// pnp1[ii].nNeigb = real(outlist[ii].size());
				// pnp1[ii].curve = curve[ii];
				pnp1[ii].pDist = vec.norm();

				if(svar.Asource == 2 && pnp1[ii].b == PartState.FREE_)
				{
					#pragma omp atomic
						cells.fNum[pnp1[ii].cellID]++;
					#pragma omp atomic
						cells.fMass[pnp1[ii].cellID] += pnp1[ii].m;

					#pragma omp critical
					{
						cells.vFn[pnp1[ii].cellID] += pn[ii].v;
						cells.vFnp1[pnp1[ii].cellID] += pnp1[ii].v;
					}	
				}

				
			}
			// else
			// {
			// 	pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;

			// 	// pnp1[pi].v = (pn[ii].v /* + pnp1[ii].vPert*/) + dt * (a * pn[ii].f + b * res[ii]);
			// 	// pnp1[ii].f = res[ii];
			// 	// pnp1[ii].Af = Af[ii];
			// 	pnp1[ii].Rrho = Rrho[ii];

			// 	pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*Rrho[ii]);
			// 	pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
			// }
			pnp1[ii].s = dp.lam_ng[ii];


			
		} /*End fluid particles*/

		/* Do the buffer particles */
		#pragma omp for schedule(static)
		for(size_t ii = 0; ii < svar.back.size(); ++ii)
		{	
			// size_t const& pID = svar.back[ii];
			// StateVecD vec = svar.Transp*(pnp1[pID].xi-svar.Start);
			// StateVecD v = pnp1[pID].v;
			for(size_t jj = 0; jj < 4; ++jj)
			{	/* Define buffer off the back particle */
				// StateVecD xi = vec;
				// xi[1] -= real(jj+1.0)*svar.dx;
				// xi = svar.Rotate*xi + svar.Start;

				size_t const& buffID = svar.buffer[ii][jj];

				// pnp1[buffID].xi = xi;
				// pnp1[buffID].v = v;
				// pnp1[buffID].rho = fvar.rhoJ;
				// pnp1[buffID].p = fvar.pPress;	
				// pnp1[buffID].Rrho = Rrho[buffID];
				// pnp1[buffID].rho = pn[buffID].rho+dt*
				// 	(a*pn[buffID].Rrho+b*Rrho[buffID]);
				// pnp1[buffID].p = B*(pow(pnp1[buffID].rho/fvar.rho0,gam)-1);

				real frac = std::min(1.0,real(jj+1)/3.0);
				// #ifdef NOALE
				// pnp1[buffID].xi = frac*xi + (1.0-frac)*(pn[buffID].xi+dt*pn[buffID].v+dt2*(c*pn[buffID].f+d*res[buffID]));
				// #else 
				// pnp1[buffID].xi = frac*xi + (1.0-frac)*(pn[buffID].xi+dt*(pn[buffID].v + pnp1[buffID].vPert)+dt2*(c*pn[buffID].f+d*res[ii]));
				// #endif
				// pnp1[buffID].v =  frac*v + (1.0-frac)*((pn[buffID].v/* + pnp1[buffID].vPert*/) +dt*(a*pn[buffID].f+b*res[buffID]));
				// pnp1[buffID].f = res[buffID];
				// // pnp1[buffID].Af = Af[buffID];
				pnp1[buffID].Rrho = Rrho[buffID];

 				pnp1[buffID].rho = frac * fvar.rhoJ + (1.0-frac) * (pn[buffID].rho+dt*(a*pn[buffID].Rrho+b*Rrho[buffID]));
				pnp1[buffID].p = frac * fvar.pPress + (1.0-frac) * (B*(pow(pnp1[buffID].rho/fvar.rho0,gam)-1));

				pnp1[buffID].xi = pn[buffID].xi + dt*pn[buffID].v;

				// pnp1[pi].v = (pn[ii].v /* + pnp1[ii].vPert*/) + dt * (a * pn[ii].f + b * res[ii]);
				// pnp1[ii].f = res[ii];
				// pnp1[ii].Af = Af[ii];
				// pnp1[buffID].Rrho = Rrho[buffID];

				// pnp1[buffID].rho = pn[buffID].rho+dt*(a*pn[buffID].Rrho+b*Rrho[buffID]);
				// pnp1[buffID].p = B*(pow(pnp1[buffID].rho/fvar.rho0,gam)-1);
				
			}	
		}

	}/*End pragma omp parallel*/

}


void Newmark_Beta(KDTREE& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, 
	real const& a, real const& b, real const& c, real const& d, real const& B, real const& gam,
	MESH& cells, vector<size_t>& cellsused, outl& outlist, DELTAP const& dp,
	real& logbase, uint& k, real& error1, real& error2, vector<StateVecD>& xih,
	State& pn, State& pnp1, StateVecD& Force, StateVecD& dropVel)
{
	uint nUnstab = 0;

	while (error1 > -5.0)
	{
		/*Previous State for error calc*/
		#pragma omp parallel for shared(pnp1)
		for (size_t  ii=start; ii < end; ++ii)
			xih[ii-start] = pnp1[ii].xi;	
		
		Do_NB_Iter(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,
			cells,cellsused,outlist,dp,pn,pnp1,Force,dropVel);

		int errstate = Check_Error(TREE,svar,fvar,start,end,error1,error2,logbase,
							cellsused,outlist,xih,pn,pnp1,k,nUnstab);

		if (errstate == 0)	/*Continue as normal*/
			k++;
		else if (errstate == 1)	/*Sub iterations exceeded*/
			break;

		error2 = error1;

		// cout << "It: " << k << " Error: " << error1 << endl;
		
	} /*End of subits*/
}

#endif