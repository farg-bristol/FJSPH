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
		State& pn, State& pnp1, uint& k)
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

	if (error1-error2 > 0.0 /*|| std::isnan(error1)*/)
	{	/*If simulation starts diverging, then reduce the timestep and try again.*/
		
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
	}	/*Check if we've exceeded the maximum iteration count*/
	else if (k > svar.subits)
	{
		return 1;
	}

	/*Otherwise, roll forwards*/
	return 0;
	
}

void Get_Resid(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, 
	real const& a, real const& b, real const& c, real const& d, real const& B, real const& gam,
	MESH& cells, vector<size_t> const& cellsused,
	vector<vector<Part>> const& neighb, outl const& outlist, DELTAP const& dp,
	State& pn, State& pnp1, State& airP, StateVecD& Force, StateVecD& dropVel)
{
	/*Update the state at time n+1*/
	#pragma omp parallel shared(pn,svar,fvar) /*reduction(+:Force,dropVel)*/
	{
		// cout << "K: " << k << endl;
		// cout << "Calculating forces" << endl;
		vector<StateVecD> res(end,StateVecD::Zero());
		vector<StateVecD> Af(end,StateVecD::Zero());
		vector<real> Rrho(end,0.0);
		vector<real> curve(end,0.0);
		vector<StateVecD> norm(end, StateVecD::Zero());
		// vector<real> curve(end,0.0);

		/*Calculate Di term for delta-SPH*/
		// vector<real> Di(end,0.0);
		// vector<real> isSurf(end,0.0);
		// dSPH_PreStep(fvar,start,end,pnp1,outlist,Di,isSurf);

		// cout << "Neighbour list size: " << neighb.size() << "  Outlist size: " << outlist.size() << endl;
		// cout << "Start: " << start << "  End: " << end << endl;
		Force = StateVecD::Zero();
		svar.AForce = StateVecD::Zero();

		Forces(svar,fvar,avar,cells,pnp1,neighb,outlist,dp,res,Rrho,Af,Force,curve); /*Guess force at time n+1*/

		const real dt = svar.dt;
		const real dt2 = dt*dt;

		if(svar.Asource == 2)
		{
			#pragma omp for schedule(static) nowait
			for(size_t const& ii : cellsused)
			{
					cells.fNum[ii] = 0;
					cells.fMass[ii] = 0.0;
					cells.vFn[ii] = StateVecD::Zero();
					cells.vFnp1[ii] = StateVecD::Zero();
			}
		}

		#pragma omp for schedule(static) nowait
		for (size_t ii=0; ii < start; ++ii)
		{	/****** BOUNDARY PARTICLES ***********/
			pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+Rrho[ii]);
			pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
			// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
			// pnp1[ii].p = fvar.pPress;
			// pnp1[ii].rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
			pnp1[ii].f = res[ii];
			pnp1[ii].Rrho = Rrho[ii];
		}


		#pragma omp for schedule(static) nowait 
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* START = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			if(pnp1[ii].b == PartState.START_ || pnp1[ii].b == PartState.BACK_)
			{   
				StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
				if(vec(1) > svar.clear)
				{	/*Tag it as clear if it's higher than the plane of the exit*/
					pnp1[ii].b=PartState.PIPE_;
				}
				else
				{	/*For the particles marked 1, perform a prescribed motion*/
					pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
					// pnp1[ii].f = res[ii];
					// pnp1[ii].Rrho = Rrho[ii];
					// pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
					// pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
				}
			}				

			if(pnp1[ii].b > PartState.START_)
			{	
				StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);

				if(pnp1[ii].b == PartState.PIPE_)
				{	/*Do a check to see if it needs to be given an aero force*/	
					if(vec(1) > 0.0)
					{
						pnp1[ii].b = PartState.FREE_;
						if(svar.Asource == 1 || svar.Asource == 2)
						{
							/*retrieve the cell it's in*/
							FirstCell(svar,end,ii,TREE.CELL, cells, pnp1, pn);
						}
					}	
				}


				/*For any other particles, intergrate as normal*/
				pnp1[ii].xi = pn[ii].xi+dt*(pn[ii].v + pnp1[ii].vPert)+dt2*(c*pn[ii].f+d*res[ii]);
				pnp1[ii].v =  (pn[ii].v + pnp1[ii].vPert) +dt*(a*pn[ii].f+b*res[ii]);
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


				if(pnp1[ii].xi(0) != pnp1[ii].xi(0))
				{
					cout << "Position is nan:  " << ii << endl;
				}

				if(pnp1[ii].v != pnp1[ii].v)
				{
					cout << "velocity is nan:  " << ii << endl;
				}

				if(pnp1[ii].f != pnp1[ii].f)
				{
					cout << "Force is nan:  " << ii << endl;
				}

// #if SIMDIM == 2
// 					if(real(outlist[ii].size()) < 5.26 * wDiff[ii] +30)
// 					{
// 						pnp1[ii].normal = norm[ii];
// 						pnp1[ii].surf = curve[ii];
// 						pnp1[ii].theta = 1.0;
// 					}
// 					else
// 					{
// 						pnp1[ii].normal = StateVecD::Zero();
// 						pnp1[ii].surf = 0.0;
// 						pnp1[ii].theta = 0.0;
// 					}
// #else
// 					if(real(outlist[ii].size()) < 297.25 * wDiff[ii] + 34)
// 					{
// 						pnp1[ii].theta = wDiff[ii];
// 						pnp1[ii].normal = norm[ii];
// 						pnp1[ii].surf = curve[ii];
// 					}
// 					else
// 					{
// 						pnp1[ii].normal = StateVecD::Zero();
// 						pnp1[ii].surf = 0.0;
// 						pnp1[ii].theta = 0.0;
// 					}
// #endif

				
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

			pnp1[ii].curve = curve[ii];
			pnp1[ii].s = outlist[ii].size();

			// pnp1[ii].theta = wDiff[ii];
			pnp1[ii].nNeigb = real(outlist[ii].size());
			pnp1[ii].bNorm = norm[ii];
			
		} /*End fluid particles*/

	}/*End pragma omp parallel*/

}


void Newmark_Beta(KDTREE& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, 
	real const& a, real const& b, real const& c, real const& d, real const& B, real const& gam,
	MESH& cells, vector<size_t>& cellsused,
	vector<vector<Part>> const& neighb, outl& outlist, DELTAP const& dp,
	real& logbase, uint& k, real& error1, real& error2, 
	State& pn, State& pnp1, State& airP, StateVecD& Force, StateVecD& dropVel)
{
	vector<StateVecD> xih(end-start);

	while (error1 > -7.0)
	{
		/*Previous State for error calc*/
		#pragma omp parallel for shared(pnp1)
		for (size_t  ii=start; ii < end; ++ii)
			xih[ii-start] = pnp1[ii].xi;	
		
		Get_Resid(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,
			cells,cellsused,neighb,outlist,dp,pn,pnp1,airP,Force,dropVel);

		int errstate = Check_Error(TREE,svar,fvar,start,end,error1,error2,logbase,
							cellsused,outlist,xih,pn,pnp1,k);

		if (errstate == 0)	/*Continue as normal*/
			k++;
		else if (errstate == 1)	/*Sub iterations exceeded*/
			break;

		error2 = error1;

		// Particle_Shift(svar,fvar,start,end,outlist,dp,pnp1);
		// cout << error1 << endl;
		
	} /*End of subits*/
}

#endif