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

	if (error1-error2 > 0.0 /*|| std::isnan(error1)*/)
	{	/*If simulation starts diverging, then reduce the timestep and try again.*/
		nUnstab++;

		if(nUnstab > 4)
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

		return 0;
	}	/*Check if we've exceeded the maximum iteration count*/
	else if (k > svar.subits)
	{
		return 1;
	}

	/*Otherwise, roll forwards*/
	nUnstab = 0;
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
	
	// cout << "Calculating forces" << endl;
	vector<StateVecD> res;
	vector<StateVecD> Af;
	vector<real> Rrho;
	vector<real> curve;

	Force = StateVecD::Zero();
	svar.AForce = StateVecD::Zero();

	Forces(svar,fvar,avar,cells,pnp1,neighb,outlist,dp,res,Rrho,Af,Force,curve); /*Guess force at time n+1*/

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

		#pragma omp for schedule(static) nowait
		for (size_t ii=0; ii < start; ++ii)
		{	/****** BOUNDARY PARTICLES ***********/
			pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+Rrho[ii]);
			pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
			// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
			// pnp1[ii].p = fvar.pPress;
			// pnp1[ii].rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
			// pnp1[ii].f = res[ii];
			pnp1[ii].Rrho = Rrho[ii];
		}


		#pragma omp for nowait 
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* START = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			if (pnp1[ii].b == PartState.BACK_)
			{
// 				StateVecD xin = svar.Transp*pn[ii].xi;
// 				StateVecD veln  = svar.Transp*(pn[ii].v /*+ pnp1[ii].vPert*/);
// 				StateVecD acc = svar.Transp*res[ii];
// 				StateVecD accn = svar.Transp*pn[ii].f;

// 				/*Integrate the pipe x-direction as normal, and allow particles to find their spacing.*/
//  				StateVecD xi = StateVecD::Zero();
//  				xi(0) = xin(0) + 0.25*(dt2*(c*accn(0)+d*acc(0)));
//  				xi(1) = xin(1) + dt*veln(1);
 				
// 				StateVecD v  = StateVecD::Zero();
// 				v(0) = veln(0) + 0.25*(dt*(a*accn(0)+b*acc(0)));
// 				v(1) = veln(1);

// #if SIMDIM == 3
// 				xi(2) = xin(2) + 0.25*(dt2*(c*accn(2)+d*acc(2)));
// 				v(2) = veln(2) + 0.25*(dt*(a*accn(2)+b*acc(2)));
// #endif

// 				xi = svar.Rotate*xi;
// 				v = svar.Rotate*v;
				
// 				pnp1[ii].xi = xi;
// 				pnp1[ii].v = v;
// #if SIMDIM == 3
// 				pnp1[ii].f = svar.Rotate*StateVecD(acc(0),0.0,acc(2));
// #else
// 				pnp1[ii].f = svar.Rotate*StateVecD(acc(0),0.0);
// #endif

				pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
				// pnp1[ii].f = res[ii];
				pnp1[ii].Rrho = Rrho[ii];
				pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
				pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
			}
			else if(pnp1[ii].b == PartState.START_ )
			{   
				StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
				if(vec(1) > svar.clear)
				{	/*Tag it as clear if it's higher than the plane of the exit*/
					pnp1[ii].b=PartState.PIPE_;
				}
				else
				{	/*For the particles marked 1, perform a prescribed motion*/

					/*Try setting a force to maintain a minimum jet-wise velocity*/

					/*Need to transform into the vertical domain...*/
// 					StateVecD xin = svar.Transp*pn[ii].xi;
// 				    StateVecD veln  = svar.Transp*(pn[ii].v /*+ pnp1[ii].vPert*/);
// 					// StateVecD vel = svar.Transp*pnp1[ii].v;
// 					StateVecD acc = svar.Transp*res[ii];
// 					StateVecD accn = svar.Transp*pn[ii].f;

// 					// if(vel(1) < (getVelocity(avar.vJet,vec(0)*vec(0),0.5*svar.Jet(0) + 2*svar.dx))(1))
// 					// {
// 					// 	acc(1) += 1e2 * pow(vel(1)-avar.vJetMag,2);
// 					// }

// 					// acc = svar.Rotate * acc;


// 					// pnp1[ii].xi = pn[ii].xi+ dt*(pn[ii].v + pnp1[ii].vPert)+dt2*(c*pn[ii].f+d*acc);
// 					// pnp1[ii].v =  (pn[ii].v /*+ pnp1[ii].vPert*/) +dt*(a*pn[ii].f+b*acc);
// 					// pnp1[ii].f = svar.Rotate*acc;

// 					// StateVecD xi = xin + dt*(veln + svar.Transp*pnp1[ii].vPert)+dt2*(c*accn+d*accnp1);
// 	 				// 	pnp1[ii].xi = svar.Rotate*(xi+svar.Start);
// 	 				// 	pnp1[ii].v = svar.Rotate*v;

// 					/*Integrate the pipe x-direction as normal, and allow particles to find their spacing.*/
// 	 				StateVecD xi = StateVecD::Zero();
// 	 				xi(0) = xin(0) + 0.25*(dt2*(c*accn(0)+d*acc(0)));
// 	 				xi(1) = xin(1) + dt*veln(1);

	 				
// 					StateVecD v  = StateVecD::Zero();
// 					v(0) = veln(0) + 0.25*(dt*(a*accn(0)+b*acc(0)));
//  					v(1) = veln(1);

// #if SIMDIM == 3
// 					xi(2) = xin(2) + 0.25*(dt2*(c*accn(2)+d*acc(2)));
// 					v(2) = veln(2) + 0.25*(dt*(a*accn(2)+b*acc(2)));
// #endif

// 					xi = svar.Rotate*xi;
//  					v = svar.Rotate*v;
					
// 					pnp1[ii].xi = xi;
// 					pnp1[ii].v = v;
// #if SIMDIM == 3
// 					pnp1[ii].f = svar.Rotate*StateVecD(acc(0),0.0,acc(2));
// #else
// 					pnp1[ii].f = svar.Rotate*StateVecD(acc(0),0.0);
// #endif
					// cout << acc(0) << "  " << acc(1) << endl;
					// pnp1[ii].xi = pn[ii].xi+dt*(pn[ii].v + pnp1[ii].vPert)+dt2*(c*pn[ii].f+d*res[ii]);
					pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
					// pnp1[ii].f = res[ii];
					pnp1[ii].Rrho = Rrho[ii];
					pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
					pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
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
				pnp1[ii].v =  (pn[ii].v/* + pnp1[ii].vPert*/) +dt*(a*pn[ii].f+b*res[ii]);
				pnp1[ii].f = res[ii];
				pnp1[ii].Af = Af[ii];
				pnp1[ii].Rrho = Rrho[ii];

				// if(res[ii].norm() > 10)
				// {
				// 	#pragma omp critical
				// 	cout << res[ii](0) << "  " << res[ii](1) << "  " << Af[ii](0) << "  " << Af[ii](1) <<  endl;
				// }


				
				pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*Rrho[ii]);
				pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
				// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
				
				Force += res[ii]*pnp1[ii].m;
				dropVel += (pnp1[ii].v + pnp1[ii].vPert);
				
				// pnp1[ii].theta = dp.lam[ii];
				// pnp1[ii].nNeigb = real(outlist[ii].size());
				// pnp1[ii].curve = curve[ii];
				pnp1[ii].pDist = vec.norm();


				// if(pnp1[ii].xi(0) != pnp1[ii].xi(0))
				// {
				// 	cout << "Position is nan:  " << ii << endl;
				// }

				// if(pnp1[ii].v != pnp1[ii].v)
				// {
				// 	cout << "velocity is nan:  " << ii << endl;
				// }

				// if(pnp1[ii].f != pnp1[ii].f)
				// {
				// 	cout << "Force is nan:  " << ii << endl;
				// }

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
			// pnp1[ii].bNorm = norm[ii];
			
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

	uint nUnstab = 0;

	while (error1 > -7.0)
	{
		/*Previous State for error calc*/
		#pragma omp parallel for shared(pnp1)
		for (size_t  ii=start; ii < end; ++ii)
			xih[ii-start] = pnp1[ii].xi;	
		
		Get_Resid(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,
			cells,cellsused,neighb,outlist,dp,pn,pnp1,airP,Force,dropVel);

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