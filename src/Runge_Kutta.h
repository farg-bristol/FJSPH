/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <chrono>
#include <unordered_set>
#include "Var.h"
#include "Resid.h"
#include "Crossing.h"

int Check_RK_Error(SIM const& svar, size_t const& start, size_t const& end, 
		real& error1, real& error2, real& logbase, State const& pn, State const& pnp1, uint& k)
{
	/****** FIND ERROR ***********/
	real errsum = 0.0;
	#pragma omp parallel for reduction(+:errsum) schedule(static)
	for (size_t ii = start; ii < end; ++ii)
	{
		StateVecD r = pnp1[ii].xi-pn[ii].xi;
		errsum += r.squaredNorm();
	}

	if(k == 1)
		logbase=log10(sqrt(errsum/(real(svar.totPts))));

	error1 = log10(sqrt(errsum/(real(svar.totPts)))) - logbase;
	// cout << RestartCount << "  " << k << "  " << error1  << "  " << svar.dt << endl;
	// cout << k << "  " << error1 << "  " << svar.dt << endl;

	cout << "Step: " << k << "  Error: " << error1 << endl;

	if (error1-error2 > 0.0 /*|| std::isnan(error1)*/)
	{	/*If simulation starts diverging, then reduce the timestep and try again.*/
		return -1;
	}

	/*Otherwise, roll forwards*/
	return 0;
}


void Step_1(real const& dt, 
		StateVecD const& xin, StateVecD const& veln, real const& rhon, /*State at time n*/
		StateVecD const& res_1, real const& Rrho_1, 	/*residuals at step 1*/
		StateVecD& xi_2, StateVecD& vel_2, real& rho_2) /*State at step 1*/
{
	xi_2 = xin + 0.5 * dt * veln;
	vel_2 = veln + 0.5 * dt * res_1;
	rho_2 = rhon + 0.5 * dt * Rrho_1;
}

void Step_2(real const& dt, 
		StateVecD const& xin, StateVecD const& veln, real const& rhon,	/*State at time n*/
		StateVecD const& vel_2, StateVecD const& res_2, real const& Rrho_2, /*residuals at step 2*/
		StateVecD& xi_3, StateVecD& vel_3, real& rho_3) /*State at step 2*/
{
	xi_3 = xin + 0.5 * dt * vel_2;
	vel_3 = veln + 0.5 * dt * res_2;
	rho_3 = rhon + 0.5 * dt * Rrho_2;
}

void Step_3(real const& dt, 
		StateVecD const& xin, StateVecD const& veln, real const& rhon,	/*State at time n*/
		StateVecD const& vel_3, StateVecD const& res_3, real const& Rrho_3, /*residuals at step 3*/
		StateVecD& xi_4, StateVecD& vel_4, real& rho_4) /*State at step 3*/
{
	xi_4 = xin + dt * vel_3;
	vel_4 = veln + dt * res_3;
	rho_4 = rhon + dt * Rrho_3;
}

void Step_4(real const& dt, 
		StateVecD const& xi_1, StateVecD const& vel_1, real const& rho_1, 	/*State at time n*/
		 StateVecD const& res_1, real const& Rrho_1,
		StateVecD const& vel_2, StateVecD const& res_2, real const& Rrho_2, /*residuals at step 2*/
		StateVecD const& vel_3, StateVecD const& res_3, real const& Rrho_3, /*residuals at step 3*/
		StateVecD const& vel_4, StateVecD const& res_4, real const& Rrho_4, /*residuals at step 4*/
		StateVecD& xi_np1, StateVecD& vel_np1, real& rho_np1)	/*State at time n+1*/
{
	xi_np1 = xi_1 +  (dt/6.0) * (vel_1 + 2.0 * vel_2 + 2.0 * vel_3 + vel_4);
	vel_np1 = vel_1 + (dt/6.0) * (res_1 + 2.0 * res_2 + 2.0 * res_3 + res_4);
	rho_np1 = rho_1 + (dt/6.0) * (Rrho_1 + 2.0 * Rrho_2 + 2.0 * Rrho_3 + Rrho_4);
}

void Get_First_RK(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, 
	MESH& cells, vector<size_t> const& cellsused,
	vector<vector<Part>> const& neighb, outl const& outlist, DELTAP const& dp, real& logbase,
	State& pn, State& st_2, real& error1)
{
	const real dt = svar.dt;
	
	// real error1 = 0.0;
	real error2 = 0.0;
	// logbase = 0.0;
	uint k = 1;

	#pragma omp parallel shared(svar,pn,st_2) /*reduction(+:Force,dropVel)*/
	{
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

		uint w = 15;

		/********************************************************************/
		/*************************  STEP 1  *********************************/
		/********************************************************************/
		#pragma omp for schedule(static) nowait
		for (size_t ii=0; ii < start; ++ii)
		{	/****** BOUNDARY PARTICLES ***********/
			st_2[ii].rho = pn[ii].rho + 0.5 * dt * pn[ii].Rrho;
			// st_2[ii].p = B*(pow(st_2[ii].rho/fvar.rho0,gam)-1);
			st_2[ii].p = fvar.Cs*fvar.Cs * (st_2[ii].rho - fvar.rho0);
		}

		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* START = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			if(st_2[ii].b == PartState.START_ || st_2[ii].b == PartState.BACK_)
			{   
				st_2[ii].xi = pn[ii].xi + 0.5*dt*pn[ii].v;
			}

			if(st_2[ii].b > PartState.START_)
			{	

				// Step_1(dt,pn[ii].xi,pn[ii].v,pn[ii].rho,pn[ii].f,pn[ii].Rrho,
				// 		st_2[ii].xi,st_2[ii].v,st_2[ii].rho);

				st_2[ii].xi = pn[ii].xi + 0.5 * dt * pn[ii].v;
				st_2[ii].v = pn[ii].v + 0.5 * dt * pn[ii].f;
				st_2[ii].rho = pn[ii].rho + 0.5 * dt * pn[ii].Rrho;
				st_2[ii].p = fvar.Cs*fvar.Cs * (st_2[ii].rho - fvar.rho0);

				#pragma omp critical
				cout << setw(6) << ii << setw(w) << st_2[ii].xi(0) << setw(w) << st_2[ii].xi(1) << 
				setw(w) << st_2[ii].v(0) << setw(w) << st_2[ii].v(1) << setw(w) << 
				setw(w) << pn[ii].f(0) << setw(w) << pn[ii].f(1) << setw(w) <<
				setw(w) << st_2[ii].rho << setw(w) << st_2[ii].m << endl;


				if(svar.Asource == 2 && st_2[ii].b == PartState.FREE_)
				{
					#pragma omp atomic
						cells.fNum[st_2[ii].cellID]++;
					#pragma omp atomic
						cells.fMass[st_2[ii].cellID] += st_2[ii].m;

					#pragma omp critical
					{
						cells.vFn[st_2[ii].cellID] += pn[ii].v;
						cells.vFnp1[st_2[ii].cellID] +=st_2[ii].v;
					}	
				}
			}
		}

	}

	void(Check_RK_Error(svar,start,end,error1,error2,logbase,pn,st_2,k));

}

int Perform_RK4(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, 
	MESH& cells, vector<size_t> const& cellsused,
	vector<vector<Part>> const& neighb, outl const& outlist, DELTAP const& dp, real& logbase,
	State& pn, State& st_2, State& pnp1, State& airP, StateVecD& Force, StateVecD& dropVel, real& error1)
{
	/*Create the vectors*/		
	vector<StateVecD> res_2(end,StateVecD::Zero());
	vector<real> Rrho_2(end,0.0);

	State st_3 = pn; /*Step 3*/
	vector<StateVecD> res_3(end,StateVecD::Zero());
	vector<real> Rrho_3(end,0.0);

	State st_4 = pn; /*Step 4*/
	vector<StateVecD> res_4(end,StateVecD::Zero());
	vector<real> Rrho_4(end,0.0);

	vector<StateVecD> Af(end,StateVecD::Zero());
	vector<real> wDiff(end,0.0);
	vector<StateVecD> norm(end, StateVecD::Zero());

	vector<real> curve;
	Force = StateVecD::Zero();
	svar.AForce = StateVecD::Zero();

	const real dt = svar.dt;
	
	// real error1 = 0.0;
	// real error2 = 0.0;
	// logbase = 0.0;
	uint k = 2;

	// #pragma omp parallel shared(svar,pn,st_2,res_2,Rrho_2) /*reduction(+:Force,dropVel)*/
	// {
	// 	if(svar.Asource == 2)
	// 	{
	// 		#pragma omp for schedule(static) nowait
	// 		for(size_t const& ii : cellsused)
	// 		{
	// 				cells.fNum[ii] = 0;
	// 				cells.fMass[ii] = 0.0;
	// 				cells.vFn[ii] = StateVecD::Zero();
	// 				cells.vFnp1[ii] = StateVecD::Zero();
	// 		}
	// 	}


	// 	/********************************************************************/
	// 	/*************************  STEP 1  *********************************/
	// 	/********************************************************************/
	// 	#pragma omp for schedule(static) nowait
	// 	for (size_t ii=0; ii < start; ++ii)
	// 	{	/****** BOUNDARY PARTICLES ***********/
	// 		st_2[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho);
	// 		// st_2[ii].p = B*(pow(st_2[ii].rho/fvar.rho0,gam)-1);
	// 		st_2[ii].p = fvar.Cs*fvar.Cs * (st_2[ii].rho - fvar.rho0);
	// 	}

	// 	#pragma omp for schedule(static) nowait
	// 	for (size_t ii=start; ii < end; ++ii)
	// 	{	/****** FLUID PARTICLES **************/
			
	// 		/* START = pipe particle receiving prescribed motion            */
	// 		/* BACK = the latest particle in that column                    */
	// 		/* PIPE = in the pipe, with free motion                         */
	// 		/* FREE = free of the pipe and receives an aero force           */

	// 		/*Check if the particle is clear of the starting area*/
	// 		if(pnp1[ii].b == PartState.START_ || pnp1[ii].b == PartState.BACK_)
	// 		{   
	// 			StateVecD vec = svar.Transp*(pn[ii].xi-svar.Start);
	// 			if(vec(1) > svar.clear)
	// 			{	/*Tag it as clear if it's higher than the plane of the exit*/
	// 				pnp1[ii].b=PartState.PIPE_;
	// 			}
	// 			else
	// 			{	/*For the particles marked 1, perform a prescribed motion*/
	// 				st_2[ii].xi = pn[ii].xi + dt*pn[ii].v;
	// 			}
	// 		}

	// 		if(pnp1[ii].b > PartState.START_)
	// 		{	
	// 			if(pnp1[ii].b == PartState.PIPE_)
	// 			{	/*Do a check to see if it needs to be given an aero force*/
	// 				StateVecD vec = svar.Transp*(pn[ii].xi-svar.Start);
	// 				if(vec(1) > 0.0)
	// 				{
	// 					pnp1[ii].b = PartState.FREE_;
	// 					if(svar.Asource == 1 || svar.Asource == 2)
	// 					{
	// 						/*retrieve the cell it's in*/
	// 						FirstCell(svar,end,ii,TREE.CELL, cells, pnp1, pn);
	// 					}
	// 				}	
	// 			}

	// 			Step_1(dt,pn[ii].xi,pn[ii].v,pn[ii].rho,pn[ii].f,pn[ii].Rrho,
	// 					st_2[ii].xi,st_2[ii].v,st_2[ii].rho);


	// 			if(svar.Asource == 2 && pnp1[ii].b == PartState.FREE_)
	// 			{
	// 				#pragma omp atomic
	// 					cells.fNum[pnp1[ii].cellID]++;
	// 				#pragma omp atomic
	// 					cells.fMass[pnp1[ii].cellID] += pnp1[ii].m;

	// 				#pragma omp critical
	// 				{
	// 					cells.vFn[pnp1[ii].cellID] += pn[ii].v;
	// 					cells.vFnp1[pnp1[ii].cellID] +=st_2[ii].v;
	// 				}	
	// 			}
	// 		}
	// 	}

	// }

	// int errstate = Check_RK_Error(svar,start,end,error1,error2,logbase,pn,st_2,k);

	// if (errstate < 0)
	// 	return -1;

	// k++;

	// uint w = 15;

	/********************************************************************/
	/*************************  STEP 2  *********************************/
	/********************************************************************/
	#pragma omp parallel shared(svar,pn,st_2,res_2,Rrho_2,st_3,res_3,Rrho_3)
	{
		Forces(svar, fvar, avar, cells, st_2, neighb, outlist, dp, res_2, Rrho_2, Af, Force, curve);
		#pragma omp for schedule(static) nowait
		for (size_t ii=0; ii < start; ++ii)
		{	/****** BOUNDARY PARTICLES ***********/
			st_3[ii].rho = pn[ii].rho+0.5*dt*(Rrho_2[ii]);
			// st_3[ii].p = B*(pow(st_3[ii].rho/fvar.rho0,gam)-1);
			st_3[ii].p = fvar.Cs*fvar.Cs * (st_3[ii].rho - fvar.rho0);
		}

		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* START = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			if(st_2[ii].b == PartState.START_ || st_2[ii].b == PartState.BACK_)
			{   
				/*For the particles marked 1, perform a prescribed motion*/
				st_3[ii].xi = pn[ii].xi + 0.5 * dt * pn[ii].v;
				
			}

			if(st_2[ii].b > PartState.START_)
			{	
				// Step_2(dt,pn[ii].xi,pn[ii].v,pn[ii].rho,
				// 		st_2[ii].v,res_2[ii],Rrho_2[ii],
				// 		st_3[ii].xi,st_3[ii].v,st_3[ii].rho);

				st_3[ii].xi = pn[ii].xi + 0.5 * dt * st_2[ii].v;
				st_3[ii].v = pn[ii].v + 0.5 * dt * res_2[ii];
				st_3[ii].rho = pn[ii].rho + 0.5 * dt * Rrho_2[ii];
				st_3[ii].p = fvar.Cs*fvar.Cs * (st_3[ii].rho - fvar.rho0);
			}

			// #pragma omp critical
			// cout << setw(6) << ii << setw(w) << st_3[ii].xi(0) << setw(w) << st_3[ii].xi(1) << 
			// setw(w) << st_3[ii].v(0) << setw(w) << st_3[ii].v(1) << setw(w) << 
			// res_2[ii](0) << setw(w) << res_2[ii](1) << setw(w) << Rrho_2[ii] <<
			// setw(w) << st_2[ii].rho << setw(w) << st_2[ii].m << endl;
		}

	}

	// if (Check_RK_Error(svar,start,end,error1,error2,logbase,st_2,st_3,k) < 0)
	// 	return -1;

	k++;

	/********************************************************************/
	/*************************  STEP 3  *********************************/
	/********************************************************************/
	#pragma omp parallel shared(svar,pn,st_3,res_3,Rrho_3,st_4,res_4,Rrho_4)
	{
		Forces(svar, fvar, avar, cells, st_3, neighb, outlist, dp, res_3, Rrho_3, Af, Force, curve);
		#pragma omp for schedule(static) nowait
		for (size_t ii=0; ii < start; ++ii)
		{	/****** BOUNDARY PARTICLES ***********/
			st_4[ii].rho = pn[ii].rho+0.5*dt*(Rrho_3[ii]);
			// st_3[ii].p = B*(pow(st_3[ii].rho/fvar.rho0,gam)-1);
			st_4[ii].p = fvar.Cs*fvar.Cs * (st_4[ii].rho - fvar.rho0);
		}

		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* START = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			if(st_3[ii].b == PartState.START_ || st_3[ii].b == PartState.BACK_)
			{   
				/*For the particles marked 1, perform a prescribed motion*/
				st_4[ii].xi = pn[ii].xi + dt*pn[ii].v;
			}

			if(st_3[ii].b > PartState.START_)
			{	
				// Step_3(dt,pn[ii].xi,pn[ii].v,pn[ii].rho,
				// 		st_3[ii].v,res_3[ii],Rrho_3[ii],
				// 		st_4[ii].xi,st_4[ii].v,st_4[ii].rho);

				st_4[ii].xi = pn[ii].xi + dt * st_3[ii].v;
				st_4[ii].v = pn[ii].v + dt * res_3[ii];
				st_4[ii].rho = pn[ii].rho + dt * Rrho_3[ii];
				st_4[ii].p = fvar.Cs*fvar.Cs * (st_4[ii].rho - fvar.rho0);
			}
		}
	}

	// if (Check_RK_Error(svar,start,end,error1,error2,logbase,st_3,st_4,k) < 0)
	// 	return -1;

	k++;
	/********************************************************************/
	/*************************  STEP 4  *********************************/
	/********************************************************************/

	#pragma omp parallel shared(svar,pn,pnp1,st_2,res_2,Rrho_2,st_3,res_3,Rrho_3,st_4,res_4,Rrho_4)
	{
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

		Forces(svar, fvar, avar, cells, st_4, neighb, outlist, dp, res_4, Rrho_4, Af, Force, curve);
		#pragma omp for schedule(static) nowait
		for (size_t ii=0; ii < start; ++ii)
		{	/****** BOUNDARY PARTICLES ***********/
			pnp1[ii].rho = pn[ii].rho + 
				(dt/6.0) * (pn[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii]);
			// st_3[ii].p = B*(pow(st_3[ii].rho/fvar.rho0,gam)-1);
			pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);

			pnp1[ii].f = res_4[ii];
			pnp1[ii].Rrho = Rrho_4[ii];
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
				}
			}

			if(pnp1[ii].b > PartState.START_)
			{	
				if(pnp1[ii].b == PartState.PIPE_)
				{	/*Do a check to see if it needs to be given an aero force*/
					StateVecD vec = svar.Transp*(st_4[ii].xi-svar.Start);
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

				// Step_4(dt,pn[ii].xi,pn[ii].v,pn[ii].rho,
				// 		pn[ii].f,pn[ii].Rrho,
				// 		st_2[ii].v,res_2[ii],Rrho_2[ii],
				// 		st_3[ii].v,res_3[ii],Rrho_3[ii],
				// 		st_4[ii].v,res_4[ii],Rrho_4[ii],
				// 		pnp1[ii].xi,pnp1[ii].v,pnp1[ii].rho);

				pnp1[ii].xi = pn[ii].xi +  
					(dt/6.0) * (pn[ii].v + 2.0 * st_2[ii].v + 2.0 * st_3[ii].v + st_4[ii].v);
				pnp1[ii].v = pn[ii].v + 
					(dt/6.0) * (pn[ii].f + 2.0 * res_2[ii] + 2.0 * res_3[ii] + res_4[ii]);
				pnp1[ii].rho = pn[ii].rho + 
					(dt/6.0) * (pn[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii]);
				pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
				pnp1[ii].f = res_4[ii];
				pnp1[ii].Rrho = Rrho_4[ii];

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

			pnp1[ii].theta = wDiff[ii];
			pnp1[ii].s = outlist[ii].size();

			// pnp1[ii].theta = wDiff[ii];
			pnp1[ii].nNeigb = real(outlist[ii].size());
			pnp1[ii].bNorm = norm[ii];
		}

	}/*End pragma omp parallel*/

	// if (Check_RK_Error(svar,start,end,error1,error2,logbase,st_4,pnp1,k) < 0)
	// 	return -1;

	return 0;
}

real Runge_Kutta4(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, 
	MESH& cells, vector<size_t> const& cellsused,
	vector<vector<Part>> const& neighb, outl const& outlist, DELTAP const& dp, real& logbase,
	State& pn, State& st_2, State& pnp1, State& airP, StateVecD& Force, StateVecD& dropVel)
{
	real error = 0.0;
	int errstate = 1;
	while(errstate != 0)
	{
		errstate = Perform_RK4(TREE,svar,fvar,avar,start,end,cells,cellsused,
								neighb,outlist,dp,logbase,pn,st_2,pnp1,airP,Force,dropVel,error);

		if(errstate < 0)
		{
			/*Timestep was unstable - reduce the timestep*/
			svar.dt *= 0.5;
		}
	}

	return error;
}

#endif