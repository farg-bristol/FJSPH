/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <chrono>
#include <unordered_set>
#include "Var.h"
#include "Resid.h"
#include "Containment.h"

real Check_RK_Error(SIM const& svar, size_t const& start, size_t const& end, 
		real& error1, real& error2, real& logbase, SPHState const& pn, SPHState const& pnp1, uint& k)
{
	/****** FIND ERROR ***********/
	real errsum = 0.0;
	#pragma omp parallel for reduction(+:errsum)
	for (size_t ii = start; ii < end; ++ii)
	{
		errsum += (pnp1[ii].xi-pn[ii].xi).squaredNorm();
	}

	if(k == 1)
		logbase=log10(sqrt(errsum/(real(svar.totPts))));

	error1 = log10(sqrt(errsum/(real(svar.totPts)))) - logbase;

	return error1;
}

/* <summary> Peform the first stage of the Runge-Kutta integration
	to get the first guess of time n+1 (regarded here as time n+1/4)
	to perform neighbour search and dissipation terms before freezing </summary */
real Get_First_RK(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, real const& B, real const& gam, MESH& cells, vector<size_t> const& cellsused,
	 OUTL const& outlist, DELTAP const& dp, real& logbase, SPHState& pn, SPHState& st_2, real& error1)
{
	const real dt = svar.dt;
	
	// real error1 = 0.0;
	real error2 = 0.0;
	// logbase = 0.0;
	uint k = 1;

	vector<StateVecD> res_1(end, StateVecD::Zero());
	vector<real> Rrho_1(end, 0.0);

	StateVecD Force = StateVecD::Zero();
	vector<StateVecD> Af;
	svar.AForce = StateVecD::Zero();

	Get_Boundary_Pressure(fvar,start,outlist,st_2);
	Forces(svar, fvar, avar, cells, st_2, outlist, dp, res_1, Rrho_1, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,st_2) /*reduction(+:Force,dropVel)*/
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

		// uint w = 15;

		/********************************************************************/
		/*************************  STEP 1  *********************************/
		/********************************************************************/
		// #pragma omp for schedule(static) nowait
		// for (size_t ii=0; ii < start; ++ii)
		// {	/****** BOUNDARY PARTICLES ***********/
		// 	st_2[ii].rho = pn[ii].rho + 0.5 * dt * Rrho_1[ii];
		// 	st_2[ii].p = B*(pow(st_2[ii].rho/fvar.rho0,gam)-1);
		// 	// st_2[ii].p = fvar.Cs*fvar.Cs * (st_2[ii].rho - fvar.rho0);
		// }

		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* BUFFER = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			if(st_2[ii].b > BUFFER)
			{	

				// Step_1(dt,pn[ii].xi,pn[ii].v,pn[ii].rho,pn[ii].acc,pn[ii].Rrho,
				// 		st_2[ii].xi,st_2[ii].v,st_2[ii].rho);
				#ifdef ALE
				st_2[ii].xi = pn[ii].xi + 0.5 * dt * (pn[ii].v + pn[ii].vPert);
				#else
				st_2[ii].xi = pn[ii].xi + 0.5 * dt * pn[ii].v;
				#endif

				st_2[ii].v = pn[ii].v + 0.5 * dt * res_1[ii];
				st_2[ii].rho = pn[ii].rho + 0.5 * dt * Rrho_1[ii];
				st_2[ii].p = B*(pow(st_2[ii].rho/fvar.rho0,gam)-1);
				// st_2[ii].p = fvar.Cs*fvar.Cs * (st_2[ii].rho - fvar.rho0);
				st_2[ii].acc = res_1[ii];
				st_2[ii].Af = Af[ii];
				st_2[ii].Rrho = Rrho_1[ii];

				if(svar.Asource == 2 && st_2[ii].b == FREE)
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
			else
			{
				st_2[ii].xi = pn[ii].xi + dt * pn[ii].v;
				st_2[ii].rho = pn[ii].rho + 0.5 * dt * Rrho_1[ii];
				st_2[ii].p = B*(pow(st_2[ii].rho/fvar.rho0,gam)-1);
			}
		}

		/* Do the buffer particles */
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < svar.back.size(); ++ii)
		{
			for (size_t jj = 0; jj < 4; ++jj)
			{ /* Define buffer off the back particle */
				size_t const &buffID = svar.buffer[ii][jj];

				real frac = std::min(1.0, real(jj + 1) / 3.0);

				st_2[buffID].rho = frac * fvar.rhoJ + (1.0 - frac) * (pn[buffID].rho + dt * Rrho_1[buffID]);
				st_2[buffID].p = frac * fvar.pPress + (1.0 - frac) * (B * (pow(st_2[buffID].rho / fvar.rho0, gam) - 1));

				st_2[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
			}
		}
	}

	return (Check_RK_Error(svar,start,end,error1,error2,logbase,pn,st_2,k));
}

/* <summary> Perform the rest of the Runge-Kutta integration, assuming frozen 
	dissipation terms. This will do step 2 to 4 (n+1/4 to n+1) </summary> */
void Perform_RK4(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, real const& B, real const& gam,
	MESH& cells, vector<size_t> const& cellsused, OUTL const& outlist, DELTAP const& dp, real& logbase,
	SPHState& pn, SPHState& st_2, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel, real& error1)
{
	/*Create the vectors*/		
	vector<StateVecD> res_2(end,StateVecD::Zero());
	vector<real> Rrho_2(end,0.0);

	SPHState st_3 = pnp1; /*Step 3*/
	vector<StateVecD> res_3(end,StateVecD::Zero());
	vector<real> Rrho_3(end,0.0);

	SPHState st_4 = pnp1; /*Step 4*/
	vector<StateVecD> res_4(end,StateVecD::Zero());
	vector<real> Rrho_4(end,0.0);

	vector<StateVecD> Af(end,StateVecD::Zero());
	vector<real> wDiff(end,0.0);
	vector<StateVecD> norm(end, StateVecD::Zero());

	Force = StateVecD::Zero();
	svar.AForce = StateVecD::Zero();

	const real dt = svar.dt;
	
	// real error1 = 0.0;
	real error2 = error1;
	// logbase = 0.0;
	uint k = 2;


	/********************************************************************/
	/*************************  STEP 2  *********************************/
	/********************************************************************/
	Get_Boundary_Pressure(fvar,start,outlist,st_2);
	Forces(svar, fvar, avar, cells, st_2, outlist, dp, res_2, Rrho_2, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,st_2,res_2,Rrho_2,st_3,res_3,Rrho_3)
	{
		
		// #pragma omp for schedule(static) nowait
		// for (size_t ii=0; ii < start; ++ii)
		// {	/****** BOUNDARY PARTICLES ***********/
		// 	st_3[ii].rho = pn[ii].rho+0.5*dt*(Rrho_2[ii]);
		// 	st_3[ii].p = B*(pow(st_3[ii].rho/fvar.rho0,gam)-1);
		// 	// st_3[ii].p = fvar.Cs*fvar.Cs * (st_3[ii].rho - fvar.rho0);
		// }

		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* BUFFER = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			if(st_3[ii].b > BUFFER)
			{	
				#ifdef ALE
				st_3[ii].xi = pn[ii].xi + 0.5 * dt * (st_2[ii].v + st_2[ii].vPert);
				#else
				st_3[ii].xi = pn[ii].xi + 0.5 * dt * st_2[ii].v;
				#endif

				st_3[ii].v = pn[ii].v + 0.5 * dt * res_2[ii];
				real const rho = pn[ii].rho + 0.5 * dt * Rrho_2[ii];
				st_3[ii].rho = rho;
				st_3[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);

			}

		}

		/* Do the buffer particles */
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < svar.back.size(); ++ii)
		{
			for (size_t jj = 0; jj < 4; ++jj)
			{ /* Define buffer off the back particle */
				size_t const &buffID = svar.buffer[ii][jj];

				real const frac = std::min(1.0, real(jj + 1) / 3.0);
				real const rho = frac * fvar.rhoJ + (1.0 - frac) * (pn[buffID].rho + dt * Rrho_2[buffID]); 
				st_3[buffID].rho = rho;
				st_3[buffID].p = frac * fvar.pPress + (1.0 - frac) * pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);

				st_3[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
			}
		}
	}

	// if (Check_RK_Error(svar,start,end,error1,error2,logbase,st_2,st_3,k) < 0)
	// 	return -1;

	k++;

	/********************************************************************/
	/*************************  STEP 3  *********************************/
	/********************************************************************/
	#ifdef ALE
		if(svar.ghost > 0)
			Particle_Shift_Ghost(svar,fvar,start,end,outlist,dp,st_3);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end,outlist,dp,st_3);
	#endif

	Get_Boundary_Pressure(fvar,start,outlist,st_3);
	Forces(svar, fvar, avar, cells, st_3, outlist, dp, res_3, Rrho_3, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,st_3,res_3,Rrho_3,st_4,res_4,Rrho_4)
	{
		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			/* BUFFER = pipe particle receiving prescribed motion            */
			/* BACK = the latest particle in that column                    */
			/* PIPE = in the pipe, with free motion                         */
			/* FREE = free of the pipe and receives an aero force           */

			/*Check if the particle is clear of the starting area*/
			// if(st_3[ii].b == START || st_3[ii].b == BACK)
			// {   
			// 	/*For the particles marked 1, perform a prescribed motion*/
			// 	st_4[ii].xi = pn[ii].xi + dt*pn[ii].v;
			// }

			if(st_4[ii].b > BUFFER)
			{	
				#ifdef ALE
				st_4[ii].xi = pn[ii].xi + dt * (st_3[ii].v + st_3[ii].vPert);
				#else
				st_4[ii].xi = pn[ii].xi + dt * st_3[ii].v;
				#endif

				st_4[ii].v = pn[ii].v + dt * res_3[ii];
				real const rho = pn[ii].rho + dt * Rrho_3[ii];
				st_4[ii].rho = rho;
				st_4[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);
				// st_4[ii].p = fvar.Cs*fvar.Cs * (st_4[ii].rho - fvar.rho0);
			}
		}

		/* Do the buffer particles */
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < svar.back.size(); ++ii)
		{
			for (size_t jj = 0; jj < 4; ++jj)
			{ /* Define buffer off the back particle */
				size_t const &buffID = svar.buffer[ii][jj];

				real frac = std::min(1.0, real(jj + 1) / 3.0);
				real const rho = frac * fvar.rhoJ + (1.0 - frac) * (pn[buffID].rho + dt * Rrho_3[buffID]);
				st_4[buffID].rho = rho;
				st_4[buffID].p = frac * fvar.pPress + (1.0 - frac) * pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);
				st_4[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
			}
		}
	}

	// if (Check_RK_Error(svar,start,end,error1,error2,logbase,st_3,st_4,k) < 0)
	// 	return -1;

	k++;
	/********************************************************************/
	/*************************  STEP 4  *********************************/
	/********************************************************************/
	#ifdef ALE
		if(svar.ghost > 0)
			Particle_Shift_Ghost(svar,fvar,start,end,outlist,dp,st_4);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end,outlist,dp,st_4);
	#endif

	Get_Boundary_Pressure(fvar,start,outlist,st_4);
	Forces(svar, fvar, avar, cells, st_4, outlist, dp, res_4, Rrho_4, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,pnp1,st_2,res_2,Rrho_2,st_3,res_3,Rrho_3,st_4,res_4,Rrho_4)
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


		#pragma omp for schedule(static) nowait
		for (size_t ii=start; ii < end; ++ii)
		{	/****** FLUID PARTICLES **************/
			
			if(pnp1[ii].b > BUFFER)
			{	
				if(pnp1[ii].b == PIPE)
				{	/*Do a check to see if it needs to be given an aero force*/
					StateVecD vec = svar.Transp*(st_4[ii].xi-svar.sim_start);
					if(vec(1) > 0.0)
					{
						pnp1[ii].b = FREE;
						if(svar.Asource == 1 || svar.Asource == 2)
						{
							/*retrieve the cell it's in*/
							uint to_del = 0;
							FirstCell(svar, TREE.CELL, cells, pnp1[ii], to_del);

							if (to_del)
							{
								#pragma omp critical
								{
									// cout << "Particle has crossed an outer boundary!" << endl;
									// cout << "Particle will be deleted." << endl;
									pnp1.erase(pnp1.begin() + ii);
									pn.erase(pn.begin() + ii);
								}
								#pragma omp atomic
								svar.totPts--;
								#pragma omp atomic
								end--;
							}
						}
					}	
				}

				#ifdef ALE
				pnp1[ii].xi = pn[ii].xi + (dt / 6.0) * 
							  ((pn[ii].v+pn[ii].vPert) +
								2.0 * (st_2[ii].v + st_2[ii].vPert) +
								2.0 * (st_3[ii].v + st_3[ii].vPert) +
								(st_4[ii].v + st_4[ii].vPert));
				#else
				pnp1[ii].xi = pn[ii].xi + (dt / 6.0) * 
						(pn[ii].v + 2.0 * st_2[ii].v + 2.0 * st_3[ii].v + st_4[ii].v);
				#endif

				pnp1[ii].v = pn[ii].v + 
					(dt/6.0) * (st_2[ii].acc + 2.0 * res_2[ii] + 2.0 * res_3[ii] + res_4[ii]);

				real const rho = pn[ii].rho + (dt / 6.0) * 
					(st_2[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii]);
				pnp1[ii].rho = rho;

				pnp1[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);

				pnp1[ii].acc = res_4[ii];
				pnp1[ii].Af = Af[ii];
				pnp1[ii].Rrho = Rrho_4[ii];

				if(svar.Asource == 2 && pnp1[ii].b == FREE)
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

			pnp1[ii].s = outlist[ii].size();
			pnp1[ii].bNorm = norm[ii];
		}

		/* Do the buffer particles */
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < svar.back.size(); ++ii)
		{
			for (size_t jj = 0; jj < 4; ++jj)
			{ /* Define buffer off the back particle */
				size_t const &buffID = svar.buffer[ii][jj];

				real frac = std::min(1.0, real(jj + 1) / 3.0);

				pnp1[buffID].Rrho = Rrho_4[buffID];

				real const rho = frac * fvar.rhoJ + (1.0 - frac) * 
								(pn[buffID].rho + (dt / 6.0) * 
								(st_2[buffID].Rrho + 2.0 * Rrho_2[buffID] + 
								 2.0 * Rrho_3[buffID] + Rrho_4[buffID]));
				pnp1[buffID].rho = rho;

				pnp1[buffID].p = frac * fvar.pPress + (1.0 - frac) * pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);

				pnp1[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
			}
		}

	}/*End pragma omp parallel*/

	Check_RK_Error(svar, start, end, error1, error2, logbase, st_4, pnp1, k);
}

real Runge_Kutta4(KDTREE const& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, real const& B, real const& gam,
	MESH& cells, vector<size_t> const& cellsused, OUTL const& outlist, DELTAP const& dp, real& logbase,
	SPHState& pn, SPHState& st_2, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel)
{
	real error = 0.0;
	// int errstate = 1;
	// while(errstate != 0)
	// {
		Perform_RK4(TREE,svar,fvar,avar,start,end,B,gam,cells,cellsused,
					outlist,dp,logbase,pn,st_2,pnp1,Force,dropVel,error);

	// 	if(errstate < 0)
	// 	{
	// 		/*Timestep was unstable - reduce the timestep*/
	// 		svar.dt *= 0.5;
	// 	}
	// }

	return error;
}

#endif