/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "Runge_Kutta.h"
#include "Kernel.h"
#include "Resid.h"
#include "Containment.h"
#include "Shifting.h"

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
real Get_First_RK(SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start,
	 size_t const& end, real const& B, real const& gam, real const& npd, MESH const& cells,
	LIMITS const& limits, OUTL const& outlist, real& logbase, SPHState& pn, SPHState& st_2, real& error1)
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

	for(size_t block = 0; block < svar.nbound; block++)
	{

		if(limits[block].nTimes != 0)
		{
			// Get the current boundary velocity
			StateVecD vel = StateVecD::Zero();
			for(size_t time = 0; time < limits[block].nTimes; time++)
			{
				if(limits[block].times[time] > svar.t)
				{
					vel = limits[block].vels[time];
				}
			}

			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_2[jj].v = vel;
			}
		}
		else
		{
			#pragma omp parallel for
			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_2[jj].v = limits[block].vels[0];
			}
		}

		if(limits[block].no_slip)
			Set_No_Slip(fvar,limits[block].index.first,limits[block].index.second,outlist,st_2);

		switch (limits[block].bound_solver)
		{
			case DBC:
			{
				Boundary_DBC(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_2,res_1);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
								pn[ii].rho + 0.5 * dt * Rrho_1[ii]));
					st_2[ii].rho = rho;
					st_2[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
					st_2[ii].Rrho = Rrho_1[ii];
				}
				break;
			}
			case pressure_G:
			{
				Get_Boundary_Pressure(svar.grav,fvar,limits[block].index.first,
								limits[block].index.second,outlist,st_2);
				break;
			}
			case ghost:
			{
				vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
				Boundary_Ghost(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_2,Rrho_1,near_inlet);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					if(!near_inlet[ii - limits[block].index.first])
					{
						real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
									pn[ii].rho + 0.5 * dt * Rrho_1[ii]));
						st_2[ii].rho = rho;
						st_2[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						st_2[ii].Rrho = Rrho_1[ii];
					}
					else
					{	// Don't allow negative pressures
						real const rho = std::max(fvar.rho0, std::min(fvar.rhoMax, 
									pn[ii].rho + 0.5 * dt * Rrho_1[ii]));
						st_2[ii].rho = rho;
						st_2[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						st_2[ii].Rrho = fmax(0.0,Rrho_1[ii]);
					}
				}
				break;
			}
		}
	}

	Forces(svar, fvar, avar, cells, st_2, outlist, npd, res_1, Rrho_1, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,st_2) /*reduction(+:Force,dropVel)*/
	{
		/********************************************************************/
		/*************************  STEP 1  *********************************/
		/********************************************************************/
		for(size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
		{
			#pragma omp for nowait 
			for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
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
                    real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
                                    pn[ii].rho + 0.5 * dt * Rrho_1[ii]));
                    st_2[ii].rho = rho;
                    st_2[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
                    // st_2[ii].p = fvar.Cs*fvar.Cs * (st_2[ii].rho - fvar.rho0);
                    st_2[ii].acc = res_1[ii];
                    st_2[ii].Af = Af[ii];
                    st_2[ii].Rrho = Rrho_1[ii];
                }
                else
                {
                    st_2[ii].xi = pn[ii].xi + dt * pn[ii].v;
                    st_2[ii].rho = pn[ii].rho + 0.5 * dt * Rrho_1[ii];
                    st_2[ii].p = B*(pow(st_2[ii].rho/fvar.rho0,gam)-1);
                }
            }
            
			/* Do the buffer particles */
			if(limits[block].block_type == inletZone)
			{
				switch(limits[block].fixed_vel_or_dynamic)
				{
				case 1:
				{
					StateVecD unorm = limits[block].insert_norm.normalized();
					#pragma omp for schedule(static) nowait
					for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
					{	
						size_t const& backID = limits[block].back[ii];
						StateVecD const& xi = st_2[backID].xi;

						for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
						{	/* Define buffer off the back particle */
							size_t const& buffID = limits[block].buffer[ii][jj];
							// Set position as related to the previous particle.
							st_2[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

							// How to set density and pressure though?
							st_2[buffID].v = st_2[backID].v;
							st_2[buffID].rho = st_2[backID].rho;
							st_2[buffID].p = st_2[backID].p; 
						}	
					}
					break;
				}
				case 0:
				{
					#pragma omp for schedule(static) nowait
					for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
					{	
						// size_t const& backID = svar.back[ii];
						for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
						{	/* Define buffer off the back particle */

							size_t const& buffID = limits[jj].buffer[ii][jj];

							st_2[buffID].Rrho = Rrho_1[buffID];
							
							real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
											pn[buffID].rho+0.5*dt*Rrho_1[buffID])); 
							st_2[buffID].rho = rho;
							st_2[buffID].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
							st_2[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
						}	
					}
					break;
				}
				}
			}
        }
	}

	return (Check_RK_Error(svar,start,end,error1,error2,logbase,pn,st_2,k));
}

/* <summary> Perform the rest of the Runge-Kutta integration, assuming frozen 
	dissipation terms. This will do step 2 to 4 (n+1/4 to n+1) </summary> */
void Perform_RK4(Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, real const& B, real const& gam, real const& npd, MESH const& cells, 
	LIMITS const& limits, OUTL const& outlist, real& logbase,
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
	for(size_t block = 0; block < svar.nbound; block++)
	{
		if(limits[block].nTimes != 0)
		{
			// Get the current boundary velocity
			StateVecD vel = StateVecD::Zero();
			for(size_t time = 0; time < limits[block].nTimes; time++)
			{
				if(limits[block].times[time] > svar.t)
				{
					vel = limits[block].vels[time];
				}
			}

			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_2[jj].v = vel;
			}
		}
		else
		{
			#pragma omp parallel for
			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_2[jj].v = limits[block].vels[0];
			}
		}

		if(limits[block].no_slip)
			Set_No_Slip(fvar,limits[block].index.first,limits[block].index.second,outlist,st_2);

		switch (limits[block].bound_solver)
		{
			case DBC:
			{
				Boundary_DBC(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_2,res_2);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
								pn[ii].rho + 0.5 * dt * Rrho_2[ii]));
					st_2[ii].rho = rho;
					st_2[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
					st_2[ii].Rrho = Rrho_2[ii];
				}
				break;
			}
			case pressure_G:
			{
				Get_Boundary_Pressure(svar.grav,fvar,limits[block].index.first,
								limits[block].index.second,outlist,st_2);
				break;
			}
			case ghost:
			{
				vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
				Boundary_Ghost(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_2,Rrho_2,near_inlet);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					if(!near_inlet[ii - limits[block].index.first])
					{
						real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
									pn[ii].rho + 0.5 * dt * Rrho_2[ii]));
						st_3[ii].rho = rho;
						st_3[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						st_3[ii].Rrho = Rrho_2[ii];
					}
					else
					{	// Don't allow negative pressures
						real const rho = std::max(fvar.rho0, std::min(fvar.rhoMax, 
									pn[ii].rho + 0.5 * dt * Rrho_2[ii]));
						st_3[ii].rho = rho;
						st_3[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						st_3[ii].Rrho = fmax(0.0,Rrho_2[ii]);
					}
				}
				break;
			}
		}
	}

	Forces(svar, fvar, avar, cells, st_2, outlist, npd, res_2, Rrho_2, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,st_2,res_2,Rrho_2,st_3,res_3,Rrho_3)
	{	
		for(size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
		{
			#pragma omp for nowait 
			for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
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
                    real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
                                    pn[ii].rho + 0.5 * dt * Rrho_2[ii]));
                    st_3[ii].rho = rho;
                    st_3[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
                }
            }
			
			/* Do the buffer particles */
			if(limits[block].block_type == inletZone)
			{
				switch (limits[block].fixed_vel_or_dynamic)
				{
					case 1:
					{
						StateVecD unorm = limits[block].insert_norm.normalized();
						#pragma omp for schedule(static) nowait
						for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
						{	
							size_t const& backID = limits[block].back[ii];
							StateVecD const& xi = st_3[backID].xi;

							for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
							{	/* Define buffer off the back particle */
								size_t const& buffID = limits[block].buffer[ii][jj];
								// Set position as related to the previous particle.
								st_3[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

								// How to set density and pressure though?
								st_3[buffID].v = st_3[backID].v;
								st_3[buffID].rho = st_3[backID].rho;
								st_3[buffID].p = st_3[backID].p; 
							}	
						}
						break;
					}
					case 0:
					{
						#pragma omp for schedule(static) nowait
						for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
						{	
							// size_t const& backID = svar.back[ii];
							for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
							{	/* Define buffer off the back particle */

								size_t const& buffID = limits[jj].buffer[ii][jj];

								st_3[buffID].Rrho = Rrho_2[buffID];
								
								real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
									pn[buffID].rho + 0.5 * dt * Rrho_2[buffID]));
								st_3[buffID].rho = rho;
								st_3[buffID].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);

								st_3[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
							}	
						}
						break;
					}
				}
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
			Particle_Shift_Ghost(svar,fvar,start,end,outlist,st_3);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end,outlist,st_3);
	#endif

	for(size_t block = 0; block < svar.nbound; block++)
	{
		if(limits[block].nTimes != 0)
		{
			// Get the current boundary velocity
			StateVecD vel = StateVecD::Zero();
			for(size_t time = 0; time < limits[block].nTimes; time++)
			{
				if(limits[block].times[time] > svar.t)
				{
					vel = limits[block].vels[time];
				}
			}

			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_3[jj].v = vel;
			}
		}
		else
		{
			#pragma omp parallel for
			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_3[jj].v = limits[block].vels[0];
			}
		}
		if(limits[block].no_slip)
			Set_No_Slip(fvar,limits[block].index.first,limits[block].index.second,outlist,st_3);

		switch (limits[block].bound_solver)
		{
			case DBC:
			{
				Boundary_DBC(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_3,res_3);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
								pn[ii].rho + 0.5 * dt * Rrho_3[ii]));
					st_4[ii].rho = rho;
					st_4[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
					st_4[ii].Rrho = Rrho_3[ii];
				}
				break;
			}
			case pressure_G:
			{
				Get_Boundary_Pressure(svar.grav,fvar,limits[block].index.first,
								limits[block].index.second,outlist,st_3);
				break;
			}
			case ghost:
			{
				vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
				Boundary_Ghost(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_3,Rrho_3,near_inlet);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					if(!near_inlet[ii - limits[block].index.first])
					{
						real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
									pn[ii].rho + 0.5 * dt * Rrho_3[ii]));
						st_4[ii].rho = rho;
						st_4[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						st_4[ii].Rrho = Rrho_3[ii];
					}
					else
					{	// Don't allow negative pressures
						real const rho = std::max(fvar.rho0, std::min(fvar.rhoMax, 
									pn[ii].rho + 0.5 * dt * Rrho_3[ii]));
						st_4[ii].rho = rho;
						st_4[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						st_4[ii].Rrho = fmax(0.0,Rrho_3[ii]);
					}
				}
				break;
			}
		}
	}

	Forces(svar, fvar, avar, cells, st_3, outlist, npd, res_3, Rrho_3, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,st_3,res_3,Rrho_3,st_4,res_4,Rrho_4)
	{
		for(size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
		{
			#pragma omp for nowait 
			for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
			{	/****** FLUID PARTICLES **************/
				if(st_4[ii].b > BUFFER)
				{	
					#ifdef ALE
					st_4[ii].xi = pn[ii].xi + dt * (st_3[ii].v + st_3[ii].vPert);
					#else
					st_4[ii].xi = pn[ii].xi + dt * st_3[ii].v;
					#endif

					st_4[ii].v = pn[ii].v + dt * res_3[ii];
					real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, pn[ii].rho + dt * Rrho_3[ii]));
					st_4[ii].rho = rho;
					st_4[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
					// st_4[ii].p = fvar.Cs*fvar.Cs * (st_4[ii].rho - fvar.rho0);
				}
			}

			/* Do the buffer particles */
			if(limits[block].block_type == inletZone)
			{
				switch (limits[block].fixed_vel_or_dynamic)
				{
					case 1:
					{
						StateVecD unorm = limits[block].insert_norm.normalized();
						#pragma omp for schedule(static) nowait
						for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
						{	
							size_t const& backID = limits[block].back[ii];
							StateVecD const& xi = st_4[backID].xi;

							for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
							{	/* Define buffer off the back particle */
								size_t const& buffID = limits[block].buffer[ii][jj];
								// Set position as related to the previous particle.
								st_4[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

								// How to set density and pressure though?
								st_4[buffID].v = st_4[backID].v;
								st_4[buffID].rho = st_4[backID].rho;
								st_4[buffID].p = st_4[backID].p; 
							}	
						}
						break;
					}
					case 0:
					{
						#pragma omp for schedule(static) nowait
						for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
						{	
							// size_t const& backID = svar.back[ii];
							for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
							{	/* Define buffer off the back particle */

								size_t const& buffID = limits[jj].buffer[ii][jj];

								st_4[buffID].Rrho = Rrho_3[buffID];
								
								real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
									pn[buffID].rho + 0.5 * dt * Rrho_3[buffID]));
								st_4[buffID].rho = rho;
								st_4[buffID].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);

								st_4[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
							}	
						}
						break;
					}
				}
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
			Particle_Shift_Ghost(svar,fvar,start,end,outlist,st_4);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end,outlist,st_4);
	#endif

	for(size_t block = 0; block < svar.nbound; block++)
	{
		if(limits[block].nTimes != 0)
		{
			// Get the current boundary velocity
			StateVecD vel = StateVecD::Zero();
			for(size_t time = 0; time < limits[block].nTimes; time++)
			{
				if(limits[block].times[time] > svar.t)
				{
					vel = limits[block].vels[time];
				}
			}

			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_4[jj].v = vel;
			}
		}
		else
		{
			#pragma omp parallel for
			for(size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
			{
				st_4[jj].v = limits[block].vels[0];
			}
		}

		if(limits[block].no_slip)
			Set_No_Slip(fvar,limits[block].index.first,limits[block].index.second,outlist,st_4);

		switch (limits[block].bound_solver)
		{
			case DBC:
			{
				Boundary_DBC(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_4,res_4);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
								pn[ii].rho + (dt / 6.0) * (st_2[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii])));
					pnp1[ii].rho = rho;
					pnp1[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
					pnp1[ii].Rrho = Rrho_2[ii];
				}
				break;
			}
			case pressure_G:
			{
				Get_Boundary_Pressure(svar.grav,fvar,limits[block].index.first,
								limits[block].index.second,outlist,st_4);
				break;
			}
			case ghost:
			{
				vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
				Boundary_Ghost(fvar,limits[block].index.first,limits[block].index.second,
								outlist,st_4,Rrho_4,near_inlet);

				#pragma omp for schedule(static) nowait
				for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
				{	/****** BOUNDARY PARTICLES ***********/
					if(!near_inlet[ii - limits[block].index.first])
					{
						real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, 
									pn[ii].rho + (dt / 6.0) * (st_2[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii])));
						pnp1[ii].rho = rho;
						pnp1[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						pnp1[ii].Rrho = Rrho_4[ii];
					}
					else
					{	// Don't allow negative pressures
						real const rho = std::max(fvar.rho0, std::min(fvar.rhoMax, 
								pn[ii].rho + (dt / 6.0) * (st_2[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii])));
						pnp1[ii].rho = rho;
						pnp1[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
						pnp1[ii].Rrho = fmax(0.0,Rrho_4[ii]);
					}
				}
				break;
			}
		}

	}

	Forces(svar, fvar, avar, cells, st_4, outlist, npd, res_4, Rrho_4, Af, Force);

	#pragma omp parallel default(shared) // shared(svar,pn,pnp1,st_2,res_2,Rrho_2,st_3,res_3,Rrho_3,st_4,res_4,Rrho_4)
	{
		for(size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
		{
			#pragma omp for nowait 
			for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
			{	/****** FLUID PARTICLES **************/
				if(pnp1[ii].b > BUFFER && pnp1[ii].b != OUTLET)
				{	
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

					real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, pn[ii].rho + (dt / 6.0) * 
						(st_2[ii].Rrho + 2.0 * Rrho_2[ii] + 2.0 * Rrho_3[ii] + Rrho_4[ii])));
					pnp1[ii].rho = rho;

					pnp1[ii].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);

					pnp1[ii].acc = res_4[ii];
					pnp1[ii].Af = Af[ii];
					pnp1[ii].Rrho = Rrho_4[ii];
				}
				else if (pnp1[ii].b == OUTLET)
				{	/* For the outlet zone, just perform euler integration of last info */
					pnp1[ii].xi = pn[ii].xi + dt*pnp1[ii].v;
				}
			}

			/* Do the buffer particles */
			if(limits[block].block_type == inletZone)
			{
				switch (limits[block].fixed_vel_or_dynamic)
				{
					case 1:
					{
						StateVecD unorm = limits[block].insert_norm.normalized();
						#pragma omp for schedule(static) nowait
						for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
						{	
							size_t const& backID = limits[block].back[ii];
							StateVecD const& xi = pnp1[backID].xi;

							for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
							{	/* Define buffer off the back particle */
								size_t const& buffID = limits[block].buffer[ii][jj];
								// Set position as related to the previous particle.
								pnp1[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

								// How to set density and pressure though?
								pnp1[buffID].v = pnp1[backID].v;
								pnp1[buffID].rho = pnp1[backID].rho;
								pnp1[buffID].p = pnp1[backID].p; 
							}	
						}
						break;
					}
					case 0:
					{
						#pragma omp for schedule(static) nowait
						for(size_t ii = 0; ii < limits[block].back.size(); ++ii)
						{	
							// size_t const& backID = svar.back[ii];
							for(size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
							{	/* Define buffer off the back particle */

								size_t const& buffID = limits[jj].buffer[ii][jj];

								st_4[buffID].Rrho = Rrho_4[buffID];
								
								real const rho = std::max(fvar.rhoMin, std::min(fvar.rhoMax, pn[buffID].rho + (dt / 6.0) * 
									(st_2[buffID].Rrho + 2.0 * Rrho_2[buffID] + 2.0 * Rrho_3[buffID] + Rrho_4[buffID])));
								pnp1[buffID].rho = rho;
								pnp1[buffID].p = pressure_equation(rho,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
								pnp1[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
							}	
						}
						break;
					}
				}
			}
		}
	}/*End pragma omp parallel*/

	Check_RK_Error(svar, start, end, error1, error2, logbase, st_4, pnp1, k);
}

real Runge_Kutta4(Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, real const& B, real const& gam, real const& npd, MESH const& cells,
	LIMITS const& limits, OUTL const& outlist, real& logbase,
	SPHState& pn, SPHState& st_2, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel)
{
	real error = 0.0;
	// int errstate = 1;
	// while(errstate != 0)
	// {
		Perform_RK4(CELL_TREE,svar,fvar,avar,start,end,B,gam,npd,cells,
					limits,outlist,logbase,pn,st_2,pnp1,Force,dropVel,error);

	// 	if(errstate < 0)
	// 	{
	// 		/*Timestep was unstable - reduce the timestep*/
	// 		svar.dt *= 0.5;
	// 	}
	// }

	return error;
}
