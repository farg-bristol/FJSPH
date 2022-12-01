
/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/*** Continuity:        Delta-SPH. Marrone et al. (2011)                                  ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)                                  ***/
/*** Smoothing Kernel:  Wendland's C2                                                     ***/
/*** Integrator:        Newmark-Beta or 4th order Runge-Kutta                             ***/
/*** Variable Timestep Criteria: CFL + Monaghan, J.J. (1989) conditions                   ***/

#include "Resid.h"
#include "Kernel.h"
#include "Aero.h"


/* Boundary pressure calculation - Adami, Hu, and Adams, 2012 - https://doi.org/10.1016/j.jcp.2012.05.005*/
void Get_Boundary_Pressure(StateVecD const& grav, FLUID const& fvar,
	size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1)
{
	std::vector<real> pressure(end-start, 0.0);
	#pragma omp parallel for schedule(static) default(shared)/*Reduction defs in Var.h*/
	for (size_t ii=start; ii < end; ++ii)
	{
		SPHPart const& pi = pnp1[ii];

		int isNearSurface = 0;
		real kernsum = 0.0;
		real pkern = 0.0;
		StateVecD acckern = StateVecD::Zero();

		for (std::pair<size_t,real> const& jj : outlist[ii])
		{	/* Neighbour list loop. */
			SPHPart const& pj = pnp1[jj.first];
			/*Check if the neighbour is a fluid particle*/
			if(pj.b > PISTON)
			{
				StateVecD const Rji = pi.xi-pj.xi;
				
				real const rr = jj.second;
				real const r = sqrt(rr);
				real const volj = pj.m/pj.rho;
				real kern = volj * Kernel(r,fvar.H,fvar.correc);

				kernsum += kern;
				pkern += pj.p * kern;
				acckern +=  kern * pj.rho * Rji;
				if(pj.surfzone)	// Check if nearby to a surface
					isNearSurface = 1;
			}
		}/*End of neighbours*/

		if(kernsum > 0.0)
		{
			real p = (pkern+(grav-pi.acc).dot(acckern))/kernsum;
			
			if(isNearSurface)	// Disallow negative pressures near surfaces. Causes problems
				pressure[ii-start] = std::max(0.0,p);
			else
				pressure[ii-start] = p;
		}
			
	} /*End of boundary parts*/

	#pragma omp parallel for default(shared)/*Reduction defs in Var.h*/
	for (size_t ii=start; ii < end; ++ii)
	{	
		pnp1[ii].p = pressure[ii-start];
		pnp1[ii].rho = density_equation(pressure[ii-start],fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
	}

}

void Boundary_DBC( FLUID const& fvar, size_t const& start, size_t const& end,
	 OUTL const& outlist, SPHState& pnp1, vector<StateVecD>& RV)
{
	real const c2 = fvar.Cs*fvar.Cs;
	/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/		
	#pragma omp for schedule(static) nowait reduction(+:RV) /*Reduction defs in Var.h*/
	for (size_t ii=start; ii < end; ++ii)
	{
		SPHPart const& pi = pnp1[ii];
		

		for (std::pair<size_t,real> const& jj : outlist[ii])
		{	/* Neighbour list loop. */
			SPHPart const& pj = pnp1[jj.first];
			if(ii == jj.first || pj.b < BUFFER)
				continue;

			StateVecD const Rji = pj.xi-pi.xi;
			StateVecD const Vji = pj.v-pi.v;
			real const rr = jj.second;
			real const r = sqrt(rr);
			real const idist2 = 1.0/(rr + 0.001*fvar.HSQ);
			// real const volj = pj.m/pj.rho;
			real const kern = Kernel(r,fvar.H,fvar.correc);
			StateVecD const gradK = /*dp.L[ii] * */GradK(Rji, r,fvar.H, fvar.correc);
			StateVecD artViscI = pj.m*ArtVisc(fvar.nu,pi,pj,fvar,Rji,Vji,idist2,gradK);
			/*Base WCSPH continuity drho/dt*/
			RV[jj.first] -= 2.0*c2* kern / pow(kern + fvar.correc,2)*gradK + artViscI;
		}/*End of neighbours*/
	} /*End of boundary parts*/
}

/* Basic ghost particles. Basically just static particles, with the same treatment as fluid particles in all equations */
/* Just don't have the position, velocity, and acceleration updated */
void Boundary_Ghost( FLUID const& fvar, size_t const& start, size_t const& end,
	 OUTL const& outlist, SPHState& pnp1, vector<real>& Rrho, vector<int>& near_inlet)
{
	/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/		
	#pragma omp for schedule(static) nowait /*Reduction defs in Var.h*/
	for (size_t ii=start; ii < end; ++ii)
	{
		SPHPart const& pi = pnp1[ii];
		real Rrhoi = 0.0;
		near_inlet[ii-start] = 0;
		for (std::pair<size_t,real> const& jj : outlist[ii])
		{	/* Neighbour list loop. */
			SPHPart const& pj = pnp1[jj.first];
			if(ii == jj.first)
				continue;
			
			if (pj.b == BACK || pj.b == BUFFER)
				near_inlet[ii-start] = 1;

			StateVecD const Rji = pj.xi-pi.xi;
			StateVecD const Vji = pj.v-pi.v;
			real const rr = jj.second;
			real const r = sqrt(rr);
			real const volj = pj.m/pj.rho;
			StateVecD const gradK = /*dp.L[ii] * */GradK(Rji, r,fvar.H, fvar.correc);

			/*Base WCSPH continuity drho/dt*/
			Rrhoi -= volj*Vji.dot(gradK);
		}/*End of neighbours*/
		
		/* No dissipation terms for the boundary */
		Rrho[ii] = Rrhoi*pi.rho;
	} /*End of boundary parts*/
}

void Set_No_Slip( FLUID const& fvar, size_t const& start, size_t const& end,
	 OUTL const& outlist, SPHState& pnp1)
{
	#pragma omp for schedule(static) nowait /*Reduction defs in Var.h*/
	for (size_t ii=start; ii < end; ++ii)
	{
		StateVecD velsum = StateVecD::Zero();
		real kernsum = 0.0;
		for (std::pair<size_t,real> const& jj : outlist[ii])
		{	/* Neighbour list loop. */
			SPHPart const& pj = pnp1[jj.first];
			if(ii == jj.first || pj.b <= PISTON)
				continue;
			real r = sqrt(jj.second);
			real kern = Kernel(r,fvar.H,fvar.correc);
			kernsum += kern;
			velsum += pj.v * kern;
		}

		if(kernsum > 0.0)
			pnp1[ii].v = 2.0 * pnp1[ii].v - velsum/kernsum;
	}
}

inline StateVecD PressureContrib(SPHPart const& pi, SPHPart const& pj, 
							real const& volj, StateVecD const& gradK)
{
	#ifdef LINEAR
		#ifndef TIC
		return pressPos(volj,pi,pj,gradK);
		#else
		if(pi.surfzone != 1)
			return pressTIC(volj,pi,pj,gradK);
		else
			return pressPos(volj,pi,pj,gradK);
		#endif
	#else
		#ifndef TIC
		return BasePos(pi,pj,gradK);
		#else
		if(pi.surfzone != 1)
			return BaseTIC(pi,pj,gradK);
		else
			return BasePos(pi,pj,gradK);
		#endif
	#endif
}

inline StateVecD SurfTenContrib(SPHPart const& pi, SPHPart const& pj, real const& volj, StateVecD const& gradK,
					StateVecD const& Rji, real const& r, real const& h, 
					real const& npdm2, real const& pi3o4)
{
		/* He et al (2014) - Diffuse Interface model */
		#ifdef HESF
		(void) Rji;
		(void) r;
		(void) h;
		(void) npdm2;
		(void) pi3o4;
		if( dp.lam[ii] < 0.75 /* pi.surf == 1 */ )
		{
			return HeST(pi.norm/* /dp.colour[ii] */,pj.norm/* /dp.colour[jj.first] */,volj*gradK);
		}
		#endif

		/* Nair & Poeschel (2017) - Pairwise force */
		#ifdef PAIRWISE
		(void) volj;
		(void) gradK;
		#ifdef ALE
		// if(pi.surfzone == 1)	
		#endif
			// surfT += SurfaceTens(Rji, r, fvar.H, fvar.sig, lam, npdm2, pi3o4, pi.b, pj.b, voli, volj);
			return SurfaceTens(Rji, r, h, npdm2, pi3o4, pi.b, pj.b);
		#endif	
}

///**************** RESID calculation **************
void Forces(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, 
	SPHState const& pnp1, OUTL const& outlist,/*  DELTAP const& dp, */ real const& npd,
	 vector<StateVecD>& RV, vector<real>& Rrho, std::vector<StateVecD>& Af, StateVecD& Force)
{

	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;

	#ifdef PAIRWISE

	/*Surface tension factor*/
	// real const lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
	// 						(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));
	#if SIMDIM==2
	real const lam = 8.0/81.0*pow(fvar.H,4)/pow(M_PI,4)*(9.0/4.0*pow(M_PI,3)-6.0*M_PI-4.0);
	
	#else
	real const lam = 3.0 / (4.0 * pow(M_PI, 4.0)) * (pow(2.0, 7) - 9.0 * pow(2.0, 4) 
					* M_PI * M_PI + 9.0 * 3.0 * pow(M_PI, 4)) * pow(fvar.H / 3.0, 5);
	#endif
	
	real const npdm2 = (0.5 * fvar.sig / lam )/(npd*npd);
	real const pi3o4 = 3.0*M_PI/4.0;
	#endif

	#pragma omp parallel default(shared)
	{
		/*Gravity Vector*/
		StateVecD const& g = svar.grav;

		/******* Loop fluid points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Force,RV,Af) schedule(static) nowait
		for (size_t ii = start; ii < end; ++ii)
		{
			SPHPart const& pi = pnp1[ii];

			StateVecD Vdiff = StateVecD::Zero();
			real Pbasei = 0.0;

			StateVecD aero = StateVecD::Zero();
			StateVecD RVi = StateVecD::Zero();
			StateVecD viscI = StateVecD::Zero();
			StateVecD surfT = StateVecD::Zero();
			
			real Rrhoi = 0.0;
			#ifdef ALE
			StateVecD RVia = StateVecD::Zero();
			real Rrhoc = 0.0;
			#endif
			
			#ifdef NOFROZEN
			StateVecD artViscI = StateVecD::Zero();
			#ifndef NODSPH
			real Rrhod = 0.0;			
			#endif
			#endif

			if( pi.lam_ng < avar.cutoff /* pi.surf == 1 */ && pi.b == FREE )
			{
				Vdiff = pi.cellV - pi.v;
								
				aero = CalcAeroAcc(avar,pi,Vdiff,pi.norm,pi.lam_ng,Pbasei,svar.dx);
				RV[ii] += aero;
				Af[ii] += aero;
				Force += aero;
			}

			for (std::pair<size_t,real> const& jj : outlist[ii])
			{	/* Neighbour list loop. */
				SPHPart const& pj = pnp1[jj.first];
				/*Check if the position is the same, and skip the particle if yes*/
				if(ii == jj.first)
					continue;

				StateVecD const Rji = pj.xi-pi.xi;
				StateVecD const Vji = pj.v-pi.v;
				real const rr = jj.second;
				real const r = sqrt(rr);
				real const idist2 = 1.0/(rr + 0.001*fvar.HSQ);
				
				// #if defined(CSF) || defined(HESF) || (!defined(NODSPH) && defined(NOFROZEN))
					real const volj = pj.m/pj.rho;
				// #endif

				// StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);/* gradK;*/	
				
				StateVecD const contrib = PressureContrib(pi,pj,volj,gradK);

				#ifdef ALE
					StateVecD const ALEcontrib = ALEMomentum(pi, pj, volj, gradK);
				#endif

				/*Laminar Viscosity - Morris (2003)*/
				StateVecD const visc = Viscosity(fvar.nu,pi,pj,Rji,Vji,idist2,gradK);
				viscI += pj.m*visc;

				/* Surface tension options */
				if(/* pi.b == FREE  && */ pj.b != GHOST)
				{
					surfT += SurfTenContrib(pi,pj,volj,gradK,Rji,r,fvar.H,npdm2,pi3o4);
				}
				
				RVi -= contrib;
				#ifdef ALE
					RVia += ALEcontrib;
				#endif

				/*Base WCSPH continuity drho/dt*/
				#ifdef ALE
					Rrhoi -= ALEContinuity(pi,pj,volj,gradK);
					Rrhoc += ALECont2ndterm(pi,pj,volj,gradK);
				#else
					// Rrhoi -= pj.m*(Vji.dot(gradK));
					Rrhoi -= volj*Vji.dot(gradK);
				#endif

				#ifdef NOFROZEN
					if(pj.b > PISTON)
					{
						#ifndef NODSPH
							Rrhod += Continuity_dSPH(Rji,idist2,fvar.HSQ,gradK,volj,
												pi.gradRho,pj.gradRho,pi,pj);
						#endif
						#ifdef LINEAR
							artViscI += Linear_ArtVisc(Rji,Vji,idist2,gradK,volj);
						#else
							artViscI += pj.m*ArtVisc(fvar.nu,pi,pj,fvar,Rji,Vji,idist2,gradK);
						#endif
					}
					else
					{
						#ifndef NODSPH
							Rrhod += Continuity_dSPH(Rji,idist2,fvar.HSQ,gradK,volj,pi,pj);
						#endif
					}
				#endif
				
				/* Smooth aero force over the neighbours */
				// RV[jj.first] += aero*Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				// if(aero.norm()!= 0)
				// cout << aero[0] << "  " << aero[1] << "  " << Kernel(r,fvar.H,fvar.correc) << "  " << dp.kernsum[ii] << endl;
				

			}/*End of neighbours*/

			#ifdef LINEAR
			RVi /= pi.rho;
			#endif

			if(pi.internal == 1)
			{	// Apply the normal boundary force
				RVi += NormalBoundaryRepulsion(fvar, cells, pi);
			}

			#ifdef CSF
			// if(dp.lam[ii] < 0.75 /* pi.surf == 1 */)
			// {
			// 	RVi += (fvar.sig/pi.rho * curve * pi.norm)/correc;
			// }
			// if(pi.surf == 1)
			// {	/* lambda tuning factor supposedly = 3 */
			// 	RVi += 3.0*(fvar.sig * pi.curve * pi.norm); 
			// }
			// surfT -= fvar.sig * pi.curve * pi.norm * dp.colourG[ii];
			#ifdef ALE
			// if(pi.surfzone == 1)	
			#endif
				RVi += 1.0/pi.rho * fvar.sig * pi.curve * pi.norm * pi.colourG;
			#endif

			#ifdef HESF
			surfT *= (0.02/2.0)*(pi.m/pi.rho);
			#endif
			
			// #pragma omp critical
			// Force += aero;	
			
			#ifdef NOFROZEN
				#ifdef LINEAR
				artViscI *= fvar.alpha * fvar.H * fvar.Cs * fvar.rho0 / pi.rho;
				#endif

				#ifdef ALE
				RV[ii] += (RVi + RVia + artViscI + viscI + surfT/pi.m /* + aero/pi.m */ + g);
				#else
				RV[ii] += (RVi + artViscI + viscI + surfT/pi.m /* + aero/pi.m */ + g);
				#endif

				#ifndef NODSPH
					#ifdef ALE
					Rrho[ii] = Rrhoi*pi.rho + Rrhoc + fvar.dCont * Rrhod;
					#else
					Rrho[ii] = Rrhoi*pi.rho + fvar.dCont * Rrhod;
					#endif 
				#else
					Rrho[ii] = Rrhoi*pi.rho;
				#endif
				
			#else
				#ifdef ALE
				RV[ii] += (RVi + RVia + pnp1[ii].aVisc + viscI + surfT/pi.m /* + aero/pi.m */ + g);
				#else
				RV[ii] += (RVi + pnp1[ii].aVisc + viscI + surfT/pi.m /* + aero/pi.m */ + g);
				#endif
				

				#ifndef NODSPH
					#ifdef ALE
					Rrho[ii] = Rrhoi*pi.rho + Rrhoc + pnp1[ii].deltaD;
					#else
					Rrho[ii] = Rrhoi*pi.rho + pnp1[ii].deltaD;
					#endif 
				#else
					Rrho[ii] = Rrhoi*pi.rho;
				#endif

			#endif
			// res_.emplace_back(aero/pi.m);
			// Rrho_.emplace_back(0.0);

		} /*End of sim parts*/		

	}	/*End of declare parallel */
 
    // cout << RV.size() << "  " << Rrho.size() << "  " << Af.size() << "  " << curv.size() << endl;
}
