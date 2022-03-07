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

#ifndef RESID_H
#define RESID_H

#include "Var.h"
#include "Kernel.h"
#include "Aero.h"
// #include "Add.h"

#if !defined(HESF) && !defined(CSF)
#define PAIRWISE
#endif


#ifdef HESF
// /* Colour field gradient in He et al (2014) method, note the bottom variable to account for no air*/
// std::vector<StateVecD> GetColourGrad(SIM const& svar, FLUID const& fvar, SPHState const& pnp1, OUTL const& outlist)
// {
// 	std::vector<StateVecD> cgrad(svar.totPts, StateVecD::Zero());
// 	uint const& start = svar.bndPts;
// 	uint const& end = svar.totPts;

// 	#pragma omp parallel for  reduction(+:cgrad)
// 	for(uint ii=start; ii < end; ++ii)
// 	{
// 		StateVecD const& pi = pnp1[ii].xi;
// 		real bottom = 0.0;
		
// 		for(auto const& jj:outlist[ii])
// 		{	/*Find the denominator to correct absence of second phase*/
// 			SPHPart const& pj = pnp1[jj.first];
// 			real const r = sqrt(jj.second);
// 			bottom += (pj.m/pj.rho)*Kernel(r,fvar.H,fvar.correc);
// 		}

// 		for(auto const& jj:outlist[ii])
// 		{	/*Find the numerator and sum*/
// 			SPHPart const& pj = pnp1[jj.first];
// 			StateVecD const Rji = pj.xi-pi;
// 			real const r = sqrt(jj.second);
// 			cgrad[ii] += (pj.m/pj.rho)*GradK(Rji,r,fvar.H,fvar.correc);
// 		}
		
// 		cgrad[ii] /= bottom;
// 	}

// 	return cgrad;
// }

/*Diffuse interface method from He et al (2014) for surface tension*/
StateVecD HeST(StateVecD const& cgradi, StateVecD const& cgradj, StateVecD const& diffK)
{
	return 0.5*((cgradi.squaredNorm()+cgradj.squaredNorm()))*diffK;
}

#endif

// /*Numerical particle density for Nair & Poeschel (2017) surface tension*/
// real GetNumpartdens(SIM const& svar, FLUID const& fvar, SPHState const& pnp1, OUTL const& outlist)
// {
// 	real npd = 0.0;
// 	uint const& end = svar.totPts;
// 	#pragma omp parallel for reduction(+:npd)
// 	for (size_t ii = 0; ii < end; ++ii)
// 	{
// 		// StateVecD const& pi = pnp1[ii].xi;
// 		for (auto const& jj:outlist[ii])
// 		{ /* Surface Tension calcs */
// 			// StateVecD const& pj = pnp1[jj.first].xi;
// 			real const r = sqrt(jj.second);
// 			// real const r = (pj-pi).norm();
// 			npd += Kernel(r,fvar.H,fvar.correc);
// 		}
// 	}
// 	return npd/real(svar.totPts);
// }

// /*Surface Tension - Nair & Poeschel (2017)*/
#ifdef PAIRWISE
inline real surface_tension_fac(int const bA, int const bB)
{
	if(bA == BOUND || bB == BOUND )
	{
		real contang = 0.5 * M_PI * 7.0/9.0;
		return (1.0 + 0.5*cos(contang));
	}
	return 1.0;
}

inline StateVecD SurfaceTens(StateVecD const& Rji, real const& r, real const& h, real const& sig,  real const& lam, 
					real const& npdm2, real const& pi3o4, int const& bA, int const& bB)
{
	real fac = surface_tension_fac(bA, bB); 

	/*npd = numerical particle density (see code above) */
	return -0.5*npdm2*(sig/lam)*fac*cos(pi3o4*r/h)*(Rji/r);
}
#endif

#ifdef CSF
inline void CSF_Curvature(FLUID const& fvar, DELTAP const& dp, SPHPart const& pi, SPHPart const& pj,
					StateVecD const& Rji, StateVecD const& gradK, real const& volj, real const& r, real& curve, real& correc)
{
	if( dp.lam[pj.partID] < 0.75 )
	{	
		curve -= (pj.norm.normalized()-pi.norm.normalized()).dot(volj*gradK);
		correc += volj * Kernel(r,fvar.H,fvar.correc)/*/dp.kernsum[ii]*/;
	}	

	/* Consider imaginary particles to improve curvature representation */
	if(pi.surf == 1 && pj.surf != 1)
	{	
		/* Host particle is on the surface, but neighbour is not. Reflect stuff then */
		StateVecD normj = 2*pi.norm.normalized() - pj.norm.normalized();
		curve += (normj.normalized()-pi.norm.normalized()).dot(volj*gradK);
		correc += volj * Kernel(r,fvar.H,fvar.correc)/*/dp.kernsum[ii]*/;
	}

	if(pj.surf == 1 && pi.surf != 1)
	{
		/* Neighbour particle is on the surface, but host is not. Not as easy */
		/* Need to check if extended particle is still in the neighbourhood */
		if(r < fvar.H)
		{	/* Assuming a support radius of 2h, if the current particle is within h */
			/* it will still be in the neighbourhood if distance is doubled */
			StateVecD gradK2 = GradK(2*Rji,2*r,fvar.H,fvar.correc);
			StateVecD normj = 2*pj.norm.normalized() - pi.norm.normalized();
			curve -= (normj.normalized()-pi.norm.normalized()).dot(volj*gradK2);
			correc += volj * Kernel(r,fvar.H,fvar.correc)/*/dp.kernsum[ii]*/;
		}
	}
}
#endif

#ifdef ALE
/* Arbitrary Lagrangian Eulerian formulation - Sun, Colagrossi, Marrone, Zhang (2018)*/
inline StateVecD ALEMomentum(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK, real const& rho0)
{
	return 	 ((pj.v * pj.vPert.transpose() + pi.v * pi.vPert.transpose()) * gradK 
			- pi.v * (pj.vPert - pi.vPert).dot(gradK)) * Vj ;
}

inline real ALEContinuity(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK)
{
	return ((pj.v+pj.vPert) - (pi.v + pi.vPert)).dot(gradK)*Vj;
}

inline real ALECont2ndterm(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK)
{
	return (pj.rho*pj.vPert + pi.rho*pi.vPert).dot(gradK)*Vj;
}
#endif

/* delta-SPH dissipation term in the continuity equation*/
#ifndef NODSPH
inline real Continuity_dSPH(StateVecD const& Rji, real const& idist2, real const& HSQ, StateVecD const& Grad, 
				real const& volj, StateVecD const& gRho_i, StateVecD const& gRho_j, SPHPart const& pi, SPHPart const& pj)
{
	return volj * ((pj.rho - pi.rho) - 0.5*(gRho_i + gRho_j).dot(Rji)) * Rji.dot(Grad)*idist2;
}
#endif


inline StateVecD BasePos(SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
	/* Positive Pressure - Monaghan (1994)*/
	return pj.m * gradK * (pj.p/(pj.rho*pj.rho)+pi.p/(pi.rho*pi.rho));
}

inline StateVecD BaseNeg(SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
	/* Negative Pressure - Monaghan (1994)*/
	return pj.m * gradK * (pj.p/(pj.rho*pj.rho)-pi.p/(pi.rho*pi.rho));
}

inline StateVecD pressPos(real const& volj, SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
	/* Positive Pressure - Marrone et al. (2011)*/
	return gradK * (pj.p + pi.p) * volj;
}

inline StateVecD pressNeg(real const& volj, SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
	/* Negative Pressure - Marrone et al. (2011)*/
	return gradK * (pj.p - pi.p) * volj;
}

// #ifdef COLEEOS
inline StateVecD ArtVisc(real const& nu, SPHPart const& pi, SPHPart const& pj, FLUID const& fvar, 
				StateVecD const& Rji, StateVecD const& Vji, real const idist2, StateVecD const& gradK)
{
	// if(pj.b != BOUND)
	// {
		real const vdotr = Vji.dot(Rji);

		if (vdotr > 0.0) 
		{
			return StateVecD::Zero();
		}
		else
		{	
			real const muij= fvar.H*vdotr*idist2;
			real const rhoij = 0.5*(pi.rho+pj.rho);
			real const cbar = 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));

			real const alphv = 2 * nu * (SIMDIM + 2)/(fvar.H*cbar);
			real const alpha = std::max(fvar.alpha,alphv);

			return gradK * alpha *cbar*muij/rhoij;
		}
	// }

	return StateVecD::Zero();
}
// #else
inline StateVecD aVisc(StateVecD const& Rji, StateVecD const& Vji, real const idist2, 
	 StateVecD const& gradK, real const& volj)
{
	real const vdotr = Vji.dot(Rji);

	if (vdotr > 0.0) 
	{
		return StateVecD::Zero();
	}
	
    return vdotr * volj * gradK * idist2;
}
// #endif



/*Laminar Viscosity - Morris (2003)*/
/*Apparently divergent at free surfaces - consider removing.*/
// StateVecD Viscosity(real const& mu, real const& HSQ, SPHPart const& pi, SPHPart const& pj, 
// 	StateVecD const& Rji, StateVecD const& Vji, real const& r, StateVecD const& gradK)
// {
// 	return Vji*(mu/(pi.rho*pj.rho))*(1.0/(r*r+0.01*HSQ))*Rji.dot(gradK);
// }

inline StateVecD Viscosity(real const& nu, SPHPart const& pi, SPHPart const& pj, 
	StateVecD const& Rji, StateVecD const& Vji, real const& idist2, StateVecD const& gradK)

{
	return nu * (pi.rho + pj.rho)/(pi.rho*pj.rho) * (Rji.dot(gradK))*idist2 * Vji;		
}


/*Repulsion for interacting with mesh surface - saves generating particles on surface*/
inline StateVecD NormalBoundaryRepulsion(FLUID const& fvar, MESH const& cells, SPHPart const& pi)
{
    real beta = 4*fvar.Cs*fvar.Cs;
    real kern = BoundaryKernel(pi.y,fvar.H,beta);
	return fvar.bndM/(fvar.bndM+fvar.simM)*kern*pi.bNorm;
}


/* Boundary pressure calculation */
void Get_Boundary_Pressure(FLUID const& fvar, size_t const& end, OUTL const& outlist, SPHState& pnp1)
{
		/*Gravity Vector*/
	#if SIMDIM == 3
		StateVecD const g(0.0,0.0,-9.81);
	#else
		StateVecD const g(0.0,-9.81);
		// StateVecD const g(0.0,0.0);
	#endif

	#pragma omp parallel for schedule(static) default(shared)/*Reduction defs in Var.h*/
	for (size_t ii=0; ii < end; ++ii)
	{
		SPHPart const& pi = pnp1[ii];

		real kernsum = 0.0;
		real pkern = 0.0;
		StateVecD acckern = StateVecD::Zero();
		// StateVecD velkern = StateVecD::Zero();

		/* Get the 'extrapolated velocity' of the boundary  */
		// for (std::pair<size_t,real> const& jj : outlist[ii])
		// {	/* Neighbour list loop. */
		// 	SPHPart const& pj = pnp1[jj.first];
		// 	/*Check if the neighbour is a fluid particle*/
		// 	if(pj.b > PISTON)
		// 	{
		// 		real const rr = jj.second;
		// 		real const r = sqrt(rr);
		// 		real kern = Kernel(r,fvar.H,fvar.correc);
		// 		kernsum += kern;
		// 		velkern += pj.v*kern;
		// 	}
		// }

		// StateVecD const ext_vel = velkern/kernsum;

		for (std::pair<size_t,real> const& jj : outlist[ii])
		{	/* Neighbour list loop. */
			SPHPart const& pj = pnp1[jj.first];

			/*Check if the neighbour is a fluid particle*/
			if(pj.b > PISTON)
			{
				StateVecD const Rji = pj.xi-pi.xi;
				
				real const rr = jj.second;
				real const r = sqrt(rr);
				real kern = Kernel(r,fvar.H,fvar.correc);
				// StateVecD const Vji = pj.v-ext_vel;
				// StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);/* gradK;*/
				// StateVecD const aVisc = ArtVisc(fvar.nu,pi,pj,fvar,Rji,Vji,rr,gradK);
				kernsum += kern;
				pkern += pj.p * kern;
				acckern +=  kern * pj.rho * Rji;
			}
		}/*End of neighbours*/

		real pressure = 0.0;
		real density = fvar.rho0;
		if(kernsum > 0.0)
		{
			pressure = (pkern-g.dot(acckern))/kernsum;
			density = density_equation(pressure,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);
		}

		pnp1[ii].p = pressure;
		pnp1[ii].rho = density;
		
	} /*End of boundary parts*/
}

///**************** RESID calculation **************
void Forces(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, SPHState const& pnp1,
	 OUTL const& outlist, DELTAP const& dp,
	 vector<StateVecD>& RV, vector<real>& Rrho, std::vector<StateVecD>& Af, StateVecD& Force)
{

	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;

	#ifdef PAIRWISE
	real const npdm2 = 1/(dp.npd*dp.npd);
	real const pi3o4 = 3.0*M_PI/4.0;

	/*Surface tension factor*/
	// #if SIMDIM==2
	// real const lam = 8.0/81.0*pow(fvar.H,4)/pow(M_PI,4)*(9.0/4.0*pow(M_PI,3)-6.0*M_PI-4.0);
	real const lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));
	// #else
	// real const lam = 3.0 / (4.0 * pow(M_PI, 4.0)) * (pow(2.0, 7) - 9.0 * pow(2.0, 4) 
	// 				* M_PI * M_PI + 9.0 * 3.0 * pow(M_PI, 4)) * pow(fvar.H / 3.0, 5);
	// #endif
	
	#endif

	/********* LOOP 0 - all points: Calculate numpartdens ************/
	// wDiff = GetWeightedDistance(svar,fvar,pnp1,outlist);
	// #ifdef HESF
	// std::vector<StateVecD> const cgrad = GetColourGrad(svar,fvar,pnp1,outlist);
	// #else
	// #ifndef CSF
	// #ifndef HESF
	// real const npd = GetNumpartdens(svar, fvar, pnp1, outlist);
	// #endif
	// #endif
	// 
	// std::vector<StateVecD> ST(svar.totPts,StateVecD::Zero()); /*Surface tension force*/		

	Rrho = vector<real>(end,0.0);
	RV = vector<StateVecD>(end,StateVecD::Zero());
	Af = vector<StateVecD>(end,StateVecD::Zero());

	#pragma omp parallel default(shared)
	{
		/*Gravity Vector*/
		StateVecD const& g = svar.grav;
		

		if(svar.bound_solver == 1)
		{
			/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/		
			#pragma omp for schedule(static) nowait /*Reduction defs in Var.h*/
			for (size_t ii=0; ii < start; ++ii)
			{
				SPHPart const& pi = pnp1[ii];
				// pnp1[ii].theta = real(size);
				// pnp1[ii].theta = wDiff[ii];
				real Rrhoi = 0.0;

				for (std::pair<size_t,real> const& jj : outlist[ii])
				{	/* Neighbour list loop. */
					SPHPart const& pj = pnp1[jj.first];
					if(ii == jj.first)
						continue;

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

/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Force,RV,Af) schedule(static) nowait
		for (size_t ii = start; ii < end; ++ii)
		{
			SPHPart const& pi = pnp1[ii];
			// size_t const size = outlist[ii].size();
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];

			real Rrhoi = 0.0;
			#ifdef ALE
			StateVecD RVia = StateVecD::Zero();
			real Rrhoc = 0.0;
			#endif
			
			#ifdef CSF
			real curve = 0.0;
			real correc = 0.0;
			#endif
			
			StateVecD Vdiff = StateVecD::Zero();
			real Pbasei = 0.0;

			StateVecD aero = StateVecD::Zero();
			StateVecD RVi = StateVecD::Zero();
			
			StateVecD viscI = StateVecD::Zero();
			StateVecD surfT = StateVecD::Zero();
			
			#ifdef NOFROZEN
			StateVecD artViscI = StateVecD::Zero();
			#ifndef NODSPH
			real Rrhod = 0.0;			
			#endif
			#endif

			if( dp.lam_ng[ii] < avar.cutoff /* pi.surf == 1 */ && pi.b == FREE )
			{
				Vdiff = pi.cellV - pi.v /* dp.avgV[ii] */;
				
				if (svar.Asource == 2)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - pi.v /* dp.avgV[ii] */;
				}
				
				aero = CalcAeroAcc(avar,pi,Vdiff,dp.norm[ii],dp.lam_ng[ii],Pbasei);
				RV[ii] += aero;
				Af[ii] += aero;
				Force += aero;
			}

			for (std::pair<size_t,real> const& jj : outlist[ii])
			{	/* Neighbour list loop. */
				SPHPart const& pj = pnp1[jj.first];
				/*Check if the position is the same, and skip the particle if yes*/
				if(ii == jj.first)
				{
					// if(/* pi.surf == 1 */ dp.lam_ng[ii] < 0.75 && pi.b == FREE)
					// {
					// 	RV[ii] += aero * fvar.correc/dp.kernsum[ii];
					// 	Af[ii] += aero * fvar.correc/dp.kernsum[ii];
					// 	Force +=  aero * fvar.correc/dp.kernsum[ii];
					// }	
					continue;
				}

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
				
				// if( dp.lam_ng[ii] < 0.75 && pi.b == FREE && 
				// 	pj.b > PISTON && pj.b != GHOST)
				// {	/* Assumed that particle masses are the same */
				// 	RV[jj.first] += aero * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				// 	Af[jj.first] += aero * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				// 	Force += aero * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				// }	

				// StateVecD const	contrib = BasePos(pi,pj,gradK);
				StateVecD const contrib = pressPos(volj,pi,pj,gradK);

				#ifdef ALE
					StateVecD const ALEcontrib = ALEMomentum(pi, pj, volj, gradK, fvar.rho0);
				#endif

				/*Laminar Viscosity - Morris (2003)*/
				StateVecD const visc = Viscosity(fvar.nu,pi,pj,Rji,Vji,idist2,gradK);
				viscI += pj.m*visc;

				/* Surface tension options */
				if(/* pi.b == FREE  && */ pj.b != GHOST)
				{
					/* He et al (2014) - Diffuse Interface model */
					#ifdef HESF
					if( dp.lam[ii] < 0.75 /* pi.surf == 1 */ )
					{
						surfT += HeST(pi.norm/* /dp.colour[ii] */,pj.norm/* /dp.colour[jj.first] */,volj*gradK);
					}
					#endif

					#ifdef CSF
					/* Ordoubadi et al. (2017) - Continuum Surface Force + reflected particles */
					if( dp.lam[ii] < 0.75 /* pi.surf == 1 */ )
					{
						CSF_Curvature(fvar,dp,pi,pj,Rji,gradK,volj,r,curve,correc);
					}
					#endif
					
					#ifdef PAIRWISE
					#ifdef ALE
					if(pi.surfzone == 1)	/* Nair & Poeschel (2017) - Pairwise force */
					#endif
						surfT += SurfaceTens(Rji,r, fvar.H, fvar.sig,lam,npdm2,pi3o4,pi.b,pj.b);
					#endif	
				}
				
				
				#ifdef ALE
					RVi -= contrib;
					RVia += ALEcontrib;
				#else
					RVi -= contrib;
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
					if(pj.b != BOUND)
					{
						#ifndef NODSPH
							Rrhod += Continuity_dSPH(Rji,idist2,fvar.HSQ,gradK,volj,
												dp.gradRho[ii],dp.gradRho[jj.first],pi,pj);
						#endif

						artViscI += pj.m*ArtVisc(fvar.nu,pi,pj,fvar,Rji,Vji,idist2,gradK);
					}
				#endif
				
				

			}/*End of neighbours*/

			if(pi.internal == 1)
			{	// Apply the normal boundary force
				RVi += NormalBoundaryRepulsion(fvar, cells, pi);
			}

			#ifdef CSF
			if(dp.lam[ii] < 0.75 /* pi.surf == 1 */)
			{
				RVi += (fvar.sig/pi.rho * curve * pi.norm)/correc;
			}
			#endif

			#ifdef HESF
			surfT *= (0.02/2.0)*(pi.m/pi.rho);
			#endif
			
			// #pragma omp critical
			// Force += aero;	
			
			#ifdef NOFROZEN
				#ifdef ALE
				RV[ii] += (RVi/pi.rho + RVia + artViscI + viscI + surfT/pi.m /* + aero/pi.m */ + g);
				#else
				RV[ii] += (RVi/pi.rho + artViscI + viscI + surfT/pi.m /* + aero/pi.m */ + g);
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
				RV[ii] += (RVi/pi.rho + RVia + pnp1[ii].aVisc + viscI + surfT/pi.m /* + aero/pi.m */ + g);
				#else
				RV[ii] += (RVi/pi.rho + pnp1[ii].aVisc + viscI + surfT/pi.m /* + aero/pi.m */ + g);
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

#endif