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

inline real pressure_equation(real const& rho, real const& B, real const& gam, real const& Cs, real const& rho0)
{
	#ifdef COLEEOS
		(void)Cs;
		return B*(pow(rho/rho0,gam)-1); /* Cole EOS */
	#else
		(void)B;
		(void)gam;
		return Cs*Cs * (rho - rho0); /* Isothermal EOS */
	#endif
}



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
// 			StateVecD const Rij = pj.xi-pi;
// 			real const r = sqrt(jj.second);
// 			cgrad[ii] += (pj.m/pj.rho)*GradK(Rij,r,fvar.H,fvar.correc);
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

#else

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
inline real surface_tension_fac(int const bA, int const bB)
{
	if(bA == BOUND || bB == BOUND )
	{
		real contang = 0.5 * M_PI * 7.0/9.0;
		return (1.0 + 0.5*cos(contang));
	}
	return 1.0;
}

inline StateVecD SurfaceTens(StateVecD const& Rij, real const& r, real const& h, real const& sig,  real const& lam, 
					real const& npdm2, real const& pi3o4, int const& bA, int const& bB)
{
	/*Surface tension factor*/
	// real const lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
	// 						(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));

	real fac = surface_tension_fac(bA, bB); 

	/*npd = numerical particle density (see code above) */
	return -0.5*npdm2*(sig/lam)*fac*cos(pi3o4*r/h)*(Rij/r);
}

#endif

#ifdef ALE
/* Arbitrary Lagrangian Eulerian formulation - Sun, Colagrossi, Marrone, Zhang (2018)*/
StateVecD ALEMomentum(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK, real const& rho0)
{
	return 	(rho0/pi.rho) * (pj.v * pj.vPert.transpose() + pi.v * pi.vPert.transpose()) * gradK * Vj 
			- pi.v * (pj.vPert - pi.vPert).dot(gradK) * Vj ;
}

real ALEContinuity(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK)
{
	return pi.rho*((pj.v+pj.vPert) - (pi.v + pi.vPert)).dot(gradK)*Vj/*  - 
			(pj.rho*pj.vPert + pi.rho*pi.vPert).dot(gradK)*Vj */ ;
}
#endif

/* delta-SPH dissipation term in the continuity equation*/
#ifndef NODSPH
real Continuity_dSPH(StateVecD const& Rij, real const& rr, real const& HSQ, StateVecD const& Grad, 
				real const& mj, StateVecD const& gRho_i, StateVecD const& gRho_j, SPHPart const& pi, SPHPart const& pj)
{
	return mj * ((pj.rho - pi.rho) + 0.5*(gRho_i + gRho_j).dot(Rij)) * Rij.dot(Grad)/(rr+0.001*HSQ);
}
#endif


// StateVecD ArtVisc(StateVecD const& Rij, StateVecD const& Vij, real const r, 
// 	 StateVecD const& gradK, real const& Vj)
// {
//     return Vij.dot(Rij) * Vj * gradK / (r*r);
// }


// StateVecD BasePos(real const& Pi, real const& Pj, real const& Vj, StateVecD const& gradK)
// {
// 	/* Positive Pressure - Sun, Colagrossi, Marrone, Zhang (2016)*/
// 	return gradK * Vj * (Pj + Pi);
// }

// StateVecD BaseNeg(real const& Pi, real const& Pj, real const& Vj, StateVecD const& gradK)
// {
// 	/* Negative Pressure - Sun, Colagrossi, Marrone, Zhang (2016)*/
// 	return gradK * Vj * (Pj - Pi);
// }

// StateVecD ArtVisc(StateVecD const& Rij, StateVecD const& Vij, real const r, 
// 	 StateVecD const& gradK, real const& Vj)
// {
//     return Vij.dot(Rij) * Vj * gradK / (r*r);
// }

StateVecD BasePos(SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
	/* Positive Pressure - Monaghan (1994)*/
	return gradK * (pj.p*pow(pj.rho,-2)+pi.p*pow(pi.rho,-2));
}

StateVecD BaseNeg(SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
	/* Negative Pressure - Monaghan (1994)*/
	return gradK * (pj.p*pow(pj.rho,-2)-pi.p*pow(pi.rho,-2));
}

StateVecD ArtVisc(real const& nu, SPHPart const& pi, SPHPart const& pj, FLUID const& fvar, StateVecD const& Rij, StateVecD const& Vij, real const rr, 
	 StateVecD const& gradK)
{
	if(pj.b != BOUND)
	{
		real const vdotr = Vij.dot(Rij);

		if (vdotr > 0.0) 
		{
			return StateVecD::Zero();
		}
		else
		{	
			real const muij= fvar.H*vdotr/(rr+0.0001*fvar.HSQ);
			real const rhoij = 0.5*(pi.rho+pj.rho);
			real const cbar = 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));

			real const alphv = 2 * nu * (SIMDIM + 2)/(fvar.H*cbar);
			real const alpha = std::max(fvar.alpha,alphv);

			return gradK * alpha *cbar*muij/rhoij;
		}
	}
	else
		return StateVecD::Zero();
}

/*Laminar Viscosity - Morris (2003)*/
/*Apparently divergent at free surfaces - consider removing.*/
// StateVecD Viscosity(real const& mu, real const& HSQ, SPHPart const& pi, SPHPart const& pj, 
// 	StateVecD const& Rij, StateVecD const& Vij, real const& r, StateVecD const& gradK)
// {
// 	return Vij*(mu/(pi.rho*pj.rho))*(1.0/(r*r+0.01*HSQ))*Rij.dot(gradK);
// }

StateVecD Viscosity(real const& nu, real const& HSQ, SPHPart const& pi, SPHPart const& pj, 
	StateVecD const& Rij, StateVecD const& Vij, real const& rr, StateVecD const& gradK)

{
	return nu * (pi.rho + pj.rho)/(pi.rho*pj.rho) * (Rij.dot(gradK))/(rr + 0.001*HSQ) * Vij;		
}


/*Repulsion for interacting with mesh surface - saves generating particles on surface*/
StateVecD NormalBoundaryRepulsion(FLUID const& fvar, MESH const& cells, SPHPart const& pi)
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
	#pragma omp parallel for schedule(static) shared(pnp1)/*Reduction defs in Var.h*/
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
				StateVecD const Rij = pj.xi-pi.xi;
				
				real const rr = jj.second;
				real const r = sqrt(rr);
				real kern = Kernel(r,fvar.H,fvar.correc);
				// StateVecD const Vij = pj.v-ext_vel;
				// StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);/* gradK;*/
				// StateVecD const aVisc = ArtVisc(fvar.nu,pi,pj,fvar,Rij,Vij,rr,gradK);
				kernsum += kern;
				pkern += pj.p * kern;
				acckern +=  kern * pj.rho * Rij;
			}
		}/*End of neighbours*/

		real pressure = 0.0;
		real density = fvar.rho0;
		if(kernsum > 0.0)
		{
			pressure = (pkern-g.dot(acckern))/kernsum;
			density = fvar.rho0*pow((pressure/fvar.B) + 1.0, 1.0/fvar.gam);
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

	real const npdm2 = 1/(dp.npd*dp.npd);
	real const pi3o4 = 3.0*M_PI/4.0;

	#if SIMDIM==2
	real const lam = 8.0/81.0*pow(fvar.H,4)/pow(M_PI,4)*(9.0/4.0*pow(M_PI,3)-6.0*M_PI-4.0);
	#else
	real const lam = 3.0 / (4.0 * pow(M_PI, 4.0)) * (pow(2.0, 7) - 9.0 * pow(2.0, 4) 
					* M_PI * M_PI + 9.0 * 3.0 * pow(M_PI, 4)) * pow(2.0 * fvar.H / 3.0, 5);
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

	#pragma omp parallel default(shared) // shared(RV,Rrho,Af)
	{
		/*Gravity Vector*/
		StateVecD const& g = svar.grav;
		

/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/		
		// #pragma omp for schedule(static) nowait /*Reduction defs in Var.h*/
		// for (size_t ii=0; ii < start; ++ii)
		// {
		// 	SPHPart const& pi = pnp1[ii];
		// 	// pnp1[ii].theta = real(size);
		// 	// pnp1[ii].theta = wDiff[ii];
		// 	real Rrhoi = 0.0;
		// 	real Rrhod = 0.0;

		// 	real kernsum = 0.0;
		// 	real pkern = 0.0;
		// 	real acckern = 0.0;

		// 	for (std::pair<size_t,real> const& jj : outlist[ii])
		// 	{	/* Neighbour list loop. */
		// 		SPHPart const& pj = pnp1[jj.first];

		// 		/*Check if the position is the same, and skip the particle if yes*/
		// 		if(pj.b > PISTON)
		// 		{
		// 			StateVecD const Rij = pj.xi-pi.xi;
		// 			StateVecD const Vij = pj.v-pi.v;
		// 			real const rr = jj.second;
		// 			real const r = sqrt(rr);
		// 			real const volj = pj.m/pj.rho;
		// 			real kern = Kernel(r,fvar.H,fvar.correc);
		// 			kernsum += kern;
		// 			pkern += pj.p * kern;
		// 			acckern += pj.rho * kern * (g  - pj.f).dot(Rij);
		// 		}

		// 		// StateVecD const Rij = pj.xi-pi.xi;
		// 		// StateVecD const Vij = pj.v-pi.v;
		// 		// real const rr = jj.second;
		// 		// real const r = sqrt(rr);
		// 		// real const volj = pj.m/pj.rho;
		// 		// StateVecD const gradK = /*dp.L[ii] * */GradK(Rij, r,fvar.H, fvar.correc);

		// 		// // StateVecD contrib = BasePos(pi.p,pj.p,volj,gradK);
		// 		// // StateVecD contrib = BasePos(pi,pj,gradK);

		// 		// /*drho/dt*/
		// 		// Rrhoi -= pj.m*(Vij.dot(gradK));
		// 		// Rrhod -= Continuity_dSPH(Rij,rr,fvar.HSQ,gradK,volj,dp.gradRho[ii],dp.gradRho[pj.partID],pi,pj);
		// 		// RV[ii] -= pj.m*contrib;
		// 	}/*End of neighbours*/
			
			 

		// 	RV[ii] = StateVecD::Zero();
		// 	Af[ii] = StateVecD::Zero();
		// 	Rrho[ii] = (Rrhoi - fvar.dCont * Rrhod);

		// } /*End of boundary parts*/

		// vector<StateVecD> gradRho(end,StateVecD::Zero());

		// #pragma omp for schedule(static) nowait
		// for (size_t ii = start; ii < end; ++ii)
		// {
		// 	SPHPart const& pi = pnp1[ii];

		// 	StateVecD gradRho_ = StateVecD::Zero();

		// 	for(std::pair<size_t,real> const& jj:outlist[ii])
		// 	{
		// 		SPHPart const& pj = pnp1[jj.first];
		// 		/*Check if the position is the same, and skip the particle if yes*/
		// 		if(ii == jj.first)
		// 			continue;
				

		// 		StateVecD const Rij = pj.xi - pi.xi;
		// 		real const r = sqrt(jj.second);
		// 		real const volj = pj.m/pj.rho;
		// 		StateVecD const Grad = GradK(-Rij,r,fvar.H,fvar.correc);

		// 		gradRho_ += volj * (pj.rho - pi.rho)*Grad;          
		// 	}

		// 	gradRho[ii] = dp.L[ii]*gradRho_;
		// }


/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Force,RV,Af) schedule(static) nowait
		for (size_t ii = start; ii < end; ++ii)
		{
			SPHPart const& pi = pnp1[ii];
			// size_t const size = outlist[ii].size();
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];

			real Rrhoi = 0.0;
			
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

			if( dp.lam_ng[ii] < 0.75 /* pi.surf == 1 */ && pi.b == FREE )
			{
				Vdiff = pi.cellV - pi.v /* dp.avgV[ii] */;
				
				if (svar.Asource == 2)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - pi.v /* dp.avgV[ii] */;
				}
				
				aero = CalcAeroAcc(avar,pi,Vdiff,dp.norm[ii],dp.lam_ng[ii],Pbasei);
			}

			for (std::pair<size_t,real> const& jj : outlist[ii])
			{	/* Neighbour list loop. */
				SPHPart const& pj = pnp1[jj.first];
				/*Check if the position is the same, and skip the particle if yes*/
				if(ii == jj.first)
				{
					if(/* pi.surf == 1 */ dp.lam_ng[ii] < 0.75 && pi.b == FREE)
					{
						RV[ii] += aero * fvar.correc/dp.kernsum[ii];
						Af[ii] += aero * fvar.correc/dp.kernsum[ii];
						Force +=  aero * fvar.correc/dp.kernsum[ii];
					}	
					continue;
				}

				StateVecD const Rij = pj.xi-pi.xi;
				StateVecD const Vij = pj.v-pi.v;
				real const rr = jj.second;
				real const r = sqrt(rr);
				
				#if defined(CSF) || defined(HESF)
					real const volj = pj.m/pj.rho;
				#endif

				// StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);/* gradK;*/

				// if(gradK == StateVecD::Zero())
				// {
				// 	cout << ii << "  " << jj.first << endl;
				// 	cout << rr << "  " << r << "  " << r/fvar.H << endl;
 				// }

				
				
				if( dp.lam_ng[ii] < 0.75 /* pi.surf == 1 */ && pi.b == FREE && 
					pj.b > PISTON && pj.b != GHOST)
				{	/* Assumed that particle masses are the same */
					RV[jj.first] += aero * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
					Af[jj.first] += aero * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
					Force += aero * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				}	

				StateVecD const	contrib = BasePos(pi,pj,gradK);

				#ifdef ALE
					StateVecD const ALEcontrib = ALEMomentum(pi, pj, volj, gradK, fvar.rho0);
				#endif

				/*Laminar Viscosity - Morris (2003)*/
				StateVecD const visc = Viscosity(fvar.nu,fvar.HSQ,pi,pj,Rij,Vij,r,gradK);

				/*Base WCSPH continuity drho/dt*/
				#ifdef ALE
					Rrhoi -= ALEContinuity(pi,pj,volj,gradK);
				#else
					Rrhoi -= pj.m*(Vij.dot(gradK));
				#endif


				if(pi.b == FREE  && pj.b != GHOST)
				{
					/* He et al (2014) - Diffuse Interface model */
					#ifdef HESF
					if( dp.lam[ii] < 0.75 /* pi.surf == 1 */ )
					{
						surfT += HeST(pi.norm/* /dp.colour[ii] */,pj.norm/* /dp.colour[jj.first] */,volj*gradK);
					}
					#else

					#ifdef CSF
					/* Ordoubadi et al. (2017) - Continuum Surface Force + reflected particles */
					if( dp.lam[ii] < 0.75 /* pi.surf == 1 */ )
					{
						if( dp.lam[pj.partID] < 0.75 /* pj.surf == 1 */ )
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
								StateVecD gradK2 = GradK(2*Rij,2*r,fvar.H,fvar.correc);
								StateVecD normj = 2*pj.norm.normalized() - pi.norm.normalized();
								curve -= (normj.normalized()-pi.norm.normalized()).dot(volj*gradK2);
								correc += volj * Kernel(r,fvar.H,fvar.correc)/*/dp.kernsum[ii]*/;
							}
						}
					}
					#else
					/* Nair & Poeschel (2017) - Pairwise force */
					surfT += SurfaceTens(Rij,r, fvar.H, fvar.sig,lam,npdm2,pi3o4,pi.b,pj.b);
					#endif
					#endif	
				}
				
				
				#ifdef ALE
					RVi -= pj.m*contrib  - ALEcontrib;
				#else
					RVi -= pj.m*contrib;
				#endif

				#ifdef NOFROZEN
				StateVecD const aVisc = ArtVisc(fvar.nu,pi,pj,fvar,Rij,Vij,rr,gradK);
				#ifndef NODSPH
				
				Rrhod -= Continuity_dSPH(Rij,rr,fvar.HSQ,gradK,pj.m,dp.gradRho[ii],dp.gradRho[jj.first],pi,pj);
				#endif

				artViscI += pj.m*aVisc;
				#endif

				viscI += pj.m*visc;

			}/*End of neighbours*/

			if(pi.internal == 1)
			{	// Apply the normal boundary force
				RVi += NormalBoundaryRepulsion(fvar, cells, pi);
			}

			#ifdef CSF
			if(pi.b == FREE && dp.lam[ii] < 0.75 /* pi.surf == 1 */)
			{
				RVi += (fvar.sig/pi.rho * curve * pi.norm)/correc;
			}
			#else
			#ifdef HESF
			surfT *= (0.02/2.0)*(pi.m/pi.rho);
			#endif
			#endif
			
			// #pragma omp critical
			// Force += aero;	
			
			#ifdef NOFROZEN
			if( pi.b == FREE )
				RV[ii] += (RVi + artViscI + viscI + surfT/pi.m /* + aero/pi.m */ + g);
			else
				RV[ii] += (RVi + artViscI + viscI /* + surfT/pi.m */ /* + aero/pi.m */);
			
			#ifndef NODSPH
			Rrho[ii] = Rrhoi - fvar.dCont * Rrhod;
			#else
			Rrho[ii] = Rrhoi;
			#endif
			#else
			if( pi.b == FREE )
				RV[ii] += (RVi + pnp1[ii].aVisc + viscI + surfT/pi.m /* + aero/pi.m */ + g);
			else
				RV[ii] += (RVi + pnp1[ii].aVisc + viscI /* + surfT/pi.m */ /* + aero/pi.m */);
			
			#ifndef NODSPH
			Rrho[ii] = Rrhoi - pnp1[ii].deltaD;
			#else
			Rrho[ii] = Rrhoi;
			#endif

			#endif
			// res_.emplace_back(aero/pi.m);
			// Rrho_.emplace_back(0.0);

		} /*End of sim parts*/		

	}	/*End of declare parallel */
 
    // cout << RV.size() << "  " << Rrho.size() << "  " << Af.size() << "  " << curv.size() << endl;
}

#endif