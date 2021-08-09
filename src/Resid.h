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
#include "Add.h"

// /*Numerical particle density for Nair & Poeschel (2017) surface tension*/
// real GetNumpartdens(SIM const& svar, FLUID const& fvar, State const& pnp1, outl const& outlist)
// {
// 	real npd = 0.0;
// 	uint const& end = svar.totPts;
// 	#pragma omp parallel for reduction(+:npd)
// 	for (uint ii=0; ii< end; ++ii)
// 	{
// 		StateVecD const& pi = pnp1[ii].xi;
// 		for (auto jj:outlist[ii])
// 		{ /* Surface Tension calcs */
// 			StateVecD const& pj = pnp1[jj].xi;
// 			real const r = (pj-pi).norm();
// 			npd += W2Kernel(r,fvar.H,fvar.correc);
// 		}
// 	}
// 	return npd/real(svar.totPts);
// }

// /*Surface Tension - Nair & Poeschel (2017)*/
// StateVecD SurfaceTens(FLUID const& fvar, Part const& pj, StateVecD const& Rij, 
// 					  real const& r, real const& npd)
// {
// 	/*Surface tension factor*/
// 	real const lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
// 							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));

// 	real fac=1.0; /*Boundary Correction Factor*/
// 	if(pj.b==PartState.BOUND_) fac=(1+0.5*cos(M_PI*(fvar.contangb/180))); 

// 	/*npd = numerical particle density (see code above) */
// 	real const sij = 0.5*pow(npd,-2.0)*(fvar.sig/lam)*fac;
// 	return -(Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*fvar.H));
// }

// /* Colour field gradient in He et al (2014) method, note the bottom variable to account for no air*/
// std::vector<StateVecD> GetColourGrad(SIM const& svar, FLUID const& fvar, State const& pnp1, outl const& outlist)
// {
// 	std::vector<StateVecD> cgrad(svar.totPts, StateVecD::Zero());
// 	const uint& start = svar.bndPts;
// 	const uint& end = svar.totPts;

// 	#pragma omp parallel for shared(outlist) reduction(+:cgrad)
// 	for(uint ii=start; ii < end; ++ii)
// 	{
// 		StateVecD const& pi = pnp1[ii].xi;
// 		real bottom = 0.0;
		
// 		for(auto jj:outlist[ii])
// 		{	/*Find the denominator to correct absence of second phase*/
// 			Part const& pj = pnp1[jj];
// 			real const r = (pj.xi-pi).norm();
// 			bottom +=(pj.m/pj.rho)*W2Kernel(r,fvar.H,fvar.correc);
// 		}

// 		for(auto jj:outlist[ii])
// 		{	/*Find the numerator and sum*/
// 			Part const& pj = pnp1[jj];
// 			StateVecD const Rij = pj.xi-pi;
// 			real const r = Rij.norm();
// 			cgrad[ii] += (pj.m/pj.rho)*W2GradK(Rij,r,fvar.H,fvar.correc);
// 		}
		
// 		cgrad[ii] = cgrad[ii]/bottom;
// 	}

// 	return cgrad;
// }

// /*Diffuse interface method from He et al (2014) for surface tension*/
// StateVecD HeST(FLUID const& fvar, Part const& pi, Part const& pj, 
// 			  StateVecD const& Rij, real const& r, StateVecD const& cgradi, StateVecD const& cgradj)
// {
// 	return (0.01/2.0)*(pi.m/pi.rho)*(pj.m/pj.rho)*((cgradi.squaredNorm()+cgradj.squaredNorm())/2.0)
// 		*W2GradK(Rij,r,fvar.H,fvar.correc);
// }

#ifndef NOALE
/* Arbitrary Lagrangian Eulerian formulation - Sun, Colagrossi, Marrone, Zhang (2018)*/
StateVecD ALEMomentum(Particle const& pi, Particle const& pj, real const& Vj, StateVecD const& gradK, real const& rho0)
{
	return 	(rho0/pi.rho) * (pj.v * pj.vPert.transpose() + pi.v * pi.vPert.transpose()) * gradK * Vj ;
			// - pi.v * (pj.vPert - pi.vPert).dot(gradK) * Vj ;
}

real ALEContinuity(Particle const& pi, Particle const& pj, real const& Vj, StateVecD const& gradK)
{
	return pi.rho*((pj.v+pj.vPert) - (pi.v + pi.vPert)).dot(gradK)*Vj - 
			(pj.rho*pj.vPert + pi.rho*pi.vPert).dot(gradK)*Vj ;
}
#endif

/* delta-SPH dissipation term in the continuity equation*/
real Continuity_dSPH(StateVecD const& Rij, real const& rr, real const& HSQ, StateVecD const& Grad, 
				real const& volj, StateVecD const& gRho_i, StateVecD const& gRho_j, Particle const& pi, Particle const& pj)
{
	real const psi_ij = ((pj.rho - pi.rho) + 0.5*(gRho_i + gRho_j).dot(Rij));
	return volj * psi_ij * Rij.dot(Grad)/(rr+0.001*HSQ);
}


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

StateVecD BasePos(Particle const& pi, Particle const& pj, StateVecD const& gradK)
{
	/* Positive Pressure - Monaghan (1994)*/
	return gradK * (pj.p*pow(pj.rho,-2)+pi.p*pow(pi.rho,-2));
}

StateVecD BaseNeg(Particle const& pi, Particle const& pj, StateVecD const& gradK)
{
	/* Negative Pressure - Monaghan (1994)*/
	return gradK * (pj.p*pow(pj.rho,-2)-pi.p*pow(pi.rho,-2));
}

StateVecD ArtVisc(real const& nu, Particle const& pi, Particle const& pj, FLUID const& fvar, StateVecD const& Rij, StateVecD const& Vij, real const rr, 
	 StateVecD const& gradK)
{
	if(pj.b != PartState.BOUND_)
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
// StateVecD Viscosity(real const& mu, real const& HSQ, Particle const& pi, Particle const& pj, 
// 	StateVecD const& Rij, StateVecD const& Vij, real const& r, StateVecD const& gradK)
// {
// 	return Vij*(mu/(pi.rho*pj.rho))*(1.0/(r*r+0.01*HSQ))*Rij.dot(gradK);
// }

StateVecD Viscosity(real const& nu, real const& HSQ, Particle const& pi, Particle const& pj, 
	StateVecD const& Rij, StateVecD const& Vij, real const& rr, StateVecD const& gradK)

{
	return nu * (pi.rho + pj.rho)/(pi.rho*pj.rho) * (Rij.dot(gradK))/(rr + 0.001*HSQ) * Vij;		
}


/*Repulsion for interacting with mesh surface - saves generating particles on surface*/
StateVecD NormalBoundaryRepulsion(FLUID const& fvar, MESH const& cells, Particle const& pi)
{
    real beta = 4*fvar.Cs*fvar.Cs;
    real kern = BoundaryKernel(pi.y,fvar.H,beta);
	return fvar.bndM/(fvar.bndM+fvar.simM)*kern*pi.bNorm;
}


/* Boundary pressure calculation */
void Get_Boundary_Pressure(FLUID const& fvar, size_t const& end, outl const& outlist, State& pnp1)
{
		/*Gravity Vector*/
	#if SIMDIM == 3
		StateVecD const g(0.0,-9.81,0.0);
	#else
		StateVecD const g(0.0,-9.81);
		// StateVecD const g(0.0,0.0);
	#endif
	#pragma omp parallel for schedule(static) shared(pnp1)/*Reduction defs in Var.h*/
	for (size_t ii=0; ii < end; ++ii)
	{
		Particle const& pi = pnp1[ii];

		real kernsum = 0.0;
		real pkern = 0.0;
		real acckern = 0.0;
		// StateVecD velkern = StateVecD::Zero();

		/* Get the 'extrapolated velocity' of the boundary  */
		// for (std::pair<size_t,real> const& jj : outlist[ii])
		// {	/* Neighbour list loop. */
		// 	Particle const& pj = pnp1[jj.first];
		// 	/*Check if the neighbour is a fluid particle*/
		// 	if(pj.b > PartState.PISTON_)
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
			Particle const& pj = pnp1[jj.first];

			/*Check if the neighbour is a fluid particle*/
			if(pj.b > PartState.PISTON_)
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
				acckern +=  kern * (pj.rho * g /* + aVisc */).dot(Rij);
			}
		}/*End of neighbours*/

		real pressure = 0.0;
		real density = fvar.rho0;
		if(kernsum > 0.0)
		{
			pressure = (pkern-acckern)/kernsum;
			density = fvar.rho0*pow((pressure/fvar.B) + 1.0, 1.0/fvar.gam);
		}

		pnp1[ii].p = pressure;
		pnp1[ii].rho = density;
		
	} /*End of boundary parts*/
}

///**************** RESID calculation **************
void Forces(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, State const& pnp1,
	 outl const& outlist, DELTAP const& dp,
	 vector<StateVecD>& RV, vector<real>& Rrho, std::vector<StateVecD>& Af, StateVecD& Force)
{
	svar.maxmu=0; 					    /* CFL Parameter */
	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;

	// #pragma omp critical
	// cout << endl << endl;

	
	// const uint piston = svar.psnPts;

	/********* LOOP 0 - all points: Calculate numpartdens ************/
	// wDiff = GetWeightedDistance(svar,fvar,pnp1,outlist);
	// real const npd = GetNumpartdens(svar, fvar, pnp1, outlist);
	// std::vector<StateVecD> const cgrad = GetColourGrad(svar,fvar,pnp1,outlist);
	// std::vector<StateVecD> ST(svar.totPts,StateVecD::Zero()); /*Surface tension force*/		

	

	// for(size_t ii = 0; ii < end; ++ii)
	// {
	// 	cout << ii << "  " << kernsum[ii] << endl;
	// }

	Rrho = vector<real>(end,0.0);
	RV = vector<StateVecD>(end,StateVecD::Zero());
	Af = vector<StateVecD>(end,StateVecD::Zero());

	#pragma omp parallel shared(RV,Rrho,Af)
	{
		/*Gravity Vector*/
		StateVecD const& g = svar.grav;
		

/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/		
		// #pragma omp for schedule(static) nowait /*Reduction defs in Var.h*/
		// for (size_t ii=0; ii < start; ++ii)
		// {
		// 	Particle const& pi = pnp1[ii];
		// 	// pnp1[ii].theta = real(size);
		// 	// pnp1[ii].theta = wDiff[ii];
		// 	real Rrhoi = 0.0;
		// 	real Rrhod = 0.0;

		// 	real kernsum = 0.0;
		// 	real pkern = 0.0;
		// 	real acckern = 0.0;

		// 	for (std::pair<size_t,real> const& jj : outlist[ii])
		// 	{	/* Neighbour list loop. */
		// 		Particle const& pj = pnp1[jj.first];

		// 		/*Check if the position is the same, and skip the particle if yes*/
		// 		if(pj.b > PartState.PISTON_)
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


/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Force,RV,Af) schedule(static) nowait
		for (size_t ii = start; ii < end; ++ii)
		{
			Particle const& pi = pnp1[ii];
			// size_t const size = outlist[ii].size();
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];

			
			real Rrhoi = 0.0;
			real Rrhod = 0.0;			
			real curve = 0.0;
			real correc = 0.0;

			StateVecD Vdiff = StateVecD::Zero();
			real Pbasei = 0.0;

			StateVecD aero = StateVecD::Zero();
			StateVecD RVi = StateVecD::Zero();
			StateVecD artViscI = StateVecD::Zero();


			if( dp.lam_ng[ii] < 0.75 /* pi.surf == 1 */ && pi.b == PartState.FREE_)
			{
				if (svar.Asource == 1)
				{
					Vdiff = (pi.cellV) - /*pi.v*/ dp.avgV[ii];
					// Pbasei = pi.cellP - avar.pRef;
				}
				else if (svar.Asource == 2)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - /*pi.v*/ dp.avgV[ii];
				}
				else 
				{
					Vdiff = avar.vInf - /* pi.v */ dp.avgV[ii];
				}

				#if SIMDIM == 3
					if(svar.Asource == 3)
					{	
						StateVecD Vel = svar.vortex.getVelocity(pi.xi);
						Vdiff = Vel - /*pi.v*/ dp.avgV[ii];
					}
				#endif
				// cout << ii << endl;
				aero = CalcAeroForce(avar,pi,Vdiff,dp.norm[ii],dp.lam_ng[ii],Pbasei);
				
			}

			for (std::pair<size_t,real> const& jj : outlist[ii])
			{	/* Neighbour list loop. */
				Particle const& pj = pnp1[jj.first];
				/*Check if the position is the same, and skip the particle if yes*/
				if(ii == jj.first)
				{
					if(/* pi.surf == 1 */ dp.lam_ng[ii] < 0.75 && pi.b == PartState.FREE_)
					{
						RV[ii] += aero / pj.rho * fvar.correc/dp.kernsum[ii];
						Af[ii] += aero / pj.rho * fvar.correc/dp.kernsum[ii];
						Force += aero / pj.rho * fvar.correc/dp.kernsum[ii];
					}
					
					continue;
				}

				StateVecD const Rij = pj.xi-pi.xi;
				StateVecD const Vij = pj.v-pi.v;
				real const rr = jj.second;
				real const r = sqrt(rr);	
				real const volj = pj.m/pj.rho;
				// StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);/* gradK;*/

				// if(gradK == StateVecD::Zero())
				// {
				// 	cout << ii << "  " << jj.first << endl;
				// 	cout << rr << "  " << r << "  " << r/fvar.H << endl;
 				// }

				if( /* dp.lam[ii] < 0.75 */ pi.surf == 1 )
				{
					// if(/* svar.iter % 10 == 0  &&*/ r > 1e-6*fvar.H)
					// 	gradLK = dp.L[ii]*gradK;

					if( /* dp.lam[pj.partID] < 0.75 */ pj.surf == 1 )
					{	
						curve -= (dp.norm[pj.partID].normalized()-dp.norm[ii].normalized()).dot(volj*gradK);
						correc += volj * Kernel(r,fvar.H,fvar.correc)/*/dp.kernsum[ii]*/;
					}	
				}
				
				if( dp.lam_ng[ii] < 0.75 /* pi.surf == 1 */ && pi.b == PartState.FREE_ && 
					pj.b > PartState.PISTON_ && pj.b != PartState.GHOST_)
				{	/* Assumed that particle masses are the same */
					RV[pj.partID] += aero / pj.rho * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
					Af[pj.partID] += aero / pj.rho * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
					Force += aero / pj.rho * Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				}	

				
				StateVecD const	contrib = BasePos(pi,pj,gradK);

				#ifndef NOALE
					StateVecD const ALEcontrib = ALEMomentum(pi, pj, volj, gradK, fvar.rho0);
				#endif

				StateVecD const aVisc = ArtVisc(fvar.nu,pi,pj,fvar,Rij,Vij,rr,gradK);

				/*Laminar Viscosity - Morris (2003)*/
				// StateVecD const visc = Viscosity(fvar.nu,fvar.HSQ,pi,pj,Rij,Vij,r,gradK);

				/*Base WCSPH continuity drho/dt*/
				#ifdef NOALE
					Rrhoi -= pj.m*(Vij.dot(gradK));
				#else
					Rrhoi -= ALEContinuity(pi,pj,volj,gradK);
				#endif
				Rrhod -= Continuity_dSPH(Rij,rr,fvar.HSQ,gradK,volj,dp.gradRho[ii],dp.gradRho[pj.partID],pi,pj);
		
				/*Surface Tension - Nair & Poeschel (2017)*/
				// SurfC = SurfaceTens(fvar,pj,Rij,r,npd);
				
				// SurfC = HeST(fvar,pi,pj,Rij,r,cgrad[ii],cgrad[pj.partID]);	

				// RVi += (-contrib + fvar.artMu * fvar.dMom * aVisc)/pi.rho/* + ALEcontrib*/ +  pj.m*visc;
				#ifdef NOALE
					RVi -= pj.m*contrib;
				#else
					RVi -= pj.m*contrib  - ALEcontrib;
				#endif

				artViscI += pj.m*aVisc;
				// viscI += pj.m*visc;

			}/*End of neighbours*/
			
			// if(ii - start < 36)
			// {
			// #pragma omp critical
			// cout << ii - start << "  " << RVi(0) << "   " << RVi(1) << "  " << artViscI(0) << "   " << artViscI(1) << "  " << viscI(0) << "   " << viscI(1) << endl; 
			// }

			if(pi.internal == 1)
			{	// Apply the normal boundary force
				RVi += NormalBoundaryRepulsion(fvar, cells, pi);
			}

			if(/* dp.lam[ii] < 0.75 */ pi.surf == 1)
			{
				RVi += (fvar.sig/pi.rho * curve * dp.norm[ii].normalized())/correc;
			}
			
			// #pragma omp critical
			// Force += aero;	

			
			if(pi.b == PartState.FREE_ )
				RV[ii] += (RVi + artViscI /* + viscI */ /* + aero/pi.m */ + g);
			else
				RV[ii] += (RVi + artViscI /* + viscI */ /* + aero/pi.m */);

			Rrho[ii] = (Rrhoi - fvar.dCont * Rrhod);

			// res_.emplace_back(aero/pi.m);
			// Rrho_.emplace_back(0.0);

		} /*End of sim parts*/		

	}	/*End of declare parallel */
 
    // cout << RV.size() << "  " << Rrho.size() << "  " << Af.size() << "  " << curv.size() << endl;
}

#endif