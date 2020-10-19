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
#include "Eigen/LU"
#include "Eigen/Eigenvalues"

/*Numerical particle density for Nair & Poeschel (2017) surface tension*/
real GetNumpartdens(SIM const& svar, FLUID const& fvar, State const& pnp1, outl const& outlist)
{
	real npd = 0.0;
	uint const& end = svar.totPts;
	#pragma omp parallel for reduction(+:npd)
	for (uint ii=0; ii< end; ++ii)
	{
		StateVecD const& pi = pnp1[ii].xi;
		for (auto jj:outlist[ii])
		{ /* Surface Tension calcs */
			StateVecD const& pj = pnp1[jj].xi;
			real const r = (pj-pi).norm();
			npd += W2Kernel(r,fvar.H,fvar.correc);
		}
	}
	return npd/real(svar.totPts);
}

/*Surface Tension - Nair & Poeschel (2017)*/
StateVecD SurfaceTens(FLUID const& fvar, Part const& pj, StateVecD const& Rij, 
					  real const& r, real const& npd)
{
	/*Surface tension factor*/
	real const lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));

	real fac=1.0; /*Boundary Correction Factor*/
	if(pj.b==PartState.BOUND_) fac=(1+0.5*cos(M_PI*(fvar.contangb/180))); 

	/*npd = numerical particle density (see code above) */
	real const sij = 0.5*pow(npd,-2.0)*(fvar.sig/lam)*fac;
	return -(Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*fvar.H));
}

/* Colour field gradient in He et al (2014) method, note the bottom variable to account for no air*/
std::vector<StateVecD> GetColourGrad(SIM const& svar, FLUID const& fvar, State const& pnp1, outl const& outlist)
{
	std::vector<StateVecD> cgrad(svar.totPts, StateVecD::Zero());
	const uint& start = svar.bndPts;
	const uint& end = svar.totPts;

	#pragma omp parallel for shared(outlist) reduction(+:cgrad)
	for(uint ii=start; ii < end; ++ii)
	{
		StateVecD const& pi = pnp1[ii].xi;
		real bottom = 0.0;
		
		for(auto jj:outlist[ii])
		{	/*Find the denominator to correct absence of second phase*/
			Part const& pj = pnp1[jj];
			real const r = (pj.xi-pi).norm();
			bottom +=(pj.m/pj.rho)*W2Kernel(r,fvar.H,fvar.correc);
		}

		for(auto jj:outlist[ii])
		{	/*Find the numerator and sum*/
			Part const& pj = pnp1[jj];
			StateVecD const Rij = pj.xi-pi;
			real const r = Rij.norm();
			cgrad[ii] += (pj.m/pj.rho)*W2GradK(Rij,r,fvar.H,fvar.correc);
		}
		
		cgrad[ii] = cgrad[ii]/bottom;
	}

	return cgrad;
}

/*Diffuse interface method from He et al (2014) for surface tension*/
StateVecD HeST(FLUID const& fvar, Part const& pi, Part const& pj, 
			  StateVecD const& Rij, real const& r, StateVecD const& cgradi, StateVecD const& cgradj)
{
	return (0.01/2.0)*(pi.m/pi.rho)*(pj.m/pj.rho)*((cgradi.squaredNorm()+cgradj.squaredNorm())/2.0)
		*W2GradK(Rij,r,fvar.H,fvar.correc);
}

/*Base WCSPH momentum equation*/
StateVecD Base(FLUID const& fvar, Part const& pi, Part const& pj, 
	 StateVecD const& Rij, StateVecD const& Vij, real const r, 
	 StateVecD const& gradK, std::vector<real> &mu)
{
	/*Pressure and artificial viscosity - Monaghan (1994) p.400*/
	real const vdotr = Vij.dot(Rij);
	real const muij= fvar.H*vdotr/(r*r+0.0001*fvar.HSQ);
	real pifac;
	if (vdotr > 0.0) 
	{
		pifac = 0.0;
	}
	else
	{
		real const rhoij = 0.5*(pi.rho+pj.rho);
		real const cbar = 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));
		pifac = fvar.alpha*cbar*muij/rhoij;
	}
	
	mu.emplace_back(muij);
	return gradK*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));
}

/*Laminar Viscosity - Morris (2003)*/
StateVecD Viscosity(FLUID const& fvar, Part const& pi, Part const& pj, 
	StateVecD const& Rij, StateVecD const& Vij, real const& r, StateVecD const& gradK)
{
	return -Vij*((fvar.mu)/(pi.rho*pj.rho))*(1.0/(r*r+0.01*fvar.HSQ))*Rij.dot(gradK);
}


/*Repulsion for interacting with mesh surface - saves generating particles on surface*/
StateVecD NormalBoundaryRepulsion(FLUID const& fvar, MESH const& cells, Part const& pi)
{
    real beta = 4*fvar.Cs*fvar.Cs;
    real kern = BoundaryKernel(pi.y,fvar.H,beta);
	return fvar.bndM/(fvar.bndM+fvar.simM)*kern*pi.bNorm;
}

/*L matrix for delta-SPH calculation*/
void dSPH_PreStep(FLUID const& fvar, size_t const& start, size_t const& end, State const& pnp1, 
		vector<vector<Part>> const& neighb, vector<StateVecD>& gradRho,  vector<real>& isSurf, vector<StateMatD>& Lmat)
{
	/************   RENORMALISATION MATRIX CALCULATION    ************/
	
	#pragma omp parallel
	{
		#pragma omp for
		for(size_t ii = 0; ii < end; ++ii)
		{
			StateMatD Lmat_ = StateMatD::Zero();
			StateVecD gradRho_ = StateVecD::Zero();
			
			Part const& pi(pnp1[ii]);

			for(Part const& pj:neighb[ii])
			{
				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				if(pj.b != PartState.GHOST_)
				{
					StateVecD const Rij = pj.xi - pi.xi;
					real const r = Rij.norm();
					real volj = pj.m/pj.rho;
					StateVecD Grad = W2GradK(Rij,r,fvar.H,fvar.correc);

					Lmat_ += volj * Rij * Grad.transpose();

					gradRho_ += volj * (pj.rho - pi.rho)*Grad;	
				}
			}

			Eigen::FullPivLU<StateMatD> lu(Lmat_);
			/*If only 1 particle is in the neighbourhood, tensile instability can occur*/
			/*In this case, L becomes singular, and invertability needs to be checked.*/
			if(lu.isInvertible())
			{	
				gradRho[ii] = lu.inverse() * gradRho_;
				Eigen::SelfAdjointEigenSolver<StateMatD> es;
		        es.computeDirect(Lmat_);
		        isSurf[ii] = (es.eigenvalues()).minCoeff(); //1 inside fluid, 0 outside fluid

		        Lmat[ii] = lu.inverse();
			}
			else
			{	/*In a sparse neighbourhood if singular, so mark as a surface*/
				isSurf[ii] = 0.0;
				Lmat[ii](0,0) = 1;
			}
		}
	}

}

/*Di term for delta-SPH calculation*/
void dSPH_Di_Term(FLUID const& fvar, size_t const& start, size_t const& end, State const& pnp1, 
		outl const& outlist, vector<real>& Di,  vector<real>& isSurf, vector<StateMatD>& Lmat)
{
	/************   RENORMALISATION MATRIX CALCULATION    ************/
	vector<StateVecD> gradRho(end, StateVecD::Zero());
	// vector<real> isSurf(end, 0.0);

	#pragma omp parallel
	{
		#pragma omp for
		for(size_t ii = 0; ii < end; ++ii)
		{
			StateMatD Lmat_ = StateMatD::Zero();
			StateVecD gradRho_ = StateVecD::Zero();
			
			Part const& pi(pnp1[ii]);

			for(size_t const& jj:outlist[ii])
			{
				Part const& pj(pnp1[jj]);
				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				if(pj.b != PartState.GHOST_)
				{
					StateVecD const Rij = pj.xi - pi.xi;
					real const r = Rij.norm();
					real volj = pj.m/pj.rho;
					StateVecD Grad = W2GradK(Rij,r,fvar.H,fvar.correc);

					Lmat_ += volj * Grad * Rij.transpose();

					gradRho_ += volj * (pj.rho - pi.rho)*Grad;	
				}
			}

			// Lmat[ii] = Lmat_;

			Eigen::FullPivLU<StateMatD> lu(Lmat_);
			/*If only 1 particle is in the neighbourhood, tensile instability can occur*/
			/*In this case, L becomes singular, and invertability needs to be checked.*/
			if(lu.isInvertible())
			{	
				gradRho[ii] = lu.inverse() * gradRho_;
				Eigen::SelfAdjointEigenSolver<StateMatD> es;
		        es.computeDirect(Lmat_);
		        isSurf[ii] = (es.eigenvalues()).minCoeff(); //1 inside fluid, 0 outside fluid

		        Lmat[ii] = lu.inverse();
			}
			else
			{	/*In a sparse neighbourhood if singular, so mark as a surface*/
				isSurf[ii] = 0.0;
				Lmat[ii](0,0) = 1;
			}
		}

		#pragma omp for
		for(size_t ii = 0; ii < end; ++ii)
		{
			real Di_ = 0.0;
			Part const& pi(pnp1[ii]);

			for(size_t const& jj:outlist[ii])
			{
				Part const& pj(pnp1[jj]);
				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				if(pj.b != PartState.GHOST_)
				{

					
					StateVecD const Rij = pj.xi - pi.xi;
					real const r = Rij.norm();
					real const rsq = Rij.squaredNorm();
					StateVecD gradK = W2GradK(Rij,r,fvar.H,fvar.correc);

					StateVecD psi_ij = ((pi.rho - pj.rho) + 0.5*(gradRho[pi.partID] + gradRho[pj.partID]).dot(Rij)) * Rij/(rsq);

					Di_+= 2*psi_ij.dot(gradK) * pj.m/pj.rho;
				}
			}
			Di[ii] = Di_;
		}

	}

}

/* delta-SPH dissipation term in the continuity equation*/
real Continuity_dSPH(StateVecD const& Rij, real const& r, StateVecD const& Grad, 
				real const& volj, vector<StateVecD> const& gradRho, Part const& pi, Part const& pj)
{
	StateVecD const psi_ij = ((pi.rho - pj.rho) + 0.5*(gradRho[pi.partID] + gradRho[pj.partID]).dot(Rij)) * Rij/(r*r);
	return volj * psi_ij.dot(Grad);
}

StateVecD Momentum_dSPH(FLUID const& fvar, StateVecD const& Rij, real const& r, StateVecD const& Grad, 
						real const& volj, Part const& pi, Part const& pj)
{
	real const pi_ij = (pj.v-pi.v).dot(Rij)/(r*r);
	return fvar.dMom * volj * pi_ij * Grad;
}


///**************** RESID calculation **************
void Forces(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, State const& pnp1/*, State& airP*/,
	 vector<vector<Part>> const& neighb, outl const& outlist, 
	 // vector<real> const& Di, vector<real> const& isSurf,
	 vector<StateVecD>& RV, vector<real>& Rrho, std::vector<StateVecD>& Af,
	 StateVecD& Force, vector<real>& kern, vector<StateVecD>& normal/*, vector<real>& curve*/)
{
	svar.maxmu=0; 					    /* CFL Parameter */
	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;

	// #pragma omp critical
	// cout << endl << endl;

	/*Gravity Vector*/
	#if SIMDIM == 3
		StateVecD const g = StateVecD(0.0,0.0,0.0);
	#else
		StateVecD const g = StateVecD(0.0,-9.81);
	#endif
	// const uint piston = svar.psnPts;

	/********* LOOP 1 - all points: Calculate numpartdens ************/
	// wDiff = GetWeightedDistance(svar,fvar,pnp1,outlist);
	// real const npd = GetNumpartdens(svar, fvar, pnp1, outlist);
	// std::vector<StateVecD> const cgrad = GetColourGrad(svar,fvar,pnp1,outlist);
	// std::vector<StateVecD> ST(svar.totPts,StateVecD::Zero()); /*Surface tension force*/		

	vector<StateVecD> gradRho(end,StateVecD::Zero());
	vector<real> isSurf(end,0.0);
	vector<StateMatD> Lmat(end,StateMatD::Zero());
	vector<StateVecD> avgV(end,StateVecD::Zero());

	dSPH_PreStep(fvar,start,end,pnp1,neighb,gradRho,isSurf,Lmat);

	#pragma omp parallel /*shared(Af, wDiff, norm, curve)*/
	{
		// real npd = 0.0;
		// vector<real> wDiff(end,0.0);
		vector<StateVecD> norm(end, StateVecD::Zero());
		vector<real> kernsum(end,0.0);
		// vector<real> curve(end,0.0);

/******** LOOP 2 - Piston points: Calculate density and pressure. **********/
		// #pragma omp for reduction(+:Rrhocontr)
		// for (uint ii=0; ii < start; ++ii)
		// {
		// 	Part pi = pnp1[ii];

		// 	for(auto j:outlist[ii])
		// 	{
		// 		Part pj = pnp1[j];
		// 		StateVecD Rij = pj.xi-pi.xi;
		// 		StateVecD Vij = pj.v-pi.v;
		// 		real r = Rij.norm();
		// 		StateVecD Grad = W2GradK(Rij, r,fvar.H,fvar.correc);
		// 		Rrhocontr[ii] -= pj.m*(Vij.dot(Grad));
		// 	}
		// 	pnp1[ii].Rrho = Rrhocontr[ii]; /*drho/dt*/
		// }



/******** LOOP 3 - Boundary points: Calculate density and pressure. **********/		
		#pragma omp for reduction(+:RV, Rrho) /*Reduction defs in Var.h*/
		for (size_t ii=0; ii < start; ++ii)
		{
			Part const pi = pnp1[ii];
			uint const size = outlist[ii].size();
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];

			std::vector<real> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.reserve(size+1);			
			
			for (Part const& pj:neighb[ii])
			{	/* Neighbour list loop. */

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				StateVecD const Rij = pj.xi-pi.xi;
				StateVecD const Vij = pj.v-pi.v;
				real const r = Rij.norm();
				StateVecD const Grad = W2GradK(Rij, r,fvar.H, fvar.correc);

				StateVecD contrib = Base(fvar,pi,pj,Rij,Vij,r,Grad,mu);
				/*drho/dt*/
				Rrho[ii] -= pj.m*(Vij.dot(Grad));
				RV[ii] += pj.m*contrib;
			}/*End of neighbours*/
			
		} /*End of boundary parts*/

		

/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Rrho, RV /*, ST*/) schedule(static) nowait
		for (size_t ii = start; ii < end; ++ii)
		{
			Part const pi = pnp1[ii];
			size_t const size = outlist[ii].size();
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];

			vector<real> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.reserve(size+1);
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/		
			// real mj = 0.0;
			// StateVecD mRij = StateVecD::Zero();

			StateVecD RVi = StateVecD::Zero();
			real Rrhoi = 0.0;
			real Rrhod = 0.0;
			StateVecD normi = StateVecD::Zero();			
			StateVecD avgV_ = StateVecD::Zero();

			for (Part const& pj:neighb[ii])
			{	/* Neighbour list loop. */

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
				{
					// cout << "PartID is the same: " << pi.partID << endl;
					continue;
				}

				StateVecD const Rij = pj.xi-pi.xi;
				StateVecD const Vij = pj.v-pi.v;
				real const r = Rij.norm();
				real const volj = pj.m/pj.rho;
				kernsum[ii] += W2Kernel(r,fvar.H,fvar.correc);

				StateVecD const gradK = W2GradK(Rij,r,fvar.H,fvar.correc);

				/*Momentum contribution - Monaghan (1994)*/
				StateVecD const contrib = Base(fvar,pi,pj,Rij,Vij,r,gradK,mu);
				StateVecD visc = StateVecD::Zero();
				StateVecD SurfC = StateVecD::Zero();

				// if (pj.b != PartState.GHOST_)
				// {
					/*Laminar Viscosity - Morris (2003)*/
				    visc = Viscosity(fvar,pi,pj,Rij,Vij,r,gradK);

					/*Surface Tension - Nair & Poeschel (2017)*/
					// SurfC = SurfaceTens(fvar,pj,Rij,r,npd);
					
					// SurfC = HeST(fvar,pi,pj,Rij,r,cgrad[ii],cgrad[pj.partID]);
					
					// Find the normal vectors
					normi += volj * gradK;

					/*Base WCSPH continuity drho/dt*/
					Rrhoi -= pj.m*(Vij.dot(gradK));
					Rrhod -= Continuity_dSPH(Rij, r, gradK, volj, gradRho, pi, pj);
				// }

				/*Calculate aerodynamic pressure correction as He et al (2014)*/
				// StateVecD const Apress = pi.cellP/pi.rho * pj.m/pj.rho * gradK;

				avgV_ += pj.v * volj * W2Kernel(r,fvar.H,fvar.correc);

				// npd += W2Kernel(r,fvar.H,fvar.correc);
				RVi += pj.m*contrib + pj.m*visc /*+ Apress*/ + SurfC/pj.m;
				
				// ST[ii] += SurfC/pj.m;
			}/*End of neighbours*/
			
			if(pi.internal == 1)
			{	// Apply the normal boundary force
				RVi += NormalBoundaryRepulsion(fvar, cells, pi);
			}


			RV[ii] = RVi + g;
			Rrho[ii] = Rrhoi + fvar.dCont * Rrhod;

			norm[ii] = Lmat[ii] * normi/**isSurf[ii]*/;
			normal[ii] =  Lmat[ii] * normi;
			// wDiff[ii] = mRij.norm()/(fvar.H*mj);

			avgV[ii] = avgV_;

			// kern[ii] = kernsum[ii];
			// kern[ii] = wDiff[ii];
			// kern[ii] = Lmat[ii].determinant();

			//CFL f_cv Calc
			real it = *max_element(mu.begin(),mu.end());
			if(it > svar.maxmu)
				svar.maxmu = it;

		} /*End of sim parts*/		
	


		// StateVecD aeroF = StateVecD::Zero();
		// Get the surface curvatures
		#pragma omp for schedule(static) nowait /*reduction(+:RV, curve)*/
		for (size_t ii = start; ii < end; ++ii) 
		{
			real woccl = 0.0;			
			Part const& pi = pnp1[ii];
			size_t size = outlist[ii].size();
			
			kern[ii] = isSurf[ii];

			// Check if interaction is between surface particles
			if(isSurf[ii] < 0.66)
			{	
				real curve = 0.0;
				real correc = 0.0;
				StateVecD aeroD = StateVecD::Zero();
				// StateVecD normi = StateVecD::Zero();

				StateVecD Vdiff;
				real Pbasei = 0.0;
				if (svar.Asource == 1)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - /*pi.v*/ avgV[ii];
					Pbasei = pi.cellP - avar.pRef;
				}
				else if (svar.Asource == 2)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - /*pi.v*/ avgV[ii];
				}
				else 
				{
					Vdiff = avar.vInf - /*pi.v*/ avgV[ii];
				}

#if SIMDIM == 3
				if(svar.Asource == 3)
				{	
					StateVecD Vel = svar.vortex.getVelocity(pi.xi);
					Vdiff = Vel- /*pi.v*/ avgV[ii];
				}
#endif


				for (size_t const& jj:outlist[ii])
				{	/* Neighbour list loop. */
					Part const& pj = pnp1[jj];
					StateVecD const Rij = pj.xi - pi.xi;
					
					real const volj = (pj.m/pj.rho);
					StateVecD const divW = volj*W2GradK(Rij, Rij.norm(), fvar.H, fvar.correc);
					/*Surface Tension - Ordoubadi (2014) CSF*/
					if(isSurf[jj] < 0.66)
					{	
						curve -= (norm[jj].normalized()-norm[ii].normalized()).dot(Lmat[ii] *divW);
						correc += volj * W2Kernel(Rij.norm(),fvar.H,fvar.correc);
					}

					/*Occlusion for Gissler Aero Method*/
					if (pi.b == PartState.FREE_ && (avar.acase == 4 || avar.acase == 1))
					{
						real const frac = -Vdiff.normalized().dot(Rij.normalized());
						
						// cout << num/denom << endl;
						if (frac > woccl)
						{
							woccl = frac;
						}
					}
				}

				// real theta = norm[ii].normalized().dot(Vdiff.normalized());

				// real Cp = 1-4*pow(sin(theta),2);

				// real aPress = 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;


				// Rrho[ii] += 1e-2*(pi.p - aPress);
				// RV[ii] += (fvar.sig/pi.rho * curve * norm[ii])/correc;
				
				// aeroD *= Pi/pi.rho;
				// aeroD *= pi.m/pi.rho* pi.m/pi.rho;

				// #pragma omp critical
				// cout << aeroD[0] << "  " << aeroD[1] << "  " << endl;
				RV[ii] += aeroD/*/correc*//*/pi.m*/;
				Force += aeroD*pi.m;
				// normal[ii] = normi;
				// aeroF += (1-woccl)*aeroD/correc;
				
	// 			/*Find the aero force*/
				if(pi.b == PartState.FREE_)
				{
			
					StateVecD aero = CalcAeroForce(avar,pi,Vdiff,norm[ii],size,woccl,Pbasei);
					
					// if (avar.acase == 6)
					// {
					// 	real theta = norm[ii].normalized().dot(Vdiff.normalized());
									
					// 	real Cp = 0.0;
					// 	if(abs(theta) < 2.4435)
					// 	{
					// 		Cp = 1.1*cos(theta*2.03)-0.1;
					// 	}
					// 	else
					// 	{
					// 		Cp = 0.075;
					// 	}

					// 	real Plocali = 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;

					// 	real pTarget = Plocali + Pbasei;

					// 	Rrho[ii] -= 1000.0 * (pi.p - pTarget);
					// }

					RV[ii] += aero;
					Force += aero*pi.m;	

				}
			}
			

		}
		// cout << aeroF << endl;
		// Force = aeroF;
	}	/*End of declare parallel */
	
}

#endif