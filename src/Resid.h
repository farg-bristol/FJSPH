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
#include "Third_Party/Eigen/LU"
#include "Third_Party/Eigen/Eigenvalues"
#include <unordered_set>

/*L matrix for delta-SPH calculation*/
void dSPH_PreStep(SIM const& svar, FLUID const& fvar,
			size_t const& start, size_t const& end, State const& pnp1, 
				 outl const& outlist, DELTAP& dp)
{
	/************   RENORMALISATION MATRIX CALCULATION    ************/
	dp.clear();
	
	vector<StateMatD> Lmat(end);
	vector<StateVecD> gradRho(end);
	vector<StateVecD> norm(end);
	vector<StateVecD> avgV(end);
	
	vector<real> lam(end);
	vector<real> lam_nb(end);
	vector<real> kernsum(end);

	#pragma omp parallel shared(Lmat,gradRho,norm,lam,lam_nb,kernsum)
	{
		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < end; ++ii)
		{
			StateMatD Lmat_ = StateMatD::Zero();
			StateMatD Lmat_nb_ = StateMatD::Zero();
			StateVecD gradRho_ = StateVecD::Zero();
			StateVecD norm_ = StateVecD::Zero();
			StateVecD avgV_ = StateVecD::Zero();
			real kernsum_ = 0.0;
		
			Part const& pi(pnp1[ii]);
	
			for(std::pair<size_t,real> const& jj:outlist[ii])
			{
				Part const& pj(pnp1[jj.first]);
				/*Check if the position is the same, and skip the particle if yes*/
				if(ii == jj.first)
				{
					kernsum_ +=  pj.m/pj.rho * fvar.correc;
					avgV_ += pj.v * pj.m/pj.rho * fvar.correc;

					continue;
				}

				StateVecD const Rij = pj.xi - pi.xi;
				real const r = sqrt(jj.second);
				real const volj = pj.m/pj.rho;
				real const kern = Kernel(r,fvar.H,fvar.correc);
				StateVecD const Grad = GradK(-Rij,r,fvar.H,fvar.correc);

				Lmat_ -= volj * Rij * Grad.transpose();

				gradRho_ += volj * (pj.rho - pi.rho)*Grad;

				if(pj.b != PartState.BOUND_)
				{
					Lmat_nb_ -= volj * Rij * Grad.transpose();

					norm_ += volj*Grad;	

					kernsum_ += volj * kern;

					avgV_ += pj.v * volj * kern;
				}
			}

			Eigen::FullPivLU<StateMatD> lu(Lmat_);
			/*If only 1 particle is in the neighbourhood, tensile instability can occur*/
			/*In this case, L becomes singular, and invertability needs to be checked.*/
			if(lu.isInvertible())
			{	
				StateMatD Linv = lu.inverse();
				Lmat[ii] = Linv;
				gradRho[ii] = Linv * gradRho_ ;

				// StateVecD ntemp = Linv * norm_;
				// norm[ii] = (norm_.norm() * ntemp.normalized());
				norm[ii] = Linv.normalized() * norm_;
				// avgV[ii] = lu.inverse() * avgV_;

				Eigen::SelfAdjointEigenSolver<StateMatD> es;
		        es.computeDirect(Lmat_);
		        lam[ii] = ((es.eigenvalues()).minCoeff()); //1 inside fluid, 0 outside fluid    
			}
			else
			{	/*In a sparse neighbourhood if singular, so mark as a surface*/
				StateMatD Ltemp = StateMatD::Zero();
				Ltemp(0,0) = 1.0;
				Lmat[ii] = (Ltemp);
				gradRho[ii] = (gradRho_);
				norm[ii] = (norm_);
				// avgV[ii] = avgV_;
				lam[ii] = (0.0);
			}

			/*Lmatrix computation for particle shifting, ignoring the boundary*/
			Eigen::FullPivLU<StateMatD> lu2(Lmat_nb_);
			if(lu2.isInvertible())
			{	
				Eigen::SelfAdjointEigenSolver<StateMatD> es;
		        es.computeDirect(Lmat_nb_);
		        lam_nb[ii] = ((es.eigenvalues()).minCoeff()); //1 inside fluid, 0 outside fluid    
			}
			else
			{	/*In a sparse neighbourhood if singular, so mark as a surface*/
				lam_nb[ii] = (0.0);
			}

			if(pi.b == PartState.BOUND_)
			{
				lam[ii] = 1.0;
				lam_nb[ii] = 1.0;
			}
	
			avgV[ii] = (avgV_/kernsum_);
			kernsum[ii] = (kernsum_);	
		}

		// #pragma omp for schedule(static) 
		// for(size_t ii = 0; ii < end; ++ii)
		// {
		// 	// StateVecD gradRho_ = StateVecD::Zero();
		// 	StateVecD norm_ = StateVecD::Zero();
		// 	StateVecD avgV_ = StateVecD::Zero();
		// 	Part const& pi(pnp1[ii]);

		// 	for(size_t const& jj:outlist[ii])
		// 	{
		// 		Part const& pj(pnp1[jj]);
		// 		/*Check if the position is the same, and skip the particle if yes*/
		// 		if(pi.partID == pj.partID)
		// 			continue;

		// 		StateVecD const Rij = pj.xi - pi.xi;
		// 		real const r = Rij.norm();
		// 		real const volj = pj.m/pj.rho;
		// 		real const kern = Kernel(r,fvar.H,fvar.correc);
				
		// 		// gradRho_ += kern * volj * gradRho[jj];
				
		// 		// Apply a shepard filter to the normal vectors
		// 		if(pj.b != PartState.BOUND_)
		// 		{
		// 			norm_ +=  kern * volj * norm[jj];

		// 			avgV_ += pj.v * volj * kern;
		// 		}
		// 	}

		// 	// gradRho[ii] = gradRho_/kernsum[ii];
		// 	// if(lam[ii] > 0.4)
		// 	// {
		// 		norm[ii] = norm_/kernsum[ii];
		// 		avgV[ii] = avgV_;
		// 	// }
		// 	// else
		// 	// {
		// 	// 	avgV[ii] = pi.v;
		// 	// }
			

		// }			

	}	/*End parallel section*/	

	dp.L = Lmat;
	dp.gradRho = gradRho;
	dp.norm = norm;
	dp.avgV = avgV;
	dp.lam = lam;
	dp.lam_nb = lam_nb;
	dp.kernsum = kernsum;
}

#ifndef NOALE
void Particle_Shift(SIM const& svar, FLUID const& fvar, size_t const& start, size_t const& end, outl const& outlist,
DELTAP const& dp, State& pnp1)
{
	/*Implement particle shifting technique - Sun, Colagrossi, Marrone, Zhang (2018)*/

	// vector<Particle>::iterator maxUi = std::max_element(pnp1.begin(),pnp1.end(),
	// 	[](Particle p1, Particle p2){return p1.v.norm()< p2.v.norm();});

	// real const maxU = fvar.maxU/**maxUi->v.norm()*/;

	// real const maxU = fvar.maxU;

	#pragma omp parallel
	{
		#pragma omp for schedule(static) nowait
		for(size_t ii = start; ii < end; ++ii)
		{	
			/*Check if the particle is too close to the surface, and ignore them if so*/
			if(dp.lam_nb[ii] < 0.55)
				continue;

			Part const& pi(pnp1[ii]);
			
			StateVecD deltaU = StateVecD::Zero();
			StateVecD gradLam = StateVecD::Zero();
			real maxUij = 0.0;

			uint f = 0;
			real woccl = 0.0;

			for (std::pair<size_t,real> const& jj:outlist[ii])
			{	/* Neighbour list loop. */
				Part const& pj(pnp1[jj.first]);

				if(pj.partID == pi.partID /*|| pj.b == PartState.BOUND_*/)
					continue;

				StateVecD const Rij = pj.xi-pi.xi;
				real const r = sqrt(jj.second);
				real const volj = pj.m/pj.rho;
				real const kern = Kernel(r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);

				deltaU += ( 1.0 + 0.2 * pow(kern/fvar.Wdx,4.0) ) * gradK * volj;
				// deltaU += 
				gradLam -= (dp.lam[jj.first]-dp.lam[ii]) * dp.L[ii] * gradK * volj;

				if(/* dp.lam[jj.first] < 0.75 || */ pj.surf == 1)
					f = 1;

				real theta = acos(dp.norm[ii].normalized().dot(dp.norm[jj.first].normalized()));
				if ( theta > woccl )
					woccl = theta;

				if ((pj.v -pi.v).norm() > maxUij )
				{
					maxUij = (pj.v-pi.v).norm();
				}
			}

			// deltaR *= -1 * fvar.sr * maxU / fvar.Cs;
			deltaU *= -2.0 * fvar.H * pi.v.norm();

			deltaU = std::min(deltaU.norm(), maxUij/2.0) * deltaU.normalized();

			if(pi.b != PartState.START_ && pi.b != PartState.BACK_)
			{
				// gradLam = gradLam.normalized();
				// cout << gradLam(0) << "  " << gradLam(1) << endl;
				// StateVecD norm = dp.norm[ii].normalized();
				StateVecD norm = gradLam.normalized();
				
				/*Apply the partial shifting*/
				if (f == 0)
				{   /*Particle in the body of fluid, so treat as normal*/
					pnp1[ii].vPert = deltaU;
				}
				else /*Already checked if eigenvalue is less than 0.55*/
				{	/*Particle in the surface region, so apply a partial shifting*/
				    if(norm.dot(deltaU) >= 0.0)
					{	/*Check if particle is heading towards or away from the surface*/
						if(woccl < M_PI/12.0)
						{
							StateMatD shift = StateMatD::Identity() - norm*norm.transpose();
							pnp1[ii].vPert = /*dp.lam[ii]*/ shift*deltaU;
						}
						else
						{
							pnp1[ii].vPert = StateVecD::Zero();
							continue;
						}
					}
					else
					{
						pnp1[ii].vPert = deltaU;
					}
				}
				/*Otherwise, particle is near a surface, and don't shift.*/
			}
			else
			{
				pnp1[ii].vPert = StateVecD::Zero();
			}	
		}
	}

	// for(size_t ii = start; ii < end; ++ii)
	// {
	// 	cout << "L Mag: " << dp.L[ii].determinant() << "  gradRho: " << dp.gradRho[ii].norm() << "  norm: " << dp.norm[ii].norm() 
	// 	<< "  avgV: " << dp.avgV[ii].norm() << "  lam: " << dp.lam[ii] << "  kernsum: " << dp.kernsum[ii] << endl;
	// }

	// cout << "L Matrix: " << dp.L[400] << endl;
}
#endif

// void Apply_XSPH(FLUID const& fvar, size_t const& start, size_t const& end, 
// 				outl const& outlist, DELTAP const& dp, State& pnp1)
// {
// 	#pragma omp parallel
// 	{
// 		#pragma omp for schedule(static) nowait
// 		for(size_t ii = start; ii < end; ++ii)
// 		{	
// 			Part const& pi(pnp1[ii]);

// 			StateVecD vPert_ = StateVecD::Zero();
// 			for (size_t const& jj:outlist[ii])
// 			{	 Neighbour list loop. 
// 				Part const& pj(pnp1[jj]);

// 				if(pj.partID == pi.partID || pj.b == PartState.BOUND_)
// 				{
// 					continue;
// 				}
				
// 				real const r = (pj.xi-pi.xi).norm();;
// 				real const kern = Kernel(r,fvar.H,fvar.correc);
// 				real rho_ij = 0.5*(pi.rho + pj.rho);

// 				vPert_ -= 0.5 * pj.m/rho_ij * (pi.v-pj.v) * kern/*/dp.kernsum[ii]*/; 

// 			}

// 			pnp1[ii].vPert = vPert_;
// 		}
// 	}
// }


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

/* Arbitrary Lagrangian Eulerian formulation - Sun, Colagrossi, Marrone, Zhang (2018)*/
StateVecD ALEMomentum(Part const& pi, Part const& pj, real const& Vj, StateVecD const& gradK, real const& rho0)
{
	return 	(rho0/pi.rho) * (pj.v * pj.vPert.transpose() + pi.v * pi.vPert.transpose()) * gradK * Vj ;
			// - pi.v * (pj.vPert - pi.vPert).dot(gradK) * Vj ;
}

real ALEContinuity(Part const& pi, Part const& pj, real const& Vj, StateVecD const& gradK)
{
	return pi.rho*((pj.v+pj.vPert) - (pi.v + pi.vPert)).dot(gradK)*Vj - 
			(pj.rho*pj.vPert + pi.rho*pi.vPert).dot(gradK)*Vj ;
}

/* delta-SPH dissipation term in the continuity equation*/
real Continuity_dSPH(StateVecD const& Rij, real const& rr, real const& HSQ, StateVecD const& Grad, 
				real const& volj, StateVecD const& gRho_i, StateVecD const& gRho_j, Part const& pi, Part const& pj)
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

StateVecD BasePos(Part const& pi, Part const& pj, StateVecD const& gradK)
{
	/* Positive Pressure - Monaghan (1994)*/
	return gradK * (pj.p*pow(pj.rho,-2)+pi.p*pow(pi.rho,-2));
}

StateVecD BaseNeg(Part const& pi, Part const& pj, StateVecD const& gradK)
{
	/* Negative Pressure - Monaghan (1994)*/
	return gradK * (pj.p*pow(pj.rho,-2)-pi.p*pow(pi.rho,-2));
}

StateVecD ArtVisc(real const& nu, Part const& pi, Part const& pj, FLUID const& fvar, StateVecD const& Rij, StateVecD const& Vij, real const rr, 
	 StateVecD const& gradK)
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

/*Laminar Viscosity - Morris (2003)*/
/*Apparently divergent at free surfaces - consider removing.*/
// StateVecD Viscosity(real const& mu, real const& HSQ, Part const& pi, Part const& pj, 
// 	StateVecD const& Rij, StateVecD const& Vij, real const& r, StateVecD const& gradK)
// {
// 	return Vij*(mu/(pi.rho*pj.rho))*(1.0/(r*r+0.01*HSQ))*Rij.dot(gradK);
// }

StateVecD Viscosity(real const& nu, real const& HSQ, Part const& pi, Part const& pj, 
	StateVecD const& Rij, StateVecD const& Vij, real const& rr, StateVecD const& gradK)

{
	return nu * (pi.rho + pj.rho)/(pi.rho*pj.rho) * (Rij.dot(gradK))/(rr + 0.001*HSQ) * Vij;		
}


/*Repulsion for interacting with mesh surface - saves generating particles on surface*/
StateVecD NormalBoundaryRepulsion(FLUID const& fvar, MESH const& cells, Part const& pi)
{
    real beta = 4*fvar.Cs*fvar.Cs;
    real kern = BoundaryKernel(pi.y,fvar.H,beta);
	return fvar.bndM/(fvar.bndM+fvar.simM)*kern*pi.bNorm;
}



///**************** RESID calculation **************
void Forces(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, State const& pnp1/*, State& airP*/,
	 vector<vector<Part>> const& neighb, outl const& outlist, DELTAP const& dp,
	 // vector<real> const& Di, vector<real> const& isSurf,
	 vector<StateVecD>& RV, vector<real>& Rrho, std::vector<StateVecD>& Af, StateVecD& Force, vector<real>& curv)
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

	Rrho = vector<real>(end);
	RV = vector<StateVecD>(end);
	Af = vector<StateVecD>(end);
	curv = vector<real>(end);

	#pragma omp parallel shared(RV,Rrho,Af,curv)
	{

/******** LOOP 1 - Piston points: Calculate density and pressure. **********/
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

		/*Gravity Vector*/
		#if SIMDIM == 3
			StateVecD const g(0.0,-9.81,0.0);
		#else
			StateVecD const g(0.0,-9.81);
			// StateVecD const g(0.0,0.0);
		#endif

/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/		
		#pragma omp for schedule(static) nowait /*Reduction defs in Var.h*/
		for (size_t ii=0; ii < start; ++ii)
		{
			Part const& pi(pnp1[ii]);
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];
			real Rrhoi = 0.0;
			real Rrhod = 0.0;

			for (Part const& pj : neighb[ii])
			{	/* Neighbour list loop. */

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				StateVecD const Rij = pj.xi-pi.xi;
				StateVecD const Vij = pj.v-pi.v;
				real const rr = pj.nDist;
				real const r = sqrt(rr);
				real const volj = pj.m/pj.rho;
				StateVecD const gradK = /*dp.L[ii] * */GradK(Rij, r,fvar.H, fvar.correc);

				// StateVecD contrib = BasePos(pi.p,pj.p,volj,gradK);
				// StateVecD contrib = BasePos(pi,pj,gradK);

				/*drho/dt*/
				Rrhoi -= pj.m*(Vij.dot(gradK));
				Rrhod -= Continuity_dSPH(Rij,rr,fvar.HSQ,gradK,volj,dp.gradRho[ii],dp.gradRho[pj.partID],pi,pj);
				// RV[ii] -= pj.m*contrib;
			}/*End of neighbours*/
			

			RV[ii] = StateVecD::Zero();
			Af[ii] = StateVecD::Zero();
			curv[ii] = 0.0;
			Rrho[ii] = (Rrhoi - fvar.dCont * Rrhod);

		} /*End of boundary parts*/


/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Force) schedule(static) nowait
		for (size_t ii = start; ii < end; ++ii)
		{
			Part const pi(pnp1[ii]);
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
			StateVecD viscI = StateVecD::Zero();


			if(dp.lam[ii] < 0.75 && pi.b == PartState.FREE_)
			{
				if (svar.Asource == 1)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - /*pi.v*/ dp.avgV[ii];
					Pbasei = pi.cellP - avar.pRef;
				}
				else if (svar.Asource == 2)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - /*pi.v*/ dp.avgV[ii];
				}
				else 
				{
					Vdiff = avar.vInf - pi.v /*dp.avgV[ii]*/;
				}

#if SIMDIM == 3
				if(svar.Asource == 3)
				{	
					StateVecD Vel = svar.vortex.getVelocity(pi.xi);
					Vdiff = Vel - /*pi.v*/ dp.avgV[ii];
				}
#endif
				// cout << ii << endl;
				aero = CalcAeroForce(avar,pi,Vdiff,dp.norm[ii],dp.lam[ii],Pbasei);
				
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

				// if (aero != aero)
				// {
				// 	cout << "Aerodynamic force has gone NaN" << endl;
				// 	cout << "Point: " << ii - start << "  normal: " << dp.norm[ii](0) << "  " << dp.norm[ii](1) << endl;
				// 	cout << "Vdiff: " << Vdiff(0) << "  " << Vdiff(1) << endl;
				// 	exit(-1);
				// }

				// cout << aero(0) << "  " << aero(1) << endl;
				// curv[ii] = (std::min(avar.nfull,real(outlist[ii].size()-1))/(avar.nfull));
			}
			else 
			{
				curv[ii] = 0.0;
			}

			for (Part const& pj : neighb[ii])
			{	/* Neighbour list loop. */

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
				{
					// if(pi.surf == 1 && pi.b == PartState.FREE_)
					// {
					// 	RVi += aero *  pj.m/pj.rho* fvar.correc/dp.kernsum[ii];
					// 	Force += aero*  pj.m/pj.rho* fvar.correc/dp.kernsum[ii];
					// }
					
					continue;
				}

				StateVecD const Rij = pj.xi-pi.xi;
				StateVecD const Vij = pj.v-pi.v;
				real const rr = pj.nDist;
				real const r = sqrt(rr);	
				real const volj = pj.m/pj.rho;
				// StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rij,r,fvar.H,fvar.correc);/* gradK;*/

				if( pi.surf == 1 /* dp.lam[ii] < 0.75  */)
				{
					// if(/* svar.iter % 10 == 0  &&*/ r > 1e-6*fvar.H)
					// 	gradLK = dp.L[ii]*gradK;

					if( pj.surf == 1 /* dp.lam[pj.partID] < 0.75  */)
					{	
						curve -= (dp.norm[pj.partID].normalized()-dp.norm[ii].normalized()).dot(volj*gradK);
						correc += volj * Kernel(r,fvar.H,fvar.correc)/*/dp.kernsum[ii]*/;
					}	
				}

				// if(pi.b == PartState.FREE_ && pi.surf == 1)
				// {
				// 	RV[pj.partID] += aero * pj.m/pj.rho* Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				// 	Force += aero * pj.m/pj.rho* Kernel(r,fvar.H,fvar.correc)/dp.kernsum[ii];
				// }

				
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

			if(/* dp.lam[ii] < 0.75 */pi.surf == 1)
			{
				RVi += (fvar.sig/pi.rho * curve * dp.norm[ii]/* .normalized() */)/correc;
			}
			
			// #pragma omp critical
			// Force += aero;	

			
			Af[ii] = aero;
			RV[ii] = (RVi + artViscI + viscI + aero/pi.m + g);
			Rrho[ii] = (Rrhoi - fvar.dCont * Rrhod);

			// res_.emplace_back(aero/pi.m);
			// Rrho_.emplace_back(0.0);

		} /*End of sim parts*/		

	}	/*End of declare parallel */
 
    // cout << RV.size() << "  " << Rrho.size() << "  " << Af.size() << "  " << curv.size() << endl;
}

#endif