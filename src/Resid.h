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


ldouble GetNumpartdens(SIM &svar, FLUID &fvar, State &pnp1, outl &outlist)
{
	ldouble npd = 0.0;
	#pragma omp parallel for reduction(+:npd)
	for (uint i=0; i< svar.totPts; ++i)
	{
		StateVecD pi = pnp1[i].xi;
		for (uint j=0; j<outlist[i].size(); ++j)
		{ /* Surface Tension calcs */
			StateVecD pj = pnp1[outlist[i][j]].xi;
			ldouble r = (pj-pi).norm();
			npd += W2Kernel(r,fvar.H,fvar.correc);
		}
	}
	return npd/ldouble(svar.totPts);
}

/* Colour field gradient in Hu et al (2014) method, note the bottom variable to account for no air*/
std::vector<StateVecD> GetColourGrad(SIM &svar, FLUID &fvar, State &pnp1, outl &outlist)
{
	std::vector<StateVecD> cgrad(svar.totPts, StateVecD::Zero());

	#pragma omp parallel for
	for(uint i=svar.bndPts; i < svar.totPts; ++i)
	{
		StateVecD pi = pnp1[i].xi;
		// StateVecD top = StateVecD::Zero();
		ldouble bottom = 0.0;
		// pnp1[i].normal = StateVecD::Zero();
		
		for(auto j:outlist[i])
		{	/*Find the denominator to correct absence of second phase*/
			Part pj = pnp1[j];
			ldouble r = (pj.xi-pi).norm();
			bottom +=(pj.m/pj.rho)*W2Kernel(r,fvar.H,fvar.correc);

		}

		for(auto j:outlist[i])
		{	/*Find the numerator and sum*/
			Part pj = pnp1[j];
			StateVecD Rij = pj.xi-pi;
			ldouble r = Rij.norm();
			cgrad[i] += (pj.m/pj.rho)*W2GradK(Rij,r,fvar.H,fvar.correc);
			// pnp1[i].normal += Rij;
		}
		
		cgrad[i] = cgrad[i]/bottom;

	}

	return cgrad;
}

StateVecD Base(FLUID &fvar, Part &pi, Part &pj, 
	StateVecD &Rij, StateVecD &Vij, ldouble &r, StateVecD &Grad, std::vector<ldouble> &mu)
{
	const static ldouble alpha = fvar.alpha; 	/* Artificial Viscosity Parameter*/

	/*Pressure and artificial viscosity - Monaghan (1994) p.400*/
	ldouble rhoij = 0.5*(pi.rho+pj.rho);
	ldouble cbar= 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));
	ldouble vdotr = Vij.dot(Rij);
	ldouble muij= fvar.H*vdotr/(r*r+0.01*fvar.HSQ);
	mu.emplace_back(muij);
	ldouble pifac = alpha*cbar*muij/rhoij;

	if (vdotr > 0) pifac = 0;
	return Grad*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));
}

/*Laminar Viscosity - Morris (2003)*/
StateVecD Viscosity(FLUID &fvar, Part &pi, Part &pj, 
	StateVecD &Rij, StateVecD &Vij, ldouble &r, StateVecD &Grad)
{
	return -Vij*(fvar.mu)/(pi.rho*pj.rho)*(1.0/(r*r+0.01*fvar.HSQ))*Rij.dot(Grad);
}

/*Surface Tension - Nair & Poeschel (2017)*/
StateVecD SurfaceTens(FLUID &fvar, Part &pj, StateVecD &Rij, ldouble &r, ldouble &npd)
{
	/*Surface tension factor*/
	const static ldouble lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));

	ldouble fac=1.0; /*Boundary Correction Factor*/
	if(pj.b==0) fac=(1+0.5*cos(M_PI*(fvar.contangb/180))); 

	/*npd = numerical particle density (see code above) */
	ldouble sij = 0.5*pow(npd,-2.0)*(fvar.sig/lam)*fac;
	return -(Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*fvar.H));
}

StateVecD HuST(FLUID &fvar, Part &pi, Part &pj, StateVecD &Rij, ldouble &r, StateVecD &cgradi, StateVecD &cgradj)
{

	StateVecD Fst;

	Fst = (0.01/2.0)*(pi.m/pi.rho)*(pj.m/pj.rho)*((cgradi.squaredNorm()+cgradj.squaredNorm())/2.0)
		*W2GradK(Rij,r,fvar.H,fvar.correc);

	return Fst;
}

///**************** RESID calculation **************
void Forces(SIM &svar, FLUID &fvar, CROSS &cvar, State &pnp1, outl &outlist)
{
	svar.maxmu=0; 					/* CFL Parameter */
	StateVecD g = StateVecD::Zero();	/*Gravity Vector*/
	#if SIMDIM == 3
		g(2) = -9.81;
	#else
		g(1) = -9.81;
	#endif
	
	/********* LOOP 1 - all points: Calculate numpartdens ************/
	// ldouble numpartdens = GetNumpartdens(svar, fvar, pnp1, outlist);
	// std::vector<StateVecD> cgrad(svar.totPts,StateVecD::Zero());
	// cgrad = GetColourGrad(svar,fvar,pnp1,outlist);

	std::vector<ldouble> Rrhocontr(svar.totPts,0.0);
	std::vector<StateVecD> RV(svar.totPts,StateVecD::Zero()); /*Residual force*/
	std::vector<StateVecD> ST(svar.totPts,StateVecD::Zero()); /*Surface tension force*/		

	#pragma omp parallel shared(g)
	{
/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/
		#pragma omp for reduction(+:Rrhocontr)
		for (uint i=0; i< svar.bndPts; ++i)
		{
			
			Part pi = pnp1[i];

			for(auto j:outlist[i])
			{
				Part pj = pnp1[j];
				StateVecD Rij = pj.xi-pi.xi;
				StateVecD Vij = pj.v-pi.v;
				ldouble r = Rij.norm();
				StateVecD Grad = W2GradK(Rij, r,fvar.H,fvar.correc);
				Rrhocontr[i] -= pj.m*(Vij.dot(Grad));
			}
			pnp1[i].Rrho = Rrhocontr[i]; /*drho/dt*/
		}

/******* LOOP 3 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Rrhocontr, RV, ST) /*Reduction defs in Aero.h*/
		for (uint i=svar.bndPts; i < svar.totPts; ++i)
		{
			Part pi = pnp1[i];
			uint size = outlist[i].size();
			pnp1[i].theta = double(size);

			std::vector<double> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.reserve(size+1);
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/
			
			for (auto j:outlist[i])
			{	/* Neighbour list loop. */
				StateVecD contrib= StateVecD::Zero();
				StateVecD visc = StateVecD::Zero();
				StateVecD SurfC= StateVecD::Zero();

				Part pj(pnp1[j]);

				/*Check if the position is the same, and skip the particle if yes*/
				if(i == j)
					continue;

				StateVecD Rij = pj.xi-pi.xi;
				StateVecD Vij = pj.v-pi.v;
				ldouble r = Rij.norm();
				StateVecD Grad = W2GradK(Rij, r,fvar.H, fvar.correc);

				contrib = Base(fvar,pi,pj,Rij,Vij,r,Grad,mu);

				/*Laminar Viscosity - Morris (2003)*/
				visc    = Viscosity(fvar,pi,pj,Rij,Vij,r,Grad);

				/*Surface Tension - Nair & Poeschel (2017)*/
				// SurfC   = SurfaceTens(fvar,pj,Rij,r,numpartdens);
				// SurfC = HuST(fvar,pi,pj,Rij,r,cgrad[i],cgrad[j]);

				/*drho/dt*/
				Rrhocontr[i] -= pj.m*(Vij.dot(Grad));

				RV[i] += pj.m*contrib + pj.m*visc + SurfC/pj.m;
				ST[i] += SurfC/pj.m;
			}
	
			//CFL f_cv Calc
			ldouble it = *max_element(mu.begin(),mu.end());
			if (it > svar.maxmu)
				svar.maxmu=it;
		} /*End of sim parts*/

		//Update the actual structure
		#pragma omp for
		for(uint i=0; i < svar.bndPts; ++i)
		{
			pnp1[i].Rrho = Rrhocontr[i]; /*drho/dt*/			
		}

		#pragma omp for
		for(uint i=svar.bndPts; i < svar.totPts; ++i)
		{
			pnp1[i].f = RV[i] + g;
			pnp1[i].Rrho = Rrhocontr[i]; /*drho/dt*/
			pnp1[i].Sf = ST[i];			
		}

	}
		/*Aerodynamic force*/
	if(svar.Bcase >= 3)
		ApplyAero(svar,fvar,cvar,pnp1,outlist);
	
}

// ///*Density Reinitialisation using Least Moving Squares as in A. Colagrossi (2003)*
// void DensityReinit(FLUID &fvar, State &pnp1, outl &outlist)
// {
// 	DensVecD one = DensVecD::Zero();
//   	one[0] = 1.0;
  	
// 	#pragma omp parallel for
// 	for(uint i=0; i< outlist.size(); ++i)
// 	{
// 		DensMatD A= DensMatD::Zero();
// 		//Find matrix A.
// 		Part pi = pnp1[i];
// 		for (auto j:outlist[i])
// 		{
// 			Particle pj = pnp1[j];
// 			StateVecD Rij = pi.xi-pj.xi;
// 			DensMatD Abar = DensMatD::Zero();
// 			// Abar << 1   , Rij(0)        , Rij(1)        , Rij(2)        ,
// 			// 	    Rij(0) , Rij(0)*Rij(0) , Rij(1)*Rij(0) , Rij(2)*Rij(0) ,
// 			// 	    Rij(1) , Rij(0)*Rij(1) , Rij(1)*Rij(1) , Rij(2)*Rij(1) ,
// 			//      Rij(2) , Rij(0)*Rij(2) , Rij(1)*Rij(2) , Rij(2)*Rij(2) ;

// 	        Abar(0,0) = 1;
// 		    for (int ii = 0; ii < Rij.cols(); ++ii)
// 		    {
// 		        Abar(ii+1,0) = Rij[ii];
// 		        Abar(0,ii+1) = Rij[ii];
// 		        for (int jj = 0; jj<=ii; ++jj)
// 		        {
// 		            Abar(ii+1,jj+1) = Rij[ii]*Rij[jj];
// 		            Abar(jj+1,ii+1) = Rij[ii]*Rij[jj];
// 		        }
// 		    }

// 			A+= W2Kernel(Rij.norm(),fvar.H,fvar.correc)*Abar*pj.m/pj.rho;
// 		}

// 		DensVecD Beta;
// 		//Check if A is invertible
// 		Eigen::FullPivLU<DensMatD> lu(A);
// 		if (lu.isInvertible())
// 			Beta = lu.inverse()*one;
// 		else
// 			Beta = (1)*one;

// 		//Find corrected kernel
// 		ldouble rho = 0.0;
// 		for (uint j=0; j< outlist[i].size(); ++j)
// 		{
// 			StateVecD Rij = pi.xi-pnp1[outlist[i][j]].xi;
// 			rho += pnp1[outlist[i][j]].m*W2Kernel(Rij.norm(),fvar.H,fvar.correc)*
// 			(Beta(0)+Beta(1)*Rij(0)+Beta(2)*Rij(1));
// 		}

// 		pnp1[i].rho = rho;
// 	}

// }

#endif