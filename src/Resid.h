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

ldouble GetNumpartdens(const SIM &svar, const FLUID &fvar, const State &pnp1, const outl &outlist)
{
	ldouble npd = 0.0;
	const uint end = svar.totPts;
	#pragma omp parallel for reduction(+:npd) shared(outlist)
	for (uint ii=0; ii< end; ++ii)
	{
		StateVecD pi = pnp1[ii].xi;
		for (auto jj:outlist[ii])
		{ /* Surface Tension calcs */
			StateVecD pj = pnp1[jj].xi;
			ldouble r = (pj-pi).norm();
			npd += W2Kernel(r,fvar.H,fvar.correc);
		}
	}
	return npd/ldouble(svar.totPts);
}

/* Colour field gradient in Hu et al (2014) method, note the bottom variable to account for no air*/
std::vector<StateVecD> GetColourGrad(const SIM &svar, const FLUID &fvar, const State &pnp1, const outl &outlist)
{
	std::vector<StateVecD> cgrad(svar.totPts, StateVecD::Zero());
	const uint start = svar.bndPts;
	const uint end = svar.totPts;

	#pragma omp parallel for shared(outlist)
	for(uint ii=start; ii < end; ++ii)
	{
		StateVecD pi = pnp1[ii].xi;
		// StateVecD top = StateVecD::Zero();
		ldouble bottom = 0.0;
		// pnp1[ii].normal = StateVecD::Zero();
		
		for(auto jj:outlist[ii])
		{	/*Find the denominator to correct absence of second phase*/
			Part pj = pnp1[jj];
			ldouble r = (pj.xi-pi).norm();
			bottom +=(pj.m/pj.rho)*W2Kernel(r,fvar.H,fvar.correc);

		}

		for(auto jj:outlist[ii])
		{	/*Find the numerator and sum*/
			Part pj = pnp1[jj];
			StateVecD Rij = pj.xi-pi;
			ldouble r = Rij.norm();
			cgrad[ii] += (pj.m/pj.rho)*W2GradK(Rij,r,fvar.H,fvar.correc);
			// pnp1[ii].normal += Rij;
		}
		
		cgrad[ii] = cgrad[ii]/bottom;

	}

	return cgrad;
}

StateVecD Base(const FLUID &fvar, const Part &pi, const Part &pj, 
	const StateVecD &Rij, const StateVecD &Vij, const ldouble r, 
	const StateVecD &Grad, std::vector<ldouble> &mu)
{
	/*Pressure and artificial viscosity - Monaghan (1994) p.400*/
	ldouble pifac;
	ldouble vdotr = Vij.dot(Rij);
	ldouble muij= fvar.H*vdotr/(r*r+0.01*fvar.HSQ);
	if (vdotr > 0.0 || pj.b == 3) 
	{
		pifac = 0.0;
	}
	else
	{
		ldouble rhoij = 0.5*(pi.rho+pj.rho);
		ldouble cbar= 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));
		pifac = fvar.alpha*cbar*muij/rhoij;
	}
	
	mu.emplace_back(muij);
	return Grad*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));
}

/*Laminar Viscosity - Morris (2003)*/
StateVecD Viscosity(const FLUID &fvar, const Part &pi, const Part &pj, 
	const StateVecD &Rij, const StateVecD &Vij, const ldouble &r, const StateVecD &Grad)
{
	return -Vij*((fvar.mu)/(pi.rho*pj.rho))*(1.0/(r*r+0.01*fvar.HSQ))*Rij.dot(Grad);
}

/*Surface Tension - Nair & Poeschel (2017)*/
StateVecD SurfaceTens(const FLUID &fvar, const Part &pj, const StateVecD &Rij, 
					  const ldouble &r, const ldouble &npd)
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

StateVecD HuST(const FLUID &fvar, const Part &pi, const Part &pj, 
			  const StateVecD &Rij, const ldouble &r, const StateVecD &cgradi, const StateVecD &cgradj)
{

	StateVecD Fst;

	Fst = (0.01/2.0)*(pi.m/pi.rho)*(pj.m/pj.rho)*((cgradi.squaredNorm()+cgradj.squaredNorm())/2.0)
		*W2GradK(Rij,r,fvar.H,fvar.correc);

	return Fst;
}

///**************** RESID calculation **************
void Forces(SIM& svar, const FLUID& fvar, const AERO& avar, const State& pnp1/*, State& airP*/,
	 const vector<vector<Part>>& neighb, const outl& outlist,
	 vector<StateVecD>& RV, vector<ldouble>& Rrho)
{
	svar.maxmu=0; 					    /* CFL Parameter */
	const uint start = svar.bndPts;
	const uint end = svar.totPts;

	/*Gravity Vector*/
	#if SIMDIM == 3
		const static StateVecD g = StateVecD(0.0,0.0,-9.81);
	#else
		const static StateVecD g = StateVecD(0.0,-9.81);
	#endif
	// const uint piston = svar.psnPts;

	/********* LOOP 1 - all points: Calculate numpartdens ************/
	// ldouble numpartdens = GetNumpartdens(svar, fvar, pnp1, outlist);
	// std::vector<StateVecD> cgrad(svar.totPts,StateVecD::Zero());
	// cgrad = GetColourGrad(svar,fvar,pnp1,outlist);

	// std::vector<StateVecD> ST(svar.totPts,StateVecD::Zero()); /*Surface tension force*/		
	#pragma omp parallel /*shared(numpartdens)*/
	{

/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/
		// #pragma omp for reduction(+:Rrhocontr)
		// for (uint ii=0; ii < start; ++ii)
		// {
		// 	Part pi = pnp1[ii];

		// 	for(auto j:outlist[ii])
		// 	{
		// 		Part pj = pnp1[j];
		// 		StateVecD Rij = pj.xi-pi.xi;
		// 		StateVecD Vij = pj.v-pi.v;
		// 		ldouble r = Rij.norm();
		// 		StateVecD Grad = W2GradK(Rij, r,fvar.H,fvar.correc);
		// 		Rrhocontr[ii] -= pj.m*(Vij.dot(Grad));
		// 	}
		// 	pnp1[ii].Rrho = Rrhocontr[ii]; /*drho/dt*/
		// }

/******** LOOP 3 - Piston points: Calculate density and pressure. **********/		
		#pragma omp for reduction(+:RV, Rrho) /*Reduction defs in Aero.h*/
		for (uint ii=0; ii < start; ++ii)
		{
			const Part pi = pnp1[ii];
			const uint size = outlist[ii].size();
			// pnp1[ii].theta = double(size);
			// pnp1[ii].theta = 0.0;

			std::vector<double> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.reserve(size+1);			
			if (ii > neighb.size())
			{
				cout << "neighb array smaller than totPts" << endl;
				exit(-1);
			}
			for (Part const& pj:neighb[ii])
			{	/* Neighbour list loop. */

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				const StateVecD Rij = pj.xi-pi.xi;
				const StateVecD Vij = pj.v-pi.v;
				const ldouble r = Rij.norm();
				const StateVecD Grad = W2GradK(Rij, r,fvar.H, fvar.correc);

				StateVecD contrib = Base(fvar,pi,pj,Rij,Vij,r,Grad,mu);
				/*drho/dt*/
				Rrho[ii] -= pj.m*(Vij.dot(Grad));
				RV[ii] += pj.m*contrib;
			}/*End of neighbours*/
			
		} /*End of boundary parts*/

		

/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		#pragma omp for reduction(+:Rrho, RV/*, ST*/)  /*Reduction defs in Aero.h*/
		for (uint ii = start; ii < end; ++ii)
		{
			Part pi = pnp1[ii];
			uint size = outlist[ii].size();
			// pnp1[ii].theta = double(size);
			// pnp1[ii].theta = 0.0;

			vector<double> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.reserve(size+1);
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/		

			if (ii > neighb.size())
			{
				cout << "neighb array smaller than totPts" << endl;
				exit(-1);
			}
			for (Part const& pj:neighb[ii])
			{	/* Neighbour list loop. */

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.partID == pj.partID)
					continue;

				const StateVecD Rij = pj.xi-pi.xi;
				const StateVecD Vij = pj.v-pi.v;
				const ldouble r = Rij.norm();
				const StateVecD Grad = W2GradK(Rij, r,fvar.H, fvar.correc);

				/*Momentum contribution - Monaghan (1994)*/
				const StateVecD contrib = Base(fvar,pi,pj,Rij,Vij,r,Grad,mu);
				
				/*Laminar Viscosity - Morris (2003)*/
				const StateVecD visc    = Viscosity(fvar,pi,pj,Rij,Vij,r,Grad);

				if (pj.b != 3)
				{
					/*Surface Tension - Nair & Poeschel (2017)*/
					// StateVecD SurfC   = SurfaceTens(fvar,pj,Rij,r,numpartdens);
					// SurfC = HuST(fvar,pi,pj,Rij,r,cgrad[ii],cgrad[jj]);
				
					/*drho/dt*/
					Rrho[ii] -= pj.m*(Vij.dot(Grad));
				}

				RV[ii] += pj.m*contrib + pj.m*visc /*+ SurfC/pj.m*/;

				// if (pj.b == 3)
				// {
				// 	ST[ii] += pj.m*contrib + pj.m*visc;
				// }

				// ST[ii] += SurfC/pj.m;
			}/*End of neighbours*/
			
			RV[ii] += g;
			//CFL f_cv Calc
			double it = *max_element(mu.begin(),mu.end());
			if(it > svar.maxmu)
				svar.maxmu = it;
		} /*End of sim parts*/

		// cout << endl;
		//Update the actual structure
		// StateVecD Force = StateVecD::Zero();
		// #pragma omp for
		// for (uint ii = 0; ii < piston; ++ii)
		// {
		// 	Force += RV[ii];
		// }
		// svar.Force = Force;

		// #pragma omp for
		// for(uint ii=0; ii < end; ++ii)
		// {
		// 	pnp1[ii].f = RV[ii] + g;
		// 	pnp1[ii].Rrho = Rrho[ii]; /*drho/dt*/
		// 	// pnp1[ii].Sf = ST[ii];			
		// }
		
	}	/*End of declare parallel */
	

		/*Aerodynamic force*/
	if(svar.Bcase > 1)
		ApplyAero(svar,fvar,avar,pnp1,outlist,RV);
	
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