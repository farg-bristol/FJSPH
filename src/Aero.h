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

#ifndef AERO_H
#define AERO_H

#include "Var.h"
#include "Kernel.h"

#pragma omp declare reduction(+: std::vector<StateVecD> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
        	[](StateVecD lhs, StateVecD rhs){return lhs + rhs;})) \
                    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(+: std::vector<ldouble> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
        	[](ldouble lhs, ldouble rhs){return lhs + rhs;})) \
                    initializer(omp_priv = omp_orig)                  

#pragma omp declare reduction(+:StateVecD : omp_out=omp_out+omp_in)\
                    initializer(omp_priv = omp_orig) 

StateVecD AeroForce(const StateVecD& Vdiff, const SIM& svar, const FLUID& fvar, const AERO& avar, const ldouble m)
{
	StateVecD Fd = StateVecD::Zero();
	// ldouble Re = Vdiff.norm()*svar.Pstep*fvar.rhog/fvar.mug;
	const ldouble Re = 2.0*fvar.rhog*Vdiff.norm()*avar.L/fvar.mug;
	ldouble Cd = 0.0;

	Cd = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));


	// Fd = Cd*Vdiff*Vdiff.norm()*fvar.rhog/(fvar.rho0*svar.Pstep);
	const ldouble Ai = M_PI*avar.L;
	Fd = 0.5*fvar.rhog*Vdiff.norm()*Vdiff*Cd*Ai/m;
	//Fd[1] = 1.0*Fd[1];
	//std::cout << "Reynolds: " << Re  << " Cd: " << Cd << " Fd: " << Fd[0] << " " << Fd[1] << std::endl;
	return Fd;
}


StateVecD CalcForce(SIM& svar, const FLUID& fvar, const AERO& avar, 
	const Part &pi, const StateVecD& Vdiff, const StateVecD& SurfC, const uint size, const ldouble woccl)
{
	StateVecD Fd= StateVecD::Zero();

	if( avar.acase == 1)
	{
		 /* All upstream particles */
		Fd = woccl*avar.Acorrect*AeroForce(Vdiff, svar, fvar, avar, pi.m);
	}
	else if(avar.acase == 2)
	{ /*Surface particles, with correction based on surface normal*/
		// cout << size << endl;
		if (size < avar.nfull)
		{
			Fd = AeroForce(Vdiff, svar, fvar, avar, pi.m);
			/*Correction based on surface tension vector*/
			ldouble correc = 1.0;
			if(size > 0.1 * avar.nfull)
			{
				ldouble num = SurfC.dot(Vdiff);
				ldouble denom = SurfC.norm()*(Vdiff).norm();
				ldouble theta = num/denom;
				
				if (theta <= 1.0  && theta >= -1.0)
				{
					correc = avar.a*W2Kernel(1-theta,avar.h1,1)+
						avar.b*W2Kernel(1-theta,avar.h2,1);
				}
			}
			Fd = /*avar.Acorrect**/correc*Fd;
		}
		
	}
	else if(avar.acase == 3)	/*Case 2, plus a correction based on number of neighbours*/
	{
		if (size < avar.nfull)
		{
			// StateVecD Vdiff = avar.vInf-pi.v;

			Fd = AeroForce(Vdiff, svar, fvar, avar, pi.m);
			ldouble correc = 1.0;
			double Acorrect = 0.0;

			/*Correction based on surface normal*/
			ldouble num = SurfC.dot(Vdiff);
			ldouble denom = SurfC.norm()*(Vdiff).norm();
			ldouble theta = num/denom;
			
			if (theta <= 1.0  && theta >= -1.0)
			{
				correc = avar.a*W2Kernel(1-theta,avar.h1,1)+
					avar.b*W2Kernel(1-theta,avar.h2,1);
			}
			
			/*Correct based on the number of neighbours*/
			Acorrect = 0.9995*exp(-0.04606*double(size))+0.0005;
			// cout << size << "  " << Acorrect << endl;

			Fd = correc*Acorrect*Fd;
		}
	}
	else if(avar.acase == 4)
	{	/*Gissler et al (2017)*/
		// if(size < nfull)
		// {	
			const double nfull = avar.nfull;
			const ldouble ymax = Vdiff.squaredNorm()*avar.ycoef;
			const ldouble Re = 2.0*fvar.rhog*Vdiff.norm()*avar.L/fvar.mug;
			ldouble Cds;

			const ldouble frac2 = std::min(nfull,double(size))
							/(nfull);
			const ldouble frac1 = (1.0 - frac2);

			// if (Re < 3500)
			//  	Cds = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));
			// else
			// 	Cds = 1.699e-05*pow(Re,1.92)*(24.0/(Re));
			
			if( Re <= 1000.0)
				Cds = (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
			else 
				Cds = 0.424;

			const ldouble Cdl = Cds*(1+2.632*ymax);
			const ldouble	Cdi = frac1*Cdl + /*0.8**/frac2;

			#if SIMDIM == 3 
				const ldouble Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
				const ldouble Aunocc = frac1*Adrop + frac2*fvar.HSQ;
			#endif
			#if SIMDIM == 2
				const ldouble Adrop = M_PI*(avar.L + avar.Cb*avar.L*ymax);
				const ldouble Aunocc = frac1*Adrop + frac2*fvar.H;
			#endif


			const ldouble Ai = (1-woccl)/*correc*/*Aunocc;

			Fd = 0.5*fvar.rhog*Vdiff.norm()*Vdiff*Cdi*Ai/pi.m;
		// }
	}

	return Fd;
}

void ApplyAero(SIM &svar, const FLUID &fvar, const AERO &avar, 
	const State &pnp1, const outl &outlist, vector<StateVecD>& res, std::vector<StateVecD>& Af)
{
	// std::vector<StateVecD> Af(svar.totPts,StateVecD::Zero());
	const uint start = svar.bndPts;
	const uint end = svar.totPts; /*DO NOT MAKE STATIC*/
	
	#pragma omp parallel for reduction(+:Af) shared(svar)
	for (uint ii=start; ii < end; ++ii)
	{
		if(pnp1[ii].b == FREE)
		{	
			uint size = outlist[ii].size();
			Part pi(pnp1[ii]);
			pi.normal = StateVecD::Zero();
			ldouble kernsum = 0.0;
			ldouble woccl = 0.0;
			StateVecD Vdiff;
#if SIMDIM == 3
			if(svar.Bcase == 4)
			{	
				StateVecD Vel = svar.vortex.getVelocity(pi.xi);
				Vdiff = Vel- pi.v;
			}
#endif
			if (svar.Bcase == 6)
			{
				Vdiff = pi.cellV - pi.v;
			}
			else 
			{
				Vdiff = avar.vInf - pi.v;
			}


			for (auto const jj:outlist[ii])
			{	/* Neighbour list loop. */
				const Part pj(pnp1[jj]);

				if(pi.partID == pj.partID)
				{
					kernsum += fvar.correc;
					continue;
				}

				StateVecD Rij = pj.xi-pi.xi;
				pi.normal += Rij;
				ldouble r = Rij.norm();

				kernsum += W2Kernel(r,fvar.H,fvar.correc);

				/*Occlusion for Gissler Aero Method*/
				if (svar.Bcase > 2 && (avar.acase == 4 || avar.acase == 1))
				{
					ldouble frac = -Vdiff.normalized().dot(Rij.normalized());
					
					// cout << num/denom << endl;
					if (frac > woccl)
					{
						woccl = frac;
					}
				}
			} /*End of neighbours*/

			/*Find the aero force*/
			StateVecD Fd = CalcForce(svar,fvar,avar,pi,Vdiff,pi.normal,size,woccl);
			// pnp1[i].theta = woccl;
			// if(size > 1.0/3.0 *avar.nfull)
			// {
			for (auto jj:outlist[ii])
			{	/* Neighbour list loop. */
				const Part pj = pnp1[jj];

				if(pj.b == BOUND)
					continue;
				
				ldouble r = (pj.xi-pi.xi).norm();
				ldouble kern = W2Kernel(r,fvar.H,fvar.correc);
				Af[jj] += Fd*kern/kernsum;
				res[ii] += Fd*kern/kernsum;
			}
			// }
			// else 
			// {
			// 	pnp1[i].Af += Fd;
			// }
		}/*End of if*/

	}/*End of ii particles*/

	// StateVecD Force = StateVecD::Zero();
	// #pragma omp parallel for/* reduction(+:Force)*/
	// for (uint ii = start; ii < end; ++ii)
	// {
		
	// 	// pnp1[ii].Af = Af[ii];

	// 	// Force += Af[ii];
	// }
	// svar.Force = Force;
	

}

#endif