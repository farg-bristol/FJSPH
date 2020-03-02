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

StateVecD AeroForce(const StateVecD& Vdiff, const SIM& svar, const FLUID& fvar, const AERO& avar, const real m)
{
	StateVecD Fd = StateVecD::Zero();
	// real Re = Vdiff.norm()*svar.Pstep*fvar.rhog/fvar.mug;
	const real Re = 2.0*fvar.rhog*Vdiff.norm()*avar.L/fvar.mug;
	real Cd = 0.0;

	Cd = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));


	// Fd = Cd*Vdiff*Vdiff.norm()*fvar.rhog/(fvar.rho0*svar.Pstep);
	const real Ai = M_PI*avar.L;
	Fd = 0.5*fvar.rhog*Vdiff.norm()*Vdiff*Cd*Ai/m;
	//Fd[1] = 1.0*Fd[1];
	//std::cout << "Reynolds: " << Re  << " Cd: " << Cd << " Fd: " << Fd[0] << " " << Fd[1] << std::endl;
	return Fd;
}


StateVecD CalcForce(SIM& svar, const FLUID& fvar, const AERO& avar, 
	const Part &pi, const StateVecD& Vdiff, const StateVecD& SurfC, const uint size, const real woccl)
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
			real correc = 1.0;
			if(size > 0.1 * avar.nfull)
			{
				real num = SurfC.dot(Vdiff);
				real denom = SurfC.norm()*(Vdiff).norm();
				real theta = num/denom;
				
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
			real correc = 1.0;
			real Acorrect = 0.0;

			/*Correction based on surface normal*/
			real num = SurfC.dot(Vdiff);
			real denom = SurfC.norm()*(Vdiff).norm();
			real theta = num/denom;
			
			if (theta <= 1.0  && theta >= -1.0)
			{
				correc = avar.a*W2Kernel(1-theta,avar.h1,1)+
					avar.b*W2Kernel(1-theta,avar.h2,1);
			}
			
			/*Correct based on the number of neighbours*/
			Acorrect = 0.9995*exp(-0.04606*real(size))+0.0005;
			// cout << size << "  " << Acorrect << endl;

			Fd = correc*Acorrect*Fd;
		}
	}
	else if(avar.acase == 4)
	{	/*Gissler et al (2017)*/
		const real nfull = avar.nfull;
		if(size < nfull)
		{	
			const real ymax = Vdiff.squaredNorm()*avar.ycoef;
			const real Re = 2.0*fvar.rhog*Vdiff.norm()*avar.L/fvar.mug;
			real Cds;

			const real frac2 = std::min(nfull,real(size))
							/(nfull);
			const real frac1 = (1.0 - frac2);

			// if (Re < 3500)
			 	Cds = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));
			// else
			// 	Cds = 1.699e-05*pow(Re,1.92)*(24.0/(Re));
			
			// if( Re <= 1000.0)
			// 	Cds = (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
			// else 
			// 	Cds = 0.424;

			const real Cdl = Cds*(1+2.632*ymax);
			const real	Cdi = frac1*Cdl + /*0.8**/frac2;

			#if SIMDIM == 3 
				const real Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
				const real Aunocc = frac1*Adrop + frac2*fvar.HSQ;
			#endif
			#if SIMDIM == 2
				const real Adrop = M_PI*(avar.L + avar.Cb*avar.L*ymax);
				const real Aunocc = frac1*Adrop + frac2*fvar.H;
			#endif


			const real Ai = (1-woccl)/*correc*/*Aunocc;

			Fd = 0.5*fvar.rhog*Vdiff.norm()*Vdiff*Cdi*Ai/pi.m;
		}
	}

	return Fd;
}

void ApplyAero(SIM &svar, const FLUID &fvar, const AERO &avar, 
	const State &pnp1, const outl &outlist,/* const vector<StateVecD>& vPert,*/
	 vector<StateVecD>& res, std::vector<StateVecD>& Af)
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
			real kernsum = 0.0;
			real woccl = 0.0;
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


			// Vdiff += pi.vPert;
			// cout << "Relative velocity: " << Vdiff(0) << "  " << Vdiff(1) << "  " << Vdiff(2) << endl;
			
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
				real r = Rij.norm();

				kernsum += W2Kernel(r,fvar.H,fvar.correc);

				/*Occlusion for Gissler Aero Method*/
				if (svar.Bcase > 2 && (avar.acase == 4 || avar.acase == 1))
				{
					real frac = -Vdiff.normalized().dot(Rij.normalized());
					
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
				
				real r = (pj.xi-pi.xi).norm();
				real kern = W2Kernel(r,fvar.H,fvar.correc);
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