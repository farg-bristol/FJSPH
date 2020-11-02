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

StateVecD AeroForce(StateVecD const& Vdiff, AERO const& avar, real const mass)
{
	real const Re = 2.0*avar.rhog*Vdiff.norm()*avar.L/avar.mug;
	real const Cdi = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.00001));

	#if SIMDIM == 3
	real const Ai = M_PI*avar.L*avar.L;
	#else
	real const Ai = M_PI*avar.L;
	#endif
		//Fd[1] = 1.0*Fd[1];
	//std::cout << "Reynolds: " << Re  << " Cd: " << Cd << " Fd: " << Fd[0] << " " << Fd[1] << std::endl;
	return (0.5*avar.rhog*Vdiff.norm()*Vdiff*Cdi*Ai/mass);
}

/*Gissler et al (2017)*/
StateVecD GisslerForce(AERO const& avar, StateVecD const& Vdiff, 
				real const& mass, real const size, real const woccl)
{
	real const nfull = avar.nfull;

	real ymax = Vdiff.squaredNorm()*avar.ycoef;
	// ymax = 0.0;
	if (ymax > 1.0)
		ymax = 1.0;

	real const Re = 2.0*avar.rhog*Vdiff.norm()*avar.L/avar.mug;
	real Cds;

	real const frac2 = std::min(nfull,real(size-1))
					/(nfull);
	real const frac1 = (1.0 - frac2);


	// if (Re < 3500)
 	Cds = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.00001));

	// else
	// 	Cds = 1.699e-05*pow(Re,1.92)*(24.0/(Re));
	
	// if( Re <= 1000.0)
	// 	Cds = (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
	// else 
	// 	Cds = 0.424;

	real const Cdl = Cds*(1+2.632*ymax);
	real const	Cdi = frac1*Cdl + /*1.37**/frac2;

#if SIMDIM == 3 
	real const Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
	// real const Adrop = M_PI*pow(avar.L,2);
#else
	real const Adrop = 2*(avar.L + avar.Cb*avar.L*ymax) * pow(avar.L,1);
	// real const Adrop = 2*avar.L;
#endif

	real const Aunocc = (frac1*Adrop + frac2*avar.aPlate);

	real const Ai = (1-woccl)*Aunocc;

	// cout << "Fractions: " << frac1 << "  " << frac2 << endl;
	// cout << "Areas: " << Ai << "  " << Adrop << "  " << avar.aPlate << endl;
	// cout << "Cds: " << Cdi << "  "  << Cdl << "  " << Cds << endl << endl;;
	// cout << avar.ycoef << endl;

	return 0.5*avar.rhog*Vdiff.norm()*Vdiff*Cdi*Ai/mass;


}

StateVecD CalcAeroForce(AERO const& avar, Part const& pi, StateVecD const& Vdiff,
		StateVecD const& norm, uint const size, real const woccl, real const Pbasei)
{
	StateVecD Fd= StateVecD::Zero();

	// cout << avar.acase << endl;
	if( avar.acase == 1)
	{
		 /* All upstream particles */
		Fd = woccl*avar.Acorrect*AeroForce(Vdiff, avar, pi.m);
	}
	else if(avar.acase == 2)
	{ /*Surface particles, with correction based on surface normal*/
		// cout << size << endl;
		if (size < avar.nfull)
		{
			Fd = AeroForce(Vdiff, avar, pi.m);
			/*Correction based on surface tension vector*/
			real correc = 1.0;
			if(size > 0.1 * avar.nfull)
			{
				real num = norm.dot(Vdiff);
				real denom = norm.norm()*(Vdiff).norm();
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
	else if(avar.acase == 3)
	{

			// real ymax = Vdiff.squaredNorm()*avar.ycoef;

			// Fd = AeroForce(Vdiff, avar, pi.m);
			// real correc = 1.0;
			// real Acorrect = 1.0-woccl;

			// /*Correction based on surface normal*/
			real theta = norm.normalized().dot(Vdiff.normalized());
			// real denom = norm.norm()*(Vdiff).norm();
			// real theta = num/denom;

			real Cp = 0.0;
			if(abs(theta) < 2.4435)
			{
				Cp = 1.1*cos(theta*2.03)-0.1;
			}
			else
			{
				Cp = 0.075;
			}

			real const frac2 = real(size-1)/(avar.nfull);
			real const frac1 = (1.0 - frac2);


#if SIMDIM == 3 
			// real const Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
			real const Adrop = M_PI*pow(avar.L,2);
#endif
#if SIMDIM == 2
			// real const Adrop = 2*(avar.L + avar.Cb*avar.L*ymax);
			real const Adrop = 2*avar.L;
#endif

			real const Aunocc = frac1*Adrop + frac2*avar.aPlate;

// #if SIMDIM == 2
// 			real area = avar.L;
// #else
// 			real area = M_PI*avar.L*avar.L;
// #endif

			real press = 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;

#if SIMDIM == 3
			Fd = -7.5*norm.normalized()*Aunocc*press/pi.m;
#else
			Fd = -norm.normalized()*Aunocc*press/pi.m;
#endif

			// if (theta <= 1.0  && theta >= -1.0)
			// {
			// 	correc = avar.a*W2Kernel(1-theta,avar.h1,1)+
			// 		avar.b*W2Kernel(1-theta,avar.h2,1);
			// }
			
			// /*Correct based on the number of neighbours*/
			// Acorrect = 0.9995*exp(-0.04606*real(size))+0.0005;
			// cout << size << "  " << Acorrect << endl;
			

			// Fd = correc*Acorrect*Fd;
		
	}
	else if(avar.acase == 4)
	{	
		Fd = GisslerForce(avar,Vdiff,pi.m,size,woccl);
	}
	else if(avar.acase == 5)
	{
							
		real theta = acos(-norm.normalized().dot(Vdiff.normalized()));
		
		real Cp = 0.0;

		// if(abs(theta) < 2.48073)
		// {
		// 	Cp = 1- 4*pow(sin(abs(theta)),3);
		// }
		// else
		// {
		// 	Cp = 0.075;
		// }

		if(abs(theta) < 2.6399)
		{
			Cp = 1.0 - 4.0*pow(sin(abs(theta)),2.0);
		}
		else
		{
			Cp = 0.075;
		}

		real Plocali = 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;
		// cout << Plocalj << endl;

		real Pi = (Plocali+Pbasei);

		// real sph_d = sqrt(pi.m/pi.rho);

		real const Re = avar.rhog*Vdiff.norm()*avar.L/avar.mug;
		real const Cdi = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));


		// cout << Cdi << endl;
		StateVecD F_drop = 0.5*avar.rhog*Vdiff.norm()*Vdiff*(M_PI*avar.L*avar.L/4)*Cdi/pi.m;
		// aeroD = -Pi * avar.aPlate * norm[ii].normalized();
		StateVecD F_kern = -1500.0*(Pi/(pi.m*pi.rho)) * norm;

		real const frac1 = std::min(real(size-1),avar.nfull)/(avar.nfull);

		// cout << F_kern.norm() << "  " << F_drop.norm() << endl;
		// cout << avar.nfull << "  " << frac1  << "  " << real(size-1) << endl;

		Fd = frac1*F_kern + (1-frac1)*F_drop;

	}

	return Fd;
}

void ApplyAero(SIM & svar, FLUID const& fvar, AERO const& avar, MESH const& cells,
	State const& pnp1, outl const& outlist, vector<StateVecD>& res, std::vector<StateVecD>& Af)
{
	// std::vector<StateVecD> Af(svar.totPts,StateVecD::Zero());
	size_t const start = svar.bndPts;
	size_t const end = svar.totPts; 
	StateVecD Force = StateVecD::Zero();

	#pragma omp parallel for shared(svar) /*reduction(+:res,Af)*/ 
	for (size_t ii = start; ii < end; ++ii)
	{
		if(pnp1[ii].b == PartState.FREE_)
		{	
			size_t size = outlist[ii].size();
			Part pi(pnp1[ii]);
			pi.normal = StateVecD::Zero();
			// real kernsum = 0.0;
			real woccl = 0.0;
			StateVecD Vdiff;

			if (svar.Asource == 1)
			{
				Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - pi.v;
			}
			else if (svar.Asource == 2)
			{
				Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - pi.v;
			}
			else 
			{
				Vdiff = avar.vInf - pi.v;
			}

#if SIMDIM == 3
			if(svar.Asource == 3)
			{	
				StateVecD Vel = svar.vortex.getVelocity(pi.xi);
				Vdiff = Vel- pi.v;
			}
#endif

			// Vdiff += pi.vPert;
			// cout << "Relative velocity: " << Vdiff(0) << "  " << Vdiff(1) << "  " << Vdiff(2) << endl;
			
			for (auto const& jj:outlist[ii])
			{	/* Neighbour list loop. */
				Part const& pj(pnp1[jj]);

				if(pi.partID == pj.partID)
				{
					// kernsum += fvar.correc;
					continue;
				}

				StateVecD const Rij = pj.xi-pi.xi;
				// real const r = Rij.norm();

				// kernsum += W2Kernel(r,fvar.H,fvar.correc);

				/*Occlusion for Gissler Aero Method*/
				if (avar.acase == 4 || avar.acase == 1)
				{
					real const frac = -Vdiff.normalized().dot(Rij.normalized());
					
					// cout << num/denom << endl;
					if (frac > woccl)
					{
						woccl = frac;
					}
				}
			} /*End of neighbours*/

			/*Find the aero force*/
			StateVecD Fd = CalcAeroForce(avar,pi,Vdiff,pi.normal,size,woccl,0.0);
			// StateVecD Fd = GisslerForce(svar,fvar,avar,Vdiff,pi.m,size,woccl);

			// cout << Fd(0) << "  " << Fd(1) << endl;

			// if(Fd.norm() > 0.0)
			// cout << Fd.norm() << endl;
			// pnp1[i].theta = woccl;
			// if(size > 1.0/3.0 *avar.nfull)
			// {
			// for (auto const& jj:outlist[ii])
			// {	/* Neighbour list loop. */
			// 	if(pnp1[jj].b == BOUND)
			// 		continue;
				
			// 	real const r = (pnp1[jj].xi-pi.xi).norm();
			// 	real const kern = W2Kernel(r,fvar.H,fvar.correc);
			// 	Af[jj] += Fd*kern/kernsum;
			// 	res[jj] += Fd*kern/kernsum;
			// }

			Af[ii] += Fd;
			res[ii] += Fd;
			Force += Fd/**pnp1[ii].m*/;
			// }
			// else 
			// {
			// 	pnp1[i].Af += Fd;
			// }



		}/*End of if*/

	}/*End of ii particles*/
}

#endif