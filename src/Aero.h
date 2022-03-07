/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/***			        + delta-SPH                                                       ***/
/*** Viscosity:      Artificial Viscosity, linked to real laminar viscosity               ***/
/*** Surface Tension:   Ordoubadi (2014)                                                  ***/
/*** delta-SPH contrib: Marrone, Colagrossi, A. and Landrini, M. (2011):                  ***/
/*** Smoothing Kernel: Wendland's C2 ***/
/*** Integrator: Newmark-Beta ****/
/*** Variable Timestep Criteria: CFL + Marrone/Sun et al. (2017) conditions ***/

#ifndef AERO_H
#define AERO_H

#include "Var.h"
#include "Kernel.h"

inline real const GetCd(real const& Re)
{
	return (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.00001));

	// if( Re <= 1000.0)
	// 	return (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
	// else 
	// 	return 0.424;
}

inline StateVecD const AeroForce(StateVecD const& Vdiff, AERO const& avar, real const mass)
{
	real const Re = 2.0*avar.rhog*Vdiff.norm()*avar.L/avar.mug;
	real const Cdi = GetCd(Re);

	#if SIMDIM == 3
	real const Ai = M_PI*avar.L*avar.L;
	#else
	real const Ai = M_PI*avar.L;
	#endif
		//acc[1] = 1.0*acc[1];
	//std::cout << "Reynolds: " << Re  << " Cd: " << Cd << " acc: " << acc[0] << " " << acc[1] << std::endl;
	return (0.5*avar.rhog*Vdiff.norm()*Vdiff*Cdi*Ai);
}


/*Sphere-Plate interpolation method - Gissler et al (2017)*/
inline StateVecD const GisslerForce(AERO const& avar, StateVecD const& Vdiff, StateVecD const& norm, 
						real const& rho, real const& press, real const& mass, real const& lam, real const& woccl)
{
	// real const nfull = avar.nfull;
	real const Re = 2.0*rho*Vdiff.norm()*avar.L/avar.mug;
	
  real const frac2 = std::min(2.0 * lam, 1.0);
	real const frac1 = (1.0 - frac2);

 	real const Cds  = GetCd(Re);
	real Cdl, Adrop;

	#if SIMDIM == 3 
		if(avar.useDef)
		{
			real ymax = Vdiff.squaredNorm()*avar.ycoef;
			if (ymax > 1.0)
				ymax = 1.0;
			Cdl = Cds*(1+2.632*ymax);
			Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
		}
		else
		{
			Cdl = Cds;
			Adrop = avar.aSphere;
		}
	#else
		if(avar.useDef)
		{
			real ymax = Vdiff.squaredNorm()*avar.ycoef;
			if (ymax > 1.0)
				ymax = 1.0;
			 Cdl = Cds*(1+2.632*ymax);
			Adrop = avar.aSphere + 2*(avar.Cb*avar.L*ymax) /** pow(avar.L,1)*/;
		}
		else
		{
			Cdl = Cds;
			Adrop = avar.aSphere;
		}
	#endif

	real const Cdi = frac1*Cdl + /* 0.5* */ /*1.37**/frac2;
	real const Aunocc = (frac1*Adrop + frac2*avar.aPlate);
	real const Ai = (1.0 - woccl)*Aunocc;

	// cout << "Fractions: " << frac1 << "  " << frac2 << endl;
	// cout << "Areas: " << Ai << "  " << Adrop << "  " << avar.aPlate << endl;
	// cout << "Re: " << Re << "  Cds: " << Cdi << "  "  << Cdl << "  " << Cds  << endl;
	// cout << "F: " << F(0) << "  " << F(1) << endl << endl;


	// return  0.5*rho*Vdiff.norm()*Vdiff*Cdi*Ai / mass;
	return  0.5 * Vdiff.norm() * Vdiff / (avar.sos*avar.sos) *
			 avar.gamma * press * Cdi * Ai / mass;

}

inline StateVecD const InducedPressure(AERO const& avar, StateVecD const& Vdiff,
		 StateVecD const& norm, real const& Pbasei, real const& lam, SPHPart const& pi )
		 
{
	real theta = abs(acos(-norm.normalized().dot(Vdiff.normalized())));
		
	real Cp_p = 0.0;
	real Cp_s = 0.0;
	real Cp_tot = 0.0;

	// if(theta < 2.48073)
	// {
	// 	Cp = 1- 4*pow(sin(theta),3);
	// }
	// else
	// {
	// 	Cp = 0.075;
	// }

	// if(theta < 2.6399)
	// {
	// 	Cp = 1.0 - 2.0*pow(sin(theta),2.0);
	// }
	// else
	// {
	// 	Cp = 0.075;
	// }

	// if(theta < 2.4877)
	// {
	// 	Cp = 1.0 - 2.5*pow(sin(theta),2.0);
	// }
	// else
	// {
	// 	Cp = 0.075;
	// }

	/*Spherical Cp*/
	if(theta < 2.4455)
	{
		Cp_s = 1.0 - (2.25) * pow(sin(theta),2.0);
	}
	else
	{
		Cp_s = 0.075;
	}

	/*Plate Cp*/
	if(theta < 1.570797)
	{
		Cp_p = cos(theta);
	}
	else if (theta < 1.9918)
	{
		Cp_p = -pow(cos(6.0*theta+0.5*M_PI),1.5);
	}
	else if (theta < 2.0838)
	{
		Cp_p = 5.5836*theta-11.5601;
	}
	else
	{
		Cp_p = 0.075;
	}


	/* Overall Cp (Ideal to machine learn at some point) */
	if(pi.curve > 200.0)
	{
		if(theta < 1.570796)
			Cp_tot = 0.5*cos(2.0*theta)+0.5;
		else
			Cp_tot = 0.0;
	}
	else if (pi.curve > 0.0)
	{
		real const frac = (pi.curve) * 5.0e-3;
		Cp_tot = frac + (1.0-frac)*Cp_p;
	}
	else if (pi.curve > -200.0)
	{
		real const frac = (pi.curve + 200.0) * 5.0e-3;
		Cp_tot = frac*Cp_p + (1.0-frac)*Cp_s;
	}
	else
	{
		Cp_tot = Cp_s;
	}

	// real Plocali = 0.5*avar.rhog*Vdiff.squaredNorm()*Cp_tot;

	/* Compressible dynamic pressure */
	real const Plocali = 0.5*Vdiff.squaredNorm()/(avar.sos*avar.sos) * avar.gamma * pi.cellP * Cp_tot;
	// cout << Plocalj << endl;

	real const Pi = (Plocali/* + Pbasei */);
	
	// #pragma omp critical
	// cout << Cp_s << "  " << Cp_p << "  " << Cp_tot << "  " << theta << "  " << pi.curve << "  " << Pi <<  endl;

	real const Re = pi.cellRho*Vdiff.norm()*avar.L/avar.mug;
	real const Cdi  = GetCd(Re);

	/* Pure droplet force */
	StateVecD const acc_drop = 0.5*Vdiff*Vdiff.norm()/(avar.sos*avar.sos) * avar.gamma * pi.cellP
					* (M_PI*avar.L*avar.L*0.25) * Cdi/pi.m;

	// aeroD = -Pi * avar.aPlate * norm[ii].normalized();

	/* Induce pressure force */
	StateVecD const acc_kern =  -/*0.5 **/ Pi * avar.aPlate * norm.normalized()/pi.m;
	// StateVecD const F_kern = -(Pi/(pi.rho)) * norm;

	/*Next, consider a skin friction force acting parallel to the surface*/
	real const Vnorm = Vdiff.dot(norm.normalized());
	StateVecD const Vpar = Vdiff - Vnorm*norm.normalized();

	real const Re_par = pi.cellRho*Vpar.norm()*avar.L/avar.mug;
	real const Cf = 0.027/pow(Re_par+1e-6,1.0/7.0); /*Prandtl seventh power law for turbulent BL*/
	// real const Cf = 0;

	// StateVecD const acc_skin = StateVecD::Zero();
	// StateVecD const acc_skin = 0.5*avar.rhog* Vpar.norm() * Cf * avar.aPlate * Vpar/pi.m;
	StateVecD const acc_skin = 0.5 * Vpar.norm() * 	Vpar / (avar.sos*avar.sos) *
			 avar.gamma * pi.cellP * Cf * avar.aPlate / pi.m;

	real const frac1 = std::min(2.0 * lam, 1.0);
	// real const frac2 = std::min(exp(pi.curve*0.001+200),1.0);

	// #pragma omp critical
	// {
	// 	cout << acc_kern[0] << "  " << acc_kern[1] << "  " << acc_kern[2] << "  " 
	// 		<< acc_skin[0] << "  " << acc_skin[1] << "  " << acc_skin[2] << "  "
	// 		<< acc_drop[0] << "  " << acc_drop[1] << "  " << acc_drop[2] << endl;
	// }

	return (frac1 * /*frac2**/ (acc_kern+acc_skin) + (1.0-frac1) * acc_drop);
}

StateVecD CalcAeroAcc(AERO const& avar, SPHPart const& pi, StateVecD const& Vdiff,
		StateVecD const& norm, real const& lam, real const& Pbasei)
{
	StateVecD acc = StateVecD::Zero();

	// cout << avar.acase << endl;
	if( avar.acase == 1)
	{	/* Original Gissler */
		acc = GisslerForce(avar,Vdiff,norm,pi.cellRho,pi.cellP,pi.m,lam,pi.woccl);
	}
	else if(avar.acase == 2)
	{	/* Induced pressure based model */	
		acc = InducedPressure(avar,Vdiff,norm,Pbasei,lam,pi);
	}
	else if (avar.acase == 3)
	{	/*Skin Friction Method*/
		/*Calculate the component of velocity in the surface normal direction (scalar projection)*/
		real Vnorm = Vdiff.dot(norm.normalized());

		if(Vnorm > 0.001)
		{
			real const Re = avar.rhog*Vdiff.norm()*avar.L/avar.mug;

			/*Consider that pressure force is the stagnation of a velocity normal to the surface*/
			/*i.e. enforcing no parallel flow pressure and Cp = 1. */
			StateVecD acc_press = 0.5*avar.rhog*Vnorm*Vnorm*avar.aPlate * norm.normalized()/pi.m;

			/*Next, consider a skin friction force acting parallel to the surface*/
			StateVecD Vpar = Vdiff - abs(Vnorm)*norm.normalized();
			real Cf = 0.027/pow(Re,1.0/7.0); /*Prandtl seventh power law for turbulent BL*/

			StateVecD acc_skin = 0.5*avar.rhog* Vpar.norm() * Cf * avar.aPlate * Vpar/pi.m;

			real const frac2 = std::min(1.5 * lam, 1.0);
			real const frac1 = (1.0 - frac2);

			/*Droplet force*/
			real const Cdi = GetCd(Re);
			StateVecD const acc_drop = 0.5*avar.rhog*Vdiff.norm()*Vdiff*(M_PI*avar.L*avar.L/4)*Cdi/pi.m;

			// cout << Fpress(0) << "  " << Fpress(1) << "  " << Fskin(0) << "  " << Fskin(1) << endl;

			acc = frac2*(acc_press + acc_skin) + frac1*acc_drop;
		}

	}
	else if(avar.acase == 4)
	{
		// real ymax = Vdiff.squaredNorm()*avar.ycoef;

		// acc = AeroForce(Vdiff, avar, pi.m);
		// real correc = 1.0;
		// real Acorrect = 1.0-woccl;

		// /*Correction based on surface normal*/
		real theta = norm.normalized().dot(Vdiff.normalized());
		// real denom = norm.norm()*(Vdiff).norm();
		// real theta = num/denom;

		real Cp = 0.0;
		if(theta < 2.4435)
		{
			Cp = 1.1*cos(theta*2.03)-0.1;
		}
		else
		{
			Cp = 0.075;
		}

		real const frac2 = std::min(1.5 * lam, 1.0);
		real const frac1 = (1.0 - frac2);


		#if SIMDIM == 3 
			// real const Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
			real const Adrop = M_PI*pow(avar.L,2);
		#else
			// real const Adrop = 2*(avar.L + avar.Cb*avar.L*ymax);
			real const Adrop = 2*avar.L;
		#endif

		real const Aunocc = frac1*Adrop + frac2*avar.aPlate;

		real press = 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;

		#if SIMDIM == 3
			acc = -7.5*norm.normalized()*Aunocc*press;
		#else
			acc = -norm.normalized()*Aunocc*press;
		#endif

	}
	return acc;
}


#endif