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

// #include "Eigen/Core"
// #include "Eigen/StdVector"
#include "Eigen/LU"
#include "Var.h"
#include "Kernel.h"
#include "Cross.h"
#include "Ghost.h"

StateVecD AeroForce(StateVecD &Vdiff, SIM &svar, FLUID &fvar, CROSS &cvar)
{
	StateVecD Fd = StateVecD::Zero();
	ldouble Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
	ldouble Cd = 0.0;

	if (Re < 3500)
	 	Cd = (1.0+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/(Re+0.000001));
	else
		Cd = (1+0.197*pow(Re,0.63)+2.6e-4*pow(Re,1.38))*(24.0/(Re+0.0001));

	// cout << "Reynolds: " << Re << " Cd: " << Cd << endl;
	//cout << fvar.Acorrect << endl;
	Fd = cvar.Acorrect*(2*svar.Pstep)*Cd*1.225*Vdiff.normalized()*Vdiff.squaredNorm()*svar.Pstep;
	Fd[1] = 0.05*Fd[1];
	// cout << "Reynolds: " << Re  << " Cd: " << Cd << " Fd: " << Fd[0] << " " << Fd[1] << endl ;
	return Fd;
}

ldouble GetNumpartdens(SIM &svar, FLUID &fvar, State &pnp1,outl &outlist)
{
	ldouble npd = 0.0;
	for (size_t i=0; i< svar.totPts; ++i)
	{
		StateVecD pi = pnp1[i].xi;
		for (size_t j=0; j<outlist[i].size(); ++j)
		{ /* Surface Tension calcs */
			StateVecD pj = pnp1[outlist[i][j]].xi;
			ldouble r = (pj-pi).norm();
			npd += W2Kernel(r,fvar.H,fvar.correc);
		}
	}
	return npd/(1.0*svar.totPts);
}

StateVecD Base(FLUID &fvar, Particle &pi, Particle &pj, 
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
	return pj.m*Grad*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));
}

/*Laminar Viscosity - Morris (2003)*/
StateVecD Viscosity(FLUID &fvar, Particle &pi,Particle &pj, 
	StateVecD &Rij, StateVecD &Vij, ldouble &r, StateVecD &Grad)
{
	return Vij*(pj.m*fvar.mu)/(pi.rho*pj.rho)*(1.0/(r*r+0.01*fvar.HSQ))*Rij.dot(Grad);
}


/*Surface Tension - Nair & Poeschel (2017)*/
StateVecD SurfaceTens(FLUID &fvar, Particle &pj, StateVecD &Rij, ldouble &r, ldouble &npd)
{
	/*Surface tension factor*/
	const static ldouble lam = (6.0/81.0*pow((2.0*fvar.H),4.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));

	ldouble fac=1.0;
	if(pj.b==true) fac=(1+0.5*cos(M_PI*(fvar.contangb/180)));
	ldouble sij = 0.5*pow(npd,-2.0)*(fvar.sig/lam)*fac;
	return (Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*fvar.H));
}

StateVecD ApplyAero(SIM &svar, FLUID &fvar, CROSS &cvar, 
	Particle &pi, StateVecD &SurfC)
{
	StateVecD Fd= StateVecD::Zero();
	if (svar.Bcase == 3 && pi.xi[1] > 2*svar.Pstep)
	{
		switch(cvar.acase)
		{
			case 0: /*No aero force*/
				break;

			case 1:	{  /*All upstream particles at start*/
				if( pi.b == 2)
				{
					StateVecD Vdiff = cvar.vInf - pi.V;
					Fd = AeroForce(Vdiff, svar, fvar, cvar);
				}

				break;
			}
			case 2:	{	/* Surface particles */
				if (SurfC.norm()*pow(svar.Pstep*100,3.7622)/pi.m > 2.5)
				{  /*				^ Need to tune this parameter... */
					StateVecD Vdiff = cvar.vInf - pi.V;
					Fd = AeroForce(Vdiff, svar, fvar, cvar);
				}
				break;
			}
			case 3: {	/* All upstream particles */
				if(pi.b == 2 /*&& */)
				{
					StateVecD Vdiff = cvar.vInf - pi.V;
					Fd = AeroForce(Vdiff, svar, fvar, cvar);
				}
				break;
			}
			case 4:	{	/* Surface particles proportional to ST*/
				/*Work in progress...*/
				break;
			}
			case 6: { /*Surface particles, with correction based on surface normal*/
				if (SurfC.norm()*fvar.sig*pow(fvar.H,2)/(pi.m*2.6E-06) > 0.8)
				{
					StateVecD Vdiff = cvar.vInf-pi.V;
					Fd = AeroForce(Vdiff, svar, fvar, cvar);
					/*Correction based on surface normal*/
					StateVecD surfNorm(SurfC[0],-SurfC[1]);
					ldouble num = surfNorm.dot(cvar.vInf);
					ldouble denom = surfNorm.norm()*cvar.vInf.norm();
					ldouble theta = acos(num/denom)/M_PI;
					pi.theta = theta;
					ldouble correc = 0.0;
					if (theta <= 1.0 )
					{
						correc = cvar.a*W2Kernel(2*theta,cvar.h1,1)+cvar.b*W2Kernel(2*theta,cvar.h2,1);
					}

					Fd = correc*Fd;
				}
				break;
			}
		}
	}
	else if (svar.Bcase == 3 && pi.xi[1] < svar.Pstep)
	{
		StateVecD Vdiff = cvar.vJet - pi.V;
		ldouble Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
		ldouble Cd = 0.1*(1+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/(Re+0.0001));

		Fd = Vdiff.normalized()*Cd*(2*svar.Pstep)*fvar.rho0*Vdiff.squaredNorm();
		Fd[1] += 9.81*pi.m;
	}
	return Fd;
}

///**************** RESID calculation **************
void Forces(Sim_Tree &NP1_INDEX, SIM &svar, FLUID &fvar, CROSS &cvar, State &pn, State &pnp1,outl &outlist)
{
	svar.maxmu=0; 					/* CFL Parameter */
	const static ldouble eps = fvar.eps; 		/* XSPH Influence Parameter*/
	
	
/********* LOOP 1 - all points: Calculate numpartdens ************/
	ldouble numpartdens = GetNumpartdens(svar, fvar, pnp1, outlist);

	// #pragma omp parallel
	{
		// #pragma omp for

/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/
		for (size_t i=0; i< svar.bndPts; ++i)
		{
			ldouble Rrhocontr = 0.0;
			Particle pi = pnp1[i];

			for(size_t j=0; j<outlist[i].size(); ++j)
			{
				Particle pj = pnp1[outlist[i][j]];
				StateVecD Rij = pj.xi-pi.xi;
				StateVecD Vij = pj.v-pi.v;
				ldouble r = Rij.norm();
				StateVecD Grad = W2GradK(Rij, r,fvar.H,fvar.correc);
				Rrhocontr -= pj.m*(Vij.dot(Grad));
			}
			pnp1[i].Rrho = Rrhocontr; /*drho/dt*/
		}


/******* LOOP 3 - only for ghost particle case: Find surface points. *********/
		if(svar.Bcase == 3 && cvar.acase == 5)
			Ghost_Particles(NP1_INDEX, numpartdens, svar, fvar, cvar, pn, pnp1,outlist);

/*		cout << "Npts: " << svar.totPts << " Air Count: " << svar.aircount << endl;
		cout << "Bound Parts: " << svar.bndPts << " Sim Points: " << svar.simPts << endl;
		cout << "Pn size: " << pn.size() << " PnP1 Size: " << pnp1.size() << endl;
		cout << "Outlist Size: " << outlist.size() << endl;*/


/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		for (size_t i=svar.bndPts; i< svar.totPts + svar.aircount; ++i)
		{
			Particle pi = pnp1[i];
			pi.V = pi.v;
			if (svar.Bcase == 3 && cvar.acase == 3)
				pi.b = 2;

			ldouble Rrhocontr = 0.0;
			StateVecD contrib= StateVecD::Zero();
			StateVecD visc = StateVecD::Zero();
			StateVecD SurfC= StateVecD::Zero();

			vector<double> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/

			for (size_t j=0; j < outlist[i].size(); ++j)
			{	/* Neighbour list loop. */
				Particle pj = pnp1[outlist[i][j]];

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.xi == pj.xi)
					continue;

				StateVecD Rij = pj.xi-pi.xi;
				StateVecD Vij = pj.v-pi.v;
				ldouble r = Rij.norm();
				ldouble Kern = W2Kernel(r,fvar.H, fvar.correc);
				StateVecD Grad = W2GradK(Rij, r,fvar.H, fvar.correc);

				contrib += Base(fvar,pi,pj,Rij,Vij,r,Grad,mu);

				/*Laminar Viscosity - Morris (2003)*/
				visc    -= Viscosity(fvar,pi,pj,Rij,Vij,r,Grad);

				/*Surface Tension - Nair & Poeschel (2017)*/
				SurfC   -= SurfaceTens(fvar,pj,Rij,r,numpartdens);

				/* XSPH Influence*/
				ldouble rhoij = 0.5*(pi.rho+pj.rho);
				pi.V += eps*(pj.m/rhoij)*Kern*Vij;

				/*drho/dt*/
				Rrhocontr -= pj.m*(Vij.dot(Grad));

				if (svar.Bcase == 3 && (cvar.acase == 3 || cvar.acase ==5))
				{
					ldouble num = -Rij.dot(cvar.vInf);
					ldouble denom = Rij.norm()*cvar.vInf.norm();
					if (num/denom > 0.98)
						pi.b = 1;
				}

			}

			/*Crossflow force*/
			StateVecD Fd = ApplyAero(svar,fvar,cvar,pi,SurfC);


			pi.Rrho = Rrhocontr; /*drho/dt*/
			pi.f= contrib + SurfC*fvar.sig/pi.m + Fd/pi.m;

			pi.Sf = SurfC*fvar.sig/pi.m;
			pi.Af = Fd/pi.m;
			pi.f[1] += 9.81; /*Add gravity*/

			pnp1[i]=pi; //Update the actual structure

			//CFL f_cv Calc
			ldouble it = *max_element(mu.begin(),mu.end());
			if (it > svar.maxmu)
				svar.maxmu=it;
		} /*End of sim parts*/
	}
}

///*Density Reinitialisation using Least Moving Squares as in A. Colagrossi (2003)*
void DensityReinit(FLUID &fvar, State &pnp1, outl &outlist)
{
	DensVecD one = DensVecD::Zero();
  one[0] = 1.0;

	for(size_t i=0; i< pnp1.size(); ++i)
	{
		DensMatD A= DensMatD::Zero();
		//Find matrix A.
		Particle pi = pnp1[i];
		for (size_t j=0; j< outlist[i].size(); ++j)
		{
			Particle pj = pnp1[outlist[i][j]];
			StateVecD Rij = pi.xi-pj.xi;
			DensMatD Abar = DensMatD::Zero();
			// Abar << 1      , Rij(0)        , Rij(1)        ,
			// 	    Rij(0) , Rij(0)*Rij(0) , Rij(1)*Rij(0) ,
			// 	    Rij(1) , Rij(1)*Rij(0) , Rij(1)*Rij(1) ;

      Abar(0,0) = 1;
      for (int ii = 0; ii < Rij.cols(); ++ii)
      {
        Abar(ii+1,0) = Rij[ii];
        Abar(0,ii+1) = Rij[ii];
        for (int jj = 0; jj<=ii; ++jj)
        {
          Abar(ii+1,jj+1) = Rij[ii]*Rij[jj];
        }
      }

			A+= W2Kernel(Rij.norm(),fvar.H,fvar.correc)*Abar*pj.m/pj.rho;
		}

		DensVecD Beta;
		//Check if A is invertible
		Eigen::FullPivLU<DensMatD> lu(A);
		if (lu.isInvertible())
			Beta = lu.inverse()*one;
		else
			Beta = (1)*one;

		//Find corrected kernel
		ldouble rho = 0.0;
		for (size_t j=0; j< outlist[i].size(); ++j)
		{
			StateVecD Rij = pi.xi-pnp1[outlist[i][j]].xi;
			rho += pnp1[outlist[i][j]].m*W2Kernel(Rij.norm(),fvar.H,fvar.correc)*
			(Beta(0)+Beta(1)*Rij(0)+Beta(2)*Rij(1));
		}

		pnp1[i].rho = rho;
	}
}



#endif
