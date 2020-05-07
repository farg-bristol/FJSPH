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

vector<real> GetWeightedDistance(SIM const& svar, FLUID const& fvar, State const& pnp1, outl const& outlist)
{
	uint const& end = svar.totPts;
	vector<real> wdiff(end,0.0);
	#pragma omp parallel for
	for (uint ii=0; ii< end; ++ii)
	{
		
	}

	return wdiff;
}

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

/* Colour field gradient in Hu et al (2014) method, note the bottom variable to account for no air*/
std::vector<StateVecD> GetColourGrad(SIM const& svar, FLUID const& fvar, State const& pnp1, outl const& outlist)
{
	std::vector<StateVecD> cgrad(svar.totPts, StateVecD::Zero());
	const uint& start = svar.bndPts;
	const uint& end = svar.totPts;

	#pragma omp parallel for shared(outlist) reduction(+:cgrad)
	for(uint ii=start; ii < end; ++ii)
	{
		StateVecD const& pi = pnp1[ii].xi;
		// StateVecD top = StateVecD::Zero();
		real bottom = 0.0;
		// pnp1[ii].normal = StateVecD::Zero();
		
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
			// pnp1[ii].normal += Rij;
		}
		
		cgrad[ii] = cgrad[ii]/bottom;

	}

	return cgrad;
}

StateVecD Base(FLUID const& fvar, Part const& pi, Part const& pj, 
	 StateVecD const& Rij, StateVecD const& Vij, real const r, 
	 StateVecD const& Grad, std::vector<real> &mu)
{
	/*Pressure and artificial viscosity - Monaghan (1994) p.400*/
	real const vdotr = Vij.dot(Rij);
	real const muij= fvar.H*vdotr/(r*r+0.001*fvar.HSQ);
	real pifac;
	if (vdotr > 0.0/* || pj.b == GHOST*/) 
	{
		pifac = 0.0;
	}
	else
	{
		real const rhoij = 0.5*(pi.rho+pj.rho);
		real const cbar= 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));
		pifac = fvar.alpha*cbar*muij/rhoij;
	}
	
	mu.emplace_back(muij);
	return Grad*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));
}

/*Laminar Viscosity - Morris (2003)*/
StateVecD Viscosity(FLUID const& fvar, Part const& pi, Part const& pj, 
	StateVecD const& Rij, StateVecD const& Vij, real const& r, StateVecD const& Grad)
{
	return -Vij*((fvar.mu)/(pi.rho*pj.rho))*(1.0/(r*r+0.01*fvar.HSQ))*Rij.dot(Grad);
}

/*Surface Tension - Nair & Poeschel (2017)*/
StateVecD SurfaceTens(FLUID const& fvar, Part const& pj, StateVecD const& Rij, 
					  real const& r, real const& npd)
{
	/*Surface tension factor*/
	real const lam = (6.0/81.0*pow((2.0*fvar.H),3.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));

	real fac=1.0; /*Boundary Correction Factor*/
	if(pj.b==BOUND) fac=(1+0.5*cos(M_PI*(fvar.contangb/180))); 

	/*npd = numerical particle density (see code above) */
	real sij = 0.5*pow(npd,-2.0)*(fvar.sig/lam)*fac;
	return -(Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*fvar.H));
}

StateVecD HuST(FLUID const& fvar, Part const& pi, Part const& pj, 
			  StateVecD const& Rij, real const& r, StateVecD const& cgradi, StateVecD const& cgradj)
{

	StateVecD Fst;

	Fst = (0.01/2.0)*(pi.m/pi.rho)*(pj.m/pj.rho)*((cgradi.squaredNorm()+cgradj.squaredNorm())/2.0)
		*W2GradK(Rij,r,fvar.H,fvar.correc);

	return Fst;
}

StateVecD NormalBoundaryRepulsion(FLUID const& fvar, MESH const& cells, Part const& pi)
{
// 	const vector<size_t> face = cells.faces[pi.faceID];
//     StateVecD norm;
// #if SIMDIM == 3
//     /*Get the face normal*/
//     StateVecD r1 = cells.verts[face[1]]- cells.verts[face[0]];
//     StateVecD r2 = cells.verts[face[2]]- cells.verts[face[0]];

//     norm = r1.cross(r2);
//     norm = norm.normalized();
// #else
//     StateVecD r1 = cells.verts[face[1]]-cells.verts[face[0]];
//     norm = StateVecD(-r1(1),r1(0)); 
//     norm = norm.normalized();
// #endif

   
//     real plane = norm.dot(cells.verts[face[1]]);
//     real dist =  (plane - pi.xi.dot(norm))/(norm.dot(norm));
    // StateVecD bpos = pi.xi + (dist-fvar.H)*norm;
    real beta = 4*fvar.Cs*fvar.Cs;

    real kern = BoundaryKernel(pi.y,fvar.H,beta);

	return fvar.bndM/(fvar.bndM+fvar.simM)*kern*pi.bNorm;
}

///**************** RESID calculation **************
void Forces(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, State const& pnp1/*, State& airP*/,
	 vector<vector<Part>> const& neighb, outl const& outlist,
	 vector<StateVecD>& RV, vector<real>& Rrho, std::vector<StateVecD>& Af/*, 
	 vector<real>& wDiff, vector<StateVecD>& norm, vector<real>& curve*/)
{
	svar.maxmu=0; 					    /* CFL Parameter */
	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;

	/*Gravity Vector*/
	#if SIMDIM == 3
		StateVecD const g = StateVecD(0.0,0.0,/*9.81*//*-9.81*/0.0);
	#else
		StateVecD const g = StateVecD(0.0,-9.81);
	#endif
	// const uint piston = svar.psnPts;

	/********* LOOP 1 - all points: Calculate numpartdens ************/
	// wDiff = GetWeightedDistance(svar,fvar,pnp1,outlist);
	// real numpartdens = GetNumpartdens(svar, fvar, pnp1, outlist);
	// std::vector<StateVecD> cgrad(svar.totPts,StateVecD::Zero());
	// cgrad = GetColourGrad(svar,fvar,pnp1,outlist);

	// std::vector<StateVecD> ST(svar.totPts,StateVecD::Zero()); /*Surface tension force*/		
	#pragma omp parallel /*shared(Af, wDiff, norm, curve)*/
	{

		vector<real> wDiff(end,0.0);
		vector<StateVecD> norm(end, StateVecD::Zero());
		vector<real> curve(end,0.0);

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
		// 		real r = Rij.norm();
		// 		StateVecD Grad = W2GradK(Rij, r,fvar.H,fvar.correc);
		// 		Rrhocontr[ii] -= pj.m*(Vij.dot(Grad));
		// 	}
		// 	pnp1[ii].Rrho = Rrhocontr[ii]; /*drho/dt*/
		// }



/******** LOOP 3 - Piston points: Calculate density and pressure. **********/		
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
			Part pi = pnp1[ii];
			size_t size = outlist[ii].size();
			// pnp1[ii].theta = real(size);
			// pnp1[ii].theta = wDiff[ii];

			vector<real> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.reserve(size+1);
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/		
			real mj = 0.0;
			StateVecD mRij = StateVecD::Zero();

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
				StateVecD const Grad = W2GradK(Rij,r,fvar.H,fvar.correc);

				/*Momentum contribution - Monaghan (1994)*/
				StateVecD const contrib = Base(fvar,pi,pj,Rij,Vij,r,Grad,mu);
				StateVecD visc = StateVecD::Zero();

				if (pj.b != GHOST)
				{
				/*Laminar Viscosity - Morris (2003)*/
				    visc = Viscosity(fvar,pi,pj,Rij,Vij,r,Grad);

				
					/*Surface Tension - Nair & Poeschel (2017)*/
					// StateVecD SurfC   = SurfaceTens(fvar,pj,Rij,r,numpartdens);
					// SurfC = HuST(fvar,pi,pj,Rij,r,cgrad[ii],cgrad[jj]);
					
					/* Surface Tension calcs */
					mRij+= (pi.m*pj.xi-pj.m*pi.xi);
					mj+= pj.m;

					// Find the normal vectors
					norm[ii] += (pj.m/pj.rho) * W2GradK(Rij,r,fvar.H,fvar.correc);

					/*drho/dt*/
					Rrho[ii] -= pj.m*(Vij.dot(Grad));
				}

				RV[ii] += pj.m*contrib + pj.m*visc /*+ SurfC/pj.m*/;


				

									
				// ST[ii] += SurfC/pj.m;
			}/*End of neighbours*/
			
			if(pi.internal == 1)
			{	// Apply the normal boundary force
				RV[ii] += NormalBoundaryRepulsion(fvar, cells, pi);
			}

			wDiff[ii] = mRij.norm()/(fvar.H*mj);

			

			RV[ii] += g;
			//CFL f_cv Calc
			real it = *max_element(mu.begin(),mu.end());
			if(it > svar.maxmu)
				svar.maxmu = it;
		} /*End of sim parts*/		
	

	// Get the surface curvatures
		#pragma omp for schedule(static) nowait /*reduction(+:RV, curve)*/
		for (size_t ii = start; ii < end; ++ii) 
		{
			real woccl = 0.0;
			// StateVecD norm = 0.0;
			// real wDiff = 0.0;
			StateVecD Vdiff;
			Part const& pi = pnp1[ii];
			size_t size = outlist[ii].size();
			
			if(pnp1[ii].b == FREE)
			{	
				// real kernsum = 0.0;
				if (svar.Asource == 1 || svar.Asource == 2)
				{
					Vdiff = (pi.cellV+cells.cPertnp1[pi.cellID]) - pi.v;
				}
				else 
				{
					Vdiff = avar.vInf - pi.v;
				}

#if  SIMDIM == 3
				if(svar.Asource == 3)
				{	
					StateVecD Vel = svar.vortex.getVelocity(pi.xi);
					Vdiff = Vel- pi.v;
				}
#endif
			}

			// Check if interaction is between surface particles
#if  SIMDIM == 2
			if(real(outlist[ii].size()) < 5.26*wDiff[ii]+30)
#else
			if(real(outlist[ii].size()) < 297.25*wDiff[ii]+34)
#endif
			{	
				
				real correc = 0.0;
				for (size_t const& jj:outlist[ii])
				{	/* Neighbour list loop. */
					Part const& pj = pnp1[jj];
					StateVecD Rij = pj.xi - pi.xi;
					if(pj.b != GHOST)
					{

#if  SIMDIM == 2
						if(real(outlist[jj].size()) < 5.26*wDiff[jj]+30)
#else
						if(real(outlist[jj].size()) < 297.25*wDiff[jj]+34)
#endif
						{	
							

							curve[ii] -= (pj.m/pj.rho)*(norm[jj].normalized()-norm[ii].normalized()).dot(
											W2GradK(Rij, Rij.norm(), fvar.H, fvar.correc));

							correc += (pj.m/pj.rho)*W2Kernel(Rij.norm(),fvar.H,fvar.correc);
						}

						/*Occlusion for Gissler Aero Method*/
						if (pi.b == FREE && (avar.acase == 4 || avar.acase == 1))
						{
							real const frac = -Vdiff.normalized().dot(Rij.normalized());
							
							// cout << num/denom << endl;
							if (frac > woccl)
							{
								woccl = frac;
							}
						}
					}
				}

				curve[ii] /= correc;
				RV[ii] += fvar.sig/pi.rho * curve[ii] * norm[ii];

				/*Find the aero force*/
				if(pi.b == FREE)
				{
	#if  SIMDIM == 2
					if(real(outlist[ii].size()) < 5.26*wDiff[ii]+30)
	#else
					if(real(outlist[ii].size()) < 297.25*wDiff[ii]+34)
	#endif
					{		
						RV[ii] += CalcAeroForce(avar,pi,Vdiff,norm[ii],size,woccl);
					}
				}
			}
			

		}

	}	/*End of declare parallel */
	/*Aerodynamic force*/
	// if(svar.Bcase > 1)
	// 	ApplyAero(svar,fvar,avar,cells,pnp1,outlist,RV,Af);
	
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
// 		real rho = 0.0;
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