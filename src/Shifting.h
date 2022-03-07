#ifndef SHIFTING_H
#define SHIFTING_H

#include "Var.h"
#include "Kernel.h"
#include "Resid.h"
#include "Third_Party/Eigen/LU"
#include "Third_Party/Eigen/Eigenvalues"


/*L matrix for delta-SPH calculation*/
void dSPH_PreStep(FLUID const& fvar, size_t const& end, SPHState const& pnp1, 
				 OUTL const& outlist, DELTAP& dp)
{
	/************   RENORMALISATION MATRIX CALCULATION    ************/
	dp.clear();
	
	vector<StateMatD> Lmat(end);
	vector<StateVecD> gradRho(end);
	vector<StateVecD> norm(end);
	vector<StateVecD> avgV(end);
	
	vector<real> lam(end);
	vector<real> lam_nb(end);
    vector<real> lam_ng(end);
	vector<real> kernsum(end);
	vector<real> colour(end);
	real npd = 0.0;
	#pragma omp parallel default(shared) reduction(+:npd) // shared(Lmat,gradRho,norm,lam,lam_nb,lam_ng,kernsum)
	{
		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < end; ++ii)
		{
			StateMatD Lmat_ = StateMatD::Zero();
			StateMatD Lmat_nb_ = StateMatD::Zero();
            StateMatD Lmat_ng_ = StateMatD::Zero();
			StateVecD gradRho_ = StateVecD::Zero();
			// StateVecD gradRho_2 = StateVecD::Zero();
			StateVecD norm_ = StateVecD::Zero();
			StateVecD avgV_ = StateVecD::Zero();
			real kernsum_ = 0.0;
			real colour_ = 0.0;
		
			SPHPart const& pi = pnp1[ii];
	
			for(std::pair<size_t,real> const& jj:outlist[ii])
			{
				/*Check if the position is the same, and skip the particle if yes*/
				if(ii == jj.first)
				{
					kernsum_ +=  fvar.correc;
					avgV_ += pi.v * fvar.correc;

					continue;
				}
				
				SPHPart const& pj = pnp1[jj.first];

				StateVecD const Rji = pj.xi - pi.xi;
				real const r = sqrt(jj.second);
				real const volj = pj.m/pj.rho;
				real const kern = Kernel(r,fvar.H,fvar.correc);
				StateVecD const Grad = GradK(-Rji,r,fvar.H,fvar.correc);

				Lmat_ -= volj * Rji * Grad.transpose();

				gradRho_ += volj * (pj.rho - pi.rho)*Grad;

				if(pj.b > PISTON)
				{
					Lmat_nb_ -= volj * Rji * Grad.transpose();
	
                    if(pj.b != GHOST)
                    {   /* Ignore ghost particles */   
                        norm_ += volj*Grad;
                        avgV_ += pj.v * kern;
                        kernsum_ += kern;
						colour_ += volj*kern;
						npd += kern;
                    }
				}

                if(pj.b != GHOST)
                {
                    Lmat_ng_ -= volj * Rji * Grad.transpose();
                }
                
			}

			Eigen::FullPivLU<StateMatD> lu(Lmat_);
			/*If only 1 particle is in the neighbourhood, tensile instability can occur*/
			/*In this case, L becomes singular, and invertability needs to be checked.*/
			StateMatD Linv = StateMatD::Zero();
			Linv(0,0) = 1;
			if(lu.isInvertible())
			{	
				Linv = lu.inverse();
			}

			/*L matrix computation for particle shifting, including all particles*/
			Eigen::SelfAdjointEigenSolver<StateMatD> es;
			es.computeDirect(Lmat_);
			lam[ii] = ((es.eigenvalues()).minCoeff()); //1 inside fluid, 0 outside fluid  

			/*L matrix computation for particle shifting, ignoring the boundary*/
			Eigen::SelfAdjointEigenSolver<StateMatD> es1;
			es1.computeDirect(Lmat_nb_);
			lam_nb[ii] = ((es1.eigenvalues()).minCoeff()); //1 inside fluid, 0 outside fluid    

            /*L matrix computation for surface identification, ignoring ghost and boundary particles*/
			Eigen::SelfAdjointEigenSolver<StateMatD> es2;
			es2.computeDirect(Lmat_ng_);
			lam_ng[ii] = ((es2.eigenvalues()).minCoeff()); //1 inside fluid, 0 outside fluid    
			
			/* Now the matrix has been found, density gradient can be calculated */
			// for(std::pair<size_t,real> const& jj:outlist[ii])
			// {
			// 	if(ii == jj.first)
			// 	{
			// 		continue;
			// 	}

			// 	SPHPart const& pj = pnp1[jj.first];

			// 	StateVecD const Rji = pj.xi - pi.xi;
			// 	real const r = sqrt(jj.second);
			// 	real const volj = pj.m/pj.rho;
			// 	StateVecD const Grad = GradK(-Rji,r,fvar.H,fvar.correc);

			// 	gradRho_2 += volj * (pj.rho - pi.rho) * Linv*Grad;
			// }

			gradRho_ = Linv*gradRho_;
			// #pragma omp critical
			// {
			// 	cout << gradRho_[0] << "  " << gradRho_[1] << "  " << gradRho_2[0] << "  " << gradRho_2[1] << endl;
			// }


			if(pi.b == BOUND)
			{
				lam[ii] = 1.0;
				lam_nb[ii] = 1.0;
                lam_ng[ii] = 1.0;
			}

			
			Lmat[ii] = Linv;
			gradRho[ii] = gradRho_ ;

			// StateVecD ntemp = Linv * norm_;
			// norm[ii] = (norm_.norm() * ntemp.normalized());
			norm[ii] = Linv.normalized() * norm_;
			// avgV[ii] = lu.inverse() * avgV_;

			

			avgV[ii] = (avgV_/kernsum_);
			kernsum[ii] = (kernsum_);
			colour[ii] = colour_;	
		}

		// #pragma omp for schedule(static) 
		// for(size_t ii = 0; ii < end; ++ii)
		// {
		// 	// StateVecD gradRho_ = StateVecD::Zero();
		// 	StateVecD norm_ = StateVecD::Zero();
		// 	StateVecD avgV_ = StateVecD::Zero();
		// 	SPHPart const& pi(pnp1[ii]);

		// 	for(size_t const& jj:outlist[ii])
		// 	{
		// 		SPHPart const& pj(pnp1[jj]);
		// 		/*Check if the position is the same, and skip the particle if yes*/
		// 		if(pi.partID == pj.partID)
		// 			continue;

		// 		StateVecD const Rji = pj.xi - pi.xi;
		// 		real const r = Rji.norm();
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

	dp.npd = npd/real(end);
	dp.L = Lmat;
	dp.gradRho = gradRho;
	dp.norm = norm;
	dp.avgV = avgV;
	dp.lam = lam;
	dp.lam_nb = lam_nb;
    dp.lam_ng = lam_ng;
	dp.kernsum = kernsum;
	dp.colour = colour;
}

/* Calculate dissipation terms before freezing. */
void dissipation_terms(FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist,
DELTAP const& dp, SPHState& pnp1)
{

	for(size_t ii = start; ii < end; ++ii )
	{
		SPHPart const& pi = pnp1[ii];
		
		StateVecD artViscI = StateVecD::Zero();

		#ifndef NODSPH
		real Rrhod = 0.0;
		#endif

		for (std::pair<size_t,real> const& jj : outlist[ii])
		{	/* Neighbour list loop. */
			
			/*Check if the position is the same, and skip the particle if yes*/
			if(ii == jj.first /* || pnp1[jj.first].b == BOUND */)
				continue;

			SPHPart const& pj = pnp1[jj.first];

			StateVecD const Rji = pj.xi-pi.xi;
			StateVecD const Vji = pj.v-pi.v;
			real const rr = jj.second;
			real const r = sqrt(rr);
			real const idist2 = 1.0/(rr + 0.001*fvar.HSQ);
			// #ifndef NODSPH	
			real const volj = pj.m/pj.rho;
			// #endif
			// StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);
			StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);/* gradK;*/


			// artViscI += pj.m*ArtVisc(fvar.nu,pi,pj,fvar,Rji,Vji,idist2,gradK);
			artViscI += aVisc(Rji,Vji,idist2,gradK,volj);

			#ifndef NODSPH
			if(pj.b != BOUND)
			Rrhod += Continuity_dSPH(Rji,idist2,fvar.HSQ,gradK,volj,dp.gradRho[ii],dp.gradRho[jj.first],pi,pj);
			#endif
		}

		pnp1[ii].aVisc = fvar.alpha * fvar.H * fvar.Cs * fvar.rho0 * artViscI / pi.rho;
		#ifndef NODSPH
		pnp1[ii].deltaD = fvar.dCont*Rrhod;
		#endif

	}
}

#ifdef ALE
void Particle_Shift_No_Ghost(SIM const& svar, FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist,
DELTAP const& dp, SPHState& pnp1)
{
	/*Implement particle shifting technique - Sun, Colagrossi, Marrone, Zhang (2018)*/

	// vector<SPHPart>::iterator maxUi = std::max_element(pnp1.begin(),pnp1.end(),
	// 	[](SPHPart p1, SPHPart p2){return p1.v.norm()< p2.v.norm();});

	// real const maxU = fvar.maxU/**maxUi->v.norm()*/;

	// real const maxU = fvar.maxU;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait
		for(size_t ii = start; ii < end; ++ii)
		{	
			/*Check if the particle is too close to the surface, and ignore them if so*/
			if(dp.lam_nb[ii] < 0.55)
				continue;

			SPHPart const& pi = pnp1[ii];
			
			StateVecD deltaU = StateVecD::Zero();
			StateVecD gradLam = StateVecD::Zero();
			real maxUij = 0.0;

			uint f = 0;
			real woccl = 0.0;

			for (std::pair<size_t,real> const& jj:outlist[ii])
			{	/* Neighbour list loop. */
				SPHPart const& pj = pnp1[jj.first];

				if(ii == jj.first /*|| pj.b == PartState.BOUND_*/)
					continue;

				StateVecD const Rji = pj.xi-pi.xi;
				real const r = sqrt(jj.second);
				real const volj = pj.m/pj.rho;
				real const kern = Kernel(r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);

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

			if(pi.b != BUFFER)
			{
				// gradLam = gradLam.normalized();
				// cout << gradLam(0) << "  " << gradLam(1) << endl;
				// StateVecD norm = dp.norm[ii].normalized();
				StateVecD norm = gradLam.normalized();
				
				/*Apply the partial shifting*/
				if (f == 0)
				{   /*SPHPart in the body of fluid, so treat as normal*/
					pnp1[ii].vPert = deltaU;
				}
				else /*Already checked if eigenvalue is less than 0.55*/
				{	/*SPHPart in the surface region, so apply a partial shifting*/
				    // if(norm.dot(deltaU) >= 0.0)
					// {	/*Check if particle is heading towards or away from the surface*/
						if(woccl < M_PI/12.0)
						{
							StateMatD shift = StateMatD::Identity() - norm*norm.transpose();
							pnp1[ii].vPert = /*dp.lam[ii]*/ shift*deltaU;
						}
						else
						{
							pnp1[ii].vPert = StateVecD::Zero();
						}
					// }
					// else
					// {
					// 	pnp1[ii].vPert = deltaU;
					// }
				}
				/*Otherwise, particle is near a surface, and don't shift.*/
			}
			else
			{
				pnp1[ii].vPert = StateVecD::Zero();
			}	
		}
	}
}

void Particle_Shift_Ghost(SIM const& svar, FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist,
DELTAP const& dp, SPHState& pnp1)
{
	/*Implement particle shifting technique - Sun, Colagrossi, Marrone, Zhang (2018)*/

	// vector<SPHPart>::iterator maxUi = std::max_element(pnp1.begin(),pnp1.end(),
	// 	[](SPHPart p1, SPHPart p2){return p1.v.norm()< p2.v.norm();});

	// real const maxU = fvar.maxU/**maxUi->v.norm()*/;

	// real const maxU = fvar.maxU;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait
		for(size_t ii = start; ii < end; ++ii)
		{	
			/*Check if the particle is too close to the surface, and ignore them if so*/
			if(dp.lam_nb[ii] < 0.55)
				continue;

			SPHPart const& pi = pnp1[ii];
			
			StateVecD deltaU = StateVecD::Zero();
			StateVecD gradLam = StateVecD::Zero();
			real maxUij = 0.0;

			real woccl = 0.0;

			for (std::pair<size_t,real> const& jj:outlist[ii])
			{	/* Neighbour list loop. */
				SPHPart const& pj = pnp1[jj.first];

				if(pj.partID == pi.partID /*|| pj.b == PartState.BOUND_*/)
					continue;

				StateVecD const Rji = pj.xi-pi.xi;
				real const r = sqrt(jj.second);
				real const volj = pj.m/pj.rho;
				real const kern = Kernel(r,fvar.H,fvar.correc);
				StateVecD const gradK = GradK(Rji,r,fvar.H,fvar.correc);

				deltaU += ( 1.0 + 0.2 * pow(kern/fvar.Wdx,4.0) ) * gradK * volj;
				// deltaU += 
				gradLam -= (dp.lam[jj.first]-dp.lam[ii]) * dp.L[ii] * gradK * volj;

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

			if(pi.b != BUFFER)
			{		
				/*Apply the partial shifting*/
                pnp1[ii].vPert = deltaU;
			}
			else
			{
				pnp1[ii].vPert = StateVecD::Zero();
			}	
		}
	}
}
#endif

// void Apply_XSPH(FLUID const& fvar, size_t const& start, size_t const& end, 
// 				OUTL const& outlist, DELTAP const& dp, SPHState& pnp1)
// {
// 	#pragma omp parallel
// 	{
// 		#pragma omp for schedule(static) nowait
// 		for(size_t ii = start; ii < end; ++ii)
// 		{	
// 			SPHPart const& pi(pnp1[ii]);

// 			StateVecD vPert_ = StateVecD::Zero();
// 			for (size_t const& jj:outlist[ii])
// 			{	 Neighbour list loop. 
// 				SPHPart const& pj(pnp1[jj]);

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
// real GetNumpartdens(SIM const& svar, FLUID const& fvar, SPHState const& pnp1, OUTL const& outlist)
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

#endif