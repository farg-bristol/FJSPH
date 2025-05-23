#ifndef SHIFTING_H
#define SHIFTING_H

#include "Helper_Functions.h"
#include "Kernel.h"
#include "Var.h"

#include "Third_Party/Eigen/Eigenvalues"
#include "Third_Party/Eigen/LU"

/*L matrix for delta-SPH calculation*/
void dSPH_PreStep(FLUID const& fvar, size_t const& end, SPHState& pnp1, OUTL const& outlist, real& npd)
{
/************   RENORMALISATION MATRIX CALCULATION    ************/
#ifdef PAIRWISE
    real npd_ = 0.0;
#pragma omp parallel default(shared) reduction(+ : npd)
#else
#pragma omp parallel default(shared)
#endif
    {
#pragma omp for schedule(static) nowait
        for (size_t ii = 0; ii < end; ++ii)
        {
            StateMatD Lmat_ = StateMatD::Zero();
            StateMatD Lmat_nb_ = StateMatD::Zero();
            StateVecD gradRho_ = StateVecD::Zero();
            // StateVecD gradRho_2 = StateVecD::Zero();
            StateVecD norm_ = StateVecD::Zero();
            // StateVecD avgV_ = StateVecD::Zero();
            real kernsum_ = 0.0;
            real colour_ = 0.0;

            SPHPart const& pi = pnp1[ii];

            for (neighbour_index const& jj : outlist[ii])
            {
                /*Check if the position is the same, and skip the particle if yes*/
                if (ii == jj.first)
                {
                    kernsum_ += fvar.W_correc;
                    // avgV_ += pi.v * fvar.W_correc;

                    continue;
                }

                SPHPart const& pj = pnp1[jj.first];
                StateVecD const Rji = pj.xi - pi.xi;
                real const r = sqrt(jj.second);
                real const volj = pj.m / pj.rho;
                real const kern = Kernel(r, fvar.H, fvar.W_correc);
                StateVecD const Grad = GradK(-Rji, r, fvar.H, fvar.W_correc);

                Lmat_ -= volj * Rji * Grad.transpose();

                gradRho_ += volj * (pj.rho - pi.rho) * Grad;

                if (pj.b > PISTON)
                {
                    Lmat_nb_ -= volj * Rji * Grad.transpose();

                    norm_ += volj * Grad;

                    // avgV_ += pj.v * kern;
                    kernsum_ += kern;
                    colour_ += volj * kern;
#ifdef PAIRWISE
                    npd_ += kern;
#endif
                }
            }

            Eigen::ColPivHouseholderQR<StateMatD> lu(Lmat_);
            /*If only 1 particle is in the neighbourhood, tensile instability can occur*/
            /*In this case, L becomes singular, and invertability needs to be checked.*/
            StateMatD Linv = StateMatD::Identity();

            if (lu.isInvertible())
            {
                Linv = lu.inverse();
            }

            gradRho_ = Linv * gradRho_;

            if (pi.b == BOUND)
            {
                pnp1[ii].lam = 1.0;
                pnp1[ii].lam_nb = 1.0;
            }
            else
            {
                /*L matrix computation for particle shifting, including all particles*/
                Eigen::SelfAdjointEigenSolver<StateMatD> es;
                es.computeDirect(Lmat_);
                pnp1[ii].lam = ((es.eigenvalues()).minCoeff()); // 1 inside fluid, 0 outside fluid

                /*L matrix computation for particle shifting, ignoring the boundary*/
                es.computeDirect(Lmat_nb_);
                pnp1[ii].lam_nb = ((es.eigenvalues()).minCoeff()); // 1 inside fluid, 0 outside fluid
            }

            pnp1[ii].L = Linv;
            pnp1[ii].gradRho = gradRho_;

            pnp1[ii].norm = Linv * norm_;

            // pnp1[ii].avgV = (avgV_/kernsum_);
            if (pnp1[ii].lam > 0.7)
                pnp1[ii].colourG = 2 /* *norm_.norm() */;
            else
                pnp1[ii].colourG = 2.0 * std::max(1.0, 1.0 / (2.0 * colour_)) /* *norm_.norm() */;

            pnp1[ii].kernsum = (kernsum_);
            pnp1[ii].colour = colour_;
        }

    } /*End parallel section*/

/* Need to make parallel */
#ifdef PAIRWISE
    npd = npd_ / real(end);
#endif
}

/* Calculate dissipation terms before freezing. */
void dissipation_terms(
    FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1
)
{

    for (size_t ii = start; ii < end; ++ii)
    {
        SPHPart const& pi = pnp1[ii];

        StateVecD artViscI = StateVecD::Zero();

#ifndef NODSPH
        real Rrhod = 0.0;
#endif

        for (neighbour_index const& jj : outlist[ii])
        { /* Neighbour list loop. */

            /*Check if the position is the same, and skip the particle if yes*/
            if (ii == jj.first /* || pnp1[jj.first].b == BOUND */)
                continue;

            SPHPart const& pj = pnp1[jj.first];

            StateVecD const Rji = pj.xi - pi.xi;
            StateVecD const Vji = pj.v - pi.v;
            real const rr = jj.second;
            real const r = sqrt(rr);
            real const idist2 = 1.0 / (rr + 0.0001 * fvar.H_sq);
            // #ifndef NODSPH
            real const volj = pj.m / pj.rho;
            // #endif
            StateVecD const gradK = GradK(Rji, r, fvar.H, fvar.W_correc);

            if (pj.b > PISTON)
#ifdef LINEAR
                artViscI += Linear_ArtVisc(Rji, Vji, idist2, gradK, volj);
#else
                artViscI += pj.m * ArtVisc(fvar, pi, pj, Rji, Vji, idist2, gradK);
#endif

#ifndef NODSPH
            if (pj.b > PISTON)
                Rrhod +=
                    Continuity_dSPH(Rji, idist2, fvar.H_sq, gradK, volj, pi.gradRho, pj.gradRho, pi, pj);
            else
                Rrhod += Continuity_dSPH(Rji, idist2, fvar.H_sq, gradK, volj, pi, pj);
#endif
        }

#ifdef LINEAR
        artViscI *= fvar.visc_alpha * fvar.H * fvar.speed_sound * fvar.rho_rest / pi.rho;
#endif

        pnp1[ii].aVisc = artViscI;

#ifndef NODSPH
        pnp1[ii].deltaD = fvar.dsph_cont * Rrhod;
#endif
    }
}

#ifdef ALE
void particle_shift(
    SIM const& svar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1
)
{
    /*Implement particle shifting technique - Sun, Colagrossi, Marrone, Zhang (2018)*/

    // vector<SPHPart>::iterator maxUi = std::max_element(pnp1.begin(),pnp1.end(),
    // 	[](SPHPart p1, SPHPart p2){return p1.v.norm()< p2.v.norm();});

    // real const maxU = svar.fluid.maxU/**maxUi->v.norm()*/;

    // real const maxU = svar.fluid.maxU;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static) nowait
        for (size_t ii = start; ii < end; ++ii)
        {
            /*Check if the particle is too close to the surface, and ignore them if so*/
            if (pnp1[ii].lam_nb < 0.55 || pnp1[ii].b == BUFFER)
            {
                pnp1[ii].vPert = StateVecD::Zero();
                continue;
            }

            SPHPart const& pi = pnp1[ii];

            StateVecD deltaU = StateVecD::Zero();
            // StateVecD gradLam = StateVecD::Zero();
            real maxUij = 0.0;

            real woccl = 0.0;

            for (neighbour_index const& jj : outlist[ii])
            { /* Neighbour list loop. */
                if (ii == jj.first /*|| pj.b == PartState.BOUND_*/)
                    continue;

                SPHPart const& pj = pnp1[jj.first];

                StateVecD const Rji = pj.xi - pi.xi;
                real const r = sqrt(jj.second);
                real const volj = pj.m / pj.rho;
                real const kern = Kernel(r, svar.fluid.H, svar.fluid.W_correc);
                StateVecD const gradK = GradK(Rji, r, svar.fluid.H, svar.fluid.W_correc);

                deltaU += (1.0 + 0.2 * pow(kern / svar.fluid.W_dx, 4.0)) * gradK * volj;
                // deltaU +=
                // gradLam -= (dp.lam[jj.first]-dp.lam[ii]) * dp.L[ii] * gradK * volj;

                if (pj.b > PISTON)
                { /* Does this want to neglect the boundary? */
                    real theta = acos(pi.norm.normalized().dot(pj.norm.normalized()));
                    if (theta > woccl)
                        woccl = theta;
                }

                if ((pj.v - pi.v).norm() > maxUij)
                {
                    maxUij = (pj.v - pi.v).norm();
                }
            }

            // deltaR *= -1 * fvar.sr * maxU / fvar.speed_sound;
            deltaU *= -2.0 * svar.fluid.H * pi.v.norm();

            deltaU = std::min(deltaU.norm(), std::min(maxUij / 2.0, svar.integrator.max_shift_vel)) *
                     deltaU.normalized();

            // gradLam = gradLam.normalized();
            // cout << gradLam(0) << "  " << gradLam(1) << endl;
            // StateVecD norm = dp.norm[ii].normalized();
            StateVecD norm = -pnp1[ii].norm.normalized();
            // StateVecD norm = gradLam.normalized();

            /*Apply the partial shifting*/
            if (pi.surfzone == 0 && pnp1[ii].lam_nb > 0.55)
            { /*SPHPart in the body of fluid, so treat as normal*/
                pnp1[ii].vPert = deltaU;
            }
            else /*Already checked if eigenvalue is less than 0.55*/
            {    /*SPHPart in the surface region, so apply a partial shifting*/
                 // if(norm.dot(deltaU) >= 0.0)
                 // {	/*Check if particle is heading towards or away from the surface*/
                if (woccl < M_PI / 12.0)
                {
                    StateMatD shift = StateMatD::Identity() - norm * norm.transpose();
                    pnp1[ii].vPert = /*dp.lam[ii]*/ shift * deltaU;
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
        }
    }
}

#endif

// void Apply_XSPH( size_t const& start, size_t const& end,
// 				OUTL const& outlist, SPHState& pnp1)
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

// 				if(pj.part_id == pi.part_id || pj.b == PartState.BOUND_)
// 				{
// 					continue;
// 				}

// 				real const r = (pj.xi-pi.xi).norm();;
// 				real const kern = Kernel(r,fvar.H,fvar.W_correc);
// 				real rho_ij = 0.5*(pi.rho + pj.rho);

// 				vPert_ -= 0.5 * pj.m/rho_ij * (pi.v-pj.v) * kern/*/pi.kernsum*/;

// 			}

// 			pnp1[ii].vPert = vPert_;
// 		}
// 	}
// }

#endif