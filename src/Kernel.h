/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef KERNEL_H
#define KERNEL_H

#include "Var.h"

#ifdef CUBIC
///******Cubic Spline Kernel*******///
inline real Kernel(real const dist, real const H, real const W_correc)
{
    real const q = dist / H;

    if (q < 1.0)
        return W_correc * (1.0 - 1.5 * q * q * (1 - 0.5 * q));
    else if (q < 2.0)
        return W_correc * 0.25 * pow(2.0 - q, 3.0);

    return 0;
}

inline StateVecD GradK(StateVecD const& Rij, real const dist, real const H, real const W_correc)
{
    real const q = dist / H;

    if (q < 1.0)
        return Rij / (H * H) * W_correc * (3.0 * (0.75 * q - 1.0));
    else if (q < 2.0)
        return Rij / (dist * H) * W_correc * -0.75 * (2.0 - q) * (2.0 - q);

    return StateVecD::Zero();
}

#else
///******Wendland's C2 Quintic Kernel*******///
inline real Kernel(real const& dist, real const& H, real const& W_correc)
{
    // if(dist/H > 2.0)
    // {
    // 	// cout << "Distance greater than 2H" << endl;
    // 	return 0.0;
    // }
    return (pow(1 - 0.5 * dist / H, 4)) * (2 * dist / H + 1) * W_correc;
}

/*Gradient*/
inline StateVecD GradK(StateVecD const& Rij, real const& dist, real const& H, real const& W_correc)
{
    if (dist / H < 1e-12)
    {
        cout << "Points are too close" << endl;
        return StateVecD::Zero();
    }
    // else if(dist/H > 2.0)
    // {
    // 	// cout << "Distance greater than 2H" << endl;
    // 	return StateVecD::Zero();
    // }
    return 5.0 * (Rij / (H * H)) * pow(1 - 0.5 * dist / H, 3) * W_correc;
}
#endif

inline real BoundaryKernel(real const dist, real const H, real const beta)
{
    real const q = dist / H;
    if (q < 2.0 / 3.0)
        return beta * 2.0 / 3.0;
    else if (2.0 / 3.0 <= q && q < 1.0)
        return beta * (2 * q - 3.0 / 2.0 * q * q);
    else if (1 <= q && q < 2)
        return 0.5 * beta * pow((2 - q), 2);

    return 0;
}

#ifdef HESF
/*Diffuse interface method from He et al (2014) for surface tension*/
inline StateVecD HeST(StateVecD const& cgradi, StateVecD const& cgradj, StateVecD const& diffK)
{
    return 0.5 * ((cgradi.squaredNorm() + cgradj.squaredNorm())) * diffK;
}

#endif

// /*Surface Tension - Nair & Poeschel (2017)*/
#ifdef PAIRWISE
inline real pairwise_ST_fac(int const bA, int const bB)
{
    if (bA == BOUND || bB == BOUND)
    {
        real contang = 0.5 * M_PI * 7.0 / 9.0;
        return (1.0 + 0.5 * cos(contang));
    }
    return 1.0;
}

// inline StateVecD SurfaceTens(StateVecD const& Rji, real const& r, real const& h, real const& sig, real
// const& lam, 					real const& npdm2, real const& pi3o4, int const& bA, int
// const& bB, real const& voli, real const& volj)
inline StateVecD pairwise_ST(
    StateVecD const& Rji, real const& r, real const& h, real const& npdm2, real const& pi3o4,
    int const& bA, int const& bB
)
{
    real const fac = pairwise_ST_fac(bA, bB);

    /*npd = numerical particle density (see code above) */
    // return -0.5*npdm2*(sig/lam)*fac*cos(pi3o4*r/h)*(Rji/r) * voli * volj;
    // (void) voli;
    // (void) volj;
    return -npdm2 * fac * cos(pi3o4 * r / h) * (Rji / r);
}
#endif

#ifdef CSF
inline void CSF_Curvature(
    DELTAP const& dp, SPHPart const& pi, SPHPart const& pj, StateVecD const& Rji, StateVecD const& gradK,
    real const& volj, real const& r, real& curve, real& W_correc
)
{
    if (dp.lam[pj.part_id] < 0.75)
    {
        curve -= (pj.norm.normalized() - pi.norm.normalized()).dot(volj * gradK);
        W_correc += volj * Kernel(r, svar.fluid.H, svar.fluid.W_correc) /*/dp.kernsum[ii]*/;
    }

    /* Consider imaginary particles to improve curvature representation */
    if (pi.surf == 1 && pj.surf != 1)
    {
        /* Host particle is on the surface, but neighbour is not. Reflect stuff then */
        StateVecD normj = 2 * pi.norm.normalized() - pj.norm.normalized();
        curve += (normj.normalized() - pi.norm.normalized()).dot(volj * gradK);
        W_correc += volj * Kernel(r, svar.fluid.H, svar.fluid.W_correc) /*/dp.kernsum[ii]*/;
    }

    if (pj.surf == 1 && pi.surf != 1)
    {
        /* Neighbour particle is on the surface, but host is not. Not as easy */
        /* Need to check if extended particle is still in the neighbourhood */
        if (r < svar.fluid.H)
        { /* Assuming a support radius of 2h, if the current particle is within h */
            /* it will still be in the neighbourhood if distance is doubled */
            StateVecD gradK2 = GradK(2 * Rji, 2 * r, svar.fluid.H, svar.fluid.W_correc);
            StateVecD normj = 2 * pj.norm.normalized() - pi.norm.normalized();
            curve -= (normj.normalized() - pi.norm.normalized()).dot(volj * gradK2);
            W_correc += volj * Kernel(r, svar.fluid.H, svar.fluid.W_correc) /*/dp.kernsum[ii]*/;
        }
    }
}
#endif

inline StateVecD BasePos(SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
    /* Positive Pressure - Monaghan (1994)*/
    return pj.m * gradK * (pi.p / (pi.rho * pi.rho) + pj.p / (pj.rho * pj.rho));
}

inline StateVecD BaseTIC(SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
    /* Negative Pressure - Monaghan (1994)*/
    return pj.m * gradK * (fabs(pi.p / (pi.rho * pi.rho)) + pj.p / (pj.rho * pj.rho));
}

inline StateVecD pressPos(real const& volj, SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
    /* Positive Pressure - Marrone et al. (2011)*/
    return gradK * (pi.p + pj.p) * volj;
}

inline StateVecD pressTIC(real const& volj, SPHPart const& pi, SPHPart const& pj, StateVecD const& gradK)
{
    /* Negative Pressure - Marrone et al. (2011)*/
    return gradK * (fabs(pi.p) + pj.p) * volj;
}

#ifdef ALE
/* Arbitrary Lagrangian Eulerian formulation - Sun, Colagrossi, Marrone, Zhang (2018)*/
inline StateVecD
ALEMomentum(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK)
{
    return ((pj.v * pj.vPert.transpose() + pi.v * pi.vPert.transpose()) * gradK -
            pi.v * (pj.vPert - pi.vPert).dot(gradK)) *
           Vj;
}

inline real ALEContinuity(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK)
{
    return ((pj.v + pj.vPert) - (pi.v + pi.vPert)).dot(gradK) * Vj;
}

inline real ALECont2ndterm(SPHPart const& pi, SPHPart const& pj, real const& Vj, StateVecD const& gradK)
{
    return (pj.rho * pj.vPert + pi.rho * pi.vPert).dot(gradK) * Vj;
}
#endif

/* delta-SPH dissipation term in the continuity equation*/
#ifndef NODSPH
inline real Continuity_dSPH(
    StateVecD const& Rji, real const& idist2, real const& H_sq, StateVecD const& Grad, real const& volj,
    StateVecD const& gRho_i, StateVecD const& gRho_j, SPHPart const& pi, SPHPart const& pj
)
{
    return volj * ((pj.rho - pi.rho) + 0.5 * (gRho_i + gRho_j).dot(Rji)) * Rji.dot(Grad) * idist2;
}

inline real Continuity_dSPH(
    StateVecD const& Rji, real const& idist2, real const& H_sq, StateVecD const& Grad, real const& volj,
    SPHPart const& pi, SPHPart const& pj
)
{
    return volj * (pj.rho - pi.rho) * Rji.dot(Grad) * idist2;
}
#endif

inline StateVecD ArtVisc(
    real const& nu, SPHPart const& pi, SPHPart const& pj, StateVecD const& Rji, StateVecD const& Vji,
    real const idist2, StateVecD const& gradK
)
{
    // if(pj.b != BOUND)
    // {
    real const vdotr = Vji.dot(Rji);

    if (vdotr > 0.0)
    {
        return StateVecD::Zero();
    }
    else
    {
        real const muij = svar.fluid.H * vdotr * idist2;
        real const rhoij = 0.5 * (pi.rho + pj.rho);
        real const cbar = 0.5 * (sqrt((svar.fluid.B * svar.fluid.gam) / pi.rho) +
                                 sqrt((svar.fluid.B * svar.fluid.gam) / pj.rho));
        // real const alphv = 2 * nu * (SIMDIM + 2)/(svar.fluid.H*cbar);
        // real const visc_alpha = std::max(svar.fluid.visc_alpha,alphv);
        real const visc_alpha = svar.fluid.visc_alpha;

        return gradK * visc_alpha * cbar * muij / rhoij;
    }
    // }
    return StateVecD::Zero();
}

inline StateVecD Linear_ArtVisc(
    StateVecD const& Rji, StateVecD const& Vji, real const idist2, StateVecD const& gradK,
    real const& volj
)
{
    return volj * idist2 * Vji.dot(Rji) * gradK;
}

/*Laminar Viscosity - Morris (2003)*/
/*Apparently divergent at free surfaces - consider removing.*/
// StateVecD Viscosity(real const& mu, real const& H_sq, SPHPart const& pi, SPHPart const& pj,
// 	StateVecD const& Rji, StateVecD const& Vji, real const& r, StateVecD const& gradK)
// {
// 	return Vji*(mu/(pi.rho*pj.rho))*(1.0/(r*r+0.01*H_sq))*Rji.dot(gradK);
// }

inline StateVecD Viscosity(
    real const& nu, SPHPart const& pi, SPHPart const& pj, StateVecD const& Rji, StateVecD const& Vji,
    real const& idist2, StateVecD const& gradK
)

{
    return nu * (pi.rho + pj.rho) / (pi.rho * pj.rho) * (Rji.dot(gradK)) * idist2 * Vji;
}

/*Repulsion for interacting with mesh surface - saves generating particles on surface*/
inline StateVecD NormalBoundaryRepulsion(MESH const& cells, SPHPart const& pi)
{
    real beta = 4 * svar.fluid.speed_sound * svar.fluid.speed_sound;
    real kern = BoundaryKernel(pi.y, svar.fluid.H, beta);
    return svar.fluid.bnd_mass / (svar.fluid.bnd_mass + svar.fluid.sim_mass) * kern * pi.bNorm;
}

#endif