/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef AERO_H
#define AERO_H

#include "Kernel.h"
#include "Var.h"

inline real GetCd(real const& Re)
{
    return (1.0 + 0.197 * pow(Re, 0.63) + 2.6e-04 * pow(Re, 1.38)) * (24.0 / (Re + 0.00001));

    // if( Re <= 1000.0)
    // 	return (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
    // else
    // 	return 0.424;
}

inline StateVecD AeroForce(StateVecD const& Vdiff, AERO const& avar, real const& L_air, real const mass)
{
    real const Re = 2.0 * avar.rho_g * Vdiff.norm() * L_air / avar.mu_g;
    real const Cdi = GetCd(Re);

#if SIMDIM == 3
    real const Ai = M_PI * L_air * L_air;
#else
    real const Ai = M_PI * L_air;
#endif
    // acc[1] = 1.0*acc[1];
    // std::cout << "Reynolds: " << Re  << " Cd: " << Cd << " acc: " << acc[0] << " " << acc[1] <<
    // std::endl;
    return (0.5 * avar.rho_g * Vdiff.norm() * Vdiff * Cdi * Ai);
}

/*Sphere-Plate interpolation method - Gissler et al (2017)*/
inline StateVecD gissler_force(
    AERO const& avar, StateVecD const& Vdiff, StateVecD const& norm, SPHPart const& pi, real const& lam,
    real const& nnieghb, real const& woccl
)
{
    // real const n_full = avar.n_full;
    real const Re = 2.0 * pi.rho * Vdiff.norm() * pi.L_air / avar.mu_g;

    real frac2 = 0.0;
    if (avar.use_lam)
        frac2 = std::min(avar.interp_fac * lam, 1.0);
    else
        frac2 = std::min(avar.interp_fac * nnieghb * avar.i_n_full, 1.0);

    real const frac1 = (1.0 - frac2);

    real const Cds = GetCd(Re);

    real Cdl, Adrop;

#if SIMDIM == 3
    if (avar.use_TAB_def)
    {
        real ymax = Vdiff.squaredNorm() * pi.ycoef;
        if (ymax > 1.0)
            ymax = 1.0;
        Cdl = Cds * (1 + 2.632 * ymax);
        Adrop = M_PI * pow((pi.L_air + avar.Cb * pi.L_air * ymax), 2);
    }
    else
    {
        Cdl = Cds;
        Adrop = pi.A_sphere;
    }
#else
    if (avar.use_TAB_def)
    {
        real ymax = Vdiff.squaredNorm() * pi.ycoef;
        if (ymax > 1.0)
            ymax = 1.0;
        Cdl = Cds * (1 + 2.632 * ymax);
        Adrop = pi.A_sphere + 2 * (avar.Cb * pi.L_air * ymax) /** pow(pi.L_air,1)*/;
    }
    else
    {
        Cdl = Cds;
        Adrop = pi.A_sphere;
    }
#endif

    real const Cdi = frac1 * Cdl + /* 0.5* */ /*1.37**/ frac2;
    real const Aunocc = (frac1 * Adrop + frac2 * pi.A_plate);
    real const Ai = (1.0 - woccl) * Aunocc;

    // cout << "Fractions: " << frac1 << "  " << frac2 << endl;
    // cout << "Areas: " << Ai << "  " << Adrop << "  " << avar.A_plate << endl;
    // cout << "Re: " << Re << "  Cds: " << Cdi << "  "  << Cdl << "  " << Cds  << endl;
    // cout << "F: " << F(0) << "  " << F(1) << endl << endl;

    // return  0.5*rho*Vdiff.norm()*Vdiff*Cdi*Ai / mass;
    return 0.5 * Vdiff.norm() * Vdiff / (avar.sos * avar.sos) * avar.gamma * pi.p * Cdi * Ai / pi.m;
}

/*
Find the aerodynamic force via the induced pressure method. This calcualtes...

Returns:
    Aerodynamic force on a particle
*/
inline StateVecD induced_pressure(
    AERO const& avar, StateVecD const& Vdiff, StateVecD const& norm, real const& Pbasei, real const& lam,
    real const& nnieghb, SPHPart const& pi
)

{
    real theta = abs(acos(-norm.normalized().dot(Vdiff.normalized())));

    real Cp_s = 0.0;
    real Cp_p = 0.0;
    real Cp_b = 0.0;
    real Cp_tot = 0.0;

    /*Spherical Cp*/
    if (theta < 2.4455)
        Cp_s = 1.0 - (2.25) * pow(sin(theta), 2.0);
    else
        Cp_s = 0.075;

    /*Plate Cp*/
    if (theta < 1.570797)
        Cp_p = cos(theta);
    else if (theta < 1.9918)
        Cp_p = -pow(cos(6.0 * theta + 0.5 * M_PI), 1.5);
    else if (theta < 2.0838)
        Cp_p = 5.5836 * theta - 11.5601;
    else
        Cp_p = 0.075;

    /* Bowl Cp */
    if (theta < 0.7854)
        Cp_b = 1.0;
    else if (theta < 1.570797)
        Cp_b = 0.5 * (cos(4.0 * theta - M_PI) + 1.0);
    else
        Cp_b = 0.0;

    /* Overall Cp (Ideal to machine learn at some point) */
    real const normCurve = pi.norm_curve;
    real const fac1 = 0.25;
    real const ifac1 = 1 / fac1;

    if (normCurve < -fac1)
    { /* Bowl Cp */
        Cp_tot = Cp_b;
    }
    else if (normCurve < 0.0)
    { /* Interpolate from bowl Cp towards plate Cp */
        real const frac = (normCurve + fac1) * ifac1;
        Cp_tot = frac * Cp_b + (1.0 - frac) * Cp_p;
    }
    else if (normCurve < fac1)
    { /* Interpolate from plate Cp to sphere Cp */
        real const frac = (normCurve)*ifac1;
        Cp_tot = frac * Cp_p + (1.0 - frac) * Cp_s;
    }
    else
        Cp_tot = Cp_s;

    // real Plocali = 0.5*avar.rho_g*Vdiff.squaredNorm()*Cp_tot;

    /* Compressible dynamic pressure */
    real const Plocali =
        0.5 * Vdiff.squaredNorm() / (avar.sos * avar.sos) * avar.gamma * pi.cellP * Cp_tot;
    // cout << Plocalj << endl;

    real const Pi = (Plocali /* + Pbasei */);

    real const Re = pi.cellRho * Vdiff.norm() * pi.L_air / avar.mu_g;
    real const Cdi = GetCd(Re);

    /* Pure droplet force */
    StateVecD const acc_drop = 0.5 * Vdiff * Vdiff.norm() / (avar.sos * avar.sos) * avar.gamma *
                               pi.cellP * (M_PI * pi.L_air * pi.L_air * 0.25) * Cdi / pi.m;

    /* Induce pressure force */
    StateVecD const acc_kern = -/*0.5 **/ Pi * pi.A_plate * norm.normalized() / pi.m;

    /*Next, consider a skin friction force acting parallel to the surface*/
    real const Vnorm = Vdiff.dot(norm.normalized());

    /*Prandtl seventh power law for turbulent BL*/
    StateVecD const Vpar = Vdiff - Vnorm * norm.normalized();
    real const Re_par = pi.cellRho * Vpar.norm() * pi.L_air / avar.mu_g;
    real const Cf = 0.027 / pow(Re_par + 1e-6, 1.0 / 7.0);
    StateVecD const acc_skin = 0.5 * Vpar.norm() * Vpar / (avar.sos * avar.sos) * avar.gamma * pi.cellP *
                               Cf * pi.A_plate / pi.m;

    real frac1 = 0.0;
    if (avar.use_lam)
        frac1 = std::min(avar.interp_fac * lam, 1.0);
    else
        frac1 = std::min(avar.interp_fac * nnieghb * avar.i_n_full, 1.0);
    // real const frac2 = std::min(exp(pi.curve*0.001+200),1.0);

    return (frac1 * /*frac2**/ (acc_kern + acc_skin) + (1.0 - frac1) * acc_drop);
}

inline StateVecD CalcAeroAcc(
    AERO const& avar, SPHPart const& pi, StateVecD const& Vdiff, StateVecD const& norm, real const& lam,
    real const& nneigh, real const& Pbasei
)
{
    StateVecD acc = StateVecD::Zero();

    // cout << avar.acase << endl;
    switch (avar.acase)
    {
    case Gissler: /* Original Gissler */
    {
        acc = gissler_force(avar, Vdiff, norm, pi, lam, nneigh, pi.woccl);
        break;
    }
    case InducedPressure: /* Induced pressure based model */
    {
        acc = induced_pressure(avar, Vdiff, norm, Pbasei, lam, nneigh, pi);
        break;
    }
    case SkinFric: /*Skin Friction Method*/
    {
        /*Calculate the component of velocity in the surface normal direction (scalar projection)*/
        real Vnorm = Vdiff.dot(norm.normalized());

        if (Vnorm > 0.001)
        {
            real const Re = avar.rho_g * Vdiff.norm() * pi.L_air / avar.mu_g;

            /*Consider that pressure force is the stagnation of a velocity normal to the surface*/
            /*i.e. enforcing no parallel flow pressure and Cp = 1. */
            StateVecD acc_press =
                0.5 * avar.rho_g * Vnorm * Vnorm * pi.A_plate * norm.normalized() / pi.m;

            /*Next, consider a skin friction force acting parallel to the surface*/
            StateVecD Vpar = Vdiff - abs(Vnorm) * norm.normalized();
            real Cf = 0.027 / pow(Re, 1.0 / 7.0); /*Prandtl seventh power law for turbulent BL*/

            StateVecD acc_skin = 0.5 * avar.rho_g * Vpar.norm() * Cf * pi.A_plate * Vpar / pi.m;

            real const frac2 = std::min(1.5 * lam, 1.0);
            real const frac1 = (1.0 - frac2);

            /*Droplet force*/
            real const Cdi = GetCd(Re);
            StateVecD const acc_drop =
                0.5 * avar.rho_g * Vdiff.norm() * Vdiff * (M_PI * pi.L_air * pi.L_air / 4) * Cdi / pi.m;

            // cout << Fpress(0) << "  " << Fpress(1) << "  " << Fskin(0) << "  " << Fskin(1) << endl;

            acc = frac2 * (acc_press + acc_skin) + frac1 * acc_drop;
        }
        break;
    }
    default:
        break;
    }

    return acc;
}

#endif