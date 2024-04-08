/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef SHIFTING_H
#define SHIFTING_H

#include "Var.h"

/*L matrix for delta-SPH calculation*/
void dSPH_PreStep(
    FLUID const& fvar, size_t const& end, SPHState& pnp1, OUTL const& outlist,
    real& npd /* , DELTAP& dp */
);

/* Calculate dissipation terms before freezing. */
void dissipation_terms(
    FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist,
    /* DELTAP const& dp, */ SPHState& pnp1
);

#ifdef ALE
void Particle_Shift_No_Ghost(
    SIM const& svar, FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist,
    /* DELTAP const& dp, */ SPHState& pnp1
);

void Particle_Shift_Ghost(
    SIM const& svar, FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist,
    /* DELTAP const& dp, */ SPHState& pnp1
);
#endif

// void Apply_XSPH(FLUID const& fvar, size_t const& start, size_t const& end,
// 				OUTL const& outlist, DELTAP const& dp, SPHState& pnp1);

#endif