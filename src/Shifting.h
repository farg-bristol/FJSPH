/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef SHIFTING_H
#define SHIFTING_H

#include "Var.h"

/*L matrix for delta-SPH calculation*/
void dSPH_PreStep(FLUID const& fvar, size_t const& end, SPHState& pnp1, OUTL const& outlist, real& npd);

/* Calculate dissipation terms before freezing. */
void dissipation_terms(
    FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1
);

#ifdef ALE
void particle_shift(
    SIM const& svar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1
);
#endif

// void Apply_XSPH( size_t const& start, size_t const& end,
// 				OUTL const& outlist, SPHState& pnp1);

#endif