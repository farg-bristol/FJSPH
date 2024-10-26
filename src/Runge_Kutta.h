/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "Var.h"

/* <summary> Peform the first stage of the Runge-Kutta integration
        to get the first guess of time n+1 (regarded here as time n+1/4)
        to perform neighbour search and dissipation terms before freezing </summary */
real Get_First_RK(
    SIM& svar, AERO const& svar.air, size_t const& start, size_t const& end, real const& npd,
    MESH const& cells, LIMITS const& limits, OUTL const& outlist, SPHState& part_n, SPHState& st_1
);

real Runge_Kutta4(
    SIM& svar, AERO const& svar.air, size_t const& start, size_t& end, real const& npd,
    MESH const& cells, LIMITS const& limits, OUTL const& outlist, real const& logbase, SPHState& part_n,
    SPHState& st_1, SPHState& part_np1
);

#endif