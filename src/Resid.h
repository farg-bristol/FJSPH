/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/*** Continuity:        Delta-SPH. Marrone et al. (2011)                                  ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)                                  ***/
/*** Smoothing Kernel:  Wendland's C2                                                     ***/
/*** Integrator:        Newmark-Beta or 4th order Runge-Kutta                             ***/
/*** Variable Timestep Criteria: CFL + Monaghan, J.J. (1989) conditions                   ***/

#ifndef RESID_H
#define RESID_H

#include "VLM.h"
#include "Var.h"

/* Boundary pressure calculation - Adami, Hu, and Adams, 2012 -
 * https://doi.org/10.1016/j.jcp.2012.05.005*/
void Get_Boundary_Pressure(
    StateVecD const& grav, FLUID const& fvar, size_t const& start, size_t const& end,
    OUTL const& outlist, SPHState& pnp1
);

void Boundary_DBC(
    FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1,
    vector<StateVecD>& RV
);

void Boundary_Ghost(
    FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1,
    vector<real>& Rrho, vector<int>& near_inlet
);

void Set_No_Slip(
    FLUID const& fvar, size_t const& start, size_t const& end, OUTL const& outlist, SPHState& pnp1
);

///**************** RESID calculation **************
void Forces(
    SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, SPHState const& pnp1,
    OUTL const& outlist, /* DELTAP const& dp, */ real const& npd, vector<StateVecD>& RV,
    vector<real>& Rrho, std::vector<StateVecD>& Af, StateVecD& Force
);

void Get_Aero_Velocity(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    MESH const& cells, VLM const& vortex, size_t const& start, size_t& end_ng, OUTL& outlist,
    LIMITS& limits, SPHState& pn, SPHState& pnp1, real& npd
);
#endif