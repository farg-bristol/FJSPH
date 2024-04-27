/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef NEWMARK_BETA_H
#define NEWMARK_BETA_H

#include "Var.h"

namespace Newmark_Beta
{
    int Check_Error(
        Sim_Tree& SPH_TREE, SIM& svar, FLUID const& fvar, size_t const& start, size_t const& end,
        real& error1, real& error2, real& logbase, OUTL& outlist, vector<StateVecD> const& xih,
        SPHState& pn, SPHState& pnp1, uint& iteration
    );

    void Do_NB_Iter(
        Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start,
        size_t& end, real const& npd, MESH const& cells, LIMITS const& limits, OUTL& outlist,
        SPHState& pn, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel
    );

    void Newmark_Beta(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
        size_t const& start, size_t& end, real const& npd, MESH const& cells, LIMITS const& limits,
        OUTL& outlist, real& logbase, uint& iteration, real& error1, real& error2,
        vector<StateVecD>& xih, SPHState& pn, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel
    );
} // namespace Newmark_Beta

#endif