/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef NEWMARK_BETA_H
#define NEWMARK_BETA_H

#include "Var.h"

namespace Newmark_Beta
{
    int Check_Error(
        Sim_Tree& SPH_TREE, SIM& svar, size_t const& start, size_t const& end, real& rms_error,
        real& logbase, OUTL& outlist, vector<StateVecD> const& xih, SPHState const& pn, SPHState& pnp1,
        uint& iteration
    );

    void Do_NB_Iter(
        Vec_Tree const& CELL_TREE, SIM& svar, AERO const& svar.air, size_t const& start, size_t& end,
        real const& npd, MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState const& pn,
        SPHState& pnp1
    );

    real Newmark_Beta(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, AERO const& svar.air,
        size_t const& start, size_t& end, real const& npd, MESH const& cells, LIMITS const& limits,
        OUTL& outlist, real& logbase, uint& iteration, vector<StateVecD>& xih, SPHState const& pn,
        SPHState& pnp1
    );
} // namespace Newmark_Beta

#endif