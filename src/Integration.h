/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "VLM.h"
#include "Var.h"

/* Time integration class to solve between each write frame.

Can currently use either Newmark-Beta (2nd order semi-implicit) or 4th order Runge-Kutta integration
methods
 */
class Integrator
{
  public:
    Integrator(int solver) : solver_method(solver){};

    real integrate(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
        VLM const& vortex, MESH const& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist,
        SPHState& pn, SPHState& pnp1, vector<IPTState>& iptdata
    );

    void first_step(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
        VLM const& vortex, MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1,
        SPHState& pn, vector<IPTState>& iptdata
    );

  private:
    void solve_prestep(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
        MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1, SPHState& pn,
        vector<StateVecD>& xih, size_t const& start, size_t& end, real& logbase, real& npd,
        uint& iteration, StateVecD& Force, StateVecD& dropVel, real& error1, real& error2
    );

    void solve_timestep(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
        MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1, SPHState& pn,
        vector<StateVecD>& xih, size_t const& start, size_t& end, real& logbase, real& npd,
        uint& iteration, StateVecD& Force, StateVecD& dropVel, real& error1, real& error2
    );

    real find_timestep(
        SIM const& svar, FLUID const& fvar, MESH const& cells, SPHState const& pnp1, size_t const& start,
        size_t const& end_ng
    );

    int solver_method = newmark_beta;
    real maxf;
    real safe_dt;
};

#endif
