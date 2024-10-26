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
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, AERO const& svar.air,
        VLM const& vortex, MESH const& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist,
        SPHState& pn, SPHState& pnp1, vector<IPTState>& iptdata
    );

    real integrate_no_update(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, AERO const& svar.air,
        VLM const& vortex, MESH const& cells, LIMITS& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1
    );

  private:
    size_t update_data(
        Sim_Tree& SPH_TREE, SIM& svar, AERO const& svar.air, MESH const& cells, LIMITS& limits,
        OUTL& outlist, SPHState& pn, SPHState& pnp1, SURFS& surf_marks, vector<IPTState>& iptdata
    );

    real solve_prestep(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, AERO const& svar.air,
        MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1,
        vector<StateVecD>& xih, size_t const& start_index, size_t& end_index, real& npd
    );

    real solve_step(
        Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, AERO const& svar.air,
        MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1,
        vector<StateVecD>& xih, size_t const& start_index, size_t& end_index, real& logbase, real& npd
    );

    real find_timestep(
        SIM const& svar, MESH const& cells, SPHState const& pnp1, size_t const& start,
        size_t const& end_ng
    );

    int solver_method = newmark_beta;
    real maxf;
    real safe_dt;

    size_t start_index, end_index;
    uint iteration;
};

#endif
