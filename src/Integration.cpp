/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Integration.h"

#include "Containment.h"
#include "Geometry.h"
#include "Helper_Functions.h"
#include "IPT.h"
#include "Neighbours.h"
#include "Newmark_Beta.h"
#include "Resid.h"
#include "Runge_Kutta.h"
#include "Shifting.h"
#include "shapes/inlet.h"
#include <chrono>
#include <set>

using namespace std::chrono;

/* Integration loop to perform the first timestep and find forces only.

Solve the forces and pressures for particles according to a calculated stable timestep (current no option
to use a fixed timestep) and do not update the new timestep. This should be used only as a prestep to
write files with starting forces.
*/
real Integrator::integrate_no_update(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, MESH const& cells, LIMITS& limits,
    OUTL& outlist, SPHState& pn, SPHState& pnp1
)
{
    start_index = svar.bound_points;
    end_index = svar.total_points;
    iteration = 0;

    real rms_error = 0.0;
    real logbase = 0.0;
    real npd = 1.0;

    // Find maximum safe timestep
    svar.integrator.delta_t = find_timestep(svar, cells, pnp1, start_index, end_index);

    outlist = update_neighbours(svar.fluid, SPH_TREE, pnp1);

    dSPH_PreStep(svar.fluid, svar.total_points, pnp1, outlist, npd);

    get_aero_velocity(
        SPH_TREE, CELL_TREE, svar, cells, start_index, end_index, outlist, limits, pn, pnp1, npd
    );

    Detect_Surface(svar, start_index, end_index, outlist, cells, pnp1);

#ifndef NOFROZEN
    dissipation_terms(svar.fluid, start_index, end_index, outlist, pnp1);
#endif

#ifdef ALE
    particle_shift(svar, start_index, end_index, outlist, pnp1);
#endif

    Check_Pipe_Outlet(CELL_TREE, svar, cells, limits, pn, pnp1, end_index);

    vector<StateVecD> xih(end_index - start_index);
#pragma omp parallel for default(shared)
    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        xih[ii - start_index] = pnp1[ii].xi;
    }

    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    logbase = solve_prestep(
        SPH_TREE, CELL_TREE, svar, cells, limits, outlist, pn, pnp1, xih, start_index, end_index, npd
    );

    // cout << "Error: " << error1 << endl;

    /****** UPDATE NEIGHBOURS AND TREE ***********/
    outlist = update_neighbours(svar.fluid, SPH_TREE, pnp1);

    dSPH_PreStep(svar.fluid, end_index, pnp1, outlist, npd);

    get_aero_velocity(
        SPH_TREE, CELL_TREE, svar, cells, start_index, end_index, outlist, limits, pn, pnp1, npd
    );

    Detect_Surface(svar, start_index, end_index, outlist, cells, pnp1);

#ifndef NOFROZEN
    dissipation_terms(svar.fluid, start_index, end_index, outlist, pnp1);
#endif

    // Apply_XSPH(svar.fluid,start_index,end_index,outlist,pnp1);

#ifdef ALE
    particle_shift(svar, start_index, end_index, outlist, pnp1);
#endif

    Check_Pipe_Outlet(CELL_TREE, svar, cells, limits, pn, pnp1, end_index);

    /*Do time integration*/
    rms_error = solve_step(
        SPH_TREE, CELL_TREE, svar, cells, limits, outlist, pn, pnp1, xih, start_index, end_index,
        logbase, npd
    );

    return rms_error;
}

size_t Integrator::update_data(
    Sim_Tree& SPH_TREE, SIM& svar, MESH const& cells, LIMITS& limits, OUTL& outlist, SPHState& pn,
    SPHState& pnp1, SURFS& surf_marks, vector<IPTState>& iptdata
)
{
    uint nAdd = update_buffer_region(svar, limits, pnp1, end_index);

    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        SPHPart& pi = pnp1[ii];
        if (pi.rho < 0.0001)
        {
            printf("Warning: Particle density is zero after updating buffer region.\n");
            pi.rho = svar.fluid.rho_rest;
        }
    }

    /* Check if any particles need to be deleted, and turned to particle tracking */
    std::set<size_t> to_del;
    vector<size_t> n_del_per_block(limits.size(), 0);
    for (size_t block = svar.n_bound_blocks; block < svar.n_bound_blocks + svar.n_fluid_blocks; block++)
    {
        // Need to check if the delete plane is active.
        if (limits[block].delconst == default_val)
            continue;

#pragma omp parallel for default(shared)
        for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
        {
            if (pnp1[ii].xi.dot(limits[block].delete_norm) > limits[block].delconst)
            {
/* SPHPart is downstream enough to convert to IPT */
#pragma omp critical
                {
                    to_del.insert(ii);
                }
#pragma omp atomic
                n_del_per_block[block]++;
            }
        }
    }

    if (svar.ipt.using_ipt && svar.Asource != constVel)
    {
        vector<IPTPart> IPT_nm1, IPT_n, IPT_np1;
        for (size_t const& ii : to_del)
        {
            IPT_nm1.emplace_back(
                IPTPart(pnp1[ii], svar.integrator.current_time, svar.ipt.ipt_diam, svar.ipt.ipt_area)
            );
        }

        /* Do particle tracking on the particles converted */
        IPT_n = IPT_nm1;
        IPT_np1 = IPT_nm1;
#pragma omp parallel for default(shared)
        for (size_t ii = 0; ii < IPT_np1.size(); ++ii)
        {
            IPT::Integrate(svar, cells, ii, IPT_nm1[ii], IPT_n[ii], IPT_np1[ii], surf_marks, iptdata);
        }
    }

    /* Delete the old particles */
    uint nDel = 0;
    if (!to_del.empty())
    {
        for (std::set<size_t>::reverse_iterator itr = to_del.rbegin(); itr != to_del.rend(); ++itr)
        {
            pnp1.erase(pnp1.begin() + *itr);
            svar.total_points--;
            svar.fluid_points--;
            end_index--;
            end_index--;
            nDel++;
            svar.delete_count++;
        }

        size_t delshift = 0;
        for (size_t block = svar.n_bound_blocks; block < svar.n_bound_blocks + svar.n_fluid_blocks;
             block++)
        {
            // Shift the blocks, which will happen cumulatively
            limits[block].index.first -= delshift;
            delshift += n_del_per_block[block];
            limits[block].index.second -= delshift;

            /* Need to shift back vector and buffer vector for correct inlet */
            for (size_t& back : limits[block].back)
                back -= delshift;

            for (vector<size_t>& buffer : limits[block].buffer)
                for (size_t& part : buffer)
                {
                    part -= delshift;
                }
        }
    }

    if (nAdd != 0 || nDel != 0)
    {
        outlist = update_neighbours(svar.fluid, SPH_TREE, pnp1);
    }

    /****** UPDATE TIME N ***********/
    if (svar.total_points != pnp1.size())
    {
        cout << "Size mismatch. Total points not equal to array size. Stopping" << endl;
        exit(-1);
    }

    // cout << "Updating. Time: " << svar.current_time << "  dt: " << svar.delta_t << endl;
    if (pn.size() != pnp1.size())
        pn.resize(pnp1.size());

    copy_omp(pnp1, pn);

    return pnp1.size();
}

/* Integration loop to perform a timestep and move forward in time.

Solve the forces and pressures for particles according to a calculated stable timestep (current no option
to use a fixed timestep) and update the data of time n.
 */
real Integrator::integrate(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, MESH const& cells, SURFS& surf_marks,
    LIMITS& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1, vector<IPTState>& ipt_data
)
{
    INTEG_SETT& integ_sett = svar.integrator;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    real step_error = integrate_no_update(SPH_TREE, CELL_TREE, svar, cells, limits, outlist, pn, pnp1);

    size_t npts = update_data(SPH_TREE, svar, cells, limits, outlist, pn, pnp1, surf_marks, ipt_data);

    // If there are no particles then exit early.
    if (npts == 0)
        return 0;

    // Process and output some useful information for each timestep
    real cfl = integ_sett.delta_t / safe_dt;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2 - t1).count();
#ifdef ALE
    printf(
        "%9.3e | %7.2e | %4.2f | %8.4f | %7.3f | %9.3e | %9.3e | %9.3e | %14ld|\n",
        integ_sett.current_time, integ_sett.delta_t, cfl, step_error, maxRho_pc, maxf, maxAf, maxShift,
        duration
    );
#else
    printf(
        "%9.3e | %7.2e | %4.2f | %9.4f | %3u | %8.3f | %9.3e | %9.3e | %14ld|\n",
        integ_sett.current_time, integ_sett.delta_t, cfl, step_error, iteration, maxRho_pc, maxf, maxAf,
        duration
    );
#endif

    /*Add time to global*/
    svar.integrator.current_time += svar.integrator.delta_t;

    /* Check the error and adjust the CFL to try and keep convergence */
    if (step_error > svar.integrator.min_residual || maxRho_pc > svar.fluid.rho_max_iter)
    {
        if (step_error > 0.6 * integ_sett.min_residual)
        { // If really unstable, immediately reduce the timestep
            integ_sett.cfl = std::max(integ_sett.cfl_min, integ_sett.cfl - integ_sett.cfl_step);
            integ_sett.n_unstable = 0;
        }
        else if (integ_sett.n_unstable > integ_sett.n_unstable_limit)
        {
            integ_sett.cfl = std::max(integ_sett.cfl_min, integ_sett.cfl - integ_sett.cfl_step);
            integ_sett.n_unstable = 0;
        }
        else
            integ_sett.n_unstable++;
    }
    else
        integ_sett.n_unstable = 0;

    if (iteration < integ_sett.subits_factor * integ_sett.max_subits && integ_sett.n_unstable == 0)
    {
        if (integ_sett.n_stable > integ_sett.n_stable_limit)
        { // Try boosting the CFL if it does converge nicely
            integ_sett.cfl = std::min(integ_sett.cfl_max, integ_sett.cfl + integ_sett.cfl_step);
            integ_sett.n_stable = 0;
        }
        else
            integ_sett.n_stable++;
    }
    else
        integ_sett.n_stable = 0;

    return step_error;
}

/* Internal pre-step for the integration methods to find base error for the timestep. */
real Integrator::solve_prestep(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, MESH const& cells, LIMITS const& limits,
    OUTL& outlist, SPHState& pn, SPHState& pnp1, vector<StateVecD>& xih, size_t const& start_index,
    size_t& end_index, real& npd
)
{
    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    real logbase = 0.0;
    real rms_error = 0.0;
    switch (solver_method)
    {
    case runge_kutta:
    {
        logbase = Get_First_RK(svar, start_index, end_index, npd, cells, limits, outlist, pn, pnp1);
        break;
    }
    case newmark_beta:
    {
        Newmark_Beta::Do_NB_Iter(
            CELL_TREE, svar, start_index, end_index, npd, cells, limits, outlist, pn, pnp1
        );

        void(Newmark_Beta::Check_Error(
            SPH_TREE, svar, start_index, end_index, rms_error, logbase, outlist, xih, pn, pnp1, iteration
        ));
        iteration++; // Update iteration count
        break;
    }
    }
    return logbase;
}

/* Internal funtion to solve the timestep given a base error from the pre-step. */
real Integrator::solve_step(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, MESH const& cells, LIMITS const& limits,
    OUTL& outlist, SPHState& pn, SPHState& pnp1, vector<StateVecD>& xih, size_t const& start_index,
    size_t& end_index, real& logbase, real& npd
)
{
    /*Do time integration*/
    real rms_error = 0.0;
    switch (solver_method)
    {
    case runge_kutta:
    {
        SPHState st_1 = pnp1;

        rms_error = Runge_Kutta4(
            svar, start_index, end_index, npd, cells, limits, outlist, logbase, pn, st_1, pnp1
        );
        break;
    }
    case newmark_beta:
    {
        rms_error = Newmark_Beta::Newmark_Beta(
            SPH_TREE, CELL_TREE, svar, start_index, end_index, npd, cells, limits, outlist, logbase,
            iteration, xih, pn, pnp1
        );
        break;
    }
    }
    return rms_error;
}

real Integrator::find_timestep(
    SIM const& svar, MESH const& cells, SPHState const& pnp1, size_t const& start_index,
    size_t const& end_index
)
{
    INTEG_SETT const& integ_sett = svar.integrator;

    // Find maximum safe timestep (avoiding a div by zero)
    maxf = MEPSILON;
    maxAf = MEPSILON;
    maxRhoi = MEPSILON;
    maxdrho = MEPSILON;
    maxU = MEPSILON;
    minST = 9999999.0;
#pragma omp parallel for reduction(max : maxf, maxU, maxdrho, maxAf, maxRhoi) reduction(min : minST)
    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        maxf = std::max(maxf, pnp1[ii].acc.norm());
        maxAf = std::max(maxAf, pnp1[ii].Af.norm());

        maxdrho = std::max(maxdrho, abs(pnp1[ii].Rrho));
        maxRhoi = std::max(maxRhoi, abs(pnp1[ii].rho - svar.fluid.rho_rest));

        minST = std::min(
            minST,
            sqrt(pnp1[ii].rho * svar.dx * svar.dx / (2.0 * M_PI * svar.fluid.sig * fabs(pnp1[ii].curve)))
        );

        maxU = std::max(maxU, pnp1[ii].v.norm());
#ifdef ALE
        maxShift = std::max(maxShift, pnp1[ii].vPert.norm());
#endif
    }
    maxRho_pc = 100 * maxRhoi / svar.fluid.rho_rest;

    vector<real> timestep_factors;
    timestep_factors.emplace_back(0.25 * sqrt(svar.fluid.H / maxf)); /* Force timestep constraint */
    timestep_factors.emplace_back(2 * svar.fluid.H / (maxU));        /* Velocity constraint */
    timestep_factors.emplace_back(
        0.125 * svar.fluid.H_sq * svar.fluid.rho_rest / svar.fluid.mu
    );                                            /* Viscosity timestep constraint */
    timestep_factors.emplace_back(0.067 * minST); /* Surface tension constraint */
    timestep_factors.emplace_back(0.5 * sqrt(svar.fluid.H / maxdrho)); /* Density gradient constraint */
    timestep_factors.emplace_back(1.5 * svar.fluid.H / svar.fluid.speed_sound);
    /* Acoustic constraint */ /* 2* can't be used without delta-SPH it seems. Divergent in tensile
                                 instability */

    // Only use if -fno-finite-math-only is on
    // if (std::isinf(maxf))
    // {
    // 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
    // 	exit(-1);
    // }

    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    safe_dt = 0.75 * *std::min_element(timestep_factors.begin(), timestep_factors.end());
    real dt = integ_sett.cfl * safe_dt;
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/

    if (dt < integ_sett.delta_t_min)
        dt = integ_sett.delta_t_min;
    else if (dt > integ_sett.delta_t_max)
        dt = integ_sett.delta_t_max;

    if (dt > integ_sett.last_frame_time + integ_sett.frame_time_interval - integ_sett.current_time)
        dt = integ_sett.last_frame_time + integ_sett.frame_time_interval - integ_sett.current_time +
             integ_sett.delta_t_min;

    return dt;
}
