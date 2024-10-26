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
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    VLM const& vortex, MESH const& cells, LIMITS& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1
)
{
    start_index = svar.bound_points;
    end_index = svar.total_points;
    iteration = 0;

    real rms_error = 0.0;
    real logbase = 0.0;
    real npd = 1.0;

    // Find maximum safe timestep
    svar.delta_t = find_timestep(svar, fvar, cells, pnp1, start_index, end_index);

    outlist = update_neighbours(SPH_TREE, fvar, pnp1);

    dSPH_PreStep(fvar, svar.total_points, pnp1, outlist, npd);

    get_aero_velocity(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, vortex, start_index, end_index, outlist, limits,
        pn, pnp1, npd
    );

    Detect_Surface(svar, fvar, avar, start_index, end_index, outlist, cells, vortex, pnp1);

#ifndef NOFROZEN
    dissipation_terms(fvar, start_index, end_index, outlist, pnp1);
#endif

#ifdef ALE
    particle_shift(svar, fvar, start_index, end_index, outlist, pnp1);
#endif

    Check_Pipe_Outlet(CELL_TREE, svar, avar, cells, limits, pn, pnp1, end_index);

    vector<StateVecD> xih(end_index - start_index);
#pragma omp parallel for default(shared)
    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        xih[ii - start_index] = pnp1[ii].xi;
    }

    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    logbase = solve_prestep(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, limits, outlist, pn, pnp1, xih, start_index,
        end_index, npd
    );

    // cout << "Error: " << error1 << endl;

    /****** UPDATE NEIGHBOURS AND TREE ***********/
    outlist = update_neighbours(SPH_TREE, fvar, pnp1);

    dSPH_PreStep(fvar, end_index, pnp1, outlist, npd);

    get_aero_velocity(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, vortex, start_index, end_index, outlist, limits,
        pn, pnp1, npd
    );

    Detect_Surface(svar, fvar, avar, start_index, end_index, outlist, cells, vortex, pnp1);

#ifndef NOFROZEN
    dissipation_terms(fvar, start_index, end_index, outlist, pnp1);
#endif

    // Apply_XSPH(fvar,start_index,end_index,outlist,pnp1);

#ifdef ALE
    particle_shift(svar, fvar, start_index, end_index, outlist, pnp1);
#endif

    Check_Pipe_Outlet(CELL_TREE, svar, avar, cells, limits, pn, pnp1, end_index);

    /*Do time integration*/
    rms_error = solve_step(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, limits, outlist, pn, pnp1, xih, start_index,
        end_index, logbase, npd
    );

    return rms_error;
}

size_t Integrator::update_data(
    Sim_Tree& SPH_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells,
    LIMITS& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1, SURFS& surf_marks,
    vector<IPTState>& iptdata
)
{
    uint nAdd = update_buffer_region(svar, limits, pnp1, end_index);

    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        SPHPart& pi = pnp1[ii];
        if (pi.rho < 0.0001)
        {
            printf("Warning: Particle density is zero after updating buffer region.\n");
            pi.rho = fvar.rho_rest;
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

    if (svar.using_ipt && svar.Asource != constVel)
    {
        vector<IPTPart> IPT_nm1, IPT_n, IPT_np1;
        for (size_t const& ii : to_del)
        {
            IPT_nm1.emplace_back(IPTPart(pnp1[ii], svar.current_time, svar.ipt_diam, svar.ipt_area));
        }

        /* Do particle tracking on the particles converted */
        IPT_n = IPT_nm1;
        IPT_np1 = IPT_nm1;
#pragma omp parallel for default(shared)
        for (size_t ii = 0; ii < IPT_np1.size(); ++ii)
        {
            IPT::Integrate(
                svar, fvar, avar, cells, ii, IPT_nm1[ii], IPT_n[ii], IPT_np1[ii], surf_marks, iptdata
            );
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
        outlist = update_neighbours(SPH_TREE, fvar, pnp1);
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
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    VLM const& vortex, MESH const& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist, SPHState& pn,
    SPHState& pnp1, vector<IPTState>& ipt_data
)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    real step_error = integrate_no_update(
        SPH_TREE, CELL_TREE, svar, fvar, avar, vortex, cells, limits, outlist, pn, pnp1
    );

    size_t npts =
        update_data(SPH_TREE, svar, fvar, avar, cells, limits, outlist, pn, pnp1, surf_marks, ipt_data);

    // If there are no particles then exit early.
    if (npts == 0)
        return 0;

    // Process and output some useful information for each timestep
    real maxRho = 0.0;
    real maxAf = 0.0;
    real maxRhoi = 0.0;
#ifdef ALE
    real maxShift = 0.0;
#endif
#pragma omp parallel for reduction(max : maxAf, maxRhoi)
    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        SPHPart const& pi = pnp1[ii];
        maxAf = std::max(maxAf, pi.Af.norm());
        maxRhoi = std::max(maxRhoi, abs(pi.rho - fvar.rho_rest));
#ifdef ALE
        maxShift = std::max(maxShift, pi.vPert.norm());
#endif
    }

    maxRho = 100 * maxRhoi / fvar.rho_rest;
    real cfl = svar.delta_t / safe_dt;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2 - t1).count();
#ifdef ALE
    printf(
        "%9.3e | %7.2e | %4.2f | %8.4f | %7.3f | %9.3e | %9.3e | %9.3e | %14ld|\n", svar.current_time,
        svar.delta_t, cfl, step_error, maxRho, maxf, maxAf, maxShift, duration
    );

#else
    printf(
        "%9.3e | %7.2e | %4.2f | %9.4f | %3u | %8.3f | %9.3e | %9.3e | %14ld|\n", svar.current_time,
        svar.delta_t, cfl, step_error, iteration, maxRho, maxf, maxAf, duration
    );
#endif

    /*Add time to global*/
    svar.current_time += svar.delta_t;

    /* Check the error and adjust the CFL to try and keep convergence */
    if (step_error > svar.min_residual || maxRho > fvar.rho_max_iter)
    {
        if (step_error > 0.6 * svar.min_residual)
        { // If really unstable, immediately reduce the timestep
            svar.cfl = std::max(svar.cfl_min, svar.cfl - svar.cfl_step);
            svar.n_unstable = 0;
        }
        else if (svar.n_unstable > svar.n_unstable_limit)
        {
            svar.cfl = std::max(svar.cfl_min, svar.cfl - svar.cfl_step);
            svar.n_unstable = 0;
        }
        else
            svar.n_unstable++;
    }
    else
        svar.n_unstable = 0;

    if (iteration < svar.subits_factor * svar.max_subits && svar.n_unstable == 0)
    {
        if (svar.n_stable > svar.n_stable_limit)
        { // Try boosting the CFL if it does converge nicely
            svar.cfl = std::min(svar.cfl_max, svar.cfl + svar.cfl_step);
            svar.n_stable = 0;
        }
        else
            svar.n_stable++;
    }
    else
        svar.n_stable = 0;

#ifdef DAMBREAK
    /* Find max x and y coordinates of the fluid */
    // Find maximum safe timestep
    vector<SPHPart>::iterator maxXi = std::max_element(
        pnp1.begin() + svar.bound_points, pnp1.end_index(),
        [](SPHPart const& p1, SPHPart const& p2) { return p1.xi[0] < p2.xi[0]; }
    );

    vector<SPHPart>::iterator maxYi = std::max_element(
        pnp1.begin() + svar.bound_points, pnp1.end_index(),
        [](SPHPart const& p1, SPHPart const& p2) { return p1.xi[1] < p2.xi[1]; }
    );

    real maxX = maxXi->xi[0];
    real maxY = maxYi->xi[1];

    dambreak << svar.current_time << "  " << maxX + 0.5 * svar.particle_step << "  "
             << maxY + 0.5 * svar.particle_step << endl;
#endif

    return step_error;
}

/* Internal pre-step for the integration methods to find base error for the timestep. */
real Integrator::solve_prestep(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1,
    vector<StateVecD>& xih, size_t const& start_index, size_t& end_index, real& npd
)
{
    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    real logbase = 0.0;
    real rms_error = 0.0;
    switch (solver_method)
    {
    case runge_kutta:
    {
        logbase = Get_First_RK(
            svar, fvar, avar, start_index, end_index, npd, cells, limits, outlist, pn, pnp1
        );
        break;
    }
    case newmark_beta:
    {
        Newmark_Beta::Do_NB_Iter(
            CELL_TREE, svar, fvar, avar, start_index, end_index, npd, cells, limits, outlist, pn, pnp1
        );

        void(Newmark_Beta::Check_Error(
            SPH_TREE, svar, fvar, start_index, end_index, rms_error, logbase, outlist, xih, pn, pnp1,
            iteration
        ));
        iteration++; // Update iteration count
        break;
    }
    }
    return logbase;
}

/* Internal funtion to solve the timestep given a base error from the pre-step. */
real Integrator::solve_step(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pn, SPHState& pnp1,
    vector<StateVecD>& xih, size_t const& start_index, size_t& end_index, real& logbase, real& npd
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
            svar, fvar, avar, start_index, end_index, npd, cells, limits, outlist, logbase, pn, st_1,
            pnp1
        );
        break;
    }
    case newmark_beta:
    {

        rms_error = Newmark_Beta::Newmark_Beta(
            SPH_TREE, CELL_TREE, svar, fvar, avar, start_index, end_index, npd, cells, limits, outlist,
            logbase, iteration, xih, pn, pnp1
        );
        break;
    }
    }
    return rms_error;
}

real Integrator::find_timestep(
    SIM const& svar, FLUID const& fvar, MESH const& cells, SPHState const& pnp1,
    size_t const& start_index, size_t const& end_index
)
{
    // Find maximum safe timestep (avoiding a div by zero)
    real maxf = MEPSILON;
    real maxdrho = MEPSILON;
    real maxU = MEPSILON;
    real minST = 9999999.0;
#pragma omp parallel for reduction(max : maxf, maxU, maxdrho) reduction(min : minST)
    for (size_t ii = start_index; ii < end_index; ++ii)
    {
        maxf = std::max(maxf, pnp1[ii].acc.norm());
        maxdrho = std::max(maxdrho, abs(pnp1[ii].Rrho));
        maxU = std::max(maxU, pnp1[ii].v.norm());
        minST = std::min(
            minST,
            sqrt(pnp1[ii].rho * svar.dx * svar.dx / (2.0 * M_PI * fvar.sig * fabs(pnp1[ii].curve)))
        );
    }

    vector<real> timestep_factors;
    timestep_factors.emplace_back(0.25 * sqrt(fvar.H / maxf)); /* Force timestep constraint */
    timestep_factors.emplace_back(2 * fvar.H / (maxU));        /* Velocity constraint */
    timestep_factors.emplace_back(
        0.125 * fvar.H_sq * fvar.rho_rest / fvar.mu
    );                                                           /* Viscosity timestep constraint */
    timestep_factors.emplace_back(0.067 * minST);                /* Surface tension constraint */
    timestep_factors.emplace_back(0.5 * sqrt(fvar.H / maxdrho)); /* Density gradient constraint */
    timestep_factors.emplace_back(1.5 * fvar.H / fvar.speed_sound);
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
    real dt = svar.cfl * safe_dt;
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/

    if (dt < svar.delta_t_min)
        dt = svar.delta_t_min;
    else if (dt > svar.delta_t_max)
        dt = svar.delta_t_max;

    if (dt > svar.last_frame_time + svar.frame_time_interval - svar.current_time)
        dt = svar.last_frame_time + svar.frame_time_interval - svar.current_time + svar.delta_t_min;

    return dt;
}
