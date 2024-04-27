/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Integration.h"

#include "Add.h"
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

/* Integration loop to perform a timestep and move forward in time.

Solve the forces and pressures for particles according to a calculated stable timestep (current no option
to use a fixed timestep) and update the data of time n.
 */
real Integrator::integrate(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    VLM const& vortex, MESH const& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist, SPHState& pn,
    SPHState& pnp1, vector<IPTState>& iptdata
)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    size_t const start = svar.bndPts;
    size_t end_ng = svar.bndPts + svar.simPts;

    real error1 = 0.0;
    real error2 = 1.0;
    real logbase = 0.0;
    real npd = 1.0;

    // Find maximum safe timestep
    svar.dt = find_timestep(svar, fvar, cells, pnp1, start, end_ng);

    outlist = update_neighbours(SPH_TREE, fvar, pnp1);

    dSPH_PreStep(fvar, svar.totPts, pnp1, outlist, npd);

    Get_Aero_Velocity(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, vortex, start, end_ng, outlist, limits, pn, pnp1,
        npd
    );

    Detect_Surface(svar, fvar, avar, start, end_ng, outlist, cells, vortex, pnp1);

    if (svar.ghost == 1)
    { /* Poisson points */
        PoissonGhost(svar, fvar, avar, cells, SPH_TREE, outlist, pn, pnp1);
        dSPH_PreStep(fvar, end_ng, pnp1, outlist, npd);
    }
    else if (svar.ghost == 2)
    { /* Lattice points */
        LatticeGhost(svar, fvar, avar, cells, SPH_TREE, outlist, pn, pnp1, limits);
        dSPH_PreStep(fvar, end_ng, pnp1, outlist, npd);
    }

    size_t end = svar.totPts;

#ifndef NOFROZEN
    dissipation_terms(fvar, start, end, outlist, pnp1);
#endif

#ifdef ALE
    if (svar.ghost > 0)
        Particle_Shift_Ghost(svar, fvar, start, end_ng, outlist, pnp1);
    else
        Particle_Shift_No_Ghost(svar, fvar, start, end_ng, outlist, pnp1);
#endif

    Check_Pipe_Outlet(CELL_TREE, svar, avar, cells, limits, pn, pnp1, end, end_ng);

    StateVecD dropVel = StateVecD::Zero();
    StateVecD Force = StateVecD::Zero();

    vector<StateVecD> xih(end - start);
#pragma omp parallel for default(shared)
    for (size_t ii = start; ii < end; ++ii)
    {
        xih[ii - start] = pnp1[ii].xi;
    }

    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    uint iteration = 0;
    solve_prestep(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, limits, outlist, pnp1, pn, xih, start, end,
        logbase, npd, iteration, Force, dropVel, error1, error2
    );

    // cout << "Error: " << error1 << endl;

    /****** UPDATE NEIGHBOURS AND TREE ***********/
    outlist = update_neighbours(SPH_TREE, fvar, pnp1);

    dSPH_PreStep(fvar, end, pnp1, outlist, npd);

    Get_Aero_Velocity(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, vortex, start, end_ng, outlist, limits, pn, pnp1,
        npd
    );

    Detect_Surface(svar, fvar, avar, start, end_ng, outlist, cells, vortex, pnp1);

#ifndef NOFROZEN
    dissipation_terms(fvar, start, end, outlist, pnp1);
#endif

    // Apply_XSPH(fvar,start,end,outlist,pnp1);

#ifdef ALE
    if (svar.ghost > 0)
        Particle_Shift_Ghost(svar, fvar, start, end_ng, outlist, pnp1);
    else
        Particle_Shift_No_Ghost(svar, fvar, start, end_ng, outlist, pnp1);
#endif

    Check_Pipe_Outlet(CELL_TREE, svar, avar, cells, limits, pn, pnp1, end, end_ng);

    dropVel = StateVecD::Zero();

    /*Do time integration*/
    solve_timestep(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, limits, outlist, pnp1, pn, xih, start, end,
        logbase, npd, iteration, Force, dropVel, error1, error2
    );

    if (svar.ghost != 0)
        Check_If_Ghost_Needs_Removing(svar, fvar, SPH_TREE, pn, pnp1);

    uint nAdd = update_buffer_region(svar, limits, pnp1, end, end_ng);

    /* Check if any particles need to be deleted, and turned to particle tracking */
    std::set<size_t> to_del;
    vector<size_t> ndelperblock(limits.size(), 0);
    for (size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
    {
        // Need to check if the delete plane is active.
        if (limits[block].delconst == default_val)
            continue;

#pragma omp parallel for default(shared)
        for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
        {
            // if(pnp1[ii].xi[0] > svar.max_x_sph - fvar.H*fvar.Hfac)
            // {
            /* Create a buffer outlet zone, to avoid truncation of the evolved physics */
            // pnp1[ii].b = OUTLET;
            if (pnp1[ii].xi.dot(limits[block].delete_norm) > limits[block].delconst)
            {
/* SPHPart is downstream enough to convert to IPT */
#pragma omp critical
                {
                    to_del.insert(ii);
                }
#pragma omp atomic
                ndelperblock[block]++;
            }
            // }

            // if(pnp1[ii].b == LOST)
            // {
            // 	to_del.emplace_back(ii);
            // }
        }
    }

    if (svar.using_ipt && svar.Asource != constVel)
    {
        vector<IPTPart> IPT_nm1, IPT_n, IPT_np1;
        for (size_t const& ii : to_del)
        {
            IPT_nm1.emplace_back(IPTPart(pnp1[ii], svar.t, svar.IPT_diam, svar.IPT_area));
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
            svar.totPts--;
            svar.simPts--;
            end_ng--;
            end--;
            nDel++;
            svar.delNum++;
        }

        size_t delshift = 0;
        for (size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
        {
            // Shift the blocks, which will happen cumulatively
            limits[block].index.first -= delshift;
            delshift += ndelperblock[block];
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

    real maxRho = 0.0;
#ifndef DEBUG
    if (svar.speedTest == 0)
#endif
    {
        real maxAf = 0.0;
        real maxRhoi = 0.0;
#ifdef ALE
        real maxShift = 0.0;
#endif
#pragma omp parallel for reduction(max : maxAf, maxRhoi)
        for (size_t ii = start; ii < end_ng; ++ii)
        {
            maxAf = std::max(maxAf, pnp1[ii].Af.norm());
            maxRhoi = std::max(maxRhoi, abs(pnp1[ii].rho - fvar.rho0));
#ifdef ALE
            maxShift = std::max(maxShift, pnp1[ii].vPert.norm());
#endif
        }
        maxRho = 100 * maxRhoi / fvar.rho0;
        real cfl = svar.dt / safe_dt;
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2 - t1).count();
#ifdef ALE
        printf(
            "%9.3e | %7.2e | %4.2f | %8.4f | %7.3f | %9.3e | %9.3e | %9.3e | %14ld|\n", svar.t, svar.dt,
            cfl, error1, maxRho, maxf, maxAf, maxShift, duration
        );

#else
        printf(
            "%9.3e | %7.2e | %4.2f | %9.4f | %3u | %8.3f | %9.3e | %9.3e | %14ld|\n", svar.t, svar.dt,
            cfl, error1, iteration, maxRho, maxf, maxAf, duration
        );
#endif
    }

    if (pnp1.size() == 0 || svar.simPts == 0)
    {
        // cout << "No more particles. Stopping." << endl;
        return 0;
    }

    if (nAdd != 0 || nDel != 0)
    {
        outlist = update_neighbours(SPH_TREE, fvar, pnp1);
    }

    /****** UPDATE TIME N ***********/
    if (svar.totPts != pnp1.size())
    {
        cout << "Size mismatch. Total points not equal to array size. Stopping" << endl;
        exit(-1);
    }

    // cout << "Updating. Time: " << svar.t << "  dt: " << svar.dt << endl;
    if (pn.size() != pnp1.size())
        pn.resize(pnp1.size());

    copy_omp(pnp1, pn);

    /*Add time to global*/
    svar.t += svar.dt;

    /* Check the error and adjust the CFL to try and keep convergence */
    if (error1 > svar.minRes || maxRho > fvar.rhoMaxIter)
    {
        if (error1 > 0.6 * svar.minRes)
        { // If really unstable, immediately reduce the timestep
            svar.cfl = std::max(svar.cfl_min, svar.cfl - svar.cfl_step);
            svar.nUnstable = 0;
        }
        else if (svar.nUnstable > svar.nUnstable_Limit)
        {
            svar.cfl = std::max(svar.cfl_min, svar.cfl - svar.cfl_step);
            svar.nUnstable = 0;
        }
        else
            svar.nUnstable++;
    }
    else
        svar.nUnstable = 0;

    if (iteration < svar.subits_factor * svar.subits && svar.nUnstable == 0)
    {
        if (svar.nStable > svar.nStable_Limit)
        { // Try boosting the CFL if it does converge nicely
            svar.cfl = std::min(svar.cfl_max, svar.cfl + svar.cfl_step);
            svar.nStable = 0;
        }
        else
            svar.nStable++;
    }
    else
        svar.nStable = 0;

#ifdef DAMBREAK
    /* Find max x and y coordinates of the fluid */
    // Find maximum safe timestep
    vector<SPHPart>::iterator maxXi = std::max_element(
        pnp1.begin() + svar.bndPts, pnp1.end(),
        [](SPHPart const& p1, SPHPart const& p2) { return p1.xi[0] < p2.xi[0]; }
    );

    vector<SPHPart>::iterator maxYi = std::max_element(
        pnp1.begin() + svar.bndPts, pnp1.end(),
        [](SPHPart const& p1, SPHPart const& p2) { return p1.xi[1] < p2.xi[1]; }
    );

    real maxX = maxXi->xi[0];
    real maxY = maxYi->xi[1];

    dambreak << svar.t << "  " << maxX + 0.5 * svar.Pstep << "  " << maxY + 0.5 * svar.Pstep << endl;
#endif

    return error1;
}

/* Integration loop to perform the first timestep and find forces only.

Solve the forces and pressures for particles according to a calculated stable timestep (current no option
to use a fixed timestep) and do not update the new timestep. This should be used only as a prestep to
write files with starting forces.
*/
void Integrator::first_step(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    VLM const& vortex, MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1,
    SPHState& pn, vector<IPTState>& iptdata
)
{
    size_t const start = svar.bndPts;
    size_t end = svar.totPts;
    size_t end_ng = svar.bndPts + svar.simPts;
    real npd = 1.0;
#if DEBUG
    fprintf(dbout, "Starting first step. Start index: %zu End index: %zu\n", start, end);
#endif

    dSPH_PreStep(fvar, end, pnp1, outlist, npd);
    if (svar.Asource == meshInfl)
    {
        FindCell(svar, avar, CELL_TREE, cells, pnp1, pn);
        if (svar.totPts != pnp1.size())
        { // Rebuild the neighbour list
            // cout << "Updating neighbour list" << endl;
            // cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
            svar.delNum += svar.totPts - pnp1.size();
            svar.totPts = pnp1.size();
            svar.simPts = svar.totPts - svar.bndPts;
            end = svar.totPts;
            end_ng = end;
            outlist = update_neighbours(SPH_TREE, fvar, pnp1);
            dSPH_PreStep(fvar, end, pnp1, outlist, npd);
        }
    }

    Detect_Surface(svar, fvar, avar, start, end_ng, outlist, cells, vortex, pnp1);

    if (svar.ghost == 1)
    { /* Poisson points */
        PoissonGhost(svar, fvar, avar, cells, SPH_TREE, outlist, pn, pnp1);
        dSPH_PreStep(fvar, svar.totPts, pnp1, outlist, npd);
        // Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,pnp1);
    }
    else if (svar.ghost == 2)
    { /* Lattice points */
        LatticeGhost(svar, fvar, avar, cells, SPH_TREE, outlist, pn, pnp1, limits);
        dSPH_PreStep(fvar, svar.totPts, pnp1, outlist, npd);
        // Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,pnp1);
    }

    end = svar.totPts;

    uint iteration = 0;
    real error1 = 1.0;
    real error2 = 0.0;
    real logbase = 0.0;

    StateVecD dropVel = StateVecD::Zero();
    StateVecD Force = StateVecD::Zero();

    /*find force at time n*/
    vector<StateVecD> res;
    vector<StateVecD> Af;
    vector<real> Rrho;
    vector<real> curve;

    // Find maximum safe timestep
    svar.dt = find_timestep(svar, fvar, cells, pnp1, start, end_ng);

#ifndef NOFROZEN
    dissipation_terms(fvar, start, end, outlist, pnp1);
#endif

    // /*Previous SPHState for error calc*/
    vector<StateVecD> xih(end - start);
#pragma omp parallel for default(shared) // shared(pnp1)
    for (size_t ii = start; ii < end; ++ii)
    {
        xih[ii - start] = pnp1[ii].xi;
    }

    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    solve_prestep(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, limits, outlist, pnp1, pn, xih, start, end,
        logbase, npd, iteration, Force, dropVel, error1, error2
    );

    outlist = update_neighbours(SPH_TREE, fvar, pnp1);

    dSPH_PreStep(fvar, end, pnp1, outlist, npd);

    if (svar.Asource == meshInfl)
    {
        FindCell(svar, avar, CELL_TREE, cells, pnp1, pn);
        if (svar.totPts != pnp1.size())
        { // Rebuild the neighbour list
            svar.delNum += svar.totPts - pnp1.size();
            svar.totPts = pnp1.size();
            svar.simPts = svar.totPts - svar.bndPts;
            end = svar.totPts;
            end_ng = end;
            outlist = update_neighbours(SPH_TREE, fvar, pnp1);
            dSPH_PreStep(fvar, end, pnp1, outlist, npd);
        }
    }

    Detect_Surface(svar, fvar, avar, start, end_ng, outlist, cells, vortex, pnp1);

#ifndef NOFROZEN
    dissipation_terms(fvar, start, end, outlist, pnp1);
#endif

    /*Do time integration*/
    solve_timestep(
        SPH_TREE, CELL_TREE, svar, fvar, avar, cells, limits, outlist, pnp1, pn, xih, start, end,
        logbase, npd, iteration, Force, dropVel, error1, error2
    );

    Check_If_Ghost_Needs_Removing(svar, fvar, SPH_TREE, pn, pnp1);

    // In the first timestep no update of the particle positions is desired.
    // It is only to establish forces and pressures prior to writing time zero files.
    // Useful to debug if starting forces are crazy
#if DEBUG
    fprintf(dbout, "Exiting first step. Error: %f\n", error1);
#endif
}

/* Internal pre-step for the integration methods to find base error for the timestep. */
void Integrator::solve_prestep(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1, SPHState& pn,
    vector<StateVecD>& xih, size_t const& start, size_t& end, real& logbase, real& npd, uint& iteration,
    StateVecD& Force, StateVecD& dropVel, real& error1, real& error2
)
{
    /*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
    switch (solver_method)
    {
    case runge_kutta:
    {
        error1 = Get_First_RK(
            svar, fvar, avar, start, end, fvar.B, fvar.gam, npd, cells, limits, outlist, logbase, pn,
            pnp1, error1
        );
        break;
    }
    case newmark_beta:
    {
        Newmark_Beta::Do_NB_Iter(
            CELL_TREE, svar, fvar, avar, start, end, npd, cells, limits, outlist, pn, pnp1, Force,
            dropVel
        );

        void(Newmark_Beta::Check_Error(
            SPH_TREE, svar, fvar, start, end, error1, error2, logbase, outlist, xih, pn, pnp1, iteration
        ));
        iteration++; // Update iteration count
        break;
    }
    }
}

/* Internal funtion to solve the timestep given a base error from the pre-step. */
void Integrator::solve_timestep(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1, SPHState& pn,
    vector<StateVecD>& xih, size_t const& start, size_t& end, real& logbase, real& npd, uint& iteration,
    StateVecD& Force, StateVecD& dropVel, real& error1, real& error2
)
{
    /*Do time integration*/
    switch (solver_method)
    {
    case runge_kutta:
    {
        SPHState st_2 = pnp1;

        error1 = Runge_Kutta4(
            CELL_TREE, svar, fvar, avar, start, end, fvar.B, fvar.gam, npd, cells, limits, outlist,
            logbase, pn, st_2, pnp1, Force, dropVel
        );
        break;
    }
    case newmark_beta:
    {

        Newmark_Beta::Newmark_Beta(
            SPH_TREE, CELL_TREE, svar, fvar, avar, start, end, npd, cells, limits, outlist, logbase,
            iteration, error1, error2, xih, pn, pnp1, Force, dropVel
        );
        break;
    }
    }
}

real Integrator::find_timestep(
    SIM const& svar, FLUID const& fvar, MESH const& cells, SPHState const& pnp1, size_t const& start,
    size_t const& end_ng
)
{
    // Find maximum safe timestep (avoiding a div by zero)
    real maxf = MEPSILON;
    real maxdrho = MEPSILON;
    real maxU = MEPSILON;
    real minST = 9999999.0;
#pragma omp parallel for reduction(max : maxf, maxU, maxdrho) reduction(min : minST)
    for (size_t ii = start; ii < end_ng; ++ii)
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
        0.125 * fvar.HSQ * fvar.rho0 / fvar.mu
    );                                                           /* Viscosity timestep constraint */
    timestep_factors.emplace_back(0.067 * minST);                /* Surface tension constraint */
    timestep_factors.emplace_back(0.5 * sqrt(fvar.H / maxdrho)); /* Density gradient constraint */
    timestep_factors.emplace_back(1.5 * fvar.H / fvar.Cs);
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
    real safe_dt = 0.75 * *std::min_element(timestep_factors.begin(), timestep_factors.end());
    real dt = svar.cfl * safe_dt;
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/

    if (dt < svar.dt_min)
        dt = svar.dt_min;
    else if (dt > svar.dt_max)
        dt = svar.dt_max;

    if (dt > svar.tframem1 + svar.framet - svar.t)
        dt = svar.tframem1 + svar.framet - svar.t + svar.dt_min;

    if (svar.speedTest)
    {
        // Bound SPH timestep by the mesh size if it's larger than that to traverse a cell, to prevent
        // skipping.
        if (dt * pnp1[0].v.norm() > cells.minlength)
        {
            dt = cells.minlength / pnp1[0].v.norm();
        }
    }

    return dt;
}
