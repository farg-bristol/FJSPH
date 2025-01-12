/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Runge_Kutta.h"
#include "Containment.h"
#include "Kernel.h"
#include "Resid.h"
#include "Shifting.h"

real Check_RK_Error(
    SIM const& svar, size_t const& start, size_t const& end, real const& logbase, SPHState const& part_n,
    SPHState const& part_np1
)
{
    /****** FIND ERROR ***********/
    real errsum = 0.0;
#pragma omp parallel for reduction(+ : errsum)
    for (size_t ii = start; ii < end; ++ii)
    {
        errsum += (part_np1[ii].xi - part_n[ii].xi).squaredNorm();
    }

    real rms_error = log10(sqrt(errsum / (real(end - start)))) - logbase;

    return rms_error;
}

SPHState do_runge_kutta_intermediate_step(
    SIM& svar, size_t const& start, size_t const& end, real const& npd, MESH const& cells,
    LIMITS const& limits, OUTL const& outlist, SPHState const& part_n, SPHState const& part_prev,
    real const& dt_inter
)
{
    SPHState part_inter = part_prev;

    for (size_t block = 0; block < svar.n_bound_blocks; block++)
    {

        if (limits[block].nTimes != 0)
        {
            // Get the current boundary velocity
            StateVecD vel = StateVecD::Zero();
            for (size_t time = 0; time < limits[block].nTimes; time++)
            {
                if (limits[block].times[time] > svar.integrator.current_time)
                {
                    vel = limits[block].vels[time];
                }
            }

            for (size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
            {
                part_inter[jj].v = vel;
            }
        }
        else
        {
#pragma omp parallel for
            for (size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
            {
                part_inter[jj].v = limits[block].vels[0];
            }
        }

        if (limits[block].no_slip)
            Set_No_Slip(
                limits[block].index.first, limits[block].index.second, outlist, svar.fluid.H,
                svar.fluid.W_correc, part_inter
            );

        switch (limits[block].bound_solver)
        {
        case DBC:
        {
            Boundary_DBC(
                svar.fluid, limits[block].index.first, limits[block].index.second, outlist, part_inter
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                real const rho = std::max(
                    svar.fluid.rho_min,
                    std::min(svar.fluid.rho_max, part_n[ii].rho + dt_inter * part_inter[ii].Rrho)
                );
                part_inter[ii].rho = rho;
                part_inter[ii].p = svar.fluid.get_pressure(rho);
            }
            break;
        }
        case pressure_G:
        {
            Get_Boundary_Pressure(
                svar.fluid, svar.grav, limits[block].index.first, limits[block].index.second, outlist,
                part_inter
            );
            break;
        }
        case ghost:
        {
            vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
            Boundary_Ghost(
                limits[block].index.first, limits[block].index.second, outlist, svar.fluid.H,
                svar.fluid.W_correc, part_inter, near_inlet
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                if (!near_inlet[ii - limits[block].index.first])
                {
                    real const rho = std::max(
                        svar.fluid.rho_min,
                        std::min(svar.fluid.rho_max, part_n[ii].rho + dt_inter * part_inter[ii].Rrho)
                    );
                    part_inter[ii].rho = rho;
                    part_inter[ii].p = svar.fluid.get_pressure(rho);
                }
                else
                { // Don't allow negative pressures
                    real const rho = std::max(
                        svar.fluid.rho_rest,
                        std::min(svar.fluid.rho_max, part_n[ii].rho + dt_inter * part_inter[ii].Rrho)
                    );
                    part_inter[ii].rho = rho;
                    part_inter[ii].p = svar.fluid.get_pressure(rho);
                    part_inter[ii].Rrho = fmax(0.0, part_inter[ii].Rrho);
                }
            }
            break;
        }
        }
    }

    get_acc_and_Rrho(svar, cells, outlist, npd, part_inter);

#pragma omp parallel default(shared) // shared(svar,part_n,st_2) /*reduction(+:Force,dropVel)*/
    {
        /********************************************************************/
        /*************************  STEP 1  *********************************/
        /********************************************************************/
        for (size_t block = svar.n_bound_blocks; block < svar.n_fluid_blocks + svar.n_bound_blocks;
             block++)
        {
#pragma omp for nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** FLUID PARTICLES **************/

                /* BUFFER = pipe particle receiving prescribed motion            */
                /* BACK = the latest particle in that column                    */
                /* PIPE = in the pipe, with free motion                         */
                /* FREE = free of the pipe and receives an aero force           */

                if (part_inter[ii].b > BUFFER)
                {

#ifdef ALE
                    part_inter[ii].xi =
                        part_n[ii].xi + dt_inter * (part_inter[ii].v + part_inter[ii].vPert);
#else
                    part_inter[ii].xi = part_n[ii].xi + dt_inter * part_inter[ii].v;
#endif

                    part_inter[ii].v = part_n[ii].v + dt_inter * part_inter[ii].acc;
                    real const rho = std::max(
                        svar.fluid.rho_min,
                        std::min(svar.fluid.rho_max, part_n[ii].rho + dt_inter * part_inter[ii].Rrho)
                    );
                    part_inter[ii].rho = rho;
                    part_inter[ii].p = svar.fluid.get_pressure(rho);
                }
            }

            /* Do the buffer particles */
            if (limits[block].block_type == inletZone)
            {
                switch (limits[block].fixed_vel_or_dynamic)
                {
                case 1:
                {
                    StateVecD unorm = limits[block].insert_norm.normalized();
#pragma omp for schedule(static) nowait
                    for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
                    {
                        size_t const& backID = limits[block].back[ii];
                        StateVecD const& xi = part_inter[backID].xi;

                        for (size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
                        { /* Define buffer off the back particle */
                            size_t const& buffID = limits[block].buffer[ii][jj];
                            // Set position as related to the previous particle.
                            part_inter[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

                            // How to set density and pressure though?
                            part_inter[buffID].v = part_inter[backID].v;
                            part_inter[buffID].rho = part_inter[backID].rho;
                            part_inter[buffID].p = part_inter[backID].p;
                        }
                    }
                    break;
                }
                case 0:
                {
#pragma omp for schedule(static) nowait
                    for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
                    {
                        // size_t const& backID = svar.back[ii];
                        for (size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
                        { /* Define buffer off the back particle */

                            size_t const& buffID = limits[jj].buffer[ii][jj];

                            real const rho = std::max(
                                svar.fluid.rho_min,
                                std::min(
                                    svar.fluid.rho_max,
                                    part_n[buffID].rho + dt_inter * part_inter[buffID].Rrho
                                )
                            );
                            part_inter[buffID].rho = rho;
                            part_inter[buffID].p = svar.fluid.get_pressure(rho);
                            part_inter[buffID].xi = part_n[buffID].xi + dt_inter * part_n[buffID].v;
                        }
                    }
                    break;
                }
                }
            }
        }
    }
    return part_inter;
}

SPHState do_runge_kutta_final_step(
    SIM& svar, size_t const& start, size_t& end, real const& npd, MESH const& cells,
    LIMITS const& limits, OUTL const& outlist, real const& dt, SPHState const& part_n,
    SPHState const& st_1, SPHState const& st_2, SPHState const& st_3
)
{
    SPHState part_np1 = st_3;

    for (size_t block = 0; block < svar.n_bound_blocks; block++)
    {
        if (limits[block].nTimes != 0)
        {
            // Get the current boundary velocity
            StateVecD vel = StateVecD::Zero();
            for (size_t time = 0; time < limits[block].nTimes; time++)
            {
                if (limits[block].times[time] > svar.integrator.current_time)
                {
                    vel = limits[block].vels[time];
                }
            }

            for (size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
            {
                part_np1[jj].v = vel;
            }
        }
        else
        {
#pragma omp parallel for
            for (size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
            {
                part_np1[jj].v = limits[block].vels[0];
            }
        }

        if (limits[block].no_slip)
            Set_No_Slip(
                limits[block].index.first, limits[block].index.second, outlist, svar.fluid.H,
                svar.fluid.W_correc, part_np1
            );

        switch (limits[block].bound_solver)
        {
        case DBC:
        {
            Boundary_DBC(
                svar.fluid, limits[block].index.first, limits[block].index.second, outlist, part_np1
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                real const rho = std::max(
                    svar.fluid.rho_min,
                    std::min(
                        svar.fluid.rho_max,
                        part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                       2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                    )
                );
                part_np1[ii].rho = rho;
                part_np1[ii].p = svar.fluid.get_pressure(rho);
                part_np1[ii].Rrho = st_3[ii].Rrho;
            }
            break;
        }
        case pressure_G:
        {
            Get_Boundary_Pressure(
                svar.fluid, svar.grav, limits[block].index.first, limits[block].index.second, outlist,
                part_np1
            );
            break;
        }
        case ghost:
        {
            vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
            Boundary_Ghost(
                limits[block].index.first, limits[block].index.second, outlist, svar.fluid.H,
                svar.fluid.W_correc, part_np1, near_inlet
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                if (!near_inlet[ii - limits[block].index.first])
                {
                    real const rho = std::max(
                        svar.fluid.rho_min,
                        std::min(
                            svar.fluid.rho_max,
                            part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                           2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                        )
                    );
                    part_np1[ii].rho = rho;
                    part_np1[ii].p = svar.fluid.get_pressure(rho);
                }
                else
                { // Don't allow negative pressures
                    real const rho = std::max(
                        svar.fluid.rho_rest,
                        std::min(
                            svar.fluid.rho_max,
                            part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                           2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                        )
                    );
                    part_np1[ii].rho = rho;
                    part_np1[ii].p = svar.fluid.get_pressure(rho);
                    part_np1[ii].Rrho = fmax(0.0, part_np1[ii].Rrho);
                }
            }
            break;
        }
        }
    }

    get_acc_and_Rrho(svar, cells, outlist, npd, part_np1);

#pragma omp parallel default(shared)
    {
        for (size_t block = svar.n_bound_blocks; block < svar.n_fluid_blocks + svar.n_bound_blocks;
             block++)
        {
#pragma omp for nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** FLUID PARTICLES **************/
                if (part_np1[ii].b > BUFFER && part_np1[ii].b != OUTLET)
                {
#ifdef ALE
                    part_np1[ii].xi = part_n[ii].xi + (dt / 6.0) * ((part_n[ii].v + part_n[ii].vPert) +
                                                                    2.0 * (st_1[ii].v + st_1[ii].vPert) +
                                                                    2.0 * (st_2[ii].v + st_2[ii].vPert) +
                                                                    (st_3[ii].v + st_3[ii].vPert));
#else
                    part_np1[ii].xi = part_n[ii].xi + (dt / 6.0) * (part_n[ii].v + 2.0 * st_1[ii].v +
                                                                    2.0 * st_2[ii].v + st_3[ii].v);
#endif

                    part_np1[ii].v = part_n[ii].v + (dt / 6.0) * (part_n[ii].acc + 2.0 * st_1[ii].acc +
                                                                  2.0 * st_2[ii].acc + st_3[ii].acc);

                    real const rho = std::max(
                        svar.fluid.rho_min,
                        std::min(
                            svar.fluid.rho_max,
                            part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                           2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                        )
                    );
                    part_np1[ii].rho = rho;

                    part_np1[ii].p = svar.fluid.get_pressure(rho);
                }
                else if (part_np1[ii].b == OUTLET)
                { /* For the outlet zone, just perform euler integration of last info */
                    part_np1[ii].xi = part_n[ii].xi + dt * part_np1[ii].v;
                }
            }

            /* Do the buffer particles */
            if (limits[block].block_type == inletZone)
            {
                switch (limits[block].fixed_vel_or_dynamic)
                {
                case 1:
                {
                    StateVecD unorm = limits[block].insert_norm.normalized();
#pragma omp for schedule(static) nowait
                    for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
                    {
                        size_t const& backID = limits[block].back[ii];
                        StateVecD const& xi = part_np1[backID].xi;

                        for (size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
                        { /* Define buffer off the back particle */
                            size_t const& buffID = limits[block].buffer[ii][jj];
                            // Set position as related to the previous particle.
                            part_np1[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

                            // How to set density and pressure though?
                            part_np1[buffID].v = part_np1[backID].v;
                            part_np1[buffID].rho = part_np1[backID].rho;
                            part_np1[buffID].p = part_np1[backID].p;
                        }
                    }
                    break;
                }
                case 0:
                {
#pragma omp for schedule(static) nowait
                    for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
                    {
                        // size_t const& backID = svar.back[ii];
                        for (size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
                        { /* Define buffer off the back particle */

                            size_t const& buffID = limits[jj].buffer[ii][jj];

                            real const rho = std::max(
                                svar.fluid.rho_min,
                                std::min(
                                    svar.fluid.rho_max,
                                    part_n[buffID].rho +
                                        (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                      2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                                )
                            );
                            part_np1[buffID].rho = rho;
                            part_np1[buffID].p = svar.fluid.get_pressure(rho);
                            part_np1[buffID].xi = part_n[buffID].xi + dt * part_n[buffID].v;
                        }
                    }
                    break;
                }
                }
            }
        }
    } /*End pragma omp parallel*/

    return part_np1;
}

/* Peform the first stage of the Runge-Kutta integration
to get the first guess of time n+1 (regarded here as time n+1/4)
to perform neighbour search and dissipation terms before freezing </summary */
real Get_First_RK(
    SIM& svar, size_t const& start, size_t const& end, real const& npd, MESH const& cells,
    LIMITS const& limits, OUTL const& outlist, SPHState& part_n, SPHState& st_1
)
{
    const real dt = 0.5 * svar.integrator.delta_t;

    st_1 = do_runge_kutta_intermediate_step(
        svar, start, end, npd, cells, limits, outlist, part_n, part_n, dt
    );

    // First iteration so error log base is zero.
    real logbase = 0.0;
    return Check_RK_Error(svar, start, end, logbase, part_n, st_1);
}

/* <summary> Perform the rest of the Runge-Kutta integration, assuming frozen
        dissipation terms. This will do step 2 to 4 (n+1/4 to n+1) </summary> */
real Runge_Kutta4(
    SIM& svar, size_t const& start, size_t& end, real const& npd, MESH const& cells,
    LIMITS const& limits, OUTL const& outlist, real const& logbase, SPHState& part_n, SPHState& st_1,
    SPHState& part_np1
)
{
    /*Create the vectors*/
    SPHState st_2 = part_np1; /*Step 2*/
    SPHState st_3 = part_np1; /*Step 3*/

    real const dt_inter = 0.5 * svar.integrator.delta_t;
    real const dt_final = svar.integrator.delta_t;

    /********************************************************************/
    /*************************  STEP 2  *********************************/
    /********************************************************************/
    st_2 = do_runge_kutta_intermediate_step(
        svar, start, end, npd, cells, limits, outlist, part_n, st_1, dt_inter
    );

    /********************************************************************/
    /*************************  STEP 3  *********************************/
    /********************************************************************/
    st_3 = do_runge_kutta_intermediate_step(
        svar, start, end, npd, cells, limits, outlist, part_n, st_2, dt_final
    );

    /********************************************************************/
    /*************************  STEP 4  *********************************/
    /********************************************************************/
    part_np1 = do_runge_kutta_final_step(
        svar, start, end, npd, cells, limits, outlist, dt_final, part_n, st_1, st_2, st_3
    );

    return Check_RK_Error(svar, start, end, logbase, st_3, part_np1);
}
