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

    real rms_error = log10(sqrt(errsum / (real(svar.totPts)))) - logbase;

    return rms_error;
}

SPHState do_runge_kutta_intermediate_step(
    SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t const& end,
    real const& npd, MESH const& cells, LIMITS const& limits, OUTL const& outlist,
    SPHState const& part_n, SPHState const& part_prev, real const& dt_inter
)
{
    SPHState part_inter = part_prev;
    vector<StateVecD> res_inter(end, StateVecD::Zero());
    vector<real> Rrho_inter(end, 0.0);

    vector<StateVecD> Af;

    for (size_t block = 0; block < svar.nbound; block++)
    {

        if (limits[block].nTimes != 0)
        {
            // Get the current boundary velocity
            StateVecD vel = StateVecD::Zero();
            for (size_t time = 0; time < limits[block].nTimes; time++)
            {
                if (limits[block].times[time] > svar.t)
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
                fvar, limits[block].index.first, limits[block].index.second, outlist, part_inter
            );

        switch (limits[block].bound_solver)
        {
        case DBC:
        {
            Boundary_DBC(
                fvar, limits[block].index.first, limits[block].index.second, outlist, part_inter,
                res_inter
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                real const rho = std::max(
                    fvar.rhoMin, std::min(fvar.rhoMax, part_n[ii].rho + dt_inter * Rrho_inter[ii])
                );
                part_inter[ii].rho = rho;
                part_inter[ii].p =
                    pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                part_inter[ii].Rrho = Rrho_inter[ii];
            }
            break;
        }
        case pressure_G:
        {
            Get_Boundary_Pressure(
                svar.grav, fvar, limits[block].index.first, limits[block].index.second, outlist,
                part_inter
            );
            break;
        }
        case ghost:
        {
            vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
            Boundary_Ghost(
                fvar, limits[block].index.first, limits[block].index.second, outlist, part_inter,
                Rrho_inter, near_inlet
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                if (!near_inlet[ii - limits[block].index.first])
                {
                    real const rho = std::max(
                        fvar.rhoMin, std::min(fvar.rhoMax, part_n[ii].rho + dt_inter * Rrho_inter[ii])
                    );
                    part_inter[ii].rho = rho;
                    part_inter[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                    part_inter[ii].Rrho = Rrho_inter[ii];
                }
                else
                { // Don't allow negative pressures
                    real const rho = std::max(
                        fvar.rho0, std::min(fvar.rhoMax, part_n[ii].rho + dt_inter * Rrho_inter[ii])
                    );
                    part_inter[ii].rho = rho;
                    part_inter[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                    part_inter[ii].Rrho = fmax(0.0, Rrho_inter[ii]);
                }
            }
            break;
        }
        }
    }

    Forces(svar, fvar, avar, cells, part_inter, outlist, npd, res_inter, Rrho_inter, Af);

#pragma omp parallel default(shared) // shared(svar,part_n,st_2) /*reduction(+:Force,dropVel)*/
    {
        /********************************************************************/
        /*************************  STEP 1  *********************************/
        /********************************************************************/
        for (size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
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

                    part_inter[ii].v = part_n[ii].v + dt_inter * res_inter[ii];
                    real const rho = std::max(
                        fvar.rhoMin, std::min(fvar.rhoMax, part_n[ii].rho + dt_inter * Rrho_inter[ii])
                    );
                    part_inter[ii].rho = rho;
                    part_inter[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                    part_inter[ii].acc = res_inter[ii];
                    part_inter[ii].Af = Af[ii];
                    part_inter[ii].Rrho = Rrho_inter[ii];
                }
                else
                {
                    part_inter[ii].xi = part_n[ii].xi + dt_inter * part_inter[ii].v;
                    part_inter[ii].rho = part_n[ii].rho + 0.5 * dt_inter * Rrho_inter[ii];
                    part_inter[ii].p = pressure_equation(
                        part_inter[ii].rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP
                    );
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

                            part_inter[buffID].Rrho = Rrho_inter[buffID];

                            real const rho = std::max(
                                fvar.rhoMin,
                                std::min(fvar.rhoMax, part_n[buffID].rho + dt_inter * Rrho_inter[buffID])
                            );
                            part_inter[buffID].rho = rho;
                            part_inter[buffID].p =
                                pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
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
    SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t& end, real const& npd,
    MESH const& cells, LIMITS const& limits, OUTL const& outlist, real const& dt, SPHState const& part_n,
    SPHState const& st_1, SPHState const& st_2, SPHState const& st_3
)
{
    SPHState part_np1 = st_3;
    vector<StateVecD> res_4(end, StateVecD::Zero());
    vector<real> Rrho_4(end, 0.0);
    vector<StateVecD> Af(end, StateVecD::Zero());

    for (size_t block = 0; block < svar.nbound; block++)
    {
        if (limits[block].nTimes != 0)
        {
            // Get the current boundary velocity
            StateVecD vel = StateVecD::Zero();
            for (size_t time = 0; time < limits[block].nTimes; time++)
            {
                if (limits[block].times[time] > svar.t)
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
            Set_No_Slip(fvar, limits[block].index.first, limits[block].index.second, outlist, part_np1);

        switch (limits[block].bound_solver)
        {
        case DBC:
        {
            Boundary_DBC(
                fvar, limits[block].index.first, limits[block].index.second, outlist, part_np1, res_4
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                real const rho = std::max(
                    fvar.rhoMin,
                    std::min(
                        fvar.rhoMax,
                        part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                       2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                    )
                );
                part_np1[ii].rho = rho;
                part_np1[ii].p =
                    pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                part_np1[ii].Rrho = st_3[ii].Rrho;
            }
            break;
        }
        case pressure_G:
        {
            Get_Boundary_Pressure(
                svar.grav, fvar, limits[block].index.first, limits[block].index.second, outlist, part_np1
            );
            break;
        }
        case ghost:
        {
            vector<int> near_inlet(limits[block].index.second - limits[block].index.first, 0);
            Boundary_Ghost(
                fvar, limits[block].index.first, limits[block].index.second, outlist, part_np1, Rrho_4,
                near_inlet
            );

#pragma omp for schedule(static) nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** BOUNDARY PARTICLES ***********/
                if (!near_inlet[ii - limits[block].index.first])
                {
                    real const rho = std::max(
                        fvar.rhoMin,
                        std::min(
                            fvar.rhoMax,
                            part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                           2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                        )
                    );
                    part_np1[ii].rho = rho;
                    part_np1[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                    part_np1[ii].Rrho = Rrho_4[ii];
                }
                else
                { // Don't allow negative pressures
                    real const rho = std::max(
                        fvar.rho0,
                        std::min(
                            fvar.rhoMax,
                            part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                           2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                        )
                    );
                    part_np1[ii].rho = rho;
                    part_np1[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                    part_np1[ii].Rrho = fmax(0.0, Rrho_4[ii]);
                }
            }
            break;
        }
        }
    }

    Forces(svar, fvar, avar, cells, part_np1, outlist, npd, res_4, Rrho_4, Af);

#pragma omp parallel default(shared)
    {
        for (size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
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
                        fvar.rhoMin,
                        std::min(
                            fvar.rhoMax,
                            part_n[ii].rho + (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                           2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                        )
                    );
                    part_np1[ii].rho = rho;

                    part_np1[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);

                    part_np1[ii].acc = res_4[ii];
                    part_np1[ii].Af = Af[ii];
                    part_np1[ii].Rrho = Rrho_4[ii];
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
                                fvar.rhoMin,
                                std::min(
                                    fvar.rhoMax,
                                    part_n[buffID].rho +
                                        (dt / 6.0) * (part_n[ii].Rrho + 2.0 * st_1[ii].Rrho +
                                                      2.0 * st_2[ii].Rrho + st_3[ii].Rrho)
                                )
                            );
                            part_np1[buffID].rho = rho;
                            part_np1[buffID].p =
                                pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
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
    SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t const& end,
    real const& npd, MESH const& cells, LIMITS const& limits, OUTL const& outlist, SPHState& part_n,
    SPHState& st_1
)
{
    const real dt = 0.5 * svar.dt;

    st_1 = do_runge_kutta_intermediate_step(
        svar, fvar, avar, start, end, npd, cells, limits, outlist, part_n, part_n, dt
    );

    // First iteration so error log base is zero.
    real logbase = 0.0;
    return Check_RK_Error(svar, start, end, logbase, part_n, st_1);
}

/* <summary> Perform the rest of the Runge-Kutta integration, assuming frozen
        dissipation terms. This will do step 2 to 4 (n+1/4 to n+1) </summary> */
real Runge_Kutta4(
    SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t& end, real const& npd,
    MESH const& cells, LIMITS const& limits, OUTL const& outlist, real const& logbase, SPHState& part_n,
    SPHState& st_1, SPHState& part_np1
)
{
    /*Create the vectors*/
    SPHState st_2 = part_np1; /*Step 2*/
    SPHState st_3 = part_np1; /*Step 3*/

    real const dt_inter = 0.5 * svar.dt;
    real const dt_final = svar.dt;

    /********************************************************************/
    /*************************  STEP 2  *********************************/
    /********************************************************************/
    st_2 = do_runge_kutta_intermediate_step(
        svar, fvar, avar, start, end, npd, cells, limits, outlist, part_n, st_1, dt_inter
    );

    /********************************************************************/
    /*************************  STEP 3  *********************************/
    /********************************************************************/
    st_3 = do_runge_kutta_intermediate_step(
        svar, fvar, avar, start, end, npd, cells, limits, outlist, part_n, st_2, dt_final
    );

    /********************************************************************/
    /*************************  STEP 4  *********************************/
    /********************************************************************/
    part_np1 = do_runge_kutta_final_step(
        svar, fvar, avar, start, end, npd, cells, limits, outlist, dt_final, part_n, st_1, st_2, st_3
    );

    return Check_RK_Error(svar, start, end, logbase, st_3, part_np1);
}
