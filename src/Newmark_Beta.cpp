/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Newmark_Beta.h"
#include "Containment.h"
#include "Kernel.h"
#include "Neighbours.h"
#include "Resid.h"

int Newmark_Beta::Check_Error(
    Sim_Tree& SPH_TREE, SIM& svar, FLUID const& fvar, size_t const& start, size_t const& end,
    real& rms_error, real& logbase, OUTL& outlist, vector<StateVecD> const& xih, SPHState const& pn,
    SPHState& pnp1, uint& iteration
)
{
    /****** FIND ERROR ***********/
    real errsum = 0.0;
#pragma omp parallel for reduction(+ : errsum) schedule(static) default(shared)
    for (size_t ii = start; ii < end; ++ii)
    {
        StateVecD r = pnp1[ii].xi - xih[ii - start];
        errsum += r.squaredNorm();
    }

    real log_error = log10(sqrt(errsum / real(svar.totPts)));

    if (iteration == 0)
        logbase = log_error;

    rms_error = log_error - logbase;

    if (iteration > svar.subits)
    {
        if (rms_error > 0.0)
        {
            pnp1 = pn;

            outlist = update_neighbours(SPH_TREE, fvar, pnp1);

            svar.dt = 0.5 * svar.dt;
            cout << "Unstable timestep. New dt: " << svar.dt << endl;
            iteration = 0;
            rms_error = 0.0;
            // RestartCount++;
            return -1;
        }
        return 1;
    }

    /*Otherwise, roll forwards*/
    return 0;
}

void Newmark_Beta::Do_NB_Iter(
    Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start,
    size_t& end, real const& npd, MESH const& cells, LIMITS const& limits, OUTL& outlist,
    SPHState const& pn, SPHState& pnp1
)
{
    /*Update the state at time n+1*/

    // cout << "Calculating forces" << endl;
    vector<StateVecD> res(end, StateVecD::Zero());
    vector<StateVecD> Af(end, StateVecD::Zero());
    vector<real> Rrho(end, 0.0);
    vector<int> near_inlet(svar.bndPts, 0);

    real const& gamma_t1 = svar.gamma;
    real const gamma_t2 = 1 - gamma_t1;
    real const& beta_t1 = svar.beta;
    real const beta_t2 = 0.5 * (1 - 2 * beta_t1);

    for (size_t block = 0; block < svar.nbound; block++)
    {
        if (limits[block].nTimes != 0)
        {
            // Get the current boundary velocity
            StateVecD vel = StateVecD::Zero();
            for (size_t time = 0; time < limits[block].nTimes; time++)
            {
                if (svar.t > limits[block].times[time])
                {
                    vel = limits[block].vels[time];
                }
            }

#pragma omp parallel for
            for (size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
            {
                pnp1[jj].v = vel;
            }
        }
        else
        {
#pragma omp parallel for
            for (size_t jj = limits[block].index.first; jj < limits[block].index.second; jj++)
            {
                pnp1[jj].v = limits[block].vels[0];
            }
        }

        if (limits[block].no_slip)
            Set_No_Slip(fvar, limits[block].index.first, limits[block].index.second, outlist, pnp1);

        switch (limits[block].bound_solver)
        {
        case DBC:
        {
            Boundary_DBC(
                fvar, limits[block].index.first, limits[block].index.second, outlist, pnp1, res
            );
            break;
        }
        case pressure_G:
        {
            Get_Boundary_Pressure(
                svar.grav, fvar, limits[block].index.first, limits[block].index.second, outlist, pnp1
            );
            break;
        }
        case ghost:
        {
            Boundary_Ghost(
                fvar, limits[block].index.first, limits[block].index.second, outlist, pnp1, Rrho,
                near_inlet
            );
            break;
        }
        default:
            break;
        }
    }

    /*Guess force at time n+1*/
    Forces(svar, fvar, avar, cells, pnp1, outlist, npd, res, Rrho, Af);

#pragma omp parallel default(shared)
    {
        const real dt = svar.dt;
        const real dt2 = dt * dt;

        for (size_t block = 0; block < svar.nbound; block++)
        {
            switch (limits[block].bound_solver)
            {
            case DBC:
            {
#pragma omp for schedule(static) nowait
                for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
                { /****** BOUNDARY PARTICLES ***********/
                    real const rho = std::max(
                        fvar.rhoMin,
                        std::min(
                            fvar.rhoMax, pn[ii].rho + dt * (gamma_t1 * Rrho[ii] + gamma_t2 * pn[ii].Rrho)
                        )
                    );
                    pnp1[ii].rho = rho;
                    pnp1[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                    pnp1[ii].Rrho = Rrho[ii];
                }
                break;
            }
            case ghost:
            {
#pragma omp for schedule(static) nowait
                for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
                { /****** BOUNDARY PARTICLES ***********/
                    if (near_inlet[ii])
                    { // Don't allow negative pressures
                        real const rho = std::max(
                            fvar.rho0, std::min(
                                           fvar.rhoMax, pn[ii].rho + dt * (gamma_t1 * Rrho[ii] +
                                                                           gamma_t2 * pn[ii].Rrho)
                                       )
                        );
                        pnp1[ii].rho = rho;
                        pnp1[ii].p =
                            pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                        pnp1[ii].Rrho = fmax(0.0, Rrho[ii]);
                    }
                    else
                    {
                        real const rho = std::max(
                            fvar.rhoMin, std::min(
                                             fvar.rhoMax, pn[ii].rho + dt * (gamma_t1 * Rrho[ii] +
                                                                             gamma_t2 * pn[ii].Rrho)
                                         )
                        );
                        pnp1[ii].rho = rho;
                        pnp1[ii].p =
                            pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                        pnp1[ii].Rrho = Rrho[ii];
                    }
                }
                break;
            }
            }
        }

        for (size_t block = svar.nbound; block < svar.nfluid + svar.nbound; block++)
        {
#pragma omp for nowait
            for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
            { /****** FLUID PARTICLES **************/

                /* BUFFER = pipe particle receiving prescribed motion            */
                /* BACK = the latest particle in that column                    */
                /* PIPE = in the pipe, with free motion                         */
                /* FREE = free of the pipe and receives an aero force           */

                /*Check if the particle is clear of the starting area*/
                if (pnp1[ii].b > BUFFER && pnp1[ii].b != OUTLET)
                {
/*For any other particles, intergrate as normal*/
#ifdef ALE
                    pnp1[ii].xi = pn[ii].xi + dt * (pn[ii].v + pnp1[ii].vPert) +
                                  dt2 * (beta_t1 * res[ii] + beta_t2 * pn[ii].acc);
#else
                    pnp1[ii].xi =
                        pn[ii].xi + dt * pn[ii].v + dt2 * (beta_t2 * pn[ii].acc + beta_t1 * res[ii]);
#endif

                    pnp1[ii].v = pn[ii].v + dt * (gamma_t1 * res[ii] + gamma_t2 * pn[ii].acc);
                    pnp1[ii].acc = res[ii];
                    pnp1[ii].Af = Af[ii];
                    pnp1[ii].Rrho = Rrho[ii];

                    real const rho = std::max(
                        fvar.rhoMin,
                        std::min(
                            fvar.rhoMax, pn[ii].rho + dt * (gamma_t1 * Rrho[ii] + gamma_t2 * pn[ii].Rrho)
                        )
                    );
                    pnp1[ii].rho = rho;
                    pnp1[ii].p =
                        pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                }
                else if (pnp1[ii].b == OUTLET)
                { /* For the outlet zone, just perform euler integration of last info */
                    pnp1[ii].xi = pn[ii].xi + dt * pnp1[ii].v;
                }

            } /*End fluid particles*/

            /* Do the buffer particles */
            if (limits[block].block_type == inletZone)
            {
                switch (limits[block].fixed_vel_or_dynamic)
                {
                case 1:
                { // Dynamic inlet
                    StateVecD unorm = limits[block].insert_norm.normalized();
#pragma omp for schedule(static) nowait
                    for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
                    {
                        size_t const& backID = limits[block].back[ii];
                        StateVecD const& xi = pnp1[backID].xi;

                        for (size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
                        { /* Define buffer off the back particle */
                            size_t const& buffID = limits[block].buffer[ii][jj];
                            // Set position as related to the previous particle.
                            pnp1[buffID].xi = xi - svar.dx * (jj + 1.0) * unorm;

                            // How to set density and pressure though?
                            pnp1[buffID].v = pnp1[backID].v;
                            pnp1[buffID].rho = pnp1[backID].rho;
                            pnp1[buffID].p = pnp1[backID].p;
                        }
                    }
                    break;
                }
                case 0:
                { // Fixed velocity
#pragma omp for schedule(static) nowait
                    for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
                    {
                        // size_t const& backID = svar.back[ii];
                        for (size_t jj = 0; jj < limits[block].buffer[ii].size(); ++jj)
                        { /* Define buffer off the back particle */

                            size_t const& buffID = limits[block].buffer[ii][jj];

                            pnp1[buffID].Rrho = Rrho[buffID];

                            real const rho = std::max(
                                fvar.rhoMin,
                                std::min(
                                    fvar.rhoMax, pn[buffID].rho + dt * (gamma_t1 * Rrho[buffID] +
                                                                        gamma_t2 * pn[buffID].Rrho)
                                )
                            );
                            pnp1[buffID].rho = rho;
                            pnp1[buffID].p =
                                pressure_equation(rho, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
                            pnp1[buffID].xi = pn[buffID].xi + dt * pn[buffID].v;
                        }
                    }
                    break;
                }
                }
            }
        } // End block count

    } /*End pragma omp parallel*/
}

real Newmark_Beta::Newmark_Beta(
    Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar,
    size_t const& start, size_t& end, real const& npd, MESH const& cells, LIMITS const& limits,
    OUTL& outlist, real& logbase, uint& iteration, vector<StateVecD>& xih, SPHState const& pn,
    SPHState& pnp1
)
{
    real rms_error = 0.0;
    while (rms_error > svar.minRes)
    {
/*Previous state for error calc*/
#pragma omp parallel for shared(pnp1)
        for (size_t ii = start; ii < end; ++ii)
            xih[ii - start] = pnp1[ii].xi;

        Do_NB_Iter(CELL_TREE, svar, fvar, avar, start, end, npd, cells, limits, outlist, pn, pnp1);

        int errstate = Check_Error(
            SPH_TREE, svar, fvar, start, end, rms_error, logbase, outlist, xih, pn, pnp1, iteration
        );

        if (errstate == 0) /*Continue as normal*/
            iteration++;
        else if (errstate == 1) /*Sub iterations exceeded*/
            break;

    } /*End of subits*/

    return rms_error;
}
