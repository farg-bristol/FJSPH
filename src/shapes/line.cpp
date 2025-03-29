/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "line.h"

#include "../Third_Party/Eigen/Geometry"

void LineShape::check_input(SIM const& svar, int& fault)
{
    // Do common input checks.
    ShapeBlock::check_input(svar, fault);

    if (!check_vector(start))
    {
        std::cout << "ERROR: Block \"" << name
                  << "\" starting position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }

    if (!check_vector(end))
    {
        std::cout << "ERROR: Block \"" << name
                  << "\" ending position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }

#if SIMDIM == 3
    if (!check_vector(right))
    {
        std::cout << "ERROR: Block \"" << name
                  << "\" ending position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }
#endif

    if (dx < 0)
    {
        if (ni < 1
#if SIMDIM == 3
            || nj < 1
#endif
        )
        {
            std::cout << "Error: Block \"" << name << "\" spacing has not been defined." << std::endl;
            fault = 1;
        }
        else
        {
#if SIMDIM == 2
            dx = (end - start).norm() / real(ni);
#else
            dx = std::min((right - start).norm() / real(ni), (end - right).norm() / real(nj));
#endif
        }
    }

    if (thickness < 0)
    {
        if (nk < 0)
        {
            std::cout << "ERROR: Block \"" << name
                      << "\" line thickness has not been correctly defined. Stopping." << std::endl;
            fault = 1;
        }
    }
    else if (nk < 0)
    { /* Use a particle count over the real thickness value */
        if (particle_order == 1)
            nk = static_cast<int>(ceil(thickness / dx / sqrt(3.0) * 2.0));
        else
            nk = static_cast<int>(ceil(thickness / dx));
    }

    nk = nk > 1 ? nk : 1;

    /* Guess how many points will be in the line */
    if (ni > 0)
    { /* Use the given count for the estimate */
#if SIMDIM == 2
        npts = ni * nk;
#else
        if (nj < 1)
        {
            std::cout << "WARNING: ni defined, but nj not." << std::endl;
            nj = static_cast<size_t>(ceil((end - right).norm() / dx));
            nj = nj > 1 ? nj : 1;
        }
        npts = ni * nj * nk;
#endif
    }
    else
    { /* Need to find it */
/* Need to estimate the number of points on the line */
#if SIMDIM == 2
        ni = static_cast<size_t>(ceil((end - start).norm() / dx));
        ni = ni > 1 ? ni : 1;
        npts = ni * nk;
#else
        ni = static_cast<size_t>(ceil((right - start).norm() / dx));
        if (particle_order == 1)
            nj = static_cast<int>(ceil((end - right).norm() / dx / sqrt(3.0) * 2.0));
        else
            nj = static_cast<size_t>(ceil((end - right).norm() / dx));

        ni = ni > 1 ? ni : 1;
        nj = nj > 1 ? nj : 1;
        npts = ni * nj * nk;
#endif
    }

    ShapeBlock::check_input_post();
}

void LineShape::generate_points()
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;
#if SIMDIM == 2
    // Create line with n thick particles or a given thickness.

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * dx - 2.0 * tol * dx);
    StateVecD delta = (end - start).normalized() * dx;
    StateVecD norm(delta[1], -delta[0]);

    for (int ii = 0; ii < ni; ++ii)
    {
        for (int jj = 0; jj < nk; ++jj)
        {
            StateVecD newPoint;
            if (particle_order == 1)
                newPoint = delta * real(ii + 0.5 * (jj % 2)) + 0.5 * norm * sqrt(3.0) * real(jj);
            else
                newPoint = delta * real(ii) + norm * real(jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += start;
            points.push_back(newPoint);
        }
    }
#else
    // Three dimensional plane with n thick particles or a given thickness.

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * dx - 2.0 * MEPSILON * dx);
    /* Find integer counts of sizes alone each dimension */
    StateVecD deltai = (right - start).normalized() * dx;
    StateVecD deltaj = (end - right).normalized() * dx;

    StateVecD norm = (deltaj.cross(deltai)).normalized() * dx;

    for (int jj = 0; jj < nj; ++jj)
    {
        for (int ii = 0; ii < ni; ++ii)
        {
            for (int kk = 0; kk < nk; ++kk)
            {
                StateVecD newPoint;
                if (particle_order == 1)
                    newPoint = 0.5 * (deltai * real(2 * ii + ((jj + kk) % 2)) +
                                      deltaj * (sqrt(3.0) * (real(jj) + real(kk % 2) / 3)) +
                                      norm * 2 * sqrt(6.0) / 3.0 * real(kk));
                else
                    newPoint = deltai * real(ii) + deltaj * real(jj) + norm * real(kk);

                newPoint += StateVecD::Constant(unif(re));
                newPoint += start;
                points.push_back(newPoint);
            }
        }
    }
#endif
    coords = points;
}
