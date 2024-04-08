/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "line.h"
#include <Eigen/Geometry>

void LineShape::check_input(shape_block& block, real& globalspacing, int& fault)
{
    if (!check_vector(block.start))
    {
        std::cout << "ERROR: Block \"" << block.name
                  << "\" starting position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }

    if (!check_vector(block.end))
    {
        std::cout << "ERROR: Block \"" << block.name
                  << "\" ending position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }

#if SIMDIM == 3
    if (!check_vector(block.right))
    {
        std::cout << "ERROR: Block \"" << block.name
                  << "\" ending position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }
#endif

    if (block.dx < 0)
    { /* Have globalspacing override counts */
        if (block.ni < 0
#if SIMDIM == 3
            || block.nj < 0
#endif
        )
        {
            std::cout << "WARNING: Neither block globalspacing or block counters have been defined."
                      << std::endl;
            std::cout << "         Using global globalspacing..." << std::endl;
            block.dx = globalspacing;
        }
    }

#if SIMDIM == 2
    if (block.ni > 0)
    {
        if (block.dx < 0)
        { /* Find the correct globalspacing */
            block.dx = (block.end - block.start).norm() / real(block.ni);
        }
    }
#else
    if (block.ni > 0)
    {
        if (block.nj > 0)
        {
            block.dx = std::min(
                (block.right - block.start).norm() / real(block.ni),
                (block.end - block.right).norm() / real(block.nj)
            );
        }
    }
#endif

    if (block.thickness < 0)
    {
        if (block.nk < 0)
        {
            std::cout << "ERROR: Block \"" << block.name
                      << "\" line thickness has not been correctly defined. Stopping." << std::endl;
            fault = 1;
            return;
        }
    }
    else if (block.nk < 0)
    { /* Use a particle count over the real thickness value */
        if (block.hcpl == 1)
            block.nk = static_cast<int>(ceil(block.thickness / block.dx / sqrt(3.0) * 2.0));
        else
            block.nk = static_cast<int>(ceil(block.thickness / block.dx));
    }

    block.nk = block.nk > 1 ? block.nk : 1;

    /* Guess how many points will be in the line */
    if (block.ni > 0)
    { /* Use the given count for the estimate */
#if SIMDIM == 2
        block.npts = block.ni * block.nk;
#else
        if (block.nj < 1)
        {
            std::cout << "WARNING: ni defined, but nj not." << std::endl;
            size_t nj = static_cast<size_t>(ceil((block.end - block.right).norm() / globalspacing));
            nj = nj > 1 ? nj : 1;
            block.nj = nj;
        }
        block.npts = block.ni * block.nj * block.nk;
#endif
    }
    else
    { /* Need to find it */
/* Need to estimate the number of points on the line */
#if SIMDIM == 2
        size_t ni = static_cast<size_t>(ceil((block.end - block.start).norm() / globalspacing));
        block.ni = ni > 1 ? ni : 1;
        block.npts = block.ni * block.nk;
#else
        size_t ni = static_cast<size_t>(ceil((block.right - block.start).norm() / globalspacing));
        size_t nj = 0;
        if (block.hcpl == 1)
            nj =
                static_cast<int>(ceil((block.end - block.right).norm() / globalspacing / sqrt(3.0) * 2.0)
                );
        else
            nj = static_cast<size_t>(ceil((block.end - block.right).norm() / globalspacing));

        block.ni = ni > 1 ? ni : 1;
        block.nj = nj > 1 ? nj : 1;
        block.npts = block.ni * block.nj * block.nk;
#endif
    }
}

#if SIMDIM == 2
// Create line with n thick particles or a given thickness.
std::vector<StateVecD> LineShape::generate_points(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);
    StateVecD delta = (block.end - block.start).normalized() * globalspacing;
    StateVecD norm(delta[1], -delta[0]);

    for (int ii = 0; ii < block.ni; ++ii)
    {
        for (int jj = 0; jj < block.nk; ++jj)
        {
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = delta * real(ii + 0.5 * (jj % 2)) + 0.5 * norm * sqrt(3.0) * real(jj);
            else
                newPoint = delta * real(ii) + norm * real(jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.push_back(newPoint);
        }
    }
    return points;
}
#else
// Three dimensional plane with n thick particles or a given thickness.
std::vector<StateVecD> LineShape::generate_points(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * MEPSILON * globalspacing);
    /* Find integer counts of sizes alone each dimension */
    StateVecD deltai = (block.right - block.start).normalized() * globalspacing;
    StateVecD deltaj = (block.end - block.right).normalized() * globalspacing;

    StateVecD norm = (deltaj.cross(deltai)).normalized() * globalspacing;

    for (int jj = 0; jj < block.nj; ++jj)
    {
        for (int ii = 0; ii < block.ni; ++ii)
        {
            for (int kk = 0; kk < block.nk; ++kk)
            {
                StateVecD newPoint;
                if (block.hcpl == 1)
                    newPoint = 0.5 * (deltai * real(2 * ii + ((jj + kk) % 2)) +
                                      deltaj * (sqrt(3.0) * (real(jj) + real(kk % 2) / 3)) +
                                      norm * 2 * sqrt(6.0) / 3.0 * real(kk));
                else
                    newPoint = deltai * real(ii) + deltaj * real(jj) + norm * real(kk);

                newPoint += StateVecD::Constant(unif(re));
                newPoint += block.start;
                points.push_back(newPoint);
            }
        }
    }
    return points;
}
#endif
