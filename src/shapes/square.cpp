/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "square.h"

void SquareShape::check_input(shape_block& block, real& globalspacing, int& fault)
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

    if (block.dx < 0)
    { /* Have globalspacing override counts */
        if (block.ni < 0 || block.nj < 0)
        {
            std::cout << "WARNING: Neither block globalspacing or block counters have been defined."
                      << std::endl;
            std::cout << "         Using global globalspacing..." << std::endl;
            block.dx = globalspacing;
        }
    }

    if (block.ni > 0)
    {
        if (block.nj > 0)
        {
#if SIMDIM == 3
            if (block.nk > 0)
            {
#endif
                if (block.dx < 0)
                {
                    /* Find the correct globalspacing */
                    StateVecD dist = block.end - block.start;
                    real di = dist[0] / real(block.ni);
                    real dj = dist[1] / real(block.nj);

/* Take the smaller of the two */
#if SIMDIM == 2
                    block.dx = std::min(di, dj);
#else
                real dk = dist[1] / real(block.nk);
                block.dx = std::min(di, std::min(dj, dk));
#endif
                }
#if SIMDIM == 3
            }
#endif
        }
    }

    /* Estimate how many points on the Block */
    int npoints;
    size_t ni, nj;
#if SIMDIM == 3
    size_t nk;
#endif
    if (block.hcpl == 1)
    {
#if SIMDIM == 2
        ni = static_cast<int>(ceil((block.end[0] - block.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((block.end[1] - block.start[1]) / globalspacing / sqrt(3.0) * 2.0));
#else
        ni = static_cast<int>(ceil((block.end[0] - block.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((block.end[1] - block.start[1]) / globalspacing / sqrt(3.0) * 2.0));
        nk = static_cast<int>(ceil((block.end[2] - block.start[2]) / globalspacing / sqrt(6.0) * 3.0));
#endif
    }
    else
    {
        ni = static_cast<int>(ceil((block.end[0] - block.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((block.end[1] - block.start[1]) / globalspacing));
#if SIMDIM == 3
        nk = static_cast<int>(ceil((block.end[2] - block.start[2]) / globalspacing));
#endif
    }
    ni = ni > 1 ? ni : 1;
    nj = nj > 1 ? nj : 1;
    block.ni = ni;
    block.nj = nj;
#if SIMDIM == 3
    nk = nk > 1 ? nk : 1;
    block.nk = nk;
    npoints = ni * nk * nj;
#else
    npoints = ni * nj;
#endif
    npoints = npoints > 1 ? npoints : 1;
    block.npts = npoints;
}

std::vector<StateVecD> SquareShape::generate_points(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);

    /* Find integer counts of sizes alone each dimension */
    int ni;
    int nj;
#if SIMDIM == 3
    int nk;
#endif

    if (block.hcpl == 1)
    {
#if SIMDIM == 2
        ni = static_cast<int>(ceil((block.end[0] - block.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((block.end[1] - block.start[1]) / globalspacing / sqrt(3.0) * 2.0));
#else
        ni = static_cast<int>(ceil((block.end[0] - block.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((block.end[1] - block.start[1]) / globalspacing / sqrt(3.0) * 2.0));
        nk = static_cast<int>(ceil((block.end[2] - block.start[2]) / globalspacing / sqrt(6.0) * 3.0));
#endif
    }
    else
    {
        ni = static_cast<int>(ceil((block.end[0] - block.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((block.end[1] - block.start[1]) / globalspacing));
#if SIMDIM == 3
        nk = static_cast<int>(ceil((block.end[2] - block.start[2]) / globalspacing));
#endif
    }
    ni = ni > 1 ? ni : 1;
    nj = nj > 1 ? nj : 1;
#if SIMDIM == 3
    nk = nk > 1 ? nk : 1;
    for (int k = 0; k < nk; ++k)
    {
#endif
        for (int j = 0; j < nj; ++j)
        {
            for (int i = 0; i < ni; ++i)
            {
                StateVecD newPoint;
#if SIMDIM == 2
                if (block.hcpl == 1)
                    newPoint =
                        0.5 *
                        StateVecD(static_cast<real>(2 * i + (j % 2)), sqrt(3.0) * static_cast<real>(j));
                else
                    newPoint = StateVecD(static_cast<real>(i), static_cast<real>(j));
#else
            if (block.hcpl == 1)
                newPoint =
                    0.5 * StateVecD(
                              real(2 * i + ((j + k) % 2)), sqrt(3.0) * (real(j) + real((k % 2)) / 3.0),
                              2.0 / 3.0 * sqrt(6.0) * real(k)
                          );
            else
                newPoint = StateVecD(real(i), real(j), real(k));
#endif

                newPoint *= globalspacing;
                newPoint += StateVecD::Constant(unif(re));
                newPoint += block.start;
                points.emplace_back(newPoint);

            } /* end k count */
        } /* end j count */
#if SIMDIM == 3
    } /* end i count */
#endif

    return points;
}