/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "square.h"

void SquareShape::check_input(SIM const& svar, int& fault)
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

    if (dx < 0)
    {
        if (ni < 1 || nj < 1
#if SIMDIM == 3
            || nk < 1
#endif
        )
        {
            std::cout << "Error: Block \"" << name << "\" spacing has not been defined." << std::endl;
            fault = 1;
        }
        else
        {
            /* Find the correct dx */
            StateVecD dist = end - start;
            real di = dist[0] / real(ni);
            real dj = dist[1] / real(nj);

/* Take the smaller of the two */
#if SIMDIM == 2
            dx = std::min(di, dj);
#else
            real dk = dist[1] / real(nk);
            dx = std::min(di, std::min(dj, dk));
#endif
        }
    }

    /* Estimate how many points on the Block */
#if SIMDIM == 3
#endif
    if (particle_order == 1)
    {
#if SIMDIM == 2
        ni = static_cast<int>(ceil((end[0] - start[0]) / dx));
        nj = static_cast<int>(ceil((end[1] - start[1]) / dx / sqrt(3.0) * 2.0));
#else
        ni = static_cast<int>(ceil((end[0] - start[0]) / dx));
        nj = static_cast<int>(ceil((end[1] - start[1]) / dx / sqrt(3.0) * 2.0));
        nk = static_cast<int>(ceil((end[2] - start[2]) / dx / sqrt(6.0) * 3.0));
#endif
    }
    else
    {
        ni = static_cast<int>(ceil((end[0] - start[0]) / dx));
        nj = static_cast<int>(ceil((end[1] - start[1]) / dx));
#if SIMDIM == 3
        nk = static_cast<int>(ceil((end[2] - start[2]) / dx));
#endif
    }
    ni = ni > 1 ? ni : 1;
    nj = nj > 1 ? nj : 1;
#if SIMDIM == 3
    nk = nk > 1 ? nk : 1;
    npts = ni * nk * nj;
#else
    npts = ni * nj;
#endif

    ShapeBlock::check_input_post();
}

void SquareShape::generate_points()
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * dx - 2.0 * tol * dx);
#if SIMDIM == 3
    for (int k = 0; k < nk; ++k)
    {
#endif
        for (int j = 0; j < nj; ++j)
        {
            for (int i = 0; i < ni; ++i)
            {
                StateVecD newPoint;
#if SIMDIM == 2
                if (particle_order == 1)
                    newPoint =
                        0.5 *
                        StateVecD(static_cast<real>(2 * i + (j % 2)), sqrt(3.0) * static_cast<real>(j));
                else
                    newPoint = StateVecD(static_cast<real>(i), static_cast<real>(j));
#else
            if (particle_order == 1)
                newPoint =
                    0.5 * StateVecD(
                              real(2 * i + ((j + k) % 2)), sqrt(3.0) * (real(j) + real((k % 2)) / 3.0),
                              2.0 / 3.0 * sqrt(6.0) * real(k)
                          );
            else
                newPoint = StateVecD(real(i), real(j), real(k));
#endif

                newPoint *= dx;
                newPoint += StateVecD::Constant(unif(re));
                newPoint += start;
                points.emplace_back(newPoint);

            } /* end k count */
        } /* end j count */
#if SIMDIM == 3
    } /* end i count */
#endif

    coords = points;
}