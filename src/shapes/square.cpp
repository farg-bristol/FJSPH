
#include "square.h"

void check_square_input(shape_block& bound, real& globalspacing, int& fault)
{
    
    if(!check_vector(bound.start))
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" starting position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }

    if(!check_vector(bound.end))
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" ending position has not been correctly defined. Stopping" << std::endl;
        fault = 1;
    }

    if(bound.dx < 0)
    {   /* Have globalspacing override counts */
        if(bound.ni < 0 || bound.nj < 0)
        {
            std::cout << "WARNING: Neither block globalspacing or block counters have been defined." << std::endl;
            std::cout << "         Using global globalspacing..." << std::endl;
            bound.dx = globalspacing;
        }
    }

    if(bound.ni > 0)
    {
        if(bound.nj > 0)
        {
            #if SIMDIM == 3
            if(bound.nk > 0)
            {
            #endif
            if(bound.dx < 0)
            {
                /* Find the correct globalspacing */
                StateVecD dist = bound.end-bound.start;
                real di = dist[0]/real(bound.ni);
                real dj = dist[1]/real(bound.nj);

                /* Take the smaller of the two */
                #if SIMDIM == 2
                bound.dx = std::min(di,dj);
                #else
                real dk = dist[1]/real(bound.nk);
                bound.dx = std::min(di,std::min(dj,dk));
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
    if (bound.hcpl == 1)
    {
        #if SIMDIM == 2
        ni = static_cast<int>(ceil((bound.end[0] - bound.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((bound.end[1] - bound.start[1]) / globalspacing / sqrt(3.0) * 2.0));
        #else
        ni = static_cast<int>(ceil((bound.end[0] - bound.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((bound.end[1] - bound.start[1]) / globalspacing / sqrt(3.0) * 2.0));
        nk = static_cast<int>(ceil((bound.end[2] - bound.start[2]) / globalspacing / sqrt(6.0) * 3.0));
        #endif
    }
    else
    {
        ni = static_cast<int>(ceil((bound.end[0] - bound.start[0]) / globalspacing));
        nj = static_cast<int>(ceil((bound.end[1] - bound.start[1]) / globalspacing));
        #if SIMDIM == 3
        nk = static_cast<int>(ceil((bound.end[2] - bound.start[2]) / globalspacing));
        #endif
        
    }
    ni = ni > 1 ? ni : 1;
    nj = nj > 1 ? nj : 1;
    bound.ni = ni;
    bound.nj = nj;
    #if SIMDIM == 3
    nk = nk > 1 ? nk : 1;
    bound.nk = nk;
    npoints = ni * nk * nj;
    #else
    npoints = ni * nj;
    #endif
    npoints = npoints > 1 ? npoints : 1;
    bound.npts = npoints;
}


std::vector<StateVecD> create_square(StateVecD const& start, StateVecD const& end, 
                        real const& globalspacing, int const& hcpl)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);

    /* Find integer counts of sizes alone each dimension */
    int ni;
    int nj;
    #if SIMDIM == 3
    int nk;
    #endif

    if (hcpl == 1)
    {
        #if SIMDIM == 2
        ni = static_cast<int>(ceil((end[0] - start[0]) / globalspacing));
        nj = static_cast<int>(ceil((end[1] - start[1]) / globalspacing / sqrt(3.0) * 2.0));
        #else
        ni = static_cast<int>(ceil((end[0] - start[0]) / globalspacing));
        nj = static_cast<int>(ceil((end[1] - start[1]) / globalspacing / sqrt(3.0) * 2.0));
        nk = static_cast<int>(ceil((end[2] - start[2]) / globalspacing / sqrt(6.0) * 3.0));
        #endif
    }
    else
    {
        ni = static_cast<int>(ceil((end[0] - start[0]) / globalspacing));
        nj = static_cast<int>(ceil((end[1] - start[1]) / globalspacing));
        #if SIMDIM == 3
        nk = static_cast<int>(ceil((end[2] - start[2]) / globalspacing));
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
            if (hcpl == 1)
                newPoint = 0.5 * StateVecD(static_cast<real>(2 * i + (j % 2)), sqrt(3.0) * static_cast<real>(j));
            else
                newPoint = StateVecD(static_cast<real>(i), static_cast<real>(j));
            #else
            if (hcpl == 1)
                newPoint = 0.5 * StateVecD(real(2 * i + ((j + k) % 2)), sqrt(3.0) * (real(j) + real((k % 2)) / 3.0), 2.0 / 3.0 * sqrt(6.0) * real(k));
            else
                newPoint = StateVecD(real(i), real(j), real(k));
            #endif

            newPoint *= globalspacing;
            newPoint += StateVecD::Constant(unif(re));
            newPoint += start;
            points.emplace_back(newPoint);

        } /* end k count */
    } /* end j count */
    #if SIMDIM == 3
    }   /* end i count */
    #endif

    return points;
}
