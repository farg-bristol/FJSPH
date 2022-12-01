#include "line.h"
#include <Eigen/Geometry>

void check_line_input(shape_block& bound, real& globalspacing)
{
    if(!check_vector(bound.start))
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" starting position has not been correctly defined. Stopping" << std::endl;
        exit(1);
    }

    if(!check_vector(bound.end))
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" ending position has not been correctly defined. Stopping" << std::endl;
        exit(1);
    }

    #if SIMDIM == 3
    if(!check_vector(bound.right))
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" ending position has not been correctly defined. Stopping" << std::endl;
        exit(1);
    }
    #endif

    if(bound.dx < 0)
    {   /* Have globalspacing override counts */
        if(bound.ni < 0 
        #if SIMDIM == 3
        || bound.nj < 0
        #endif
        )
        {
            std::cout << "WARNING: Neither block globalspacing or block counters have been defined." << std::endl;
            std::cout << "         Using global globalspacing..." << std::endl;
            bound.dx = globalspacing;
        }
    }

    #if SIMDIM == 2
    if(bound.ni > 0)
    {
        if(bound.dx < 0)
        {   /* Find the correct globalspacing */
            bound.dx = (bound.end - bound.start).norm() / real(bound.ni);
        }
    }
    #else
    if(bound.ni > 0)
    {
        if(bound.nj > 0)
        {
            bound.dx = std::min((bound.right - bound.start).norm() / real(bound.ni), 
                                (bound.end - bound.right).norm() / real(bound.nj));
        }
    }
    #endif

    if(bound.thickness < 0)
    {
        if(bound.nk < 0)
        {
            std::cout << "ERROR: Block \"" << bound.name << "\" line thickness has not been correctly defined. Stopping." << std::endl;
            exit(1);
        }
    }
    else if (bound.nk < 0)
    {   /* Use a particle count over the real thickness value */
        if(bound.hcpl == 1)
            bound.nk = static_cast<int>(ceil(bound.thickness / bound.dx / sqrt(3.0) * 2.0));
        else
            bound.nk = static_cast<int>(ceil(bound.thickness / bound.dx));
    }

    bound.nk = bound.nk > 1 ? bound.nk : 1;

    /* Guess how many points will be in the line */
    if (bound.ni > 0)
    {   /* Use the given count for the estimate */
        #if SIMDIM == 2 
        bound.npts = bound.ni*bound.nk;
        #else
        if(bound.nj < 1)
        {
            std::cout << "WARNING: ni defined, but nj not." << std::endl;
            size_t nj = static_cast<size_t>(ceil((bound.end - bound.right).norm() / globalspacing));
            nj = nj > 1 ? nj : 1;
            bound.nj = nj;
        }
        bound.npts = bound.ni*bound.nj*bound.nk;
        #endif
    }
    else
    {   /* Need to find it */
        /* Need to estimate the number of points on the line */
        #if SIMDIM == 2
        size_t ni = static_cast<size_t>(ceil((bound.end - bound.start).norm() / globalspacing));
        bound.ni = ni > 1 ? ni : 1;
        bound.npts = bound.ni* bound.nk;
        #else
        size_t ni = static_cast<size_t>(ceil((bound.right - bound.start).norm() / globalspacing));
        size_t nj = 0;
        if (bound.hcpl == 1)
            nj = static_cast<int>(ceil((bound.end - bound.right).norm() / globalspacing / sqrt(3.0) * 2.0));        
        else
            nj = static_cast<size_t>(ceil((bound.end - bound.right).norm() / globalspacing));

        bound.ni = ni > 1 ? ni : 1;
        bound.nj = nj > 1 ? nj : 1;
        bound.npts = bound.ni * bound.nj * bound.nk;
        #endif
    }
}

// Fine line of one particle thick, used in particular for fibres
std::vector<StateVecD> create_fine_line(StateVecD const& start, StateVecD const& end, 
                         real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);
    int ni = static_cast<int>(ceil((end - start).norm() / globalspacing));
    ni = ni > 1 ? ni : 1;
    StateVecD delta = (end - start) / static_cast<real>(ni);

    for(int ii = 0; ii < ni; ++ii)
    {
        StateVecD newPoint = delta * real(ii);
        newPoint += StateVecD::Constant(unif(re));
        newPoint += start;
        points.push_back(newPoint);
    }
    return points;
}

#if SIMDIM == 2
// Create line with n thick particles or a given thickness.
std::vector<StateVecD> create_line(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);
    StateVecD delta = (block.end - block.start).normalized() * globalspacing;
    StateVecD norm(delta[1], -delta[0]);

    for(int ii = 0; ii < block.ni; ++ii)
    {
        for(int jj = 0; jj < block.nk; ++jj)
        {
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = delta * real(ii + 0.5*(jj % 2)) + 0.5 * norm * sqrt(3.0) * real(jj);
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
std::vector<StateVecD> create_plane(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    /* To check if any points end up too close to each other */
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * MEPSILON * globalspacing);
    /* Find integer counts of sizes alone each dimension */
    StateVecD deltai = (block.right - block.start).normalized() * globalspacing;
    StateVecD deltaj = (block.end - block.right).normalized() * globalspacing;

    StateVecD norm = (deltaj.cross(deltai)).normalized() * globalspacing;

    for(int jj = 0; jj < block.nj; ++jj)
    {
        for(int ii = 0; ii < block.ni; ++ii)
        {
            for(int kk = 0; kk < block.nk; ++kk)
            {
                StateVecD newPoint;
                if (block.hcpl == 1)
                    newPoint =  0.5*(deltai * real(2*ii + ((jj+kk)%2)) + 
                                deltaj * (sqrt(3.0) * (real(jj) + real(kk%2)/3)) + 
                                norm * 2 * sqrt(6.0) / 3.0 * real(kk)); 
                else
                    newPoint = deltai * real(ii) + deltaj * real(jj) + norm*real(kk);

                newPoint += StateVecD::Constant(unif(re));
                newPoint += block.start;
                points.push_back(newPoint);
            }
        }
    }
    return points;
}

#endif

