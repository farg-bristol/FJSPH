
#include "circle.h"

void check_circle_input(shape_block& bound, real& globalspacing)
{
    
    if(!check_vector(bound.centre))
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" centre position has not been correctly defined. Stopping" << std::endl;
        exit(1);
    }

    if(bound.radius < 0)
    {
        std::cout << "ERROR: Block \"" << bound.name << "\" radius has not been correctly defined. Stopping." << std::endl;
        exit(1);
    }

    if(bound.dx < 0)
    {
        if(bound.ni < 0)
        {
            std::cout << "WARNING: Block globalspacing has not been defined. Using global globalspacing." << std::endl;
            bound.dx = globalspacing;
        }
        else
        {
            /* Use ni as a number along the diameter, so define globalspacing off that */
            real di = (2.0*bound.radius)/real(bound.ni);
            bound.dx = di;
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
    npoints = ni * nj;
    #if SIMDIM == 3
    nk = nk > 1 ? nk : 1;
    bound.nk = nk;
    npoints *= nk;
    #endif
    npoints = npoints > 1 ? npoints : 1;
    bound.npts = npoints;
    
    /* Multiply by the ratio of circle to square area */
    #if SIMDIM == 2
    bound.npts = ceil(real(bound.npts)*M_PI/4.0);
    #else
    bound.npts = ceil(real(bound.npts)*M_PI/6.0);
    #endif
}


std::vector<StateVecD> create_circle(StateVecD const& centre, real const& radius, 
                real const& globalspacing, int const& hcpl )
{
    std::vector<StateVecD> points;

    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);

    /* Start by creating a lattice rectangle, then test if point lies inside the circle radius */
    StateVecD start = centre - StateVecD::Constant(radius);
    StateVecD end = centre + StateVecD::Constant(radius);

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

            if((centre - newPoint).squaredNorm() > radius)
            {
                continue;
            }

            points.push_back(newPoint);

        } /* end k count */
    } /* end j count */
    #if SIMDIM == 3
    }   /* end i count */
    #endif

    return points;
}