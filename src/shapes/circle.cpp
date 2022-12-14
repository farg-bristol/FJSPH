
#include "circle.h"

void check_circle_input(shape_block& block, real& globalspacing)
{
    
    int fault = 0;
    if(!check_vector(block.centre))
    {
        printf("ERROR: Block \"%s\" centre position has not been correctly defined.\n", block.name.c_str());
        fault = 1;
    }
    

    if(block.radius < 0)
    {
        printf("ERROR: Block \"%s\" radius has not been correctly defined.\n", block.name.c_str());
        fault = 1;
    }

    block.start = block.centre;
    block.start -= StateVecD::Constant(block.radius);
    block.end = block.centre;
    block.end += StateVecD::Constant(block.radius);


    if(block.dx < 0)
    {
        if(block.ni < 0)
        {
            std::cout << "WARNING: Block globalspacing has not been defined. Using global globalspacing." << std::endl;
            block.dx = globalspacing;
        }
        else
        {
            /* Use ni as a number along the diameter, so define globalspacing off that */
            real di = (2.0*block.radius)/real(block.ni);
            block.dx = di;
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
    npoints = ni * nj;
    #if SIMDIM == 3
    nk = nk > 1 ? nk : 1;
    block.nk = nk;
    npoints *= nk;
    #endif
    npoints = npoints > 1 ? npoints : 1;
    block.npts = npoints;
    
    /* Multiply by the ratio of circle to square area */
    #if SIMDIM == 2
    block.npts = ceil(real(block.npts)*M_PI/4.0);
    #else
    block.npts = ceil(real(block.npts)*M_PI/6.0);
    #endif

    if(fault)
    {
        printf("Check of circle or sphere geometry finished with errors. Stopping.\n");
        exit(-1);
    }
}

#if SIMDIM == 2
std::vector<StateVecD> create_circle(shape_block& block, real const& globalspacing)
{
    std::vector<StateVecD> points;

    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re;

    StateVecD const& centre = block.centre;
    real const& radius = block.radius;    
    StateMatD const& rotmat = block.rotmat;
    for(real rad = radius; rad > 0.99*globalspacing; rad -= globalspacing)
	{
		// % Find spacing to have a well defined surface.
		real dtheta = atan(globalspacing / rad);
		int ncirc = floor(abs(2.0 * M_PI / dtheta));
		dtheta = 2.0 * M_PI / real(ncirc);
		
		for(real theta = 0.0; theta < 2*M_PI-0.5*dtheta; theta += dtheta)
		{	/* Create a ring of points */
			real x = rad * sin(theta);
			real y = rad * cos(theta);

			StateVecD newPoint(x,y);
			newPoint += StateVecD(unif(re),unif(re));
            newPoint  = rotmat * newPoint;
            newPoint += centre;
            points.emplace_back(newPoint);
		}
	}
	
	/* Create centre point */
	StateVecD newPoint(0.0,0.0);
    newPoint += StateVecD(unif(re),unif(re));
    newPoint  = rotmat * newPoint;
    newPoint += centre;
    points.emplace_back(newPoint);
    
    return points;
}
#endif

#if SIMDIM == 3
std::vector<StateVecD> create_circle(shape_block& block, real const& globalspacing)
{
    std::vector<StateVecD> points;

    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 
    // real searchDist = pow2(0.5 * globalspacing - 2.0 * tol * globalspacing);

    StateVecD const& centre = block.centre;
    real const& radius = block.radius * block.radius;
    int const& hcpl = block.hcpl;

    /* Start by creating a lattice rectangle, then test if point lies inside the circle radius */
    // StateVecD start = centre - StateVecD::Constant(radius);
    // StateVecD end = centre + StateVecD::Constant(radius);
    StateVecD start = block.start;
    StateVecD end = block.end;

    /* Find integer counts of sizes alone each dimension */
    int ni;
    int nj;
    int nk;
    if (hcpl == 1)
    {
        ni = static_cast<int>(ceil((end[0] - start[0]) / globalspacing));
        nj = static_cast<int>(ceil((end[1] - start[1]) / globalspacing / sqrt(3.0) * 2.0));
        nk = static_cast<int>(ceil((end[2] - start[2]) / globalspacing / sqrt(6.0) * 3.0));
    }
    else
    {
        ni = static_cast<int>(ceil((end[0] - start[0]) / globalspacing));
        nj = static_cast<int>(ceil((end[1] - start[1]) / globalspacing));
        nk = static_cast<int>(ceil((end[2] - start[2]) / globalspacing));        
    }
    ni = ni > 1 ? ni : 1;
    nj = nj > 1 ? nj : 1;
    nk = nk > 1 ? nk : 1;
    for (int k = 0; k < nk; ++k)
    {
        for (int j = 0; j < nj; ++j)
        {
            for (int i = 0; i < ni; ++i)
            {
                StateVecD newPoint;
                if (hcpl == 1)
                    newPoint = 0.5 * StateVecD(real(2 * i + ((j + k) % 2)), sqrt(3.0) * (real(j) + real((k % 2)) / 3.0), 2.0 / 3.0 * sqrt(6.0) * real(k));
                else
                    newPoint = StateVecD(real(i), real(j), real(k));

                newPoint *= globalspacing;
                newPoint  = block.rotmat * newPoint;
                newPoint += StateVecD(unif(re),unif(re),unif(re));
                newPoint += start;

                if((newPoint - centre).squaredNorm() > radius)
                {
                    continue;
                }

                points.push_back(newPoint);

            } /* end k count */
        } /* end j count */
    }   /* end i count */

    return points;
}
#endif