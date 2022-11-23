
#include "circle.h"

void check_cylinder_input(shape_block& bound, real& globalspacing)
{
    int has_config = 0;
    int arc_defined = 0;

    #if SIMDIM == 2
    // These checks are only valid for 2D
    if(check_vector(bound.centre))
    {   /* Default to these values first */
        if(bound.arc_start >= 0 && bound.arc_end >= 0 && bound.radius > 0)
        {
            has_config = 1;
        }
        else if(bound.arc_start >= 0 && bound.ni > 0 && bound.radius > 0)
        {
            has_config = 1;
        }
        else if(bound.arc_start >= 0 && bound.arclength != default_val && bound.radius > 0)
        {
            has_config = 1;
            arc_defined = 1;
        }
    }
    #endif
    
    

    // #if SIMDIM == 3
    // if(!has_config)
    // {
    //     // check arch geometry first. Want to use an arclength first
    //     // Midpoint is required to define the plane appropriately, but not for the circle
    //     if(check_vector(bound.centre) && check_vector(bound.start) && check_vector(bound.mid)
    //         && bound.arclength != default_val)
    //     {
    //         has_config = 1;
    //         arc_defined = 1;
    //     }
    // }
    // #else
    // if(!has_config)
    // {
    //     if(check_vector(bound.start) && check_vector(bound.end) && check_vector(bound.mid))
    //     {   // No arclength, but just the end acquired
    //         has_config = 1;
    //         #if SIMDIM == 2
    //         get_arclength_midpoint(bound.start,bound.end,bound.mid,
    //                     bound.centre,bound.radius,bound.arc_start,bound.arc_end);
    //         #else
    //         get_arclength_midpoint(bound.start,bound.end,bound.mid,bound.centre,
    //                     bound.right,bound.radius,bound.arc_start,bound.arc_end);
    //         #endif
    //     }
    // }
    // #endif

    if(!has_config)
    {   // Check if the arclength is defined
        if(check_vector(bound.centre) && check_vector(bound.start) && bound.arclength != default_val)
        {  
            #if SIMDIM == 3
            if(check_vector(bound.right))
            {   // Need to check if normal is defined to know the plane
                has_config = 1;
                arc_defined = 1;
            }
            #else
            has_config = 1;
            arc_defined = 1;
            #endif
        }
    }

    if(!has_config)
    {
        std::cout << "ERROR: Block \"" << bound.name << 
            "\" arc geometry has not been sufficiently defined. Stopping." << std::endl;
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
    

    if(bound.thickness < 0)
    {
        if(bound.nk < 0)
        {
            std::cout << "ERROR: Block \"" << bound.name << "\" arc thickness has not been correctly defined. Stopping." << std::endl;
            exit(1);
        }
    }
    else if (bound.nk < 0)
    {   /* Use a particle count over the real thickness value */
        if(bound.hcpl == 1)
            bound.nk = static_cast<int>(ceil(bound.thickness / globalspacing / sqrt(3.0) * 2.0));
        else
            bound.nk = static_cast<int>(ceil(bound.thickness / globalspacing));
    }

    // Check that start coordinate exists
    // #if SIMDIM == 2
    // // Need to find the start vector first if it's not defined
    // if(!check_vector(bound.start))
    // {
    //     bound.start(bound.radius*cos(bound.arc_start), bound.radius*sin(bound.arc_start));
    //     bound.start += bound.centre;
    // }
    // #endif

    real dtheta = globalspacing/bound.radius;
    if(arc_defined)
    {

        #if SIMDIM == 3
        // #else
            
            // get_arc_end(bound.start,bound.centre,bound.arclength,bound.end,bound.radius,bound.arc_start);
        #endif
        // bound.arclength *= M_PI/180.0;
  
    }
    else
    {
        /* Estimate how many points on the Block */     
        if(bound.arc_end < 0 && bound.ni > 0)
        {   /* Can define arc using only starting angle and number of particles, to generate an end */
            bound.arc_end = 180.0/M_PI * (bound.arc_start*M_PI/180 + real(bound.ni) * dtheta);
        }
        bound.arclength = (bound.arc_end-bound.arc_start);
    }

    int ni = static_cast<int>(ceil((fabs(bound.arclength)*M_PI/180)/dtheta));
    ni = ni > 1 ? ni : 1;
    bound.ni = ni;
    // Need to add straights in to the count
    
    int smax = bound.sstraight > 0 ? ceil(bound.sstraight/globalspacing) : 0;
    int emax = bound.estraight > 0 ? ceil(bound.estraight/globalspacing) : 0;

    bound.npts = (ni + smax + emax) * bound.nk;
    #if SIMDIM == 3
    if(bound.length < 0)
    {
        if(bound.nj < 0)
        {
            std::cout << "ERROR: Block \"" << bound.name << "\" arch length has not been correctly defined. Stopping." << std::endl;
            exit(1);
        }
    }
    else if (bound.nj < 0)
    {
        if(bound.hcpl == 1)
            bound.nj = static_cast<int>(ceil((bound.length) / globalspacing / sqrt(6.0) * 3.0));
        else
            bound.nj = static_cast<int>(ceil(bound.length / globalspacing));
        bound.nj = bound.nj > 1 ? bound.nj : 1;
    }
    bound.npts *= bound.nj;
    #endif
}

#if SIMDIM == 2
std::vector<StateVecD> create_cylinder(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 


    StateVecD norm = block.normal.normalized();
    StateVecD left(norm[1],-norm[0]);
    real r = block.radius;
    for(int ii = 0; ii < block.ni; ii++)
    {   // Depth
        for(int kk = 0; kk < block.nk; kk++)
        {   //Thickness
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = globalspacing * (norm * real(ii + 0.5*(kk % 2)) + 0.5 * left * sqrt(3.0) * real(kk))
                    + r * left;
            else
                newPoint = globalspacing * (norm * real(ii) + left * real(kk)) + r * left;

            newPoint += StateVecD(unif(re),unif(re));
            newPoint += block.centre;
            points.push_back(newPoint);
        }
    }
    // Other wall, but minus radius
    for(int ii = 0; ii < block.ni; ii++)
    {   // Depth
        for(int kk = 0; kk < block.nk; kk++)
        {   //Thickness
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = globalspacing * (norm * real(ii + 0.5*(kk % 2)) + 0.5 * left * sqrt(3.0) * real(kk))
                    + r * left;
            else
                newPoint = globalspacing * (norm * real(ii) + left * real(kk)) + r * left;

            newPoint += StateVecD(unif(re),unif(re));
            newPoint += block.centre;
            points.push_back(newPoint);
        }
    }
    return points;   
}
#endif

#if SIMDIM == 3
std::vector<StateVecD> create_cylinder(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    StateVecD norm = block.normal.normalized();
    StateVecD left = (block.right - block.centre).normalized();
    StateVecD right = left.cross(norm).normalized();

    real dtheta = globalspacing / block.radius;
    size_t nj = ceil(2*M_PI/dtheta);
        // Create local coordinate system (u,v,w)
    dtheta = 2*M_PI/real(nj); // Redefine dtheta so that it makes a perfect circle

    for(int ii = 0; ii < block.ni; ii++)
    {   // Depth
        for(int kk = 0; kk < block.nk; kk++)
        {   //Thickness
            real r = block.radius;
            real l = real(kk)*globalspacing;
            real doffset = 0;
            if(block.hcpl == 1)
            {
                r += 1.0 / 3.0 * sqrt(6.0) * real(kk) * globalspacing;
                l = 0.5*sqrt(3)*(real(ii) + real(kk%2)/3.0)*globalspacing;
                doffset = 0.5*dtheta*((kk+ii)%2);
            }
            else
                r += real(kk)*globalspacing;

            for(size_t jj = 0; jj < nj; jj++)
            {   // Radius
                real theta = static_cast<real>(jj)*dtheta;
                real l = real(ii)*globalspacing;
                if(block.hcpl == 1)
                {   /* Shift by half a globalspacing if jj is even */
                    theta += doffset;
                }

                real a = cos(theta)*r; real b = sin(theta)*r;
                StateVecD newPoint = a * left + b * right + l*norm;

                newPoint += StateVecD(unif(re),unif(re),unif(re));
                newPoint += block.centre;
                points.push_back(newPoint);
            }
        }
    }
    return points;   
}
#endif