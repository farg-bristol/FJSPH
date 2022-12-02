
#include "cylinder.h"
#include <Eigen/Geometry>

void check_cylinder_input(shape_block& block, real& globalspacing)
{
    int has_config = 0;
    size_t fault = 0;

    if(block.subshape == "Hollow")
    {
        block.sub_bound_type = hollow;

        // Need wall thickness defining, and length of cylinder. Both 2D and 3D
        if(block.thickness < 0 && block.nk < 0)
        {
            printf("ERROR: Cylinder block \"%s\" thickness or wall count has not been defined\n",block.name.c_str());
            fault = 1;
        }

        if(block.length < 0 && block.nj < 0)
        {
            printf("ERROR: Cylinder block \"%s\" legnth or j-count has not been defined\n",block.name.c_str());
            fault = 1;
        }

        // Now need to define the geometry of the plane
        if(check_vector(block.centre))
        {   /* Default to these values first */
            if(block.radius > 0)
            {
                has_config = 3;
                block.start = block.centre;
                block.start[1] -= block.radius;

                block.end = block.centre;
                block.end[1] += block.radius;

                #if SIMDIM == 3
                block.start[2] -= block.radius;
                block.end[2] += block.radius;
                block.right = block.centre;
                block.right[1] += block.radius;
                block.right[2] -= block.radius;
                #endif
            }
            #if SIMDIM == 3 //Only valid for 3D
            else if(check_vector(block.right) && check_vector(block.end))
            {
                has_config = 4;
                
                StateVecD v = block.right - block.centre;
                StateVecD u = block.end - block.right;
                block.start = block.centre - v - u;
                block.right = block.start + 2.0 * v;
                block.end = block.right + 2.0 * u;
                block.radius = v.norm();
            }
            #endif
        }
        else if (check_vector(block.start) && check_vector(block.end) 
            #if SIMDIM == 3
            && check_vector(block.right) 
            #endif
            )
        {
            has_config = 1;
            block.centre = 0.5 * (block.start + block.end);
            StateVecD v = block.right - block.start;
            block.radius = 0.5 * v.norm();
            block.right = block.centre + 0.5 * v;
        }
        
    }
    else if (block.subshape == "Solid")
    {
        block.sub_bound_type = solid;

        if(block.length < 0 || block.nk < 0)
        {
            printf("ERROR: Cylinder block \"%s\" length or k-count has not been defined\n",block.name.c_str());
            fault = 1;
        }

        if (check_vector(block.start) && check_vector(block.right) && check_vector(block.end))
        {
            has_config = 1;
            block.centre = 0.5 * (block.start + block.end);
            StateVecD v = block.right - block.start;
            block.radius = 0.5 * v.norm();
        }
        else if(check_vector(block.centre))
        {   /* Default to these values first */
            if(block.radius > 0)
            {
                has_config = 3;
                block.start = block.centre;
                block.start[1] -= block.radius;

                block.end = block.centre;
                block.end[1] += block.radius;

                #if SIMDIM == 3
                block.start[2] -= block.radius;
                block.end[2] += block.radius;
                block.right = block.centre;
                block.right[1] += block.radius;
                block.right[2] -= block.radius;
                #endif
            }

            if(check_vector(block.right) && check_vector(block.end))
            {
                has_config = 4;
                // Need to find the start
                StateVecD v = block.right - block.centre;
                StateVecD u = block.end - block.right;
                block.start = block.centre - v - u;
                block.right = block.start + 2.0 * v;
                block.end = block.right + 2.0 * u;
                block.radius = v.norm();
            }
        }
        else if (check_vector(block.start) && block.ni > 0 && block.nj > 0)
        {
            has_config = 2;
        }
        
    }
    else
    {
        printf("ERROR: Cylinder block \"%s\" subtype not defined appropriately. Choose from the following options:\n",
            block.name.c_str());
        printf("\t1. Hollow\n\t2. Solid\n");
        fault = 1;
    }
    
    if(has_config == 0)
    {
        printf("ERROR: Cylinder block \"%s\" geometry has not been sufficiently defined.\n",block.name.c_str());
        fault = 1;
    }

    // Generate the rotation matrix if it exists
    if(block.angles.norm() != 0)
    {
        block.angles *= M_PI/180.0; // Convert to radians

        // Find the rotation matrix
        #if SIMDIM == 3
            StateMatD rotx, roty, rotz;
            rotx << 1.0 , 0.0                   , 0.0                  ,
                    0.0 , cos(block.angles(0))  , sin(block.angles(0)) ,
                    0.0 , -sin(block.angles(0)) , cos(block.angles(0)) ;

            roty << cos(block.angles(1)) , 0.0 , -sin(block.angles(1)) ,
                    0.0                  , 1.0 , 0.0                   ,
                    sin(block.angles(1)) , 0.0 , cos(block.angles(1))  ;

            rotz << cos(block.angles(2)) , -sin(block.angles(2)) , 0.0 ,
                    sin(block.angles(2)) , cos(block.angles(2))  , 0.0 ,
                    0.0                  , 0.0                   , 1.0 ;

            block.rotmat = rotx*roty*rotz;
        #else
            StateMatD rotmat;
            rotmat << cos(block.angles(0)), -sin(block.angles(0)),
                      sin(block.angles(0)),  cos(block.angles(0));

            block.rotmat = rotmat;
        #endif
        block.normal = StateVecD::Zero();
        block.normal[0] = 1.0;
        block.normal = block.rotmat * block.normal;
    }   
    else if(block.normal != default_norm)
    {
        StateMatD rotmat = StateMatD::Identity();
        // Need to find the rotation matrix now 
        // (will lose some information, as no knowledge of rotation around normal)
        block.normal.normalize();

        #if SIMDIM == 3
        StateVecD origin = default_norm;
        StateVecD v = origin.cross(block.normal);
        // If magnitude of the cross is small, normal is essentially parallel
        if(v.norm() > 1e-10)
        {
            real mag = origin.dot(block.normal);
            real s = v.norm();

            StateMatD k;
            k <<  0.0 , -v[2],  v[1],
                  v[2],  0   , -v[0],
                 -v[1],  v[0],  0;

            rotmat += k + k * k * (1.0 - mag) / (s * s);
        }
        #else
        // Find rotation matrix
        block.angles[0] = atan2(block.normal[1],block.normal[0]);
        rotmat << cos(block.angles(0)), sin(block.angles(0)),
            -sin(block.angles(0)),  cos(block.angles(0));
        #endif
        block.rotmat = rotmat;
    } 

    if(has_config == 1)
    {
        //Have three points to define the plane. Override any normal and rotation matrix def
        StateVecD ab = (block.end - block.start).normalized();
        StateMatD rotmat;
        #if SIMDIM == 3
        StateVecD ac = (block.right - block.start).normalized();
        block.normal = (ab.cross(ac)).normalized();
        rotmat << ab[0]          , ab[1]          , ab[2]          ,
                  ac[0]          , ac[1]          , ac[2]          , 
                  block.normal[0], block.normal[1], block.normal[2];
        #else
        block.angles[0] = atan2(ab[1],ab[0]);
        rotmat << cos(block.angles(0)), sin(block.angles(0)),
            -sin(block.angles(0)),  cos(block.angles(0));
        block.normal = rotmat * StateVecD(1.0,0.0);
        #endif
        block.rotmat = rotmat;

    }
    if(has_config == 2)
    {
        // Start and ni and nj counts have been defined
        block.centre = block.start;
        block.radius = 0.5*globalspacing*block.ni;
        
        block.centre[1] += block.radius;
        #if SIMDIM == 3
        block.right = block.start;
        block.right[1] += globalspacing * block.ni;
        block.end = block.right;
        block.end[2] += globalspacing * block.nj;
        block.centre[2] += block.radius;
        #else
        block.end = block.start;
        block.end[1] += globalspacing * block.ni;
        #endif

        if(block.rotmat != StateMatD::Identity())
        {
            block.centre = block.rotmat * (block.centre - block.start) + block.start;
            block.end = block.rotmat * (block.end - block.start) + block.start;
            #if SIMDIM == 3
            block.right = block.rotmat * (block.right - block.start) + block.start;
            #endif
        }
    }
    else if (has_config == 3 || has_config == 4)
    {
        if(block.rotmat != StateMatD::Identity())
        {
            block.start = block.rotmat * (block.start - block.centre) + block.centre;
            block.end = block.rotmat * (block.end - block.centre) + block.centre;
            #if SIMDIM == 3
            block.right = block.rotmat * (block.right - block.centre) + block.centre;
            #endif
        }
    }
    
    
    block.dx = globalspacing; // Potential to allow different size particles, but not right now.
    
    // if(block.dx < 0)
    // {
    //     if(block.ni < 0)
    //     {
    //         printf("WARNING: Cylinder block \"%s\" globalspacing has not been defined. Using global globalspacing.\n",block.name.c_str());
    //         block.dx = globalspacing;
    //     }
    //     else
    //     {
    //         /* Use ni as a number along the diameter, so define globalspacing off that */
    //         real di = (2.0*block.radius)/real(block.ni);
    //         block.dx = di;
    //     }
    // }
    
    if(block.sub_bound_type == hollow)
    {
        if(block.thickness < 0)
        {
            if(block.hcpl == 1)
                block.thickness = static_cast<double>(block.nk) * globalspacing * sqrt(3.0) * 2.0;
            else
                block.thickness = static_cast<double>(block.nk) * globalspacing;
        }
        else if (block.nk < 0)
        {   /* Use a particle count over the real thickness value */
            if(block.hcpl == 1)
                block.nk = static_cast<int>(ceil(block.thickness / (globalspacing * sqrt(3.0) * 2.0)));
            else
                block.nk = static_cast<int>(ceil(block.thickness / globalspacing));
            
            block.nk = block.nk > 1 ? block.nk : 1;
        }

        #if SIMDIM == 3
        real dtheta = globalspacing/block.radius;

        int ni = static_cast<int>(ceil((2*M_PI)/dtheta));
        ni = ni > 1 ? ni : 1;
        block.ni = ni;
        #endif

        int nj = ceil(block.length/globalspacing)+1;
        block.nj = nj > 1 ? nj : 1; 

        block.npts = block.nj * block.nk;
        #if SIMDIM == 3
        block.npts *= block.ni;
        #else
        block.npts *= 2;
        #endif
    }
    else if (block.sub_bound_type == solid)
    {
        block.ni = ceil(2.0*block.radius/globalspacing);
        block.nj = ceil(block.length/globalspacing);
        #if SIMDIM == 3
        block.nk = ceil(2.0*block.radius/globalspacing);
        #endif
        block.ni = block.ni > 1 ? block.ni : 1;
        block.nj = block.nj > 1 ? block.nj : 1;
        block.nk = block.nk > 1 ? block.nk : 1;
        
        #if SIMDIM == 3
        block.npts = block.nj * block.nk * M_PI_4 * block.ni; 
        block.npts *= block.ni;
        #else
        block.npts = block.ni * block.nj;
        #endif
    }

    

    if(fault)
    {
        printf("Check of cylinder block \"%s\" finished with errors. Stopping.\n",block.name.c_str());
        exit(-1);
    }
}

#if SIMDIM == 2
std::vector<StateVecD> create_hollow_cylinder(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 


    StateVecD norm = block.normal.normalized();
    StateVecD left(norm[1],-norm[0]);
    real r = block.radius;
    for(int ii = 0; ii < block.nj; ii++)
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
    for(int ii = 0; ii < block.nj; ii++)
    {   // Depth
        for(int kk = 0; kk < block.nk; kk++)
        {   //Thickness
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = globalspacing * (norm * real(ii + 0.5*(kk % 2)) - 0.5 * left * sqrt(3.0) * real(kk))
                    - r * left;
            else
                newPoint = globalspacing * (norm * real(ii) - left * real(kk)) - r * left;

            newPoint += StateVecD(unif(re),unif(re));
            newPoint += block.centre;
            points.push_back(newPoint);
        }
    }
    return points;   
}

std::vector<StateVecD> create_solid_cylinder(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    StateVecD norm = -block.normal.normalized();
    StateVecD left(-norm[1],norm[0]);

    for (int jj = 0; jj < block.nj; ++jj)
    {
        for (int ii = 0; ii > block.ni; ++ii)
        {
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = norm * real(ii + 0.5*(jj % 2)) + left * sqrt(3.0) * real(jj);
            else
                newPoint = norm * real(jj) + left * real(ii);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.emplace_back(newPoint);

        } /* end i count */
    } /* end k count */

    return points;
}
#endif

#if SIMDIM == 3
std::vector<StateVecD> create_hollow_cylinder(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    // Create local coordinate system (u,v,w)
    double dtheta = 2*M_PI/real(block.ni); // Redefine dtheta so that it makes a perfect circle

    for(int jj = 0; jj < block.nj; jj++)
    {   // Depth
        for(int kk = 0; kk < block.nk; kk++)
        {   //Thickness
            real r = block.radius;
            real l = -real(jj)*globalspacing;
            real doffset = 0;
            if(block.hcpl == 1)
            {
                r += 1.0 / 3.0 * sqrt(6.0) * real(kk) * globalspacing;
                l = 0.5*sqrt(3)*(real(jj) + real(kk%2)/3.0)*globalspacing;
                doffset = 0.5*dtheta*((kk+jj)%2);
            }
            else
                r += real(kk)*globalspacing;

            for(int ii = 0; ii < block.ni; ii++)
            {   // Radius
                real theta = static_cast<real>(ii)*dtheta;
                // real l = real(ii)*globalspacing;
                if(block.hcpl == 1)
                {   /* Shift by half a globalspacing if jj is even */
                    theta += doffset;
                }

                real a = cos(theta)*r; real b = sin(theta)*r;
                StateVecD newPoint(l,a,b);
                newPoint += StateVecD(unif(re),unif(re),unif(re));
                newPoint = block.rotmat * newPoint;
                newPoint += block.centre;
                points.push_back(newPoint);
            }
        }
    }
    return points;   
}

std::vector<StateVecD> create_solid_cylinder(shape_block const& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    
    real const& radius = block.radius;
    real const rsq = radius*radius;
    // StateVecD const& centre = block.centre;
    for (int jj = 0; jj < block.nj; ++jj)
    {
        for (int kk = 0; kk > block.nk; ++kk)
        {
            for (int ii = 0; ii < block.ni; ++ii)
            {
                StateVecD newPoint;
                if (block.hcpl == 1)
                    newPoint = 0.5 * StateVecD(sqrt(3.0) * (real(kk) + real((jj % 2)) / 3.0),
                         real(2 * ii + ((kk + jj) % 2)),                         
                        2.0 / 3.0 * sqrt(6.0) * real(kk));
                else
                    newPoint = StateVecD(real(ii), real(kk), real(jj));

                newPoint *= globalspacing;
                newPoint += StateVecD(unif(re),unif(re),unif(re));

                real a = newPoint[1] - radius;
                real b = newPoint[2] - radius;

                if((a*a + b*b) > rsq)
                    continue;

                newPoint  = block.rotmat * newPoint;
                newPoint += block.start;
                points.emplace_back(newPoint);
            }   /* end i count */
        }   /* end j count */
    }   /* end k count */
    return points;
}
#endif


std::vector<StateVecD> create_cylinder(shape_block const& block, real const& globalspacing)
{
    if(block.sub_bound_type == hollow)
    {
        return create_hollow_cylinder(block,globalspacing);
    }
    else if (block.sub_bound_type == solid)
    {
        return create_solid_cylinder(block,globalspacing);
    }
    return std::vector<StateVecD>();
}