/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "cylinder.h"
#include "../Geometry.h"

enum geometry_config
{
    NO_CONFIG = 0,
    CentreRadius,
    CentreEnd,
    ThreePoints,
    StartCount
};

void CylinderShape::check_input(SIM const& svar, int& fault)
{
    // Do common input checks.
    ShapeBlock::check_input(svar, fault);

    int has_config = NO_CONFIG;
    if (subshape == "Hollow")
    {
        sub_bound_type = hollow;

        // Need wall thickness defining, and length of cylinder. Both 2D and 3D
        if (thickness < 0 && nk < 0)
        {
            printf(
                "ERROR: Cylinder block \"%s\" thickness or wall count has not been defined\n",
                name.c_str()
            );
            fault = 1;
        }

        if (length < 0 && nj < 0)
        {
            printf(
                "ERROR: Cylinder block \"%s\" legnth or j-count has not been defined\n", name.c_str()
            );
            fault = 1;
        }

        // Now need to define the geometry of the plane
        if (check_vector(centre))
        { /* Default to these values first */
            if (radius > 0)
            {
                has_config = CentreRadius;
                start = centre;
                start[1] -= radius;

                end = centre;
                end[1] += radius;

#if SIMDIM == 3
                start[2] -= radius;
                end[2] += radius;
                right = centre;
                right[1] += radius;
                right[2] -= radius;
#endif
            }
#if SIMDIM == 3 // Only valid for 3D
            else if (check_vector(right) && check_vector(end))
            {
                has_config = CentreEnd;

                StateVecD v = right - centre;
                StateVecD u = end - right;
                start = centre - v - u;
                right = start + 2.0 * v;
                end = right + 2.0 * u;
                radius = v.norm();
            }
#endif
        }
        else if (check_vector(start) && check_vector(end)
#if SIMDIM == 3
                 && check_vector(right)
#endif
        )
        {
            has_config = ThreePoints;
            centre = 0.5 * (start + end);
            StateVecD v = right - start;
            radius = 0.5 * v.norm();
            right = centre + 0.5 * v;
        }
    }
    else if (subshape == "Solid")
    {
        sub_bound_type = solid;

        if (length < 0 || nk < 0)
        {
            printf(
                "ERROR: Cylinder block \"%s\" length or k-count has not been defined\n", name.c_str()
            );
            fault = 1;
        }

        if (check_vector(start) && check_vector(right) && check_vector(end))
        {
            has_config = ThreePoints;
            centre = 0.5 * (start + end);
            StateVecD v = right - start;
            radius = 0.5 * v.norm();
        }
        else if (check_vector(centre))
        { /* Default to these values first */
            if (radius > 0)
            {
                has_config = CentreRadius;
                start = centre;
                start[1] -= radius;

                end = centre;
                end[1] += radius;

#if SIMDIM == 3
                start[2] -= radius;
                end[2] += radius;
                right = centre;
                right[1] += radius;
                right[2] -= radius;
#endif
            }

            if (check_vector(right) && check_vector(end))
            {
                has_config = CentreEnd;
                // Need to find the start
                StateVecD v = right - centre;
                StateVecD u = end - right;
                start = centre - v - u;
                right = start + 2.0 * v;
                end = right + 2.0 * u;
                radius = v.norm();
            }
        }
        else if (check_vector(start) && ni > 0 && nj > 0)
        {
            has_config = StartCount;
        }
    }
    else
    {
        printf(
            "ERROR: Cylinder block \"%s\" subtype not defined appropriately. Choose from the "
            "following "
            "options:\n",
            name.c_str()
        );
        printf("\t1. Hollow\n\t2. Solid\n");
        fault = 1;
    }

    // TODO - Need to revisit as this is almost certainly not adequate for all configs.
    if (dx < 0)
    {
        if (ni < 0)
        {
            std::cout << "Error: Block \"" << name << "\" spacing has not been defined." << std::endl;
            fault = 1;
        }
        else
        {
            /* Use ni as a number along the diameter, so define dx off that */
            real di = (2.0 * radius) / real(ni);
            dx = di;
        }
    }

    if (has_config == NO_CONFIG)
    {
        printf(
            "ERROR: Cylinder block \"%s\" geometry has not been sufficiently defined.\n", name.c_str()
        );
        fault = 1;
    }

    // Generate the rotation matrix if it exists
    if (angles.norm() != 0)
    {
        angles *= M_PI / 180.0; // Convert to radians

        // Find the rotation matrix
        // #if SIMDIM == 3
        //     StateMatD rotx, roty, rotz;
        //     rotx << 1.0 , 0.0                   , 0.0                  ,
        //             0.0 , cos(angles(0))  , sin(angles(0)) ,
        //             0.0 , -sin(angles(0)) , cos(angles(0)) ;

        //     roty << cos(angles(1)) , 0.0 , -sin(angles(1)) ,
        //             0.0                  , 1.0 , 0.0                   ,
        //             sin(angles(1)) , 0.0 , cos(angles(1))  ;

        //     rotz << cos(angles(2)) , -sin(angles(2)) , 0.0 ,
        //             sin(angles(2)) , cos(angles(2))  , 0.0 ,
        //             0.0                  , 0.0                   , 1.0 ;

        //     rotmat = rotx*roty*rotz;
        // #else
        //     StateMatD rotmat;
        //     rotmat << cos(angles(0)), -sin(angles(0)),
        //               sin(angles(0)),  cos(angles(0));

        //     rotmat = rotmat;
        // #endif
        rotmat = GetRotationMat(angles);
        normal = StateVecD::UnitX();
        normal = rotmat * normal;
    }
    else if (normal != StateVecD::UnitX())
    {
        StateMatD rotmat = StateMatD::Identity();
        // Need to find the rotation matrix now
        // (will lose some information, as no knowledge of rotation around normal)
        normal.normalize();

#if SIMDIM == 3
        StateVecD origin = StateVecD::UnitX();
        StateVecD v = origin.cross(normal);
        // If magnitude of the cross is small, normal is essentially parallel
        if (v.norm() > 1e-10)
        {
            real mag = origin.dot(normal);
            real s = v.norm();

            StateMatD k;
            k << 0.0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

            rotmat += k + k * k * (1.0 - mag) / (s * s);
        }
#else
        // Find rotation matrix
        angles[0] = atan2(normal[1], normal[0]);
        rotmat << cos(angles(0)), sin(angles(0)), -sin(angles(0)), cos(angles(0));
#endif
        rotmat = rotmat;
    }

    if (has_config == ThreePoints)
    {
        // Have three points to define the plane. Override any normal and rotation matrix def
        StateVecD ab = (end - start).normalized();
        StateMatD rotmat;
#if SIMDIM == 3
        StateVecD ac = (right - start).normalized();
        normal = (ab.cross(ac)).normalized();
        rotmat << ab[0], ab[1], ab[2], ac[0], ac[1], ac[2], normal[0], normal[1], normal[2];
#else
        angles[0] = atan2(ab[1], ab[0]);
        rotmat << cos(angles(0)), sin(angles(0)), -sin(angles(0)), cos(angles(0));
        normal = rotmat * StateVecD(1.0, 0.0);
#endif
        rotmat = rotmat;
    }
    if (has_config == StartCount)
    {
        // Start and ni and nj counts have been defined
        centre = start;
        radius = 0.5 * dx * ni;

        centre[1] += radius;
#if SIMDIM == 3
        right = start;
        right[1] += dx * ni;
        end = right;
        end[2] += dx * nj;
        centre[2] += radius;
#else
        end = start;
        end[1] += dx * ni;
#endif

        if (rotmat != StateMatD::Identity())
        {
            centre = rotmat * (centre - start) + start;
            end = rotmat * (end - start) + start;
#if SIMDIM == 3
            right = rotmat * (right - start) + start;
#endif
        }
    }
    else if (has_config == CentreRadius || has_config == CentreEnd)
    {
        if (rotmat != StateMatD::Identity())
        {
            start = rotmat * (start - centre) + centre;
            end = rotmat * (end - centre) + centre;
#if SIMDIM == 3
            right = rotmat * (right - centre) + centre;
#endif
        }
    }

    if (sub_bound_type == hollow)
    {
        if (thickness < 0)
        {
            if (particle_order == 1)
                thickness = static_cast<double>(nk) * dx * sqrt(3.0) * 2.0;
            else
                thickness = static_cast<double>(nk) * dx;
        }
        else if (nk < 0)
        { /* Use a particle count over the real thickness value */
            if (particle_order == 1)
                nk = static_cast<int>(ceil(thickness / (dx * sqrt(3.0) * 2.0)));
            else
                nk = static_cast<int>(ceil(thickness / dx));

            nk = nk > 1 ? nk : 1;
        }

#if SIMDIM == 3
        real dtheta = dx / radius;

        ni = static_cast<int>(ceil((2 * M_PI) / dtheta));
        ni = ni > 1 ? ni : 1;
#endif

        nj = ceil(length / dx) + 1;
        nj = nj > 1 ? nj : 1;

        npts = nj * nk;
#if SIMDIM == 3
        npts *= ni;
#else
        npts *= 2;
#endif
    }
    else if (sub_bound_type == solid)
    {
        ni = ceil(2.0 * radius / dx);
        nj = ceil(length / dx);
#if SIMDIM == 3
        nk = ceil(2.0 * radius / dx);
#endif
        ni = ni > 1 ? ni : 1;
        nj = nj > 1 ? nj : 1;
        nk = nk > 1 ? nk : 1;

#if SIMDIM == 3
        npts = nj * nk * M_PI_4 * ni;
        npts *= ni;
#else
        npts = ni * nj;
#endif
    }

    ShapeBlock::check_input_post();
}

#if SIMDIM == 2
std::vector<StateVecD> create_hollow_cylinder(ShapeBlock const& block, real const& dx)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    StateVecD norm = block.normal.normalized();
    StateVecD left(norm[1], -norm[0]);
    real r = block.radius;
    for (int ii = 0; ii < block.nj; ii++)
    { // Depth
        for (int kk = 0; kk < block.nk; kk++)
        { // Thickness
            StateVecD newPoint;
            if (block.particle_order == 1)
                newPoint = dx * (norm * real(-ii + 0.5 * (kk % 2)) + 0.5 * left * sqrt(3.0) * real(kk)) +
                           r * left;
            else
                newPoint = dx * (norm * real(-ii) + left * real(kk)) + r * left;

            newPoint += StateVecD(unif(re), unif(re));
            newPoint += block.centre;
            points.push_back(newPoint);
        }
    }
    // Other wall, but minus radius
    for (int ii = 0; ii < block.nj; ii++)
    { // Depth
        for (int kk = 0; kk < block.nk; kk++)
        { // Thickness
            StateVecD newPoint;
            if (block.particle_order == 1)
                newPoint = dx * (norm * real(-ii + 0.5 * (kk % 2)) - 0.5 * left * sqrt(3.0) * real(kk)) -
                           r * left;
            else
                newPoint = dx * (norm * real(-ii) - left * real(kk)) - r * left;

            newPoint += StateVecD(unif(re), unif(re));
            newPoint += block.centre;
            points.push_back(newPoint);
        }
    }
    return points;
}

std::vector<StateVecD> create_solid_cylinder(ShapeBlock const& block, real const& dx)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    StateVecD norm = -block.normal.normalized();
    StateVecD left(-norm[1], norm[0]);

    for (int jj = 0; jj < block.nj; ++jj)
    {
        for (int ii = 0; ii > block.ni; ++ii)
        {
            StateVecD newPoint;
            if (block.particle_order == 1)
                newPoint = norm * real(-ii + 0.5 * (jj % 2)) + left * sqrt(3.0) * real(jj);
            else
                newPoint = norm * real(-ii) + left * real(jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.emplace_back(newPoint);

        } /* end i count */
    } /* end k count */

    return points;
}
#endif

#if SIMDIM == 3
std::vector<StateVecD> create_hollow_cylinder(ShapeBlock const& block, real const& dx)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    // Create local coordinate system (u,v,w)
    double dtheta = 2 * M_PI / real(block.ni); // Redefine dtheta so that it makes a perfect circle

    for (int jj = 0; jj < block.nj; jj++)
    { // Depth
        for (int kk = 0; kk < block.nk; kk++)
        { // Thickness
            real r = block.radius;
            real l = real(jj) * dx;
            real doffset = 0;
            if (block.particle_order == 1)
            {
                r += 1.0 / 3.0 * sqrt(6.0) * real(kk) * dx;
                l = 0.5 * sqrt(3) * (real(jj) + real(kk % 2) / 3.0) * dx;
                doffset = 0.5 * dtheta * ((kk + jj) % 2);
            }
            else
                r += real(kk) * dx;

            for (int ii = 0; ii < block.ni; ii++)
            { // Radius
                real theta = static_cast<real>(ii) * dtheta;
                // real l = real(ii)*dx;
                if (block.particle_order == 1)
                { /* Shift by half a dx if jj is even */
                    theta += doffset;
                }

                real a = cos(theta) * r;
                real b = sin(theta) * r;
                StateVecD newPoint(-l, a, b);
                newPoint += StateVecD(unif(re), unif(re), unif(re));
                newPoint = block.rotmat * newPoint;
                newPoint += block.centre;
                points.push_back(newPoint);
            }
        }
    }
    return points;
}

std::vector<StateVecD> create_solid_cylinder(ShapeBlock const& block, real const& dx)
{
    std::vector<StateVecD> points;
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    real const& radius = block.radius;
    real const rsq = radius * radius;
    // StateVecD const& centre = block.centre;
    for (int jj = 0; jj < block.nj; ++jj)
    {
        for (int kk = 0; kk > block.nk; ++kk)
        {
            for (int ii = 0; ii < block.ni; ++ii)
            {
                StateVecD newPoint;
                if (block.particle_order == 1)
                    newPoint = 0.5 * StateVecD(
                                         -sqrt(3.0) * (real(kk) + real((jj % 2)) / 3.0),
                                         real(2 * ii + ((kk + jj) % 2)), 2.0 / 3.0 * sqrt(6.0) * real(kk)
                                     );
                else
                    newPoint = StateVecD(-real(ii), real(kk), real(jj));

                newPoint *= dx;
                newPoint += StateVecD(unif(re), unif(re), unif(re));

                real a = newPoint[1] - radius;
                real b = newPoint[2] - radius;

                if ((a * a + b * b) > rsq)
                    continue;

                newPoint = block.rotmat * newPoint;
                newPoint += block.start;
                points.emplace_back(newPoint);
            } /* end i count */
        } /* end j count */
    } /* end k count */
    return points;
}
#endif

void CylinderShape::generate_points()
{
    if (sub_bound_type == hollow)
    {
        coords = create_hollow_cylinder(*this, dx);
    }
    else if (sub_bound_type == solid)
    {
        coords = create_solid_cylinder(*this, dx);
    }
}