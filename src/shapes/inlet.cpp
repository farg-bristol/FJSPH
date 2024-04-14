/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "inlet.h"
#include "../Geometry.h"

void InletShape::check_input(shape_block& block, real& globalspacing, int& fault)
{
    int has_config = 0;

#if SIMDIM == 3
    if (block.subshape == "Square")
    {
        block.sub_bound_type = squareCube;
        if (check_vector(block.start) && check_vector(block.right) && check_vector(block.end))
        {
            has_config = 1;
        }

        if (check_vector(block.start) && block.ni > 0 && block.nj > 0)
            has_config = 2;
    }
    else if (block.subshape == "Circle")
    {
        block.sub_bound_type = circleSphere;
        if (check_vector(block.centre))
        { /* Default to these values first */
            if (check_vector(block.right) && check_vector(block.end))
            {
                has_config = 4;
                StateVecD v = block.right - block.centre;
                StateVecD u = block.end - block.right;
                block.start = block.centre - v - u;
                block.right = block.start + 2.0 * v;
                block.end = block.right + 2.0 * u;
            }

            if (block.radius > 0)
            {
                has_config = 3;
                block.start = block.centre;
                block.start[1] -= block.radius;
                block.start[2] -= block.radius;

                block.end = block.centre;
                block.end[1] += block.radius;
                block.end[2] += block.radius;

                block.right = block.centre;
                block.right[1] += block.radius;
                block.right[2] -= block.radius;
            }
        }
        else if (check_vector(block.start) && check_vector(block.right) && check_vector(block.end))
        {
            has_config = 1;
            block.centre = 0.5 * (block.start + block.end);
            StateVecD v = block.right - block.start;
            block.radius = 0.5 * v.norm();
            block.right = block.centre + 0.5 * v;
        }
    }
    else
    {
        printf(
            "ERROR: Inlet block \"%s\" sub shape type has not been correctly defined.\n",
            block.name.c_str()
        );
        printf("Please choose from:\n\t1. Square\n\t2. Circle\n");
        fault = 1;
    }
#else
    if (check_vector(block.start) && check_vector(block.end))
    {
        has_config = 1;
        block.centre = 0.5 * (block.start + block.end);
        block.radius = (block.centre - block.start).norm();
    }
    else if (check_vector(block.centre) && block.radius > 0)
    {
        has_config = 3;
        block.start = block.centre;
        block.start[1] -= block.radius;
        block.end = block.centre;
        block.end[1] += block.radius;
    }
#endif

    if (has_config == 0)
    {
        printf("ERROR: Inlet block \"%s\" geometry not sufficiently defined\n", block.name.c_str());
        fault = 1;
    }

    block.dx = globalspacing; // Potential to allow different size particles, but not right now.

    // Generate the rotation matrix if it exists
    if (block.angles.norm() != 0)
    {
        block.angles *= M_PI / 180.0; // Convert to radians

        // Find the rotation matrix
        // #if SIMDIM == 3
        // StateMatD rotx, roty, rotz;
        // rotx << 1.0, 0.0            , 0.0           ,
        //         0.0, cos(block.angles(0)) , sin(block.angles(0)),
        //         0.0, -sin(block.angles(0)), cos(block.angles(0));

        // roty << cos(block.angles(1)) , 0.0 , -sin(block.angles(1)),
        //         0.0            , 1.0 , 0.0            ,
        //         sin(block.angles(1)) , 0.0 , cos(block.angles(1));

        // rotz << cos(block.angles(2)) , -sin(block.angles(2)) , 0.0 ,
        //         sin(block.angles(2)), cos(block.angles(2)) , 0.0 ,
        //         0.0            , 0.0            , 1.0 ;

        // block.rotmat = rotx*roty*rotz;
        // #else
        //     StateMatD rotmat;
        //     rotmat << cos(block.angles(0)), -sin(block.angles(0)),
        //               sin(block.angles(0)),  cos(block.angles(0));

        //     block.rotmat = rotmat;
        // #endif
        block.rotmat = GetRotationMat(block.angles);
        block.normal = StateVecD::UnitX();
        block.normal = block.rotmat * block.normal;
        block.insert_norm = block.normal;
    }
    else if (block.normal != StateVecD::UnitX())
    {
        StateMatD rotmat = StateMatD::Identity();
        // Need to find the rotation matrix now
        // (will lose some information, as no knowledge of rotation around normal)
        block.normal.normalize();

#if SIMDIM == 3
        StateVecD origin = StateVecD::UnitX();
        StateVecD v = origin.cross(block.normal);
        // If magnitude of the cross is small, normal is essentially parallel
        if (v.norm() > 1e-10)
        {
            real mag = origin.dot(block.normal);
            real s = v.norm();

            StateMatD k;
            k << 0.0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

            rotmat += k + k * k * (1.0 - mag) / (s * s);
        }
#else
        // Find rotation matrix
        block.angles[0] = atan2(block.normal[1], block.normal[0]);
        rotmat << cos(block.angles(0)), sin(block.angles(0)), -sin(block.angles(0)),
            cos(block.angles(0));
#endif
        block.rotmat = rotmat;
        block.insert_norm = block.normal;
    }
    else if (block.insert_norm != StateVecD::UnitX())
    {
        StateMatD rotmat = StateMatD::Identity();
        // Need to find the rotation matrix now
        // (will lose some information, as no knowledge of rotation around normal)
        block.normal = block.insert_norm;
        block.normal.normalize();

#if SIMDIM == 3
        StateVecD origin = StateVecD::UnitX();
        StateVecD v = origin.cross(block.normal);
        // If magnitude of the cross is small, normal is essentially parallel
        if (v.norm() > 1e-10)
        {
            real mag = origin.dot(block.normal);
            real s = v.norm();

            StateMatD k;
            k << 0.0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

            rotmat += k + k * k * (1.0 - mag) / (s * s);
        }
#else
        // Find rotation matrix
        block.angles[0] = atan2(block.normal[1], block.normal[0]) - M_PI / 2.0;
        rotmat << cos(block.angles(0)), sin(block.angles(0)), -sin(block.angles(0)),
            cos(block.angles(0));
#endif
        block.rotmat = rotmat;
    }

    if (has_config == 1)
    {
        // Have three points to define the plane. Override any normal and rotation matrix def
        StateVecD ab = (block.end - block.start).normalized();
        StateMatD rotmat;
#if SIMDIM == 3
        StateVecD ac = (block.right - block.start).normalized();
        block.normal = (ab.cross(ac)).normalized();
        rotmat << ab[0], ab[1], ab[2], ac[0], ac[1], ac[2], block.normal[0], block.normal[1],
            block.normal[2];
#else
        block.angles[0] = atan2(ab[1], ab[0]);
        rotmat << cos(block.angles(0)), sin(block.angles(0)), -sin(block.angles(0)),
            cos(block.angles(0));
        block.normal = rotmat * StateVecD(0.0, 1.0);
#endif
        block.rotmat = rotmat;

        block.insert_norm = block.normal;
        // Define the insertion constant using the length
        StateVecD test = block.start - (block.length - 0.01 * globalspacing) * block.insert_norm;
        block.insconst = block.insert_norm.dot(test);
    }
    else if (has_config == 2)
    {
// Start and ni and nj counts have been defined (nk doesn't need to be defined for the inlet)
#if SIMDIM == 3
        block.right = block.start;
        block.right[1] += globalspacing * block.ni;
        block.end = block.right;
        block.end[2] += globalspacing * block.nj;
#else
        block.end = block.start;
        block.end[1] += globalspacing * block.ni;
#endif

        if (block.rotmat != StateMatD::Identity())
        {
            block.end = block.rotmat * (block.end - block.start) + block.start;
#if SIMDIM == 3
            block.right = block.rotmat * (block.right - block.start) + block.start;
#endif
        }
    }
    else if (has_config == 3 || has_config == 4)
    {
        // Centre and radius defined, so start, right and end derived.
        if (block.rotmat != StateMatD::Identity())
        {
            block.start = block.rotmat * (block.start - block.centre) + block.centre;
            block.end = block.rotmat * (block.end - block.centre) + block.centre;
#if SIMDIM == 3
            block.right = block.rotmat * (block.right - block.centre) + block.centre;
#endif
        }
    }

/* Estimate how many points on the Block */
// Need to have lengths along the axis to accurately have the counts
#if SIMDIM == 3
    real xlength = (block.right - block.start).norm();
    real ylength = (block.end - block.right).norm();
#else
    real xlength = (block.end - block.start).norm();
#endif

    size_t ni, nk;
#if SIMDIM == 3
    size_t nj;
#endif

    ni = static_cast<int>(ceil(xlength / globalspacing));
#if SIMDIM == 3
    nj = static_cast<int>(ceil(ylength / globalspacing));
#endif
    nk = static_cast<int>(ceil(block.length / globalspacing));

    ni = ni > 1 ? ni : 1;
    nk = nk > 1 ? nk : 1;
    block.nk = nk;
    block.ni = ni;
    int nBuff = 4;
#if SIMDIM == 3
    nj = nj > 1 ? nj : 1;
    block.nj = nj;
#endif

#if SIMDIM == 3
    if (block.sub_bound_type == circleSphere)
    {
        // Find the centre if it's not already defined
        if (!check_vector(block.centre))
        {
            block.centre = 0.5 * (block.start + block.end);
        }

        if (block.radius < 0)
        {
            block.radius = (block.centre - block.start).norm();
            block.radius *= cos(M_PI_4);
        }
        block.npts = block.ni * block.nj * M_PI_4 * (nk + nBuff + 1);
    }
    else
    {
        block.npts = block.ni * block.nj * (nk + nBuff + 1);
    }
#endif

#if SIMDIM == 2
    block.npts = block.ni * (block.nk + nBuff);
#endif

    if (block.vmag != 0)
    {
        block.vel = block.vmag * block.insert_norm;
    }
}

#if SIMDIM == 2
std::vector<StateVecD> InletShape::generate_points(shape_block& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    StateVecD delta = (block.end - block.start) / static_cast<real>(block.ni);
    StateVecD norm(delta[1], -delta[0]);
    int jj = 0;
    for (jj = 0; jj < block.nk; ++jj)
    {
        for (int ii = 0; ii < block.ni; ++ii)
        {
            StateVecD newPoint;

            newPoint = delta * real(ii) + norm * real(-jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.emplace_back(newPoint);
            block.bc.emplace_back(PIPE);

        } /* end i count */
    } /* end k count */

    // Create the back particles
    for (int ii = 0; ii < block.ni; ++ii)
    {
        StateVecD newPoint;
        newPoint = delta * real(ii) + norm * real(-jj);

        newPoint += StateVecD::Constant(unif(re));
        newPoint += block.start;
        points.emplace_back(newPoint);
        block.bc.emplace_back(BACK);
        block.back.emplace_back(points.size() - 1);
    } /* end i count */

    // Buffer zone
    int nBuff = 4;
    block.buffer = vector<vector<size_t>>(block.ni, vector<size_t>(nBuff));
    size_t buff = 0;
    for (jj = block.nk + 1; jj <= block.nk + nBuff; ++jj)
    {
        for (int ii = 0; ii < block.ni; ++ii)
        {
            StateVecD newPoint;

            newPoint = delta * real(ii) + norm * real(-jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.emplace_back(newPoint);
            block.bc.emplace_back(BUFFER);
            block.buffer[ii][buff] = points.size() - 1;
        } /* end i count */
        buff++;
    } /* end k count */
    return points;
}

#else

std::vector<StateVecD>
create_radial_disk(shape_block const& block, real const& dx, real const& radius, int const& kk)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    // find the depth of the layer
    real x = -dx * kk;

    for (real rad = radius; rad > 0.99 * dx; rad -= dx)
    {
        // % Find spacing to have a well defined surface.
        real dtheta = atan(dx / rad);
        int ncirc = floor(abs(2.0 * M_PI / dtheta));
        dtheta = 2.0 * M_PI / real(ncirc);

        for (real theta = 0.0; theta < 2 * M_PI - 0.5 * dtheta; theta += dtheta)
        { /* Create a ring of points */
            real y = rad * sin(theta);
            real z = rad * cos(theta);

            StateVecD newPoint(x, y, z);
            newPoint += StateVecD(unif(re), unif(re), unif(re));
            newPoint = block.rotmat * newPoint;
            newPoint += block.centre;
            points.emplace_back(newPoint);
        }
    }

    /* Create centre point */
    StateVecD newPoint(x, 0.0, 0.0);
    newPoint += StateVecD(unif(re), unif(re), unif(re));
    newPoint = block.rotmat * newPoint;
    newPoint += block.centre;
    points.emplace_back(newPoint);
    return points;
}

std::vector<StateVecD>
create_lattice_disk(shape_block const& block, real const& dx, real const& radius, int const& kk)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * dx);
    std::default_random_engine re;

    real const rsq = radius * radius;
    for (int jj = 0; jj < block.nj; ++jj)
    {
        for (int ii = 0; ii < block.ni; ++ii)
        {
            StateVecD newPoint = StateVecD(-real(kk), real(ii), real(jj));

            newPoint *= dx;

            real a = newPoint[1] - radius;
            real b = newPoint[2] - radius;
            if ((a * a + b * b) > rsq)
                continue;

            newPoint += StateVecD(unif(re), unif(re), unif(re));
            newPoint = block.rotmat * newPoint;
            newPoint += block.start;
            points.emplace_back(newPoint);

        } /* end i count */
    } /* end j count */
    return points;
}

std::vector<StateVecD>
create_disk(shape_block const& block, real const& dx, real const& radius, int const& kk)
{
    if (block.particle_order)
        return create_radial_disk(block, dx, radius, kk);
    else
        return create_lattice_disk(block, dx, radius, kk);
}

std::vector<StateVecD> InletShape::generate_points(shape_block& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    // Which geometry to create?
    if (block.sub_bound_type == squareCube)
    {
        // Depth
        int kk = 0;
        for (kk = 0; kk < block.nk; ++kk)
        { // Width
            for (int jj = 0; jj < block.nj; ++jj)
            { // Height
                for (int ii = 0; ii < block.ni; ++ii)
                {
                    StateVecD newPoint = StateVecD(real(-kk), real(ii), real(jj));

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re), unif(re), unif(re));
                    newPoint = block.rotmat * newPoint;
                    newPoint += block.start;
                    points.emplace_back(newPoint);
                    block.bc.emplace_back(PIPE);
                } /* end i count */
            } /* end j count */
        } /* end k count */

        // Create the back particles
        for (int jj = 0; jj < block.nj; ++jj)
        {
            for (int ii = 0; ii < block.ni; ++ii)
            {
                StateVecD newPoint(real(-kk), real(jj), real(ii));

                newPoint *= globalspacing;
                newPoint += StateVecD(unif(re), unif(re), unif(re));
                newPoint = block.rotmat * newPoint;
                newPoint += block.start;
                points.emplace_back(newPoint);
                block.bc.emplace_back(BACK);
                block.back.emplace_back(points.size() - 1);

            } /* end i count */
        } /* end j count */

        // Create buffer zone
        int nBuff = 4;
        block.buffer = vector<vector<size_t>>(block.back.size(), vector<size_t>(nBuff));
        size_t buff = 0;
        for (kk = block.nk + 1; kk <= block.nk + nBuff; ++kk)
        {
            size_t partID = 0;
            for (int jj = 0; jj < block.nj; ++jj)
            {
                for (int ii = 0; ii < block.ni; ++ii)
                {
                    StateVecD newPoint = StateVecD(real(-kk), real(ii), real(jj));

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re), unif(re), unif(re));
                    newPoint = block.rotmat * newPoint;
                    newPoint += block.start;
                    points.emplace_back(newPoint);
                    block.bc.emplace_back(BUFFER);
                    block.buffer[partID][buff] = points.size() - 1;
                    partID++;
                } /* end i count */
            } /* end j count */
            buff++;
        } /* end k count */
    }
    else if (block.sub_bound_type == circleSphere)
    {
        real const& radius = block.radius;

        // StateVecD const& centre = block.centre;
        int kk = 0;
        for (kk = 0; kk < block.nk; ++kk)
        {
            std::vector<StateVecD> temp = create_disk(block, globalspacing, radius, kk);

            for (StateVecD const& pos : temp)
            {
                points.emplace_back(pos);
                block.bc.emplace_back(PIPE);
            }
        } /* end k count */

        // Create the back particles
        std::vector<StateVecD> temp = create_disk(block, globalspacing, radius, kk);

        for (StateVecD const& pos : temp)
        {
            block.back.emplace_back(points.size());
            points.emplace_back(pos);
            block.bc.emplace_back(BACK);
        }

        // Create buffer zone
        int nBuff = 4;
        block.buffer = vector<vector<size_t>>(block.back.size(), vector<size_t>(nBuff));
        size_t buff = 0;
        for (kk = block.nk + 1; kk <= block.nk + nBuff; ++kk)
        {
            std::vector<StateVecD> temp = create_disk(block, globalspacing, radius, kk);

            for (size_t ii = 0; ii < temp.size(); ii++)
            {
                block.buffer[ii][buff] = points.size();
                points.emplace_back(temp[ii]);
                block.bc.emplace_back(BUFFER);
            }

            buff++;
        } /* end k count */
    }
    return points;
}

#endif

uint update_buffer_region(SIM& svar, LIMITS& limits, SPHState& pnp1, size_t& end, size_t& end_ng)
{
    uint nAdd = 0;
    /*Check if more particles need to be created*/
    for (size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
    {
        if (limits[block].block_type == inletZone)
        {
            size_t& partID = svar.partID;
            /* Check if any buffer particles have become pipe particles */
            for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
            {
                size_t const& pID = limits[block].back[ii];
                StateVecD const& pos = pnp1[pID].xi;
                /*Check that the starting area is clear first...*/
                if (pos.dot(limits[block].insert_norm) > limits[block].insconst)
                {
                    /* particle has cleared the back zone */
                    pnp1[pID].b = PIPE;

                    /* Update the back vector */
                    pnp1[limits[block].buffer[ii][0]].b = BACK;
                    limits[block].back[ii] = limits[block].buffer[ii][0];

                    /* Update the buffer vector */
                    for (size_t jj = 0; jj < limits[block].buffer[0].size() - 1; ++jj)
                        limits[block].buffer[ii][jj] = limits[block].buffer[ii][jj + 1];

                    /* Create a new particle */
                    if (svar.totPts < svar.finPts)
                    {
                        StateVecD xi = pnp1[limits[block].buffer[ii].back()].xi -
                                       svar.dx * limits[block].insert_norm;

                        pnp1.insert(
                            pnp1.begin() + limits[block].index.second,
                            SPHPart(xi, pnp1[limits[block].buffer[ii].back()], BUFFER, partID)
                        );
                        limits[block].buffer[ii].back() = limits[block].index.second;

                        limits[block].index.second++;
                        // Also need to adjust all blocks after the current
                        for (size_t jj = block + 1; jj < svar.nbound + svar.nfluid; jj++)
                        {
                            limits[jj].index.first++;
                            limits[jj].index.second++;
                        }

                        svar.simPts++;
                        svar.totPts++;
                        end_ng++;
                        end++;
                        partID++;
                        nAdd++;
                    }
                }
            }
        }
    }
    return nAdd;
}
