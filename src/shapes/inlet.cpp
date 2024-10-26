/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "inlet.h"
#include "../Geometry.h"

void InletShape::check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault)
{
    // Do common input checks.
    ShapeBlock::check_input(svar, fvar, globalspacing, fault);

    int has_config = 0;

#if SIMDIM == 3
    if (subshape == "Square")
    {
        sub_bound_type = squareCube;
        if (check_vector(start) && check_vector(right) && check_vector(end))
        {
            has_config = 1;
        }

        if (check_vector(start) && ni > 0 && nj > 0)
            has_config = 2;
    }
    else if (subshape == "Circle")
    {
        sub_bound_type = circleSphere;
        if (check_vector(centre))
        { /* Default to these values first */
            if (check_vector(right) && check_vector(end))
            {
                has_config = 4;
                StateVecD v = right - centre;
                StateVecD u = end - right;
                start = centre - v - u;
                right = start + 2.0 * v;
                end = right + 2.0 * u;
            }

            if (radius > 0)
            {
                has_config = 3;
                start = centre;
                start[1] -= radius;
                start[2] -= radius;

                end = centre;
                end[1] += radius;
                end[2] += radius;

                right = centre;
                right[1] += radius;
                right[2] -= radius;
            }
        }
        else if (check_vector(start) && check_vector(right) && check_vector(end))
        {
            has_config = 1;
            centre = 0.5 * (start + end);
            StateVecD v = right - start;
            radius = 0.5 * v.norm();
            right = centre + 0.5 * v;
        }
    }
    else
    {
        printf(
            "ERROR: Inlet block \"%s\" sub shape type has not been correctly defined.\n", name.c_str()
        );
        printf("Please choose from:\n\t1. Square\n\t2. Circle\n");
        fault = 1;
    }
#else
    if (check_vector(start) && check_vector(end))
    {
        has_config = 1;
        centre = 0.5 * (start + end);
        radius = (centre - start).norm();
    }
    else if (check_vector(centre) && radius > 0)
    {
        has_config = 3;
        start = centre;
        start[1] -= radius;
        end = centre;
        end[1] += radius;
    }
#endif

    if (has_config == 0)
    {
        printf("ERROR: Inlet block \"%s\" geometry not sufficiently defined\n", name.c_str());
        fault = 1;
    }

    dx = globalspacing; // Potential to allow different size particles, but not right now.

    // Generate the rotation matrix if it exists
    if (angles.norm() != 0)
    {
        angles *= M_PI / 180.0; // Convert to radians

        // Find the rotation matrix
        // #if SIMDIM == 3
        // StateMatD rotx, roty, rotz;
        // rotx << 1.0, 0.0            , 0.0           ,
        //         0.0, cos(angles(0)) , sin(angles(0)),
        //         0.0, -sin(angles(0)), cos(angles(0));

        // roty << cos(angles(1)) , 0.0 , -sin(angles(1)),
        //         0.0            , 1.0 , 0.0            ,
        //         sin(angles(1)) , 0.0 , cos(angles(1));

        // rotz << cos(angles(2)) , -sin(angles(2)) , 0.0 ,
        //         sin(angles(2)), cos(angles(2)) , 0.0 ,
        //         0.0            , 0.0            , 1.0 ;

        // rotmat = rotx*roty*rotz;
        // #else
        //     StateMatD rotmat;
        //     rotmat << cos(angles(0)), -sin(angles(0)),
        //               sin(angles(0)),  cos(angles(0));

        //     rotmat = rotmat;
        // #endif
        rotmat = GetRotationMat(angles);
        normal = StateVecD::UnitX();
        normal = rotmat * normal;
        insert_norm = normal;
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
        insert_norm = normal;
    }
    else if (insert_norm != StateVecD::UnitX())
    {
        StateMatD rotmat = StateMatD::Identity();
        // Need to find the rotation matrix now
        // (will lose some information, as no knowledge of rotation around normal)
        normal = insert_norm;
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
        angles[0] = atan2(normal[1], normal[0]) - M_PI / 2.0;
        rotmat << cos(angles(0)), sin(angles(0)), -sin(angles(0)), cos(angles(0));
#endif
        rotmat = rotmat;
    }

    if (has_config == 1)
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
        normal = rotmat * StateVecD(0.0, 1.0);
#endif
        rotmat = rotmat;

        insert_norm = normal;
        // Define the insertion constant using the length
        StateVecD test = start - (length - 0.01 * globalspacing) * insert_norm;
        insconst = insert_norm.dot(test);
    }
    else if (has_config == 2)
    {
// Start and ni and nj counts have been defined (nk doesn't need to be defined for the inlet)
#if SIMDIM == 3
        right = start;
        right[1] += globalspacing * ni;
        end = right;
        end[2] += globalspacing * nj;
#else
        end = start;
        end[1] += globalspacing * ni;
#endif

        if (rotmat != StateMatD::Identity())
        {
            end = rotmat * (end - start) + start;
#if SIMDIM == 3
            right = rotmat * (right - start) + start;
#endif
        }
    }
    else if (has_config == 3 || has_config == 4)
    {
        // Centre and radius defined, so start, right and end derived.
        if (rotmat != StateMatD::Identity())
        {
            start = rotmat * (start - centre) + centre;
            end = rotmat * (end - centre) + centre;
#if SIMDIM == 3
            right = rotmat * (right - centre) + centre;
#endif
        }
    }

/* Estimate how many points on the Block */
// Need to have lengths along the axis to accurately have the counts
#if SIMDIM == 3
    real xlength = (right - start).norm();
    real ylength = (end - right).norm();
#else
    real xlength = (end - start).norm();
#endif

#if SIMDIM == 3
#endif

    ni = static_cast<int>(ceil(xlength / globalspacing));
#if SIMDIM == 3
    nj = static_cast<int>(ceil(ylength / globalspacing));
#endif
    nk = static_cast<int>(ceil(length / globalspacing));

    ni = ni > 1 ? ni : 1;
    nk = nk > 1 ? nk : 1;
    nk = nk;
    ni = ni;
    int nBuff = 4;
#if SIMDIM == 3
    nj = nj > 1 ? nj : 1;
    nj = nj;
#endif

#if SIMDIM == 3
    if (sub_bound_type == circleSphere)
    {
        // Find the centre if it's not already defined
        if (!check_vector(centre))
        {
            centre = 0.5 * (start + end);
        }

        if (radius < 0)
        {
            radius = (centre - start).norm();
            radius *= cos(M_PI_4);
        }
        npts = ni * nj * M_PI_4 * (nk + nBuff + 1);
    }
    else
    {
        npts = ni * nj * (nk + nBuff + 1);
    }
#endif

#if SIMDIM == 2
    npts = ni * (nk + nBuff);
#endif

    if (vmag != 0)
    {
        vel = vmag * insert_norm;
    }

    ShapeBlock::check_input_post(globalspacing);
}

#if SIMDIM == 2
void InletShape::generate_points(real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    StateVecD delta = (end - start) / static_cast<real>(ni);
    StateVecD norm(delta[1], -delta[0]);
    int jj = 0;
    for (jj = 0; jj < nk; ++jj)
    {
        for (int ii = 0; ii < ni; ++ii)
        {
            StateVecD newPoint;

            newPoint = delta * real(ii) + norm * real(-jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += start;
            points.emplace_back(newPoint);
            bc.emplace_back(PIPE);

        } /* end i count */
    } /* end k count */

    // Create the back particles
    for (int ii = 0; ii < ni; ++ii)
    {
        StateVecD newPoint;
        newPoint = delta * real(ii) + norm * real(-jj);

        newPoint += StateVecD::Constant(unif(re));
        newPoint += start;
        points.emplace_back(newPoint);
        bc.emplace_back(BACK);
        back.emplace_back(points.size() - 1);
    } /* end i count */

    // Buffer zone
    int nBuff = 4;
    buffer = vector<vector<size_t>>(ni, vector<size_t>(nBuff));
    size_t buff = 0;
    for (jj = nk + 1; jj <= nk + nBuff; ++jj)
    {
        for (int ii = 0; ii < ni; ++ii)
        {
            StateVecD newPoint;

            newPoint = delta * real(ii) + norm * real(-jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += start;
            points.emplace_back(newPoint);
            bc.emplace_back(BUFFER);
            buffer[ii][buff] = points.size() - 1;
        } /* end i count */
        buff++;
    } /* end k count */

    coords = points;
}

#else

std::vector<StateVecD>
create_radial_disk(ShapeBlock const& block, real const& dx, real const& radius, int const& kk)
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
create_lattice_disk(ShapeBlock const& block, real const& dx, real const& radius, int const& kk)
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
create_disk(ShapeBlock const& block, real const& dx, real const& radius, int const& kk)
{
    if (block.particle_order)
        return create_radial_disk(block, dx, radius, kk);
    else
        return create_lattice_disk(block, dx, radius, kk);
}

void InletShape::generate_points(real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON * globalspacing);
    std::default_random_engine re;

    // Which geometry to create?
    if (sub_bound_type == squareCube)
    {
        // Depth
        int kk = 0;
        for (kk = 0; kk < nk; ++kk)
        { // Width
            for (int jj = 0; jj < nj; ++jj)
            { // Height
                for (int ii = 0; ii < ni; ++ii)
                {
                    StateVecD newPoint = StateVecD(real(-kk), real(ii), real(jj));

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re), unif(re), unif(re));
                    newPoint = rotmat * newPoint;
                    newPoint += start;
                    points.emplace_back(newPoint);
                    bc.emplace_back(PIPE);
                } /* end i count */
            } /* end j count */
        } /* end k count */

        // Create the back particles
        for (int jj = 0; jj < nj; ++jj)
        {
            for (int ii = 0; ii < ni; ++ii)
            {
                StateVecD newPoint(real(-kk), real(jj), real(ii));

                newPoint *= globalspacing;
                newPoint += StateVecD(unif(re), unif(re), unif(re));
                newPoint = rotmat * newPoint;
                newPoint += start;
                points.emplace_back(newPoint);
                bc.emplace_back(BACK);
                back.emplace_back(points.size() - 1);

            } /* end i count */
        } /* end j count */

        // Create buffer zone
        int nBuff = 4;
        buffer = vector<vector<size_t>>(back.size(), vector<size_t>(nBuff));
        size_t buff = 0;
        for (kk = nk + 1; kk <= nk + nBuff; ++kk)
        {
            size_t part_id = 0;
            for (int jj = 0; jj < nj; ++jj)
            {
                for (int ii = 0; ii < ni; ++ii)
                {
                    StateVecD newPoint = StateVecD(real(-kk), real(ii), real(jj));

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re), unif(re), unif(re));
                    newPoint = rotmat * newPoint;
                    newPoint += start;
                    points.emplace_back(newPoint);
                    bc.emplace_back(BUFFER);
                    buffer[part_id][buff] = points.size() - 1;
                    part_id++;
                } /* end i count */
            } /* end j count */
            buff++;
        } /* end k count */
    }
    else if (sub_bound_type == circleSphere)
    {
        // StateVecD const& centre = centre;
        int kk = 0;
        for (kk = 0; kk < nk; ++kk)
        {
            std::vector<StateVecD> temp = create_disk(*this, globalspacing, radius, kk);

            for (StateVecD const& pos : temp)
            {
                points.emplace_back(pos);
                bc.emplace_back(PIPE);
            }
        } /* end k count */

        // Create the back particles
        std::vector<StateVecD> temp = create_disk(*this, globalspacing, radius, kk);

        for (StateVecD const& pos : temp)
        {
            back.emplace_back(points.size());
            points.emplace_back(pos);
            bc.emplace_back(BACK);
        }

        // Create buffer zone
        int nBuff = 4;
        buffer = vector<vector<size_t>>(back.size(), vector<size_t>(nBuff));
        size_t buff = 0;
        for (kk = nk + 1; kk <= nk + nBuff; ++kk)
        {
            std::vector<StateVecD> temp = create_disk(*this, globalspacing, radius, kk);

            for (size_t ii = 0; ii < temp.size(); ii++)
            {
                buffer[ii][buff] = points.size();
                points.emplace_back(temp[ii]);
                bc.emplace_back(BUFFER);
            }

            buff++;
        } /* end k count */
    }

    coords = points;
}

#endif

uint update_buffer_region(SIM& svar, LIMITS& limits, SPHState& pnp1, size_t& end)
{
    uint nAdd = 0;
    /*Check if more particles need to be created*/
    for (size_t block_id = svar.n_bound_blocks; block_id < svar.n_bound_blocks + svar.n_fluid_blocks;
         block_id++)
    {
        bound_block& block = limits[block_id];
        if (block.block_type == inletZone)
        {
            size_t& part_id = svar.part_id;
            /* Check if any buffer particles have become pipe particles */
            for (size_t ii = 0; ii < block.back.size(); ++ii)
            {
                size_t const& pID = block.back[ii];
                StateVecD const& pos = pnp1[pID].xi;
                if (pos.dot(block.insert_norm) > block.insconst)
                {
                    /* Old particle has cleared the back zone, a new particle can be added. */
                    pnp1[pID].b = PIPE;

                    /* Update the back vector */
                    pnp1[block.buffer[ii][0]].b = BACK;
                    block.back[ii] = block.buffer[ii][0];

                    /* Update the buffer vector */
                    for (size_t jj = 0; jj < block.buffer[0].size() - 1; ++jj)
                        block.buffer[ii][jj] = block.buffer[ii][jj + 1];

                    /* Create a new particle */
                    if (svar.total_points < svar.max_points)
                    {
                        StateVecD xi = pnp1[block.buffer[ii].back()].xi - svar.dx * block.insert_norm;

                        // Insert the new particle, using the properties of the old particle.
                        // Use insert to keep particles closer together in memory.
                        SPHPart new_part = SPHPart(xi, pnp1[block.buffer[ii].back()], BUFFER, part_id);
                        pnp1.insert(pnp1.begin() + block.index.second, new_part);

                        block.buffer[ii].back() = block.index.second;

                        // Adjust block indexes
                        block.index.second++;
                        // Also need to adjust all blocks after the current
                        for (size_t jj = block_id + 1; jj < svar.n_bound_blocks + svar.n_fluid_blocks;
                             jj++)
                        {
                            limits[jj].index.first++;
                            limits[jj].index.second++;
                        }

                        svar.fluid_points++;
                        svar.total_points++;
                        end++;
                        part_id++;
                        nAdd++;
                    }
                }
            }
        }
    }
    return nAdd;
}
