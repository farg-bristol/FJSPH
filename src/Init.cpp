/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Init.h"
#include "IOFunctions.h"
#include "Kernel.h"
#include "Neighbours.h"

#include "shapes/shapes.h"

#include <cctype>
#include <filesystem>

// Case insensitive comparison for single character.
bool ichar_equals(char a, char b)
{
    return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b));
}

// Cast insensitive comparison for strings.
bool iequals(const std::string& a, const std::string& b)
{
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), ichar_equals);
}

void get_boundary_velocity(ShapeBlock* boundvar)
{
    size_t ntimes = boundvar->ntimes;
    if (ntimes != 0)
    {
        boundvar->vels.resize(ntimes - 1);
        for (size_t jj = 0; jj < ntimes - 1; jj++)
        {
            boundvar->vels[jj] = (boundvar->pos[jj + 1] - boundvar->pos[jj]) /
                                 (boundvar->times[jj + 1] - boundvar->times[jj]);
        }
    }
}

std::vector<StateVecD> Read_Geom_File(std::string const& filename)
{
    std::ifstream input_file(filename, std::ios::in);

    int nPts;
    input_file >> nPts;
    std::vector<StateVecD> points(nPts);

    for (int ii = 0; ii < nPts; ++ii)
    {
        for (int dim = 0; dim < SIMDIM; ++dim)
        {
            input_file >> points[ii][dim];
        }
    }

    input_file.close();

    return points;
}

void Check_Intersection(SIM const& svar, Shapes& boundvar, Shapes& fluvar)
{
    // Search for overlap with particles

    // Check boundary blocks first
    for (size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        boundvar.block[blockID]->intersect.assign(boundvar.block[blockID]->npts, 0);
    }

    for (size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        Vec_Tree tree(SIMDIM, boundvar.block[blockID]->coords, 50);
        tree.index->buildIndex();
        double searchDist = 0.9 * boundvar.block[blockID]->dx;
        searchDist = searchDist * searchDist;

        // Search for overlap with boundary particles (do the upper diagonal only)
        for (size_t ii = blockID; ii < boundvar.nblocks; ++ii)
        {
            /*if (ii == blockID)
                continue;*/

            for (size_t jj = 0; jj < boundvar.block[ii]->coords.size(); ++jj)
            {
                if (boundvar.block[ii]->intersect[jj] == 0)
                {
                    std::vector<neighbour_index> matches =
                        radius_search(tree, boundvar.block[ii]->coords[jj], searchDist);

                    // std::printf(ll << "  " << blockID << "  " << matches.size() <<
                    // std::endl;
                    for (auto const& match : matches)
                    {
                        if (match.first != jj)
                            boundvar.block[blockID]->intersect[match.first] = 1;
                    }
                }
            }
        }
    }

    // Now check fluid blocks
    for (size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        fluvar.block[blockID]->intersect.assign(fluvar.block[blockID]->npts, 0);
    }

    for (size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        Vec_Tree tree(SIMDIM, fluvar.block[blockID]->coords, 10);
        tree.index->buildIndex();
        double searchDist = 0.9 * fluvar.block[blockID]->dx;
        searchDist = searchDist * searchDist;

        // Search for overlap with boundary particles
        for (size_t ii = 0; ii < boundvar.nblocks; ++ii)
        {
            for (size_t jj = 0; jj < boundvar.block[ii]->coords.size(); ++jj)
            {
                if (boundvar.block[ii]->intersect[jj] == 0)
                {
                    std::vector<neighbour_index> matches =
                        radius_search(tree, boundvar.block[ii]->coords[jj], searchDist);

                    for (auto const& match : matches)
                        fluvar.block[blockID]->intersect[match.first] = 1;
                }
            }
        }

        // Check for intersection with other fluid blocks (do the upper diagonal
        // only)
        for (size_t ii = blockID; ii < fluvar.nblocks; ++ii)
        {
            if (ii == blockID)
                continue; /* Ignore points within its own block. Should not self
                             intersect */

            for (size_t jj = 0; jj < fluvar.block[ii]->npts; ++jj)
            {
                if (fluvar.block[ii]->intersect[jj] == 0)
                {
                    std::vector<neighbour_index> matches =
                        radius_search(tree, fluvar.block[ii]->coords[jj], searchDist);

                    // std::printf(ii << "  " << blockID << "  " << matches.size() <<
                    // std::endl;
                    for (size_t kk = 0; kk < matches.size(); kk++)
                    {
                        fluvar.block[blockID]->intersect[matches[kk].first] = 1;
                    }
                }
            }
        }
    }

    for (size_t ii = 0; ii < fluvar.nblocks; ++ii)
    {
        for (size_t bID = 0; bID < fluvar.block[ii]->back.size(); bID++)
        {
            size_t backID = fluvar.block[ii]->back[bID];
            int does_int = 0;
            if (fluvar.block[ii]->intersect[backID])
            { // If the back particle interstects, don't
              // place any of the inlet
                does_int = 1;
            }
            else
            { // Check buffer of that back particle now
                for (size_t const buffID : fluvar.block[ii]->buffer[bID])
                {
                    if (fluvar.block[ii]->intersect[buffID])
                        does_int = 1;
                }
            }

            if (does_int)
            { // Set the intersection of all buffer and back particles
              // to have correct numbers
                fluvar.block[ii]->intersect[backID] = 1;

                for (size_t const buffID : fluvar.block[ii]->buffer[bID])
                    fluvar.block[ii]->intersect[buffID] = 1;
            }
        }
    }

    size_t nBoundPts = 0, nFluidPts = 0;
    for (size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        // Update the points to be accurate to that retained.
        boundvar.block[blockID]->npts = std::count(
            boundvar.block[blockID]->intersect.begin(), boundvar.block[blockID]->intersect.end(), 0
        );

        nBoundPts += boundvar.block[blockID]->npts;
    }

    for (size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        // Update the points to be accurate to that retained.
        fluvar.block[blockID]->npts = std::count(
            fluvar.block[blockID]->intersect.begin(), fluvar.block[blockID]->intersect.end(), 0
        );

        nFluidPts += fluvar.block[blockID]->npts;
    }

    // Get number of fluid particles to keep
    printf(
        "Keeping %zu out of %zu boundary points. %.2f%% intersected.\n", nBoundPts,
        boundvar.total_points, 100.0 * (boundvar.total_points - nBoundPts) / boundvar.total_points
    );

    printf(
        "Keeping %zu out of %zu fluid points. %.2f%% intersected.\n", nFluidPts, fluvar.total_points,
        100.0 * (fluvar.total_points - nFluidPts) / fluvar.total_points
    );

    printf(
        "Keeping %zu out of %zu total particles.\n", nBoundPts + nFluidPts,
        fluvar.total_points + boundvar.total_points
    );
    boundvar.total_points = nBoundPts;
    fluvar.total_points = nFluidPts;
}

size_t Generate_Points(SIM const& svar, Shapes& var)
{
    size_t total_points = 0;
    size_t diff = 0;
    for (ShapeBlock* bound : var.block)
    {
        printf("Creating boundary block: %s\t...\t", bound->name.c_str());
        bound->generate_points();

        printf("npts: %zu\n", bound->coords.size());

        if (bound->coords.size() != bound->npts)
        {
            printf(
                "Number of boundary points generated for block \"%s\" differs "
                "from expected amount "
                "by %d\n",
                bound->name.c_str(),
                static_cast<int>(bound->npts) - static_cast<int>(bound->coords.size())
            );
            diff += bound->coords.size() - bound->npts;
        }

        bound->npts = bound->coords.size();
        total_points += bound->npts;

        if (svar.use_global_gas_law)
        { // Set block gas law properties as the global
          // properties
            bound->rho_rest = svar.fluid.rho_rest;
            bound->gamma = svar.fluid.gam;
            bound->speedOfSound = svar.fluid.speed_sound;
            bound->backgroundP = svar.fluid.press_pipe;
        }
    }
    var.total_points = total_points;
    // var.total_points = boundary_intersection(globalspacing, var);

    // get_boundary_velocity(boundvar);

    return diff;
}

void Init_Particles(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits)
{
    // Read boundary blocks
    Shapes boundvar;
    printf("Reading boundary settings...\n");
    if (iequals(std::filesystem::path(svar.io.input_bound_file).extension(), ".json"))
        boundvar = read_shapes_JSON(svar.io.input_bound_file, svar);
    else
        boundvar = read_shapes_bmap(svar.io.input_bound_file, svar);

    // Read fluid
    Shapes fluvar;
    printf("Reading fluid settings...\n");
    if (iequals(std::filesystem::path(svar.io.input_fluid_file).extension(), ".json"))
        fluvar = read_shapes_JSON(svar.io.input_fluid_file, svar);
    else
        fluvar = read_shapes_bmap(svar.io.input_fluid_file, svar);

    // Now generate points and add to indexes
    Generate_Points(svar, boundvar);
    Generate_Points(svar, fluvar);

    // Now check for intersections
    Check_Intersection(svar, boundvar, fluvar);

    // Insert points into the state vectors. Change to SoA?
    size_t part_id = 0;
    for (size_t block = 0; block < boundvar.nblocks; block++)
    {
        limits.emplace_back(part_id, 0);

        for (size_t ii = 0; ii < boundvar.block[block]->coords.size(); ii++)
        {
            if (!boundvar.block[block]->intersect[ii])
            {
                pn.emplace_back(SPHPart(
                    boundvar.block[block]->coords[ii], boundvar.block[block]->vel,
                    boundvar.block[block]->dens, boundvar.block[block]->mass,
                    boundvar.block[block]->press, BOUND, part_id
                ));
                part_id++;
            }
        }
        if (!boundvar.block[block]->times.empty())
        {
            if (!boundvar.block[block]->pos.empty() && boundvar.block[block]->vels.empty())
                get_boundary_velocity(boundvar.block[block]);
            else if (boundvar.block[block]->vels.empty())
            {
                printf(
                    "No velocity or position data available for boundary block %zu "
                    "even though times "
                    "were defined.\n",
                    block
                );
                exit(-1);
            }

            limits.back().times = boundvar.block[block]->times;
            limits.back().vels = boundvar.block[block]->vels;
            limits.back().nTimes = boundvar.block[block]->ntimes;
        }
        else
        {
            limits.back().nTimes = 0;
            limits.back().vels.emplace_back(boundvar.block[block]->vel);
        }

        limits.back().index.second = part_id;

        limits.back().name = boundvar.block[block]->name;
        limits.back().fixed_vel_or_dynamic = boundvar.block[block]->fixed_vel_or_dynamic;
        limits.back().particle_order = boundvar.block[block]->particle_order;
        limits.back().no_slip = boundvar.block[block]->no_slip;
        limits.back().bound_solver = boundvar.block[block]->bound_solver;
        limits.back().block_type = boundvar.block[block]->bound_type;
        limits.back().dx = boundvar.block[block]->dx;
    }

    svar.n_bound_blocks = boundvar.nblocks;
    svar.bound_points = part_id;

    // Fluid points
    for (size_t block = 0; block < fluvar.nblocks; block++)
    {
        limits.emplace_back(part_id, 0);
        if (fluvar.block[block]->bound_type == inletZone)
        {
            size_t ii = 0;
            size_t nBuff = fluvar.block[block]->particle_order == 1 ? 5 : 4;

            while (fluvar.block[block]->bc[ii] == PIPE)
            {
                if (!fluvar.block[block]->intersect[ii])
                {
                    pn.emplace_back(SPHPart(
                        fluvar.block[block]->coords[ii], fluvar.block[block]->vel,
                        fluvar.block[block]->dens, fluvar.block[block]->mass, fluvar.block[block]->press,
                        PIPE, part_id
                    ));
                    part_id++;
                }
                ii++;
            }

            // Check if the back line or buffer ends up intersecting
            vector<size_t> back_intersect(fluvar.block[block]->back.size(), 0);

            for (size_t bID = 0; bID < fluvar.block[block]->back.size(); bID++)
            {
                size_t does_int = 0;
                size_t backID = fluvar.block[block]->back[bID];

                // If the back particle interstects, don't place any of the inlet
                if (fluvar.block[block]->intersect[backID])
                    does_int = 1;

                back_intersect[bID] = does_int;
            }

            // Insert the back particles
            for (size_t bID = 0; bID < fluvar.block[block]->back.size(); bID++)
            {
                if (!back_intersect[bID])
                {
                    size_t pointID = fluvar.block[block]->back[bID];
                    pn.emplace_back(SPHPart(
                        fluvar.block[block]->coords[pointID], fluvar.block[block]->vel,
                        fluvar.block[block]->dens, fluvar.block[block]->mass, fluvar.block[block]->press,
                        BACK, part_id
                    ));

                    limits.back().back.emplace_back(part_id);
                    limits.back().buffer.emplace_back(nBuff, 0);

                    part_id++;
                }
            }
            // Insert the buffer particles
            for (size_t buffID = 0; buffID < nBuff; buffID++)
            {
                size_t bIndex = 0;
                for (size_t bID = 0; bID < fluvar.block[block]->back.size(); bID++)
                {
                    if (!back_intersect[bID])
                    {
                        size_t pointID = fluvar.block[block]->buffer[bID][buffID];
                        pn.emplace_back(SPHPart(
                            fluvar.block[block]->coords[pointID], fluvar.block[block]->vel,
                            fluvar.block[block]->dens, fluvar.block[block]->mass,
                            fluvar.block[block]->press, BUFFER, part_id
                        ));

                        limits.back().buffer[bIndex][buffID] = part_id;

                        part_id++;
                        bIndex++;
                    }
                }
            }
        }
        else
        {
            for (size_t ii = 0; ii < fluvar.block[block]->coords.size(); ii++)
            {
                if (!fluvar.block[block]->intersect[ii])
                {
                    pn.emplace_back(SPHPart(
                        fluvar.block[block]->coords[ii], fluvar.block[block]->vel,
                        fluvar.block[block]->dens, fluvar.block[block]->mass, fluvar.block[block]->press,
                        FREE, part_id
                    ));
                    part_id++;
                }
            }
        }

        limits.back().index.second = part_id;
        if (!fluvar.block[block]->times.empty())
        {
            limits.back().times = fluvar.block[block]->times;
            limits.back().vels = fluvar.block[block]->vels;
            limits.back().nTimes = fluvar.block[block]->ntimes;
        }
        else
        {
            limits.back().nTimes = 0;
            limits.back().vels.emplace_back(StateVecD::Zero());
        }

        limits.back().name = fluvar.block[block]->name;
        limits.back().insert_norm = fluvar.block[block]->insert_norm;
        limits.back().delete_norm = fluvar.block[block]->delete_norm;
        limits.back().aero_norm = fluvar.block[block]->aero_norm;
        limits.back().insconst = fluvar.block[block]->insconst;
        limits.back().delconst = fluvar.block[block]->delconst;
        limits.back().aeroconst = fluvar.block[block]->aeroconst;

        limits.back().fixed_vel_or_dynamic = fluvar.block[block]->fixed_vel_or_dynamic;
        limits.back().particle_order = fluvar.block[block]->particle_order;
        limits.back().no_slip = fluvar.block[block]->no_slip;
        limits.back().bound_solver = fluvar.block[block]->bound_solver;
        limits.back().block_type = fluvar.block[block]->bound_type;
        limits.back().dx = fluvar.block[block]->dx;
    }

    svar.n_fluid_blocks = fluvar.nblocks;
    svar.total_points = part_id;
    svar.fluid_points = svar.total_points - svar.bound_points;
    svar.part_id = part_id;

    if (svar.init_hydro_pressure)
    { /* If using hydrostatic initialisation, set
         pressure */
        printf("Initialising hydrostatic pressure...\n");
#pragma omp parallel for
        for (size_t ii = 0; ii < svar.total_points; ++ii)
        {
            real press =
                std::max(0.0, -svar.fluid.rho_rest * svar.grav[1] * (svar.hydro_height - pn[ii].xi[1]));
            real dens = svar.fluid.get_density(press);
            pn[ii].p = press;
            pn[ii].rho = dens;
        }
    }

    pnp1 = pn;
}

/* Initialise the limits data for a restart, basically everything except actual
 * points */
void Init_Particles_Restart(SIM& svar, LIMITS& limits)
{
    // Read boundary blocks
    Shapes boundvar;

    // Get the extension of the files. If it's JSON, then use the new functions.
    if (iequals(std::filesystem::path(svar.io.input_bound_file).extension(), ".json"))
        boundvar = read_shapes_JSON(svar.io.input_bound_file, svar);
    else
        boundvar = read_shapes_bmap(svar.io.input_bound_file, svar);

    // Read fluid
    Shapes fluvar;
    printf("Reading fluid settings...\n");
    if (iequals(std::filesystem::path(svar.io.input_fluid_file).extension(), ".json"))
        fluvar = read_shapes_JSON(svar.io.input_fluid_file, svar);
    else
        fluvar = read_shapes_bmap(svar.io.input_fluid_file, svar);

    limits.resize(boundvar.nblocks + fluvar.nblocks, bound_block(0));

    for (size_t block = 0; block < boundvar.nblocks; block++)
    {
        if (!boundvar.block[block]->times.empty())
        {
            if (!boundvar.block[block]->pos.empty() && boundvar.block[block]->vels.empty())
                get_boundary_velocity(boundvar.block[block]);
            else if (boundvar.block[block]->vels.empty())
            {
                printf(
                    "No velocity or position data available for boundary block %zu "
                    "even though times "
                    "were defined.\n",
                    block
                );
                exit(-1);
            }

            limits[block].times = boundvar.block[block]->times;
            limits[block].vels = boundvar.block[block]->vels;
            limits[block].nTimes = boundvar.block[block]->ntimes;
        }
        else
        {
            limits[block].nTimes = 0;
            limits[block].vels.emplace_back(boundvar.block[block]->vel);
        }

        limits[block].name = boundvar.block[block]->name;
        limits[block].fixed_vel_or_dynamic = boundvar.block[block]->fixed_vel_or_dynamic;
        limits[block].particle_order = boundvar.block[block]->particle_order;
        limits[block].no_slip = boundvar.block[block]->no_slip;
        limits[block].bound_solver = boundvar.block[block]->bound_solver;
        limits[block].block_type = boundvar.block[block]->bound_type;
    }

    size_t offset = boundvar.nblocks;
    for (size_t block = 0; block < fluvar.nblocks; block++)
    {
        if (!fluvar.block[block]->times.empty())
        {
            limits[block + offset].times = fluvar.block[block]->times;
            limits[block + offset].vels = fluvar.block[block]->vels;
            limits[block + offset].nTimes = fluvar.block[block]->ntimes;
        }
        else
        {
            limits[block + offset].nTimes = 0;
            limits[block + offset].vels.emplace_back(StateVecD::Zero());
        }

        limits[block + offset].name = fluvar.block[block]->name;
        limits[block + offset].insert_norm = fluvar.block[block]->insert_norm;
        limits[block + offset].delete_norm = fluvar.block[block]->delete_norm;
        limits[block + offset].aero_norm = fluvar.block[block]->aero_norm;
        limits[block + offset].insconst = fluvar.block[block]->insconst;
        limits[block + offset].delconst = fluvar.block[block]->delconst;
        limits[block + offset].aeroconst = fluvar.block[block]->aeroconst;

        limits[block + offset].fixed_vel_or_dynamic = fluvar.block[block]->fixed_vel_or_dynamic;
        limits[block + offset].particle_order = fluvar.block[block]->particle_order;
        limits[block + offset].no_slip = fluvar.block[block]->no_slip;
        limits[block + offset].bound_solver = fluvar.block[block]->bound_solver;
        limits[block + offset].block_type = fluvar.block[block]->bound_type;
    }

    svar.n_bound_blocks = boundvar.nblocks;
    svar.n_fluid_blocks = fluvar.nblocks;
}

void Init_Surface(SIM const& svar, MESH const& cells, vector<SURF>& surf_marks)
{
    surf_marks = vector<SURF>(svar.io.tau_markers.size());

    vector<vector<size_t>> faceIDs(svar.io.tau_markers.size());
    vector<vector<int>> tau_markers(svar.io.tau_markers.size());
    /* for each surface, find how many faces are in it. */
    for (std::pair<size_t, int> const& marker : cells.smarkers)
    {
        auto index = find(svar.io.tau_markers.begin(), svar.io.tau_markers.end(), marker.second);
        if (index != svar.io.tau_markers.end())
        {
            size_t mark = index - svar.io.tau_markers.begin();
            faceIDs[mark].emplace_back(marker.first);
            tau_markers[mark].emplace_back(marker.second);

            // surf_faces[mark].back().faceID  = marker.first;
            // if()
            // printf(mark << "  " << tau_markers.first << endl;
            // printf(surf_faces[mark].back().faceID << "  " << tau_markers.first << endl;
            // surf_faces[mark].back().marker  = marker.second;
        }
        else
        {
            printf("Couldn't find the marker in the index\n");
        }
    }

    for (size_t ii = 0; ii < svar.io.tau_markers.size(); ii++)
    {
        surf_marks[ii].name = svar.io.tau_bnames[ii];
        surf_marks[ii].marker = svar.io.tau_markers[ii];
        surf_marks[ii].output = svar.io.tau_bwrite[ii];

        size_t nFaces = faceIDs[ii].size();
        surf_marks[ii].face_IDs = faceIDs[ii];
        surf_marks[ii].face_count = vector<uint>(nFaces, 0);
        surf_marks[ii].face_beta = vector<real>(nFaces, 0.0);
        surf_marks[ii].face_area = vector<real>(nFaces, 0.0);
    }
}
