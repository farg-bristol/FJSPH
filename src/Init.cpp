/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Init.h"
#include "Add.h"
#include "IOFunctions.h"
#include "Kernel.h"
#include "Neighbours.h"

#include "shapes/arc.h"
#include "shapes/circle.h"
#include "shapes/cylinder.h"
#include "shapes/inlet.h"
#include "shapes/line.h"
#include "shapes/square.h"

#include <filesystem>

void get_boundary_velocity(shape_block& boundvar)
{
    size_t ntimes = boundvar.ntimes;
    if (ntimes != 0)
    {
        boundvar.vels.resize(ntimes - 1);
        for (size_t jj = 0; jj < ntimes - 1; jj++)
        {
            boundvar.vels[jj] = (boundvar.pos[jj + 1] - boundvar.pos[jj]) /
                                (boundvar.times[jj + 1] - boundvar.times[jj]);
        }
    }
}

void Read_Shapes(
    Shapes& var, real& globalspacing, SIM const& svar, FLUID const& fvar, std::string const& filename
)
{
    /* Read shapes from file */
    std::ifstream inFile(filename);
    if (!inFile.is_open())
    {
        std::cerr << filename << " file missing\n";
        exit(1);
    }

    std::string line;
    size_t nblocks = 0;

    /* Find number of blocks to initialise */
    while (getline(inFile, line))
    {
        size_t end = line.find_first_of('#');
        if (end != std::string::npos)
            line = line.substr(0, end + 1);

        if (line.find("block end") != std::string::npos)
        {
            nblocks++;
        }
    }

    inFile.clear();
    inFile.seekg(0);

    std::vector<shape_block> shapes(nblocks);

    size_t block = 0;
    while (getline(inFile, line))
    {
        size_t end = line.find_first_of('#');
        if (end != std::string::npos)
            line = line.substr(0, end);

        Get_String(line, "Name", shapes[block].name);
        Get_String(line, "Shape", shapes[block].shape);
        Get_String(line, "Sub-shape", shapes[block].subshape);
        Get_String(line, "Boundary solver", shapes[block].solver_name);

        Get_Bool(line, "Write surface data (0/1)", shapes[block].write_data);
        Get_Number(line, "Fixed velocity or dynamic inlet BC (0/1)", shapes[block].fixed_vel_or_dynamic);

        // Pipe exit plane.
        Get_Vector(line, "Aerodynamic entry normal", shapes[block].aero_norm);
        Get_Vector(line, "Deletion normal", shapes[block].delete_norm);
        Get_Vector(line, "Insertion normal", shapes[block].insert_norm);
        Get_Number(line, "Aerodynamic entry plane constant", shapes[block].aeroconst);
        Get_Number(line, "Deletion plane constant", shapes[block].delconst);
        Get_Number(line, "Insertion plane constant", shapes[block].insconst);
        Get_Number(line, "Pipe depth", shapes[block].thickness);

        Get_Number(line, "i-direction count", shapes[block].ni);
        Get_Number(line, "j-direction count", shapes[block].nj);
        Get_Number(line, "k-direction count", shapes[block].nk);

        Get_Vector(line, "Stretching factor", shapes[block].stretch);

        // Rotation definitions
        Get_Vector(line, "Normal vector", shapes[block].normal);
        Get_Vector(line, "Rotation angles", shapes[block].angles);
        Get_Number(line, "Rotation angle", shapes[block].angles[0]);

        // Square block definitions
        Get_Vector(line, "Start coordinate", shapes[block].start);
        Get_Vector(line, "End coordinate", shapes[block].end);
        Get_Vector(line, "Right coordinate", shapes[block].right);

        // Circle/arc definitions
        Get_Vector(line, "Midpoint coordinate", shapes[block].mid);
        Get_Vector(line, "Centre coordinate", shapes[block].centre);
        Get_Vector(line, "Arch normal", shapes[block].right);
        Get_Number(line, "Radius", shapes[block].radius);
        Get_Number(line, "Length", shapes[block].length);
        Get_Number(line, "Arc start (degree)", shapes[block].arc_start);
        Get_Number(line, "Arc end (degree)", shapes[block].arc_end);
        Get_Number(line, "Arc length (degree)", shapes[block].arclength);
        Get_Number(line, "Starting straight length", shapes[block].sstraight);
        Get_Number(line, "Ending straight length", shapes[block].estraight);

        Get_Number(line, "Particle spacing", shapes[block].dx);
        Get_Bool(line, "Particle ordering (0=grid,1=HCP)", shapes[block].particle_order);
        Get_Number(line, "Wall thickness", shapes[block].thickness);
        Get_Number(line, "Wall radial particle count", shapes[block].nk);
        Get_Bool(line, "Wall is no slip (0/1)", shapes[block].no_slip);

        Get_String(line, "Coordinate filename", shapes[block].filename);

        if (line.find("Coordinate data:") != std::string::npos)
        {
            std::string tmp;
            getline(inFile, tmp); /* Go to the next line and get data */
            size_t npts;
            inFile >> npts;
            shapes[block].npts = npts;
            shapes[block].coords = std::vector<StateVecD>(npts);

            for (size_t ii = 0; ii < shapes[block].npts; ++ii)
            {
                for (size_t dim = 0; dim < SIMDIM; ++dim)
                {
                    inFile >> shapes[block].coords[ii][dim];
                }
            }
        }

        Get_Vector(line, "Starting velocity", shapes[block].vel);
        Get_Number(line, "Starting jet velocity", shapes[block].vmag);
        Get_Number(line, "Starting pressure", shapes[block].press);
        Get_Number(line, "Starting density", shapes[block].dens);
        // Get_Number(line, "Starting mass", shapes[block].mass);

        Get_Number(line, "Cole EOS gamma", shapes[block].gamma);
        Get_Number(line, "Speed of sound", shapes[block].speedOfSound);
        Get_Number(line, "Resting density", shapes[block].rho0);
        Get_Number(line, "Volume to target", shapes[block].renorm_vol);

        Get_String(line, "Time data filename", shapes[block].position_filename);

        if (line.find("Time position data") != std::string::npos)
        {
            std::string tmp;
            getline(inFile, tmp); /* Go to the next line and get data */
            size_t ntimes;
            std::istringstream iss(tmp);
            iss >> ntimes;
            // inFile >> ntimes;
            shapes[block].ntimes = ntimes;
            shapes[block].times = std::vector<real>(ntimes);
            shapes[block].pos = std::vector<StateVecD>(ntimes);

            for (size_t ii = 0; ii < ntimes; ++ii)
            {
                getline(inFile, tmp); /* Go to the next line and get data */
                std::istringstream iss2(tmp);
                iss2 >> shapes[block].times[ii];
                for (size_t dim = 0; dim < SIMDIM; ++dim)
                {
                    iss2 >> shapes[block].pos[ii][dim];
                }
            }
        }

        if (line.find("Time velocity data") != std::string::npos)
        {
            std::string tmp;
            getline(inFile, tmp); /* Go to the next line and get data */
            size_t ntimes;
            std::istringstream iss(tmp);
            iss >> ntimes;
            // inFile >> ntimes;
            shapes[block].ntimes = ntimes;
            shapes[block].times = std::vector<real>(ntimes);
            shapes[block].vels = std::vector<StateVecD>(ntimes);

            for (size_t ii = 0; ii < ntimes; ++ii)
            {
                getline(inFile, tmp); /* Go to the next line and get data */
                std::istringstream iss2(tmp);
                iss2 >> shapes[block].times[ii];
                for (size_t dim = 0; dim < SIMDIM; ++dim)
                {
                    iss2 >> shapes[block].vels[ii][dim];
                }
            }
        }

        if (line.find("block end") != std::string::npos)
        {
            block++;
            if (block == nblocks)
                break;
        }
    }

    inFile.close();

    /* Check enough information has been provided */
    int fault = 0;
    for (shape_block& bound : shapes)
    {
        bound.check_input(svar, fvar, globalspacing, fault);
    }

    if (fault)
    {
        printf("Check of parameters for geometry blocks finished with errors.\n");
        printf("See output for details.\n");
        exit(-1);
    }

    var.nblocks = nblocks;
    var.block = shapes;

    // Estimate the number of particles in each Block
    for (size_t ii = 0; ii < var.nblocks; ++ii)
    {
        var.totPts += var.block[ii].npts;
    }

    return;
}

std::vector<StateVecD> Read_Geom_File(std::string const& filename)
{
    std::ifstream inFile(filename, std::ios::in);

    int nPts;
    inFile >> nPts;
    std::vector<StateVecD> points(nPts);

    for (int ii = 0; ii < nPts; ++ii)
    {
        for (int dim = 0; dim < SIMDIM; ++dim)
        {
            inFile >> points[ii][dim];
        }
    }

    inFile.close();

    return points;
}

void Check_Intersection(SIM const& svar, Shapes& boundvar, Shapes& fluvar)
{
    // Search for overlap with particles
    double searchDist = 0.9 * svar.dx;
    searchDist = searchDist * searchDist;

    // Check boundary blocks first
    for (size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        boundvar.block[blockID].intersect.assign(boundvar.block[blockID].npts, 0);
    }

    for (size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        Vec_Tree tree(SIMDIM, boundvar.block[blockID].coords, 50);
        tree.index->buildIndex();

        // Search for overlap with boundary particles (do the upper diagonal only)
        for (size_t ii = blockID; ii < boundvar.nblocks; ++ii)
        {
            /*if (ii == blockID)
                continue;*/

            for (size_t jj = 0; jj < boundvar.block[ii].coords.size(); ++jj)
            {
                if (boundvar.block[ii].intersect[jj] == 0)
                {
                    std::vector<neighbour_index> matches =
                        radius_search(tree, boundvar.block[ii].coords[jj], searchDist);

                    // std::printf(ll << "  " << blockID << "  " << matches.size() <<
                    // std::endl;
                    for (auto const& match : matches)
                    {
                        if (match.first != jj)
                            boundvar.block[blockID].intersect[match.first] = 1;
                    }
                }
            }
        }
    }

    // Now check fluid blocks
    for (size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        fluvar.block[blockID].intersect.assign(fluvar.block[blockID].npts, 0);
    }

    for (size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        Vec_Tree tree(SIMDIM, fluvar.block[blockID].coords, 10);
        tree.index->buildIndex();

        // Search for overlap with boundary particles
        for (size_t ii = 0; ii < boundvar.nblocks; ++ii)
        {
            for (size_t jj = 0; jj < boundvar.block[ii].coords.size(); ++jj)
            {
                if (boundvar.block[ii].intersect[jj] == 0)
                {
                    std::vector<neighbour_index> matches =
                        radius_search(tree, boundvar.block[ii].coords[jj], searchDist);

                    for (auto const& match : matches)
                        fluvar.block[blockID].intersect[match.first] = 1;
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

            for (size_t jj = 0; jj < fluvar.block[ii].npts; ++jj)
            {
                if (fluvar.block[ii].intersect[jj] == 0)
                {
                    std::vector<neighbour_index> matches =
                        radius_search(tree, fluvar.block[ii].coords[jj], searchDist);

                    // std::printf(ii << "  " << blockID << "  " << matches.size() <<
                    // std::endl;
                    for (size_t kk = 0; kk < matches.size(); kk++)
                    {
                        fluvar.block[blockID].intersect[matches[kk].first] = 1;
                    }
                }
            }
        }
    }

    for (size_t ii = 0; ii < fluvar.nblocks; ++ii)
    {
        for (size_t bID = 0; bID < fluvar.block[ii].back.size(); bID++)
        {
            size_t backID = fluvar.block[ii].back[bID];
            int does_int = 0;
            if (fluvar.block[ii].intersect[backID])
            { // If the back particle interstects, don't
              // place any of the inlet
                does_int = 1;
            }
            else
            { // Check buffer of that back particle now
                for (size_t const buffID : fluvar.block[ii].buffer[bID])
                {
                    if (fluvar.block[ii].intersect[buffID])
                        does_int = 1;
                }
            }

            if (does_int)
            { // Set the intersection of all buffer and back particles
              // to have correct numbers
                fluvar.block[ii].intersect[backID] = 1;

                for (size_t const buffID : fluvar.block[ii].buffer[bID])
                    fluvar.block[ii].intersect[buffID] = 1;
            }
        }
    }

    size_t nBoundPts = 0, nFluidPts = 0;
    for (size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        // Update the points to be accurate to that retained.
        boundvar.block[blockID].npts = std::count(
            boundvar.block[blockID].intersect.begin(), boundvar.block[blockID].intersect.end(), 0
        );

        nBoundPts += boundvar.block[blockID].npts;
    }

    for (size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        // Update the points to be accurate to that retained.
        fluvar.block[blockID].npts = std::count(
            fluvar.block[blockID].intersect.begin(), fluvar.block[blockID].intersect.end(), 0
        );

        nFluidPts += fluvar.block[blockID].npts;
    }

    // Get number of fluid particles to keep
    printf(
        "Keeping %zu out of %zu boundary points. %.2f%% intersected.\n", nBoundPts, boundvar.totPts,
        100.0 * (boundvar.totPts - nBoundPts) / boundvar.totPts
    );

    printf(
        "Keeping %zu out of %zu fluid points. %.2f%% intersected.\n", nFluidPts, fluvar.totPts,
        100.0 * (fluvar.totPts - nFluidPts) / fluvar.totPts
    );

    printf(
        "Keeping %zu out of %zu total particles.\n", nBoundPts + nFluidPts,
        fluvar.totPts + boundvar.totPts
    );
    boundvar.totPts = nBoundPts;
    fluvar.totPts = nFluidPts;
}

size_t Generate_Points(SIM const& svar, FLUID const& fvar, double const& globalspacing, Shapes& var)
{
    size_t totPts = 0;
    size_t diff = 0;
    for (shape_block& bound : var.block)
    {
        printf("Creating boundary block: %s\t...\t", bound.name.c_str());
        bound.generate_points(globalspacing);

        printf("npts: %zu\n", bound.coords.size());

        if (bound.coords.size() != bound.npts)
        {
            printf(
                "Number of boundary points generated for block \"%s\" differs "
                "from expected amount "
                "by %d\n",
                bound.name.c_str(), static_cast<int>(bound.npts) - static_cast<int>(bound.coords.size())
            );
            diff += bound.coords.size() - bound.npts;
        }

        bound.npts = bound.coords.size();
        totPts += bound.npts;

        if (svar.use_global_gas_law)
        { // Set block gas law properties as the global
          // properties
            bound.rho0 = fvar.rho0;
            bound.gamma = fvar.gam;
            bound.speedOfSound = fvar.Cs;
            bound.backgroundP = fvar.pPress;
        }
    }
    var.totPts = totPts;
    // var.totPts = boundary_intersection(globalspacing, var);

    // get_boundary_velocity(boundvar);

    return diff;
}

void Init_Particles(SIM& svar, FLUID& fvar, AERO& avar, SPHState& pn, SPHState& pnp1, LIMITS& limits)
{
    real dx = svar.dx;
    // Read boundary blocks
    Shapes boundvar;
    printf("Reading boundary settings...\n");

    // Get the extension of the files. If it's JSON, then use the new functions.
    if (std::filesystem::path(svar.boundfile).extension() == ".json" ||
        std::filesystem::path(svar.boundfile).extension() == ".JSON")
        boundvar = read_shapes_JSON(svar.boundfile, svar, fvar, dx);
    else
        Read_Shapes(boundvar, dx, svar, fvar, svar.boundfile);

    // Read fluid
    Shapes fluvar;
    printf("Reading fluid settings...\n");
    if (std::filesystem::path(svar.boundfile).extension() == ".json" ||
        std::filesystem::path(svar.boundfile).extension() == ".JSON")
        fluvar = read_shapes_JSON(svar.fluidfile, svar, fvar, dx);
    else
        Read_Shapes(fluvar, dx, svar, fvar, svar.fluidfile);

    // Now generate points and add to indexes
    Generate_Points(svar, fvar, svar.dx, boundvar);
    Generate_Points(svar, fvar, svar.dx, fluvar);

    // Now check for intersections
    Check_Intersection(svar, boundvar, fluvar);

    // Insert points into the state vectors. Change to SoA?
    size_t pID = 0;
    for (size_t block = 0; block < boundvar.nblocks; block++)
    {
        limits.emplace_back(pID, 0);

        for (size_t ii = 0; ii < boundvar.block[block].coords.size(); ii++)
        {
            if (!boundvar.block[block].intersect[ii])
            {
                pn.emplace_back(SPHPart(
                    boundvar.block[block].coords[ii], boundvar.block[block].vel,
                    boundvar.block[block].dens, fvar.bndM, boundvar.block[block].press, BOUND, pID
                ));
                pID++;
            }
        }
        if (!boundvar.block[block].times.empty())
        {
            if (!boundvar.block[block].pos.empty() && boundvar.block[block].vels.empty())
                get_boundary_velocity(boundvar.block[block]);
            else if (boundvar.block[block].vels.empty())
            {
                printf(
                    "No velocity or position data available for boundary block %zu "
                    "even though times "
                    "were defined.\n",
                    block
                );
                exit(-1);
            }

            limits.back().times = boundvar.block[block].times;
            limits.back().vels = boundvar.block[block].vels;
            limits.back().nTimes = boundvar.block[block].ntimes;
        }
        else
        {
            limits.back().nTimes = 0;
            limits.back().vels.emplace_back(boundvar.block[block].vel);
        }

        limits.back().index.second = pID;

        limits.back().name = boundvar.block[block].name;
        limits.back().fixed_vel_or_dynamic = boundvar.block[block].fixed_vel_or_dynamic;
        limits.back().particle_order = boundvar.block[block].particle_order;
        limits.back().no_slip = boundvar.block[block].no_slip;
        limits.back().bound_solver = boundvar.block[block].bound_solver;
        limits.back().block_type = boundvar.block[block].bound_type;
    }

    svar.nbound = boundvar.nblocks;
    svar.bndPts = pID;

    // Fluid points
    for (size_t block = 0; block < fluvar.nblocks; block++)
    {
        limits.emplace_back(pID, 0);
        if (fluvar.block[block].bound_type == inletZone)
        {
            size_t ii = 0;
            size_t nBuff = fluvar.block[block].particle_order == 1 ? 5 : 4;

            while (fluvar.block[block].bc[ii] == PIPE)
            {
                if (!fluvar.block[block].intersect[ii])
                {
                    pn.emplace_back(SPHPart(
                        fluvar.block[block].coords[ii], fluvar.block[block].vel,
                        fluvar.block[block].dens, fvar.bndM, fluvar.block[block].press, PIPE, pID
                    ));
                    pID++;
                }
                ii++;
            }

            // Check if the back line or buffer ends up intersecting
            vector<size_t> back_intersect(fluvar.block[block].back.size(), 0);

            for (size_t bID = 0; bID < fluvar.block[block].back.size(); bID++)
            {
                size_t does_int = 0;
                size_t backID = fluvar.block[block].back[bID];

                // If the back particle interstects, don't place any of the inlet
                if (fluvar.block[block].intersect[backID])
                    does_int = 1;

                back_intersect[bID] = does_int;
            }

            // Insert the back particles
            for (size_t bID = 0; bID < fluvar.block[block].back.size(); bID++)
            {
                if (!back_intersect[bID])
                {
                    size_t pointID = fluvar.block[block].back[bID];
                    pn.emplace_back(SPHPart(
                        fluvar.block[block].coords[pointID], fluvar.block[block].vel,
                        fluvar.block[block].dens, fvar.simM, fluvar.block[block].press, BACK, pID
                    ));

                    limits.back().back.emplace_back(pID);
                    limits.back().buffer.emplace_back(nBuff, 0);

                    pID++;
                }
            }
            // Insert the buffer particles
            for (size_t buffID = 0; buffID < nBuff; buffID++)
            {
                size_t bIndex = 0;
                for (size_t bID = 0; bID < fluvar.block[block].back.size(); bID++)
                {
                    if (!back_intersect[bID])
                    {
                        size_t pointID = fluvar.block[block].buffer[bID][buffID];
                        pn.emplace_back(SPHPart(
                            fluvar.block[block].coords[pointID], fluvar.block[block].vel,
                            fluvar.block[block].dens, fvar.simM, fluvar.block[block].press, BUFFER, pID
                        ));

                        limits.back().buffer[bIndex][buffID] = pID;

                        pID++;
                        bIndex++;
                    }
                }
            }
        }
        else
        {
            for (size_t ii = 0; ii < fluvar.block[block].coords.size(); ii++)
            {
                if (!fluvar.block[block].intersect[ii])
                {
                    pn.emplace_back(SPHPart(
                        fluvar.block[block].coords[ii], fluvar.block[block].vel,
                        fluvar.block[block].dens, fvar.simM, fluvar.block[block].press, FREE, pID
                    ));
                    pID++;
                }
            }
        }

        limits.back().index.second = pID;
        if (!fluvar.block[block].times.empty())
        {
            limits.back().times = fluvar.block[block].times;
            limits.back().vels = fluvar.block[block].vels;
            limits.back().nTimes = fluvar.block[block].ntimes;
        }
        else
        {
            limits.back().nTimes = 0;
            limits.back().vels.emplace_back(StateVecD::Zero());
        }

        limits.back().name = fluvar.block[block].name;
        limits.back().insert_norm = fluvar.block[block].insert_norm;
        limits.back().delete_norm = fluvar.block[block].delete_norm;
        limits.back().aero_norm = fluvar.block[block].aero_norm;
        limits.back().insconst = fluvar.block[block].insconst;
        limits.back().delconst = fluvar.block[block].delconst;
        limits.back().aeroconst = fluvar.block[block].aeroconst;

        limits.back().fixed_vel_or_dynamic = fluvar.block[block].fixed_vel_or_dynamic;
        limits.back().particle_order = fluvar.block[block].particle_order;
        limits.back().no_slip = fluvar.block[block].no_slip;
        limits.back().bound_solver = fluvar.block[block].bound_solver;
        limits.back().block_type = fluvar.block[block].bound_type;
    }

    svar.nfluid = fluvar.nblocks;
    svar.totPts = pID;
    svar.simPts = svar.totPts - svar.bndPts;
    svar.partID = pID;

    if (svar.init_hydro_pressure)
    { /* If using hydrostatic initialisation, set
         pressure */
        printf("Initialising hydrostatic pressure...\n");
#pragma omp parallel for
        for (size_t ii = 0; ii < svar.totPts; ++ii)
        {
            real press = std::max(0.0, -fvar.rho0 * svar.grav[1] * (svar.hydro_height - pn[ii].xi[1]));
            real dens = density_equation(press, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
            pn[ii].p = press;
            pn[ii].rho = dens;
        }
    }

    pnp1 = pn;
}

/* Initialise the limits data for a restart, basically everything except actual
 * points */
void Init_Particles_Restart(SIM& svar, FLUID& fvar, LIMITS& limits)
{
    real dx = svar.dx;
    // Read boundary blocks
    Shapes boundvar;
    printf("Reading boundary settings...\n");
    Read_Shapes(boundvar, dx, svar, fvar, svar.boundfile);

    // Read fluid
    Shapes fluvar;
    printf("Reading fluid settings...\n");
    Read_Shapes(fluvar, dx, svar, fvar, svar.fluidfile);

    limits.resize(boundvar.nblocks + fluvar.nblocks, bound_block(0));

    for (size_t block = 0; block < boundvar.nblocks; block++)
    {
        if (!boundvar.block[block].times.empty())
        {
            if (!boundvar.block[block].pos.empty() && boundvar.block[block].vels.empty())
                get_boundary_velocity(boundvar.block[block]);
            else if (boundvar.block[block].vels.empty())
            {
                printf(
                    "No velocity or position data available for boundary block %zu "
                    "even though times "
                    "were defined.\n",
                    block
                );
                exit(-1);
            }

            limits[block].times = boundvar.block[block].times;
            limits[block].vels = boundvar.block[block].vels;
            limits[block].nTimes = boundvar.block[block].ntimes;
        }
        else
        {
            limits[block].nTimes = 0;
            limits[block].vels.emplace_back(boundvar.block[block].vel);
        }

        limits[block].name = boundvar.block[block].name;
        limits[block].fixed_vel_or_dynamic = boundvar.block[block].fixed_vel_or_dynamic;
        limits[block].particle_order = boundvar.block[block].particle_order;
        limits[block].no_slip = boundvar.block[block].no_slip;
        limits[block].bound_solver = boundvar.block[block].bound_solver;
        limits[block].block_type = boundvar.block[block].bound_type;
    }

    size_t offset = boundvar.nblocks;
    for (size_t block = 0; block < fluvar.nblocks; block++)
    {
        if (!fluvar.block[block].times.empty())
        {
            limits[block + offset].times = fluvar.block[block].times;
            limits[block + offset].vels = fluvar.block[block].vels;
            limits[block + offset].nTimes = fluvar.block[block].ntimes;
        }
        else
        {
            limits[block + offset].nTimes = 0;
            limits[block + offset].vels.emplace_back(StateVecD::Zero());
        }

        limits[block + offset].name = fluvar.block[block].name;
        limits[block + offset].insert_norm = fluvar.block[block].insert_norm;
        limits[block + offset].delete_norm = fluvar.block[block].delete_norm;
        limits[block + offset].aero_norm = fluvar.block[block].aero_norm;
        limits[block + offset].insconst = fluvar.block[block].insconst;
        limits[block + offset].delconst = fluvar.block[block].delconst;
        limits[block + offset].aeroconst = fluvar.block[block].aeroconst;

        limits[block + offset].fixed_vel_or_dynamic = fluvar.block[block].fixed_vel_or_dynamic;
        limits[block + offset].particle_order = fluvar.block[block].particle_order;
        limits[block + offset].no_slip = fluvar.block[block].no_slip;
        limits[block + offset].bound_solver = fluvar.block[block].bound_solver;
        limits[block + offset].block_type = fluvar.block[block].bound_type;
    }

    svar.nbound = boundvar.nblocks;
    svar.nfluid = fluvar.nblocks;
}

void Init_Surface(SIM const& svar, MESH const& cells, vector<SURF>& surf_marks)
{
    surf_marks = vector<SURF>(svar.markers.size());

    vector<vector<size_t>> faceIDs(svar.markers.size());
    vector<vector<int>> markers(svar.markers.size());
    /* for each surface, find how many faces are in it. */
    for (std::pair<size_t, int> const& marker : cells.smarkers)
    {
        auto index = find(svar.markers.begin(), svar.markers.end(), marker.second);
        if (index != svar.markers.end())
        {
            size_t mark = index - svar.markers.begin();
            faceIDs[mark].emplace_back(marker.first);
            markers[mark].emplace_back(marker.second);

            // surf_faces[mark].back().faceID  = marker.first;
            // if()
            // printf(mark << "  " << markers.first << endl;
            // printf(surf_faces[mark].back().faceID << "  " << markers.first << endl;
            // surf_faces[mark].back().marker  = marker.second;
        }
        else
        {
            printf("Couldn't find the marker in the index\n");
        }
    }

    for (size_t ii = 0; ii < svar.markers.size(); ii++)
    {
        surf_marks[ii].name = svar.bnames[ii];
        surf_marks[ii].marker = svar.markers[ii];
        surf_marks[ii].output = svar.bwrite[ii];

        size_t nFaces = faceIDs[ii].size();
        surf_marks[ii].faceIDs = faceIDs[ii];
        surf_marks[ii].face_count = vector<uint>(nFaces, 0);
        surf_marks[ii].face_beta = vector<real>(nFaces, 0.0);
        surf_marks[ii].face_area = vector<real>(nFaces, 0.0);
    }
}
