/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "shapes.h"

#include "arc.h"
#include "circle.h"
#include "cylinder.h"
#include "inlet.h"
#include "line.h"
#include "square.h"

#include "../Third_Party/nlohmann/json.hpp"

using json = nlohmann::json;

template <typename T> void get_var(json const& file_data, string const& param, T& value)
{
    // Try-catch to give some more information about where the error occurred.
    try
    {
        // See if the parameter exists before trying to retrieve data.
        if (file_data.find(param) != file_data.end())
        {
            // parameter exists, now retrieve it in the format of the data it will be placed in.
            value = file_data[param].template get<T>();
        }
    }
    catch (const std::exception& e)
    {
        printf("An error occured trying to read a parameter from the JSON file.\n");
        printf("File: %s\n", file_data["filename"].get<string>().c_str());
        printf("Parameter: %s\n", param.c_str());
        printf("Error: %s\n", e.what());
        exit(-1);
    }
}

// Create a separate function for Eigen vectors, as it doesn't have the right setup.
void get_var(json const& file_data, string const& param, StateVecD& value)
{
    // Use a std vector to allow use of the base get_var
    std::vector<real> temp_vec;
    get_var(file_data, param, temp_vec);

    // If the vector has the right size, then use it to initialise the Eigen vector.
    // If it didn't find the parameter, then the vector will be size zero.
    if (temp_vec.size() == SIMDIM)
        value(temp_vec.data());
}

// Create a separate function for Eigen vectors, as it doesn't have the right setup.
void get_var(json const& file_data, string const& param, std::vector<StateVecD>& value)
{
    // Use a std vector to allow use of the base get_var
    std::vector<std::vector<real>> temp_vec;
    get_var(file_data, param, temp_vec);

    // If it didn't find the parameter, then the vector will be size zero.
    if (!temp_vec.empty())
    {
        size_t vec_length = temp_vec.size();
        value = std::vector<StateVecD>(vec_length);
        // Need to loop through to initialise each vector in the vector.
        for (size_t ii = 0; ii < vec_length; ii++)
            value[ii] = StateVecD(temp_vec[ii].data());
    }
}

void read_shape_JSON(
    json const& input_block, SIM const& svar, FLUID const& fvar, real& globalspacing,
    shape_block& output_block, int& fault
)
{
    get_var(input_block, "filename", output_block.filename);
    get_var(input_block, "Shape", output_block.shape);
    get_var(input_block, "Sub-shape", output_block.subshape);
    get_var(input_block, "Boundary solver", output_block.solver_name);
    get_var(input_block, "Aerodynamic entry normal", output_block.aero_norm);

    get_var(input_block, "Fixed velocity or dynamic inlet BC", output_block.inlet_bc_type);

    // Pipe exit plane.
    get_var(input_block, "Aerodynamic entry normal", output_block.aero_norm);
    get_var(input_block, "Deletion normal", output_block.delete_norm);
    get_var(input_block, "Insertion normal", output_block.insert_norm);
    get_var(input_block, "Aerodynamic entry plane constant", output_block.aeroconst);
    get_var(input_block, "Deletion plane constant", output_block.delconst);
    get_var(input_block, "Insertion plane constant", output_block.insconst);
    get_var(input_block, "Pipe depth", output_block.thickness);

    get_var(input_block, "i-direction count", output_block.ni);
    get_var(input_block, "j-direction count", output_block.nj);
    get_var(input_block, "k-direction count", output_block.nk);

    get_var(input_block, "Stretching factor", output_block.stretch);

    // Rotation definitions
    get_var(input_block, "Normal vector", output_block.normal);
    get_var(input_block, "Rotation angles (degree)", output_block.angles);
    get_var(input_block, "Rotation angle (degree)", output_block.angles[0]);

    // Square input_block definitions
    get_var(input_block, "Start coordinate", output_block.start);
    get_var(input_block, "End coordinate", output_block.end);
    get_var(input_block, "Right coordinate", output_block.right);

    // Circle/arc definitions
    get_var(input_block, "Midpoint coordinate", output_block.mid);
    get_var(input_block, "Centre coordinate", output_block.centre);
    get_var(input_block, "Arch normal", output_block.right);
    get_var(input_block, "Radius", output_block.radius);
    get_var(input_block, "Length", output_block.length);
    get_var(input_block, "Arc start (degree)", output_block.arc_start);
    get_var(input_block, "Arc end (degree)", output_block.arc_end);
    get_var(input_block, "Arc length (degree)", output_block.arclength);
    get_var(input_block, "Starting straight length", output_block.sstraight);
    get_var(input_block, "Ending straight length", output_block.estraight);

    get_var(input_block, "Particle spacing", output_block.dx);
    get_var(input_block, "Particle ordering (Grid/HCP)", output_block.particle_order);
    get_var(input_block, "Wall thickness", output_block.thickness);
    get_var(input_block, "Wall radial particle count", output_block.nk);
    get_var(input_block, "Wall is no-slip", output_block.no_slip);

    get_var(input_block, "Starting velocity", output_block.vel);
    get_var(input_block, "Starting jet velocity", output_block.vmag);
    get_var(input_block, "Starting pressure", output_block.press);
    get_var(input_block, "Starting density", output_block.dens);
    // get_var(input_block, "Starting mass", output_block.mass);

    get_var(input_block, "Cole EOS gamma", output_block.gamma);
    get_var(input_block, "Speed of sound", output_block.speedOfSound);
    get_var(input_block, "Resting density", output_block.rho0);
    get_var(input_block, "Volume to target", output_block.renorm_vol);

    get_var(input_block, "Coordinate filename", output_block.filename);

    get_var(input_block, "Coordinate data", output_block.coords);

    get_var(input_block, "Time data filename", output_block.position_filename);
    get_var(input_block, "Time data", output_block.times);
    get_var(input_block, "Position data", output_block.pos);

    output_block.check_input(svar, fvar, globalspacing, fault);
}

// Check the inputs for a block to make sure they are valid.
void shape_block::check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault)
{
    // Scale things first before doing checks and distances
    if (svar.scale != 1.0)
    {
        if (check_vector(start))
            start *= svar.scale;

        if (check_vector(end))
            end *= svar.scale;

        if (check_vector(right))
            right *= svar.scale;

        if (check_vector(mid))
            mid *= svar.scale;

        if (check_vector(centre))
            centre *= svar.scale;

        if (aeroconst != default_val)
            aeroconst *= svar.scale;

        if (insconst != default_val)
            insconst *= svar.scale;

        if (delconst != default_val)
            delconst *= svar.scale;
    }

#if SIMDIM == 2
    if (shape == "Line")
#else
    if (shape == "Plane")
#endif
    {
        bound_type = linePlane;
        LineShape::check_input(*this, globalspacing, fault);
    }
#if SIMDIM == 2
    else if (shape == "Square")
#else
    else if (shape == "Cube")
#endif
    {
        bound_type = squareCube;
        SquareShape::check_input(*this, globalspacing, fault);
    }
#if SIMDIM == 2
    else if (shape == "Circle")
#else
    else if (shape == "Sphere")
#endif
    {
        bound_type = circleSphere;
        CircleShape::check_input(*this, globalspacing, fault);
    }
#if SIMDIM == 2
    else if (shape == "Arc")
#else
    else if (shape == "Arch")
#endif
    {
        bound_type = arcSection;
        ArcShape::check_input(*this, globalspacing, fault);
    }
    else if (shape == "Cylinder")
    {
        bound_type = cylinder;
        CylinderShape::check_input(*this, globalspacing, fault);
    }
    else if (shape == "Inlet")
    {
        bound_type = inletZone;
        InletShape::check_input(*this, globalspacing, fault);
    }
    else if (shape == "Coordinates")
    {
        bound_type = coordDef;
        if (filename.empty())
        {
            if (coords.empty())
            {
                printf(
                    "ERROR: Block \"%s\" coordinates have not been ingested "
                    "properly. Stopping.\n",
                    name.c_str()
                );
                fault = 1;
            }
            else
            {
                npts = coords.size();
            }
        }
    }

    if (position_filename.empty())
    {
        if (ntimes != 0)
        {
            if (times.empty())
            {
                printf("ERROR: Block \"%s\" has no time information\n", name.c_str());
                fault = 1;
            }
            if (pos.empty() && vels.empty())
            {
                printf(
                    "ERROR: Block \"%s\" position or velocity data has not been "
                    "ingested properly.\n",
                    name.c_str()
                );
                fault = 1;
            }
        }
    }

    if (bound_type == -1)
    {
        printf("ERROR: Block \"%s\" shape has not been correctly defined.\n", name.c_str());
        printf("File: %s\n", filename.c_str());
        fault = 1;
    }

    if (!solver_name.empty())
    {
        if (solver_name == "DBC")
        {
            bound_solver = DBC;
        }
        else if (solver_name == "Pressure-Gradient")
        {
            bound_solver = pressure_G;
        }
        else if (solver_name == "Ghost")
        {
            bound_solver = ghost;
        }
        else
        {
            printf("WARNING: Unrecognised boundary solver defined. Continuing to "
                   "use the default, Pressure-Gradient\n");
            printf("         Choose from the following options for a correct "
                   "definition: \n");
            printf("         \t1. DBC\n         \t2. Pressure-Gradient\n         "
                   "\t3. Ghost\n");
        }
    }

    if (!inlet_bc_type.empty())
    {
        if (inlet_bc_type == "Fixed Velocity")
        {
            fixed_vel_or_dynamic = fixedVel;
        }
        else if (inlet_bc_type == "Dynamic")
        {
            fixed_vel_or_dynamic = dynamicVel;
        }
        else
        {
            printf("WARNING: Unrecognised inlet velocity type defined. Continuing to "
                   "use the default, Fixed Velocity\n");
            printf("         Choose from the following options for a correct "
                   "definition: \n");
            printf("         \t1. Fixed Velocity\n         \t2. Dynamic\n");
        }
    }

    if (!particle_order_str.empty())
    {
        if (particle_order_str == "Grid")
        {
            particle_order = grid;
        }
        else if (particle_order_str == "HCP")
        {
            particle_order = hcp;
        }
        else
        {
            printf("WARNING: Unrecognised particle packing type defined. Continuing to "
                   "use the default, Grid\n");
            printf("         Choose from the following options for a correct "
                   "definition: \n");
            printf("         \t1. Grid\n         \t2. HCP\n");
        }
    }

    globalspacing = std::max(dx, globalspacing);
    npts = npts > 1 ? npts : 1;
    if (bound_type != coordDef)
        coords.reserve(npts);

    if (press != 0)
    {
        real Bconst = fvar.rho0 * (fvar.Cs * fvar.Cs) / fvar.gam;
        dens = pow(((press - fvar.pPress) / Bconst + 1.0), fvar.gam) * fvar.rho0;
    }
    else
    {
        dens = fvar.rho0;
    }

    if (nu < 0)
    {
        nu = fvar.nu;
    }

    if (fault)
    {
        printf("Check of parameters for block \"%s\" finished with errors.\n\n", name.c_str());
    }
}

Shapes
read_shapes_JSON(std::string const& filename, SIM const& svar, FLUID const& fvar, real& globalspacing)
{

    /* Read shapes from file */
    std::ifstream inFile(filename);
    if (!inFile.is_open())
    {
        std::cerr << filename << " file missing\n";
        exit(-1);
    }

    json data = json::parse(inFile);

    // Store the filename for error reporting
    data["filename"] = filename;

    // Create an empty vector of shapes. This will be filled as the file is read.
    Shapes shapes;
    int fault = 0;
    for (auto& input_block : data.items())
    {
        shapes.emplace_back();
        shape_block& output_block = shapes.back();

        // JSON key is the name of the particle block.
        output_block.name = input_block.key();

        read_shape_JSON(input_block, svar, fvar, globalspacing, output_block, fault);
    }

    if (fault)
    {
        printf("Check of parameters for geometry blocks finished with errors.\n");
        printf("See output for details.\n");
        exit(-1);
    }

    return shapes;
}

void shape_block::generate_points(real const& globalspacing)
{
    switch (bound_type)
    {
    case linePlane:
        coords = LineShape::generate_points(*this, globalspacing);
        break;
    case squareCube:
        coords = SquareShape::generate_points(*this, globalspacing);
        break;
    case circleSphere:
        coords = CircleShape::generate_points(*this, globalspacing);
        break;
    case arcSection:
        coords = ArcShape::generate_points(*this, globalspacing);
        break;
    case cylinder:
        coords = CircleShape::generate_points(*this, globalspacing);
        break;
    case inletZone:
        coords = InletShape::generate_points(*this, globalspacing);
        break;
    case coordDef:
        break;
    default:
        printf("Unrecognised boundary type");
        break;
    }
}