/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "FileIO.h"
#include "BinaryIO.h"
#include "H5IO.h"
#include "Helper_Functions.h"
#include "IO.h"
#include "IOFunctions.h"
#include "VarDefs.h"

#include <hdf5.h>

FileReadWriter::FileReadWriter(file_extension const& _file_ext, SIM const& svar)
    : file_ext(_file_ext), scale(svar.scale), single_file(svar.io.single_file),
      var_names(svar.io.var_names), var_types(svar.io.var_types)
{
    // Get the correct file extension string
    switch (file_ext)
    {
    case dat:
        ext_str = ".dat";
        break;
    case h5:
        ext_str = ".h5";
        break;
    case szplt:
        ext_str = ".szplt";
        break;
    case h5part:
        ext_str = ".h5part";
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        exit(-1);
    }
}

FileReadWriter::~FileReadWriter()
{
    // Ensure file is closed
    if (get_file_is_open())
    {
        close_file();
    }
}

int FileReadWriter::open_file(std::string const& filename, std::string const& zone_name = "")
{
    /*
     * Create a new file. If file exists its contents will be overwritten.
     */
    std::filesystem::path _filename{filename};
    _filename.replace_extension(ext_str);
    // if (std::filesystem::exists(_filename))
    // { // Remove the file
    //     std::filesystem::remove(_filename);
    // }

    switch (file_ext)
    {
    case dat:
    {
        /* code */
        break;
    }
    case h5:
    case h5part:
    {
        h5_file_handle = HDF5::Open_File(_filename.string());
        break;
    }
    case szplt:
    {
        if (SZLIO::open_file(_filename.string(), zone_name, var_names, tec_file_handle))
        {
            cout << "Failed to open " << _filename.string() << endl;
            exit(-1);
        }
        break;
    }
    }
}

void FileReadWriter::close_file()
{
    if (!get_file_is_open())
        return; // Already closed

    switch (file_ext)
    {
    case dat:
    {
        /* code */
        break;
    }
    case h5:
    case h5part:
    {
        HDF5::Close_Handle(h5_file_handle);
        break;
    }
    case szplt:
    {
        SZLIO::close_file(tec_file_handle);
        break;
    }
    }
}

void FileReadWriter::write_file_attributes(SIM const& svar)
{
    // Write crucial simulation auxiliary data. Potential to expand
    write_file_attribute("Simulation current time", svar.integrator.current_time);
    write_file_attribute("Particle index to add", svar.part_id);
    write_file_attribute("Number of boundary blocks", svar.n_bound_blocks);
    write_file_attribute("Number of fluid blocks", svar.n_fluid_blocks);
    write_file_attribute("Number of boundary points", svar.bound_points);
    write_file_attribute("Number of fluid points", svar.fluid_points);
    write_file_attribute("Current frame", svar.integrator.current_frame);

    write_file_attribute("Input para filename", svar.io.input_file);
    write_file_attribute("Input fluid definition filename", svar.io.input_fluid_file);
    write_file_attribute("Input boundary definition filename", svar.io.input_bound_file);

    write_file_attribute("Primary grid face filename", svar.io.tau_mesh);
    write_file_attribute("Boundary mapping filename", svar.io.tau_bmap);
    write_file_attribute("Restart-data prefix", svar.io.tau_sol);
    write_file_attribute("Dimension offset vector", svar.io.offset_axis);
    write_file_attribute("Solution angle of attack", svar.io.angle_alpha);
    write_file_attribute("Grid scale", svar.scale);

    write_file_attribute("OpenFOAM input directory", svar.io.foam_dir);
    write_file_attribute("OpenFOAM solution directory", svar.io.foam_sol);
    write_file_attribute("OpenFOAM binary", svar.io.foam_is_binary);
    write_file_attribute("Label size", svar.io.foam_label_size);
    write_file_attribute("Scalar size", svar.io.foam_scalar_size);
    write_file_attribute("OpenFOAM is buoyant", svar.io.foam_buoyant_sim);
    write_file_attribute("OpenFOAM is foam_is_incompressible", svar.io.foam_is_incomp);

    write_file_attribute("VLM definition filename", svar.io.vlm_file);

    write_file_attribute("Single file for output", svar.io.single_file);
    write_file_attribute("Output files prefix", svar.io.output_prefix);
    write_file_attribute("SPH restart prefix", svar.io.restart_prefix);
    write_file_attribute("SPH frame time interval", svar.integrator.frame_time_interval);
    write_file_attribute("SPH frame count", svar.integrator.max_frames);
    write_file_attribute("SPH previous frame time", svar.integrator.last_frame_time);
    write_file_attribute("SPH output encoding", svar.io.out_encoding);
    write_file_attribute("Variable list", svar.io.output_names);

    /* Fluid data */
    write_file_attribute("Reference density", svar.air.rho_g);
    write_file_attribute("Init hydrostatic pressure", svar.init_hydro_pressure);
    write_file_attribute("Hydrostatic height", svar.hydro_height);

    /* Aerodynamic data */
    write_file_attribute("Reference velocity", svar.air.v_ref);
    write_file_attribute("Reference pressure", svar.air.p_ref);
    write_file_attribute("Reference Mach number", svar.air.M_ref);
    write_file_attribute("Reference temperature", svar.air.temp_g);
    write_file_attribute("Gas constant gamma", svar.air.gamma);

    /* Simulation settings */
    write_file_attribute("SPH integration solver", svar.integrator.solver_name);
    write_file_attribute("SPH boundary solver", svar.integrator.bound_solver);
    write_file_attribute("SPH solver minimum residual", svar.integrator.min_residual);
    write_file_attribute("SPH maximum timestep", svar.integrator.delta_t_max);
    write_file_attribute("SPH minimum timestep", svar.integrator.delta_t_min);
    write_file_attribute("SPH maximum CFL", svar.integrator.cfl_max);
    write_file_attribute("SPH minimum CFL", svar.integrator.cfl_min);
    write_file_attribute("SPH CFL condition", svar.integrator.cfl);
    write_file_attribute("SPH unstable CFL step", svar.integrator.cfl_step);
    write_file_attribute("SPH unstable CFL count limit", svar.integrator.n_unstable_limit);
    write_file_attribute("SPH stable CFL count limit", svar.integrator.n_stable_limit);
    write_file_attribute("SPH stable CFL count iteration factor", svar.integrator.subits_factor);
    write_file_attribute("SPH maximum shifting velocity", svar.integrator.max_shift_vel);
    write_file_attribute("SPH stable CFL count", svar.integrator.n_stable);
    write_file_attribute("SPH unstable CFL count", svar.integrator.n_unstable);

    write_file_attribute("SPH density variation", svar.rho_var);
    write_file_attribute("SPH maximum density", svar.rho_max);
    write_file_attribute("SPH minimum density", svar.rho_min);
    write_file_attribute("SPH Newmark Beta iteration limit", svar.integrator.max_subits);
    write_file_attribute("SPH gravity vector", svar.grav);

    write_file_attribute("SPH smoothing length factor", svar.H_fac);
    write_file_attribute("SPH aerodynamic case", svar.air.aero_case);
    write_file_attribute("SPH SP diameter definition", svar.air.use_dx);
    write_file_attribute("SPH use TAB deformation", svar.air.use_TAB_def);
    write_file_attribute("SPH global offset coordinate", svar.offset_vec);
    write_file_attribute("SPH maximum particle count", svar.max_points);
    write_file_attribute("SPH aerodynamic cutoff value", svar.air.lam_cutoff);
    write_file_attribute("SPH freestream velocity", svar.air.v_inf);
    write_file_attribute("SPH restart fit tolerance", svar.io.restart_tol);

    /* Particle tracking settings */
    write_file_attribute("Transition to IPT", svar.ipt.using_ipt);
    write_file_attribute("Velocity equation order", svar.ipt.ipt_eq_order);
    write_file_attribute("SPH tracking conversion x coordinate", svar.ipt.max_x_sph);
    write_file_attribute("Maximum x trajectory coordinate", svar.ipt.max_x);
    write_file_attribute("Particle scatter output", svar.ipt.part_out);
    write_file_attribute("Particle streak output", svar.ipt.streak_out);
    write_file_attribute("Particle cell intersection output", svar.ipt.cells_out);
}

void FileReadWriter::write_zone_attributes(bound_block const& limits)
{
    write_zone_attribute("Block name", limits.name);
    // Add information about the insertion, aero, and deletion planes.
    write_zone_attribute("Insertion normal x", limits.insert_norm[0]);
    write_zone_attribute("Insertion normal y", limits.insert_norm[1]);
    write_zone_attribute("Deletion normal x", limits.delete_norm[0]);
    write_zone_attribute("Deletion normal y", limits.delete_norm[1]);
    write_zone_attribute("Aerodynamic normal x", limits.aero_norm[0]);
    write_zone_attribute("Aerodynamic normal y", limits.aero_norm[1]);
#if SIMDIM == 3
    write_zone_attribute("Insertion normal z", limits.insert_norm[2]);
    write_zone_attribute("Deletion normal z", limits.delete_norm[2]);
    write_zone_attribute("Aerodynamic normal z", limits.aero_norm[2]);
#endif

    write_zone_attribute("Insertion plane constant", limits.insconst);
    write_zone_attribute("Deletion plane constant", limits.delconst);
    write_zone_attribute("Aerodynamic plane constant", limits.aeroconst);

    // Add boundary time information
    if (limits.nTimes != 0)
    {
        write_zone_attribute("Times count", limits.nTimes);
        for (size_t time = 0; time < limits.nTimes; time++)
        {
            write_zone_attribute("Time " + std::to_string(time), limits.times[time]);
            write_zone_attribute("Velocity " + std::to_string(time) + " x", limits.vels[time][0]);
            write_zone_attribute("Velocity " + std::to_string(time) + " y", limits.vels[time][1]);
#if SIMDIM == 3
            write_zone_attribute("Velocity " + std::to_string(time) + " z", limits.vels[time][2]);
#endif
        }
    }
    else
    {
        write_zone_attribute("Times count", 0);
        write_zone_attribute("Velocity x", limits.vels[0][0]);
        write_zone_attribute("Velocity y", limits.vels[0][1]);
#if SIMDIM == 3
        write_zone_attribute("Velocity z", limits.vels[0][2]);
#endif
    }

    write_zone_attribute("Fixed velocity or dynamic inlet", limits.fixed_vel_or_dynamic);
    write_zone_attribute("Lattice or HCPL packing", limits.particle_order);
    write_zone_attribute("Boundary solver type", limits.bound_solver);
    write_zone_attribute("Boundary is no slip", limits.no_slip);
    write_zone_attribute("Block type", limits.block_type);

    // Fluid properties
    write_zone_attribute("Cole EOS gamma", limits.fluid.gam);
    write_zone_attribute("Reference dispersed density", limits.fluid.rho_rest);
    write_zone_attribute("Reference dispersed viscosity", limits.fluid.mu);
    write_zone_attribute("Reference surface tension", limits.fluid.sig);
    write_zone_attribute("Surface tension contact angle", limits.fluid.contangb);
    write_zone_attribute("Equation of state", limits.fluid.pressure_rel);
    write_zone_attribute("Background pressure", limits.fluid.press_back);
    write_zone_attribute("Starting pressure", limits.fluid.press_start);
    write_zone_attribute("SPH delta coefficient", limits.fluid.dsph_delta);
    write_zone_attribute("Artificial viscosity factor", limits.fluid.visc_alpha);
    write_zone_attribute("Speed of sound", limits.fluid.speed_sound);
}

template <typename T>
void FileReadWriter::write_file_attribute(std::string const& attr_name, T const& attr_data)
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        HDF5::Write_Attribute(h5_file_handle, attr_name, attr_data);
        break;
    case szplt:
        SZLIO::add_file_aux_data(tec_file_handle, attr_name, attr_data);
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        break;
    }
}

template <typename T>
void FileReadWriter::write_zone_attribute(std::string const& attr_name, T const& attr_data)
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        HDF5::Write_Attribute(h5_zone_handle, attr_name, attr_data);
        break;
    case szplt:
        SZLIO::add_zone_aux_data(tec_file_handle, tec_output_zone, attr_name, attr_data);
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        break;
    }
}

template <typename T>
void FileReadWriter::write_scalar_var(
    std::string const& var_name, std::vector<T> const& var_data, int32_t& var_count
)
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        hsize_t dims[] = {var_data.size()};
        HDF5::Write_Variable(h5_zone_handle, var_name, dims, var_data);
        break;
    case szplt:
        SZLIO::Write_Vector(tec_file_handle, tec_output_zone, var_count, 0, imax, var_data, var_name);
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        break;
    }
}

void FileReadWriter::create_zone(
    std::string const& zone_header, real const& solution_time = 0, int64_t size = 0,
    int32_t const strandID = 0
)
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        h5_zone_handle = HDF5::Open_Zone(h5_file_handle, zone_header);
        break;
    case szplt:
        tec_output_zone =
            SZLIO::create_zone(tec_file_handle, strandID, zone_header, var_types, size, solution_time);
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        break;
    }
}

void FileReadWriter::close_zone()
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        HDF5::Close_Handle(h5_zone_handle);
        break;
    case szplt:
        tec_output_zone = 0;
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        break;
    }
}

bool FileReadWriter::get_file_is_open()
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        return h5_file_handle != H5I_INVALID_HID;
        break;
    case szplt:
        return tec_file_handle != NULL;
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        return false;
        break;
    }
}

bool FileReadWriter::get_zone_is_open()
{
    switch (file_ext)
    {
    case dat:
        /* code */
        break;
    case h5:
    case h5part:
        return h5_zone_handle != H5I_INVALID_HID;
        break;
    case szplt:
        return tec_output_zone != 0;
        break;
    default:
        cerr << "Incorrect file extension defined" << endl;
        return false;
        break;
    }
}

RestartFile::RestartFile(SIM const& svar) : FileReadWriter(h5, svar){};

/**
 * @brief Write a .h5 file to provide restart data.
 *
 * @param svar Structure of simulation variables
 * @param pnp1 SoA containing SPH particles at time n+1 (at the time of writing, n = n+1)
 * @param limits array of block indexes for particles in the global array.
 * @return void
 */
void RestartFile::write_restart_file(SIM& svar, SPHState const& pnp1, LIMITS const& limits)
{
    assert(file_ext == h5);

    std::string filename = svar.io.output_prefix + "_particles.h5";
    open_file(filename);

    write_file_attributes(svar);

    /* Write the actual data */
    size_t blockID = 0;

    // Boundary particles
    if (svar.n_bound_blocks > 0)
    {
        for (size_t ii = 0; ii < svar.n_bound_blocks; ++ii)
        {
            std::string zoneHeader = "Boundary " + std::to_string(ii);

            // Create the group
            h5_zone_handle = HDF5::Open_Zone(h5_file_handle, zoneHeader);

            write_zone_attributes(limits[ii]);
            write_restart_zone(svar.scale, pnp1, limits[ii].index.first, limits[ii].index.second);

            close_zone();

            blockID++;
        }
    }

    // Fluid particles
    if (svar.n_fluid_blocks > 0)
    {
        size_t flublock = 0;
        for (size_t ii = svar.n_bound_blocks; ii < svar.n_bound_blocks + svar.n_fluid_blocks; ++ii)
        {
            std::string zoneHeader = "Fluid " + std::to_string(flublock);

            // Create the group
            h5_zone_handle = HDF5::Open_Zone(h5_file_handle, zoneHeader);

            write_zone_attributes(limits[ii]);
            write_restart_zone(svar.scale, pnp1, limits[ii].index.first, limits[ii].index.second);
            if (limits[ii].block_type == inletZone)
                HDF5::Write_Inlet_Data(h5_zone_handle, limits[ii]);

            close_zone();

            blockID++;
            flublock++;
        }
    }

    // Close the files
    close_file();

    // If the output_prefix differs from the restart prefix, write a new restart prefix using
    // that
    if (svar.io.restart && svar.io.output_prefix == svar.io.restart_prefix)
    {
        svar.io.restart_header_written = 1;
    }

    if (!svar.io.restart_header_written)
    {
        /* Append a restart prefix to the end of the para file */
        FILE* para = fopen(svar.io.input_file.c_str(), "a");
        if (para == NULL)
        {
            printf("Failed to reopen para file\n");
            exit(-1);
        }
        std::time_t now = std::time(0);
        char* time = ctime(&now);
        fprintf(para, "\n");
        fprintf(para, "    solver at %s", time);
        fprintf(para, "                                          ");
        fprintf(para, "SPH restart prefix: %s\n", svar.io.output_prefix.c_str());
        fclose(para);
        svar.io.restart_header_written = 1;
    }
    printf("Finished writing restart HDF5 file.\n");
}

int RestartFile::write_restart_zone(
    real const& scale, SPHState const& pnp1, size_t const& start, size_t const& end
)
{
    // Write out all the parameters for a timestep

    /* Need to write position, velocity, acceleration, pressure, density gradient, mass, boundary
     * flag */
    /* Write as a 2D array, or 1D arrays for vectors? 1D for now */
    size_t length = end - start;
    // hsize_t dims[] = {length,2};
    hsize_t dims[] = {length};
    vector<real> vec(length);
    int32_t var_count = 0;

    /* Position */
    vector<StateVecD> vvec(length);
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vvec[ii - start] = pnp1[ii].xi / scale;

    write_scalar_var("Position coordinate", vvec, var_count);

    /* Velocity */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vvec[ii - start] = pnp1[ii].v;

    write_scalar_var("Velocity", vvec, var_count);

    /* Acceleration */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vvec[ii - start] = pnp1[ii].acc;

    write_scalar_var("Acceleration", vvec, var_count);
    vvec.clear();

/* Pressure */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].p;

    write_scalar_var("Pressure", vec, var_count);

/* Density */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].rho;

    write_scalar_var("Density", vec, var_count);

/* Density gradient */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].Rrho;

    write_scalar_var("Density gradient", vec, var_count);

/* Mass */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].m;

    write_scalar_var("Mass", vec, var_count);

    vec.clear();

/* Particle spacing */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].dx;

    write_scalar_var("Particle spacing", vec, var_count);

    vec.clear();

    /* Boundary flag */
    vector<uint> uvec(length);
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        uvec[ii - start] = pnp1[ii].b;

    write_scalar_var("Boundary condition", uvec, var_count);
    uvec.clear();

    /* Particle ID */
    vector<size_t> svec(length);
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        svec[ii - start] = pnp1[ii].part_id;

    write_scalar_var("Particle ID", svec, var_count);
    svec.clear();

    /* Fluid ID */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        svec[ii - start] = pnp1[ii].fluid_id;

    write_scalar_var("Fluid ID", svec, var_count);
    svec.clear();

    /* Cell ID */
    vector<long> lvec(length);
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        lvec[ii - start] = pnp1[ii].cellID;

    write_scalar_var("Cell ID", lvec, var_count);

    /* Cell velocity */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vvec[ii - start] = pnp1[ii].cellV;

    write_scalar_var("Cell velocity", vvec, var_count);

/* Cell density */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].cellRho;

    write_scalar_var("Cell density", vec, var_count);

/* Cell pressure */
#pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
        vec[ii - start] = pnp1[ii].cellP;

    write_scalar_var("Cell pressure", vec, var_count);

    lvec.clear();
}

DataFile::DataFile(file_extension const& _file_ext, SIM const& svar, std::string const& _data_type)
    : FileReadWriter(_file_ext, svar), output_var(svar.io.output_variables),
      output_prefix(svar.io.output_prefix), data_type(_data_type)
{
    // Default file name.
    data_lower = data_type;
    std::transform(data_lower.begin(), data_lower.end(), data_lower.begin(), ::tolower);
};

string DataFile::get_filename(real time = -1)
{
    if (time < 0)
        return output_prefix + "_" + data_lower + ext_str;
    else
    {
        std::ostringstream oss;
        oss << std::setprecision(5) << time;
        std::string time_str = oss.str();
        return output_prefix + "_time_" + time_str + "_" + data_lower + ext_str;
    }
}

void DataFile::write_zone(
    SPHState const& pnp1, real const& rho_rest, size_t const& start, size_t const& end
)
{

    int32_t var_count = 1;
    int32_t imax = end - start;

    vector<real> vec(imax);
    // Restart essential variables. These should not be zero, but could possibly in future
    if (output_var.at("pos-vec").write)
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].xi / scale;

        string name;
        switch (file_ext)
        {
        case szplt:
        {
            name = "Position coordinate";
            break;
        }
        default:
        {
            // For h5part want just an empty string and it will fill in x/y/z for paraview.
            name = "";
            break;
        }
        }

        write_scalar_var(name, vvec, var_count);
    }

    if (output_var.at("vel-vec").write)
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].v;

        write_scalar_var("Velocity", vvec, var_count);
    }

    if (output_var.at("acc-vec").write)
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].acc;

        write_scalar_var("Acceleration", vvec, var_count);
    }

    if (output_var.at("press").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].p;
        write_scalar_var("Pressure", vec, var_count);
    }

    if (output_var.at("dRho").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].Rrho;
        write_scalar_var("Density gradient", vec, var_count);
    }

    if (output_var.at("part_id").write)
    {
        vector<int> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].part_id;
        write_scalar_var("Particle ID", uvec, var_count);
    }

    if (output_var.at("cellID").write)
    {
        vector<int> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].cellID;
        write_scalar_var("Cell ID", uvec, var_count);
    }

    if (output_var.at("bound").write)
    {
        vector<uint8_t> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].b;
        write_scalar_var("Particle status condition", uvec, var_count);
    }

    // Non essential variables, but to provide further information
    if (output_var.at("dens").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].rho;
        write_scalar_var("Density", vec, var_count);
    }

    if (output_var.at("densVar").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = 100.0 * (pnp1[ii].rho / rho_rest - 1.0);
        write_scalar_var("Density Variation", vec, var_count);
    }

    if (output_var.at("vmag").write)
    {
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].v.norm();
        write_scalar_var("Velocity magnitude", vec, var_count);
    }

    if (output_var.at("surf").write)
    {
        vector<uint8_t> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].surf;
        write_scalar_var("Surface flag", uvec, var_count);
    }

    if (output_var.at("surfZ").write)
    {
        vector<uint8_t> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].surfzone;
        write_scalar_var("Surface zone flag", uvec, var_count);
    }

    if (output_var.at("aero-mag").write)
    {
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].Af.norm();
        write_scalar_var("Aerodynamic force magnitude", vec, var_count);
    }

    if (output_var.at("aero-vec").write)
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].Af;

        write_scalar_var("Aerodynamic force", vvec, var_count);
    }

    if (output_var.at("curv").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].curve;
        write_scalar_var("Curvature", vec, var_count);
    }

    if (output_var.at("occl").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].woccl;
        write_scalar_var("Occlusion factor", vec, var_count);
    }

    if (output_var.at("cellP").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellP;
        write_scalar_var("Cell pressure", vec, var_count);
    }

    if (output_var.at("cellRho").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellP;
        write_scalar_var("Cell density", vec, var_count);
    }

    if (output_var.at("cellV-mag").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellV.norm();
        write_scalar_var("Cell velocity magnitude", vec, var_count);
    }

    if (output_var.at("cellV-vec").write)
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].cellV;

        write_scalar_var("Cell velocity", vvec, var_count);
    }

    if (output_var.at("dsphG-vec").write) // Delta SPH density gradient
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].gradRho;

        write_scalar_var("d-SPH density gradient", vvec, var_count);
    }

    if (output_var.at("lam").write) // Lambda eigenvalue of L matrix
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].lam;
        write_scalar_var("Lmabda eigenvalue", vec, var_count);
    }

    if (output_var.at("lam-nb").write) // Lambda eigenvalue of L matrix without boundary
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].lam_nb;
        write_scalar_var("Lambda eigenvalue without boundary", vec, var_count);
    }

    if (output_var.at("colour").write) // Surface colour function
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].colour;
        write_scalar_var("Colour function", vec, var_count);
    }

    if (output_var.at("colour-G").write) // Gradient of the colour function for CSF
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].colourG;
        write_scalar_var("Colour function gradient", vec, var_count);
    }

    if (output_var.at("norm-vec").write) // Surface normal
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].norm;

        write_scalar_var("Surface normal", vvec, var_count);
    }

    if (output_var.at("shiftV-mag").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].vPert.norm();
        write_scalar_var("Shifting velocity magnitude", vec, var_count);
    }

    if (output_var.at("shiftV-vec").write)
    {
        vector<StateVecD> vvec(imax);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vvec[ii - start] = pnp1[ii].vPert;

        write_scalar_var("Shifting velocity", vvec, var_count);
    }
}

H5PartFile::H5PartFile(SIM const& svar, std::string const& _data_type)
    : DataFile(h5part, svar, _data_type){};

/**
 * @brief Write a timestep of the simulation to .h5part files.
 *
 * @param svar Structure of simulation variables
 * @param pnp1 SoA containing SPH particles
 * @return void
 */
void H5PartFile::write_timestep(
    SPHState const& pnp1, LIMITS const& limits, size_t start, size_t end, uint frame, real const& time
)
{
    // Open the file for each timestep. Will append if the file exists.
    string _file_name = get_filename();
    if (!single_file)
        _file_name = get_filename(time);

    open_file(_file_name);

    // Ensure that files are open
    assert(get_file_is_open());

    std::string zone_header = "Step#" + std::to_string(frame);
    create_zone(zone_header);
    assert(get_zone_is_open());

    // Write attributes for the zone
    write_zone_attribute("TimeValue", time);

    // Write the arrays to the zone
    write_zone(pnp1, limits[0].fluid.rho_rest, start, end);

    // close the zone and return
    close_zone();

    // Close the file for h5part as it will be reopened and appended to.
    close_file();
}

SZLFile::SZLFile(
    SIM const& svar, std::string const& _data_type, size_t const& _start_block, size_t const& _end_block
)
    : DataFile(szplt, svar, _data_type), start_block(_start_block), end_block(_end_block)
{
    if (single_file)
    {
        // Open the file now. Default file name is correct
        string zone_name = data_type + " Particles";
        open_file(get_filename(), zone_name);
    }
};

void SZLFile::write_timestep(SPHState const& pnp1, LIMITS const& limits, real const& time)
{
    if (!single_file)
    {
        // Open the file for this time.
        string _file_name = get_filename(time);
        string zone_name = data_type + " Particles";
        open_file(_file_name, zone_name);
    }

    write_timestep_data(pnp1, limits, time);

    if (!single_file)
    {
        close_file();
    }
}

void SZLFile::flush_file() { SZLIO::flush_file(tec_file_handle); }

void SZLFile::combine_file()
{
    if (!single_file)
        Combine_SZPLT(get_filename());
};

void SZLFile::write_timestep_data(SPHState const& pnp1, LIMITS const& limits, real const& time)
{
    // Ensure that files are open
    assert(get_file_is_open());

    std::ostringstream oss;
    oss << std::setprecision(5) << time;
    std::string time_str = oss.str();

    for (size_t block = start_block; block < end_block; block++)
    {
        bound_block const& limit = limits[block];
        size_t start = limit.index.first;
        size_t end = limit.index.second;
        int32_t const size = end - start;
        string title = data_type + "_" + std::to_string(block) + "_" + limits[block].name + "_time_" +
                       time_str + "s";

        int32_t strand_id = static_cast<int32_t>(block + 1);

        int32_t outputZone =
            SZLIO::create_zone(tec_file_handle, strand_id, title, var_types, size, time);

        write_zone_attributes(limit);
        write_zone(pnp1, limit.fluid.rho_rest, start, end);
        flush_file();
        close_zone();
    }
}
