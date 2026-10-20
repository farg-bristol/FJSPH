/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "BinaryIO.h"
#include "H5IO.h"
#include "Helper_Functions.h"
#include "IO.h"
#include "IOFunctions.h"
#include "VarDefs.h"

#include <hdf5.h>

/* Create a wrapper class to write to files, so that one location is needed for updating files. */

class FileWriter
{
  public:
    // Need to handle Tecplot, and HDF5, which have different types of things
    // Create with a flag to give an ID for what kind of file to interact with?
    FileWriter(
        file_extension const& _file_ext, OutputMap const& _outvar, real const& _scale, real const& _rho
    )
        : file_ext(_file_ext), output_var(_outvar), scale(_scale), rho_rest(_rho){};

    int open_file(std::string const& filename, std::string const& zone_name, SIM const& svar)
    {
        // Get the correct file extension string
        string ext = "";
        switch (file_ext)
        {
        case dat:
            ext = ".dat";
            break;
        case h5:
            ext = ".h5";
            break;
        case szplt:
            ext = ".szplt";
            break;
        case h5part:
            ext = ".h5part";
            break;
        default:
            break;
        }

        /*
         * Create a new file. If file exists its contents will be overwritten.
         */
        std::filesystem::path _filename{filename};
        _filename.replace_extension(ext);
        if (std::filesystem::exists(_filename))
        { // Remove the file
            std::filesystem::remove(_filename);
        }

        switch (file_ext)
        {
        case dat:
        {
            /* code */
            break;
        }
        case h5:
        {
            h5_file_handle = HDF5::Open_HDF5_File(_filename.string(), svar);
            break;
        }
        case szplt:
        {
            if (Open_PLT_File(_filename.string(), zone_name, svar, tec_file_handle))
            {
                cout << "Failed to open " << _filename.string() << endl;
                exit(-1);
            }
            break;
        }
        case h5part:
        {
            h5_file_handle = HDF5::Open_HDF5_File(_filename.string(), svar);
            break;
        }
        }
    }

    void write_file_attributes(SIM const& svar)
    {
        // Write crucial simulation auxiliary data. Potential to expand
        write_attribute("Simulation current time", svar.integrator.current_time);
        write_attribute("Particle index to add", svar.part_id);
        write_attribute("Number of boundary blocks", svar.n_bound_blocks);
        write_attribute("Number of fluid blocks", svar.n_fluid_blocks);
        write_attribute("Number of boundary points", svar.bound_points);
        write_attribute("Number of fluid points", svar.fluid_points);
        write_attribute("Current frame", svar.integrator.current_frame);

        write_attribute("Input para filename", svar.io.input_file);
        write_attribute("Input fluid definition filename", svar.io.input_fluid_file);
        write_attribute("Input boundary definition filename", svar.io.input_bound_file);

        write_attribute("Primary grid face filename", svar.io.tau_mesh);
        write_attribute("Boundary mapping filename", svar.io.tau_bmap);
        write_attribute("Restart-data prefix", svar.io.tau_sol);
        write_attribute("Dimension offset vector", svar.io.offset_axis);
        write_attribute("Solution angle of attack", svar.io.angle_alpha);
        write_attribute("Grid scale", svar.scale);

        write_attribute("OpenFOAM input directory", svar.io.foam_dir);
        write_attribute("OpenFOAM solution directory", svar.io.foam_sol);
        write_attribute("OpenFOAM binary", svar.io.foam_is_binary);
        write_attribute("Label size", svar.io.foam_label_size);
        write_attribute("Scalar size", svar.io.foam_scalar_size);
        write_attribute("OpenFOAM is buoyant", svar.io.foam_buoyant_sim);
        write_attribute("OpenFOAM is foam_is_incompressible", svar.io.foam_is_incomp);

        write_attribute("VLM definition filename", svar.io.vlm_file);

        write_attribute("Single file for output", svar.io.single_file);
        write_attribute("Output files prefix", svar.io.output_prefix);
        write_attribute("SPH restart prefix", svar.io.restart_prefix);
        write_attribute("SPH frame time interval", svar.integrator.frame_time_interval);
        write_attribute("SPH frame count", svar.integrator.max_frames);
        write_attribute("SPH previous frame time", svar.integrator.last_frame_time);
        write_attribute("SPH output encoding", svar.io.out_encoding);
        write_attribute("Variable list", svar.io.output_names);

        /* Fluid data */
        write_attribute("Reference density", svar.air.rho_g);
        write_attribute("Init hydrostatic pressure", svar.init_hydro_pressure);
        write_attribute("Hydrostatic height", svar.hydro_height);

        /* Aerodynamic data */
        write_attribute("Reference velocity", svar.air.v_ref);
        write_attribute("Reference pressure", svar.air.p_ref);
        write_attribute("Reference Mach number", svar.air.M_ref);
        write_attribute("Reference temperature", svar.air.temp_g);
        write_attribute("Gas constant gamma", svar.air.gamma);

        /* Simulation settings */
        write_attribute("SPH integration solver", svar.integrator.solver_name);
        write_attribute("SPH boundary solver", svar.integrator.bound_solver);
        write_attribute("SPH solver minimum residual", svar.integrator.min_residual);
        write_attribute("SPH maximum timestep", svar.integrator.delta_t_max);
        write_attribute("SPH minimum timestep", svar.integrator.delta_t_min);
        write_attribute("SPH maximum CFL", svar.integrator.cfl_max);
        write_attribute("SPH minimum CFL", svar.integrator.cfl_min);
        write_attribute("SPH CFL condition", svar.integrator.cfl);
        write_attribute("SPH unstable CFL step", svar.integrator.cfl_step);
        write_attribute("SPH unstable CFL count limit", svar.integrator.n_unstable_limit);
        write_attribute("SPH stable CFL count limit", svar.integrator.n_stable_limit);
        write_attribute("SPH stable CFL count iteration factor", svar.integrator.subits_factor);
        write_attribute("SPH maximum shifting velocity", svar.integrator.max_shift_vel);
        write_attribute("SPH stable CFL count", svar.integrator.n_stable);
        write_attribute("SPH unstable CFL count", svar.integrator.n_unstable);

        write_attribute("SPH density variation", svar.rho_var);
        write_attribute("SPH maximum density", svar.rho_max);
        write_attribute("SPH minimum density", svar.rho_min);
        write_attribute("SPH Newmark Beta iteration limit", svar.integrator.max_subits);
        write_attribute("SPH gravity vector", svar.grav);

        write_attribute("SPH smoothing length factor", svar.H_fac);
        write_attribute("SPH aerodynamic case", svar.air.aero_case);
        write_attribute("SPH SP diameter definition", svar.air.use_dx);
        write_attribute("SPH use TAB deformation", svar.air.use_TAB_def);
        write_attribute("SPH global offset coordinate", svar.offset_vec);
        write_attribute("SPH maximum particle count", svar.max_points);
        write_attribute("SPH aerodynamic cutoff value", svar.air.lam_cutoff);
        write_attribute("SPH freestream velocity", svar.air.v_inf);
        write_attribute("SPH restart fit tolerance", svar.io.restart_tol);

        /* Particle tracking settings */
        write_attribute("Transition to IPT", svar.ipt.using_ipt);
        write_attribute("Velocity equation order", svar.ipt.ipt_eq_order);
        write_attribute("SPH tracking conversion x coordinate", svar.ipt.max_x_sph);
        write_attribute("Maximum x trajectory coordinate", svar.ipt.max_x);
        write_attribute("Particle scatter output", svar.ipt.part_out);
        write_attribute("Particle streak output", svar.ipt.streak_out);
        write_attribute("Particle cell intersection output", svar.ipt.cells_out);
    }

    int
    write_restart_zone(real const& scale, SPHState const& pnp1, size_t const& start, size_t const& end)
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

    void write_data_zone(SPHState const& pnp1, size_t const& start, size_t const& end)
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

  private:
    file_extension file_ext;
    OutputMap output_var;
    real scale = 1.0;
    real rho_rest = -1;

    hid_t h5_file_handle;
    void* tec_file_handle;
    int32_t tec_output_zone;

    inline string get_dim(uint const& dim)
    {
        switch (dim)
        {
        case 0:
            return "x";
            break;
        case 1:
            return "y";
            break;
        case 2:
            return "z";
            break;
        default:
            return "x";
            break;
        }
        return "x";
    }

    template <typename T> void write_attribute(std::string const& attr_name, T const& attr_data)
    {
        switch (file_ext)
        {
        case dat:
            /* code */
            break;
        case h5:
            HDF5::Write_Attribute(h5_file_handle, var_name, var_data);
            break;
        case szplt:
            SZLIO::Write_Scalar_Vector(
                tec_file_handle, tec_output_zone, var_count, 0, imax, var_data, var_name
            );
            break;
        default:
            cerr << "Incorrect file extension defined" << endl;
            break;
        }
    }

    template <typename T>
    void
    write_scalar_var(std::string const& var_name, std::vector<T> const& var_data, int32_t& var_count)
    {
        switch (file_ext)
        {
        case dat:
            /* code */
            break;
        case h5:
            hsize_t dims[] = {var_data.size()};
            HDF5::Write_Variable(h5_file_handle, var_name, dims, var_data);
            break;
        case szplt:
            SZLIO::Write_Vector(
                tec_file_handle, tec_output_zone, var_count, 0, imax, var_data, var_name
            );
            break;
        default:
            cerr << "Incorrect file extension defined" << endl;
            break;
        }
    }
};

class FileReader
{
  public:
    // Need to handle Tecplot, and HDF5, which have different types of things
    // Create with a flag to give an ID for what kind of file to interact with?
    FileReader(file_extension const& _file_ext) { file_ext = _file_ext; };

    int open_file(std::string const& filename, SIM const& svar)
    {
        switch (file_ext)
        {
        case dat:
            /* code */
            break;
        case h5:
            h5_file_handle = HDF5::Open_HDF5_File(filename, svar);
            break;
        case szplt:
            break;
        default:
            break;
        }
    }

    int read_zone()
    {
        // Write out all the parameters for a timestep
    }

  private:
    file_extension file_ext;
    hid_t h5_file_handle;
    void* tec_file_handle;
};
