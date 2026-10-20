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
        case h5part:
            break;
        default:
            break;
        }
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
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].xi(dim) / scale;

            string name = "Position coordinate " + get_dim(dim);
            write_scalar_var(name, vec, var_count);
        }

        /* Velocity */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].v(dim);

            string name = "Velocity " + get_dim(dim);
            write_scalar_var(name, vec, var_count);
        }

        /* Acceleration */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].acc(dim);

            string name = "Acceleration " + get_dim(dim);
            write_scalar_var(name, vec, var_count);
        }

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
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellV(dim);

            string name = "Cell velocity " + get_dim(dim);
            write_scalar_var(name, vec, var_count);
        }

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
                vvec[ii] = pnp1[ii].xi / scale;

            write_vector_var("Position coordinate", vvec, var_count);
        }

        if (output_var.at("vel-vec").write)
        {
            vector<StateVecD> vvec(imax);
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vvec[ii] = pnp1[ii].v;

            write_vector_var("Velocity", vvec, var_count);
        }

        if (output_var.at("acc-vec").write)
        {
            vector<StateVecD> vvec(imax);
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vvec[ii] = pnp1[ii].acc;

            write_vector_var("Acceleration", vvec, var_count);
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
                vvec[ii] = pnp1[ii].Af;

            write_vector_var("Aerodynamic force", vvec, var_count);
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
                vvec[ii] = pnp1[ii].cellV;

            write_vector_var("Cell velocity", vvec, var_count);
        }

        if (output_var.at("dsphG-vec").write) // Delta SPH density gradient
        {
            vector<StateVecD> vvec(imax);
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vvec[ii] = pnp1[ii].gradRho;

            write_vector_var("d-SPH density gradient", vvec, var_count);
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
                vvec[ii] = pnp1[ii].norm;

            write_vector_var("Surface normal", vvec, var_count);
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
                vvec[ii] = pnp1[ii].vPert;

            write_vector_var("Shifting velocity", vvec, var_count);
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
            HDF5::Write_Variable_Scalar(h5_file_handle, var_name, dims, var_data);
            break;
        case szplt:
            SZLIO::Write_Scalar_Vector(
                tec_file_handle, tec_output_zone, var_count, 0, imax, var_data, var_name
            );
            break;
        default:
            break;
        }
    }

    void write_vector_var(
        std::string const& var_name, std::vector<StateVecD> const& var_data, int32_t& var_count
    )
    {
        switch (file_ext)
        {
        case dat:
            /* code */
            break;
        case h5:
            HDF5::Write_Vector_Data(h5_file_handle, var_name, var_data);
            break;
        case szplt:
            SZLIO::Write_State_Vector(tec_file_handle, tec_output_zone, var_count, var_data, var_name);
            break;
        default:
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
