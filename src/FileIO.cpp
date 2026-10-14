/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

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
    FileWriter(file_extension const& _file_ext) { file_ext = _file_ext; };

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

    int write_zone(real const& scale, SPHState const& pnp1, size_t const& start, size_t const& end)
    {
        // Write out all the parameters for a timestep

        /* Need to write position, velocity, acceleration, pressure, density gradient, mass, boundary
         * flag */
        /* Write as a 2D array, or 1D arrays for vectors? 1D for now */
        size_t length = end - start;
        // hsize_t dims[] = {length,2};
        hsize_t dims[] = {length};
        vector<real> vec(length);

        /* Position */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].xi(dim) / scale;

            string name = "Position coordinate " + get_dim(dim);
            write_scalar_var(name, vec);
        }

        /* Velocity */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].v(dim);

            string name = "Velocity " + get_dim(dim);
            write_scalar_var(name, vec);
        }

        /* Acceleration */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].acc(dim);

            string name = "Acceleration " + get_dim(dim);
            write_scalar_var(name, vec);
        }

/* Pressure */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].p;

        write_scalar_var("Pressure", vec);

/* Density */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].rho;

        write_scalar_var("Density", vec);

/* Density gradient */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].Rrho;

        write_scalar_var("Density gradient", vec);

/* Mass */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].m;

        write_scalar_var("Mass", vec);

        vec.clear();

/* Particle spacing */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].dx;

        write_scalar_var("Particle spacing", vec);

        vec.clear();

        /* Boundary flag */
        vector<uint> uvec(length);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].b;

        write_scalar_var("Boundary condition", uvec);
        uvec.clear();

        /* Particle ID */
        vector<size_t> svec(length);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            svec[ii - start] = pnp1[ii].part_id;

        write_scalar_var("Particle ID", svec);
        svec.clear();

        /* Fluid ID */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            svec[ii - start] = pnp1[ii].fluid_id;

        write_scalar_var("Fluid ID", svec);
        svec.clear();

        /* Cell ID */
        vector<long> lvec(length);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            lvec[ii - start] = pnp1[ii].cellID;

        write_scalar_var("Cell ID", lvec);

        /* Cell velocity */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellV(dim);

            string name = "Cell velocity " + get_dim(dim);
            write_scalar_var(name, vec);
        }

/* Cell density */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellRho;

        write_scalar_var("Cell density", vec);

/* Cell pressure */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellP;

        write_scalar_var("Cell pressure", vec);

        lvec.clear();
    }

  private:
    file_extension file_ext;
    hid_t h5_file_handle;
    void* tec_file_handle;

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
    void write_scalar_var(std::string const& var_name, std::vector<T> const& var_data)
    {
        switch (file_ext)
        {
        case dat:
            /* code */
            break;
        case h5:
            hsize_t dims[] = {var_data.size()};
            HDF5::Write_Variable_Scalar(h5_file_handle, var_name, dims, 0, var_data);
            break;
        case szplt:
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
