/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "H5IO.h"
#include "Helper_Functions.h"
#include "IO.h"
#include "IOFunctions.h"
#include "Kernel.h"

#include <filesystem>
#include <hdf5.h>

#define asize(array) (sizeof(array) / sizeof(array[0]))

namespace HDF5
{
    void Write_Uint_Attribute(int64_t const& fout, std::string const& attr_name, size_t const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate(fout, (attr_name).c_str(), H5T_NATIVE_ULONG, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_ULONG, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void Write_Uint_Attribute(int64_t const& fout, std::string const& attr_name, uint const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate(fout, (attr_name).c_str(), H5T_NATIVE_UINT, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_UINT, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void Write_Int_Attribute(int64_t const& fout, std::string const& attr_name, int const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate(fout, (attr_name).c_str(), H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_INT, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void Write_Real_Attribute(int64_t const& fout, std::string const& attr_name, double const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate(fout, (attr_name).c_str(), H5T_IEEE_F64LE, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_DOUBLE, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void
    Write_Vector_Attribute(int64_t const& fout, std::string const& attr_name, StateVecD const& attr_data)
    {
        std::string name_ = attr_name /* + " x" */;
        hsize_t dims[1] = {SIMDIM};
        herr_t retval;
        int64_t aid = H5Screate_simple(1, dims, NULL);
        int64_t attr = H5Acreate(fout, (name_).c_str(), H5T_IEEE_F64LE, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_DOUBLE, &attr_data[0])))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void Write_String_Attribute(
        int64_t const& fout, std::string const& attr_name, std::string const& attr_data
    )
    {
        if (attr_data.length() > 0)
        {
            char data[attr_data.length() + 1];
            strcpy(data, attr_data.c_str());

            int64_t attr_type = H5Tcopy(H5T_FORTRAN_S1);
            herr_t retval = H5Tset_size(attr_type, attr_data.length());
            int64_t memtype = H5Tcopy(H5T_C_S1);
            retval = H5Tset_size(memtype, attr_data.length() + 1);

            int64_t attr_space = H5Screate(H5S_SCALAR);

            int64_t attr =
                H5Acreate(fout, attr_name.c_str(), attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
            if ((retval = H5Awrite(attr, memtype, data)))
            {
                printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
                exit(-1);
            }
            // retval = H5Sclose(attr_type);
            retval = H5Sclose(attr_space);
            retval = H5Aclose(attr);
        }
    }

    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<size_t>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start][0])
            ))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<int>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start][0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<float>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start][0])
            ))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<double>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start][0])
            ))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    // #if UINTPTR_MAX == 0xffffffffffffffff
    /* 64-bit pointers */
    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<size_t> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }
    // #elif UINTPTR_MAX == 0xffffffff
    /* 32-bit pointers */
    // void Write_Variable_Scalar(int64_t const& fout, std::string const& var_name,
    //         hsize_t const* dims, size_t const& start, std::vector<size_t> const& var_data)
    // {
    //     herr_t retval;
    //     int64_t  dataspace = H5Screate(H5S_SIMPLE);
    //     retval = H5Sset_extent_simple(dataspace,1,dims,NULL);
    //     int64_t dataset = H5Dcreate2(fout, var_name.c_str(), H5T_NATIVE_UINT, dataspace, H5P_DEFAULT,
    //     H5P_DEFAULT, H5P_DEFAULT); if ((retval = H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
    //     H5P_DEFAULT, &var_data[start])))
    //     {
    //         printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
    //         exit(-1);
    //     }
    //     retval = H5Sclose(dataspace);
    //     retval = H5Dclose(dataset);
    // }

    // #endif

    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<uint> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<int> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<long> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_LONG, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<float> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<double> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Vector_Data(
        int64_t const& fout, std::string const& name, std::vector<StateVecD> const& data,
        size_t const& start, size_t const& end
    )
    {
        hsize_t length = end - start;
        std::vector<double> vect(length);
        hsize_t dim[] = {length};
        std::string name_;

#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vect[ii - start] = data[ii][0];

        name_ = name.empty() ? "x" : name + " x";
        Write_Variable_Scalar(fout, name_.c_str(), dim, 0, vect);

#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vect[ii - start] = data[ii][1];

        name_ = name.empty() ? "y" : name + " y";
        Write_Variable_Scalar(fout, name_.c_str(), dim, 0, vect);

        name_ = name.empty() ? "z" : name + " z";
#if SIMDIM == 3
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vect[ii - start] = data[ii][2];

        Write_Variable_Scalar(fout, name_.c_str(), dim, 0, vect);
#else
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vect[ii - start] = 0.0;

        Write_Variable_Scalar(fout, name_.c_str(), dim, 0, vect);
#endif
    }

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

    void Write_Zone_Data(
        int64_t const& fout, real const& scale, SPHState const& pnp1, size_t const& start,
        size_t const& end
    )
    {
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
            Write_Variable_Scalar(fout, name, dims, 0, vec);
        }

        /* Velocity */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].v(dim);

            string name = "Velocity " + get_dim(dim);
            Write_Variable_Scalar(fout, name, dims, 0, vec);
        }

        /* Acceleration */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].acc(dim);

            string name = "Acceleration " + get_dim(dim);
            Write_Variable_Scalar(fout, name, dims, 0, vec);
        }

/* Pressure */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].p;

        Write_Variable_Scalar(fout, "Pressure", dims, 0, vec);

/* Density */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].rho;

        Write_Variable_Scalar(fout, "Density", dims, 0, vec);

/* Density gradient */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].Rrho;

        Write_Variable_Scalar(fout, "Density gradient", dims, 0, vec);

/* Mass */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].m;

        Write_Variable_Scalar(fout, "Mass", dims, 0, vec);

        vec.clear();

        /* Boundary flag */
        vector<uint> uvec(length);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].b;

        Write_Variable_Scalar(fout, "Boundary condition", dims, 0, uvec);
        uvec.clear();

        /* Particle ID */
        vector<size_t> svec(length);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            svec[ii - start] = pnp1[ii].part_id;

        Write_Variable_Scalar(fout, "Particle ID", dims, 0, svec);
        svec.clear();

        /* Cell ID */
        vector<long> lvec(length);
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            lvec[ii - start] = pnp1[ii].cellID;

        Write_Variable_Scalar(fout, "Cell ID", dims, 0, lvec);

        /* Cell velocity */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellV(dim);

            string name = "Cell velocity " + get_dim(dim);
            Write_Variable_Scalar(fout, name, dims, 0, vec);
        }

/* Cell density */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellRho;

        Write_Variable_Scalar(fout, "Cell density", dims, 0, vec);

/* Cell pressure */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellP;

        Write_Variable_Scalar(fout, "Cell pressure", dims, 0, vec);

        lvec.clear();
    }

    void Write_Inlet_Data(int64_t const& fout, bound_block const& limits)
    {
        size_t length = limits.back.size();
        hsize_t dim[] = {length};
        hsize_t dims[] = {length * 4};

        Write_Variable_Scalar(fout, "Back index", dim, 0, limits.back);

        vector<size_t> buff(length * 4);
        size_t accm(0);
        for (size_t ii = 0; ii < length; ii++)
            for (size_t jj = 0; jj < 4; jj++)
                buff[accm++] = limits.buffer[ii][jj];

        Write_Variable_Scalar(fout, "Buffer index", dims, 0, buff);
    }

    /*************************************************************************/
    /*************************** READ FUNCTIONS ******************************/
    /*************************************************************************/
    void Read_String_Attribute(int64_t const& fin, std::string const& attr_name, std::string& attr_data)
    {
        herr_t retval;

        if (H5Aexists(fin, attr_name.c_str()))
        {
            int64_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
            int64_t attr_type = H5Aget_type(attr);

            size_t sdim = H5Tget_size(attr_type) + 1;
            int64_t memtype = H5Tcopy(H5T_C_S1);
            retval = H5Tset_size(memtype, sdim);

            char* data;
            data = (char*)malloc(sdim * sizeof(char));
            if ((retval = H5Aread(attr, memtype, data)))
            {
                printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
                exit(-1);
            }
            H5Aclose(attr);
            attr_data = data;
        }
    }

    void Read_Uint_Attribute(int64_t const& fin, std::string const& attr_name, size_t& attr_data)
    {
        herr_t retval;
        int64_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        size_t data;
        if ((retval = H5Aread(attr, H5T_NATIVE_ULONG, &data)))
        {
            printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Aclose(attr);
        attr_data = data;
    }

    void Read_Uint_Attribute(int64_t const& fin, std::string const& attr_name, uint& attr_data)
    {
        herr_t retval;
        int64_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        uint data;
        if ((retval = H5Aread(attr, H5T_NATIVE_UINT, &data)))
        {
            printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Aclose(attr);
        attr_data = data;
    }

    void Read_Int_Attribute(int64_t const& fin, std::string const& attr_name, int& attr_data)
    {
        herr_t retval;
        int64_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        int data;
        if ((retval = H5Aread(attr, H5T_NATIVE_INT, &data)))
        {
            printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Aclose(attr);
        attr_data = data;
    }

    void Read_Real_Attribute(int64_t const& fin, std::string const& attr_name, double& attr_data)
    {
        herr_t retval;
        int64_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        double data;
        if ((retval = H5Aread(attr, H5T_NATIVE_DOUBLE, &data)))
        {
            printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Aclose(attr);
        attr_data = data;
    }

    void Read_Vector_Attribute(int64_t const& fin, std::string const& attr_name, StateVecD& attr_data)
    {
        herr_t retval;
        int64_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Aget_space(attr);
        hsize_t dims[1];

        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);

        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        double* data = (double*)malloc(dims[0] * sizeof(double));
        if ((retval = H5Aread(attr, H5T_NATIVE_DOUBLE, data)))
        {
            printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Aclose(attr);

        if (SIMDIM > dims[0])
        {
            printf("ERROR: Vector dimension is smaller than simulation dimension. Will cause a fault so "
                   "stopping.\n");
            exit(-1);
        }
        else if (dims[0] != SIMDIM)
        {
            printf("WARNING: dimension of vector does not match simulation dimension\n");
        }

        for (size_t ii = 0; ii < SIMDIM; ii++)
            attr_data[ii] = data[ii];
    }

    void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<size_t>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("WARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<size_t>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<int>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("WARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<int>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<float>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("WARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<float>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<double>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("WARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<double>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<size_t>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<uint>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<int>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<long>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<float>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<double>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void Read_Vector_Data(int64_t const& fin, std::string const& name, std::vector<StateVecD>& var_data)
    {
        std::string name_;

        std::vector<double> vectx, vecty;

        name_ = name.empty() ? "x" : name + " x";
        Read_Variable_Scalar(fin, name_.c_str(), vectx);

        name_ = name.empty() ? "y" : name + " y";
        Read_Variable_Scalar(fin, name_.c_str(), vecty);

#if SIMDIM == 3
        std::vector<double> vectz;
        name_ = name.empty() ? "z" : name + " z";
        Read_Variable_Scalar(fin, name_.c_str(), vectz);

        if (vectx.size() == vecty.size() && vectx.size() == vectz.size())
        {
            var_data.resize(vectx.size());
#pragma omp parallel for
            for (size_t ii = 0; ii < vectx.size(); ii++)
            {
                var_data[ii][0] = vectx[ii];
                var_data[ii][1] = vecty[ii];
                var_data[ii][2] = vectz[ii];
            }
        }
#else
        if (vectx.size() == vecty.size())
        {
            var_data.resize(vectx.size());
#pragma omp parallel for
            for (size_t ii = 0; ii < vectx.size(); ii++)
            {
                var_data[ii][0] = vectx[ii];
                var_data[ii][1] = vecty[ii];
            }
        }
#endif
        else
        {
            printf("Mismatch of vector sizes. Cannot continue.\n");
            exit(-1);
        }
    }

    void Read_Zone_Data(int64_t const& fin, real const& scale, SPHState& pnp1, size_t& length)
    {
        std::vector<StateVecD> avect;
        std::vector<double> dvect;
        std::vector<int> ivect;
        std::vector<size_t> svect;
        std::vector<long> lvect;

        SPHState temp;
        /* Position */
        Read_Vector_Data(fin, "Position coordinate", avect);
        length = avect.size();
        temp.resize(length);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].xi = avect[ii] * scale;

        /* Velocity */
        Read_Vector_Data(fin, "Velocity", avect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].v = avect[ii];

        /* Acceleration */
        Read_Vector_Data(fin, "Acceleration", avect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].acc = avect[ii];

        /* Pressure */
        Read_Variable_Scalar(fin, "Pressure", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].p = dvect[ii];

        /* Density */
        Read_Variable_Scalar(fin, "Density", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].rho = dvect[ii];

        /* Density gradient */
        Read_Variable_Scalar(fin, "Density gradient", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].Rrho = dvect[ii];

        /* Mass */
        Read_Variable_Scalar(fin, "Mass", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].m = dvect[ii];

        /* Boundary flag */
        Read_Variable_Scalar(fin, "Boundary condition", ivect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].b = ivect[ii];

        /* Particle ID */
        Read_Variable_Scalar(fin, "Particle ID", svect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].part_id = svect[ii];

        /* Cell ID */
        Read_Variable_Scalar(fin, "Cell ID", lvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].cellID = lvect[ii];

        /* Cell density */
        Read_Variable_Scalar(fin, "Cell density", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].cellRho = dvect[ii];

        /* Cell pressure */
        Read_Variable_Scalar(fin, "Cell pressure", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].cellP = dvect[ii];

        /* Velocity */
        Read_Vector_Data(fin, "Cell velocity", avect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].cellV = avect[ii];

        // Check each vector is the right size.
        if (avect.size() != dvect.size() || avect.size() != ivect.size())
        {
            printf("ERROR: Mismatch of ingested particle data lengths. Stopping\n");
            exit(-1);
        }
        pnp1.resize(pnp1.size() + length);
        copy_omp(temp.begin(), temp.end(), pnp1.end() - length);
    }

    void Read_Inlet_Data(int64_t const& fin, bound_block& limits)
    {
        Read_Variable_Scalar(fin, "Back index", limits.back);

        vector<size_t> buff;
        Read_Variable_Scalar(fin, "Buffer index", buff);

        limits.buffer.resize(limits.back.size(), vector<size_t>(4, 0));
        size_t accm(0);
        for (size_t ii = 0; ii < limits.back.size(); ii++)
            for (size_t jj = 0; jj < 4; jj++)
                limits.buffer[ii][jj] = buff[accm++];
    }
} // namespace HDF5

void Write_HDF5_Attributes(int64_t const& file, SIM const& svar, FLUID const& fvar, AERO const& avar)
{

    // Write crucial simulation auxiliary data. Potential to expand
    HDF5::Write_Real_Attribute(file, "Simulation current time", svar.current_time);
    HDF5::Write_Uint_Attribute(file, "Particle index to add", svar.part_id);
    HDF5::Write_Uint_Attribute(file, "Number of boundary blocks", svar.n_bound_blocks);
    HDF5::Write_Uint_Attribute(file, "Number of fluid blocks", svar.n_fluid_blocks);
    HDF5::Write_Uint_Attribute(file, "Number of boundary points", svar.bound_points);
    HDF5::Write_Uint_Attribute(file, "Number of fluid points", svar.fluid_points);
    HDF5::Write_Uint_Attribute(file, "Current frame", svar.current_frame);

    HDF5::Write_String_Attribute(file, "Input para filename", svar.input_file);
    HDF5::Write_String_Attribute(file, "Input fluid definition filename", svar.input_fluid_file);
    HDF5::Write_String_Attribute(file, "Input boundary definition filename", svar.input_bound_file);

    HDF5::Write_String_Attribute(file, "Primary grid face filename", svar.tau_mesh);
    HDF5::Write_String_Attribute(file, "Boundary mapping filename", svar.tau_bmap);
    HDF5::Write_String_Attribute(file, "Restart-data prefix", svar.tau_sol);
    HDF5::Write_Int_Attribute(file, "Dimension offset vector", svar.offset_axis);
    HDF5::Write_Real_Attribute(file, "Solution angle of attack", svar.angle_alpha);
    HDF5::Write_Real_Attribute(file, "Grid scale", svar.scale);

    HDF5::Write_String_Attribute(file, "OpenFOAM input directory", svar.foam_dir);
    HDF5::Write_String_Attribute(file, "OpenFOAM solution directory", svar.foam_sol);
    HDF5::Write_Int_Attribute(file, "OpenFOAM binary", svar.foam_is_binary);
    HDF5::Write_Int_Attribute(file, "Label size", svar.foam_label_size);
    HDF5::Write_Int_Attribute(file, "Scalar size", svar.foam_scalar_size);
    HDF5::Write_Int_Attribute(file, "OpenFOAM is buoyant", svar.foam_buoyant_sim);
    HDF5::Write_Int_Attribute(file, "OpenFOAM is foam_is_incompressible", svar.foam_is_incomp);

    HDF5::Write_String_Attribute(file, "VLM definition filename", svar.vlm_file);

    HDF5::Write_Uint_Attribute(file, "Single file for output", svar.single_file);
    HDF5::Write_String_Attribute(file, "Output files prefix", svar.output_prefix);
    HDF5::Write_String_Attribute(file, "SPH restart prefix", svar.restart_prefix);
    HDF5::Write_Real_Attribute(file, "SPH frame time interval", svar.frame_time_interval);
    HDF5::Write_Uint_Attribute(file, "SPH frame count", svar.max_frames);
    HDF5::Write_Real_Attribute(file, "SPH previous frame time", svar.last_frame_time);
    HDF5::Write_Uint_Attribute(file, "SPH output encoding", svar.out_encoding);
    HDF5::Write_String_Attribute(file, "Variable list", svar.output_names);

    /* Fluid data */
    HDF5::Write_Real_Attribute(file, "Reference density", avar.rhog);
    HDF5::Write_Real_Attribute(file, "Reference dispersed density", fvar.rho0);
    HDF5::Write_Real_Attribute(file, "Sutherland reference viscosity", avar.mug);
    HDF5::Write_Real_Attribute(file, "Reference dispersed viscosity", fvar.mu);
    HDF5::Write_Real_Attribute(file, "Reference surface tension", fvar.sig);
    HDF5::Write_Real_Attribute(file, "SPH surface tension contact angle", fvar.contangb);
    HDF5::Write_Int_Attribute(file, "Init hydrostatic pressure", svar.init_hydro_pressure);
    HDF5::Write_Real_Attribute(file, "Hydrostatic height", svar.hydro_height);

    /* Aerodynamic data */
    HDF5::Write_Real_Attribute(file, "Reference velocity", avar.vRef);
    HDF5::Write_Real_Attribute(file, "Reference pressure", avar.pRef);
    HDF5::Write_Real_Attribute(file, "Reference Mach number", avar.MRef);
    HDF5::Write_Real_Attribute(file, "Reference temperature", avar.T);
    HDF5::Write_Real_Attribute(file, "Gas constant gamma", avar.gamma);

    /* Simulation settings */
    HDF5::Write_String_Attribute(file, "SPH integration solver", svar.solver_name);
    HDF5::Write_Uint_Attribute(file, "SPH boundary solver", svar.bound_solver);
    HDF5::Write_Real_Attribute(file, "SPH solver minimum residual", svar.min_residual);
    HDF5::Write_Real_Attribute(file, "SPH maximum timestep", svar.delta_t_max);
    HDF5::Write_Real_Attribute(file, "SPH minimum timestep", svar.delta_t_min);
    HDF5::Write_Real_Attribute(file, "SPH maximum CFL", svar.cfl_max);
    HDF5::Write_Real_Attribute(file, "SPH minimum CFL", svar.cfl_min);
    HDF5::Write_Real_Attribute(file, "SPH CFL condition", svar.cfl);
    HDF5::Write_Real_Attribute(file, "SPH unstable CFL step", svar.cfl_step);
    HDF5::Write_Uint_Attribute(file, "SPH unstable CFL count limit", svar.n_unstable_limit);
    HDF5::Write_Uint_Attribute(file, "SPH stable CFL count limit", svar.n_stable_limit);
    HDF5::Write_Real_Attribute(file, "SPH stable CFL count iteration factor", svar.subits_factor);
    HDF5::Write_Real_Attribute(file, "SPH maximum shifting velocity", svar.max_shift_vel);
    HDF5::Write_Uint_Attribute(file, "SPH stable CFL count", svar.n_stable);
    HDF5::Write_Uint_Attribute(file, "SPH unstable CFL count", svar.n_unstable);

    HDF5::Write_Real_Attribute(file, "SPH background pressure", fvar.backP);
    HDF5::Write_Real_Attribute(file, "SPH starting pressure", fvar.pPress);
    HDF5::Write_Real_Attribute(file, "SPH density variation", fvar.rhoVar);
    HDF5::Write_Real_Attribute(file, "SPH maximum density", fvar.rhoMax);
    HDF5::Write_Real_Attribute(file, "SPH minimum density", fvar.rhoMin);
    HDF5::Write_Real_Attribute(file, "SPH delta coefficient", fvar.delta);

    HDF5::Write_Real_Attribute(file, "SPH artificial viscosity factor", fvar.alpha);
    HDF5::Write_Real_Attribute(file, "SPH speed of sound", fvar.Cs);
    HDF5::Write_Uint_Attribute(file, "SPH Newmark Beta iteration limit", svar.max_subits);
    HDF5::Write_Vector_Attribute(file, "SPH gravity vector", svar.grav);

    HDF5::Write_Real_Attribute(file, "SPH initial spacing", svar.particle_step);
    HDF5::Write_Real_Attribute(file, "SPH boundary spacing factor", svar.bound_step_factor);
    HDF5::Write_Real_Attribute(file, "SPH smoothing length factor", fvar.Hfac);
    HDF5::Write_String_Attribute(file, "SPH aerodynamic case", avar.aero_case);
    HDF5::Write_Int_Attribute(file, "SPH SP diameter definition", avar.use_dx);
    HDF5::Write_Int_Attribute(file, "SPH use TAB deformation", avar.useDef);
    HDF5::Write_Vector_Attribute(file, "SPH global offset coordinate", svar.offset_vec);
    HDF5::Write_Uint_Attribute(file, "SPH maximum particle count", svar.max_points);
    HDF5::Write_Real_Attribute(file, "SPH aerodynamic cutoff value", avar.cutoff);
    HDF5::Write_Vector_Attribute(file, "SPH freestream velocity", avar.vInf);
    HDF5::Write_Real_Attribute(file, "SPH restart fit tolerance", svar.restart_tol);

    /* Particle tracking settings */
    HDF5::Write_Int_Attribute(file, "Transition to IPT", svar.using_ipt);
    HDF5::Write_Int_Attribute(file, "Velocity equation order", svar.ipt_eq_order);
    HDF5::Write_Real_Attribute(file, "SPH tracking conversion x coordinate", svar.max_x_sph);
    HDF5::Write_Real_Attribute(file, "Maximum x trajectory coordinate", svar.max_x);
    HDF5::Write_Uint_Attribute(file, "Particle scatter output", svar.part_out);
    HDF5::Write_Uint_Attribute(file, "Particle streak output", svar.streak_out);
    HDF5::Write_Uint_Attribute(file, "Particle cell intersection output", svar.cells_out);
}

void Write_Zone_Attributes(int64_t const& zone, bound_block const& limits)
{
    HDF5::Write_String_Attribute(zone, "Block name", limits.name);
    // Add information about the insertion, aero, and deletion planes.
    HDF5::Write_Real_Attribute(zone, "Insertion normal x", limits.insert_norm[0]);
    HDF5::Write_Real_Attribute(zone, "Insertion normal y", limits.insert_norm[1]);
    HDF5::Write_Real_Attribute(zone, "Deletion normal x", limits.delete_norm[0]);
    HDF5::Write_Real_Attribute(zone, "Deletion normal y", limits.delete_norm[1]);
    HDF5::Write_Real_Attribute(zone, "Aerodynamic normal x", limits.aero_norm[0]);
    HDF5::Write_Real_Attribute(zone, "Aerodynamic normal y", limits.aero_norm[1]);
#if SIMDIM == 3
    HDF5::Write_Real_Attribute(zone, "Insertion normal z", limits.insert_norm[2]);
    HDF5::Write_Real_Attribute(zone, "Deletion normal z", limits.delete_norm[2]);
    HDF5::Write_Real_Attribute(zone, "Aerodynamic normal z", limits.aero_norm[2]);
#endif

    HDF5::Write_Real_Attribute(zone, "Insertion plane constant", limits.insconst);
    HDF5::Write_Real_Attribute(zone, "Deletion plane constant", limits.delconst);
    HDF5::Write_Real_Attribute(zone, "Aerodynamic plane constant", limits.aeroconst);

    // Add boundary time information
    if (limits.nTimes != 0)
    {
        HDF5::Write_Real_Attribute(zone, "Times count", limits.nTimes);
        for (size_t time = 0; time < limits.nTimes; time++)
        {
            HDF5::Write_Real_Attribute(zone, "Time " + std::to_string(time), limits.times[time]);
            HDF5::Write_Real_Attribute(
                zone, "Velocity " + std::to_string(time) + " x", limits.vels[time][0]
            );
            HDF5::Write_Real_Attribute(
                zone, "Velocity " + std::to_string(time) + " y", limits.vels[time][1]
            );
#if SIMDIM == 3
            HDF5::Write_Real_Attribute(
                zone, "Velocity " + std::to_string(time) + " z", limits.vels[time][2]
            );
#endif
        }
    }
    else
    {
        HDF5::Write_Real_Attribute(zone, "Times count", 0);
        HDF5::Write_Real_Attribute(zone, "Velocity x", limits.vels[0][0]);
        HDF5::Write_Real_Attribute(zone, "Velocity y", limits.vels[0][1]);
#if SIMDIM == 3
        HDF5::Write_Real_Attribute(zone, "Velocity z", limits.vels[0][2]);
#endif
    }

    HDF5::Write_Real_Attribute(zone, "Fixed velocity or dynamic inlet", limits.fixed_vel_or_dynamic);
    HDF5::Write_Real_Attribute(zone, "Lattice or HCPL packing", limits.particle_order);
    HDF5::Write_Real_Attribute(zone, "Boundary solver type", limits.bound_solver);
    HDF5::Write_Real_Attribute(zone, "Boundary is no slip", limits.no_slip);
    HDF5::Write_Real_Attribute(zone, "Block type", limits.block_type);
}

/**
 * @brief Initialise and add the boundary particles to the SPH particles
 *
 * @param svar Structure of simulation variables
 * @param parts SoA containing SPH particles
 * @param gas SoA containing gas particles and settings
 * @param fibvar SoA containing fibre particles and settings
 * @param act_time Current time for the simulation
 * @param nFluid number of fluid blocks
 * @param nBound number of boundary blocks
 * @param nFib number of fibre blocks
 * @param nGas number of gas blocks
 * @return void
 */
void Write_HDF5(
    SIM& svar, FLUID const& fvar, AERO const& avar, SPHState const& pnp1, LIMITS const& limits
)
{
    // printf("Starting writing restart HDF5 file...\n");
    std::string filename = svar.output_prefix + "_particles.h5";

    /*
     * Create a new file. If file exists its contents will be overwritten.
     */
    if (std::filesystem::exists(filename))
    { // Remove the file
        std::filesystem::remove(filename);
    }

    int64_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // int64_t group = H5Gcreate2(file,"Particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Write_HDF5_Attributes(file, svar, fvar, avar);

    /* Write the actual data */
    size_t blockID = 0;

    // Boundary particles
    if (svar.n_bound_blocks > 0)
    {
        for (size_t ii = 0; ii < svar.n_bound_blocks; ++ii)
        {
            std::string zoneHeader = "Boundary " + std::to_string(ii);

            // Create the group
            int64_t bdrzone =
                H5Gcreate2(file, (zoneHeader).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            HDF5::Write_Zone_Data(
                bdrzone, svar.scale, pnp1, limits[ii].index.first, limits[ii].index.second
            );

            // Write attribute data
            Write_Zone_Attributes(bdrzone, limits[ii]);

            if (H5Gclose(bdrzone))
            {
                printf("Failed to close HDF group\n");
                exit(-1);
            }

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
            int64_t fluzone =
                H5Gcreate2(file, (zoneHeader).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            HDF5::Write_Zone_Data(
                fluzone, svar.scale, pnp1, limits[ii].index.first, limits[ii].index.second
            );

            if (limits[ii].block_type == inletZone)
                HDF5::Write_Inlet_Data(fluzone, limits[ii]);

            // Need to write auxiliary data for the zone
            Write_Zone_Attributes(fluzone, limits[ii]);

            if (H5Gclose(fluzone))
            {
                printf("Failed to close HDF group\n");
                exit(-1);
            }

            blockID++;
            flublock++;
        }
    }

    // Close the files
    if (H5Fclose(file))
    {
        printf("Failed to close HDF file\n");
        exit(-1);
    }

    // If the output_prefix differs from the restart prefix, write a new restart prefix using that
    if (svar.restart && svar.output_prefix == svar.restart_prefix)
    {
        svar.restart_header_written = 1;
    }

    if (!svar.restart_header_written)
    {
        /* Append a restart prefix to the end of the para file */
        FILE* para = fopen(svar.input_file.c_str(), "a");
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
        fprintf(para, "SPH restart prefix: %s\n", svar.output_prefix.c_str());
        fclose(para);
        svar.restart_header_written = 1;
    }
    printf("Finished writing restart HDF5 file.\n");
}

void Read_HDF5(
    SIM& svar, FLUID& fvar, AERO& avar, VLM& vortex, SPHState& pn, SPHState& pnp1, LIMITS& limits
)
{
    printf("Starting reading restart HDF5 file...\n");

    std::string filename = svar.restart_prefix + "_particles.h5";
    int64_t file;
    if (std::filesystem::exists(filename))
        file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    else
    {
        printf("ERROR: Attempting to read file that doesn't exist: %s\n", filename.c_str());
        exit(-1);
    }

    // Read attributes, but only essential ones to not cause discrepancy
    HDF5::Read_Real_Attribute(file, "Simulation current time", svar.current_time);
    HDF5::Read_Uint_Attribute(file, "Current frame", svar.current_frame);
    HDF5::Read_Real_Attribute(file, "SPH initial spacing", svar.particle_step);
    HDF5::Read_Uint_Attribute(file, "Single file for output", svar.single_file);
    HDF5::Read_Real_Attribute(file, "SPH frame time interval", svar.frame_time_interval);
    // HDF5::Read_Real_Attribute(file,   "SPH previous frame time", svar.last_frame_time);

    HDF5::Read_String_Attribute(file, "Variable list", svar.output_names);
    HDF5::Read_Real_Attribute(file, "Grid scale", svar.scale);
    HDF5::Read_Real_Attribute(file, "SPH CFL condition", svar.cfl);

    HDF5::Read_Uint_Attribute(file, "Number of boundary blocks", svar.n_bound_blocks);
    HDF5::Read_Uint_Attribute(file, "Number of fluid blocks", svar.n_fluid_blocks);
    HDF5::Read_Uint_Attribute(file, "Number of boundary points", svar.bound_points);
    HDF5::Read_Uint_Attribute(file, "Number of fluid points", svar.fluid_points);
    HDF5::Read_Uint_Attribute(file, "Current frame", svar.current_frame);
    HDF5::Read_Uint_Attribute(file, "Particle index to add", svar.part_id);
    HDF5::Read_Uint_Attribute(file, "SPH stable CFL count", svar.n_stable);
    HDF5::Read_Uint_Attribute(file, "SPH unstable CFL count", svar.n_unstable);
    svar.last_frame_time = svar.current_time;

    Set_Values(svar, fvar, avar, vortex);

    size_t start = 0;
    size_t block = 0;
    if (svar.n_bound_blocks > 0)
    {
        // Will already be set from reading the boundary file
        // parts.ntimes.assign(nBound,0);
        // parts.b_times.resize(nBound);
        // parts.b_vels.resize(nBound);
        svar.bound_points = 0;
        for (size_t ii = 0; ii < svar.n_bound_blocks; ii++)
        {
            std::string zoneHeader = "Boundary " + std::to_string(ii);
            int64_t group = H5Gopen(file, zoneHeader.c_str(), H5P_DEFAULT);
            // size_t n0 = parts.grdnp1.size();
            limits[block].index.first = start;
            size_t length = 0;
            HDF5::Read_Zone_Data(group, svar.scale, pn, length);
            start += length;
            limits[block].index.second = start;

            // std::string name;
            // HDF5::Read_String_Attribute(group, "Block name", name);

            H5Gclose(group);
            block++;
            svar.bound_points += length;
        }
    }

    if (svar.n_fluid_blocks > 0)
    {
        svar.fluid_points = 0;
        for (size_t ii = 0; ii < svar.n_fluid_blocks; ii++)
        {
            // Assume the block definition files have been read.
            std::string zoneHeader = "Fluid " + std::to_string(ii);
            int64_t group = H5Gopen(file, zoneHeader.c_str(), H5P_DEFAULT);
            // size_t n0 = parts.grdnp1.size();

            limits[block].index.first = start;
            size_t length = 0;
            HDF5::Read_Zone_Data(group, svar.scale, pn, length);
            start += length;
            limits[block].index.second = start;

            // read the inlet data
            if (limits[block].block_type == inletZone)
                HDF5::Read_Inlet_Data(group, limits[block]);

            // Read any attributes?
            // std::string name;
            // HDF5::Read_String_Attribute(group, "Block name", name);
            // HDF5::Read_String_Attribute(group, "Block type", name);

            H5Gclose(group);
            block++;
            svar.fluid_points += length;
        }
    }

    H5Fclose(file);
    svar.total_points = svar.bound_points + svar.fluid_points;
    pnp1.resize(pn.size());
    copy_omp(pn.begin(), pn.end(), pnp1.begin());

    printf("Finished reading restart HDF5 file.\n");
}

namespace h5part
{

    void Write_Zone_Data(
        int64_t const& fout, real const& scale, double const& time, SPHState const& pnp1,
        size_t const& start, size_t const& end, OutputMap const& output_variables, double const& rho0
    )
    {
        /* Need to write position, velocity, acceleration, pressure, density gradient, mass, boundary
         * flag */
        /* Write as a 2D array, or 1D arrays for vectors? 1D for now */
        size_t length = end - start;
        hsize_t dims[] = {length};

        HDF5::Write_Real_Attribute(fout, "TimeValue", time);
        std::vector<double> vec(length);

        /* Position */
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].xi(dim) / scale;

            string name = "" + HDF5::get_dim(dim);
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
        }
#if SIMDIM == 2
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = 0;

        string name = "z";
        HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif

/* Density */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].rho;

        HDF5::Write_Variable_Scalar(fout, "Density", dims, 0, vec);

/* Mass */
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].m;

        HDF5::Write_Variable_Scalar(fout, "Mass", dims, 0, vec);

        /* Kinda optional outputs */

        /* Velocity */
        if (output_variables.at("vel-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].v(dim);

                string name = "Vel " + HDF5::get_dim(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "Vel z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }

        /* Acceleration */
        if (output_variables.at("acc-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].acc(dim);

                string name = "Acc " + HDF5::get_dim(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "Acc z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }

        /* Pressure */
        if (output_variables.at("press").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].p;
            HDF5::Write_Variable_Scalar(fout, "Pressure", dims, 0, vec);
        }

        if (output_variables.at("dRho").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].Rrho;
            HDF5::Write_Variable_Scalar(fout, "Density gradient", dims, 0, vec);
        }

        if (output_variables.at("part_id").write)
        {
            vector<int> uvec(length);
            for (size_t ii = start; ii < end; ++ii)
                uvec[ii - start] = pnp1[ii].part_id;
            HDF5::Write_Variable_Scalar(fout, "Particle ID", dims, 0, uvec);
        }

        if (output_variables.at("cellID").write)
        {
            vector<int> uvec(length);
            for (size_t ii = start; ii < end; ++ii)
                uvec[ii - start] = pnp1[ii].cellID;
            HDF5::Write_Variable_Scalar(fout, "Cell ID", dims, 0, uvec);
        }

        if (output_variables.at("bound").write)
        {
            vector<uint> uvec(length);
            for (size_t ii = start; ii < end; ++ii)
                uvec[ii - start] = pnp1[ii].b;
            HDF5::Write_Variable_Scalar(fout, "Particle status condition", dims, 0, uvec);
        }

        // Non essential variables, but to provide further information
        if (output_variables.at("dens").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].rho;
            HDF5::Write_Variable_Scalar(fout, "Density", dims, 0, vec);
        }

        if (output_variables.at("densVar").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 100.0 * (pnp1[ii].rho / rho0 - 1.0);
            HDF5::Write_Variable_Scalar(fout, "Density Variation", dims, 0, vec);
        }

        if (output_variables.at("vmag").write)
        {
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].v.norm();
            HDF5::Write_Variable_Scalar(fout, "Velocity magnitude", dims, 0, vec);
        }

        if (output_variables.at("surf").write)
        {
            vector<uint> uvec(length);
            for (size_t ii = start; ii < end; ++ii)
                uvec[ii - start] = pnp1[ii].surf;
            HDF5::Write_Variable_Scalar(fout, "Surface flag", dims, 0, uvec);
        }

        if (output_variables.at("surfZ").write)
        {
            vector<uint> uvec(length);
            for (size_t ii = start; ii < end; ++ii)
                uvec[ii - start] = pnp1[ii].surfzone;
            HDF5::Write_Variable_Scalar(fout, "Surface zone flag", dims, 0, uvec);
        }

        if (output_variables.at("aero-mag").write)
        {
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].Af.norm();
            HDF5::Write_Variable_Scalar(fout, "Aerodynamic force magnitude", dims, 0, vec);
        }

        if (output_variables.at("aero-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].Af(dim);

                string name = "Aerodynamic force " + std::to_string(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "Aerodynamic force z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }

        if (output_variables.at("curv").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].curve;
            HDF5::Write_Variable_Scalar(fout, "Curvature", dims, 0, vec);
        }

        if (output_variables.at("occl").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].woccl;
            HDF5::Write_Variable_Scalar(fout, "Occlusion factor", dims, 0, vec);
        }

        if (output_variables.at("cellP").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellP;
            HDF5::Write_Variable_Scalar(fout, "Cell pressure", dims, 0, vec);
        }

        if (output_variables.at("cellRho").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellP;
            HDF5::Write_Variable_Scalar(fout, "Cell density", dims, 0, vec);
        }

        if (output_variables.at("cellV-mag").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellV.norm();
            HDF5::Write_Variable_Scalar(fout, "Cell velocity magnitude", dims, 0, vec);
        }

        if (output_variables.at("cellV-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].cellV(dim);

                string name = "Cell velocity " + std::to_string(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "Cell velocity z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }

        if (output_variables.at("dsphG-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].gradRho(dim);

                string name = "dSPH density gradient" + std::to_string(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "dSPH densit gradient z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }

        if (output_variables.at("lam").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].lam;
            HDF5::Write_Variable_Scalar(fout, "Lambda eigenvalue", dims, 0, vec);
        }

        if (output_variables.at("lam-nb").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].lam_nb;
            HDF5::Write_Variable_Scalar(fout, "Lambda eigenvalue without boundary", dims, 0, vec);
        }

        if (output_variables.at("colour").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].colour;
            HDF5::Write_Variable_Scalar(fout, "Colour function", dims, 0, vec);
        }

        if (output_variables.at("colour-G").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].colourG;
            HDF5::Write_Variable_Scalar(fout, "Colour function gradient", dims, 0, vec);
        }

        if (output_variables.at("norm-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].norm(dim);

                string name = "Surface normal" + std::to_string(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "Surface normal z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }

        if (output_variables.at("shiftV-mag").write)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].vPert.norm();
            HDF5::Write_Variable_Scalar(fout, "Shifting velocity magnitude", dims, 0, vec);
        }

        if (output_variables.at("shiftV-vec").write)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (size_t ii = start; ii < end; ++ii)
                    vec[ii - start] = pnp1[ii].vPert(dim);

                string name = "Shifting velocity " + std::to_string(dim);
                HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
            }
#if SIMDIM == 2
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = 0;

            string name = "Shifting velocity z";
            HDF5::Write_Variable_Scalar(fout, name, dims, 0, vec);
#endif
        }
    }
} // namespace h5part

void open_h5part_files(
    SIM const& svar, FLUID const& fvar, AERO const& avar, string const& prefix,
    int64_t& h5part_fluid_file, int64_t& h5part_bound_file
)
{

    if (svar.fluid_points > 0)
    {
        std::string filename = prefix + "_fluid.h5part";

        if (std::filesystem::exists(filename))
        {
            h5part_fluid_file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        else
        {
            h5part_fluid_file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            Write_HDF5_Attributes(h5part_fluid_file, svar, fvar, avar);
        }
    }

    if (svar.bound_points > 0)
    {
        std::string filename = prefix + "_boundary.h5part";

        if (std::filesystem::exists(filename))
        {
            h5part_bound_file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        else
        {
            h5part_bound_file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            Write_HDF5_Attributes(h5part_bound_file, svar, fvar, avar);
        }
    }
}

/**
 * @brief Initialise and add the boundary particles to the SPH particles
 *
 * @param h5part_fluid_file fluid file identifier
 * @param h5part_bound_file boundary file identifier
 * @param sfile structure file identifier
 * @param svar Structure of simulation variables
 * @param parts SoA containing SPH particles
 * @param gas SoA containing gas particles and settings
 * @param fibvar SoA containing fibre particles and settings
 * @param act_time Current time for the simulation
 * @return void
 */
void write_h5part_data(
    int64_t& h5part_fluid_file, int64_t& h5part_bound_file, SIM const& svar, FLUID const& fvar,
    SPHState const& pnp1
)
{
    std::string zoneHeader = "Step#" + std::to_string(svar.current_frame);
    if (svar.n_fluid_blocks > 0)
    {
        // Create a group of the current timestep
        int64_t fluzone =
            H5Gcreate2(h5part_fluid_file, zoneHeader.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        h5part::Write_Zone_Data(
            fluzone, svar.scale, svar.current_time, pnp1, svar.bound_points, svar.total_points,
            svar.output_variables, fvar.rho0
        );

        if (H5Gclose(fluzone))
        {
            printf("Failed to close fluid HDF group\n");
            exit(-1);
        }

        // Close the file
        if (H5Fclose(h5part_fluid_file))
        {
            printf("Failed to close fluid HDF file\n");
            exit(-1);
        }
    }

    if (svar.n_bound_blocks > 0)
    {
        // Create a group of the current timestep
        int64_t bndzone =
            H5Gcreate2(h5part_bound_file, zoneHeader.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        h5part::Write_Zone_Data(
            bndzone, svar.scale, svar.current_time, pnp1, 0, svar.bound_points, svar.output_variables,
            fvar.rho0
        );

        if (H5Gclose(bndzone))
        {
            printf("Failed to close boundary HDF group\n");
            exit(-1);
        }

        // Close the file
        if (H5Fclose(h5part_bound_file))
        {
            printf("Failed to close boundary HDF file\n");
            exit(-1);
        }
    }
}
