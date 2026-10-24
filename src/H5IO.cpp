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
    hid_t Open_File(std::string const& filename)
    {
        hid_t file_handle;
        if (std::filesystem::exists(filename))
        {
            file_handle = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

            if (file_handle == H5I_INVALID_HID)
            {
                cerr << "Failed to open new instance of HDF5 file: " << filename.c_str() << endl;
                exit(-1);
            }
        }
        else
        {
            file_handle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            if (file_handle == H5I_INVALID_HID)
            {
                cerr << "Failed to open existing instance of HDF5 file: " << filename.c_str() << endl;
                exit(-1);
            }
        }

        return file_handle;
    }

    hid_t Open_Zone(hid_t const& file_handle, string const& zone_header)
    {
        assert(file_handle != H5I_INVALID_HID);

        hid_t bdrzone =
            H5Gcreate2(file_handle, (zone_header).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    void Close_File(hid_t handle)
    {
        if (H5Fclose(handle))
        {
            printf("Failed to close fluid HDF file\n");
            exit(-1);
        }
    }

    void Close_Handle(hid_t handle)
    {
        assert(handle != H5I_INVALID_HID);

        if (H5Gclose(handle))
        {
            cerr << "Failed to close HDF group" << endl;
            exit(-1);
        }
    }

    // Do not call directly. Needs to be called from a typed attribute writer.
    template <typename T>
    void _Write_Attribute(
        hid_t const& fout, std::string const& attr_name, T const& attr_data, hid_t const& var_type,
        hid_t const& write_type
    )
    {
        herr_t retval;
        hid_t aid = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate(fout, (attr_name).c_str(), var_type, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, write_type, &attr_data)))
        {
            printf(
                "\n\tError: Failed to write attribute \"%s\" to handle %ld\n", attr_name.c_str(), fout
            );
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, size_t const& attr_data)
    {
        _Write_Attribute(fout, attr_name, attr_data, H5T_NATIVE_ULONG, H5T_NATIVE_ULONG);
    }

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, uint const& attr_data)
    {
        _Write_Attribute(fout, attr_name, attr_data, H5T_NATIVE_UINT, H5T_NATIVE_UINT);
    }

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, int const& attr_data)
    {
        _Write_Attribute(fout, attr_name, attr_data, H5T_NATIVE_UINT, H5T_NATIVE_UINT);
    }

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, double const& attr_data)
    {
        _Write_Attribute(fout, attr_name, attr_data, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE);
    }

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, StateVecD const& attr_data)
    {
        std::string name_ = attr_name /* + " x" */;
        hsize_t dims[1] = {SIMDIM};
        herr_t retval;
        hid_t aid = H5Screate_simple(1, dims, NULL);
        hid_t attr = H5Acreate(fout, (name_).c_str(), H5T_IEEE_F64LE, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_DOUBLE, &attr_data[0])))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, std::string const& attr_data)
    {
        if (attr_data.length() > 0)
        {
            char data[attr_data.length() + 1];
            strcpy(data, attr_data.c_str());

            hid_t attr_type = H5Tcopy(H5T_FORTRAN_S1);
            herr_t retval = H5Tset_size(attr_type, attr_data.length());
            hid_t memtype = H5Tcopy(H5T_C_S1);
            retval = H5Tset_size(memtype, attr_data.length() + 1);

            hid_t attr_space = H5Screate(H5S_SCALAR);

            hid_t attr =
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

    // Do not call directly
    template <typename T>
    void Write_Var_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<T>> const& var_data, hid_t const& var_type, hid_t const& write_type
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate_simple(2, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), var_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, write_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<size_t>> const& var_data
    )
    {
        Write_Var_Array(fout, var_name, dims, var_data, H5T_NATIVE_ULONG, H5T_NATIVE_ULONG);
    }

    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<int>> const& var_data
    )
    {
        Write_Var_Array(fout, var_name, dims, var_data, H5T_NATIVE_INT, H5T_NATIVE_INT);
    }

    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<float>> const& var_data
    )
    {
        Write_Var_Array(fout, var_name, dims, var_data, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT);
    }

    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<double>> const& var_data
    )
    {
        Write_Var_Array(fout, var_name, dims, var_data, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE);
    }

    // #if UINTPTR_MAX == 0xffffffffffffffff
    /* 64-bit pointers */
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<size_t> const& var_data
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }
    // #elif UINTPTR_MAX == 0xffffffff
    /* 32-bit pointers */
    // void Write_Variable_Scalar(hid_t const& fout, std::string const& var_name,
    //         hsize_t const* dims, size_t const& start, std::vector<size_t> const& var_data)
    // {
    //     herr_t retval;
    //     hid_t  dataspace = H5Screate(H5S_SIMPLE);
    //     retval = H5Sset_extent_simple(dataspace,1,dims,NULL);
    //     hid_t dataset = H5Dcreate2(fout, var_name.c_str(), H5T_NATIVE_UINT, dataspace, H5P_DEFAULT,
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

    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<uint> const& var_data
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<int> const& var_data
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<long> const& var_data
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_LONG, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<float> const& var_data
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<double> const& var_data
    )
    {
        herr_t retval;
        hid_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        hid_t dataset = H5Dcreate(
            fout, var_name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    void Write_Variable(
        hid_t const& fout, std::string const& name, hsize_t const* dims,
        std::vector<StateVecD> const& data
    )
    {
        size_t length = data.size();
        std::vector<double> vect(data.size());
        hsize_t dim[] = {data.size()};
        std::string name_;

#pragma omp parallel for
        for (size_t ii = 0; ii < length; ++ii)
            vect[ii] = data[ii][0];

        name_ = name.empty() ? "x" : name + " x";
        Write_Variable(fout, name_.c_str(), dims, vect);

#pragma omp parallel for
        for (size_t ii = 0; ii < length; ++ii)
            vect[ii] = data[ii][1];

        name_ = name.empty() ? "y" : name + " y";
        Write_Variable(fout, name_.c_str(), dims, vect);

        // Always need the z-dimension written for HDF5/h5part
        name_ = name.empty() ? "z" : name + " z";
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ++ii)
#if SIMDIM == 3
            vect[ii] = data[ii][2];
#else
            vect[ii] = 0.0;
#endif

        Write_Variable(fout, name_.c_str(), dims, vect);
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

    void Write_Inlet_Data(hid_t const& fout, bound_block const& limits)
    {
        size_t length = limits.back.size();
        hsize_t dim[] = {length};
        hsize_t dims[] = {length * 4};

        Write_Variable(fout, "Back index", dim, limits.back);

        vector<size_t> buff(length * 4);
        size_t accm(0);
        for (size_t ii = 0; ii < length; ii++)
            for (size_t jj = 0; jj < 4; jj++)
                buff[accm++] = limits.buffer[ii][jj];

        Write_Variable(fout, "Buffer index", dims, buff);
    }

    /*************************************************************************/
    /*************************** READ FUNCTIONS ******************************/
    /*************************************************************************/
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, std::string& attr_data)
    {
        herr_t retval;

        if (H5Aexists(fin, attr_name.c_str()))
        {
            hid_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
            hid_t attr_type = H5Aget_type(attr);

            size_t sdim = H5Tget_size(attr_type) + 1;
            hid_t memtype = H5Tcopy(H5T_C_S1);
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

    // Do not call directly.
    template <typename T>
    void
    _Read_Attribute(hid_t const& fin, std::string const& attr_name, T& attr_data, hid_t const& read_type)
    {
        herr_t retval;
        hid_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        size_t data;
        if ((retval = H5Aread(attr, read_type, &data)))
        {
            printf("\n\tError: Failed to read attribute: \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Aclose(attr);
        attr_data = data;
    }

    void Read_Attribute(hid_t const& fin, std::string const& attr_name, size_t& attr_data)
    {
        _Read_Attribute(fin, attr_name, attr_data, H5T_NATIVE_ULONG);
    }

    void Read_Attribute(hid_t const& fin, std::string const& attr_name, uint& attr_data)
    {
        _Read_Attribute(fin, attr_name, attr_data, H5T_NATIVE_UINT);
    }

    void Read_Attribute(hid_t const& fin, std::string const& attr_name, int& attr_data)
    {
        _Read_Attribute(fin, attr_name, attr_data, H5T_NATIVE_INT);
    }

    void Read_Attribute(hid_t const& fin, std::string const& attr_name, double& attr_data)
    {
        _Read_Attribute(fin, attr_name, attr_data, H5T_NATIVE_DOUBLE);
    }

    void Read_Attribute(hid_t const& fin, std::string const& attr_name, StateVecD& attr_data)
    {
        herr_t retval;
        hid_t attr = H5Aopen(fin, attr_name.c_str(), H5P_DEFAULT);
        hid_t space = H5Aget_space(attr);
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

    template <typename T>
    void Read_Var_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<T>>& var_data,
        hid_t const& read_type
    )
    {
        herr_t retval;
        hid_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("WARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<T>(dims[1]));
        if ((retval = H5Dread(dset, read_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<size_t>>& var_data
    )
    {
        Read_Var_Array(fin, var_name, var_data, H5T_NATIVE_ULONG);
    }

    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<int>>& var_data
    )
    {
        Read_Var_Array(fin, var_name, var_data, H5T_NATIVE_INT);
    }

    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<float>>& var_data
    )
    {
        Read_Var_Array(fin, var_name, var_data, H5T_NATIVE_FLOAT);
    }

    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<double>>& var_data
    )
    {
        Read_Var_Array(fin, var_name, var_data, H5T_NATIVE_DOUBLE);
    }

    // Do not call directly.
    template <typename T>
    void Read_Var_Scalar(
        hid_t const& fin, std::string const& var_name, std::vector<T>& var_data, hid_t const& read_type
    )
    {
        herr_t retval;
        hid_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("WARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, read_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<size_t>& var_data)
    {
        Read_Var_Scalar(fin, var_name, var_data, H5T_NATIVE_ULONG);
    }

    void Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<uint>& var_data)
    {
        Read_Var_Scalar(fin, var_name, var_data, H5T_NATIVE_UINT);
    }

    void Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<int>& var_data)
    {
        Read_Var_Scalar(fin, var_name, var_data, H5T_NATIVE_INT);
    }

    void Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<long>& var_data)
    {
        Read_Var_Scalar(fin, var_name, var_data, H5T_NATIVE_LONG);
    }

    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<float>& var_data)
    {
        Read_Var_Scalar(fin, var_name, var_data, H5T_NATIVE_FLOAT);
    }

    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<double>& var_data)
    {
        Read_Var_Scalar(fin, var_name, var_data, H5T_NATIVE_DOUBLE);
    }

    void Read_Vector_Data(hid_t const& fin, std::string const& name, std::vector<StateVecD>& var_data)
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

    void Read_Zone_Data(hid_t const& fin, real const& scale, SPHState& pnp1, size_t& length)
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

        /* Particle spacing */
        Read_Variable_Scalar(fin, "Particle spacing", dvect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].dx = dvect[ii];

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

        /* Fluid ID */
        Read_Variable_Scalar(fin, "Fluid ID", svect);
#pragma omp parallel for
        for (size_t ii = 0; ii < length; ii++)
            temp[ii].fluid_id = svect[ii];

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

    void Read_Inlet_Data(hid_t const& fin, bound_block& limits)
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

void Read_Zone_Attributes(hid_t const& zone, bound_block& limits)
{
    HDF5::Read_Attribute(zone, "Block name", limits.name);
    // Add information about the insertion, aero, and deletion planes.
    HDF5::Read_Attribute(zone, "Insertion normal x", limits.insert_norm[0]);
    HDF5::Read_Attribute(zone, "Insertion normal y", limits.insert_norm[1]);
    HDF5::Read_Attribute(zone, "Deletion normal x", limits.delete_norm[0]);
    HDF5::Read_Attribute(zone, "Deletion normal y", limits.delete_norm[1]);
    HDF5::Read_Attribute(zone, "Aerodynamic normal x", limits.aero_norm[0]);
    HDF5::Read_Attribute(zone, "Aerodynamic normal y", limits.aero_norm[1]);
#if SIMDIM == 3
    HDF5::Read_Attribute(zone, "Insertion normal z", limits.insert_norm[2]);
    HDF5::Read_Attribute(zone, "Deletion normal z", limits.delete_norm[2]);
    HDF5::Read_Attribute(zone, "Aerodynamic normal z", limits.aero_norm[2]);
#endif

    HDF5::Read_Attribute(zone, "Insertion plane constant", limits.insconst);
    HDF5::Read_Attribute(zone, "Deletion plane constant", limits.delconst);
    HDF5::Read_Attribute(zone, "Aerodynamic plane constant", limits.aeroconst);

    HDF5::Read_Attribute(zone, "Fixed velocity or dynamic inlet", limits.fixed_vel_or_dynamic);
    HDF5::Read_Attribute(zone, "Lattice or HCPL packing", limits.particle_order);
    HDF5::Read_Attribute(zone, "Boundary solver type", limits.bound_solver);
    HDF5::Read_Attribute(zone, "Boundary is no slip", limits.no_slip);
    HDF5::Read_Attribute(zone, "Block type", limits.block_type);

    // Fluid properties
    HDF5::Read_Attribute(zone, "Cole EOS gamma", limits.fluid.gam);
    HDF5::Read_Attribute(zone, "Reference dispersed density", limits.fluid.rho_rest);
    HDF5::Read_Attribute(zone, "Reference dispersed viscosity", limits.fluid.mu);
    HDF5::Read_Attribute(zone, "Reference surface tension", limits.fluid.sig);
    HDF5::Read_Attribute(zone, "Surface tension contact angle", limits.fluid.contangb);
    HDF5::Read_Attribute(zone, "Equation of state", limits.fluid.pressure_rel);
    HDF5::Read_Attribute(zone, "Background pressure", limits.fluid.press_back);
    HDF5::Read_Attribute(zone, "Starting pressure", limits.fluid.press_start);
    HDF5::Read_Attribute(zone, "SPH delta coefficient", limits.fluid.dsph_delta);
    HDF5::Read_Attribute(zone, "Artificial viscosity factor", limits.fluid.visc_alpha);
    HDF5::Read_Attribute(zone, "Speed of sound", limits.fluid.speed_sound);
}

/**
 * @brief Read a .h5 file to ingest restart data.
 *
 * @param svar Structure of simulation variables
 * @param pn SoA of containing SPH particles at time n (at time of reading, n = n+1, but need to
 * initialise the data)
 * @param pnp1 SoA containing SPH particles at time n+1
 * @param limits array of block indexes for particles in the global array.
 * @return void
 */
void Read_HDF5(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits)
{
    printf("Starting reading restart HDF5 file...\n");

    std::string filename = svar.io.restart_prefix + "_particles.h5";
    hid_t file;
    if (std::filesystem::exists(filename))
        file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    else
    {
        printf("ERROR: Attempting to read file that doesn't exist: %s\n", filename.c_str());
        exit(-1);
    }

    // Read attributes, but only essential ones to not cause discrepancy
    HDF5::Read_Attribute(file, "Simulation current time", svar.integrator.current_time);
    HDF5::Read_Attribute(file, "Current frame", svar.integrator.current_frame);
    HDF5::Read_Attribute(file, "Single file for output", svar.io.single_file);
    HDF5::Read_Attribute(file, "SPH frame time interval", svar.integrator.frame_time_interval);
    // HDF5::Read_Attribute(file,   "SPH previous frame time", svar.last_frame_time);

    HDF5::Read_Attribute(file, "Variable list", svar.io.output_names);
    HDF5::Read_Attribute(file, "Grid scale", svar.scale);
    HDF5::Read_Attribute(file, "SPH CFL condition", svar.integrator.cfl);

    HDF5::Read_Attribute(file, "Number of boundary blocks", svar.n_bound_blocks);
    HDF5::Read_Attribute(file, "Number of fluid blocks", svar.n_fluid_blocks);
    HDF5::Read_Attribute(file, "Number of boundary points", svar.bound_points);
    HDF5::Read_Attribute(file, "Number of fluid points", svar.fluid_points);
    HDF5::Read_Attribute(file, "Particle index to add", svar.part_id);
    HDF5::Read_Attribute(file, "SPH stable CFL count", svar.integrator.n_stable);
    HDF5::Read_Attribute(file, "SPH unstable CFL count", svar.integrator.n_unstable);
    svar.integrator.last_frame_time = svar.integrator.current_time;

    Set_Values(svar);

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
            hid_t group = H5Gopen(file, zoneHeader.c_str(), H5P_DEFAULT);
            // size_t n0 = parts.grdnp1.size();
            limits[block].index.first = start;
            size_t length = 0;
            HDF5::Read_Zone_Data(group, svar.scale, pn, length);
            start += length;
            limits[block].index.second = start;

            // std::string name;
            // HDF5::Read_Attribute(group, "Block name", name);

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
            hid_t group = H5Gopen(file, zoneHeader.c_str(), H5P_DEFAULT);
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
            // HDF5::Read_Attribute(group, "Block name", name);
            // HDF5::Read_Attribute(group, "Block type", name);

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
