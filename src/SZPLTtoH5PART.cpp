/*******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/*******          Created by Jamie MacLeod, University of Bristol        ***********/
/*******         Convert FJSPH Tecplot szplt files to h5part files       ***********/
/*******          Currently the other way is not supported (Sorry)       ***********/

/*******      Usage: FJSPH_conv infile.szplt -o outfile.h5part           ***********/
/*******      If no output file is specified, the input filename is      ***********/
/*******      used with the extension replaced.                          ***********/

#include "TECIO.h"
#include <hdf5.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <vector>

inline std::string ltrim(const std::string& s)
{
    size_t start = s.find_first_not_of(" \n\r\t\f\v");
    return (start == std::string::npos) ? "" : s.substr(start);
}

inline std::string rtrim(const std::string& s)
{
    size_t end = s.find_last_not_of(" \n\r\t\f\v");
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

struct ProgressBar
{
    ProgressBar() : percent(0), bar_length(40)
    { // Start with an empty bar
        std::string spaces = std::string((bar_length), ' ');
        printf("\rPercent: [%s] %d%%", spaces.c_str(), int(round(percent * 100)));
        fflush(stdout);
    };

    void reset()
    {
        percent = 0;
        std::string spaces = std::string((bar_length), ' ');
        printf("\rPercent: [%s] %d%%", spaces.c_str(), int(round(percent * 100)));
        fflush(stdout);
    }

    void update(float percent_)
    {
        int new_percent = round(percent_ * 100);
        if (new_percent > round(percent * 100))
        {
            percent = percent_;
            int bar_fill = round(percent * bar_length);
            std::string hashes = std::string(bar_fill, '#');
            std::string spaces = std::string((bar_length - hashes.size()), ' ');
            printf("\rPercent: [%s%s] %d%%", hashes.c_str(), spaces.c_str(), int(round(percent * 100)));
            fflush(stdout);
            // std::cout << "\rPercent: [" << hashes << spaces << "] " << int(round(percent * 100))
            //     << "%" << std::flush;
        }
    }

    float percent;
    int bar_length;
};

namespace HDF5
{
    inline void
    Write_Attribute(int64_t const& fout, std::string const& attr_name, size_t const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate2(fout, (attr_name).c_str(), H5T_NATIVE_ULONG, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_ULONG, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    inline void Write_Attribute(int64_t const& fout, std::string const& attr_name, uint const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate2(fout, (attr_name).c_str(), H5T_NATIVE_UINT, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_UINT, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    inline void Write_Attribute(int64_t const& fout, std::string const& attr_name, int const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate2(fout, (attr_name).c_str(), H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_INT, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    inline void
    Write_Attribute(int64_t const& fout, std::string const& attr_name, double const& attr_data)
    {
        herr_t retval;
        int64_t aid = H5Screate(H5S_SCALAR);
        int64_t attr =
            H5Acreate2(fout, (attr_name).c_str(), H5T_IEEE_F64LE, aid, H5P_DEFAULT, H5P_DEFAULT);
        if ((retval = H5Awrite(attr, H5T_NATIVE_DOUBLE, &attr_data)))
        {
            printf("\n\tError: Failed to write attribute \"%s\"\n", attr_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(aid);
        retval = H5Aclose(attr);
    }

    inline void
    Write_Attribute(int64_t const& fout, std::string const& attr_name, std::string const& attr_data)
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
                H5Acreate2(fout, attr_name.c_str(), attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
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

    inline void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<int>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate2(
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

    inline void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<int16_t>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate2(
            fout, var_name.c_str(), H5T_NATIVE_SHORT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start][0])
            ))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    inline void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<uint8_t>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate2(
            fout, var_name.c_str(), H5T_NATIVE_UCHAR, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start][0])
            ))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    inline void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<float>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate2(
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

    inline void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<std::vector<double>> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate_simple(2, dims, NULL);
        int64_t dataset = H5Dcreate2(
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
    inline void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<size_t> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate2(
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

    inline void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<int> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate2(
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

    inline void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<int16_t> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate2(
            fout, var_name.c_str(), H5T_NATIVE_SHORT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    inline void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<uint8_t> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate2(
            fout, var_name.c_str(), H5T_NATIVE_UCHAR, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        if ((retval =
                 H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[start])))
        {
            printf("\n\tError: Failed to write \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(dataspace);
        retval = H5Dclose(dataset);
    }

    inline void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<float> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate2(
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

    inline void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims, size_t const& start,
        std::vector<double> const& var_data
    )
    {
        herr_t retval;
        int64_t dataspace = H5Screate(H5S_SIMPLE);
        retval = H5Sset_extent_simple(dataspace, 1, dims, NULL);
        int64_t dataset = H5Dcreate2(
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

    /*************************************************************************/
    /*************************** READ FUNCTIONS ******************************/
    /*************************************************************************/
    inline void Read_Attribute(int64_t const& fin, std::string const& attr_name, std::string& attr_data)
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

    inline void Read_Attribute(int64_t const& fin, std::string const& attr_name, size_t& attr_data)
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

    inline void Read_Attribute(int64_t const& fin, std::string const& attr_name, uint& attr_data)
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

    inline void Read_Attribute(int64_t const& fin, std::string const& attr_name, int& attr_data)
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

    inline void Read_Attribute(int64_t const& fin, std::string const& attr_name, double& attr_data)
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

    inline void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<size_t>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("\nWARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<size_t>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<int>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("\nWARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<int>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<float>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("\nWARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<float>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void Read_Variable_Array(
        int64_t const& fin, std::string const& var_name, std::vector<std::vector<double>>& var_data
    )
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[2];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 2)
            printf("\nWARNING: Not expected dimensions to array\n");

        var_data.resize(dims[0], std::vector<double>(dims[1]));
        if ((retval = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0][0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }

        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<size_t>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("\nWARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void
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

    inline void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<int16_t>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("\nWARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<uint8_t>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("\nWARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<float>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("\nWARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }

    inline void
    Read_Variable_Scalar(int64_t const& fin, std::string const& var_name, std::vector<double>& var_data)
    {
        herr_t retval;
        int64_t dset = H5Dopen(fin, var_name.c_str(), H5P_DEFAULT);
        int64_t space = H5Dget_space(dset);
        hsize_t dims[1];
        int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
        if (ndims != 1)
            printf("\nWARNING: More than one dimension to array\n");

        var_data.resize(dims[0]);
        if ((retval = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_data[0])))
        {
            printf("\n\tError: Failed to read \"%s\" data\n", var_name.c_str());
            exit(-1);
        }
        retval = H5Sclose(space);
        retval = H5Dclose(dset);
    }
} // namespace HDF5

void trans_var_tec_to_hdf5(
    void* const& inFile, int32_t const& inputZone, std::string const& varname,
    std::vector<int32_t> const& varTypes, int32_t const& var, int64_t const& hFile
)
{
    int64_t numValues;
    int32_t i = tecZoneVarGetNumValues(inFile, inputZone, var, &numValues);
    hsize_t dim[] = {static_cast<hsize_t>(numValues)};

    // For large zones, could "chunk" this input/output--read/write the var in pieces instead of all at
    // once
    switch ((FieldDataType_e)varTypes[var - 1])
    {
    case FieldDataType_Float:
    {
        std::vector<float> values(numValues);
        i = tecZoneVarGetFloatValues(inFile, inputZone, var, 1, numValues, &values[0]);
        HDF5::Write_Variable_Scalar(hFile, varname, dim, 0, values);
        break;
    }
    case FieldDataType_Double:
    {
        std::vector<double> values(numValues);
        i = tecZoneVarGetDoubleValues(inFile, inputZone, var, 1, numValues, &values[0]);
        HDF5::Write_Variable_Scalar(hFile, varname, dim, 0, values);
        break;
    }
    case FieldDataType_Int32:
    {
        std::vector<int32_t> values(numValues);
        i = tecZoneVarGetInt32Values(inFile, inputZone, var, 1, numValues, &values[0]);
        HDF5::Write_Variable_Scalar(hFile, varname, dim, 0, values);
        break;
    }
    case FieldDataType_Int16:
    {
        std::vector<int16_t> values(numValues);
        i = tecZoneVarGetInt16Values(inFile, inputZone, var, 1, numValues, &values[0]);
        HDF5::Write_Variable_Scalar(hFile, varname, dim, 0, values);
        break;
    }
    case FieldDataType_Byte:
    {
        std::vector<uint8_t> values(numValues);
        i = tecZoneVarGetUInt8Values(inFile, inputZone, var, 1, numValues, &values[0]);
        HDF5::Write_Variable_Scalar(hFile, varname, dim, 0, values);
        break;
    }
    default:
        i = -1;
        break;
    }

    if (i)
    {
        printf("\nAn error occured transcribing a variable\n");
        exit(-1);
    }
}

void trans_zone_tec_to_hdf5(
    void* const& inFile, int32_t const& inputZone, int64_t const& hFile,
    std::vector<std::string> const& varnames, double const& rho_rest, double const& sos, double& dx
)
{
    int32_t i = 0;
    std::string zoneHeader = "Step#" + std::to_string(inputZone - 1);
    int64_t outputZone = H5Gcreate2(hFile, zoneHeader.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    char* zoneTitle = NULL;
    i = tecZoneGetTitle(inFile, inputZone, &zoneTitle);
    HDF5::Write_Attribute(outputZone, "Zone title", zoneTitle);
    tecStringFree(&zoneTitle);

    double solutionTime;
    i = tecZoneGetSolutionTime(inFile, inputZone, &solutionTime);
    HDF5::Write_Attribute(outputZone, "TimeValue", solutionTime);

    int32_t numVars;
    i = tecDataSetGetNumVars(inFile, &numVars);
    std::vector<int32_t> varTypes(numVars);
    for (int32_t var = 1; var <= numVars; ++var)
        i = tecZoneVarGetType(inFile, inputZone, var, &varTypes[var - 1]);

    for (int32_t var = 1; var <= numVars; ++var)
        trans_var_tec_to_hdf5(inFile, inputZone, varnames[var - 1], varTypes, var, outputZone);

    // Check if mess and density vectors exists
    int hasmass = 0;
    int hasdens = 0;
    int is3d = 0;
    int hasx(0), hasy(0), hasz(0);
    for (int32_t var = 0; var < numVars; ++var)
    {
        if (varnames[var] == "mass" || varnames[var] == "m")
            hasmass = 1;

        if (varnames[var] == "density" || varnames[var] == "rho")
            hasdens = 1;

        if (varnames[var] == "X")
            hasx = 1;

        if (varnames[var] == "Y")
            hasy = 1;

        if (varnames[var] == "Z")
            hasz = 1;
    }

    int offset = 0;
    if (hasx && hasy && hasz)
    {
        is3d = 1;
    }
    else
    {
        if (hasx && hasy)
            offset = 3;
        if (hasx && hasz)
            offset = 2;
        if (hasy && hasz)
            offset = 1;
    }

    if (is3d == 0)
    {
        // Need particle count
        int64_t numValues;
        i = tecZoneVarGetNumValues(inFile, inputZone, 1, &numValues);
        std::vector<double> val(numValues, 0.0);
        hsize_t dim[] = {static_cast<hsize_t>(numValues)};

        if (offset == 1)
            HDF5::Write_Variable_Scalar(outputZone, "X", dim, 0, val);
        else if (offset == 2)
            HDF5::Write_Variable_Scalar(outputZone, "Y", dim, 0, val);
        else if (offset == 3)
            HDF5::Write_Variable_Scalar(outputZone, "Z", dim, 0, val);

        for (int32_t var = 0; var < numVars; ++var)
        {
            std::string temp = varnames[var].substr(0, varnames[var].length() - 1);
            // Add third dimension of vectors

            if (temp == "A" || temp == "V" || temp == "Af" || temp == "cellV" || temp == "dsphG" ||
                temp == "surf-norm" || temp == "shiftV" || temp == "a_" || temp == "v_")
            {
                if (offset == 1)
                    HDF5::Write_Variable_Scalar(outputZone, temp + "x", dim, 0, val);
                else if (offset == 2)
                    HDF5::Write_Variable_Scalar(outputZone, temp + "y", dim, 0, val);
                else if (offset == 3)
                    HDF5::Write_Variable_Scalar(outputZone, temp + "z", dim, 0, val);
                var++;
            }
        }
    }

    // Zone auxiliary data
    double mass = 0.0;
    int32_t numItems;
    i = tecZoneAuxDataGetNumItems(inFile, inputZone, &numItems);
    for (int32_t whichItem = 1; whichItem <= numItems; ++whichItem)
    {
        char* name = NULL;
        char* value = NULL;
        i = tecZoneAuxDataGetItem(inFile, inputZone, whichItem, &name, &value);
        if (!H5Aexists(outputZone, name))
            HDF5::Write_Attribute(outputZone, std::string(name), std::string(value));

        if (strcmp(name, "Particle_mass") == 0)
        {
            std::istringstream iss(value);
            iss >> mass;
        }

        tecStringFree(&name);
        tecStringFree(&value);
    }

    if (hasmass == 0)
    {
        // Need particle count
        int64_t numValues;
        i = tecZoneVarGetNumValues(inFile, inputZone, 1, &numValues);
        if (mass == 0.0)
            mass = rho_rest * pow(dx, 2 + is3d);

        if (mass == 0.0)
        {
            printf("\nERROR: Mass is zero, so no information was found.\n");
            exit(-1);
        }
        std::vector<double> m(numValues, mass);
        hsize_t dim[] = {static_cast<hsize_t>(numValues)};
        HDF5::Write_Variable_Scalar(outputZone, "Mass", dim, 0, m);
    }

    if (hasdens == 0)
    {
        // Need particle count
        int64_t numValues;
        i = tecZoneVarGetNumValues(inFile, inputZone, 1, &numValues);
        hsize_t dim[] = {static_cast<hsize_t>(numValues)};

        // Retrieve the pressure values
        std::vector<double> p(numValues);
        for (int32_t var = 1; var <= numVars; ++var)
            if (varnames[var - 1] == "Pressure")
                i = tecZoneVarGetDoubleValues(inFile, inputZone, var, 1, numValues, &p[0]);

        std::vector<double> den(numValues, 0.0);
        double Bconst(rho_rest * sos * sos / 7.0);
        for (int64_t ii = 0; ii < numValues; ii++)
            den[ii] = pow((p[ii] / Bconst + 1.0), 1.0 / 7.0) * rho_rest;

        if (rho_rest == 0.0)
        {
            printf("\nERROR: Density is zero, so no information was found.\n");
            exit(-1);
        }

        HDF5::Write_Variable_Scalar(outputZone, "Density", dim, 0, den);
    }

    if (H5Gclose(outputZone))
    {
        printf("\nFailed to close HDF group\n");
        exit(-1);
    }

    if (i)
    {
        printf("An error occured transcribing a zone\n");
        exit(-1);
    }
}

void trans_file_tec_to_hdf5(
    void* const& inFile, int64_t const& hFile, double& rho_rest, double& sos, double& dx
)
{
    int32_t i = 0;

    // double rho_rest = 0.0;
    // double sos = 0.0;

    // Write aux data
    int32_t numItems;
    i = tecDataSetAuxDataGetNumItems(inFile, &numItems);
    for (int32_t whichItem = 1; whichItem <= numItems; ++whichItem)
    {
        char* name = NULL;
        char* value = NULL;
        i = tecDataSetAuxDataGetItem(inFile, whichItem, &name, &value);
        if (!H5Aexists(hFile, name))
            HDF5::Write_Attribute(hFile, std::string(name), std::string(value));

        if (strcmp(name, "Reference_dispersed_density") == 0)
        {
            std::istringstream iss(value);
            iss >> rho_rest;
        }

        if (strcmp(name, "SPH_speed_of_sound") == 0)
        {
            std::istringstream iss(value);
            iss >> sos;
        }

        tecStringFree(&name);
        tecStringFree(&value);
    }

    // Number of variables
    int32_t numVars;
    if (tecDataSetGetNumVars(inFile, &numVars))
    {
        printf("Failed to get number of variables.\n");
        exit(-1);
    }

    std::ostringstream outputStream;
    std::vector<std::string> varnames;
    for (int32_t var = 1; var <= numVars; ++var)
    {
        char* name = NULL;
        if (tecVarGetName(inFile, var, &name))
        {
            printf("Failed to get variable name.\n");
            exit(-1);
        }
        outputStream << name;
        varnames.emplace_back(name);
        if (var < numVars)
            outputStream << ',';
        tecStringFree(&name);
    }

    // Read the zones.
    int32_t numZones;
    if (tecDataSetGetNumZones(inFile, &numZones))
    {
        printf("Failed to get number of zones.\n");
        exit(-1);
    }

    ProgressBar progress;
    for (int32_t inputZone = 1; inputZone <= numZones; inputZone++)
    {
        trans_zone_tec_to_hdf5(inFile, inputZone, hFile, varnames, rho_rest, sos, dx);
        progress.update(float(inputZone - 1.0) / float(numZones));
    }

    // Close the file
    if (H5Fclose(hFile))
    {
        printf("\nFailed to close HDF file\n");
        exit(-1);
    }

    if (i)
    {
        printf("\nAn error occured transcribing a file\n");
        exit(-1);
    }
}

void convert_tecplot(
    std::string const& szfile, std::string const& hfile, double& rho_rest, double& sos, double& dx
)
{
    // Open the last file to get the variables and aux data
    void* inFile = NULL;
    if (tecFileReaderOpen(szfile.c_str(), &inFile))
    {
        printf("Failed to open %s\n", szfile.c_str());
        exit(-1);
    }

    // Number of variables
    int32_t numVars;
    if (tecDataSetGetNumVars(inFile, &numVars))
    {
        printf("Failed to get number of variables.\n");
        exit(-1);
    }

    std::ostringstream outputStream;
    for (int32_t var = 1; var <= numVars; ++var)
    {
        char* name = NULL;
        if (tecVarGetName(inFile, var, &name))
        {
            printf("Failed to get variable name.\n");
            exit(-1);
        }

        outputStream << name;
        if (var < numVars)
            outputStream << ',';
        tecStringFree(&name);
    }

    int64_t hFile = H5Fcreate(hfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    trans_file_tec_to_hdf5(inFile, hFile, rho_rest, sos, dx);
}

int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        printf("No inputs provdied. Stopping.");
        exit(0);
    }

    // int index;
    int c;

    opterr = 0;
    std::string outprefix;
    std::string infile;

    while ((c = getopt(argc, argv, "o:")) != -1)
        switch (c)
        {
        case 'o':
            outprefix = optarg;
            break;
        // case 'h':
        //     h5part = true;
        //     break;
        case '?':
            if (optopt == 'o')
                fprintf(stderr, "Option -%o requires an argument.\n", optopt);
            else if (isprint(optopt))
                fprintf(stderr, "Unknown option `-%o'.\n", optopt);
            else
                fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            return 1;
        default:
            abort();
        }

    infile = argv[optind];
    double rho_rest(0.0), sos(0.0), dx(0.0);
    if (argc - optind > 2)
    {
        std::istringstream iss(argv[optind + 1]);
        iss >> rho_rest;

        std::istringstream iss2(argv[optind + 2]);
        iss2 >> sos;

        std::istringstream iss3(argv[optind + 3]);
        iss3 >> dx;
    }
    else
    {
        printf(
            "WARNING: No info on density and spacing provided, going to try and find it in the file.\n"
        );
    }

    // Determining mass and density of particles...

    // for(index = optind; index < argc; index++)
    // {
    //     infile = argv[index];
    // }

    int szorh5 = 0;
    std::filesystem::path fent(infile);
    if (std::filesystem::exists(infile))
    {
        if (fent.extension() == ".szplt")
            szorh5 = 0;
        else if (fent.extension() == ".h5part")
            szorh5 = 1;
    }
    else
    {
        printf("Input file not found. Exiting.\n");
        exit(-1);
    }

    if (outprefix.empty())
        outprefix = fent.stem();

    outprefix = ltrim(rtrim(outprefix));
    std::string outfile = outprefix;
    if (szorh5 == 0)
        outfile += ".h5part";
    else
        outfile += ".szplt";

    printf(" Input file: %s\nOutput file: %s\n", infile.c_str(), outfile.c_str());

    if (szorh5 == 0)
        convert_tecplot(infile, outfile, rho_rest, sos, dx);

    printf("\nConversion complete!\n");

    return 0;
}