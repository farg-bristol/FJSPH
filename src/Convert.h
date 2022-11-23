#ifndef CONVERT_H
#define CONVERT_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <netcdf.h>

#define NC_ERR 2
#define ERR(e)                                     \
{                                                  \
    printf("\n\tError: %s\n", nc_strerror(e)); \
    exit(-1);                                      \
}

using std::cout;
using std::endl;
using std::setw;
using std::string;
using std::vector;

/* Define data type. */
/* Define data type. */
#ifndef FOD
#define FOD 1 /*0 = float, 1 = double*/
#endif

#if FOD == 1
typedef double real;
#else
typedef float real;
#endif

typedef unsigned int uint;

#ifdef DEBUG
	/*Open debug file to write to*/
	std::ofstream dbout("Cell_Convert.log",std::ios::out);
#endif

uint index(uint ii, uint jj, uint nPts)
{
	return(ii*nPts + jj);
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(uint ii=0; ii < num - 1; ++ii){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

namespace netCDF
{
    // Functions for reading
    void Get_Dimension(int &fin, string const& variable, int &dimID, size_t &dimVal)
    {
        int retval;
        if ((retval = nc_inq_dimid(fin, variable.c_str(), &dimID)))
        {
            cout << "\n\tError: failed getting the ID of dimension \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed getting the ID of dimension \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        if ((retval = nc_inq_dimlen(fin, dimID, &dimVal)))
        {
            cout << "\n\tError: failed getting the length of dimension \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed getting the length of dimension \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }
    }

    double* Get_Real_Scalar(int &fin, string const& variable, size_t const &nPnts)
    {
        int retval, varID;
        if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
        {
            cout << "\n\tError: failed to get ID for variable \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed to get ID for variable \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        double *var = new double[nPnts];

        if ((retval = nc_get_var_double(fin, varID, &var[0])))
        {
            cout << "\n\tError: failed to get data for variable  \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed to get data for variable  \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        return var;
    }

    int* Get_Int_Scalar(int &fin, string const& variable, size_t const &nPnts)
    {
        int retval, varID;
        if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
        {
            cout << "\n\tError: failed to get ID for variable \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed to get ID for variable \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        int* var = new int[nPnts];

        if ((retval = nc_get_var_int(fin, varID, &var[0])))
        {
            cout << "\n\tError: failed to get data for variable  \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed to get data for variable  \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        return var;
    }

    /*To run on the mesh file*/
    vector<vector<size_t>> Get_Element(int &fin, string const& variable, size_t &nElem, size_t &nPnts)
    {
        #ifdef DEBUG
            dbout << "Reading Element: " << variable << endl;
        #endif
        int retval, varID;

        if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
        {
            cout << "\n\tError: failed to get ID for variable \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed to get ID for variable \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        // cout << nElem << "  " << nPoints << endl;
        #ifdef DEBUG
            dbout << "Allocating array of: " << nElem << " by " << nPnts << endl;
        #endif

        /*Allocate on the heap (can be big datasets)*/
        int *elemArray = new int[nElem * nPnts];

        /*Get the actual data from the file*/
        size_t startp[] = {0, 0};
        size_t countp[] = {nElem, nPnts};

        #ifdef DEBUG
            dbout << "Attempting to read NetCDF elements." << endl;
        #endif

        if ((retval = nc_get_vara_int(fin, varID, startp, countp, &elemArray[0])))
        {
            cout << "\n\tError: failed to get data for variable  \"" << variable << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: failed to get data for variable  \"" << variable << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }

        // cout << "Successfully read: " << variable << endl;
        // cout << "Number of elements: " << nElem << endl;

        #ifdef DEBUG
            dbout << "Successfully read elements" << endl;
        #endif

        /*Convert it to a vector to store*/
        size_t ii, jj;
        vector<vector<size_t>> elemVec(nElem, vector<size_t>(nPnts));
        for (ii = 0; ii < nElem; ++ii)
        {
            for (jj = 0; jj < nPnts; ++jj)
                elemVec[ii][jj] = static_cast<size_t>(elemArray[index(ii, jj, nPnts)]);
        }

        #ifdef DEBUG
            dbout << "Returning vector" << endl;
        #endif
        return elemVec;
    }
        
    /*To run on the mesh file*/
    vector<std::array<real,3>> Get_Coordinate_Vector(int &fin, size_t const &nPnts)
    {
        #ifdef DEBUG
            dbout << "Reading coordinates." << endl;
        #endif

        /*Get the coordinate data*/
        double *coordX = Get_Real_Scalar(fin, "points_xc", nPnts);
        double *coordY = Get_Real_Scalar(fin, "points_yc", nPnts);
        double *coordZ = Get_Real_Scalar(fin, "points_zc", nPnts);

        /*Convert it to a vector to store*/
        vector<std::array<real,3>> coordVec(nPnts);
        for (uint ii = 0; ii < nPnts; ++ii)
        {
            coordVec[ii] = std::array<real,3>{coordX[ii], coordY[ii], coordZ[ii]};
        }
        #ifdef DEBUG
            dbout << "Returning coordinates." << endl;
        #endif
        return coordVec;
    }

    // Functions for writing
    void Define_Dimension(int& fout, string const& dim_name, int const& dim_val, int* dim_id)
    {
        int retval;
        if ((retval = nc_def_dim(fout, dim_name.c_str(), dim_val, dim_id)))
        {
            cout << "\n\tError: Failed to define \"" << dim_name << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to define \"" << dim_name << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }
    }

    void Define_Variable(int& fout, string const& var_name, int const& var_type,
                        int const& ndims, int* dim_id, int* var_id)
    {
        int retval;
        if ((retval = nc_def_var(fout, var_name.c_str(), var_type, ndims, dim_id, var_id)))
        {
            cout << "\n\tError: Failed to define \"" << var_name << "\"" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to define \"" << var_name << "\"" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }
    }

    void Write_Variable_Array(int& fout, int const& var_id, size_t* const start,
                        size_t* const end, int* const var_data, string const& var_name)
    {
        int retval;
        if ((retval = nc_put_vara_int(fout, var_id, start, end, &var_data[0])))
        {
            cout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }	
    }

    void Write_Variable_Array(int& fout, int const& var_id, size_t* const start,
                        size_t* const end, float* const var_data, string const& var_name)
    {
        int retval;
        if ((retval = nc_put_vara_float(fout, var_id, start, end, &var_data[0])))
        {
            cout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }	
    }

    void Write_Variable_Array(int& fout, int const& var_id, size_t* const start,
                        size_t* const end, double* const var_data, string const& var_name)
    {
        int retval;
        if ((retval = nc_put_vara_double(fout, var_id, start, end, &var_data[0])))
        {
            cout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }	
    }

    void Write_Variable_Scalar(int& fout, int const& var_id,
                    int* const var_data, string const& var_name)
    {
        int retval;
        if ((retval = nc_put_var_int(fout, var_id, &var_data[0])))
        {
            cout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }	
    }

    void Write_Variable_Scalar(int& fout, int const& var_id,
                    float* const var_data, string const& var_name)
    {
        int retval;
        if ((retval = nc_put_var_float(fout, var_id, &var_data[0])))
        {
            cout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }	
    }

    void Write_Variable_Scalar(int& fout, int const& var_id, 
                    double* const var_data, string const& var_name)
    {
        int retval;
        if ((retval = nc_put_var_double(fout, var_id, &var_data[0])))
        {
            cout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #ifdef DEBUG
                dbout << "\n\tError: Failed to write \"" << var_name << "\" data" << endl;
            #endif
            ERR(retval);
            exit(-1);
        }	
    }
}

#endif