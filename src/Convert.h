#ifndef CONVERT_H
#define CONVERT_H

#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <netcdf>

#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"

#define NC_ERR 2
#define ERR(e)                             \
{                                          \
    printf("Error: %s\n", nc_strerror(e)); \
    exit(-1);                              \
}

using std::cout;
using std::endl;
using std::setw;
using std::string;
using std::vector;

/* Define data type. */
typedef double real;
typedef unsigned int uint;

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

void Get_Dimension(int &fin, string variable, int &dimID, size_t &dimVal)
{
    int retval;
    if ((retval = nc_inq_dimid(fin, variable.c_str(), &dimID)))
    {
        cout << "Error: failed getting the ID of dimension \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }

    if ((retval = nc_inq_dimlen(fin, dimID, &dimVal)))
    {
        cout << "Error: failed getting the length of dimension \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }
}

double* Get_Real_Scalar(int &fin, string variable, size_t const &nPnts)
{
    int retval, varID;
    if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
    {
        cout << "Error: failed to get ID for variable \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }

    double *var = new double[nPnts];

    if ((retval = nc_get_var_double(fin, varID, &var[0])))
    {
        cout << "Error: failed to get data for variable  \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }

    return var;
}

int* Get_Int_Scalar(int &fin, string variable, size_t const &nPnts)
{
    int retval, varID;
    if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
    {
        cout << "Error: failed to get ID for variable \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }

    int* var = new int[nPnts];

    if ((retval = nc_get_var_int(fin, varID, &var[0])))
    {
        cout << "Error: failed to get data for variable  \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }

    return var;
}

/*To run on the mesh file*/
vector<vector<size_t>> Get_Element(int &fin, string variable, size_t &nElem, size_t &nPnts)
{
    #ifdef DEBUG
        dbout << "Reading Element: " << variable << endl;
    #endif
    int retval, varID;

    if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
    {
        cout << "Error: failed to get ID for variable \"" << variable << "\"" << endl;
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
        cout << "Error: failed to get data for variable  \"" << variable << "\"" << endl;
        ERR(retval);
        exit(-1);
    }

    cout << "Successfully read: " << variable << endl;
    cout << "Number of elements: " << nElem << endl;

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
vector<Eigen::Matrix<real,3,1>> Get_Coordinate_Vector(int &fin, size_t const &nPnts)
{
    #ifdef DEBUG
        dbout << "Reading coordinates." << endl;
    #endif

    /*Get the coordinate data*/
    double *coordX = Get_Real_Scalar(fin, "points_xc", nPnts);
    double *coordY = Get_Real_Scalar(fin, "points_yc", nPnts);
    double *coordZ = Get_Real_Scalar(fin, "points_zc", nPnts);

    /*Convert it to a vector to store*/
    vector<Eigen::Matrix<real,3,1>> coordVec(nPnts);
    for (uint ii = 0; ii < nPnts; ++ii)
    {
        coordVec[ii] = Eigen::Matrix<real,3,1>(coordX[ii], coordY[ii], coordZ[ii]);
    }
    #ifdef DEBUG
        dbout << "Returning coordinates." << endl;
    #endif
    return coordVec;
}

#endif