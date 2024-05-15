/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Geometry.h"
#include "IOFunctions.h"
#include "Var.h"
#include <netcdf.h>
#include <set>
#define NC_ERR 2
#define ERR(e)                                                                                          \
    {                                                                                                   \
        printf("Error: %s\n", nc_strerror(e));                                                          \
        exit(-1);                                                                                       \
    }

void Average_Point_to_Cell(
    vector<StateVecD> const& pData, vector<StateVecD>& cData, vector<vector<size_t>> const& cFaces,
    vector<vector<size_t>> const& faces
)
{
    cData = vector<StateVecD>(cFaces.size());
#pragma omp parallel for default(shared)
    for (uint ii = 0; ii < cFaces.size(); ++ii)
    {
        std::set<size_t> elem;
        for (auto const& faceID : cFaces[ii])
        {
            vector<size_t> const face = faces[faceID];
            for (auto const& vert : face)
            { /*use set to ensure vertex doesn't exist in the elems vector yet.*/
                elem.insert(vert);
            }
        }

        StateVecD sum_ = StateVecD::Zero();
        for (auto const& jj : elem)
        {
            sum_ += pData[jj];
        }
        cData[ii] = sum_ / elem.size();
    }
}

void Average_Point_to_Cell(
    vector<real> const& pData, vector<real>& cData, vector<vector<size_t>> const& cFaces,
    vector<vector<size_t>> const& faces
)
{
    cData = vector<real>(cFaces.size(), 0.0);
#pragma omp parallel for default(shared)
    for (uint ii = 0; ii < cFaces.size(); ++ii)
    {
        std::set<size_t> elem;
        for (auto const& faceID : cFaces[ii])
        {
            vector<size_t> const face = faces[faceID];
            for (auto const& vert : face)
            { /*use set to ensure vertex doesn't exist in the elems vector yet.*/
                elem.insert(vert);
            }
        }

        real sum_ = 0.0;
        for (auto const& jj : elem)
        {
            sum_ += pData[jj];
        }
        cData[ii] = sum_ / elem.size();
    }
}

real KahanSum(vector<real> const& values, std::set<size_t> const& index)
{
    real sum = 0.0; // Prepare the accumulator.
    real c = 0.0;   // A running compensation for lost low-order bits.

    for (size_t const& ii : index)
    {
        real y = values[ii] - c; // c is zero the first time around.
        real t = sum + y;        // Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y; // (t - sum) cancels the high-order part of y; subtracting y recovers negative
                           // (low part of y)
        sum = t;           // Algebraically, c should always be zero. Beware overly-aggressive optimizing
                 // compilers! Next time around, the lost low part will be added to y in a fresh attempt.
    }
    return sum;
}

StateVecD KahanSum(vector<StateVecD> const& values, std::set<size_t> const& index)
{
    StateVecD sumVec;
    for (size_t dim = 0; dim < SIMDIM; dim++)
    {
        real sum = 0.0; // Prepare the accumulator.
        real c = 0.0;   // A running compensation for lost low-order bits.

        for (size_t const& ii : index)
        {
            real y = values[ii][dim] - c; // c is zero the first time around.
            real t = sum + y;  // Alas, sum is big, y small, so low-order digits of y are lost.
            c = (t - sum) - y; // (t - sum) cancels the high-order part of y; subtracting y recovers
                               // negative (low part of y)
            sum = t; // Algebraically, c should always be zero. Beware overly-aggressive optimizing
                     // compilers! Next time around, the lost low part will be added to y in a fresh
                     // attempt.
        }
        sumVec[dim] = sum;
    }
    return sumVec;
}

void Average_Point_Data_to_Cell(
    SIM const& svar, vector<StateVecD> const& vel, vector<real> const& press, vector<real> const& dens,
    MESH& cells
)
{
#pragma omp parallel for default(shared)
    for (uint ii = 0; ii < cells.cFaces.size(); ++ii)
    {
        std::set<size_t> elem;
        for (auto const& faceID : cells.cFaces[ii])
        {
            vector<size_t> const face = cells.faces[faceID];
            for (auto const& vert : face)
            { /*use set to ensure vertex doesn't exist in the elems vector yet.*/
                elem.insert(vert);
            }
        }

        StateVecD cCentre_ = StateVecD::Zero();
        StateVecD cVel_ = StateVecD::Zero();
        real cPress_ = 0.0;
        real cRho_ = 0.0;

        real const nVerts = real(elem.size());

        // for (auto const& jj : elem)
        // {
        // 	cCentre_ += cells.verts[jj];
        // 	cVel_ += vel[jj];
        // 	cPress_ += press[jj];
        // 	cRho_ += dens[jj];
        // }

        cCentre_ = KahanSum(cells.verts, elem);
        cVel_ = KahanSum(vel, elem);
        cPress_ = KahanSum(press, elem);
        cRho_ = KahanSum(dens, elem);

        cells.cCentre[ii] = cCentre_ / nVerts;
        cells.cVel[ii] = cVel_ / nVerts;
        cells.cP[ii] = cPress_ / nVerts;
        cells.cRho[ii] = cRho_ / nVerts;
    }
}

/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
namespace TAU
{
    void Get_Dim_ID(int& meshID, string const& dim, int& dimID)
    {
        int retval = 0;
        int dimID_;
        if ((retval = nc_inq_dimid(meshID, dim.c_str(), &dimID_)))
        {
            printf("Failed to get the dimension ID of \"%s\"\n", dim.c_str());
            exit(-1);
        }
        dimID = dimID_;
    }

    void Get_Dim_Length(int& meshID, int const& dimID, size_t& dimLen)
    {
        int retval = 0;
        size_t dimLen_;
        if ((retval = nc_inq_dimlen(meshID, dimID, &dimLen_)))
        {
            printf("Failed to get the dimension length of ID %d\n", dimID);
            exit(-1);
        }
        dimLen = dimLen_;
    }

    void Get_Var_ID(int& meshID, string const& var, int& varID)
    {
        int retval = 0;
        int varID_;
        if ((retval = nc_inq_varid(meshID, var.c_str(), &varID_)))
        {
            printf("Failed to get the variable ID of \"%s\"\n", var.c_str());
            exit(-1);
        }
        varID = varID_;
    }

    void Find_Angle_Alpha(SIM& svar)
    {
#ifdef DEBUG
        fprintf(dbout, "Entering Find_Bmap_Markers...\n");
#endif
        string const& bmapIn = svar.tau_bmap;
        std::ifstream fin(bmapIn, std::ios::in);

        if (!fin.is_open())
        {
            cout << "Couldn't open the boundary map file." << endl;
            cout << "Attempted path: " << bmapIn << endl;
            exit(-1);
        }

        string line;
        while (getline(fin, line))
        {
            Get_Number(line, "Angle alpha (degree)", svar.angle_alpha);
        }

        svar.angle_alpha *= M_PI / 180.0;

        StateVecD oldgrav = svar.grav;
#if SIMDIM == 3
        svar.grav[2] = oldgrav[2] * cos(svar.angle_alpha);
#else
        svar.grav[1] = oldgrav[1] * cos(svar.angle_alpha);
#endif
        svar.grav[0] = -oldgrav[1] * sin(svar.angle_alpha);

#ifdef DEBUG
        fprintf(dbout, "Exiting Find_Angle_Alpha...\n");
#endif
    }

    void Read_BMAP(SIM& svar)
    {
#ifdef DEBUG
        fprintf(dbout, "Entering Find_Bmap_Markers...\n");
#endif
        string const& bmapIn = svar.tau_bmap;
        std::ifstream fin(bmapIn, std::ios::in);

        if (!fin.is_open())
        {
            cout << "Couldn't open the boundary map file." << endl;
            cout << "Attempted path: " << bmapIn << endl;
            exit(-1);
        }

        /* Find out how many blocks, i.e. how many boundaries to expect. */
        uint nBlocks = 0;
        string line;
        while (getline(fin, line))
        {
            if (line.find("block end") != string::npos)
            { /*We're on a new block*/
                nBlocks++;
            }
        }
        fin.clear();
        fin.seekg(0);

        uint blockno = 0;

        svar.tau_markers = vector<int>(nBlocks, 0);
        svar.tau_bnames = vector<string>(nBlocks);
        svar.tau_bwrite = vector<int>(nBlocks, 0);

        vector<string> types(nBlocks);

        while (getline(fin, line))
        {
            /* Remove whitespace */
            line = ltrim(line);
            /* Check if line is commented out */
            if (line[0] == '#')
                continue;

            Get_Number(line, "Angle alpha (degree)", svar.angle_alpha);
            Get_Number(line, "tau_markers", svar.tau_markers[blockno]);
            Get_String(line, "Name", svar.tau_bnames[blockno]);
            Get_String(line, "Type", types[blockno]);
            Get_Number(line, "Write surface data (0/1)", svar.tau_bwrite[blockno]);

            if (line.find("block end") != string::npos)
            { /*We're on a new block*/
                blockno++;
            }
        }

        fin.close();

        /* Check if the boundary has a name. If not, then give it the type name*/
        for (size_t ii = 0; ii < nBlocks; ++ii)
        {
            if (svar.tau_bnames[ii].empty())
            {
                svar.tau_bnames[ii] = types[ii];
            }
        }

        svar.angle_alpha *= M_PI / 180.0;

        StateVecD oldgrav = svar.grav;
#if SIMDIM == 3
        svar.grav[2] = oldgrav[2] * cos(svar.angle_alpha);
#else
        svar.grav[1] = oldgrav[1] * cos(svar.angle_alpha);
#endif
        svar.grav[0] = -oldgrav[1] * sin(svar.angle_alpha);

#ifdef DEBUG
        fprintf(dbout, "Exiting Find_Bmap_Markers...\n");
#endif
    }

    /*To run on the solution file*/
    vector<real> Get_Scalar_Property_real(int& fin, string const& variable, size_t const& nPts)
    {
#ifdef DEBUG
        fprintf(dbout, "Reading variable: %s\n", variable.c_str());
#endif
        int retval;
        int varID;

        Get_Var_ID(fin, variable.c_str(), varID);

#ifdef DEBUG
        fprintf(dbout, "Allocating array of: %zu", nPts);
#endif

        double* array = new double[nPts];

#ifdef DEBUG
        fprintf(dbout, "Attempting to read NetCDF variable: %s\n", variable.c_str());
#endif

        if ((retval = nc_get_var_double(fin, varID, &array[0])))
        {
            printf("Failed to get variable data for \"%s\"\n", variable.c_str());
            exit(-1);
        }

        /*Convert it to a vector to store*/
        vector<double> propVec;
        propVec.insert(propVec.end(), &array[0], &array[nPts]);
        vector<real> var = propVec;

#ifdef DEBUG
        fprintf(dbout, "Returning vector\n");
#endif
        return var;
    }

    vector<int> Get_Scalar_Property_int(int& fin, string const& variable, int const& nPts)
    {
#ifdef DEBUG
        fprintf(dbout, "Reading variable: %s\n", variable.c_str());
#endif
        int retval;
        int varID;

        Get_Var_ID(fin, variable.c_str(), varID);

#ifdef DEBUG
        fprintf(dbout, "Allocating array of: %d\n", nPts);
#endif

        int* array = new int[nPts];

#ifdef DEBUG
        fprintf(dbout, "Attempting to read NetCDF variable: %s\n", variable.c_str());
#endif

        if ((retval = nc_get_var_int(fin, varID, &array[0])))
        {
            printf("Failed to get the variable data of \"%s\"\n", variable.c_str());
            exit(-1);
        }

        /*Convert it to a vector to store*/
        vector<int> propVec;
        propVec.insert(propVec.end(), &array[0], &array[nPts]);

#ifdef DEBUG
        fprintf(dbout, "Returning vector\n");
#endif
        return propVec;
    }

    /*To run on the mesh file*/
    void Get_Element(
        int& fin, string const& variable, size_t const& nElem, size_t const& nPpEd,
        vector<vector<size_t>>& vec
    )
    {
#ifdef DEBUG
        fprintf(dbout, "Reading variable: %s\n", variable.c_str());
#endif
        int retval;
        int varID;

        if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
        {
            printf("Failed to get the variable ID of \"%s\"\n", variable.c_str());
            exit(-1);
        }

        // cout << nElem << "  " << nPoints << endl;
#ifdef DEBUG
        fprintf(dbout, "Allocating array of: %zu by %zu\n", nElem, nPpEd);
#endif

        vector<vector<size_t>> elemVec(nElem, vector<size_t>(nPpEd));

        try
        {
            /*Allocate on the heap (can be big datasets)*/
            int* elemArray = new int[nElem * nPpEd];

            size_t start[] = {0, 0};
            size_t end[] = {nElem, nPpEd};

            /*Get the actual data from the file*/
            if ((retval = nc_get_vara_int(fin, varID, start, end, &elemArray[0])))
            {
                printf("Failed to get the variable data of \"%s\"\n", variable.c_str());
                exit(-1);
            }

#ifdef DEBUG
            fprintf(dbout, "Attempting to read NetCDF elements.\n");
#endif

            cout << "Successfully read \"" << variable << "\"" << endl;
            cout << "Number of cells: " << nElem << endl;

#ifdef DEBUG
            fprintf(dbout, "Successfully read variable: %s\n", variable.c_str());
#endif

            /*Convert it to a vector to store*/
            for (size_t ii = 0; ii < nElem; ++ii)
                for (size_t jj = 0; jj < nPpEd; ++jj)
                    elemVec[ii][jj] = static_cast<size_t>(elemArray[index(jj, ii, nPpEd)]);
        }
        catch (std::bad_alloc& ba)
        {

            std::cerr << "Bad alloc caught. Failed to allocate \"" << variable << "\"" << endl;
            exit(-1);
        }

#ifdef DEBUG
        fprintf(dbout, "Returning vector\n");
#endif
        vec.insert(vec.end(), elemVec.begin(), elemVec.end());

        // return elemVec;
    }

    /*To run on the mesh file*/
    uint Get_Coordinates(int& fin, size_t const& nPnts, vector<StateVecD>& coordVec)
    {
#ifdef DEBUG
        fprintf(dbout, "Reading coordinates.\n");
#endif

        int retval;
        int xcID, ycID, zcID;

#if SIMDIM == 3

        if ((retval = nc_inq_varid(fin, "points_xc", &xcID)))
        {
            cout << "Failed to get the variable ID of: " << "points_xc" << endl;
            cout << "Stopping" << endl;
            exit(-1);
        }

        if ((retval = nc_inq_varid(fin, "points_yc", &ycID)))
        {
            cout << "Failed to get the variable ID of: " << "points_yc" << endl;
            cout << "Stopping" << endl;
            exit(-1);
        }

        if ((retval = nc_inq_varid(fin, "points_zc", &zcID)))
        {
            cout << "Failed to get the variable ID of: " << "points_zc" << endl;
            cout << "Stopping" << endl;
            exit(-1);
        }

#else
        uint ignore = 0;
        uint nFail = 0;
        if ((retval = nc_inq_varid(fin, "points_xc", &xcID)))
        {
            cout << "Failed to get the variable ID of: " << "points_xc" << endl;
            cout << "Assuming this is the ignored dimension" << endl;
            ignore = 1;
            nFail++;
        }

        if ((retval = nc_inq_varid(fin, "points_yc", &ycID)))
        {
            cout << "Failed to get the variable ID of: " << "points_yc" << endl;
            cout << "Assuming this is the ignored dimension" << endl;
            ignore = 2;
            nFail++;
        }

        if ((retval = nc_inq_varid(fin, "points_zc", &zcID)))
        {
            cout << "Failed to get the variable ID of: " << "points_zc" << endl;
            cout << "Assuming this is the ignored dimension" << endl;
            ignore = 3;
            nFail++;
        }

        if (nFail > 1)
        {
            cout << "More than one dimension was not aquired, meaning something went wrong.\n\tStopping"
                 << endl;
            exit(-1);
        }

#ifdef DEBUG
        fprintf(dbout, "Ignored dimension: %u\n", ignore);
#endif
        cout << "Ignored dimension: " << ignore << endl;
#endif

#ifdef DEBUG
        fprintf(dbout, "Number of points: %zu\n", nPnts);
#endif
        /*Allocate on the heap (can be big datasets)*/
        real* coordX = new real[nPnts];
        real* coordY = new real[nPnts];

        coordVec = vector<StateVecD>(nPnts);
#if SIMDIM == 2
        /*Get the actual data from the file*/
        if (ignore == 1)
        {
            if ((retval = nc_get_var_double(fin, ycID, &coordX[0])))
            {
                printf("Failed to get the variable data of y-coordinate\n");
                exit(-1);
            }
            if ((retval = nc_get_var_double(fin, zcID, &coordY[0])))
            {
                printf("Failed to get the variable data of z-coordinate\n");
                exit(-1);
            }
        }
        else if (ignore == 2)
        {
            if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
            {
                printf("Failed to get the variable data of x-coordinate\n");
                exit(-1);
            }
            if ((retval = nc_get_var_double(fin, zcID, &coordY[0])))
            {
                printf("Failed to get the variable data of z-coordinate\n");
                exit(-1);
            }
        }
        else if (ignore == 3)
        {
            if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
            {
                printf("Failed to get the variable data of x-coordinate\n");
                exit(-1);
            }
            if ((retval = nc_get_var_double(fin, ycID, &coordY[0])))
            {
                printf("Failed to get the variable data of y-coordinate\n");
                exit(-1);
            }
        }
        else
        {
            cout << "The ignored dimension was not found." << endl;
            exit(-1);
        }

        for (uint ii = 0; ii < nPnts; ++ii)
        { /*Convert it to a vector to store*/
            coordVec[ii] = StateVecD(coordX[ii], coordY[ii]);
        }

#else
        real* coordZ = new real[nPnts];
        if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
        {
            printf("Failed to get the variable data of x-coordinate\n");
            exit(-1);
        }
        if ((retval = nc_get_var_double(fin, ycID, &coordY[0])))
        {
            printf("Failed to get the variable data of y-coordinate\n");
            exit(-1);
        }
        if ((retval = nc_get_var_double(fin, zcID, &coordZ[0])))
        {
            printf("Failed to get the variable data of z-coordinate\n");
            exit(-1);
        }

#pragma omp parallel for
        for (uint ii = 0; ii < nPnts; ++ii)
        { /*Convert it to a vector to store*/
            coordVec[ii] = StateVecD(coordX[ii], coordY[ii], coordZ[ii]);
        }
#endif

#ifdef DEBUG
        fprintf(dbout, "Returning coordinates.\n");
#endif

#if SIMDIM == 3
        return 0;
#else
        return ignore;
#endif
    }

    void Get_Cell_Faces(
        vector<vector<uint>> const& cell, vector<vector<uint>> const& facenum,
        std::vector<std::vector<std::vector<uint>>>& cFaces
    )
    {
        for (uint ii = 0; ii < cell.size(); ++ii)
        {
            uint jj = 0;
            cFaces.emplace_back();
            for (auto faces : facenum)
            {
                cFaces[ii].emplace_back();
                for (auto vert : faces)
                {
                    cFaces[ii][jj].emplace_back(cell[ii][vert]);
                }
                ++jj;
            }
        }
    }

    /*****************************************************************************/
    /***************** READING NETCDF SOLUTION DATA FUNCTIONS ********************/
    /*****************************************************************************/

    void Read_SOLUTION(
        SIM const& svar, FLUID const& fvar, AERO const& avar, uint const ignored, MESH& cells,
        vector<uint> const& usedVerts
    )
    {
        cout << "Reading solution file..." << endl;

#ifdef DEBUG
        fprintf(dbout, "Opening solultion file.\n");
#endif
        string const& solIn = svar.tau_sol;
        int retval;
        int solID;

        if ((retval = nc_open(solIn.c_str(), NC_NOWRITE, &solID)))
        {
            cout << "A netCDF error occured whilst trying to open the solution file:" << endl;
            cout << "\t" << solIn << endl << endl;
            ERR(retval);
            exit(-1);
        }

        cout << "Solution file opened, reading contents..." << endl;

        int ptDimID;
        size_t solPts;

        Get_Dim_ID(solID, "no_of_points", ptDimID);
        Get_Dim_Length(solID, ptDimID, solPts);

#ifdef DEBUG
        fprintf(dbout, "Solution points: %zu\n", solPts);
#endif

#if SIMDIM == 3
        if (solPts != cells.nPnts)
        {
            cout << "Solution file does not have the same number of vertices as the mesh." << endl;
            cout << "Please check again." << endl;
            exit(-1);
        }
#else
        if (solPts / 2 != usedVerts.size())
        {
            cout << "Solution file size does not match size of mesh. Please check inputs." << endl;
        }
#endif

        vector<real> realDens = Get_Scalar_Property_real(solID, "density", solPts);

        /*Get the velocities*/
        vector<real> xvel, yvel, zvel;
#if SIMDIM == 3
        xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
        yvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
        zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);
#else
        if (ignored == 1)
        {
            xvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
            zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);
        }
        else if (ignored == 2)
        {
            xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
            zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);
        }
        else if (ignored == 3)
        {
            xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
            zvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
        }
#endif

        vector<StateVecD> vel = vector<StateVecD>(solPts);
        /*Test for size*/
        if (xvel.size() == solPts)
        { /*Turn the arrays into a state vector*/
#pragma omp parallel for default(shared)
            for (uint ii = 0; ii < solPts; ++ii)
            {
#if SIMDIM == 3
                vel[ii] = StateVecD(xvel[ii], yvel[ii], zvel[ii]);
#else
                vel[ii] = StateVecD(xvel[ii], zvel[ii]);
#endif
            }

            if (usedVerts.size() != 0)
            { /*Get the used vertices*/
                vector<StateVecD> newVel;
                for (auto index : usedVerts)
                {
                    newVel.emplace_back(vel[index]);
                }
                vel = newVel;
            }
            else if (cells.verts.size() == solPts)
            { // Do nothing
              // vel = vel;
            }
            else
            {
                cout << "Don't know which values to use!" << endl;
                exit(-1);
            }
        }
        else
        {
            cout << "velocities do not have the same number of vertices as the mesh." << endl;
            cout << xvel.size() << "  " << solPts << endl;
            cout << "Please check again." << endl;
            exit(-1);
        }

        /*Get other scalar data that may or may not exist*/
        vector<real> press = Get_Scalar_Property_real(solID, "pressure", solPts);
        vector<real> dens = Get_Scalar_Property_real(solID, "density", solPts);

        if ((retval = nc_close(solID)))
        {
            cout << "A netCDF error occured whilst trying to close the solution file:" << endl;
            cout << "\t" << solIn << endl << endl;
            ERR(retval);
            exit(-1);
        }
        // vector<real> dens(solPts);

        // #pragma omp parallel for
        // for (uint ii = 0; ii < press.size(); ++ii)
        // {
        // 	press[ii] -= avar.pRef; /* Want value to be gauge pressure */
        // 	// dens[ii] = fvar.rho0 * pow((press[ii] / fvar.B + 1), 1 / fvar.gam);
        // }

        if (press.size() == 0)
        {
            cout << "The solution data has not been correctly ingested. Please check the solution file."
                 << endl;
        }

        /*Update the point arrays with the actual stuff*/
        if (usedVerts.size() != 0)
        { /*Get the used vertices*/
            vector<real> newP;
            vector<real> newR;
            newP.reserve(usedVerts.size());
            newR.reserve(usedVerts.size());

            for (auto index : usedVerts)
            {
                newP.emplace_back(press[index]);
                newR.emplace_back(dens[index]);
            }
            press = newP;
            dens = newR;
        }
        else if (cells.verts.size() == solPts)
        { // Do nothing
          //  cells.pointP = press;
          //  cells.pointRho = dens;
        }
        else
        {
            cout << "Don't know which values to use!" << endl;
            exit(-1);
        }

        Average_Point_Data_to_Cell(svar, vel, press, dens, cells);
    }

    /*****************************************************************************/
    /*************** READING NETCDF EDGE BASED DATA FUNCTIONS ********************/
    /*****************************************************************************/

#if SIMDIM == 2

    void
    Place_Edges(int& fin, size_t const& nElem, size_t const& nPnts, size_t const& nEdge, MESH& cells)
    {
#ifdef DEBUG
        fprintf(dbout, "Reading element left/right and placing faces\n");
#endif

        vector<int> left = Get_Scalar_Property_int(fin, "left_element_of_edges", nEdge);
        vector<int> right = Get_Scalar_Property_int(fin, "right_element_of_edges", nEdge);

        // vector<std::pair<int, int>> leftright(nEdge);

        // #pragma omp parallel default(shared) //shared(leftright)
        // {
        // 	/*Create local of */

        // 	#pragma omp for schedule(static) nowait
        // 	for (size_t ii = 0; ii < nEdge; ++ii)
        // 	{
        // 		int lindex = left[ii];
        // 		int rindex = right[ii];
        // 		leftright[ii] = std::pair<int, int>(lindex, rindex);
        // 		#pragma omp critical
        // 		{
        // 			cells.cFaces[lindex].emplace_back(ii);
        // 			if (rindex >= 0)
        // 				cells.cFaces[rindex].emplace_back(ii);
        // 		}
        // 	}
        // }

        vector<int> tau_markers =
            Get_Scalar_Property_int(fin, "boundarymarker_of_surfaces", cells.nSurf);

        vector<std::pair<int, int>> leftright(nEdge);
        cells.smarkers = vector<std::pair<size_t, int>>(cells.nSurf);

        real& maxedge = cells.maxlength;
        real& minedge = cells.minlength;
        maxedge = 0.0;

        vector<size_t> surf_IDs;
#pragma omp parallel default(shared) /* shared(leftright) */ reduction(max : maxedge)                   \
    reduction(min : minedge)
        {
            vector<size_t> local;
#pragma omp for schedule(static) nowait
            for (size_t ii = 0; ii < nEdge; ++ii)
            {
                int lindex = left[ii];
                int rindex = right[ii];
                leftright[ii] = std::pair<int, int>(lindex, rindex);
#pragma omp critical
                {
                    cells.cFaces[lindex].emplace_back(ii);
                    if (rindex >= 0)
                        cells.cFaces[rindex].emplace_back(ii);
                    else
                    {
                        local.emplace_back(ii);
                    }
                }

                /* Find the longest edge */
                real max_e = 0.0;

                vector<size_t> const& face = cells.faces[ii];
                max_e = (cells.verts[face[0]] - cells.verts[face[1]]).norm();

                maxedge = std::max(maxedge, max_e); // Provides a sanity reference for particle tracking
                minedge = std::min(minedge, max_e);
            }

#pragma omp for schedule(static) ordered
            for (int ii = 0; ii < omp_get_num_threads(); ++ii)
            {
#pragma omp ordered
                {
                    surf_IDs.insert(surf_IDs.end(), local.begin(), local.end());
                }
            }

#pragma omp single
            {
                if (surf_IDs.size() != cells.nSurf)
                {
                    cout << "Mismatch of number of surface faces identified, and the number given"
                         << endl;
                    cout << "Identified: " << surf_IDs.size() << "  Given: " << cells.nSurf << endl;
                    exit(-1);
                }
            }

#pragma omp for schedule(static) nowait
            for (size_t ii = 0; ii < cells.nSurf; ++ii)
            {
                cells.smarkers[ii].first = surf_IDs[ii];
                cells.smarkers[ii].second = tau_markers[ii];
            }
        } /* End parallel */

        maxedge *= 4.0;

        cells.leftright = leftright;

#ifdef DEBUG
        fprintf(dbout, "End of placing edges in elements.\n\n");
#endif

        /*Now go through the faces and see which vertices are unique, to get element data*/
        // cout << "Building elements..." << endl;
        // #pragma omp parallel default(shared)
        // {
        // 	// vector<vector<uint>> local;
        // 	#pragma omp for schedule(static) nowait
        // 	for (uint ii = 0; ii < cells.cFaces.size(); ++ii)
        // 	{
        // 		for (auto const &index : cells.cFaces[ii])
        // 		{
        // 			vector<size_t> const face = cells.faces[index];
        // 			for (auto const &vert : face)
        // 			{
        // 				if (std::find(cells.elems[ii].begin(), cells.elems[ii].end(),
        // vert) == cells.elems[ii].end()) 				{ /*Vertex doesn't exist in the
        // elems vector yet.*/ 					#pragma omp critical
        // 					{
        // 						cells.elems[ii].emplace_back(vert);
        // 					}
        // 				}
        // 			}
        // 		}
        // 	}

        // 	/*Find cell centres*/
        // 	#pragma omp single
        // 	{
        // 		cout << "Finding cell centres..." << endl;
        // 	}

        // Average_Point_to_Cell(cells.verts, cells.cCentre, cells.elems);

        // 	/*Find cell centres*/
        // 	#pragma omp single
        // 	{
        // 		cout << "Finding cell volumes..." << endl;
        // 	}

        // 	// Find cell volumes
        // 	#pragma omp for schedule(static) nowait
        // 	for (size_t ii = 0; ii < cells.elems.size(); ++ii)
        // 	{
        // 		cells.cVol[ii] = Cell_Volume(cells.verts, cells.faces, cells.elems[ii],
        // cells.cFaces[ii], cells.cCentre[ii]);
        // 	}
        // }
        // cout << "Elements built." << endl;

#ifdef DEBUG
        fprintf(dbout, "All elements defined.\n");
#endif
    }

    void
    Read_tau_mesh_EDGE(SIM& svar, MESH& cells, FLUID const& fvar, AERO const& avar, vector<uint>& uVerts)
    {
        string meshIn = svar.tau_mesh;
        string solIn = svar.tau_sol;

#ifdef DEBUG
        fprintf(dbout, "Attempting read of NetCDF file.\n");
        fprintf(dbout, "Mesh file: %s\n", meshIn.c_str());
        fprintf(dbout, "Solution file: %s\n", solIn.c_str());
#endif

        /*Read the mesh data*/
        int retval;
        int meshID;

        if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
        {
            cout << "A netCDF error occured whilst trying to open the mesh file:" << endl;
            cout << "\t" << meshIn << endl << endl;
            ERR(retval);
            exit(-1);
        }

        cout << "Mesh file open. Reading edge data..." << endl;

        int ptDimID, elemDimID, edgeDimID, nPpEDimID, surfDimID;
        size_t nPnts, nElem, nEdge, nPpEd, nSurf;

        // Retrieve how many elements there are.
        Get_Dim_ID(meshID, "no_of_elements", elemDimID);
        Get_Dim_Length(meshID, elemDimID, nElem);

        // Retrieve edge number
        Get_Dim_ID(meshID, "no_of_edges", edgeDimID);
        Get_Dim_Length(meshID, edgeDimID, nEdge);

        // Retrieve points per edge
        Get_Dim_ID(meshID, "points_per_edge", nPpEDimID);
        Get_Dim_Length(meshID, nPpEDimID, nPpEd);

        Get_Dim_ID(meshID, "no_of_surfaceelements", surfDimID);
        Get_Dim_Length(meshID, surfDimID, nSurf);

        Get_Dim_ID(meshID, "no_of_points", ptDimID);
        Get_Dim_Length(meshID, ptDimID, nPnts);

#ifdef DEBUG
        fprintf(dbout, "nElem: %zu nPnts: %zu nEdge: %zu nSurf: %zu\n", nElem, nPnts, nEdge, nSurf);
#endif

        cout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << " nSurf: " << nSurf
             << endl;

        cells.alloc(nPnts, nElem, nEdge, nSurf);

        // cells.nPnts = nPnts;
        // cells.nElem = nElem;
        // cells.nFace = nEdge;
        // cells.nSurf = nSurf;

        // cells.elems = vector<vector<size_t>>(nElem);
        // cells.cFaces = vector<vector<size_t>>(nElem);
        // cells.verts = vector<StateVecD>(nPnts);
        // cells.cVol = vector<real>(nElem);
        // cells.leftright = vector<std::pair<int,int>>(nFace);

        /*Get the faces of the mesh*/

        Get_Element(meshID, "points_of_element_edges", nEdge, nPpEd, cells.faces);

        /*Get the vertices that are in use to take the values from the solution file*/
        vector<int> usedVerts = Get_Scalar_Property_int(meshID, "vertices_in_use", nPnts);

#ifdef DEBUG
        fprintf(dbout, "Successfully read vertices_in_use data.\n");
#endif

        uVerts.resize(nPnts);
        for (size_t ii = 0; ii < nPnts; ++ii)
            uVerts[ii] = static_cast<uint>(usedVerts[ii]);

        /*Get the coordinates of the mesh*/
        (void)Get_Coordinates(meshID, nPnts, cells.verts);
        if (cells.verts.size() != nPnts)
        {
            cout << "Some data has been missed.\nPlease check how many points." << endl;
        }

        /*Adjust the scale*/
        cells.scale = svar.scale;
        for (auto& vert : cells.verts)
        {
            vert *= cells.scale;
        }
        // svar.Jet/=cells.scale;

        /*Get face left and right, and put the faces in the elements*/
        Place_Edges(meshID, nElem, nPnts, nEdge, cells);

#ifdef DEBUG
        fprintf(dbout, "End of interaction with mesh file and ingested data.\n\n");

#endif
    }

#endif

/*****************************************************************************/
/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
#if SIMDIM == 3
    void Place_Faces(int& fin, size_t const& nFace, SIM const& svar, MESH& cells)
    {
#ifdef DEBUG
        fprintf(dbout, "Reading element left/right and placing faces\n");
#endif

        vector<int> left = Get_Scalar_Property_int(fin, "left_element_of_faces", nFace);
        vector<int> right = Get_Scalar_Property_int(fin, "right_element_of_faces", nFace);
        vector<int> tau_markers =
            Get_Scalar_Property_int(fin, "boundarymarker_of_surfaces", cells.nSurf);

        vector<std::pair<int, int>> leftright(nFace);
        cells.smarkers = vector<std::pair<size_t, int>>(cells.nSurf);

        real& maxedge = cells.maxlength;
        real& minedge = cells.minlength;
        vector<size_t> surf_IDs;
#pragma omp parallel shared(leftright) reduction(max : maxedge) reduction(min : minedge)
        {
            vector<size_t> local;
#pragma omp for schedule(static) nowait
            for (size_t ii = 0; ii < nFace; ++ii)
            {
                int lindex = left[ii];
                int rindex = right[ii];
                leftright[ii] = std::pair<int, int>(lindex, rindex);
#pragma omp critical
                {
                    cells.cFaces[lindex].emplace_back(ii);
                    if (rindex >= 0)
                        cells.cFaces[rindex].emplace_back(ii);
                    else
                        local.emplace_back(ii);
                }

                /* Find the longest edge */
                real max_e = 0.0;
                real min_e = 0.0;
                if (cells.faces[ii].size() == 3)
                {
                    vector<size_t> const& face = cells.faces[ii];
                    max_e = std::max(
                        (cells.verts[face[0]] - cells.verts[face[1]]).norm(),
                        std::max(
                            (cells.verts[face[0]] - cells.verts[face[2]]).norm(),
                            (cells.verts[face[1]] - cells.verts[face[2]]).norm()
                        )
                    );

                    min_e = std::min(
                        (cells.verts[face[0]] - cells.verts[face[1]]).norm(),
                        std::min(
                            (cells.verts[face[0]] - cells.verts[face[2]]).norm(),
                            (cells.verts[face[1]] - cells.verts[face[2]]).norm()
                        )
                    );
                }
                else
                {
                    /* Assume the diagonals are the longest edges*/
                    /* Untrue for a trapezium, but not needed to be perfectly accurate */
                    vector<size_t> const& face = cells.faces[ii];
                    max_e = std::max(
                        (cells.verts[face[0]] - cells.verts[face[2]]).norm(),
                        (cells.verts[face[1]] - cells.verts[face[3]]).norm()
                    );

                    min_e = std::min(
                        (cells.verts[face[0]] - cells.verts[face[1]]).norm(),
                        std::min(
                            (cells.verts[face[1]] - cells.verts[face[2]]).norm(),
                            std::min(
                                (cells.verts[face[2]] - cells.verts[face[3]]).norm(),
                                (cells.verts[face[3]] - cells.verts[face[0]]).norm()
                            )
                        )
                    );
                }

                maxedge = std::max(maxedge, max_e);
                minedge = std::min(minedge, min_e);
            }

#pragma omp for schedule(static) ordered
            for (int ii = 0; ii < omp_get_num_threads(); ++ii)
            {
#pragma omp ordered
                {
                    surf_IDs.insert(surf_IDs.end(), local.begin(), local.end());
                }
            }

#pragma omp single
            {
                if (surf_IDs.size() != cells.nSurf)
                {
                    cout << "Mismatch of number of surface faces identified, and the number given"
                         << endl;
                    cout << "Identified: " << surf_IDs.size() << "  Given: " << cells.nSurf << endl;
                    exit(-1);
                }
            }

#pragma omp for schedule(static) nowait
            for (size_t ii = 0; ii < cells.nSurf; ++ii)
            {
                cells.smarkers[ii].first = surf_IDs[ii];
                cells.smarkers[ii].second = tau_markers[ii];
            }
        } /* End parallel */

        maxedge *= 5.0;

#ifdef DEBUG
        fprintf(dbout, "End of placing faces in elements.\n");
        fprintf(dbout, "Building boundary indexes.\n");
#endif

        cells.leftright = leftright;

#ifdef DEBUG
        fprintf(dbout, "All elements defined.\n\n");
#endif
    }

    void Read_tau_mesh_FACE(SIM& svar, MESH& cells, FLUID const& fvar, AERO const& avar)
    {
        string meshIn = svar.tau_mesh;
        string solIn = svar.tau_sol;

#ifdef DEBUG
        fprintf(dbout, "Attempting read of NetCDF file.\n");
        fprintf(dbout, "Mesh file: %s\n", meshIn.c_str());
        fprintf(dbout, "Solution file: %s\n", solIn.c_str());
#endif

        /*Read the mesh data*/
        int retval;
        int meshID;

        if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
        {
            cout << "A netCDF error occured whilst trying to open the mesh file:" << endl;
            cout << "\t" << meshIn << endl << endl;
            ERR(retval);
            exit(-1);
        }

        cout << "Mesh file open. Reading face data..." << endl;

        int hasTrig, hasQuad;
        int ptDimID, elemDimID, faceDimID, surfDimID, faceTDimID, faceQDimID, nPpTFDimID, nPpQFDimID;
        size_t nPnts, nElem, nFace, nSurf, nTFace, nQFace, nPpTFc, nPpQFc;

        // Retrieve how many elements there are.
        Get_Dim_ID(meshID, "no_of_elements", elemDimID);
        Get_Dim_Length(meshID, elemDimID, nElem);

        // Retrieve face number
        Get_Dim_ID(meshID, "no_of_faces", faceDimID);
        Get_Dim_Length(meshID, faceDimID, nFace);

        /* Retrieve how many points there are */
        Get_Dim_ID(meshID, "no_of_points", ptDimID);
        Get_Dim_Length(meshID, ptDimID, nPnts);

        Get_Dim_ID(meshID, "no_of_surfaceelements", surfDimID);
        Get_Dim_Length(meshID, surfDimID, nSurf);

        // Retrieve triangle face dimensions
        if ((retval = nc_inq_dimid(meshID, "no_of_triangles", &faceTDimID)))
        {
            cout << "No triangle faces" << endl;
            hasTrig = 0;
        }
        else
        {
            Get_Dim_Length(meshID, faceTDimID, nTFace);
            Get_Dim_ID(meshID, "points_per_triangle", nPpTFDimID);
            Get_Dim_Length(meshID, nPpTFDimID, nPpTFc);

            hasTrig = 1;
        }

        // Retrieve quadrilateral face dimensions
        if ((retval = nc_inq_dimid(meshID, "no_of_quadrilaterals", &faceQDimID)))
        {
            cout << "No quadrilateral faces" << endl;
            hasQuad = 0;
        }
        else
        {
            Get_Dim_Length(meshID, faceQDimID, nQFace);
            Get_Dim_ID(meshID, "points_per_quadrilateral", nPpQFDimID);
            Get_Dim_Length(meshID, nPpQFDimID, nPpQFc);

            hasQuad = 1;
        }

#ifdef DEBUG
        fprintf(dbout, "nElem: %zu nPnts: %zu nFace: %zu nSurf: %zu\n", nElem, nPnts, nFace, nSurf);
#endif

        cout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << " nSurf: " << nSurf
             << endl;

        cells.alloc(nPnts, nElem, nFace, nSurf);

        /*Get the faces of the mesh*/
        if (hasTrig)
        {
            Get_Element(meshID, "points_of_triangles", nTFace, nPpTFc, cells.faces);
        }

        if (hasQuad)
        {
            Get_Element(meshID, "points_of_quadrilaterals", nQFace, nPpQFc, cells.faces);
        }

        if (cells.faces.size() != cells.nFace)
        {
            cout << "Mismatch of number of faces to that defined." << endl;
            cout << "number of faces: " << cells.nFace << " faces size: " << cells.faces.size() << endl;
        }

        /*Get the coordinates of the mesh*/
        Get_Coordinates(meshID, nPnts, cells.verts);
        if (cells.verts.size() != nPnts)
        {
            cout << "Some data has been missed.\nPlease check how many points." << endl;
        }

        /*Adjust the scale*/
        cells.scale = svar.scale;
        for (auto& vert : cells.verts)
        {
            vert *= cells.scale;
        }

        /*Get face left and right, and put the faces in the elements*/
        Place_Faces(meshID, nFace, svar, cells);

        if ((retval = nc_close(meshID)))
        {
            cout << "A netCDF error occured whilst trying to close the mesh file:" << endl;
            cout << "\t" << meshIn << endl << endl;
            ERR(retval);
            exit(-1);
        }

#ifdef DEBUG
        fprintf(dbout, "End of interaction with mesh file and ingested data.\n\n");
#endif
    }
#endif
} // namespace TAU