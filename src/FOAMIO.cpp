/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Var.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
using std::fstream;
using std::string;

/* CURRENTLY 3D ONLY */
#if SIMDIM == 3
namespace FOAM
{
    /****************************************************/
    /**************** ASCII Functions *******************/
    /****************************************************/
    namespace ascii
    {
        /* General function to read the label data from an OpenFOAM ascii file */
        void Read_Label_Data(ifstream& fin, size_t const& nFaces, vector<int>& label, size_t& nCells)
        {
            nCells = 0;
            string line;
            getline(fin, line);
            label = vector<int>(nFaces);
            for (size_t ii = 0; ii < nFaces; ++ii)
            {
                getline(fin, line);

                std::istringstream sline(line);

                sline >> label[ii];

                if (label[ii] + 1 > int(nCells + 1))
                {
                    nCells = label[ii] + 1;
                }
            }
        }

        /* General function to read the scalar data from an OpenFOAM ascii file */
        void Read_Scalar_Data(ifstream& fin, size_t const& nPnts, vector<real>& field)
        {
            string line;
            getline(fin, line);
            field = vector<real>(nPnts);
            for (size_t ii = 0; ii < nPnts; ++ii)
            {
                getline(fin, line);
                std::istringstream iss(line);

                iss >> field[ii];
            }
        }

        /* General function to read the vector data from an OpenFOAM ascii file */
        void
        Read_Vector_Data(ifstream& fin, size_t const& nPnts, vector<Eigen::Matrix<real, 3, 1>>& field)
        {
            string line;
            getline(fin, line);
            field = vector<Eigen::Matrix<real, 3, 1>>(nPnts);
            for (size_t ii = 0; ii < nPnts; ++ii)
            {
                getline(fin, line);
                line.erase(std::remove(line.begin(), line.end(), '('), line.end());
                line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

                std::istringstream iss(line);
                Eigen::Matrix<real, 3, 1> vec;
                iss >> vec(0);
                iss >> vec(1);
                iss >> vec(2);

                field[ii] = vec;
            }
        }

        /* General function to read the face data from an OpenFOAM ascii file */
        void Read_Face_Data(ifstream& fin, size_t const& nFaces, vector<vector<size_t>>& faces)
        {
            string line;
            getline(fin, line);
            faces = vector<vector<size_t>>(nFaces, vector<size_t>());
            for (size_t ii = 0; ii < nFaces; ++ii)
            {
                getline(fin, line);

                std::istringstream sline(line);

                size_t nPnts;
                sline >> nPnts;
                size_t left = line.find("(");
                size_t right = line.find(")");

                string temp = line.substr(left + 1, right - left);

                std::istringstream iss(temp);
                vector<size_t> face(nPnts);
                for (size_t jj = 0; jj < nPnts; jj++)
                {
                    iss >> face[jj];
                }
                faces[ii] = face;
            }
        }

    } // namespace ascii

    // Binary functions
    namespace binary
    {
        /* General function to read the label data from an OpenFOAM binary file */
        void Read_Label_Data(
            std::ifstream& fin, int const& foam_label_size, int const& foam_scalar_size,
            size_t const& nFaces, vector<int>& field, size_t& nCells
        )
        {
            /* Read binary stream */
            field = vector<int>(nFaces);
            if (foam_label_size == 32)
            {
                for (size_t ii = 0; ii < nFaces; ++ii)
                {
                    int32_t a;

                    fin.read(reinterpret_cast<char*>(&a), sizeof(a));

                    field[ii] = a;
                    if (static_cast<size_t>(a + 1) > nCells + 1)
                    {
                        nCells = static_cast<size_t>(a + 1);
                    }
                }
            }
            else if (foam_label_size == 64)
            {
                for (size_t ii = 0; ii < nFaces; ++ii)
                {
                    int64_t a;

                    fin.read(reinterpret_cast<char*>(&a), sizeof(a));

                    field[ii] = a;

                    if (static_cast<size_t>(a + 1) > nCells + 1)
                    {
                        nCells = static_cast<size_t>(a + 1);
                    }
                }
            }
        }

        /* General function to read the scalar data from an OpenFOAM binary file */
        void Read_Scalar_Data(
            std::ifstream& fin, int const& foam_label_size, int const& foam_scalar_size,
            size_t const& nPnts, vector<real>& field
        )
        {
            /* Read binary stream */
            field = vector<real>(nPnts);

            if (foam_scalar_size == 32)
            {
                for (size_t ii = 0; ii < nPnts; ++ii)
                {
                    float a;

                    fin.read(reinterpret_cast<char*>(&a), sizeof(a));

                    field[ii] = a;
                }
            }
            else if (foam_scalar_size == 64)
            {
                for (size_t ii = 0; ii < nPnts; ++ii)
                {
                    double a;

                    fin.read(reinterpret_cast<char*>(&a), sizeof(a));

                    field[ii] = a;
                }
            }
        }

        /* General function to read the vector data from an OpenFOAM binary file */
        void Read_Vector_Data(
            std::ifstream& fin, int const& foam_label_size, int const& foam_scalar_size,
            size_t const& nPnts, vector<Eigen::Matrix<real, 3, 1>>& field
        )
        {
            field = vector<Eigen::Matrix<real, 3, 1>>(nPnts);

            if (foam_scalar_size == 32)
            {
                for (size_t ii = 0; ii < nPnts; ++ii)
                {
                    float a;

                    for (size_t jj = 0; jj < 3; ++jj)
                    {
                        fin.read(reinterpret_cast<char*>(&a), sizeof(a));
                        field[ii](jj) = static_cast<real>(a);
                    }
                }
            }
            else if (foam_scalar_size == 64)
            {
                for (size_t ii = 0; ii < nPnts; ++ii)
                {
                    double a;

                    for (size_t jj = 0; jj < 3; ++jj)
                    {
                        fin.read(reinterpret_cast<char*>(&a), sizeof(a));
                        field[ii](jj) = static_cast<real>(a);
                    }
                }
            }
        }

        /* General function to read the face data from an OpenFOAM binary file */
        void Read_Face_Data(
            std::ifstream& fin, int const& foam_label_size, int const& foam_scalar_size,
            size_t const& nFaces, vector<vector<size_t>>& faces
        )
        {
            faces = vector<vector<size_t>>(nFaces, vector<size_t>());

            if (foam_label_size == 32)
            {
                /* Read indexes */
                vector<int32_t> index(nFaces + 1);
                for (size_t ii = 0; ii < nFaces + 1; ++ii)
                {
                    int32_t a;

                    fin.read(reinterpret_cast<char*>(&a), sizeof(a));

                    index[ii] = a;
                }

                char temp = '0';
                int nPnts = 0;
                char bracket = '(';
                vector<char> interim;

                /* Read the interim data between the arrays */
                while (temp != bracket)
                {
                    fin.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    interim.emplace_back(temp);
                }

                /* Put the vector into a string to extract the number */
                string str(interim.begin(), interim.end());
                /* Clean out characters other than the number */
                str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
                str.erase(std::remove(str.begin(), str.end(), '('), str.end());
                str.erase(std::remove(str.begin(), str.end(), ')'), str.end());

                nPnts = std::stoi(str);

                int count = 0;
                for (size_t ii = 0; ii < nFaces; ++ii)
                {
                    int32_t length = index[ii + 1] - index[ii];
                    int32_t ind;

                    for (int jj = 0; jj < length; ++jj)
                    {
                        fin.read(reinterpret_cast<char*>(&ind), sizeof(ind));
                        faces[ii].emplace_back(static_cast<size_t>(ind));

                        if (count > nPnts)
                        {
                            cout << "Warning: Exceeded number of specified for face array" << endl;
                        }
                        count++;
                    }
                }
            }
            else if (foam_label_size == 64)
            {
                /* Read indexes */
                vector<int64_t> index(nFaces + 1);
                for (size_t ii = 0; ii < nFaces + 1; ++ii)
                {
                    int64_t a;

                    fin.read(reinterpret_cast<char*>(&a), sizeof(a));

                    index[ii] = a;
                }

                char temp = '0';
                int nPnts = 0;
                char bracket = '(';
                vector<char> interim;

                /* Read the interim data between the arrays */
                while (temp != bracket)
                {
                    fin.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    interim.emplace_back(temp);
                }

                /* Put the vector into a string to extract the number */
                string str(interim.begin(), interim.end());
                /* Clean out characters other than the number */
                str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
                str.erase(std::remove(str.begin(), str.end(), '('), str.end());
                str.erase(std::remove(str.begin(), str.end(), ')'), str.end());

                nPnts = std::stoi(str);

                int count = 0;
                for (size_t ii = 0; ii < nFaces; ++ii)
                {
                    int64_t length = index[ii + 1] - index[ii];

                    int64_t ind;

                    for (int jj = 0; jj < length; ++jj)
                    {
                        fin.read(reinterpret_cast<char*>(&ind), sizeof(ind));
                        faces[ii].emplace_back(static_cast<size_t>(ind));

                        if (count > nPnts)
                        {
                            cout << "Warning: Exceeded number of specified for face array" << endl;
                        }
                        count++;
                    }
                }
            }
        }

    } // namespace binary

    /* Read the header section of an openfoam file, and identify the file class and binary information.
     */
    void Read_Header(
        std::ifstream& fin, string& class_, int& binary, int& foam_label_size, int& foam_scalar_size
    )
    {
        /* Read the header info */
        /* Information that is wanted: format and class */
        string line;

        while (line.find("}") == string::npos)
        {
            std::getline(fin, line);

            if (line.find("format") != string::npos)
            {
                if (line.find("ascii") != string::npos)
                { /* File format doesn't matter for the boundary file, but suggests  */
                    // cout << "File format is in ascii..." << endl;
                    binary = 0;
                }
                else
                {
                    // cout << "File format is in binary..." << endl;
                    binary = 1;
                }
            }

            if (line.find("class") != string::npos)
            {
                std::istringstream iss(line);
                string temp;
                iss >> temp;
                iss >> class_;

                string b = &class_.back();
                if (b == ";")
                {
                    class_.pop_back();
                }
            }

            if (line.find("arch") != string::npos)
            {
                size_t p1 = line.find("label");

                string lsize = line.substr(p1 + 6, 2);

                std::istringstream iss(lsize);
                iss >> foam_label_size;

                size_t p2 = line.find("scalar");

                string ssize = line.substr(p2 + 7, 2);

                iss = std::istringstream(ssize);
                iss >> foam_scalar_size;
            }
        }
    }

    /* Read information for a patch in the boundary file, and whether it should be marked as a wall */
    void Read_Patch(ifstream& fin, size_t& nFaces, size_t& startFace, int& wall)
    {
        /* Optionally can get the name? */
        string name;
        getline(fin, name);
        string line;

        while (line.find("}") == string::npos)
        {
            getline(fin, line);

            if (line.find("type") != string::npos)
            {
                if (line.find("wall") != string::npos)
                    wall = 1;
                else
                    wall = 0;
            }
            else if (line.find("nFaces") != string::npos)
            {
                std::istringstream iss(line);
                string temp;
                iss >> temp;
                iss >> nFaces;
            }
            else if (line.find("startFace") != string::npos)
            {
                std::istringstream iss(line);
                string temp;
                iss >> temp;
                iss >> startFace;
            }
        }
    }

    /* Read the standard preamble for an openfoam file, including header, and number of data points in
     * array */
    void Read_Preamble(
        ifstream& fin, string const& file, string const exp_class, int& binary, int& foam_label_size,
        int& foam_scalar_size, size_t& nVals
    )
    {
        string line;

        while (line.find("FoamFile") == string::npos)
        { /* Ignore banner. Not important */
            getline(fin, line);
        }

        string class_;
        Read_Header(fin, class_, binary, foam_label_size, foam_scalar_size);

        if (class_ != exp_class)
        {
            cout << "File " << file << " is the expected class." << endl;
            cout << "File is class \"" << class_ << "\" and should be \"" << exp_class << "\"" << endl;
            cout << "Please check your input" << endl;
            exit(-1);
        }

        /* Skip to the size declaration */
        getline(fin, line);
        while (line.find("//") == 0 || line == "" || line.find("dimensions") != string::npos)
        {
            getline(fin, line);
        }

        /* Skip line if this is not a polyMesh file */
        if (line.find("internalField") != string::npos)
        {
            getline(fin, line);
        }

        std::istringstream iss(line);
        iss >> nVals;
    }

    /* Read the boundary file of the mesh, and find out where the boundary faces begin. */
    /* Boundary is in ascii no matter what. Can use to determine if system is in binary or not */
    void Read_Boundary(
        SIM& svar, vector<std::pair<size_t, size_t>>& surfInfo, vector<int>& walls, int& binary
    )
    {
        // Open folder
        cout << "Reading boundary file..." << endl;
        string file = svar.foam_dir;
        file.append("/constant/polyMesh/boundary");
        std::ifstream fin(file);

        if (!fin.is_open())
        {
            cout << "Failed to open boundary file" << endl;
            exit(-1);
        }

        string line;

        while (line.find("FoamFile") == string::npos)
        { /* Ignore banner. Not important */
            std::getline(fin, line);
        }

        string class_;
        int foam_label_size, foam_scalar_size;
        Read_Header(fin, class_, binary, foam_label_size, foam_scalar_size);

        getline(fin, line);
        while (line == "" || line.find("//") != string::npos)
        { /* Skip through to things begin */
            getline(fin, line);
        }

        /* Get how many surfaces to expect */
        size_t nSurfaces;
        std::istringstream iss(line);
        iss >> nSurfaces;

        surfInfo = vector<std::pair<size_t, size_t>>(nSurfaces);
        walls = vector<int>(nSurfaces);
        for (size_t ii = 0; ii < nSurfaces; ++ii)
        {
            // cout << "Reading patch: " << ii << endl;
            size_t nFaces, startFace;
            int isWall;
            Read_Patch(fin, nFaces, startFace, isWall);
            surfInfo[ii] = (std::pair<size_t, size_t>(nFaces, startFace));
            walls[ii] = isWall;
        }

        fin.close();
    }

    void Post_Process(
        vector<std::pair<size_t, size_t>> const& surfBounds, vector<int> const& walls,
        vector<StateVecD> const& pnts, vector<vector<size_t>> const& faces_, vector<int> const& left_,
        vector<int>& right_, size_t const& nCells, MESH& cells
    )
    {
        /* Now check, and fill in faces */
        if (right_.size() != left_.size())
        {
            // /* Sum the wall surfaces, and see if that adds up to the missing faces */
            // size_t internalSum = 0;
            // size_t externalSum = 0;
            // for(size_t ii = 0; ii < walls.size(); ++ii)
            // {
            //     if(walls[ii] == 1)
            //     {
            //         /* Means it's an internal wall */
            //         internalSum += surfBounds[ii].first;
            //     }
            //     else
            //     {
            //         /* Means it's an external patch */
            //         externalSum += surfBounds[ii].first;

            //     }
            // }

            // size_t diff = left_.size() - right_.size();
            // cout << "Face difference: " << diff << "  internal faces: " << internalSum << "  external
            // faces: " << externalSum << endl;

            for (size_t ii = 0; ii < walls.size(); ++ii)
            {
                if (walls[ii] == 1)
                {
                    /* Means it's an internal wall */
                    vector<int> temp(surfBounds[ii].first, -1);
                    right_.insert(right_.end(), temp.begin(), temp.end());
                }
                else
                {
                    /* Means it's an external patch */
                    vector<int> temp(surfBounds[ii].first, -2);
                    right_.insert(right_.end(), temp.begin(), temp.end());
                }
            }
        }

        if (left_.size() != faces_.size())
        {
            cout << "Mismatch of number of faces and face owner sizes." << endl;
            cout << "Number of faces: " << faces_.size() << "Owner size: " << left_.size() << endl;
        }

        if (right_.size() != faces_.size())
        {
            cout << "Mismatch of number of faces and face neighbour sizes." << endl;
            cout << "Number of faces: " << faces_.size() << "  Neighbour size: " << right_.size()
                 << endl;
        }

        if (right_.size() != left_.size())
        {
            cout << "Mismatch of cell owner and neighbour sizes." << endl;
            cout << "Owner size: " << left_.size() << "  Neighbour size: " << right_.size() << endl;
        }

        /* need to break faces of four points*/
        cout << "Splitting faces... " << endl;
        vector<vector<size_t>> faces;
        vector<std::pair<int, int>> leftright;
        for (size_t ii = 0; ii < faces_.size(); ++ii)
        {
            if (faces_[ii].size() > 3)
            {
                for (size_t jj = 0; jj < faces_[ii].size() - 2; ++jj)
                {
                    vector<size_t> face = {faces_[ii][0], faces_[ii][jj + 1], faces_[ii][jj + 2]};

                    faces.emplace_back(face);
                    leftright.emplace_back(std::pair<int, int>(left_[ii], right_[ii]));
                }
            }
            else
            {
                faces.emplace_back(faces_[ii]);
                leftright.emplace_back(std::pair<int, int>(left_[ii], right_[ii]));
            }
        }

        size_t nFace = faces.size();

        cout << "Placing faces..." << endl;
        vector<vector<size_t>> cFaces(nCells);

        for (size_t ii = 0; ii < nFace; ++ii)
        {

            cFaces[leftright[ii].first].emplace_back(ii);
            if (leftright[ii].second >= 0)
                cFaces[leftright[ii].second].emplace_back(ii);
        }

        /* Find cell centres  */
        cout << "Finding cell centres..." << endl;
        vector<StateVecD> cCentre(nCells);
        // vector<vector<size_t>> elems(nCells);
        for (size_t ii = 0; ii < nCells; ++ii)
        {
            vector<size_t> verts;
            for (auto const& face : cFaces[ii])
            {
                verts.insert(verts.end(), faces[face].begin(), faces[face].end());
            }

            /* Sort and delete repetitions to average values */
            std::sort(verts.begin(), verts.end());
            std::unique(verts.begin(), verts.end());

            /* Average the point values */
            StateVecD sum = StateVecD::Zero();
            for (auto jj : verts)
            {
                sum += pnts[jj];
            }
            sum /= verts.size();
            cCentre[ii] = sum;
            // elems[ii] = verts;
        }

        cout << "Placing data in cell structure..." << endl;
        /* Place into the mesh structure. */
        cells.verts = pnts; /* Point data */

        cells.faces = faces; /* Face data */
        cells.leftright = leftright;

        // cells.elems = elems; /* Cell data */
        cells.cCentre = cCentre;
        cells.cFaces = cFaces;

        cells.nElem = nCells;
        cells.nPnts = pnts.size();
        cells.nFace = faces.size();
    }

    /* General function to read the scalar data from an OpenFOAM ascii file */
    void Read_Label_Field(string& file, vector<int>& field, size_t& nCells)
    {
        std::ifstream fin(file);

        int binary, foam_label_size, foam_scalar_size;
        size_t nFaces;

        Read_Preamble(fin, file, "labelList", binary, foam_label_size, foam_scalar_size, nFaces);

        cout << "Number of faces: " << nFaces << endl;

        if (binary == 0)
        {
            ascii::Read_Label_Data(fin, nFaces, field, nCells);
        }
        else
        {
            int pos = fin.tellg();
            fin.close();
            fin.open(file, std::ifstream::binary);
            fin.seekg(pos + 1);
            binary::Read_Label_Data(fin, foam_label_size, foam_scalar_size, nFaces, field, nCells);
        }

        fin.close();
    }

    /* General function to read the scalar data from an OpenFOAM ascii file */
    void Read_Solution_Scalar(string& file, vector<real>& field)
    {
        std::ifstream fin(file, std::ios::in);

        int binary, foam_label_size, foam_scalar_size;
        size_t nInternal;

        Read_Preamble(fin, file, "volScalarField", binary, foam_label_size, foam_scalar_size, nInternal);

        cout << "Number of cells: " << nInternal << endl;

        if (binary == 0)
        {
            ascii::Read_Scalar_Data(fin, nInternal, field);
        }
        else
        {
            int pos = fin.tellg();
            fin.close();
            fin.open(file, std::ifstream::binary);
            fin.seekg(pos + 1);
            binary::Read_Scalar_Data(fin, foam_label_size, foam_scalar_size, nInternal, field);
        }

        fin.close();
    }

    /* General function to read the vector data from an OpenFOAM ascii file */
    void Read_Solution_Vector(string& file, vector<Eigen::Matrix<real, 3, 1>>& field)
    {
        std::ifstream fin(file, std::ios::in);

        int binary, foam_label_size, foam_scalar_size;
        size_t nInternal;

        Read_Preamble(fin, file, "volVectorField", binary, foam_label_size, foam_scalar_size, nInternal);

        cout << "Number of cells: " << nInternal << endl;

        if (binary == 0)
        {
            ascii::Read_Vector_Data(fin, nInternal, field);
        }
        else
        {
            int pos = fin.tellg();
            fin.close();
            fin.open(file, std::ifstream::binary);
            fin.seekg(pos + 1);
            binary::Read_Vector_Data(fin, foam_label_size, foam_scalar_size, nInternal, field);
        }

        fin.close();
    }

    /* General function to read the vector data from an OpenFOAM ascii file */
    void Read_Points(string& file, vector<Eigen::Matrix<real, 3, 1>>& pnts)
    {
        std::ifstream fin(file);

        int binary, foam_label_size, foam_scalar_size;
        size_t nPnts;

        Read_Preamble(fin, file, "vectorField", binary, foam_label_size, foam_scalar_size, nPnts);

        cout << "Number of points: " << nPnts << endl;

        if (binary == 0)
        {
            ascii::Read_Vector_Data(fin, nPnts, pnts);
        }
        else
        {
            int pos = fin.tellg();
            fin.close();
            fin.open(file, std::ifstream::binary);
            fin.seekg(pos + 1);
            binary::Read_Vector_Data(fin, foam_label_size, foam_scalar_size, nPnts, pnts);
        }

        fin.close();
    }

    /* Function to read the face data from an OpenFOAM binary file */
    void Read_Faces(string& file, vector<vector<size_t>>& faces)
    {
        std::ifstream fin(file, std::ios::in);

        string line;
        getline(fin, line); /* Get first line */

        while (line.find("FoamFile") == string::npos)
        { /* Ignore banner. Not important */
            getline(fin, line);
        }

        string class_;
        int binary, foam_label_size, foam_scalar_size;

        Read_Header(fin, class_, binary, foam_label_size, foam_scalar_size);

        if (binary == 0)
        {
            if (class_ != "faceList")
            {
                cout << "File " << file << " is not the correct class." << endl;
                cout << "File is class \"" << class_ << "\" and should be \"faceList\"" << endl;
                cout << "Please check your input" << endl;
                exit(-1);
            }
        }
        else
        {
            if (class_ != "faceCompactList")
            {
                cout << "File " << file << " is not the correct class." << endl;
                cout << "File is class \"" << class_ << "\" and should be \"faceCompactList\"" << endl;
                cout << "Please check your input" << endl;
                exit(-1);
            }
        }

        getline(fin, line);
        while (line == "" || line.find("//") == 0)
        { /* Skip through to things begin */
            getline(fin, line);
        }

        /* Get how many faces to expect */
        size_t nFaces;
        std::istringstream iss(line);
        iss >> nFaces;

        if (binary == 0)
        {
            cout << "Number of faces: " << nFaces << endl;
            ascii::Read_Face_Data(fin, nFaces, faces);
        }
        else
        {
            /* Will be the number of faces plus 1. */
            nFaces--;
            cout << "Number of faces: " << nFaces << endl;

            int pos = fin.tellg();
            fin.close();
            fin.open(file, std::ifstream::binary);
            fin.seekg(pos + 1);
            binary::Read_Face_Data(fin, foam_label_size, foam_scalar_size, nFaces, faces);
        }

        fin.close();
    }

    void Read_polyMesh(
        SIM const& svar, vector<std::pair<size_t, size_t>> const& surfBounds, vector<int> const& walls,
        MESH& cells
    )
    {
        string file = svar.foam_dir;

        /* Get point values */
        cout << "Reading point data..." << endl;
        file.append("/constant/polyMesh/points");
        vector<StateVecD> pnts;
        Read_Points(file, pnts);

        /* Get face vertex data */
        cout << "Reading face data..." << endl;
        vector<vector<size_t>> faces;
        file = svar.foam_dir;
        file.append("/constant/polyMesh/faces");
        Read_Faces(file, faces);

        /* Get owner data. Considered as cell left (normal direction unimportant) */
        cout << "Reading owner data..." << endl;
        file = svar.foam_dir;
        file.append("/constant/polyMesh/owner");
        vector<int> fOwner;
        size_t nCells = 0;
        Read_Label_Field(file, fOwner, nCells);

        /* Read neighbour data. This may well not have the same number of faces. */
        cout << "Reading neighbour data..." << endl;
        file = svar.foam_dir;
        file.append("/constant/polyMesh/neighbour");
        vector<int> fNeigh; /* Neighbour file */
        Read_Label_Field(file, fNeigh, nCells);

        Post_Process(surfBounds, walls, pnts, faces, fOwner, fNeigh, nCells, cells);
    }

    void Read_Solution(SIM& svar, MESH& cells)
    {
        string timef = svar.foam_dir;
        timef.append("/");
        timef.append(svar.foam_sol);

        cout << "Reading pressure data..." << endl;
        string file = timef;
        if (svar.foam_buoyant_sim == 0)
            file.append("/p");
        else
            file.append("/p_rgh");

        vector<real> press;
        Read_Solution_Scalar(file, press);

        cout << "Reading velocity data..." << endl;
        file = timef;
        file.append("/U");
        vector<StateVecD> vel;
        Read_Solution_Vector(file, vel);

        if (press.size() != cells.nElem)
        {
            cout << "Mismatch between solution size and mesh size" << endl;
            cout << "Pressure vector size: " << press.size() << "  Mesh size: " << cells.nElem << endl;
        }

        if (vel.size() != cells.nElem)
        {
            cout << "Mismatch between solution size and mesh size" << endl;
            cout << "Velocity vector size: " << vel.size() << "  Mesh size: " << cells.nElem << endl;
        }

        cells.cVel = vel;
        cells.cP = press;
    }

    void Read_FOAM(SIM& svar, MESH& cells)
    {
        /* Read boundary values, and find where wall faces begin.*/
        vector<std::pair<size_t, size_t>> surfBounds;
        vector<int> walls;
        int binary;
        Read_Boundary(svar, surfBounds, walls, binary);

        Read_polyMesh(svar, surfBounds, walls, cells);

        Read_Solution(svar, cells);
    }
} // namespace FOAM
#endif