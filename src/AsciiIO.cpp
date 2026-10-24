/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "AsciiIO.h"
#include "Var.h"
#include <regex>

#ifdef __STDC_ALLOC_LIB__
#define __STDC_WANT_LIB_EXT2__ 1
#else
#define _POSIX_C_SOURCE 200809L
#endif

inline std::string space2underscore(std::string text)
{
    std::replace(text.begin(), text.end(), ' ', '_');
    return text;
}
/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/

namespace ASCII
{
    void add_file_aux_data(FILE* fp, string const& name, string const& var)
    {
        fprintf(fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), var.c_str());
    }

    void add_file_aux_data(FILE* fp, std::string const& name, size_t const& var)
    {
        std::ostringstream stream;
        stream << var;
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str()
        );
    }

    void add_file_aux_data(FILE* fp, std::string const& name, int const& var)
    {
        std::ostringstream stream;
        stream << var;
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str()
        );
    }

    void add_file_aux_data(FILE* fp, std::string const& name, uint const& var)
    {
        std::ostringstream stream;
        stream << var;
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str()
        );
    }

    void add_file_aux_data(FILE* fp, std::string const& name, double const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var;
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str()
        );
    }

    void add_file_aux_data(FILE* fp, std::string const& name, float const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var;
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str()
        );
    }

    void add_file_aux_data(FILE* fp, std::string const& name, StateVecD const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var[0];
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name + " x").c_str(),
            stream.str().c_str()
        );

        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[1];
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name + " y").c_str(),
            stream.str().c_str()
        );

#if SIMDIM == 3
        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[2];
        fprintf(
            fp, "DATASETAUXDATA %s = \"%s\"\n", space2underscore(name + " z").c_str(),
            stream.str().c_str()
        );

#endif
    }

    void add_zone_aux_data(FILE* fp, string const& name, string const& var)
    {
        fprintf(fp, "AUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), var.c_str());
    }

    void add_zone_aux_data(FILE* fp, std::string const& name, size_t const& var)
    {
        std::ostringstream stream;
        stream << var;
        fprintf(fp, "AUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str());
    }

    void add_zone_aux_data(FILE* fp, std::string const& name, int const& var)
    {
        std::ostringstream stream;
        stream << var;
        fprintf(fp, "AUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str());
    }

    void add_zone_aux_data(FILE* fp, std::string const& name, uint const& var)
    {
        std::ostringstream stream;
        stream << var;
        fprintf(fp, "AUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str());
    }

    void add_zone_aux_data(FILE* fp, std::string const& name, double const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var;
        fprintf(fp, "AUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str());
    }

    void add_zone_aux_data(FILE* fp, std::string const& name, float const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var;
        fprintf(fp, "AUXDATA %s = \"%s\"\n", space2underscore(name).c_str(), stream.str().c_str());
    }

    void add_zone_aux_data(FILE* fp, std::string const& name, StateVecD const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var[0];
        fprintf(
            fp, "AUXDATA %s = \"%s\"\n", space2underscore(name + " x").c_str(), stream.str().c_str()
        );

        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[1];
        fprintf(
            fp, "AUXDATA %s = \"%s\"\n", space2underscore(name + " y").c_str(), stream.str().c_str()
        );

#if SIMDIM == 3
        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[2];
        fprintf(
            fp, "AUXDATA %s = \"%s\"\n", space2underscore(name + " z").c_str(), stream.str().c_str()
        );

#endif
    }

    void Write_Vector(FILE* fp, vector<real> const& vec)
    {
        size_t newl = 0;
        size_t const end = vec.size();
        for (size_t ii = 0; ii < end; ++ii)
        {
            fprintf(fp, "%3.7e", vec[ii]);
            if (newl > 4)
            {
                fprintf(fp, "\n");
                newl = 0;
            }
            else
                newl++;
        }
        fprintf(fp, "\n");
    }

    void Write_Vector(FILE* fp, vector<uint> const& vec)
    {
        size_t newl = 0;
        size_t const end = vec.size();
        for (size_t ii = 0; ii < end; ++ii)
        {
            fprintf(fp, "%u", vec[ii]);
            if (newl > 6)
            {
                fprintf(fp, "\n");
                newl = 0;
            }
            else
                newl++;
        }
        fprintf(fp, "\n");
    }

    void Write_Vector(FILE* fp, vector<size_t> const& vec)
    {
        size_t newl = 0;
        size_t const end = vec.size();
        for (size_t ii = 0; ii < end; ++ii)
        {
            fprintf(fp, "%zu", vec[ii]);
            if (newl > 6)
            {
                fprintf(fp, "\n");
                newl = 0;
            }
            else
                newl++;
        }
        fprintf(fp, "\n");
    }

    void Write_Vector(FILE* fp, vector<long> const& vec)
    {
        size_t newl = 0;
        size_t const end = vec.size();
        for (size_t ii = 0; ii < end; ++ii)
        {
            fprintf(fp, "%zu", vec[ii]);
            if (newl > 6)
            {
                fprintf(fp, "\n");
                newl = 0;
            }
            else
                newl++;
        }
        fprintf(fp, "\n");
    }

    void Write_Vector(FILE* fp, vector<int> const& vec)
    {
        size_t newl = 0;
        size_t const end = vec.size();
        for (size_t ii = 0; ii < end; ++ii)
        {
            fprintf(fp, "%d", vec[ii]);
            if (newl > 6)
            {
                fprintf(fp, "\n");
                newl = 0;
            }
            else
                newl++;
        }
        fprintf(fp, "\n");
    }

    void Write_Vector(FILE* fp, vector<StateVecD> const& vec)
    {
        size_t const end = vec.size();
        for (size_t ii = 0; ii < end; ++ii)
        {
            for (uint dim = 0; dim < SIMDIM; ++dim)
                fprintf(fp, "%3.7e", vec[ii][dim]);
            fprintf(fp, "\n");
        }
    }

    void Write_File_header(FILE* fp, string const& var_names, string const& title)
    {
        std::regex r(",");
        string var = std::regex_replace(var_names, r, "\",\"");
        fprintf(fp, "TITLE=\"%s\"\n", title.c_str());
        fprintf(fp, "VARIABLES=\"%s\"", var.c_str());
    }

    void Write_Zone_Header(
        FILE* fp, std::string zone_header, int64_t size, int32_t strand_id, real solution_time
    )
    {
        fprintf(
            fp, "ZONE T=\"%s\", I=%zu, DATAPACKING=BLOCK, STRANDID=%d, SOLUTIONTIME=%.7g\n", zone_header,
            size, strand_id, solution_time
        );
    }

} // namespace ASCII

void Write_Cell_Data(MESH const& cdata)
{
#ifdef DEBUG
    fprintf(dbout, "Entering Write_Cell_Data...\n");
#endif

    cout << "Writing cell based data." << endl;

    std::ofstream fout("Cell.dat", std::ios::out);
    if (!fout.is_open())
    {
        cout << "Failed to open data file for writing mesh." << endl;
        exit(-1);
    }

    fout << "TITLE = \"3D Mesh Solution\"\n";
    fout << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\"\n";
    fout << "ZONE T=\"Cell Data\"" << endl;
    fout << "N=" << cdata.nPnts << ", E=" << cdata.nElem << ", F=FEBLOCK, ET=BRICK" << endl << endl;

    /*Write vertices*/
    fout << std::left << std::scientific << std::setprecision(6);
    fout << std::setw(1);
    for (uint ii = 0; ii < SIMDIM; ++ii)
    {
        uint kk = 0;
        for (uint jj = 0; jj < cdata.verts.size(); ++jj)
        {
            fout << std::setw(15) << cdata.verts[jj][ii];
            kk++;

            if (kk == 5)
            {
                fout << endl;
                fout << std::setw(1);
                kk = 0;
            }
        }

        if (kk % 5 != 0)
            fout << "\n";
    }

    /*Write element indexing*/
    fout << std::fixed;
    for (uint index = 0; index < cdata.cFaces.size(); ++index)
    {
        /* Build the element unique index list */
        vector<size_t> elem;
        for (auto const& faceID : cdata.cFaces[index])
        {
            vector<size_t> const face = cdata.faces[faceID];
            for (auto const& vert : face)
            {
                if (std::find(elem.begin(), elem.end(), vert) == elem.end())
                { /*Vertex doesn't exist in the elems vector yet.*/
#pragma omp critical
                    {
                        elem.emplace_back(vert);
                    }
                }
            }
        }

        for (auto const& vert : elem)
        {
            fout << std::setw(6) << vert + 1;
        }
        fout << "\n";
    }

    fout.close();

#ifdef DEBUG
    fprintf(dbout, "Exiting Write_Cell_Data...\n");
#endif
}

void Write_Face_Data(MESH const& cells)
{

    std::ofstream f1("foam_mesh.dat", std::ios::out);
    f1 << "VARIABLES= \"X\", \"Y\", \"Z\"" << endl;
    uint w = 17;
    f1 << std::left << std::scientific << std::setprecision(8);

    /*Write zone header information*/
    f1 << "ZONE T=\"OpenFOAM MESH\"" << endl;
    f1 << "ZONETYPE=FEPOLYHEDRON" << endl;
    f1 << "NODES=" << cells.nPnts << " ELEMENTS=" << cells.nElem << " FACES=" << cells.nFace << endl;
    size_t TotalNumFaceNodes = cells.nFace * 3;
    f1 << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
    f1 << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

    /*Write vertices in block format*/
    size_t n = 0;

    for (size_t dim = 0; dim < SIMDIM; ++dim)
    {
        for (auto const& pnt : cells.verts)
        {
            f1 << std::setw(w) << pnt[dim];
            n++;

            if (n == 6)
            {
                f1 << endl;
                n = 0;
            }
        }
        f1 << endl;
    }

    /*Write how many vertices per face*/
    n = 0;
    for (size_t ii = 0; ii < cells.nFace; ++ii)
    {
        f1 << std::setw(5) << 3;
        n++;

        if (n == 6)
        {
            f1 << endl;
            n = 0;
        }
    }
    f1 << endl;

    /*Write the face vertex list*/
    n = 0;
    w = 10;
    for (auto const& face : cells.faces)
    {
        for (auto const& vertex : face)
        {
            f1 << std::setw(w) << vertex + 1;
            n++;

            if (n == 6)
            {
                f1 << endl;
                n = 0;
            }
        }
    }
    f1 << endl;

    /*Write left elements*/
    n = 0;
    for (auto const& lr : cells.leftright)
    {
        f1 << std::setw(w) << lr.first + 1;
        n++;

        if (n == 6)
        {
            f1 << endl;
            n = 0;
        }
    }
    f1 << endl;

    /*Write right elements*/
    n = 0;
    for (auto const& lr : cells.leftright)
    {
        if (lr.second < 0)
        {
            f1 << std::setw(w) << 0;
        }
        else
        {
            f1 << std::setw(w) << lr.second + 1;
        }
        n++;

        if (n == 6)
        {
            f1 << endl;
            n = 0;
        }
    }
    f1 << endl << endl;
    f1.close();
    // exit(0);
}
