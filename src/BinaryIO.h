/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef BINARYIO_H
#define BINARYIO_H

#include "Var.h"

namespace SZLIO
{
    bool open_file(
        string const& file_name, string const& zone_name, string const& var_names, void*& file_handle
    );

    int32_t create_zone(
        void* const& file_handle, int32_t const& strandID, std::string const& zone_header,
        vector<int32_t> const& var_types, int64_t const& size, double const& solution_time
    );

    void close_file(void* handle);

    void flush_file(void* handle);

    void add_file_aux_data(void* const& file, std::string const& name, std::string const& var);
    void add_file_aux_data(void* const& file, std::string const& name, size_t const& var);
    void add_file_aux_data(void* const& file, std::string const& name, int const& var);
    void add_file_aux_data(void* const& file, std::string const& name, uint const& var);
    void add_file_aux_data(void* const& file, std::string const& name, double const& var);
    void add_file_aux_data(void* const& file, std::string const& name, StateVecD const& var);

    void add_zone_aux_data(
        void* const& file, int32_t const& zone, std::string const& name, std::string const& var
    );
    void add_zone_aux_data(
        void* const& file, int32_t const& zone, std::string const& name, size_t const& var
    );
    void
    add_zone_aux_data(void* const& file, int32_t const& zone, std::string const& name, int const& var);
    void add_zone_aux_data(
        void* const& file, int32_t const& zone, std::string const& name, double const& var
    );
    void add_zone_aux_data(
        void* const& file, int32_t const& zone, std::string const& name, StateVecD const& var
    );

    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<real> const& varVec, string const& varName
    );

    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<int32_t> const& varVec, string const& varName
    );
    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<uint8_t> const& varVec, string const& varName
    );

    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<StateVecD> const& varVec, string const& varName
    );

}; // namespace SZLIO

inline int Combine_SZPLT(string const& file)
{
    string cmd = "exec szcombine \"";
    cmd.append(file);
    cmd.append("\"");

#ifdef DEBUG
    fprintf(dbout, "Attempting to combine szplt.\n");
    fprintf(dbout, "Command: %s\n", cmd.c_str());
#endif

    printf("Combining szplt: %s\n", file.c_str());
    if (system(cmd.c_str()))
    {
        printf("Could not combine szplt file.\n");
        printf("Command: %s\n", cmd.c_str());
        return -1;
    }
    return 0;
}

// IPT Binary functions (Allows more explicit library includes, so that tecio is only in BinaryIO.cpp.)
namespace IPT
{
    namespace BINARY
    {
        void Init_IPT_Files(SIM& svar);

        void Write_Point(SIM const& svar, IPTPart const& pnp1);

        void
        Write_State(SIM const& svar, IPTState const& pnp1, string const& zoneName, void*& fileHandle);

        void Write_Cells(
            SIM& svar, MESH const& cells, IPTState const& pnp1, int32_t const& totalNumFaceNodes,
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left,
            vector<int32_t> const& right
        );
    } // namespace BINARY
} // namespace IPT
#endif