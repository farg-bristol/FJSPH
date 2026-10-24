/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef ASCIIIO_H
#define ASCIIIO_H

#include "Var.h"

namespace ASCII
{
    void add_file_aux_data(FILE* fp, string const& name, string const& var);
    void add_file_aux_data(FILE* fp, string const& name, size_t const& var);
    void add_file_aux_data(FILE* fp, string const& name, int const& var);
    void add_file_aux_data(FILE* fp, string const& name, uint const& var);
    void add_file_aux_data(FILE* fp, string const& name, double const& var);
    void add_file_aux_data(FILE* fp, string const& name, float const& var);
    void add_file_aux_data(FILE* fp, string const& name, StateVecD const& var);

    void add_zone_aux_data(FILE* fp, string const& name, string const& var);
    void add_zone_aux_data(FILE* fp, string const& name, size_t const& var);
    void add_zone_aux_data(FILE* fp, string const& name, int const& var);
    void add_zone_aux_data(FILE* fp, string const& name, uint const& var);
    void add_zone_aux_data(FILE* fp, string const& name, double const& var);
    void add_zone_aux_data(FILE* fp, string const& name, float const& var);
    void add_zone_aux_data(FILE* fp, string const& name, StateVecD const& var);

    void Write_Vector(FILE* fp, vector<real> const& vec);
    void Write_Vector(FILE* fp, vector<uint> const& vec);
    void Write_Vector(FILE* fp, vector<long> const& vec);
    void Write_Vector(FILE* fp, vector<size_t> const& vec);
    void Write_Vector(FILE* fp, vector<int> const& vec);
    void Write_Vector(FILE* fp, vector<StateVecD> const& vec);

    void Write_File_header(FILE* fp, string const& var_names, string const& title);
    void Write_Zone_Header(
        FILE* fp, std::string zone_header, int64_t size, int32_t strand_id, real solution_time
    );
} // namespace ASCII

/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/

void Write_Cell_Data(MESH const& cdata);

void Write_Face_Data(MESH const& cells);

#endif