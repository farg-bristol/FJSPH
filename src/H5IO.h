#ifndef H5IO_H
#define H5IO_H

#include "VLM.h"
#include "Var.h"

#include <hdf5.h>

namespace HDF5
{
    hid_t Open_HDF5_File(std::string const& filename, SIM const& svar);

    template <typename T>
    void Write_Attribute(
        int64_t const& fout, std::string const& attr_name, T const& attr_data, hid_t const& var_type,
        hid_t const& write_type
    );

    void
    Write_Uint_Attribute(int64_t const& fout, std::string const& attr_name, size_t const& attr_data);
    void Write_Uint_Attribute(int64_t const& fout, std::string const& attr_name, uint const& attr_data);
    void Write_Int_Attribute(int64_t const& fout, std::string const& attr_name, int const& attr_data);
    void
    Write_Real_Attribute(int64_t const& fout, std::string const& attr_name, double const& attr_data);

    void Write_Vector_Attribute(
        int64_t const& fout, std::string const& attr_name, StateVecD const& attr_data
    );
    void Write_String_Attribute(
        int64_t const& fout, std::string const& attr_name, std::string const& attr_data
    );

    template <typename T>
    void Write_Var_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<T>> const& var_data, hid_t const& var_type, hid_t const& write_type
    );

    template <typename T>
    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<T>> const& var_data
    );

    template <typename T>
    void Write_Variable_Scalar(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<T> const& var_data
    );

    void
    Write_Vector_Data(int64_t const& fout, std::string const& name, std::vector<StateVecD> const& data);

} // namespace HDF5

void Write_HDF5(SIM& svar, SPHState const& pnp1, LIMITS const& limits);

void Read_HDF5(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

void open_h5part_files(SIM& svar, string const& prefix);

void close_h5part_files(SIM& svar);

void write_h5part_data(SIM& svar, LIMITS const& limits, SPHState const& pnp1);

#endif