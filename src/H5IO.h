#ifndef H5IO_H
#define H5IO_H

#include "VLM.h"
#include "Var.h"

#include <hdf5.h>

namespace HDF5
{
    hid_t Open_HDF5_File(std::string const& filename, SIM const& svar);

    void Write_Attribute(int64_t const& fout, std::string const& attr_name, size_t const& attr_data);
    void Write_Attribute(int64_t const& fout, std::string const& attr_name, uint const& attr_data);
    void Write_Attribute(int64_t const& fout, std::string const& attr_name, int const& attr_data);
    void Write_Attribute(int64_t const& fout, std::string const& attr_name, double const& attr_data);
    void Write_Attribute(int64_t const& fout, std::string const& attr_name, StateVecD const& attr_data);
    void
    Write_Attribute(int64_t const& fout, std::string const& attr_name, std::string const& attr_data);

    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<size_t>> const& var_data
    );
    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<int>> const& var_data
    );
    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<float>> const& var_data
    );
    void Write_Variable_Array(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<double>> const& var_data
    );

    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<size_t> const& var_data
    );
    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<uint> const& var_data
    );
    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<int> const& var_data
    );
    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<long> const& var_data
    );
    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<float> const& var_data
    );
    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<double> const& var_data
    );
    void Write_Variable(
        int64_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<StateVecD> const& var_data
    );

} // namespace HDF5

void Write_HDF5(SIM& svar, SPHState const& pnp1, LIMITS const& limits);

void Read_HDF5(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

void open_h5part_files(SIM& svar, string const& prefix);

void close_h5part_files(SIM& svar);

void write_h5part_data(SIM& svar, LIMITS const& limits, SPHState const& pnp1);

#endif