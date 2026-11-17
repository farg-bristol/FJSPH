#ifndef H5IO_H
#define H5IO_H

#include "VLM.h"
#include "Var.h"

#include <hdf5.h>

namespace HDF5
{
    hid_t Open_File(std::string const& filename, bool const& reading = false);
    hid_t Open_Zone(hid_t const& file_handle, string const& zone_header);
    hid_t Create_Zone(hid_t const& file_handle, string const& zone_header);
    void Close_Handle(hid_t handle);

    void Write_Attribute(hid_t const& fout, std::string const& attr_name, size_t const& attr_data);
    void Write_Attribute(hid_t const& fout, std::string const& attr_name, uint const& attr_data);
    void Write_Attribute(hid_t const& fout, std::string const& attr_name, int const& attr_data);
    void Write_Attribute(hid_t const& fout, std::string const& attr_name, double const& attr_data);
    void Write_Attribute(hid_t const& fout, std::string const& attr_name, StateVecD const& attr_data);
    void Write_Attribute(hid_t const& fout, std::string const& attr_name, std::string const& attr_data);

    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<size_t>> const& var_data
    );
    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<int>> const& var_data
    );
    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<float>> const& var_data
    );
    void Write_Variable_Array(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<std::vector<double>> const& var_data
    );

    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<size_t> const& var_data
    );
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<uint> const& var_data
    );
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<int> const& var_data
    );
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<long> const& var_data
    );
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<float> const& var_data
    );
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<double> const& var_data
    );
    void Write_Variable(
        hid_t const& fout, std::string const& var_name, hsize_t const* dims,
        std::vector<StateVecD> const& var_data
    );

    void Write_Inlet_Data(hid_t const& fout, bound_block const& limits);

    void Write_Particle_Array(
        hid_t const& fout, std::string var_name, SPHState const& pnp1, size_t const& start,
        size_t const& end, hid_t restart_type
    );

    // Reading
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, std::string& attr_data);
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, size_t& attr_data);
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, uint& attr_data);
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, int& attr_data);
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, double& attr_data);
    void Read_Attribute(hid_t const& fin, std::string const& attr_name, StateVecD& attr_data);

    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<size_t>>& var_data
    );
    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<int>>& var_data
    );
    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<float>>& var_data
    );
    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<double>>& var_data
    );
    void Read_Variable_Array(
        hid_t const& fin, std::string const& var_name, std::vector<std::vector<size_t>>& var_data
    );

    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<size_t>& var_data);
    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<uint>& var_data);
    void Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<int>& var_data);
    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<long>& var_data);
    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<float>& var_data);
    void
    Read_Variable_Scalar(hid_t const& fin, std::string const& var_name, std::vector<double>& var_data);
    void Read_Variable_Scalar(
        hid_t const& fin, std::string const& var_name, std::vector<StateVecD>& var_data
    );

    void
    Read_Particle_Array(hid_t const& fin, std::string var_name, SPHState& var_data, hid_t restart_type);

    hid_t generate_compound_particle();
} // namespace HDF5

void Read_HDF5(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

#endif