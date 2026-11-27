/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#pragma once

#include "Var.h"

#include <hdf5.h>

class FileReadWriter
{
  public:
    void close_file();
    bool file_is_open();
    bool zone_is_open();

  protected:
    FileReadWriter(){};
    FileReadWriter(file_extension const& _file_ext, SIM const& svar);

    ~FileReadWriter();

    void open_file(std::string const& filename, std::string const& zone_name = "");
    virtual void rw_file_attributes(SIM& svar);
    virtual void rw_zone_attributes(bound_block& limits);

    inline string get_dim(uint const& dim)
    {
        switch (dim)
        {
        case 0:
            return "x";
        case 1:
            return "y";
        case 2:
            return "z";
        default:
            return "x";
        }
    }

    template <typename T> void write_file_attribute(std::string const& attr_name, T const& attr_data);
    template <typename T> void write_zone_attribute(std::string const& attr_name, T const& attr_data);
    template <typename T>
    void
    write_scalar_var(std::string const& var_name, std::vector<T> const& var_data, int32_t& var_count);

    template <typename T> void read_file_attribute(std::string const& attr_name, T& attr_data);
    template <typename T> void read_zone_attribute(std::string const& attr_name, T& attr_data);
    template <typename T>
    void read_scalar_var(std::string const& var_name, std::vector<T>& var_data, int32_t& var_count);

    // Functions that do the appropriate thing based on reading or writing.
    template <typename T> void rw_file_attribute(std::string const& attr_name, T& attr_data);
    template <typename T> void rw_zone_attribute(std::string const& attr_name, T& attr_data);
    template <typename T>
    void rw_scalar_var(std::string const& var_name, std::vector<T>& var_data, int32_t& var_count);

    void create_zone(
        std::string const& zone_header, real const& solution_time = 0, int64_t const size = 0,
        int32_t const strandID = 0
    );
    void open_zone(std::string const& zone_header);
    void close_zone();

    file_extension file_ext;
    std::string ext_str;
    real scale = 1.0;
    ReadOrWrite read_or_write;

    hid_t h5_file_handle;
    hid_t h5_zone_handle;

    void* tec_file_handle;
    int32_t tec_output_zone;
    bool single_file;
    string var_names;
    vector<int32_t> var_types;

    FILE* tec_dat_handle;
};

class RestartFile : public FileReadWriter
{
  public:
    RestartFile(SIM const& svar);

    void write_restart_file(SIM& svar, SPHState& pnp1, LIMITS& limits);
    void read_restart_file(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

  private:
    void write_restart_zone(
        std::string const& zone_header, SPHState const& pnp1, LIMITS& limits, size_t const& block_index
    );
    void read_restart_zone(
        std::string const& zone_header, SPHState& pn, LIMITS& limits, size_t const& block_index
    );

    hid_t restart_type;
};

class DataFile : public FileReadWriter
{
  public:
    DataFile(){};
    DataFile(file_extension const& _file_ext, SIM const& svar, std::string const& _data_type);

    void rw_file_attributes(SIM& svar) override;
    void rw_zone_attributes(bound_block& limits) override;

  protected:
    OutputMap output_var;
    std::string output_prefix;
    std::string data_type;  // Fluid or Boundary particles
    std::string data_lower; // lower case of the same
    std::string file_name;

    string get_filename(real time = -1);

    void write_zone(SPHState const& pnp1, real const& rho_rest, size_t const& start, size_t const& end);
};

class H5PartFile : public DataFile
{
  public:
    H5PartFile(){};
    H5PartFile(SIM const& svar, std::string const& _data_type);

    void write_timestep(
        SPHState const& pnp1, LIMITS const& limits, size_t start, size_t end, uint frame,
        real const& time
    );
};

class SZLFile : public DataFile
{

  public:
    SZLFile(){};
    SZLFile(
        SIM& svar, std::string const& _data_type, size_t const& _start_block, size_t const& _end_block
    );

    void write_timestep(SIM& svar, SPHState const& pnp1, LIMITS& limits, real const& time);

    void combine_file();
    void flush_file();

  private:
    void write_timestep_data(SPHState const& pnp1, LIMITS& limits, real const& time);
    size_t start_block, end_block;
};

class ASCIIFile : public DataFile
{
  public:
    ASCIIFile(){};
    ASCIIFile(
        SIM& svar, std::string const& _data_type, size_t const& _start_block, size_t const& _end_block
    );

    void write_timestep(SIM& svar, SPHState const& pnp1, LIMITS& limits, real const& time);

  private:
    void write_timestep_data(SPHState const& pnp1, LIMITS& limits, real const& time);
    size_t start_block, end_block;
};