#ifndef H5IO_H
#define H5IO_H

#include "VLM.h"
#include "Var.h"

void Write_HDF5(SIM& svar, SPHState const& pnp1, LIMITS const& limits);

void Read_HDF5(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

void open_h5part_files(
    SIM const& svar, string const& prefix, int64_t& h5part_fluid_file, int64_t& h5part_bound_file
);

void write_h5part_data(
    int64_t& h5part_fluid_file, int64_t& h5part_bound_file, SIM const& svar, SPHState const& pnp1
);

#endif