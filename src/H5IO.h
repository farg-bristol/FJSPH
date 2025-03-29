#ifndef H5IO_H
#define H5IO_H

#include "VLM.h"
#include "Var.h"

void Write_HDF5(SIM& svar, SPHState const& pnp1, LIMITS const& limits);

void Read_HDF5(SIM& svar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

void open_h5part_files(SIM& svar, string const& prefix);

void close_h5part_files(SIM& svar);

void write_h5part_data(SIM& svar, SPHState const& pnp1);

#endif