#ifndef H5IO_H
#define H5IO_H

#include "Var.h"
#include "VLM.h"

void Write_HDF5(SIM& svar, FLUID const& fvar, AERO const& avar, 
        SPHState const& pnp1, LIMITS const& limits);

void Read_HDF5(SIM& svar, FLUID& fvar, AERO& avar, VLM& vortex, SPHState& pn, SPHState& pnp1, LIMITS& limits);

void open_h5part_files(SIM const& svar, FLUID const& fvar, AERO const& avar, 
                     int64_t& ffile, int64_t& bfile);

void write_h5part_data(int64_t& ffile, int64_t& bfile,
    SIM const& svar, FLUID const& fvar, SPHState const& pnp1);

#endif