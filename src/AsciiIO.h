/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef ASCIIIO_H
#define ASCIIIO_H

#include "Var.h"

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void ASCII_Restart(SIM& svar, FLUID const& fvar, SPHState& pn);

/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/
void Write_ASCII_header(FILE* fp, SIM const& svar, char const* title);

void Write_ASCII_Timestep(
    SIM& svar, real const& rho_rest, SPHState const& pnp1, size_t const& start, size_t const& end,
    char const* name, uint const& strandID, FILE* fp
);

void Write_Cell_Data(MESH const& cdata);

void Write_Face_Data(MESH const& cells);

#endif