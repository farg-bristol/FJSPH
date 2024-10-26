/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef IO_H
#define IO_H

#include "VLM.h"
#include "Var.h"

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void Set_Values(SIM& svar, FLUID& svar.fluid, AERO& svar.air, VLM& vortex);

void GetInput(int argc, char** argv, SIM& svar, FLUID& svar.fluid, AERO& svar.air, VLM& vortex);

void Write_Tec_Headers(
    FILE* ff, FILE* fb, FILE* fg, SIM& svar, AERO const& svar.air, std::string const& prefix
);

void Write_h5part_Headers(SIM& svar, AERO const& svar.air, std::string const& prefix);

void Write_Timestep(
    FILE* ff, FILE* fb, FILE* fg, SIM& svar, AERO const& svar.air, LIMITS const& limits,
    SPHState const& pnp1
);

void Remove_Old_Files(SIM const& svar);

void Check_Output_Variables(SIM& svar);

#endif