/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef IO_H
#define IO_H

#include "VLM.h"
#include "Var.h"

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void Set_Values(SIM& svar);

void GetInput(int argc, char** argv, SIM& svar);

void init_writers(SIM& svar);
void close_writers();
void combine_szl_files();

void Write_Timestep(SIM& svar, LIMITS& limits, SPHState& pnp1);

void Remove_Old_Files(SIM const& svar);

void Check_Output_Variables(SIM& svar);

#endif