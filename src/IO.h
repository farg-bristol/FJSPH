/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H

#include "Var.h"

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void GetInput(int argc, char **argv, SIM& svar, FLUID& fvar, AERO& avar);

void Write_Headers(FILE* f1, FILE* fb, FILE* fg, SIM& svar);

void Write_Timestep(FILE* f1, FILE* fb, FILE* fg, SIM& svar, LIMITS const& limits,
				 SPHState const& pnp1/* , DELTAP const& dp *//* , SPHState const& airP */);

void Append_Restart_Prefix(SIM const& svar);

void Check_Output_Variables(SIM& svar);

void Restart_Simulation(SIM& svar, FLUID const& fvar, AERO const& avar, 
		MESH const& cells, SPHState& pn, SPHState& pnp1, LIMITS& limits);


#endif