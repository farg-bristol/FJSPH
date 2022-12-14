/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"

void Init_Particles(SIM& svar, FLUID& fvar, AERO& avar, SPHState& pn, SPHState& pnp1, LIMITS& limits);

void Init_Surface(SIM const& svar, MESH const& cells, vector<SURF>& surf_marks);

#endif