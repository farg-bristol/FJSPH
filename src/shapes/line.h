/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef LINE_H
#define LINE_H

#include "shapes.h"

namespace LineShape
{
void check_input(shape_block& block, real& globalspacing, int& fault);

// Create line with n thick particles or a given thickness.
std::vector<StateVecD> generate_points(shape_block const& block, real const& globalspacing);
} // namespace LineShape

#endif