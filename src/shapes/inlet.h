/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef INLET_H
#define INLET_H

#include "shapes.h"

class InletShape : public ShapeBlock
{
  public:
    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault);

    void generate_points(real const& globalspacing);
};

uint update_buffer_region(SIM& svar, LIMITS& limits, SPHState& pnp1, size_t& end, size_t& end_ng);

#endif