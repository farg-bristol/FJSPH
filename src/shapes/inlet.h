/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef INLET_H
#define INLET_H

#include "shapes.h"

class InletShape : public ShapeBlock
{
  public:
    InletShape() { bound_type = inletZone; };

    void check_input(SIM const& svar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

uint update_buffer_region(SIM& svar, LIMITS& limits, SPHState& pnp1, size_t& end);

#endif