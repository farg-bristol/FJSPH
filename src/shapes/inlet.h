#ifndef INLET_H
#define INLET_H

#include "shapes.h"

void check_inlet_input(shape_block& bound, real& globalspacing);

std::vector<StateVecD> create_inlet_zone(shape_block& block, real const& globalspacing);

uint update_buffer_region(SIM& svar, LIMITS& limits, SPHState& pnp1, size_t& end, size_t& end_ng);
#endif