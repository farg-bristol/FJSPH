#ifndef INLET_H
#define INLET_H

#include "shapes.h"

namespace InletShape
{
    void check_input(shape_block &block, real &globalspacing, int &fault);

    std::vector<StateVecD> generate_points(shape_block &block, real const &globalspacing);
}

uint update_buffer_region(SIM &svar, LIMITS &limits, SPHState &pnp1, size_t &end, size_t &end_ng);


#endif