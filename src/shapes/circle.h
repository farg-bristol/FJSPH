#ifndef CIRCLE_H
#define CIRCLE_H

#include "shapes.h"

void check_circle_input(shape_block& bound, real& globalspacing, int& fault);

std::vector<StateVecD> create_circle(shape_block& block, real const& globalspacing);

#endif