#ifndef CYLINDER_H
#define CYLINDER_H

#include "shapes.h"

void check_cylinder_input(shape_block& bound, real& globalspacing);

std::vector<StateVecD> create_cylinder(shape_block const& block, real const& globalspacing);

#endif