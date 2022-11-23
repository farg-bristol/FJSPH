#ifndef CIRCLE_H
#define CIRCLE_H

#include "shapes.h"

void check_circle_input(shape_block& bound, real& globalspacing);

std::vector<StateVecD> create_circle(StateVecD const& centre, real const& radius, 
                real const& globalspacing, int const& hcpl );

#endif