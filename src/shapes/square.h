#ifndef SQUARE_H
#define SQUARE_H

#include "shapes.h"

void check_square_input(shape_block& bound, real& globalspacing);

std::vector<StateVecD> create_square(StateVecD const& start, StateVecD const& end, 
                        real const& globalspacing, int const& hcpl);

#endif