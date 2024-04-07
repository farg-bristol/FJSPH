#ifndef CIRCLE_H
#define CIRCLE_H

#include "shapes.h"

namespace CircleShape
{
    void check_input(shape_block &block, real &globalspacing, int &fault);

    std::vector<StateVecD> generate_points(shape_block const &block, real const &globalspacing);
}

#endif