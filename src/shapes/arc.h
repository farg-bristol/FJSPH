#ifndef ARC_H
#define ARC_H

#include "shapes.h"

namespace ArcShape
{
    void check_input(shape_block &block, real &globalspacing, int &fault);

    std::vector<StateVecD> generate_points(shape_block const &block, real const &globalspacing);
}

#endif