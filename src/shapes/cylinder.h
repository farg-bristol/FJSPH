#ifndef CYLINDER_H
#define CYLINDER_H

#include "shapes.h"

namespace CylinderShape
{
    void check_input(shape_block &block, real &globalspacing, int &fault);

    std::vector<StateVecD> generate_points(shape_block const &block, real const &globalspacing);
}

#endif