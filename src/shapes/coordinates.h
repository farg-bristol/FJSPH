#ifndef COORDINATES_H
#define COORDINATES_H

#include "shapes.h"

class CoordShape : public ShapeBlock
{
  public:
    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault);

    void generate_points(real const& globalspacing);
};

#endif