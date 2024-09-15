#ifndef COORDINATES_H
#define COORDINATES_H

#include "shapes.h"

class CoordShape : public ShapeBlock
{
  public:
    CoordShape() { bound_type = coordDef; };

    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

#endif