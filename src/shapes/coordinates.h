#ifndef COORDINATES_H
#define COORDINATES_H

#include "shapes.h"

class CoordShape : public ShapeBlock
{
  public:
    CoordShape() { bound_type = coordDef; };

    void check_input(SIM const& svar, int& fault) override;

    void generate_points() override;
};

#endif