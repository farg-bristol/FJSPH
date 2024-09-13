#ifndef CIRCLE_H
#define CIRCLE_H

#include "shapes.h"

class CircleShape : public ShapeBlock
{
  public:
    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault);

    void generate_points(real const& globalspacing);
};

#endif