#ifndef CIRCLE_H
#define CIRCLE_H

#include "shapes.h"

class CircleShape : public ShapeBlock
{
  public:
    CircleShape() { bound_type = circleSphere; };

    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

#endif