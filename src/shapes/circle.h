#ifndef CIRCLE_H
#define CIRCLE_H

#include "shapes.h"

class CircleShape : public ShapeBlock
{
  public:
    CircleShape() { bound_type = circleSphere; };

    void check_input(SIM const& svar, int& fault) override;

    void generate_points() override;
};

#endif