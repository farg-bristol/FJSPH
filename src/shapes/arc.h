#ifndef ARC_H
#define ARC_H

#include "shapes.h"

class ArcShape : public ShapeBlock
{
  public:
    ArcShape() { bound_type = arcSection; };

    void check_input(SIM const& svar, int& fault) override;

    void generate_points() override;
};

#endif