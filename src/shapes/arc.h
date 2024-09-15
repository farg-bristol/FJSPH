#ifndef ARC_H
#define ARC_H

#include "shapes.h"

class ArcShape : public ShapeBlock
{
  public:
    ArcShape() { bound_type = arcSection; };

    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

#endif