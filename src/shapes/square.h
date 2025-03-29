/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef SQUARE_H
#define SQUARE_H

#include "shapes.h"

class SquareShape : public ShapeBlock
{
  public:
    SquareShape() { bound_type = squareCube; };

    void check_input(SIM const& svar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

#endif