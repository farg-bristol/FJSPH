/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef SQUARE_H
#define SQUARE_H

#include "shapes.h"

class SquareShape : public ShapeBlock
{
  public:
    SquareShape() { bound_type = squareCube; };

    virtual void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault);

    virtual void generate_points(real const& globalspacing);

    virtual ~SquareShape(){};
};

#endif