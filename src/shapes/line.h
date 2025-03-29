/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef LINE_H
#define LINE_H

#include "shapes.h"

class LineShape : public ShapeBlock
{
  public:
    LineShape() { bound_type = linePlane; };

    void check_input(SIM const& svar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

#endif