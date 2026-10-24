/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef LINE_H
#define LINE_H

#include "shapes.h"

class LineShape : public ShapeBlock
{
  public:
    LineShape() { bound_type = linePlane; };

    void check_input(SIM const& svar, int& fault) override;

    void generate_points() override;
};

#endif