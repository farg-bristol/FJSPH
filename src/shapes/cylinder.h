/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef CYLINDER_H
#define CYLINDER_H

#include "shapes.h"

class CylinderShape : public ShapeBlock
{
  public:
    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault);

    void generate_points(real const& globalspacing);
};

#endif