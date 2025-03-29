/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef CYLINDER_H
#define CYLINDER_H

#include "shapes.h"

class CylinderShape : public ShapeBlock
{
  public:
    CylinderShape() { bound_type = cylinder; };

    void check_input(SIM const& svar, real& globalspacing, int& fault) override;

    void generate_points(real const& globalspacing) override;
};

#endif