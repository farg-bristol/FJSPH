/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef FOAMIO_H
#define FOAMIO_H

#include "Var.h"

/* CURRENTLY 3D ONLY */
#if SIMDIM == 3
namespace FOAM
{
    void Read_FOAM(SIM& svar, MESH& cells);
}
#endif
#endif