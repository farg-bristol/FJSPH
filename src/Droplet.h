/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef DROPLET_H
#define DROPLET_H

/* Library to perform drag assessment for a range of particle resolutions */
#include "Var.h"

void Droplet_Drag_Sweep(SIM& svar, FLUID& fvar, AERO& avar);

#endif