/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef SPEEDTEST_H
#define SPEEDTEST_H
/* Library to perform drag assessment for a range of particle resolutions */
#include "Var.h"

void Speed_Test_Sweep(SIM& svar, FLUID& fvar, AERO& avar);

#endif