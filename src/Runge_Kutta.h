/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "Var.h"

/* <summary> Peform the first stage of the Runge-Kutta integration
	to get the first guess of time n+1 (regarded here as time n+1/4)
	to perform neighbour search and dissipation terms before freezing </summary */
real Get_First_RK(SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t const& end, real const& B, real const& gam, real const& npd, MESH& cells,
	LIMITS const& limits, OUTL const& outlist, real& logbase, SPHState& pn, SPHState& st_2, real& error1);


real Runge_Kutta4(Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, real const& B, real const& gam, real const& npd, MESH& cells, 
	LIMITS const& limits, OUTL const& outlist, real& logbase,
	SPHState& pn, SPHState& st_2, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel);

#endif