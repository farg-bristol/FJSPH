/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEWMARK_BETA_H
#define NEWMARK_BETA_H

#include "Var.h"

int Check_Error(Sim_Tree& SPH_TREE, SIM& svar, FLUID const& fvar, size_t const& start, size_t const& end, 
		real& error1, real& error2, real& logbase, OUTL& outlist, vector<StateVecD> const& xih,
		SPHState& pn, SPHState& pnp1, uint& k, uint& nUnstab);

void Do_NB_Iter(Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, 
	real const& a, real const& b, real const& c, real const& d, real const& B, real const& gam,
	real const& npd, MESH& cells, LIMITS const& limits, OUTL& outlist, /* DELTAP const& dp, */ 
	SPHState& pn, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel);

void Newmark_Beta(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	size_t const& start, size_t& end, real const& a, real const& b, real const& c, real const& d,
	real const& B, real const& gam, real const& npd, MESH& cells, LIMITS const& limits, OUTL& outlist, /* DELTAP const& dp, */
	real& logbase, uint& k, real& error1, real& error2, vector<StateVecD>& xih,
	SPHState& pn, SPHState& pnp1, StateVecD& Force, StateVecD& dropVel);

#endif