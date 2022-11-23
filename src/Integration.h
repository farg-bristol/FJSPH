/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "Var.h"

///**************** Integration loop **************///
real Integrate(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	MESH& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist, /* DELTAP& dp, */ SPHState& pn, SPHState& pnp1, 
	vector<IPTState>& iptdata);

void First_Step(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	MESH& cells, LIMITS const& limits, OUTL& outlist, /* DELTAP& dp, */ SPHState& pnp1, SPHState& pn, vector<IPTState>& iptdata);

#endif
