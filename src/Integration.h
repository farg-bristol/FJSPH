/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "Var.h"
#include "VLM.h"
///**************** Integration loop **************///
real Integrate(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	VLM const& vortex, MESH const& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist,
	SPHState& pn, SPHState& pnp1, vector<IPTState>& iptdata);

void First_Step(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	VLM const& vortex, MESH const& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1, SPHState& pn, vector<IPTState>& iptdata);

#endif
