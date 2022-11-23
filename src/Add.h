/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef ADD_H
#define ADD_H

#include "Var.h"

namespace PoissonSample
{
	/**
		Return a vector of generated points
		sampleLimit - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
	**/
	SPHState generatePoissonPoints(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells,
	 uint const& host, SPHState const& pnp1, OUTL const& outlist/* , StateVecD const& norm, StateVecD const& avgV */);
}

void PoissonGhost(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, 
			Sim_Tree& NP1_INDEX, OUTL& outlist, SPHState& pn, SPHState& pnp1);


void LatticeGhost(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, 
			Sim_Tree& SPH_TREE, OUTL& outlist, SPHState& pn, SPHState& pnp1);

/* Function to check if a ghost particle needs to be removed from the simulation, because */
/* it's left the support of the fluid */
void Check_If_Ghost_Needs_Removing(SIM& svar, FLUID const& fvar, 
				Sim_Tree& NP1_INDEX, SPHState& pn, SPHState& pnp1);

#endif
