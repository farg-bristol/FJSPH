/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "Var.h"

///**************** Update neighbour list **************
void FindNeighbours(Sim_Tree const& NP1_INDEX, FLUID const& fvar, SPHState const& pnp1, OUTL& outlist);

void FindCellNeighbours(Vec_Tree const& CELL_INDEX, vector<StateVecD> const& cells, celll& outlist);

#endif
