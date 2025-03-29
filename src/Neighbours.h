/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "Var.h"

OUTL update_neighbours(FLUID const& fvar, Sim_Tree const& tree, SPHState const& pnp1);

OUTL find_neighbours(FLUID const& fvar, Sim_Tree const& NP1_INDEX, SPHState const& pnp1);

std::vector<neighbour_index>
radius_search(Sim_Tree const& tree, StateVecD const& test_point, real const& search_radius);

std::vector<neighbour_index>
radius_search(Vec_Tree const& tree, StateVecD const& test_point, real const& search_radius);

#endif
