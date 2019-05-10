/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "Var.h"

///**************** Update neighbour list **************
void FindNeighbours(Sim_Tree &NP1_INDEX, FLUID &fvar, State &pnp1,outl &outlist)
{
	nanoflann::SearchParams params;
	outlist.erase(outlist.begin(),outlist.end());
	ldouble search_radius = fvar.sr;
	/*Find neighbour list*/
	for(unsigned int i=0; i<pnp1.size(); ++i)
	{
		std::vector<std::pair<size_t, ldouble>> matches; /* Nearest Neighbour Search*/

		NP1_INDEX.index->radiusSearch(&pnp1[i].xi[0], search_radius, matches, params);

		std::vector<unsigned int> temp;
		for (auto &j:matches)
		{
			temp.emplace_back(j.first);
		}
		outlist.emplace_back(temp);
	}
}

#endif
