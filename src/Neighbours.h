/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "Var.h"

///**************** Update neighbour list **************
void FindNeighbours(const Sim_Tree& NP1_INDEX, const FLUID& fvar, const State& pnp1, outl& outlist)
{
	const nanoflann::SearchParams params;
	const ldouble search_radius = fvar.sr;
	outlist.clear();
	
	#pragma omp parallel
	{	/*Find neighbour list*/
		outl local; /*Local processor copy*/
		#pragma omp for schedule(static) nowait 
		for(uint i=0; i < pnp1.size(); ++i)
		{
			// std::cout << pnp1[i].list.size();
			std::vector<std::pair<size_t, ldouble>> matches; /* Nearest Neighbour Search*/
			NP1_INDEX.index->radiusSearch(&pnp1[i].xi[0], search_radius, matches, params);
			
			std::vector<uint> temp;
			for (auto &j:matches)
			{
				temp.emplace_back(static_cast<uint>(j.first));
			}
			local.emplace_back(temp);
			// std::cout << "  " << pnp1[i].list.size() << std::endl;
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<NTHREADS; i++)
    	{
    		#pragma omp ordered
    		outlist.insert(outlist.end(),local.begin(),local.end());
    	}
	}
}

#endif
