/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "Var.h"

///**************** Update neighbour list **************
void FindNeighbours(Sim_Tree const& NP1_INDEX, FLUID const& fvar, State const& pnp1, outl& outlist)
{
	const nanoflann::SearchParams params(0,0,false);
	const real search_radius = fvar.sr;
	outlist.clear();
	outlist.reserve(pnp1.size());
	
	#pragma omp parallel
	{	/*Find neighbour list*/
		outl local; /*Local processor copy*/
		local.reserve(pnp1.size());
		#pragma omp for schedule(static) nowait 
		for(size_t ii=0; ii < pnp1.size(); ++ii)
		{
			// std::cout << pnp1[i].list.size();
			std::vector<std::pair<size_t, real>> matches; /* Nearest Neighbour Search*/
			#if SIMDIM == 3
			matches.reserve(250);
			#else
			matches.reserve(20);
			#endif

			NP1_INDEX.index->radiusSearch(&pnp1[ii].xi[0], search_radius, matches, params);
			
			local.emplace_back(std::vector<size_t>(matches.size()));
			for (size_t jj = 0; jj < matches.size(); ++jj)
			{
				local.back()[jj] = matches[jj].first;
			}
			
			// std::cout << "  " << pnp1[i].list.size() << std::endl;
		}

		#pragma omp for schedule(static) ordered
    	for(int ii=0; ii<omp_get_num_threads(); ii++)
    	{
    		#pragma omp ordered
    		outlist.insert(outlist.end(),local.begin(),local.end());
    	}
	}
}

void FindCellNeighbours(Vec_Tree const& CELL_INDEX, vector<StateVecD> const& cells, outl& outlist)
{
	
	#if SIMDIM == 3
	const size_t num_results = 80;
	#else
	const size_t num_results = 80;
	#endif

	outlist = vector<vector<size_t>>(cells.size(),vector<size_t>(num_results));

	#pragma omp parallel
	{	/*Find neighbour list*/
		
		#pragma omp for schedule(static) nowait 
		for(size_t ii=0; ii < cells.size(); ++ii)
		{
			StateVecD testp = cells[ii];
			vector<size_t> ret_indexes(num_results);
			vector<real> out_dists_sqr(num_results);

			nanoflann::KNNResultSet<real> resultSet(num_results);
			resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			
			CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));
			
			outlist[ii] = ret_indexes;			
		}
	}
}


#endif
