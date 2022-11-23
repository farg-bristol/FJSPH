/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "Neighbours.h"


///**************** Update neighbour list **************
void FindNeighbours(Sim_Tree const& NP1_INDEX, FLUID const& fvar, SPHState const& pnp1, OUTL& outlist)
{
	const nanoflann::SearchParams params(0,0,false);
	const real search_radius = fvar.sr;
	outlist.clear();
	outlist.reserve(pnp1.size());
	
	// vector<vector<std::pair<size_t,real>>> plist;
	// plist.reserve(pnp1.size());

	#pragma omp parallel default(shared)
	{	/*Find neighbour list*/
		OUTL local; /*Local processor copy*/
		// vector<vector<std::pair<size_t,real>>> plocal;
		local.reserve(pnp1.size());
		// plocal.reserve(pnp1.size());
		#pragma omp for schedule(static) nowait 
		for(size_t ii=0; ii < pnp1.size(); ++ii)
		{
			// std::cout << pnp1[i].list.size();
			std::vector<std::pair<size_t, real>> matches; /* Nearest Neighbour Search*/
			#if SIMDIM == 3
				matches.reserve(250);
			#else
				matches.reserve(47);
			#endif

			NP1_INDEX.index->radiusSearch(&pnp1[ii].xi[0], search_radius, matches, params);
			
			local.emplace_back(matches);
			// plocal.emplace_back(vector<std::pair<size_t,real>>(matches.size()-1));
			// for (size_t jj = 1; jj < matches.size(); ++jj)
			// {
			// 	local.back()[jj-1] = matches[jj].first;
			// 	// plocal.back()[jj-1] = matches[jj];
			// }
			
			// std::cout << "  " << pnp1[i].list.size() << std::endl;
		}

		#pragma omp for schedule(static) ordered
    	for(int ii=0; ii<omp_get_num_threads(); ii++)
    	{
    		#pragma omp ordered
    		{
    			outlist.insert(outlist.end(),local.begin(),local.end());
    			// plist.insert(plist.end(),plocal.begin(),plocal.end());
			}
    	}
	}


	// std::fstream fn("Neighbours",std::ios::out);
	// for(size_t ii = 0; ii < outlist.size(); ++ii)
	// {
	// 	fn << "Particle: " << ii << endl;

	// 	for (size_t jj = 0; jj < outlist[ii].size(); ++jj)
	// 	{
	// 		fn << "\tNeighbour: " << outlist[ii][jj].first << "  Distance: " << sqrt(outlist[ii][jj].second) << endl;
	// 	}
	// }
	// fn.close();

}

void FindCellNeighbours(Vec_Tree const& CELL_INDEX, vector<StateVecD> const& cells, celll& outlist)
{
	
	#if SIMDIM == 3
	const size_t num_results = 80;
	#else
	const size_t num_results = 80;
	#endif

	outlist = vector<vector<size_t>>(cells.size(),vector<size_t>(num_results));

	#pragma omp parallel default(shared)
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

