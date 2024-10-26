/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Neighbours.h"

/* Update tree index and return a new neighbour list */
OUTL update_neighbours(Sim_Tree const& tree, SPHState const& pnp1)
{
    // Combine to make it impossible to forgot to update the tree index before finding neighbours.
    tree.index->buildIndex();
    return find_neighbours(tree, svar.fluid, pnp1);
}

/* Find the list of neighbours for each particle and return as a 2D vector */
OUTL find_neighbours(Sim_Tree const& NP1_INDEX, SPHState const& pnp1)
{
    const real search_radius = svar.fluid.sr;
    OUTL neighbour_list(pnp1.size());

#pragma omp parallel default(shared)
    { /*Find neighbour list*/
#pragma omp for schedule(static) nowait
        for (size_t ii = 0; ii < pnp1.size(); ++ii)
        {
            neighbour_list[ii] =
                radius_search(NP1_INDEX, pnp1[ii].xi, search_radius); /* Nearest Neighbour Search */
        }
    }

    return neighbour_list;
}

/** Perform a radius search using NanoFLANN */
std::vector<neighbour_index>
radius_search(Sim_Tree const& tree, StateVecD const& test_point, real const& search_radius)
{
    std::vector<neighbour_index> matches; /* Nearest Neighbour Search*/
#if SIMDIM == 3
    matches.reserve(250);
#else
    matches.reserve(47);
#endif

    tree.index->radiusSearch(&test_point[0], search_radius, matches, flann_params);

    return matches;
}

/** Perform a radius search using NanoFLANN */
std::vector<neighbour_index>
radius_search(Vec_Tree const& tree, StateVecD const& test_point, real const& search_radius)
{
    std::vector<neighbour_index> matches; /* Nearest Neighbour Search*/
#if SIMDIM == 3
    matches.reserve(250);
#else
    matches.reserve(47);
#endif

    tree.index->radiusSearch(&test_point[0], search_radius, matches, flann_params);

    return matches;
}