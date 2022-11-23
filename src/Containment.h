/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef CONTAINMENT_H
#define CONTAINMENT_H

#include "Var.h"

void Write_Containment(SIM const& svar, vector<size_t> const& ret_indexes, MESH const& cells, StateVecD const& testp);

uint CheckCell(size_t const& cell, MESH const& cells, StateVecD const& testp);


/* <summary> Find the first cell for particles transitioning from the pipe to being in free air.  */
/* It is assumed that the particle does not have a previous cell defined, and so goes */
/* immediately to the KD Tree to find the nearest cells to iterate through. </summary> */
void FirstCell(SIM& svar, Vec_Tree const& CELL_INDEX, MESH const& cells, SPHPart& pi, uint& to_del);

/* <summary> Find the cells for all non-boundary particles. Checks for whether a paricle is free.  */
/* It is assumed that the particle does have a previous cell defined, and checks this cell */
/* before going to the KD Tree to find the nearest cells to iterate through. </summary> */
vector<size_t> FindCell(SIM& svar, AERO const& avar, Vec_Tree const& CELL_TREE, MESH const& cells,
     /* DELTAP const& dp, */ SPHState& pn, SPHState& pnp1);

void Check_Pipe_Outlet(Vec_Tree const& CELL_TREE, SIM& svar, AERO const& avar, MESH const& cells, 
            LIMITS& limits, /* DELTAP& dp, */ SPHState& pn, SPHState& pnp1, size_t& end, size_t& end_ng);

/* IMPLICIT PARTICLE TRACKING FUNCTIONS */


real FindFace(SIM const& svar, MESH const& cells, IPTPart const& pn, IPTPart& pnp1);


#endif