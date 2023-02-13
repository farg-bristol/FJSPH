#ifndef ARC_H
#define ARC_H

#include "shapes.h"

void check_arc_input(shape_block& bound, real& globalspacing, int& fault);

#if SIMDIM == 2
void get_arc_end(StateVecD const& start, StateVecD const& centre, real const& arclength, 
                StateVecD& end, real& radius, real& theta0);
                
void get_arclength_centrepoint(StateVecD const& start, StateVecD const& end, StateVecD const& centre,
                        real& radius, real& theta0, real& theta1);

void get_arclength_midpoint(StateVecD const& start, StateVecD const& end, StateVecD const& midpoint, 
                StateVecD& centre, real& radius, real& theta0, real& theta1);
#else
void get_arc_end(StateVecD const& start, StateVecD const& centre, StateVecD const& normal, 
                 real const& arclength, StateVecD& end, real& radius);

void get_arclength_centrepoint(StateVecD const& start, StateVecD const& end, StateVecD const& centre, StateVecD& right,
                        real& radius, real& theta0, real& theta1);

void get_arclength_midpoint(StateVecD const& start, StateVecD const& end, StateVecD const& midpoint, 
                StateVecD& centre, StateVecD& right, real& radius, real& theta0, real& theta1);
#endif

std::vector<StateVecD> create_arc_segment(shape_block const& block, real const& globalspacing);

std::vector<StateVecD> create_fibre_arch(shape_block& var, real const globalspacing);
#endif