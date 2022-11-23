#ifndef LINE_H
#define LINE_H

#include "shapes.h"

void check_line_input(shape_block& bound, real& globalspacing);

#if SIMDIM == 2
// Fine line of one particle thick, used in particular for fibres
std::vector<StateVecD> create_fine_line(StateVecD const& start, StateVecD const& end, 
                         real const& globalspacing);

// Create line with n thick particles or a given thickness.
std::vector<StateVecD> create_line(shape_block const& block, real const& globalspacing);

#else 
// Three dimensional plane with n thick particles or a given thickness.
std::vector<StateVecD> create_plane(shape_block const& block, real const& globalspacing);

// Three dimensional plane with i and j indexing for fibre sheets (only 1 particle thick)
std::vector<StateVecD> create_fibre_plane(shape_block& var, real const globalspacing);

#endif        

#endif