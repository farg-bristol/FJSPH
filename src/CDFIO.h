/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef CDFIO_H
#define CDFIO_H

#include "Var.h"

void Average_Point_to_Cell(
    vector<StateVecD> const& pData, vector<StateVecD>& cData, vector<vector<size_t>> const& cFaces,
    vector<vector<size_t>> const& faces
);

void Average_Point_to_Cell(
    vector<real> const& pData, vector<real>& cData, vector<vector<size_t>> const& cFaces,
    vector<vector<size_t>> const& faces
);

/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
namespace TAU
{
    void Read_BMAP(SIM& svar);

    /*****************************************************************************/
    /***************** READING NETCDF SOLUTION DATA FUNCTIONS ********************/
    /*****************************************************************************/
    void Read_SOLUTION(
        SIM const& svar, AERO const& svar.air, uint const ignored, MESH& cells,
        vector<uint> const& usedVerts
    );

/*****************************************************************************/
/*************** READING NETCDF EDGE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
#if SIMDIM == 2
    void Read_tau_mesh_EDGE(SIM& svar, MESH& cells, AERO const& svar.air, vector<uint>& uVerts);
#endif

/*****************************************************************************/
/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
#if SIMDIM == 3
    void Read_tau_mesh_FACE(SIM& svar, MESH& cells, AERO const& svar.air);
#endif
} // namespace TAU
#endif