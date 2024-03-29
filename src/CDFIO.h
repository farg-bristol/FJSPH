#ifndef CDFIO_H
#define CDFIO_H

#include "Var.h"	


void Average_Point_to_Cell(vector<StateVecD> const& pData, vector<StateVecD> &cData,
						   vector<vector<size_t>> const& cFaces, vector<vector<size_t>> const& faces);
						   
void Average_Point_to_Cell(vector<real> const& pData, vector<real> &cData,
						vector<vector<size_t>> const& cFaces, vector<vector<size_t>> const& faces);
						
/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
namespace TAU
{
	void Read_BMAP(SIM& svar);

	/*****************************************************************************/
	/***************** READING NETCDF SOLUTION DATA FUNCTIONS ********************/
	/*****************************************************************************/
	void Read_SOLUTION(SIM const& svar, FLUID const& fvar, AERO const& avar,
					uint const ignored, MESH& cells, vector<uint> const& usedVerts);

	/*****************************************************************************/
	/*************** READING NETCDF EDGE BASED DATA FUNCTIONS ********************/
	/*****************************************************************************/
	#if SIMDIM == 2
	void Read_TAUMESH_EDGE(SIM &svar, MESH &cells, FLUID const &fvar, 
							AERO const &avar, vector<uint>& uVerts);
	#endif

	/*****************************************************************************/
	/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
	/*****************************************************************************/
	#if SIMDIM == 3
	void Read_TAUMESH_FACE(SIM &svar, MESH &cells, FLUID const &fvar, AERO const &avar);
	#endif
}
#endif