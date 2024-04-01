/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include "Var.h"

inline int Combine_SZPLT(string& file)
{
	string cmd = "exec szcombine \"";
	cmd.append(file);
	cmd.append("\"");

	#ifdef DEBUG
		fprintf(dbout,"Attempting to combine szplt.\n");
		fprintf(dbout,"Command: %s\n", cmd.c_str());
	#endif

	printf("Combining szplt: %s\n",file.c_str());
	if(system(cmd.c_str()))
	{
    	printf("Could not combine szplt file.\n");
    	printf("Command: %s\n",cmd.c_str());
    	return -1;
	}
	return 0;
}

void Write_Binary_Timestep(SIM const& svar, real const& rho0, SPHState const& pnp1,
	bound_block const& limits, char const* group, int32_t const& strandID, void* const& fileHandle);

void Init_Binary_PLT(SIM &svar, FLUID const& fvar, AERO const& avar, string const& prefix,
			 string const& filename, string const& zoneName, void* &fileHandle);

void close_file(void* handle);

void flush_file(void* handle);

// IPT Binary functions (Allows more explicit library includes, so that tecio is only in BinaryIO.cpp.)
namespace IPT
{
	namespace BINARY
	{
        void Init_IPT_Files(SIM& svar);

		void Write_Point(SIM const& svar, IPTPart const& pnp1);

        void Write_State(SIM const& svar, IPTState const& pnp1, string const& zoneName, void* &fileHandle);

		void Write_Cells(SIM& svar, MESH const& cells, IPTState const& pnp1, int32_t const& totalNumFaceNodes, 
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left, vector<int32_t> const& right);
	}
}
#endif