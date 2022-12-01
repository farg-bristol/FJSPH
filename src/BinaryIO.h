/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include "Var.h"
#include <TECIO.h>

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

/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/
inline void Write_Real_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, 
								int32_t const& start, int32_t const& size,
								vector<real> const& varVec, string const& varName)
{
	int retval;
	#if FOD == 1
		retval =  tecZoneVarWriteDoubleValues(fileHandle, outputZone, varCount, 0, size, &varVec[start]);
	#else
		retval =  tecZoneVarWriteFloatValues(fileHandle, outputZone, varCount, 0, size, &varVec[start]);
	#endif	

	if(retval)
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}


inline void Write_State_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
								int32_t const& start, int32_t const& size,
								vector<StateVecD> const& varVec, string const& varName)
{
	vector<real> vec(size);
	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(int ii = start; ii < start+size; ++ii)
			vec[ii-start] = varVec[ii][dim];

		string name = varName + " " + std::to_string(dim);
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, size, vec, name);
	}
}

inline void Write_Int_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, 
								int32_t const& start, int32_t const& size,
								vector<int32_t> const& varVec, string const& varName)
{
	if(tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &varVec[start]))
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

inline void Write_UInt_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, 
								int32_t const& start, int32_t const& size,
								vector<uint8_t> const& varVec, string const& varName)
{
	if(tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &varVec[start]))
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

inline void Write_Real_Value(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
						real const& value, string const& varName)
{
	int retval;
	vector<real> rvec(size,0.0);
	rvec[0] = value;
	#if FOD == 1
		retval = tecZoneVarWriteDoubleValues(fileHandle, outputZone, varCount, 0, size, &rvec[0]);
	#else
		retval = tecZoneVarWriteFloatValues(fileHandle, outputZone, varCount, 0, size, &rvec[0]);
	#endif	

	if(retval)
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

inline void Write_Int_Value(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
						int32_t const& value, string const& varName)
{
	int retval;
	vector<int32_t> rvec(size,0.0);
	rvec[0] = value;

	retval = tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &rvec[0]);

	if(retval)
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

inline void Write_UInt_Value(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
						uint8_t const& value, string const& varName)
{
	int retval;
	vector<uint8_t> rvec(size,0.0);
	rvec[0] = value;

	retval = tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &rvec[0]);

	if(retval)
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_Binary_Timestep(SIM const& svar, real const& rho0, SPHState const& pnp1,
	bound_block const& limits, char const* group, int32_t const& strandID, void* const& fileHandle);

void Init_Binary_PLT(SIM &svar, string const& filename, string const& zoneName, void* &fileHandle);

/*************************************************************************/
/**************************** BINARY INPUTS ******************************/
/*************************************************************************/
void Restart_Binary(SIM& svar, FLUID const& fvar, SPHState& pn, LIMITS& limits);

void Write_Cell_Centres(MESH const& cells);

#endif