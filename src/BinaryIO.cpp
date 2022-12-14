/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "BinaryIO.h"
#include <filesystem> //Requires C++17 standard, needed for restart functionality
#include <chrono>
#include <thread>

using std::string;
using std::filesystem::directory_iterator;

inline std::string space2underscore(std::string text)
{
    std::replace(text.begin(), text.end(), ' ', '_');
    return text;
}

/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/
inline void add_file_aux_data(void* const& file, std::string const& name, 
				std::string const& var)
{
	if(!var.empty())
	{
		if(tecDataSetAddAuxData(file,  space2underscore(name).c_str(), space2underscore(var).c_str()))
		{
			std::cout << "Failed to write file auxiliary data: " << name << std::endl;
			std::cout << "Stopping." << std::endl; exit(-1);
		}
	}
}

inline void add_file_aux_data(void* const& file, std::string const& name, 
				size_t const& var)
{
	std::ostringstream stream;
	stream << var;
	if(tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_file_aux_data(void* const& file, std::string const& name, 
				int const& var)
{
	std::ostringstream stream;
	stream << var;
	if(tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_file_aux_data(void* const& file, std::string const& name, 
				uint const& var)
{
	std::ostringstream stream;
	stream << var;
	if(tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_file_aux_data(void* const& file, std::string const& name, 
				double const& var)
{
	std::ostringstream stream;
	stream << std::setprecision(9) << var;
	if(tecDataSetAddAuxData(file,  space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_file_aux_data(void* const& file, std::string const& name, 
				StateVecD const& var)
{
	string name_ = name + " x";
	std::ostringstream stream;
	stream << std::setprecision(9) << var[0];
	if(tecDataSetAddAuxData(file,  space2underscore(name_).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}

	name_ = name + " y";
	stream.str(""); stream.clear();
	stream << std::setprecision(9) << var[1];
	if(tecDataSetAddAuxData(file,  space2underscore(name_).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
	
	#if SIMDIM == 3
	name_ = name + " z";
	stream.str(""); stream.clear();
	stream << std::setprecision(9) << var[2];
	if(tecDataSetAddAuxData(file,  space2underscore(name_).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write file auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
	#endif
}

inline void add_zone_aux_data(void* const& file, int32_t const& zone,
	 std::string const& name, std::string const& var)
{
	if(!var.empty())
	{
		if(tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), space2underscore(var).c_str()))
		{
			std::cout << "Failed to write zone " << zone << 
					" auxiliary data: " << name << std::endl;
			std::cout << "Stopping." << std::endl; exit(-1);
		}
	}
}

inline void add_zone_aux_data(void* const& file, int32_t const& zone,
	 std::string const& name, size_t const& var)
{
	std::ostringstream stream;
	stream << var;
	if(tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write zone " << zone << 
				" auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_zone_aux_data(void* const& file, int32_t const& zone,
	 std::string const& name, int const& var)
{
	std::ostringstream stream;
	stream << var;
	if(tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write zone " << zone << 
				" auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_zone_aux_data(void* const& file, int32_t const& zone,
	 std::string const& name, double const& var)
{
	std::ostringstream stream;
	stream << std::setprecision(9) << var;
	if(tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write zone " << zone << 
				" auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
}

inline void add_zone_aux_data(void* const& file, int32_t const& zone,
	 std::string const& name, StateVecD const& var)
{
	string name_ = name + " x";
	std::ostringstream stream;
	stream << std::setprecision(9) << var[0];
	if(tecZoneAddAuxData(file, zone, space2underscore(name_).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write zone " << zone << 
				" auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}

	name_ = name + " y";
	stream.str(""); stream.clear();
	stream << var[1];
	if(tecZoneAddAuxData(file, zone, space2underscore(name_).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write zone " << zone << 
				" auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
	
	#if SIMDIM == 3
	name_ = name + " z";
	stream.str(""); stream.clear();
	stream << std::setprecision(9) << var[2];
	if(tecZoneAddAuxData(file, zone, space2underscore(name_).c_str(), stream.str().c_str()))
	{
		std::cout << "Failed to write zone " << zone << 
				" auxiliary data: " << name << std::endl;
		std::cout << "Stopping." << std::endl; exit(-1);
	}
	#endif
}


inline void Write_Zone(SPHState const& pnp1, real const& rho0, real const& scale, std::vector<uint> const& outvar, 
					size_t const& start, size_t const& end, void* const& fileHandle, int32_t const& outputZone)
{
	int32_t varCount = 1;
	int32_t imax = end - start;
	
	vector<real> vec(imax);
	// Restart essential variables. These should not be zero, but could possibly in future
	if(outvar[0])
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].xi(dim)/scale;

			string name = "Position coordinate " + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}

	if(outvar[1])
	{	
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].acc(dim);

			string name = "Acceleration " + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}

	if(outvar[2])
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].v(dim);

			string name = "Velocity " + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}

	if(outvar[3])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].p;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Pressure");
	}

	if(outvar[4])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].Rrho;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Density gradient");
	}
	
	if(outvar[5])
	{
		vector<int> uvec(imax);
		for (size_t ii = start; ii < end; ++ii)
			uvec[ii - start] = pnp1[ii].partID;
		Write_Int_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Particle ID");
	}

	if(outvar[6])
	{
		vector<int> uvec(imax);
		for (size_t ii = start; ii < end; ++ii)
			uvec[ii - start] = pnp1[ii].cellID;
		Write_Int_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Cell ID" );
	}

	if(outvar[7])
	{
		vector<uint8_t> uvec(imax);
		for (size_t ii = start; ii < end; ++ii)
			uvec[ii - start] = pnp1[ii].b;
		Write_UInt_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Particle status condition" );
	}

	// if(outvar[5])
	// 	Write_Real_Vector(sphFile, outputZone, varCount, imax, start, parts.m, "Mass" );
		
	// if(outvar[6])
	// 	Write_Int_Vector(sphFile, outputZone, varCount, imax, start, parts.marker, "Marker" );

	// Non essential variables, but to provide further information
	if(outvar[8])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].rho;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Density");
	}

	if(outvar[9])
	{
		#pragma omp parallel for
		for (size_t ii = start; ii < end; ++ii)
			vec[ii - start] = 100.0*(pnp1[ii].rho/rho0-1.0);
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Density Variation" );
	}
		
	if(outvar[10])
	{
		for (size_t ii = start; ii < end; ++ii)
			vec[ii - start] = pnp1[ii].v.norm();
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Velocity magnitude" );
	}

	if(outvar[11])
	{
		vector<uint8_t> uvec(imax);
		for (size_t ii = start; ii < end; ++ii)
			uvec[ii - start] = pnp1[ii].surf;
		Write_UInt_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Surface flag" );
	}
	
	if(outvar[12])
	{
		vector<uint8_t> uvec(imax);
		for (size_t ii = start; ii < end; ++ii)
			uvec[ii - start] = pnp1[ii].surfzone;
		Write_UInt_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Surface zone flag" );
	}
		
	if(outvar[13])
	{
		for (size_t ii = start; ii < end; ++ii)
			vec[ii - start] = pnp1[ii].Af.norm();
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Aerodynamic force magnitude" );
	}
		
	if(outvar[14])
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].Af(dim);

			string name = "Aerodynamic force " + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}

	if(outvar[15])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].curve;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Curvature");
	}
	
	if(outvar[16])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].woccl;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Occlusion factor");
	}
	
	if(outvar[17])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].cellP;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Cell pressure");
	}
	
	if(outvar[18])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].cellP;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Cell density");
	}
	
	if(outvar[19])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].cellV.norm();
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Cell velocity magnitude");
	}
		
	if(outvar[20])
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].cellV(dim);

			string name = "Cell velocity " + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}
	

	if(outvar[21]) // Delta SPH density gradient
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].gradRho(dim);

			string name = "dSPH density gradient" + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}

	if(outvar[22]) // Lambda eigenvalue of L matrix
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].lam;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Lambda eigenvalue");
	}

	if(outvar[23]) // Lambda eigenvalue of L matrix without boundary
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].lam_nb;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Lambda eigenvalue without boundary");
	}

	if(outvar[24]) // Surface colour function
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].colour;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Colour function");
	}

	if(outvar[25]) // Gradient of the colour function for CSF
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].colourG;
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Colour function gradient");
	}


	if(outvar[26]) // Surface normal
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].norm(dim);

			string name = "Surface normal" + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}

	if(outvar[27])
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end; ++ii)
			vec[ii-start] = pnp1[ii].vPert.norm();
		Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Shifting velocity magnitude");
	}

	if(outvar[28])
	{
		for(uint dim = 0; dim < SIMDIM; ++dim)
		{
			#pragma omp parallel for
			for(size_t ii = start; ii < end; ++ii)
				vec[ii-start] = pnp1[ii].vPert(dim);

			string name = "Shifting velocity " + std::to_string(dim);
			Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
		}
	}
}


void Write_Binary_Timestep(SIM const& svar, real const& rho0, SPHState const& pnp1,
	bound_block const& limits, char const* group, int32_t const& strandID, void* const& fileHandle)
{
	size_t start = limits.index.first;
	size_t end = limits.index.second;
	int32_t const size = end - start;

	double solTime = svar.t;     
	int32_t outputZone;

	// Get zone data types
	vector<int32_t> varTypes = svar.var_types;
	vector<int32_t> shareVarFromZone(varTypes.size(),0);
    vector<int32_t> valueLocation(varTypes.size(),1);
    vector<int32_t> passiveVarList(varTypes.size(),0);

	if(tecZoneCreateIJK(fileHandle,group,size,1,1,&varTypes[0],
		&shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&outputZone))
	{
		cerr << "Failed to create IJK zone." << endl;
		exit(-1);
	}

	if(strandID != 0)
		if(tecZoneSetUnsteadyOptions(fileHandle, outputZone, solTime, strandID))
		{
			cerr << "Failed to add unsteady options." << endl;
			exit(-1);
		}

	Write_Zone(pnp1,rho0,svar.scale,svar.outvar,start,end,fileHandle,outputZone);
	
	// Write some auxiliary data, such as mass
	add_zone_aux_data(fileHandle,outputZone,"Particle mass",pnp1[start].m);
	// Add information about the insertion, aero, and deletion planes.
	add_zone_aux_data(fileHandle,outputZone,"Insertion normal x",limits.insert_norm[0]);
	add_zone_aux_data(fileHandle,outputZone,"Insertion normal y",limits.insert_norm[1]);
	add_zone_aux_data(fileHandle,outputZone,"Deletion normal x",limits.delete_norm[0]);
	add_zone_aux_data(fileHandle,outputZone,"Deletion normal y",limits.delete_norm[1]);
	add_zone_aux_data(fileHandle,outputZone,"Aerodynamic normal x",limits.aero_norm[0]);
	add_zone_aux_data(fileHandle,outputZone,"Aerodynamic normal y",limits.aero_norm[1]);
	#if SIMDIM == 3
	add_zone_aux_data(fileHandle,outputZone,"Insertion normal z",limits.insert_norm[2]);
	add_zone_aux_data(fileHandle,outputZone,"Deletion normal z",limits.delete_norm[2]);
	add_zone_aux_data(fileHandle,outputZone,"Aerodynamic normal z",limits.aero_norm[2]);
	#endif
	
	add_zone_aux_data(fileHandle,outputZone,"Insertion plane constant",limits.insconst);
	add_zone_aux_data(fileHandle,outputZone,"Deletion plane constant",limits.delconst);
	add_zone_aux_data(fileHandle,outputZone,"Aerodynamic plane constant",limits.aeroconst);

	// Add boundary time information
	if(limits.nTimes != 0)
	{
		add_zone_aux_data(fileHandle,outputZone,"Times count",limits.nTimes);
		for(size_t time = 0; time < limits.nTimes; time++)
		{
			add_zone_aux_data(fileHandle,outputZone,"Time " + std::to_string(time),limits.times[time]);
			add_zone_aux_data(fileHandle,outputZone,"Velocity " + std::to_string(time) + " x",limits.vels[time]);
			add_zone_aux_data(fileHandle,outputZone,"Velocity " + std::to_string(time) + " y",limits.vels[time]);
			#if SIMDIM == 3
			add_zone_aux_data(fileHandle,outputZone,"Velocity " + std::to_string(time) + " z",limits.vels[time]);
			#endif
		}
	}
	else
	{
		add_zone_aux_data(fileHandle,outputZone,"Times count",0);
		add_zone_aux_data(fileHandle,outputZone,"Velocity x",limits.vels[0]);
		add_zone_aux_data(fileHandle,outputZone,"Velocity y",limits.vels[1]);
		#if SIMDIM == 3
		add_zone_aux_data(fileHandle,outputZone,"Velocity z",limits.vels[2]);
		#endif

	}

	add_zone_aux_data(fileHandle,outputZone,"Fixed velocity or dynamic inlet",limits.fixed_vel_or_dynamic);
	add_zone_aux_data(fileHandle,outputZone,"Lattice or HCPL packing",limits.hcpl);
	add_zone_aux_data(fileHandle,outputZone,"Boundary solver type",limits.bound_solver);
	add_zone_aux_data(fileHandle,outputZone,"Boundary is no slip",limits.no_slip);
	add_zone_aux_data(fileHandle,outputZone,"Block type",limits.block_type);

    // cout << "Flushing results." << endl;
    // INTEGER4 numZonesToRetain = 0;
    if(tecFileWriterFlush(fileHandle,0,NULL))
    {
    	cout << "Failed to flush data. Retrying..." << endl;
    	std::this_thread::sleep_for(std::chrono::milliseconds(10));
    	// retry the flush
    	if(tecFileWriterFlush(fileHandle,0,NULL))
    	{
	    	cerr << "Failed to flush data to file: " << fileHandle << endl;
	    	exit(-1);
	   	}
    }
}

void Init_Binary_PLT(SIM &svar, string const& filename, string const& zoneName, void* &fileHandle)
{
    int32_t FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
    int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT

    string file = svar.output_prefix + filename;

	if(tecFileWriterOpen(file.c_str(),zoneName.c_str(),svar.var_names.c_str(),fileFormat,FileType,1,NULL,&fileHandle))
    {
    	cout << "Failed to open " << file << endl;
    	exit(-1);
    }

	#ifdef DEBUG
	if(tecFileSetDiagnosticsLevel(fileHandle, 1))
	{
		cerr << "Failed to set debug option for output file: " << file << endl;
		exit(-1);
	}
	#endif

	// Write crucial simulation auxiliary data. Potential to expand
	add_file_aux_data(fileHandle,"Number of boundary blocks",svar.nbound);
	add_file_aux_data(fileHandle,"Number of fluid blocks",svar.nfluid);
	add_file_aux_data(fileHandle,"Number of boundary points",svar.bndPts);
	add_file_aux_data(fileHandle,"Number of fluid points",svar.simPts);
	add_file_aux_data(fileHandle,"Current frame", svar.frame);
	
}

/*************************************************************************/
/**************************** BINARY INPUTS ******************************/
/*************************************************************************/
inline void Read_Real_Value(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
					int64_t const& iMax, real& var, string const& varName)
{
	int retval;
	vector<real> varvec(iMax);
	#if FOD == 1
		retval = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &varvec[0]);
	#else
		retval = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &varvec[0]);
	#endif

	if(retval)
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	var = varvec[0];

	++varCount;
}

inline void Read_Int_Value(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
					int64_t const& iMax, int& var, string const& varName)
{
	vector<int32_t> varvec(iMax);
	if(tecZoneVarGetInt32Values(inputHandle, frame, varCount, 1, iMax, &varvec[0]))
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	var = static_cast<int>(varvec[0]);
	++varCount;
}

inline void Read_UInt_Value(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
					int64_t const& iMax, uint& var, string const& varName)
{
	vector<uint8_t> varvec(iMax);
	if(tecZoneVarGetUInt8Values(inputHandle, frame, varCount, 1, iMax, &varvec[0]))
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	var = static_cast<uint>(varvec[0]);
	++varCount;
}

inline void Read_Real_Vector(void* const& inputHandle, int32_t const& frame, int32_t const& varCount, 
					int64_t const& iMax, vector<real>& var, string const& varName)
{
	int retval;
	var = vector<real>(iMax);
	#if FOD == 1
		retval = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &var[0]);
	#else
		retval = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &var[0]);
	#endif

	if(retval)
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
}

inline void Read_Int_Vector(void* const& inputHandle, int32_t const& frame, int32_t const& varCount, 
					int64_t const& iMax, vector<int>& var, string const& varName)
{
	vector<int32_t> varvec(iMax);
	if(tecZoneVarGetInt32Values(inputHandle, frame, varCount, 1, iMax, &varvec[0]))
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	var = vector<int>(iMax);
	for(int64_t ii = 0; ii < iMax; ++ii)
	{
		var[ii] = static_cast<int>(varvec[ii]);
	}
}

inline void Read_UInt_Vector(void* const& inputHandle, int32_t const& frame, int32_t const& varCount, 
					int64_t const& iMax, vector<uint>& var, string const& varName)
{
	vector<uint8_t> varvec(iMax);
	if(tecZoneVarGetUInt8Values(inputHandle, frame, varCount, 1, iMax, &varvec[0]))
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	var = vector<uint>(iMax);
	for(int64_t ii = 0; ii < iMax; ++ii)
	{
		var[ii] = static_cast<uint>(varvec[ii]);
	}
}

inline vector<StateVecD> Read_Binary_Vector(void* inputHandle, int32_t& frame, 
		int32_t& varCount, int64_t& iMax, string const& varName)
{
	vector<real> xVar(iMax,0.0);
	vector<real> yVar(iMax,0.0);

    // cout << "Trying to get vector x-component. Var: " << varCount << endl;
	string name = "x-component of ";
	name.append(varName);
	Read_Real_Vector(inputHandle, frame, varCount, iMax, xVar, name);
	
    // cout << "Trying to get vector y-component. Var: " << varCount << endl;
	name = "y-component of ";
	name.append(varName);
	Read_Real_Vector(inputHandle, frame, varCount, iMax, yVar, name);

	#if SIMDIM == 3
	vector<real> zVar(iMax,0.0);
	// cout << "Trying to get vector z-component. Var: " << varCount << endl;
	name = "z-component of ";
	name.append(varName);
	Read_Real_Vector(inputHandle, frame, varCount, iMax, zVar, name);
	#endif

	vector<StateVecD> vec(iMax);
	#pragma omp parallel for
	for(uint ii = 0; ii < iMax; ++ii)
	{
		vec[ii](0) = xVar[ii];
		vec[ii](1) = yVar[ii];
		
	#if SIMDIM == 3
		vec[ii](2) = zVar[ii];
	#endif

		// cout << tsData.verts[ii](0) << "  " << tsData.verts[ii](1) << endl;
	}

	return vec;
}

void Read_Binary_Timestep(void* inputHandle, SIM& svar, FLUID const& fvar, int32_t frame, 
			vector<int32_t> const& varIndex, SPHState& pn, bound_block& limits)
{
	real mass = 0.0;

	cout << "Reading zone number: " << frame << endl;

	int64_t iMax, jMax, kMax;
	if(tecZoneGetIJK(inputHandle, frame, &iMax, &jMax, &kMax))
	{
		cout << "Failed to read frame " << frame << " IJK info. Stopping." << endl;
		exit(-1);
	}

	cout << "Number of particles in zone: " << iMax << endl << endl;


	if(tecZoneGetSolutionTime(inputHandle, frame, &svar.t))
	{
		cerr << "Failed to read frame " << frame << " time" << endl;
		exit(-1);
	}

	// Read auxiliary data to find the mass. For now it's the only aux data.
	char* auxname = NULL;
	char* auxvalue = NULL;
	std::istringstream stream;
	int32_t numAux;
	if(tecZoneAuxDataGetNumItems(inputHandle,frame,&numAux))
	{
		std::cout << "Failed to get the number of auxiliary items, stopping." << std::endl;
		exit(-1);
	}

	for(int32_t aux = 1; aux <= numAux; aux++)
	{
		real data;
		if(tecZoneAuxDataGetItem(inputHandle, frame, aux, &auxname, &auxvalue))
		{
			std::cout << "Failed to read boundary auxiliary data number " << aux << ". Stopping" << std::endl;
			exit(-1); 
		}
		
		stream = std::istringstream(auxvalue);
		stream >> data;

		if(auxname == space2underscore("Particle mass"))
			mass = data;

		if(auxname == space2underscore("Insertion normal x"))
			limits.insert_norm[0] = data;
			
		if(auxname == space2underscore("Insertion normal y"))
			limits.insert_norm[1] = data;
			
		if(auxname == space2underscore("Insertion normal z"))
			limits.insert_norm[2] = data;
			
		if(auxname == space2underscore("Deletion normal x"))
			limits.delete_norm[0] = data;
			
		if(auxname == space2underscore("Deletion normal y"))
			limits.delete_norm[1] = data;
			
		if(auxname == space2underscore("Deletion normal z"))
			limits.delete_norm[2] = data;
			
		if(auxname == space2underscore("Aerodynamic normal x"))
			limits.aero_norm[0] = data;
			
		if(auxname == space2underscore("Aerodynamic normal y"))
			limits.aero_norm[1] = data;
			
		if(auxname == space2underscore("Aerodynamic normal z"))
			limits.aero_norm[2] = data;
			
		if(auxname == space2underscore("Insertion plane constant"))
			limits.insconst = data;
			
		if(auxname == space2underscore("Deletion plane constant"))
			limits.delconst = data;
			
		if(auxname == space2underscore("Aerodynamic plane constant"))
			limits.aeroconst = data;

		if(auxname == space2underscore("Times count"))
			limits.nTimes = data;

		if(auxname == space2underscore("Fixed velocity or dynamic inlet"))
			limits.fixed_vel_or_dynamic = data;

		if(auxname == space2underscore("Lattice or HCPL packing"))
			limits.hcpl = data;

		if(auxname == space2underscore("Boundary solver type"))
			limits.bound_solver = data;

		if(auxname == space2underscore("Boundary is no slip"))
			limits.no_slip = data;

		if(auxname == space2underscore("Block type"))
			limits.block_type = data;

		tecStringFree(&auxname);
		tecStringFree(&auxvalue);
	}

	// If the number of times is non zero, then resize to have the right size
	if(limits.nTimes != 0)
	{
		limits.times.resize(limits.nTimes);
		limits.vels.resize(limits.nTimes);
	}
	else
		limits.vels.resize(1);
		
	// Get velocity data this time
	for(int32_t aux = 1; aux <= numAux; aux++)
	{
		real data;
		if(tecZoneAuxDataGetItem(inputHandle, frame, aux, &auxname, &auxvalue))
		{
			std::cout << "Failed to read boundary auxiliary data number " << aux << ". Stopping" << std::endl;
			exit(-1); 
		}
		std::string name(auxname);
		
		stream = std::istringstream(auxvalue);
		stream >> data;

		if(name.find(space2underscore("Time ")) != string::npos)
		{
			if(name == "Times count")
			{
				tecStringFree(&auxname);
				tecStringFree(&auxvalue);
				continue;
			}

			std::istringstream str(name.substr(5));
			size_t jj;
			str >> jj;

			limits.times[jj] = data;
		}

		if(name.find(space2underscore("Velocity ")) != std::string::npos)
		{
			// Find the count, and know the fixed width of the preamble
			if(name.size() == 9)
			{	// Avoid the false positive.
				if(name.back() == 'x')
					limits.vels[0][0] = data;
				else if (name.back() == 'y')
					limits.vels[0][1] = data;
				#if SIMDIM == 3
				else if (name.back() == 'z')
					limits.vels[0][2] = data;
				#endif	

				tecStringFree(&auxname);
				tecStringFree(&auxvalue);
				continue;
			}

			std::istringstream str(name.substr(9));
			size_t jj;
			str >> jj;

			// Now find which dimension it refers to.
			if(name.back() == 'x')
				limits.vels[jj][0] = data;
			else if (name.back() == 'y')
				limits.vels[jj][1] = data;
			#if SIMDIM == 3
			else if (name.back() == 'z')
				limits.vels[jj][2] = data;
			#endif		
		}

		tecStringFree(&auxname);
		tecStringFree(&auxvalue);
	}


	if(mass == 0.0)
	{
		cout << "Failed to retrieve the particle mass from the auxiliary data. Stopping" << endl;
		exit(-1);
	}

	vector<real> x(iMax);
	vector<real> y(iMax);
	vector<real> ax(iMax);
	vector<real> ay(iMax);
	vector<real> vx(iMax);
	vector<real> vy(iMax);
	#if SIMDIM == 3
	vector<real> z(iMax);
	vector<real> az(iMax);
	vector<real> vz(iMax);
	#endif
	vector<real> press(iMax,0.0);
	vector<real> Rrho(iMax,0.0);
	vector<int>  partID(iMax);
	vector<int>  cellID(iMax);
	vector<uint> b(iMax,0);

	Read_Real_Vector(inputHandle, frame, varIndex[0], iMax, x, "position x");
	Read_Real_Vector(inputHandle, frame, varIndex[1], iMax, y, "position y");
	#if SIMDIM == 3
	Read_Real_Vector(inputHandle, frame, varIndex[11], iMax, z, "position z");
	#endif

	Read_Real_Vector(inputHandle, frame, varIndex[2], iMax, ax, "acceleration x");
	Read_Real_Vector(inputHandle, frame, varIndex[3], iMax, ay, "acceleration y");
	#if SIMDIM == 3
	Read_Real_Vector(inputHandle, frame, varIndex[12], iMax, az, "acceleration z");
	#endif
	Read_Real_Vector(inputHandle, frame, varIndex[4], iMax, vx, "velocity x");
	Read_Real_Vector(inputHandle, frame, varIndex[5], iMax, vy, "velocity y");
	#if SIMDIM == 3
	Read_Real_Vector(inputHandle, frame, varIndex[13], iMax, vz, "velocity z");
	#endif

	Read_Real_Vector(inputHandle, frame, varIndex[6], iMax, press, "pressure");
	Read_Real_Vector(inputHandle, frame, varIndex[7], iMax, Rrho, "density gradient");
	Read_Int_Vector(inputHandle, frame, varIndex[8], iMax, partID, "particle ID");
	Read_Int_Vector(inputHandle, frame, varIndex[9], iMax, cellID, "cell ID");
	Read_UInt_Vector(inputHandle, frame, varIndex[10], iMax, b, "particle type");

	pn = vector<SPHPart>(iMax);
	/*Now put it into the state vector*/
	for(int64_t ii = 0; ii < iMax; ++ii)
	{
		pn[ii].xi[0]= x[ii];
		pn[ii].xi[1]= y[ii];
		pn[ii].acc[0] = ax[ii];
		pn[ii].acc[1] = ay[ii];
		pn[ii].v[0] = vx[ii];
		pn[ii].v[1] = vy[ii];

		#if SIMDIM == 3
		pn[ii].xi[2]= z[ii];
		pn[ii].acc[2] = az[ii];
		pn[ii].v[2] = vz[ii];
		#endif
		
		pn[ii].p = press[ii];
		// pn[ii].rho = density_equation(press[ii],fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
		pn[ii].Rrho = Rrho[ii];
		pn[ii].m = mass;
		pn[ii].b = b[ii]; // static_cast<uint>(b[ii]);
		pn[ii].partID = partID[ii];
		pn[ii].cellID = static_cast<size_t>(cellID[ii]);	
	} 
}


void CheckContents_Binary(void* const& inputHandle, SIM& svar, int32_t& numZones, double& time, 
			vector<int32_t>& vecIndex)
{
	int32_t numVars;
	if(tecDataSetGetNumVars(inputHandle, &numVars))
	{
		cout << "Couldn't obtain number of variables. Stopping" << endl;
		exit(-1);
	}
	
	cout << "Number of variables in output file: " << numVars << endl;

	string axis1 = "X"; string axis1s = "x";
	string axis2 = "Y"; string axis2s = "y";
	#if SIMDIM == 2
	if(svar.offset_axis != 0)
	{
		if(svar.offset_axis == 1)
		{
			axis1 = "Y"; axis1s = "y";
			axis2 = "Z"; axis2s = "z";
		}
		else if(svar.offset_axis == 2)
		{
			axis1 = "X"; axis1s = "x";
			axis2 = "Z"; axis2s = "z";
		}
		else if(svar.offset_axis == 3)
		{
			axis1 = "X"; axis1s = "x";
			axis2 = "Y"; axis2s = "y";
		}
	}
	#endif

	// Check the restart data exists
	vecIndex = vector<int32_t>(3*SIMDIM + 5, 0);
	vector<uint> hasVar(3*SIMDIM + 5, 0);
	for (int32_t var = 1; var <= numVars; ++var)
    {
        char* name = NULL;
        if(tecVarGetName(inputHandle, var, &name))
		{
			cout << "Failed to get variable name" << endl;
			exit(-1);
		}
		string sname(name);

		if(sname == axis1)
		{
			vecIndex[0] = var;
			hasVar[0] = realType;
		}
		if(sname == axis2)
		{
			vecIndex[1] = var;
			hasVar[1] = realType;
		}
		if(sname == "A" + axis1s)
		{
			vecIndex[2] = var;
			hasVar[2] = realType;
		}
		if(sname == "A" + axis2s)
		{
			vecIndex[3] = var;
			hasVar[3] = realType;
		}
		if(sname == "V" + axis1s)
		{
			vecIndex[4] = var;
			hasVar[4] = realType;
		}
		if(sname == "V" + axis2s)
		{
			vecIndex[5] = var;
			hasVar[5] = realType;
		}
		if(sname == "Pressure")
		{
			vecIndex[6] = var;
			hasVar[6] = realType;
		}
		if(sname == "dRho")
		{
			vecIndex[7] = var;
			hasVar[7] = realType;
		}
		if(sname == "partID")
		{
			vecIndex[8] = var;
			hasVar[8] = 3;
		}
		if(sname == "cellID")
		{
			vecIndex[9] = var;
			hasVar[9] = 3;
		}
		if(sname == "bound")
		{
			vecIndex[10] = var;
			hasVar[10] = 5;
		}
		#if SIMDIM == 3
		if(sname == "Z")
		{
			vecIndex[11] = var;
			hasVar[11] = realType;
		}
		if(sname == "Az")
		{
			vecIndex[12] = var;
			hasVar[12] = realType;
		}
		if(sname == "Vz")
		{
			vecIndex[13] = var;
			hasVar[13] = realType;
		}
		#endif
		
        tecStringFree(&name);
	}

	for(uint const& var:hasVar)
	{
		if(var == 0)
		{
			cout << "Failed to find a variable required for restart in the file. Stopping." << endl;
			exit(-1);
		}
	}

    if(tecDataSetGetNumZones(inputHandle, &numZones))
	{
		cout << "Failed to get number of zones in file" << endl;
		exit(-1);
	}
    cout << "Number of zones in file: " << numZones << endl;

    if(tecZoneGetSolutionTime(inputHandle, numZones, &time))
	{
		cout << "Failed to get solution time of zone: " << numZones << endl;
		exit(-1);
	}
    cout << "Latest zone time: " << time << endl << endl;

	// Get auxiliary data
	char* auxname = NULL;
	char* auxvalue = NULL;
	std::istringstream stream;
	int32_t numAux;
	if(tecDataSetAuxDataGetNumItems(inputHandle,&numAux))
	{
		std::cout << "Failed to get the number of auxiliary items, stopping." << std::endl;
		exit(-1);
	}

	for(int32_t aux = 1; aux <= numAux; aux++)
	{
		size_t data;
		if(tecDataSetAuxDataGetItem(inputHandle, aux, &auxname, &auxvalue))
		{
			std::cout << "Failed to read dataset auxiliary data. Stopping" << std::endl;
			exit(-1); 
		}
		stream = std::istringstream(auxvalue);
		stream >> data;

		if(auxname == space2underscore("Number of boundary blocks"))
			svar.nbound = data;
		
		if(auxname == space2underscore("Number of fluid blocks"))
			svar.nfluid = data;

		if(auxname == space2underscore("Number of boundary points"))
			svar.bndPts = data;

		if(auxname == space2underscore("Number of fluid points"))
			svar.simPts = data;
			
		tecStringFree(&auxname);
		tecStringFree(&auxvalue);
	}
}

void Restart_Binary(SIM& svar, FLUID const& fvar, SPHState& pn, LIMITS& limits)
{	
	// Read the values from the solution folder, then check. 
	#ifdef DEBUG
		fprintf(dbout,"Reading binary files for restart information.\n");
	#endif

	void* fluidHandle = NULL;
	void* boundHandle = NULL;

	string dir;
	string prefix = svar.restart_prefix;
	if(svar.restart_prefix.find_last_of("/\\") != string::npos)
		dir = svar.restart_prefix.substr(0,svar.restart_prefix.find_last_of("/\\"));
	else
	{
		dir = ".";
		prefix = "./" + svar.restart_prefix;
	}
	
	string pltext = ".szplt";
	string datext = ".szdat";

	std::vector<std::filesystem::directory_entry> pltfiles;
	std::vector<std::filesystem::directory_entry> datfiles;
	for(auto const& file : directory_iterator(dir))
	{
		// cout << file.path().string() << endl;
		// Check it has the right prefix
		if(file.path().string().find(prefix) == string::npos)
			continue;
		
		// Now check if a szplt already exists, and if the data exists
		if(file.path().extension() == pltext)
			pltfiles.emplace_back(file);
		else if(file.path().extension() == datext)
			datfiles.emplace_back(file);
	}

	float latestPltTime = 0.0;
	float latestDatTime = 0.0;
	float latestBPltTime = 0.0;
	float latestBDatTime = 0.0;
	string single_plt_file, multi_plt_file, bound_plt_file, bound_multi_plt_file;
	string single_dat_file, multi_dat_file, bound_dat_file, bound_multi_dat_file;

	std::filesystem::file_time_type plt_write_time_single, plt_write_time_multi, 
								plt_write_time_bsing, plt_write_time_bmulti;

	std::filesystem::file_time_type dat_write_time_single, dat_write_time_multi,
							 	dat_write_time_bsing, dat_write_time_bmulti;

	if(!pltfiles.empty())
	{	// Find the most recent timestep. Can do it off of write time or frame number?
		for(std::filesystem::directory_entry const& file:pltfiles)
		{
			string name = file.path().string();
			// cout << "Name: " << name << endl;
			if(name.find("time") != string::npos)
			{	//Multiple files are being used, so extract from the most recent
				size_t i1 = name.find("time") + 5;
				size_t i2 = name.find("_fluid");

				if(i2 != string::npos)
				{
					string substr = name.substr(i1,i2);
					std::istringstream timestr(substr);
					float time;
					timestr >> time;
					if(time > latestPltTime)
					{
						latestPltTime = time;
						plt_write_time_multi = file.last_write_time();
						multi_plt_file = name;
					}
				}
				else
				{
					i2 = name.find("_boundary");
					if(i2 != string::npos)
					{
						string substr = name.substr(i1,i2);
						std::istringstream timestr(substr);
						float time;
						timestr >> time;
						if(time > latestPltTime)
						{
							latestBPltTime = time;
							plt_write_time_bmulti = file.last_write_time();
							bound_multi_plt_file = name;
						}
					}
				}
			}
			else if(name == (prefix + "_fluid.szplt"))
			{
				single_plt_file = name;
				plt_write_time_single = file.last_write_time();
			}
			else if(name == (prefix + "_boundary.szplt"))
			{
				bound_plt_file = name;
				dat_write_time_bsing = file.last_write_time();
			}
		}
	}

	if(!datfiles.empty())
	{	// Find the most recent timestep. Can do it off of write time or frame number?
		for(std::filesystem::directory_entry const& file:datfiles)
		{
			string name = file.path().string();
			
			if(name.find("time") != string::npos)
			{	//Multiple files are being used, so extract from the most recent
				size_t i1 = name.find("time") + 5;
				size_t i2 = name.find("_fluid");

				if(i2 != string::npos)
				{
					string substr = name.substr(i1,i2);
					std::istringstream timestr(substr);
					double time;
					timestr >> time;
					if(time > latestDatTime)
					{
						latestDatTime = time;
						dat_write_time_multi = file.last_write_time();
						multi_dat_file = name;
					}
				}
				else
				{
					i2 = name.find("_boundary");
					if(i2 != string::npos)
					{
						string substr = name.substr(i1,i2);
						std::istringstream timestr(substr);
						double time;
						timestr >> time;
						if(time > latestBDatTime)
						{
							latestBDatTime = time;
							dat_write_time_bmulti = file.last_write_time();
							bound_multi_dat_file = name;
						}
					}
				}
			}
			else if(name == (prefix + "_fluid.szplt.szdat"))
			{
				single_dat_file = name;
				dat_write_time_single = file.last_write_time();
			}
			else if(name == (prefix + "_boundary.szplt.szdat"))
			{
				bound_dat_file = name;
				dat_write_time_bsing = file.last_write_time();
			}
				
		}
	}

	string fluid_restart_file;
	string bound_restart_file;

	if(svar.single_file)
	{
		// See if the plt file needs combining
		if(!single_plt_file.empty() && !single_dat_file.empty())
		{
			if(dat_write_time_single > plt_write_time_single)
			{	// More recent, so combine it
				string name = single_dat_file.substr(0,single_dat_file.size()-6);
				if(Combine_SZPLT(name) == -1)
					exit(-1);
				fluid_restart_file = name;
			}
			else
				fluid_restart_file = single_plt_file;
		}
		else if (!single_plt_file.empty())
		{	// Data files don't exist, but the szplt file does
			fluid_restart_file = single_plt_file;
			std::cout << "WARNING: No data files exist, so time history will be lost." << std::endl;
		}
		else if (!single_dat_file.empty())
		{	// szplt file doesnt exist yet, so create it
			string name = single_dat_file.substr(0,single_dat_file.size()-6);
			if(Combine_SZPLT(name) == -1)
				exit(-1);
			fluid_restart_file = name;
		}
		else
		{
			std::cout << "No files capable of restarting from exist. Stopping." << std::endl;
			exit(-1);
		}

		// Do the same for the boundary file
		if(!bound_plt_file.empty() && !bound_dat_file.empty())
		{
			if(dat_write_time_bsing > plt_write_time_bsing)
			{	// More recent, so combine it
				string name = bound_dat_file.substr(0,bound_dat_file.size()-6);
				if(Combine_SZPLT(name) == -1)
					exit(-1);
				bound_restart_file = name;
			}
			else
				bound_restart_file = bound_plt_file;
		}
		else if (!bound_plt_file.empty())
		{	// Data files don't exist, but the szplt file does
			bound_restart_file = bound_plt_file;
			std::cout << "WARNING: No data files exist, so time history will be lost." << std::endl;
		}
		else if (!bound_dat_file.empty())
		{	// szplt file doesnt exist yet, so create it
			string name = bound_dat_file.substr(0,bound_dat_file.size()-6);
			if(Combine_SZPLT(name) == -1)
				exit(-1);
			bound_restart_file = name;
		}
		else
		{
			std::cout << "No boundary files exist." << std::endl;
		}

		int32_t fluidFrames, boundFrames;
		double fluidTime, boundTime;

		if(tecFileReaderOpen(fluid_restart_file.c_str(),&fluidHandle))
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << fluid_restart_file << endl;
			exit(-1);
		}

		// Check how many frames are in the fluid file.
		cout << "Checking Fuel file..." << endl;
		vector<int32_t> varIndex_fluid, varIndex_bound;
		CheckContents_Binary(fluidHandle,svar,fluidFrames,fluidTime,varIndex_fluid);

		if(fluidFrames % svar.nfluid != 0)
		{
			cout << "WARNING: Number of fluid frames is not a factor of the defined number of fluid blocks." << endl;
		}

		limits.resize(svar.nbound+svar.nfluid,0);

		if(svar.nbound != 0)
		{
			if(!std::filesystem::exists(bound_restart_file))
			{
				cout << "ERROR: Simulation fluid file states existence of boundary blocks, but no boundary file exists." << endl;
				exit(-1);
			}
			
			if(tecFileReaderOpen(bound_restart_file.c_str(),&boundHandle))
			{
				cout << "Error opening szplt file. Path:" << endl;
				cout << bound_restart_file << endl;
				exit(-1);
			}

			cout << "Checking Boundary file..." << endl;
			CheckContents_Binary(boundHandle,svar,boundFrames,boundTime,varIndex_bound);
			
			if(boundFrames % svar.nbound != 0)
			{
				cout << "WARNING: Number of boundary frames is not a factor of the defined number of boundary blocks." << endl;
			}

			if(fluidTime != boundTime)
			{
				cout << "WARNING: Frame times are not consistent between fluid and boundary files." << endl;

				if(fluidTime > boundTime)
				{
					double time = 0.0;
					for(int32_t frame = fluidFrames-1; frame > 1; frame--)
					{
						if(tecZoneGetSolutionTime(fluidHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from fluid file." << endl;
							continue;
						}
						
						if(time == boundTime)
						{
							cout << "Found the correct frame" << endl;
							fluidFrames = frame;
							fluidTime = time;
							break;
						}
					}
				}
				else
				{
					double time = 0.0;
					for(int32_t frame = boundFrames-1; frame >= 1; frame--)
					{
						if(tecZoneGetSolutionTime(boundHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from boundary file." << endl;
							continue;
						}

						if(time == fluidTime)
						{
							cout << "Found the correct frame" << endl;
							boundFrames = frame;
							boundTime = time;
							break;
						}
					}
				}

				if(fluidTime != boundTime)
				{
					cout << "Could not find a consistent time in each file. Stopping." << endl;
					exit(-1);
				}
			}

			cout << "Attempting to read the boundary..." << endl;
			size_t nBnd = 0;
			for(size_t ii = 0; ii < svar.nbound; ++ii)
			{
				int32_t frame = boundFrames - svar.nbound + ii + 1;
				SPHState block;
				Read_Binary_Timestep(boundHandle,svar,fvar,frame,varIndex_bound,block,limits[ii]);
				pn.insert(pn.end(),block.begin(),block.end());
				nBnd += block.size();
			}
			svar.bndPts = nBnd;
		}

		svar.frame = fluidFrames/svar.nfluid - 1;

		// Read the actual data.
		cout  << "Attempting to read the fluid..." << endl;
		size_t nFlu = 0;
		for(size_t ii = 0; ii < svar.nfluid; ++ii)
		{
			int32_t frame = fluidFrames - svar.nfluid + ii + 1;
			SPHState block;
			Read_Binary_Timestep(fluidHandle,svar,fvar,frame,varIndex_fluid,block,limits[ii+svar.nbound]);
			pn.insert(pn.end(),block.begin(),block.end());
			nFlu += block.size();
		}
		svar.simPts = nFlu;
	}	// End single file
	else
	{
		if(!multi_plt_file.empty() && !multi_dat_file.empty())
		{	// Both file exist, so decide if the data needs combining
			if(latestDatTime > latestPltTime)
			{
				// Try combine the file, but may not be fully written, since it should be combined
				string name = multi_dat_file.substr(0,multi_dat_file.size()-6);
				if(Combine_SZPLT(name) == -1)
				{
					std::cout << "Failed to combine multi file szlplt data. Using latest combined data." << std::endl;
					fluid_restart_file = multi_plt_file;
				}	
				else
					fluid_restart_file = name;
			}
		}
		else if (!multi_plt_file.empty())
		{	// Data files don't exist, but the szplt file does
			fluid_restart_file = multi_plt_file;
		}
		else if (!multi_dat_file.empty())
		{	// szplt file doesnt exist yet, so create it
			string name = multi_dat_file.substr(0,multi_dat_file.size()-6);
			if(Combine_SZPLT(name) == -1)
				exit(-1);
			fluid_restart_file = name;
		}
		else
		{
			std::cout << "No files capable of restarting from exist. Stopping." << std::endl;
			exit(-1);
		}

		// Do the same for boundary files
		if(!bound_multi_plt_file.empty() && !bound_multi_dat_file.empty())
		{	// Both file exist, so decide if the data needs combining
			if(latestBDatTime > latestBPltTime)
			{
				// Try combine the file, but may not be fully written, since it should be combined
				string name = bound_multi_dat_file.substr(0,bound_multi_dat_file.size()-6);
				if(Combine_SZPLT(name) == -1)
				{
					std::cout << "Failed to combine multi file szlplt data. Using latest combined data." << std::endl;
					fluid_restart_file = bound_multi_plt_file;
				}	
				else
					fluid_restart_file = name;
			}
		}
		else if (!bound_multi_plt_file.empty())
		{	// Data files don't exist, but the szplt file does
			fluid_restart_file = bound_multi_plt_file;
		}
		else if (!multi_dat_file.empty())
		{	// szplt file doesnt exist yet, so create it
			string name = multi_dat_file.substr(0,bound_multi_dat_file.size()-6);
			if(Combine_SZPLT(name) == -1)
				exit(-1);
			fluid_restart_file = name;
		}
		else
		{
			std::cout << "No boundary files exist." << std::endl;
		}
		
		int32_t fluidFrames, boundFrames;
		double fluidTime, boundTime;

		if(tecFileReaderOpen(fluid_restart_file.c_str(),&fluidHandle))
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << fluid_restart_file << endl;
			exit(-1);
		}

		// Check how many frames are in the fluid file.
		cout << "Checking Fuel file..." << endl;
		vector<int32_t> varIndex_fluid, varIndex_bound;
		CheckContents_Binary(fluidHandle,svar,fluidFrames,fluidTime,varIndex_fluid);

		if(static_cast<int>(svar.nfluid) != fluidFrames)
		{
			cout << "ERROR: Mismatch of number of defined fluid blocks, and number of defined zones in file." << endl;
			exit(-1);
		}

		limits.resize(svar.nbound+svar.nfluid,0);

		if(svar.nbound != 0)
		{
			if(!std::filesystem::exists(bound_restart_file))
			{
				cout << "ERROR: Simulation fluid file states existence of boundary blocks, but no boundary file exists." << endl;
				exit(-1);
			}
			
			if(tecFileReaderOpen(bound_restart_file.c_str(),&boundHandle))
			{
				cout << "Error opening szplt file. Path:" << endl;
				cout << bound_restart_file << endl;
				exit(-1);
			}

			cout << "Checking Boundary file..." << endl;
			CheckContents_Binary(boundHandle,svar,boundFrames,boundTime,varIndex_bound);
			
			if(static_cast<int>(svar.nbound) != boundFrames)
			{
				cout << "ERROR: Number of boundary frames is not a factor of the defined number of boundary blocks." << endl;
			}

			if(fluidTime != boundTime)
			{
				cout << "ERROR: Frame times are not consistent between fluid and boundary files." << endl;
				exit(-1);
			}

			cout << "Attempting to read the boundary..." << endl;
			size_t nBnd = 0;
			for(size_t ii = 0; ii < svar.nbound; ++ii)
			{
				int32_t frame = ii + 1;
				SPHState block;
				Read_Binary_Timestep(boundHandle,svar,fvar,frame,varIndex_bound,block,limits[ii]);
				pn.insert(pn.end(),block.begin(),block.end());
				nBnd += block.size();
			}
			svar.bndPts = nBnd;
		}

		
		// Read the actual data.
		cout  << "Attempting to read the fluid..." << endl;
		size_t nFlu = 0;
		for(size_t ii = 0; ii < svar.nfluid; ++ii)
		{
			int32_t frame = ii + 1;
			SPHState block;
			Read_Binary_Timestep(fluidHandle,svar,fvar,frame,varIndex_fluid,block,limits[ii+svar.nbound]);
			pn.insert(pn.end(),block.begin(),block.end());
			nFlu += block.size();
		}
		svar.simPts = nFlu;

	}	// End multi file

	/* Get time */
	svar.tframem1 = svar.t;

	svar.totPts = pn.size();
	
	if(svar.simPts + svar.bndPts != svar.totPts)
	{
		cout << "Mismatch of array sizes. Total array is not the sum of sim and boundary arrays" << endl;
		exit(-1);
	}	
}

void Write_Cell_Centres(MESH const& cells)
{
	void* fileHandle;
	int32_t FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
    int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT

    string file = "Cell_Centres.szplt";
	string zoneName = "Cell Centres";
	
	string variables = "X,Y,Z"; 
	vector<int32_t> varTypes;
	#if FOD == 1
		int32_t realType = 2;
	#else
		int32_t realType = 1;
	#endif

	varTypes.emplace_back(realType);
	varTypes.emplace_back(realType);
	varTypes.emplace_back(realType);

	if(tecFileWriterOpen(file.c_str(),zoneName.c_str(),variables.c_str(),
							fileFormat,FileType,1,NULL,&fileHandle))
    {
    	cout << "Failed to open " << file << endl;
    	exit(-1);
    }

	#ifdef DEBUG
	if(tecFileSetDiagnosticsLevel(fileHandle, 1))
	{
		cerr << "Failed to set debug option for output file: " << file << endl;
		exit(-1);
	}
	#endif

	int32_t outputZone;

	// Get zone data types
	vector<int32_t> shareVarFromZone(varTypes.size(),0);
    vector<int32_t> valueLocation(varTypes.size(),1);
    vector<int32_t> passiveVarList(varTypes.size(),0);

	if(tecZoneCreateIJK(fileHandle,zoneName.c_str(),cells.nElem,1,1,&varTypes[0],
		&shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&outputZone))
	{
		cerr << "Failed to create IJK zone." << endl;
		exit(-1);
	}


	/*Write the basic position data present for all outputs*/
	vector<real> x(cells.nElem);
	int32_t var = 1;

	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(uint ii = 0; ii < cells.nElem; ++ii)
			x[ii] = cells.cCentre[ii][dim];

		string name = "position coordinate ";
		name.append(std::to_string(dim));
		Write_Real_Vector(fileHandle, outputZone, var, 0, cells.nElem, x, name);
	}

	if(tecZoneCreateIJK(fileHandle,"Node coordinates",cells.nPnts,1,1,&varTypes[0],
		&shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&outputZone))
	{
		cerr << "Failed to create IJK zone." << endl;
		exit(-1);
	}


	/*Write the basic position data present for all outputs*/
	x = vector<real>(cells.nPnts);
	var = 1;

	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(uint ii = 0; ii < cells.nPnts; ++ii)
			x[ii] = cells.verts[ii][dim];

		string name = "position coordinate ";
		name.append(std::to_string(dim));
		Write_Real_Vector(fileHandle, outputZone, var, 0, cells.nPnts, x, name);
	}

	
    if(tecFileWriterClose(&fileHandle))
	{
		cout << "Failed to close file" << endl;
		exit(-1);
	}

}