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

	if(tecZoneCreateIJK(fileHandle,group, size,1,1,&varTypes[0],
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

	Write_Zone(pnp1,rho0, svar.scale, svar.outvar, start,end,fileHandle,outputZone);
	
	// Write some auxiliary data, such as mass
	add_zone_aux_data(fileHandle,outputZone, "Particle mass",pnp1[start].m);
	// Add information about the insertion, aero, and deletion planes.
	add_zone_aux_data(fileHandle,outputZone, "Insertion normal x",limits.insert_norm[0]);
	add_zone_aux_data(fileHandle,outputZone, "Insertion normal y",limits.insert_norm[1]);
	add_zone_aux_data(fileHandle,outputZone, "Deletion normal x",limits.delete_norm[0]);
	add_zone_aux_data(fileHandle,outputZone, "Deletion normal y",limits.delete_norm[1]);
	add_zone_aux_data(fileHandle,outputZone, "Aerodynamic normal x",limits.aero_norm[0]);
	add_zone_aux_data(fileHandle,outputZone, "Aerodynamic normal y",limits.aero_norm[1]);
	#if SIMDIM == 3
	add_zone_aux_data(fileHandle,outputZone, "Insertion normal z",limits.insert_norm[2]);
	add_zone_aux_data(fileHandle,outputZone, "Deletion normal z",limits.delete_norm[2]);
	add_zone_aux_data(fileHandle,outputZone, "Aerodynamic normal z",limits.aero_norm[2]);
	#endif
	
	add_zone_aux_data(fileHandle,outputZone, "Insertion plane constant",limits.insconst);
	add_zone_aux_data(fileHandle,outputZone, "Deletion plane constant",limits.delconst);
	add_zone_aux_data(fileHandle,outputZone, "Aerodynamic plane constant",limits.aeroconst);

	// Add boundary time information
	if(limits.nTimes != 0)
	{
		add_zone_aux_data(fileHandle,outputZone, "Times count",limits.nTimes);
		for(size_t time = 0; time < limits.nTimes; time++)
		{
			add_zone_aux_data(fileHandle,outputZone, "Time " + std::to_string(time),limits.times[time]);
			add_zone_aux_data(fileHandle,outputZone, "Velocity " + std::to_string(time) + " x",limits.vels[time][0]);
			add_zone_aux_data(fileHandle,outputZone, "Velocity " + std::to_string(time) + " y",limits.vels[time][1]);
			#if SIMDIM == 3
			add_zone_aux_data(fileHandle,outputZone, "Velocity " + std::to_string(time) + " z",limits.vels[time][2]);
			#endif
		}
	}
	else
	{
		add_zone_aux_data(fileHandle,outputZone, "Times count",0);
		add_zone_aux_data(fileHandle,outputZone, "Velocity x",limits.vels[0][0]);
		add_zone_aux_data(fileHandle,outputZone, "Velocity y",limits.vels[0][1]);
		#if SIMDIM == 3
		add_zone_aux_data(fileHandle,outputZone, "Velocity z",limits.vels[0][2]);
		#endif

	}

	add_zone_aux_data(fileHandle,outputZone, "Fixed velocity or dynamic inlet",limits.fixed_vel_or_dynamic);
	add_zone_aux_data(fileHandle,outputZone, "Lattice or HCPL packing",limits.hcpl);
	add_zone_aux_data(fileHandle,outputZone, "Boundary solver type",limits.bound_solver);
	add_zone_aux_data(fileHandle,outputZone, "Boundary is no slip",limits.no_slip);
	add_zone_aux_data(fileHandle,outputZone, "Block type",limits.block_type);

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

void Init_Binary_PLT(SIM &svar, FLUID const& fvar, AERO const& avar, string const& prefix,
			 string const& filename, string const& zoneName, void* &fileHandle)
{
    int32_t FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
    int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT

    string file = prefix + filename;

	if(tecFileWriterOpen(file.c_str(),zoneName.c_str(), svar.var_names.c_str(),fileFormat,FileType,1,NULL,&fileHandle))
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
	add_file_aux_data(fileHandle, "Number of boundary blocks", svar.nbound);
	add_file_aux_data(fileHandle, "Number of fluid blocks", svar.nfluid);
	add_file_aux_data(fileHandle, "Number of boundary points", svar.bndPts);
	add_file_aux_data(fileHandle, "Number of fluid points", svar.simPts);
	add_file_aux_data(fileHandle, "Current frame", svar.frame);

	add_file_aux_data(fileHandle, "Input para filename", svar.infile);
	add_file_aux_data(fileHandle, "Input fluid definition filename", svar.fluidfile);
	add_file_aux_data(fileHandle, "Input boundary definition filename", svar.boundfile);

	// add_file_aux_data(fileHandle, "Primary grid face filename", svar.taumesh);
	// add_file_aux_data(fileHandle, "Boundary mapping filename", svar.taubmap);
	// add_file_aux_data(fileHandle, "Restart-data prefix", svar.tausol);
	add_file_aux_data(fileHandle, "Dimension offset vector", svar.offset_axis);
	add_file_aux_data(fileHandle, "Solution angle of attack", svar.angle_alpha);
	add_file_aux_data(fileHandle, "Grid scale", svar.scale);

	add_file_aux_data(fileHandle, "OpenFOAM input directory", svar.foamdir);
	add_file_aux_data(fileHandle, "OpenFOAM solution directory", svar.foamsol);
	add_file_aux_data(fileHandle, "OpenFOAM binary", svar.isBinary);
	add_file_aux_data(fileHandle, "Label size", svar.labelSize);
	add_file_aux_data(fileHandle, "Scalar size", svar.scalarSize);
	add_file_aux_data(fileHandle, "OpenFOAM buoyant", svar.buoyantSim);

	add_file_aux_data(fileHandle, "VLM definition filename", svar.vlm_file); 

	add_file_aux_data(fileHandle, "Single file for output", svar.single_file);
	add_file_aux_data(fileHandle, "Output files prefix", svar.output_prefix);
	add_file_aux_data(fileHandle, "SPH restart prefix", svar.restart_prefix);
	add_file_aux_data(fileHandle, "SPH frame time interval", svar.framet);
	add_file_aux_data(fileHandle, "SPH frame count", svar.Nframe);
	add_file_aux_data(fileHandle, "SPH output encoding", svar.out_encoding);
	add_file_aux_data(fileHandle, "SPH ghost output", svar.gout);
	add_file_aux_data(fileHandle, "Variable list", svar.output_names);

	
	/* Fluid data */
	add_file_aux_data(fileHandle, "Reference density", avar.rhog);
	add_file_aux_data(fileHandle, "Reference dispersed density", fvar.rho0);
	add_file_aux_data(fileHandle, "Sutherland reference viscosity", avar.mug);
	add_file_aux_data(fileHandle, "Reference dispersed viscosity", fvar.mu);
	add_file_aux_data(fileHandle, "Reference surface tension", fvar.sig);
	add_file_aux_data(fileHandle, "SPH surface tension contact angle", fvar.contangb);
	add_file_aux_data(fileHandle, "Init hydrostatic pressure", svar.init_hydro_pressure);
	add_file_aux_data(fileHandle, "Hydrostatic height", svar.hydro_height);

	/* Aerodynamic data */
	add_file_aux_data(fileHandle, "Reference velocity", avar.vRef);
	add_file_aux_data(fileHandle, "Reference pressure", avar.pRef);
	add_file_aux_data(fileHandle, "Reference Mach number", avar.MRef);
	add_file_aux_data(fileHandle, "Reference temperature", avar.T);
	add_file_aux_data(fileHandle, "Gas constant gamma", avar.gamma);

	/* Simulation settings */
	add_file_aux_data(fileHandle, "SPH integration solver", svar.solver_name);
	add_file_aux_data(fileHandle, "SPH boundary solver", svar.bound_solver);
	add_file_aux_data(fileHandle, "SPH solver minimum residual", svar.minRes);
	add_file_aux_data(fileHandle, "SPH maximum timestep", svar.dt_max);
	add_file_aux_data(fileHandle, "SPH minimum timestep", svar.dt_min);
	add_file_aux_data(fileHandle, "SPH maximum CFL", svar.cfl_max);
	add_file_aux_data(fileHandle, "SPH minimum CFL", svar.cfl_min);
	add_file_aux_data(fileHandle, "SPH CFL condition", svar.cfl);
	add_file_aux_data(fileHandle, "SPH unstable CFL step", svar.cfl_step);
	add_file_aux_data(fileHandle, "SPH unstable CFL count limit", svar.nUnstable_Limit);
	add_file_aux_data(fileHandle, "SPH stable CFL count limit", svar.nStable_Limit);
	add_file_aux_data(fileHandle, "SPH stable CFL count iteration factor", svar.subits_factor);
	add_file_aux_data(fileHandle, "SPH maximum shifting velocity", svar.maxshift);

	add_file_aux_data(fileHandle, "SPH background pressure", fvar.backP);
	add_file_aux_data(fileHandle, "SPH starting pressure", fvar.pPress);
	add_file_aux_data(fileHandle, "SPH density variation", fvar.rhoVar);
	add_file_aux_data(fileHandle, "SPH maximum density", fvar.rhoMax);
	add_file_aux_data(fileHandle, "SPH minimum density", fvar.rhoMin);
	add_file_aux_data(fileHandle, "SPH delta coefficient", fvar.delta);

	add_file_aux_data(fileHandle, "SPH artificial viscosity factor", fvar.alpha);
	add_file_aux_data(fileHandle, "SPH speed of sound", fvar.Cs);
	add_file_aux_data(fileHandle, "SPH Newmark Beta iteration limit", svar.subits);
	add_file_aux_data(fileHandle, "SPH gravity vector", svar.grav);

	add_file_aux_data(fileHandle, "SPH initial spacing", svar.Pstep);
	add_file_aux_data(fileHandle, "SPH boundary spacing factor", svar.Bstep);
	add_file_aux_data(fileHandle, "SPH smoothing length factor", fvar.Hfac);
	add_file_aux_data(fileHandle, "SPH aerodynamic case", avar.aero_case);
	add_file_aux_data(fileHandle, "SPH SP diameter definition", avar.use_dx);
	add_file_aux_data(fileHandle, "SPH use TAB deformation", avar.useDef);
	add_file_aux_data(fileHandle, "SPH use ghost particles", svar.ghost);
	add_file_aux_data(fileHandle, "SPH global offset coordinate", svar.offset_vec);
	add_file_aux_data(fileHandle, "SPH maximum particle count",svar.finPts);
	add_file_aux_data(fileHandle, "SPH aerodynamic cutoff value", avar.cutoff);
	add_file_aux_data(fileHandle, "SPH freestream velocity", avar.vInf);
	add_file_aux_data(fileHandle, "SPH restart fit tolerance", svar.restart_tol);

	/* Particle tracking settings */
	add_file_aux_data(fileHandle, "Transition to IPT", svar.using_ipt);
	add_file_aux_data(fileHandle, "Velocity equation order", svar.eqOrder);
	add_file_aux_data(fileHandle, "SPH tracking conversion x coordinate", svar.max_x_sph);
	add_file_aux_data(fileHandle, "Maximum x trajectory coordinate", svar.max_x);
	add_file_aux_data(fileHandle, "Particle scatter output", svar.partout);
	add_file_aux_data(fileHandle, "Particle streak output", svar.streakout);
	add_file_aux_data(fileHandle, "Particle cell intersection output", svar.cellsout);

	/* Droplet drag sweep settings */
	add_file_aux_data(fileHandle, "Do droplet drag sweep", svar.dropDragSweep);
	add_file_aux_data(fileHandle, "Do speed test", svar.speedTest);
	add_file_aux_data(fileHandle, "Speed test run count", svar.nRuns);

	string str;
	if(!svar.nacross.empty())
	{
		for(auto const& x:svar.nacross)
			str.append(std::to_string(x) + ",");
		str.pop_back(); //Remove the last comma
		add_file_aux_data(fileHandle, "Droplet resolutions",str);
		str.clear();
	}

	if(!svar.diameters.empty())
	{
		for(auto const& x:svar.diameters)
			str.append(std::to_string(x) + ",");
		str.pop_back(); //Remove the last comma
		add_file_aux_data(fileHandle, "Droplet diameters", str);
		str.clear();
	}
	
	if(!svar.velocities.empty())
	{
		for(auto const& x:svar.velocities)
			str.append(std::to_string(x) + ",");
		str.pop_back(); //Remove the last comma
		add_file_aux_data(fileHandle, "Droplet velocities", str);
		str.clear();
	}

	
	if(!svar.Reynolds.empty())
	{
		for(auto const& x:svar.Reynolds)
			str.append(std::to_string(x) + ",");
		str.pop_back(); //Remove the last comma
		add_file_aux_data(fileHandle, "Droplet Reynolds numbers", str);
		str.clear();
	}
	
}
