/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include <TECIO.h>
#include "Var.h"
#include "Neighbours.h"
#include "Kernel.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <chrono>
#include <thread>
using std::string;

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

/*************************************************************************/
/**************************** BINARY INPUTS ******************************/
/*************************************************************************/
int Read_Real_Value(void* const& inputHandle, int32_t const& frame, int32_t const& varCount, 
					int64_t const& iMax, vector<real>& var)
{
#if FOD == 1
	return tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &var[0]);
#else
	return tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &var[0]);
#endif
}

vector<StateVecD> Read_Binary_Vector(void* inputHandle, INTEGER4& frame, 
		INTEGER4& varCount, INTEGER8& iMax)
{
	vector<real> xVar(iMax,0.0);
	vector<real> yVar(iMax,0.0);

    // cout << "Trying to get vector x-component. Var: " << varCount << endl;
	if(Read_Real_Value(inputHandle, frame, varCount, iMax, xVar))
	{
		cerr << "Failed to read x-component of vector. frame: " << frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	++varCount;
    // cout << "Trying to get vector y-component. Var: " << varCount << endl;
	if(Read_Real_Value(inputHandle, frame, varCount, iMax, yVar))
	{
		cerr << "Failed to read y-component of vector. frame: " << frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	++varCount;

#if SIMDIM == 3
	vector<real> zVar(iMax,0.0);
	// cout << "Trying to get vector z-component. Var: " << varCount << endl;
	if(Read_Real_Value(inputHandle, frame, varCount, iMax, zVar))
	{
		cerr << "Failed to read z-component of vector. frame: " << frame << ". varCount: " << varCount << endl;
		exit(-1);
	}
	++varCount;
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


void Read_Binary_Timestep(void* inputHandle, SIM& svar, int32_t frame, State& pn)
{
// variables = "x y z rho P m v_x v_y v_z a_x a_y a_z Cell_ID Cell_Vx Cell_Vy Cell_Vz Cell_P Cell_Rho";
	
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

	vector<StateVecD> xi(iMax);
	vector<real> rho(iMax,0.0);
	vector<real> Rrho(iMax,0.0);
	vector<real> m(iMax,0.0);
	vector<StateVecD> vel(iMax);
	vector<StateVecD> acc(iMax);
	vector<uint8_t> b(iMax,0);
	vector<StateVecD> cellV;
	vector<real> cellP;
	vector<int32_t> cellID;
	
	int32_t varCount = 1;

	xi = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	/*Get density and mass*/
	// cout << "Trying to get variable: " << varCount << endl;
	if(Read_Real_Value(inputHandle, frame, varCount, iMax, rho))
	{
		cout << "Failed to read density. frame: " << frame << endl;
		exit(-1);
	}	
	++varCount;
	// cout << "Trying to get variable: " << varCount << endl;
	if(Read_Real_Value(inputHandle, frame, varCount, iMax, Rrho))
	{
		cout << "Failed to read Rrho. frame: " << frame << endl;
		exit(-1);
	}
	++varCount;  
	// cout << "Trying to get variable: " << varCount << endl;
	if(Read_Real_Value(inputHandle, frame, varCount, iMax, m))
	{
		cout << "Failed to read mass. frame: " << frame << endl;
		exit(-1);
	}
	++varCount;  

	vel = Read_Binary_Vector(inputHandle, frame, varCount, iMax); 


	acc = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	if(tecZoneVarGetUInt8Values(inputHandle, frame, varCount, 1, iMax, &b[0]))
	{
		cout << "Failed to read particle type. frame: " << frame << endl;
		exit(-1);
	}
	++varCount; 

	if(svar.outform == 3)
	{	/*Get the cell information for the points*/
		cellV = vector<StateVecD>(iMax);
		cellP = vector<real>(iMax,0.0);
		cellID = vector<int32_t>(iMax,0);

		cellV = Read_Binary_Vector(inputHandle, frame, varCount, iMax);
		
		/*Get pressure*/
		// cout << "Trying to get variable: " << varCount << endl;
		if(Read_Real_Value(inputHandle, frame, varCount, iMax, cellP))
		{
			cerr << "Failed to read cell pressure. frame: " << frame << endl;
			exit(-1);
		}
		++varCount;
		// cout << "Trying to get variable: " << varCount << endl;
		if(tecZoneVarGetInt32Values(inputHandle, frame, varCount, 1, iMax, &cellID[0]))
		{
			cerr << "Failed to read cell ID. frame: " << frame << endl;
			exit(-1);
		}
		++varCount;
	}
	else if (svar.outform == 5)
	{
		// cout << "Trying to get variable: " << varCount << endl;
		cellID = vector<int32_t>(iMax,0);
		if(tecZoneVarGetInt32Values(inputHandle, frame, varCount, 1, iMax, &cellID[0]))
		{
			cerr << "Failed to read cell ID. frame: " << frame << endl;
			exit(-1);
		}
		++varCount;
	}
	
	pn = vector<Particle>(iMax);
	/*Now put it into the state vector*/
	for(size_t ii= 0; ii < static_cast<size_t>(iMax); ++ii)
	{
		pn[ii].xi = xi[ii];
		pn[ii].rho = rho[ii];
		pn[ii].Rrho = Rrho[ii];
		pn[ii].m = m[ii];
		pn[ii].v = vel[ii];
		pn[ii].f = acc[ii];
		pn[ii].b = b[ii]; // static_cast<uint>(b[ii]);

		if(svar.outform == 5 || svar.outform == 3)
		{
			pn[ii].cellID = static_cast<size_t>(cellID[ii]);
			if(svar.outform == 3)
			{
				pn[ii].cellV = cellV[ii];
				pn[ii].cellP = cellP[ii];
			}
		}		
	} 
}


/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/
int Write_Real_Vector(void* const& fileHandle, int32_t& outputZone, int32_t const& var, int64_t const& size,
						vector<real> const& varVec)
{
#if FOD == 1
	return tecZoneVarWriteDoubleValues(fileHandle, outputZone, var, 0, size, &varVec[0]);
#else
	return tecZoneVarWriteFloatValues(fileHandle, outputZone, var, 0, size, &varVec[0]);
#endif	
}

int Write_Real_Value(void* const& fileHandle, int32_t& outputZone, int32_t const& var, int64_t const& size,
						real const& value)
{
	vector<real> rvec(size,0.0);
	rvec[0] = value;
#if FOD == 1
	return tecZoneVarWriteDoubleValues(fileHandle, outputZone, var, 0, size, &rvec[0]);
#else
	return tecZoneVarWriteFloatValues(fileHandle, outputZone, var, 0, size, &rvec[0]);
#endif	
}

void Write_Binary_Timestep(SIM const& svar, State const& pnp1, 
	uint const start, uint const end, char const* group, int32_t const& strandID, void* const& fileHandle)
{
	int64_t const size = end - start;

	double solTime = svar.t;     
	int32_t outputZone;

	// Get zone data types
	vector<int32_t> varTypes = svar.varTypes;
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

	/*Write the basic position data present for all outputs*/
	vector<real> x(size);
	int32_t var = 1;

	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(uint ii = start; ii < end; ++ii)
			x[ii-start] = pnp1[ii].xi(dim)/svar.scale;

		if(Write_Real_Vector(fileHandle, outputZone, var, size, x))
		{
			cerr << "Failed to write position coordinate " << dim << endl;
			exit(-1);
		}
		var++;
	}

	if(svar.outform != 0)
	{
		vector<real> rho(size);
		vector<real> Rrho(size);
		vector<real> m(size);

		#pragma omp parallel for
	  	for(uint ii = start; ii < end; ++ii)
		{
			rho[ii-start] = pnp1[ii].rho;
			Rrho[ii-start] = pnp1[ii].Rrho;

			// cout << Rrho[ii-start] << "  ";
			m[ii-start] = pnp1[ii].m;
		}

			if(Write_Real_Vector(fileHandle, outputZone, var, size, rho))
			{
				cerr << "Failed to write density" << endl;
				exit(-1);
			}
			var++;
			if(Write_Real_Vector(fileHandle, outputZone, var, size, Rrho))
			{
				cerr << "Failed to write Rrho" << endl;
				exit(-1);
			}
			var++;
			if(Write_Real_Vector(fileHandle, outputZone, var, size, m))
			{
				cerr << "Failed to write mass" << endl;
				exit(-1);
			}
			var++;
	}

	if(svar.outform == 1 || svar.outform == 4)
	{
		vector<real> v(size);
		vector<real> a(size);
		
		#pragma omp parallel for
		for(uint ii = start; ii < end; ++ii)
		{
			v[ii-start] = pnp1[ii].v.norm();
			a[ii-start] = pnp1[ii].f.norm();
		}

		if(Write_Real_Vector(fileHandle, outputZone, var, size, v))
		{
			cerr << "Failed to write velocity magnitude" << endl;
			exit(-1);
		}
		var++;
		if(Write_Real_Vector(fileHandle, outputZone, var, size, a))
		{
			cerr << "Failed to write acceleration magnitude" << endl;
			exit(-1);
		}
		var++;

		if(svar.outform == 4)
		{
			vector<real> nNb(size);
			vector<uint8_t> aF(size);

			#pragma omp parallel for
		  	for(uint ii = start; ii < end; ++ii)
		  	{
		  		nNb[ii-start] = pnp1[ii].curve;
		  		aF[ii-start] = pnp1[ii].surf;
	  		}

			if(Write_Real_Vector(fileHandle, outputZone, var, size, nNb))
			{
				cerr << "Failed to write real value" << endl;
				exit(-1);
			}
			var++;
			if(tecZoneVarWriteUInt8Values(fileHandle, outputZone, var, 0, size, &aF[0]))
			{
				cerr << "Failed to write number of Neighbours" << endl;
				exit(-1);
			}
			var++;
		}
	}
	else if (svar.outform == 2 || svar.outform == 3 || svar.outform == 5)
	{

		vector<real> vx(size);
		vector<real> vy(size);
		vector<real> ax(size);
		vector<real> ay(size);
#if SIMDIM == 3
		vector<real> vz(size);
		vector<real> az(size);
#endif
		vector<uint8_t> b(size);

		#pragma omp parallel for
		for(uint ii = start; ii < end; ++ii)
		{
			vx[ii-start] = pnp1[ii].v(0);
			vy[ii-start] = pnp1[ii].v(1);
			ax[ii-start] = pnp1[ii].f(0);
			ay[ii-start] = pnp1[ii].f(1);
			
			#if SIMDIM == 3
				vz[ii-start] = pnp1[ii].v(2);
				az[ii-start] = pnp1[ii].f(2);
			#endif
  			b[ii-start] = static_cast<uint8_t>(pnp1[ii].b);
  	}

  		if(Write_Real_Vector(fileHandle, outputZone, var, size, vx))
  		{
			cerr << "Failed to write velocity x-component" << endl;
			exit(-1);
		}
		var++;

		if(Write_Real_Vector(fileHandle, outputZone, var, size, vy))
  		{
			cerr << "Failed to write velocity y-component" << endl;
			exit(-1);
		}
		var++;

		#if SIMDIM == 3
			if(Write_Real_Vector(fileHandle, outputZone, var, size, vz))
			{
				cerr << "Failed to write velocity z-component" << endl;
				exit(-1);
			}
			var++;
		#endif

  		if(Write_Real_Vector(fileHandle, outputZone, var, size, ax))
  		{
			cerr << "Failed to write acceleration x-component" << endl;
			exit(-1);
		}
		var++;

		if(Write_Real_Vector(fileHandle, outputZone, var, size, ay))
  		{
			cerr << "Failed to write acceleration y-component" << endl;
			exit(-1);
		}
		var++;

		#if SIMDIM == 3
			if(Write_Real_Vector(fileHandle, outputZone, var, size, az))
			{
				cerr << "Failed to write acceleration z-component" << endl;
				exit(-1);
			}
			var++;
		#endif
		if(tecZoneVarWriteUInt8Values(fileHandle, outputZone, var, 0, size, &b[0]))
		{
			cerr << "Failed to write particle type" << endl;
			exit(-1);
		}
		var++;

		if (svar.outform == 3)
		{
			vector<int32_t> cID(size);
			vector<real> cVx(size);
			vector<real> cVy(size);
			vector<real> cP(size);

			#pragma omp parallel for
			for(uint ii = start; ii < end; ++ii)
			{
				cID[ii-start] = static_cast<int32_t>(pnp1[ii].cellID);
				cVx[ii-start] = pnp1[ii].cellV(0);
				cVy[ii-start] = pnp1[ii].cellV(1);
				cP[ii-start] = pnp1[ii].cellP;
			}

			if(Write_Real_Vector(fileHandle, outputZone, var, size, cVx))
			{
				cerr << "Failed to write cell velocity x-component" << endl;
				exit(-1);
			}
			var++;
			if(Write_Real_Vector(fileHandle, outputZone, var, size, cVy))
			{
				cerr << "Failed to write cell velocity y-component" << endl;
				exit(-1);
			}
			var++;

			#if SIMDIM == 3
				vector<real> cVz(size);
				#pragma omp parallel for
				for(uint ii = start; ii < end; ++ii)
					cVz[ii-start] = pnp1[ii].cellV(2);
				if(Write_Real_Vector(fileHandle, outputZone, var, size, cVz))
				{
					cerr << "Failed to write cell velocity z-component" << endl;
					exit(-1);
				}
				var++;
			#endif

			if(Write_Real_Vector(fileHandle, outputZone, var, size, cP))
			{
				cerr << "Failed to write cell pressure" << endl;
				exit(-1);
			}
			var++;
			if(tecZoneVarWriteInt32Values(fileHandle, outputZone, var, 0, size, &cID[0]))
			{
				cerr << "Failed to write cell ID" << endl;
				exit(-1);
			}
			var++;
		}
		else if (svar.outform == 5)
		{
			vector<int32_t> cID(size);
			#pragma omp parallel for
			for(uint ii = start; ii < end; ++ii)
			{
				cID[ii-start] = static_cast<int32_t>(pnp1[ii].cellID);
			}

			if(tecZoneVarWriteInt32Values(fileHandle, outputZone, var, 0, size, &cID[0]))
			{
				cerr << "Failed to write cell ID" << endl;
				exit(-1);
			}
			var++;
		}
	    
	}
	
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

    string file = svar.outfolder;
    file.append(filename);

    // cout << file << endl;

    /* Open the file and write the tecplot datafile header information  */
#if SIMDIM == 2
    	std::string variables = "X,Z";	
		if (svar.outform == 1)
		{
			variables = "X,Z,rho,Rrho,m,v,a";
		}
		else if (svar.outform == 2)
		{
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b";
		}
		else if (svar.outform == 3)
		{
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_Vx,Cell_Vz,Cell_P,Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "X,Z,rho,Rrho,m,v,a,Neighbours,Aero";
		}
		else if (svar.outform == 5)
		{
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_ID";
		}
#endif

#if SIMDIM == 3
		std::string variables = "X,Y,Z";  
		if (svar.outform == 1)
		{
			variables = "X,Y,Z,rho,Rrho,m,v,a";
		}
		else if (svar.outform == 2)
		{
			variables = "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b";
		}
		else if (svar.outform == 3)
		{
			variables = 
		"X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_Vx,Cell_Vy,Cell_Vz,Cell_P,Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "X,Y,Z,rho,Rrho,m,v,a,Neighbours,Aero";
		}
		else if (svar.outform == 5)
		{
			variables = "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_ID";
		}
#endif

	vector<int32_t> varTypes;
#if FOD == 1
	int32_t realType = 2;
#else
	int32_t realType = 1;
#endif

	varTypes.emplace_back(realType);
	varTypes.emplace_back(realType);
#if SIMDIM == 3
	varTypes.emplace_back(realType);	//x,y,z
#endif
	if(svar.outform != 0)
	{
		varTypes.emplace_back(realType);  //rho
		varTypes.emplace_back(realType);  //Rrho
		varTypes.emplace_back(realType);  //m

		if(svar.outform == 1 || svar.outform == 4)
		{
			varTypes.emplace_back(realType);  //vmag
			varTypes.emplace_back(realType);  //amag

			if(svar.outform == 4)
			{
				varTypes.emplace_back(realType);  //Nneighb
				varTypes.emplace_back(realType);  //vy
			}
		}
		else if (svar.outform == 2 || svar.outform == 3 || svar.outform == 5)
		{
			varTypes.emplace_back(realType);  //vx
			varTypes.emplace_back(realType);  //vy
			varTypes.emplace_back(realType);  //ax
			varTypes.emplace_back(realType);  //ay
#if SIMDIM == 3
			varTypes.emplace_back(realType);  //vz
			varTypes.emplace_back(realType);  //az
#endif
			varTypes.emplace_back(5);  //b

			if(svar.outform == 3)
			{
				varTypes.emplace_back(realType);  //cell_Vx
				varTypes.emplace_back(realType);  //cell_Vy
				#if SIMDIM == 3
				varTypes.emplace_back(realType);  //cell_Vz
				#endif
				varTypes.emplace_back(realType);  //cell_P
				varTypes.emplace_back(3);  //cell_ID
			}

			if(svar.outform == 5)
			{
				varTypes.emplace_back(3);  //cell_ID
			}

		}
	}
	svar.varTypes = varTypes;

	if(tecFileWriterOpen(file.c_str(),zoneName.c_str(),variables.c_str(),fileFormat,FileType,1,NULL,&fileHandle))
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

    
}



#endif