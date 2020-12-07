/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include <tecio/TECIO.h>
#include "Var.h"
#include "Neighbours.h"
#include "Kernel.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
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
	INTEGER4 I;
	vector<real> xVar(iMax);
	vector<real> yVar(iMax);

    // cout << "Trying to get vector x-component. Var: " << varCount << endl;
	I = Read_Real_Value(inputHandle, frame, varCount, iMax, xVar);
	++varCount;
    // cout << "Trying to get vector y-component. Var: " << varCount << endl;
	I = Read_Real_Value(inputHandle, frame, varCount, iMax, yVar);
	++varCount;

#if SIMDIM == 3
	vector<real> zVar(iMax);
	// cout << "Trying to get vector z-component. Var: " << varCount << endl;
	I = Read_Real_Value(inputHandle, frame, varCount, iMax, zVar);
	++varCount;
#endif

	vector<StateVecD> vec(iMax);
	#pragma omp parallel for
	for(uint ii = 0; ii < xVar.size(); ++ii)
	{
		vec[ii](0) = xVar[ii];
		vec[ii](1) = yVar[ii];
		
#if SIMDIM == 3
		vec[ii](2) = zVar[ii];
#endif

		// cout << tsData.verts[ii](0) << "  " << tsData.verts[ii](1) << endl;
	}

	if(I == -1)
	{
		cout << "Errors occured obtaining vector. Variable count: " << varCount << endl;
	}

	return vec;
}


void Read_Binary_Timestep(void* inputHandle, SIM& svar, INTEGER4 frame, State& pn)
{
// variables = "x y z rho P m v_x v_y v_z a_x a_y a_z Cell_ID Cell_Vx Cell_Vy Cell_Vz Cell_P Cell_Rho";
	
	cout << "Reading zone number: " << frame << endl;

	INTEGER4 I;
	INTEGER8 iMax, jMax, kMax;
	I = tecZoneGetIJK(inputHandle, frame, &iMax, &jMax, &kMax);
	if(I == -1)
	{
		cout << "Reading zone data caused a tecplot error. Stopping." << endl;
		exit(-1);
	}

	cout << "Number of particles in zone: " << iMax << endl;



	I = tecZoneGetSolutionTime(inputHandle, frame, &svar.t);
	vector<StateVecD> xi(iMax);
	vector<real> rho(iMax);
	vector<real> Rrho(iMax);
	vector<real> m(iMax);
	vector<StateVecD> vel(iMax);
	vector<StateVecD> acc(iMax);
	vector<uint8_t> b(iMax);
	vector<StateVecD> cellV;
	vector<real> cellP;
	// vector<real> cellRho;
	vector<int32_t> cellID;
	INTEGER4 varCount = 1;


	xi = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	/*Get density and mass*/
	// cout << "Trying to get variable: " << varCount << endl;
	I = Read_Real_Value(inputHandle, frame, varCount, iMax, rho);
	++varCount;
	// cout << "Trying to get variable: " << varCount << endl;
	I = Read_Real_Value(inputHandle, frame, varCount, iMax, Rrho);
	++varCount;  
	// cout << "Trying to get variable: " << varCount << endl;
	I = Read_Real_Value(inputHandle, frame, varCount, iMax, m);
	++varCount;  

	vel = Read_Binary_Vector(inputHandle, frame, varCount, iMax); 

	acc = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	I = tecZoneVarGetUInt8Values(inputHandle, frame, varCount, 1, iMax, &b[0]);
	++varCount; 

	if(svar.outform == 3)
	{	/*Get the cell information for the points*/
		cellV = vector<StateVecD>(iMax);
		cellP = vector<real>(iMax);
		cellID = vector<int32_t>(iMax);

		cellV = Read_Binary_Vector(inputHandle, frame, varCount, iMax);
		
		/*Get pressure*/
		// cout << "Trying to get variable: " << varCount << endl;
		I = Read_Real_Value(inputHandle, frame, varCount, iMax, cellP);
		++varCount;
		// cout << "Trying to get variable: " << varCount << endl;
		I = tecZoneVarGetInt32Values(inputHandle, frame, varCount, 1, iMax, &cellID[0]);
		++varCount;
	}
	else if (svar.outform == 5)
	{
		// cout << "Trying to get variable: " << varCount << endl;
		cellID = vector<int32_t>(iMax);
		I = tecZoneVarGetInt32Values(inputHandle, frame, varCount, 1, iMax, &cellID[0]);
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
		pn[ii].b = static_cast<uint>(std::round(b[ii]));

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
int Write_Real_Value(void* const& fileHandle, int32_t& outputZone, int32_t const& var, int64_t const& size,
						vector<real>& varVec)
{
#if FOD == 1
	return tecZoneVarWriteDoubleValues(fileHandle, outputZone, var, 0, size, &varVec[0]);
#else
	return tecZoneVarWriteFloatValues(fileHandle, outputZone, var, 0, size, &varVec[0]);
#endif	
}


void Write_Binary_Timestep(SIM const& svar, State const& pnp1, 
	uint const start, uint const end, char const* group, int32_t const& strandID, void* const& fileHandle)
{
	int64_t const size = end - start;

	double solTime = svar.t;     
	int32_t outputZone;
	int I = 0;

	// Get zone data types
	vector<int32_t> varTypes = svar.varTypes;
	vector<int32_t> shareVarFromZone(varTypes.size(),0);
    vector<int32_t> valueLocation(varTypes.size(),1);
    vector<int32_t> passiveVarList(varTypes.size(),0);

	I = tecZoneCreateIJK(fileHandle,group,size,1,1,&varTypes[0],
		&shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&outputZone);

	if(I == -1)
		exit(-1);

	if(strandID != 0)
		I = tecZoneSetUnsteadyOptions(fileHandle, outputZone, solTime, strandID);

	/*Write the basic position data present for all outputs*/
	vector<real> x(size);
	int32_t var = 1;

	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(uint ii = start; ii < end; ++ii)
			x[ii-start] = pnp1[ii].xi(dim)/svar.scale;

		I = Write_Real_Value(fileHandle, outputZone, var, size, x);
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

			I = Write_Real_Value(fileHandle, outputZone, var, size, rho);
			var++;
			I = Write_Real_Value(fileHandle, outputZone, var, size, Rrho);
			var++;
			I = Write_Real_Value(fileHandle, outputZone, var, size, m);
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

		I = Write_Real_Value(fileHandle, outputZone, var, size, v);
		var++;
		I = Write_Real_Value(fileHandle, outputZone, var, size, a);
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

			I = Write_Real_Value(fileHandle, outputZone, var, size, nNb);
			var++;
			I = tecZoneVarWriteUInt8Values(fileHandle, outputZone, var, 0, size, &aF[0]);
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

  		I = Write_Real_Value(fileHandle, outputZone, var, size, vx);
		var++;
		I = Write_Real_Value(fileHandle, outputZone, var, size, vy);
		var++;
#if SIMDIM == 3
		I = Write_Real_Value(fileHandle, outputZone, var, size, vz);
		var++;
#endif
  		I = Write_Real_Value(fileHandle, outputZone, var, size, ax);
		var++;
		I = Write_Real_Value(fileHandle, outputZone, var, size, ay);
		var++;
#if SIMDIM == 3
		I = Write_Real_Value(fileHandle, outputZone, var, size, az);
		var++;
#endif
		I = tecZoneVarWriteUInt8Values(fileHandle, outputZone, var, 0, size, &b[0]);
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

			I = Write_Real_Value(fileHandle, outputZone, var, size, cVx);
			var++;
			I = Write_Real_Value(fileHandle, outputZone, var, size, cVy);
			var++;

			#if SIMDIM == 3
				vector<real> cVz(size);
				#pragma omp parallel for
				for(uint ii = start; ii < end; ++ii)
					cVz[ii-start] = pnp1[ii].cellV(2);
				I = Write_Real_Value(fileHandle, outputZone, var, size, cVz);
				var++;
			#endif

			I = Write_Real_Value(fileHandle, outputZone, var, size, cP);
			var++;
			I = tecZoneVarWriteInt32Values(fileHandle, outputZone, var, 0, size, &cID[0]);
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

			I = tecZoneVarWriteInt32Values(fileHandle, outputZone, var, 0, size, &cID[0]);
			var++;
		}
	    
	}
	
    // cout << "Flushing results." << endl;
    // INTEGER4 numZonesToRetain = 0;
    I = tecFileWriterFlush(fileHandle,0,NULL);
    if(I == -1)
    {
    	exit(-1);
    }
}

void Init_Binary_PLT(SIM &svar, string const& filename, string const& zoneName, void* &fileHandle)
{
    int32_t FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
    int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT
    int I          = 0;     /* Used to track return codes */

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

	I = tecFileWriterOpen(file.c_str(),zoneName.c_str(),variables.c_str(),fileFormat,FileType,1,NULL,&fileHandle);

#ifdef DEBUG
	I = tecFileSetDiagnosticsLevel(fileHandle, 1);
#endif

    if(I == -1)
    {
    	cout << "Failed to open " << file << endl;
    	exit(-1);
    }
}

#endif