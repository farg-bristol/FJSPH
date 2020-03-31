/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include <tecio/TECIO.h>
#include "Var.h"
#include "Neighbours.h"
#include "Kernel.h"
#include "IO.h"
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

vector<StateVecD> Read_Binary_Vector(void* inputHandle, INTEGER4& frame, 
		INTEGER4& varCount, INTEGER8& iMax)
{
	INTEGER4 I;
	vector<real> xVar(iMax);
	vector<real> yVar(iMax);
#if FOD == 1
    // cout << "Trying to get vector x-component. Var: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &xVar[0]);
	++varCount;
    // cout << "Trying to get vector y-component. Var: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &yVar[0]);
	++varCount;

#if SIMDIM == 3
	vector<real> zVar(iMax);
	// cout << "Trying to get vector z-component. Var: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &zVar[0]);
	++varCount;
#endif
#else
	I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &xVar[0]);
	++varCount;
	I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &yVar[0]);
	++varCount;

#if SIMDIM == 3
	vector<real> zVar(iMax);
	I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &zVar[0]);
	++varCount;
#endif
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
	vector<real> b(iMax);
	vector<StateVecD> cellV;
	vector<real> cellP;
	vector<real> cellRho;
	vector<real> cellID;
	INTEGER4 varCount = 1;

#if FOD == 1
	xi = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	/*Get density and mass*/
	// cout << "Trying to get variable: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &rho[0]);
	++varCount;
	// cout << "Trying to get variable: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &Rrho[0]);
	++varCount;  
	// cout << "Trying to get variable: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &m[0]);
	++varCount;  

	vel = Read_Binary_Vector(inputHandle, frame, varCount, iMax); 

	acc = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &b[0]);
	++varCount; 

	if(svar.outform == 3)
	{	/*Get the cell information for the points*/
		
		cellV = Read_Binary_Vector(inputHandle, frame, varCount, iMax);
		cellRho = vector<real>(iMax);
		cellP = vector<real>(iMax);
		cellID = vector<real>(iMax);
		/*Get density and pressure*/
		// cout << "Trying to get variable: " << varCount << endl;
		I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &cellRho[0]);
		++varCount;
		// cout << "Trying to get variable: " << varCount << endl;
		I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &cellP[0]);
		++varCount;
		// cout << "Trying to get variable: " << varCount << endl;
		I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &cellID[0]);
		++varCount;
	}
	else if (svar.outform == 5)
	{
		// cout << "Trying to get variable: " << varCount << endl;
		cellID = vector<real>(iMax);
		I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &cellID[0]);
		++varCount;
	}
	
#else
	
	xi = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	/*Get density and mass*/
	I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &rho[0]);
	++varCount;
	I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &Rrho[0]);
	++varCount;
	I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &m[0]);
	++varCount;  

	vel = Read_Binary_Vector(inputHandle, frame, varCount, iMax); 

	acc = Read_Binary_Vector(inputHandle, frame, varCount, iMax);

	if(svar.outform == 3)
	{	/*Get the cell information for the points*/
		
		cellV = Read_Binary_Vector(inputHandle, frame, varCount, iMax);
		cellRho = vector<real>(iMax);
		cellP = vector<real>(iMax);
		/*Get density and pressure*/
		I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &cellRho[0]);
		++varCount;
		I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &cellP[0]);
		++varCount;
		I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &cellID[0]);
		++varCount;
	}
	else if (svar.outform == 5)
	{
		I = tecZoneVarGetFloatValues(inputHandle, frame, varCount, 1, iMax, &cellID[0]);
		++varCount;
	}
#endif

	pn = vector<Particle>(iMax);
	/*Now put it into the state vector*/
	for(size_t ii= 0; ii < static_cast<size_t>(iMax); ++ii)
	{
		pn[ii].xi = xi[ii];
		pn[ii].v = vel[ii];
		pn[ii].f = acc[ii];
		pn[ii].xi = xi[ii];
		pn[ii].Rrho = Rrho[ii];
		pn[ii].rho = rho[ii];
		pn[ii].m = m[ii];
		pn[ii].b = static_cast<uint>(b[ii]);
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

void Write_Binary_Timestep(const SIM& svar, const State& pnp1, 
	const uint start, const uint end, const char* group, const INTEGER4 StrandID)
{
	if(StrandID == 0)
	{
		INTEGER4 File = 2;
		TECFIL142(&File);
	}
	else
		TECFIL142(&StrandID);

	uint size = end - start;

	#ifdef DEBUG
		
	#endif

	const INTEGER4 IMax = size;
	const INTEGER4 JMax = 1;
	const INTEGER4 KMax = 1;
	const INTEGER4 ICellMax                 = 0;
	const INTEGER4 JCellMax                 = 0;
	const INTEGER4 KCellMax                 = 0;
#if FOD == 0
	const INTEGER4 DIsDouble                = 0;
#else
	const INTEGER4 DIsDouble                = 1;
#endif
	const double   SolTime                  = svar.t;     
	const INTEGER4 ParentZn                 = 0;
	const INTEGER4 IsBlock                  = 1;      /* Block */
	const INTEGER4 NFConns                  = 0;
	const INTEGER4 FNMode                   = 0;
	const INTEGER4 TotalNumFaceNodes        = 0;
	const INTEGER4 TotalNumBndryFaces       = 0;
	const INTEGER4 TotalNumBndryConnections = 0;
	const INTEGER4 ShrConn                  = 0;
	const INTEGER4 ZoneType = 0;
	INTEGER4 I = 0;


	I = TECZNE142(group,
	              &ZoneType,
	              &IMax,
	              &JMax,
	              &KMax,
	              &ICellMax,
	              &JCellMax,
	              &KCellMax,
	              &SolTime,
	              &StrandID,
	              &ParentZn,
	              &IsBlock,
	              &NFConns,
	              &FNMode,
	              &TotalNumFaceNodes,
	              &TotalNumBndryFaces,
	              &TotalNumBndryConnections,
	              NULL,
	              NULL,
	              NULL,
	              &ShrConn);
	if(I == -1)
		exit(-1);

	/*Write the basic position data present for all outputs*/
	real* x = new real[size];

	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(uint ii = start; ii < end; ++ii)
			x[ii-start] = pnp1[ii].xi(dim)/svar.scale;

		I   = TECDAT142(&IMax, x, &DIsDouble);
	}

	if(svar.outform != 0)
	{
		real* rho = new real[size];
		real* Rrho = new real[size];
		real* m = new real[size];

		#pragma omp parallel for
	  	for(uint ii = start; ii < end; ++ii)
		{
			rho[ii-start] = pnp1[ii].rho;
			Rrho[ii-start] = pnp1[ii].Rrho;
			m[ii-start] = pnp1[ii].m;
		}

			I   = TECDAT142(&IMax, rho, &DIsDouble);
			I   = TECDAT142(&IMax, Rrho, &DIsDouble);
			I   = TECDAT142(&IMax, m, &DIsDouble);
	}

	if(svar.outform == 1)
	{
		real* v = new real[size];
		real* a = new real[size];
		
		#pragma omp parallel for
		for(uint ii = start; ii < end; ++ii)
		{
			v[ii-start] = pnp1[ii].v.norm();
			a[ii-start] = pnp1[ii].f.norm();
		}

		I   = TECDAT142(&IMax, v, &DIsDouble);
	  I   = TECDAT142(&IMax, a, &DIsDouble);
	}
	else if (svar.outform == 2 || svar.outform == 3 || svar.outform == 5)
	{

		real* vx = new real[size];
		real* vy = new real[size];
		real* ax = new real[size];
		real* ay = new real[size];
#if SIMDIM == 3
		real* az = new real[size];
		real* vz = new real[size];
#endif
		real* b = new real[size];

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
  		b[ii-start] = static_cast<real>(pnp1[ii].b);
  	}

		I   = TECDAT142(&IMax, vx, &DIsDouble);
		I   = TECDAT142(&IMax, vy, &DIsDouble);
#if SIMDIM == 3
		I   = TECDAT142(&IMax, vz, &DIsDouble);
#endif
		I   = TECDAT142(&IMax, ax, &DIsDouble);
		I   = TECDAT142(&IMax, ay, &DIsDouble);
#if SIMDIM == 3
		I   = TECDAT142(&IMax, az, &DIsDouble);
#endif
		I   = TECDAT142(&IMax, b, &DIsDouble);

		if (svar.outform == 3)
		{
			real* cID = new real[size];
			real* cVx = new real[size];
			real* cVy = new real[size];
			real* cP = new real[size];

			#pragma omp parallel for
			for(uint ii = start; ii < end; ++ii)
			{
				cID[ii-start] = static_cast<real>(pnp1[ii].cellID);
				cVx[ii-start] = pnp1[ii].cellV(0);
				cVy[ii-start] = pnp1[ii].cellV(1);
				cP[ii-start] = pnp1[ii].cellP;
			}

			I   = TECDAT142(&IMax, cVx, &DIsDouble);
			I   = TECDAT142(&IMax, cVy, &DIsDouble);

			#if SIMDIM == 3
				real* cVz = new real[size];
				#pragma omp parallel for
				for(uint ii = start; ii < end; ++ii)
					cVz[ii-start] = pnp1[ii].cellV(2);
				I   = TECDAT142(&IMax, cVz, &DIsDouble);
			#endif

			I   = TECDAT142(&IMax, cP, &DIsDouble);
			I   = TECDAT142(&IMax, cID, &DIsDouble);
		}
		else if (svar.outform == 5)
		{
			real* cID = new real[size];
			#pragma omp parallel for
			for(uint ii = start; ii < end; ++ii)
			{
				cID[ii-start] = static_cast<real>(pnp1[ii].cellID);
			}

			I   = TECDAT142(&IMax, cID, &DIsDouble);
		}
	    
	}
	else if (svar.outform == 4)
		{
			// variables = "x y z rho P m v a b Neighbours Aero";

			real* v = new real[size];
			real* a = new real[size];
			float* nNb = new float[size];
			real* aF = new real[size];

			#pragma omp parallel for
		  	for(uint ii = start; ii < end; ++ii)
		  	{
		  		v[ii-start] = pnp1[ii].v.norm();
		  		a[ii-start] = pnp1[ii].f.norm();
		  		nNb[ii-start] = pnp1[ii].theta;
		  		aF[ii-start] = pnp1[ii].vPert.norm();
	  		}

	  		INTEGER4 IsFloat = 0;
			I   = TECDAT142(&IMax, v, &DIsDouble);
			I   = TECDAT142(&IMax, a, &DIsDouble);
			I   = TECDAT142(&IMax, nNb, &IsFloat);
			I   = TECDAT142(&IMax, aF, &DIsDouble);
		}
    // cout << "Flushing results." << endl;
    INTEGER4 numZonesToRetain = 0;
    I = TECFLUSH142(&numZonesToRetain,NULL);
}

void Init_Binary_PLT(const SIM &svar, const string& filename, const string& zoneName)
{
#if DEBUG
		const INTEGER4 Debug = 1;
#else 
		const INTEGER4 Debug = 0;
#endif
#if FOD == 0
    const INTEGER4 VIsDouble  = 0;
#else
	const INTEGER4 VIsDouble  = 1;
#endif
    
    const INTEGER4 FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
    const INTEGER4 fileFormat = 1; // 0 == PLT, 1 == SZPLT
    INTEGER4 I          = 0; /* Used to track return codes */

    string file = svar.outfolder;
    file.append(filename);
    /*
     * Open the file and write the tecplot datafile
     * header information
     */
#if SIMDIM == 2
    	std::string variables = "X Z";	
		if (svar.outform == 1)
		{
			variables = "X Z rho Rrho m v a";
		}
		else if (svar.outform == 2)
		{
			variables = "X Z rho Rrho m v_x v_z a_x a_z b";
		}
		else if (svar.outform == 3)
		{
			variables = "X Z rho Rrho m v_x v_z a_x a_z b Cell_Vx Cell_Vz Cell_P Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "X Z rho Rrho m v a Neighbours Aero";
		}
		else if (svar.outform == 5)
		{
			variables = "X Z rho Rrho m v_x v_z a_x a_z b Cell_ID";
		}
#endif

#if SIMDIM == 3
		std::string variables = "X Y Z";  
		if (svar.outform == 1)
		{
			variables = "X Y Z rho Rrho m v a";
		}
		else if (svar.outform == 2)
		{
			variables = "X Y Z rho Rrho m v_x v_y v_z a_x a_y a_z b";
		}
		else if (svar.outform == 3)
		{
			variables = 
		"X Y Z rho Rrho m v_x v_y v_z a_x a_y a_z b Cell_Vx Cell_Vy Cell_Vz Cell_P Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "X Y Z rho Rrho m v a Neighbours Aero";
		}
		else if (svar.outform == 5)
		{
			variables = "X Y Z rho Rrho m v_x v_y v_z a_x a_y a_z b Cell_ID";
		}
#endif

	I = TECINI142(  zoneName.c_str(), 
					variables.c_str(),  
		            file.c_str(),
		            svar.outfolder.c_str(),      /* Scratch Directory */
		            &fileFormat,
		            &FileType,
		            &Debug,
		            &VIsDouble);

    if(I == -1)
    {
    	cout << "Failed to open Fuel.szplt" << endl;
    	exit(-1);
    }
}

#endif