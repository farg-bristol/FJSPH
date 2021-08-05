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
void Read_Real_Value(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
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

void Read_Int_Value(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
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

void Read_UInt_Value(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
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

void Read_Real_Vector(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
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

	++varCount;
}

void Read_Int_Vector(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
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
	++varCount;
}

void Read_UInt_Vector(void* const& inputHandle, int32_t const& frame, int32_t& varCount, 
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
	++varCount;
}

vector<StateVecD> Read_Binary_Vector(void* inputHandle, INTEGER4& frame, 
		INTEGER4& varCount, INTEGER8& iMax, string const& varName)
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
	vector<uint> b(iMax,0);
	vector<StateVecD> cellV;
	vector<real> cellP;
	vector<real> cellRho;
	vector<int> cellID;
	
	int32_t varCount = 1;

	xi = Read_Binary_Vector(inputHandle, frame, varCount, iMax, "position");

	/*Get density and mass*/
	// cout << "Trying to get variable: " << varCount << endl;
	Read_Real_Vector(inputHandle, frame, varCount, iMax, rho, "density");

	// cout << "Trying to get variable: " << varCount << endl;
	Read_Real_Vector(inputHandle, frame, varCount, iMax, Rrho, "density gradient");

	// cout << "Trying to get variable: " << varCount << endl;
	Read_Real_Vector(inputHandle, frame, varCount, iMax, m, "mass");

	vel = Read_Binary_Vector(inputHandle, frame, varCount, iMax, "velocity"); 

	acc = Read_Binary_Vector(inputHandle, frame, varCount, iMax, "acceleration");

	Read_UInt_Vector(inputHandle, frame, varCount, iMax, b, "particle type");


	if(svar.outform == 3)
	{	/*Get the cell information for the points*/
		cellV = vector<StateVecD>(iMax);
		cellP = vector<real>(iMax,0.0);
		cellID = vector<int>(iMax,0);

		cellV = Read_Binary_Vector(inputHandle, frame, varCount, iMax, "cell velocity");
		
		/*Get pressure*/
		// cout << "Trying to get variable: " << varCount << endl;
		Read_Real_Vector(inputHandle, frame, varCount, iMax, cellP, "cell pressure");

		Read_Real_Vector(inputHandle, frame, varCount, iMax, cellRho, "cell density");

		// cout << "Trying to get variable: " << varCount << endl;
		Read_Int_Vector(inputHandle, frame, varCount, iMax, cellID, "cell ID");
	}
	else if (svar.outform == 5)
	{
		// cout << "Trying to get variable: " << varCount << endl;
		cellID = vector<int>(iMax,0);
		Read_Int_Vector(inputHandle, frame, varCount, iMax, cellID, "cell ID");
	}
	
	pn = vector<Particle>(iMax);
	/*Now put it into the state vector*/
	for(int64_t ii= 0; ii < iMax; ++ii)
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
				pn[ii].cellRho = cellRho[ii];
			}
		}		
	} 
}


/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/
void Write_Real_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int64_t const& size,
						vector<real> const& varVec, string const& varName)
{
	int retval;
#if FOD == 1
	retval =  tecZoneVarWriteDoubleValues(fileHandle, outputZone, varCount, 0, size, &varVec[0]);
#else
	retval =  tecZoneVarWriteFloatValues(fileHandle, outputZone, varCount, 0, size, &varVec[0]);
#endif	

	if(retval)
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_Int_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int64_t const& size,
						vector<int32_t> const& varVec, string const& varName)
{
	if(tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &varVec[0]))
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_UInt_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int64_t const& size,
						vector<uint8_t> const& varVec, string const& varName)
{
	if(tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &varVec[0]))
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_Real_Value(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int64_t const& size,
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
		cout << "Failed to read \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_Int_Value(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int64_t const& size,
						int32_t const& value, string const& varName)
{
	int retval;
	vector<int32_t> rvec(size,0.0);
	rvec[0] = value;

	retval = tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &rvec[0]);

	if(retval)
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_UInt_Value(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int64_t const& size,
						uint8_t const& value, string const& varName)
{
	int retval;
	vector<uint8_t> rvec(size,0.0);
	rvec[0] = value;

	retval = tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &rvec[0]);

	if(retval)
	{
		cout << "Failed to read \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
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

		string name = "position coordinate ";
		name.append(std::to_string(dim));
		Write_Real_Vector(fileHandle, outputZone, var, size, x, name);
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

		Write_Real_Vector(fileHandle, outputZone, var, size, rho, "density");

		Write_Real_Vector(fileHandle, outputZone, var, size, Rrho, "density gradient");

		Write_Real_Vector(fileHandle, outputZone, var, size, m, "mass");

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

		Write_Real_Vector(fileHandle, outputZone, var, size, v, "velocity magnitude");

		Write_Real_Vector(fileHandle, outputZone, var, size, a, "acceleration magnitude");

		if(svar.outform == 4)
		{
			vector<real> nNb(size);
			vector<uint8_t> aF(size);

			#pragma omp parallel for
		  	for(uint ii = start; ii < end; ++ii)
		  	{
		  		nNb[ii-start] = pnp1[ii].s;
		  		aF[ii-start] = pnp1[ii].surf;
	  		}

			Write_Real_Vector(fileHandle, outputZone, var, size, nNb, "real value");
			
			Write_UInt_Vector(fileHandle, outputZone, var, size, aF, "neighbour count");
		}
	}
	else if (svar.outform == 2 || svar.outform == 3 || svar.outform == 5 
			|| svar.outform == 6 || svar.outform == 7)
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

  		Write_Real_Vector(fileHandle, outputZone, var, size, vx, "velocity x-component");

		Write_Real_Vector(fileHandle, outputZone, var, size, vy, "velocity y-component");

		#if SIMDIM == 3
			Write_Real_Vector(fileHandle, outputZone, var, size, vz, "velocity z-component");
		#endif

  		Write_Real_Vector(fileHandle, outputZone, var, size, ax, "acceleration x-component");

		Write_Real_Vector(fileHandle, outputZone, var, size, ay, "acceleration y-component");

		#if SIMDIM == 3
			Write_Real_Vector(fileHandle, outputZone, var, size, az, "acceleration z-component");
		#endif

		Write_UInt_Vector(fileHandle, outputZone, var, size, b, "particle type");


		if (svar.outform == 3)
		{
			vector<int32_t> cID(size);
			vector<real> cVx(size);
			vector<real> cVy(size);
			#if SIMDIM == 3
				vector<real> cVz(size);
			#endif			
			vector<real> cP(size);
			vector<real> cRho(size);

			#pragma omp parallel for
			for(uint ii = start; ii < end; ++ii)
			{
				cID[ii-start] = static_cast<int32_t>(pnp1[ii].cellID);
				cVx[ii-start] = pnp1[ii].cellV(0);
				cVy[ii-start] = pnp1[ii].cellV(1);
				#if SIMDIM == 3
				cVz[ii-start] = pnp1[ii].cellV(2);
				#endif
				cP[ii-start] = pnp1[ii].cellP;
				cRho[ii-start] = pnp1[ii].cellRho;
			}


			Write_Real_Vector(fileHandle, outputZone, var, size, cVx, "cell velocity x-component");

			Write_Real_Vector(fileHandle, outputZone, var, size, cVy, "cell velocity y-component");

			#if SIMDIM == 3
				Write_Real_Vector(fileHandle, outputZone, var, size, cVz, "cell velocity z-component");
			#endif

			Write_Real_Vector(fileHandle, outputZone, var, size, cP, "cell pressure");

			Write_Real_Vector(fileHandle, outputZone, var, size, cRho, "cell density");

			Write_Int_Vector(fileHandle, outputZone, var, size, cID, "cell ID");
		}
		else if (svar.outform == 5)
		{
			vector<int32_t> cID(size);
			#pragma omp parallel for
			for(uint ii = start; ii < end; ++ii)
			{
				cID[ii-start] = static_cast<int32_t>(pnp1[ii].cellID);
			}

			Write_Int_Vector(fileHandle, outputZone, var, size, cID, "cell ID");
		}
		else if(svar.outform == 6 || svar.outform == 7)
		{
			vector<real> nNb(size);
			vector<real> aF(size);

			#pragma omp parallel for
		  	for(uint ii = start; ii < end; ++ii)
		  	{
		  		nNb[ii-start] = pnp1[ii].s;
		  		aF[ii-start] = pnp1[ii].surf;
	  		}

			Write_Real_Vector(fileHandle, outputZone, var, size, nNb, "real value");
			
			Write_Real_Vector(fileHandle, outputZone, var, size, aF, "aero force");

			if(svar.outform == 7)
			{
				vector<real> aFx(size);
				vector<real> aFy(size);
				#if SIMDIM == 3
					vector<real> aFz(size);
				#endif			

				#pragma omp parallel for
				for(uint ii = start; ii < end; ++ii)
				{
					aFx[ii-start] = pnp1[ii].Af(0);
					aFy[ii-start] = pnp1[ii].Af(1);
					#if SIMDIM == 3
					aFz[ii-start] = pnp1[ii].Af(2);
					#endif
				}

				Write_Real_Vector(fileHandle, outputZone, var, size, aFx, "aero force x-component");

				Write_Real_Vector(fileHandle, outputZone, var, size, aFy, "aero force y-component");

				#if SIMDIM == 3
					Write_Real_Vector(fileHandle, outputZone, var, size, aFz, "aero force z-component");
				#endif
			}
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
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_Vx,Cell_Vz,Cell_P,Cell_Rho,Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "X,Z,rho,Rrho,m,v,a,Neighbours,Aero";
		}
		else if (svar.outform == 5)
		{
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_ID";
		}
		else if (svar.outform == 6)
		{
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,lambda,surface";
		}
		else if (svar.outform == 7)
		{
			variables = "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,lambda,surface,a_aero_x,a_aero_z";
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
		"X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_Vx,Cell_Vy,Cell_Vz,Cell_P,Cell_Rho,Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "X,Y,Z,rho,Rrho,m,v,a,Neighbours,Aero";
		}
		else if (svar.outform == 5)
		{
			variables = "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_ID";
		}
		else if (svar.outform == 6)
		{
			variables = "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,lambda,surface";
		}
		else if (svar.outform == 7)
		{
			variables = "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,lambda,surface,a_aero_x,a_aero_z";
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
		else if (svar.outform == 2 || svar.outform == 3 || 
			svar.outform == 5 || svar.outform == 6 || svar.outform == 7)
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
				varTypes.emplace_back(realType);  //cell_Rho
				varTypes.emplace_back(3);  //cell_ID
			}

			if(svar.outform == 5)
			{
				varTypes.emplace_back(3);  //cell_ID
			}

			if(svar.outform == 6 || svar.outform == 7)
			{
				varTypes.emplace_back(realType); /* Lambda */
				varTypes.emplace_back(realType);	/* surface */
				if(svar.outform == 7)
				{
					varTypes.emplace_back(realType);
					varTypes.emplace_back(realType);
					#if SIMDIM == 3
					varTypes.emplace_back(realType);  //cell_Vz
					#endif
				}
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

void Write_Cell_Data(MESH const& cdata)
{
#ifdef DEBUG
	dbout << "Entering Write_Cell_Data..." << endl;
#endif

	cout << "Writing cell based data." << endl;

	std::ofstream fout("Cell.dat", std::ios::out);
	if (!fout.is_open())
	{
		cout << "Failed to open data file for writing mesh." << endl;
		exit(-1);
	}

	fout << "TITLE = \"3D Mesh Solution\"\n";
	fout << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\"\n";
	fout << "ZONE T=\"Cell Data\"" << endl;
	fout << "N=" << cdata.numPoint << ", E=" << cdata.numElem << ", F=FEBLOCK, ET=BRICK" << endl
		 << endl;

	/*Write vertices*/
	fout << std::left << std::scientific << std::setprecision(6);
	fout << std::setw(1);
	for (uint ii = 0; ii < SIMDIM; ++ii)
	{
		uint kk = 0;
		for (uint jj = 0; jj < cdata.verts.size(); ++jj)
		{
			fout << std::setw(15) << cdata.verts[jj][ii];
			kk++;

			if (kk == 5)
			{
				fout << endl;
				fout << std::setw(1);
				kk = 0;
			}
		}

		if (kk % 5 != 0)
			fout << "\n";
	}

	/*Write element indexing*/
	fout << std::fixed;
	for (uint ii = 0; ii < cdata.elems.size(); ++ii)
	{
		for (auto elem : cdata.elems[ii])
		{
			fout << std::setw(6) << elem + 1;
		}
		fout << "\n";
	}

	fout.close();

#ifdef DEBUG
	dbout << "Exiting Write_Cell_Data..." << endl;
#endif
}

void Write_Face_Data(MESH const& cells)
{
	
	std::ofstream f1("foam_mesh.dat", std::ios::out);
	f1 << "VARIABLES= \"X\", \"Y\", \"Z\"" << endl;
	uint w = 17;
	f1 << std::left << std::scientific << std::setprecision(8);

	/*Write zone header information*/
	f1 << "ZONE T=\"OpenFOAM MESH\"" << endl;
	f1 << "ZONETYPE=FEPOLYHEDRON" << endl;
	f1 << "NODES=" << cells.numPoint << " ELEMENTS=" << cells.numElem << " FACES=" << cells.numFace << endl;
	size_t TotalNumFaceNodes = cells.numFace * 3;
	f1 << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
	f1 << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

	
	/*Write vertices in block format*/
	size_t n = 0;

	for (size_t dim = 0; dim < SIMDIM; ++dim)
	{
		for (auto const &pnt : cells.verts)
		{
			f1 << std::setw(w) << pnt[dim];
			n++;

			if(n == 6)
			{
				f1 << endl;
				n = 0;
			}
		}
		f1 << endl;
	}

	/*Write how many vertices per face*/
	n = 0;
	for (size_t ii = 0; ii < cells.numFace; ++ii)
	{
		f1 << std::setw(5) << 3;
		n++;

		if(n == 6)
		{
			f1 << endl;
			n = 0;
		}
	}
	f1 << endl;

	/*Write the face vertex list*/
	n = 0;
	w = 10;
	for (auto const &face : cells.faces)
	{
		for (auto const &vertex : face)
		{ 
			f1 << std::setw(w) << vertex + 1;
			n++;

			if (n == 6)
			{
				f1 << endl;
				n = 0;
			}
		}
	}
	f1 << endl;

	/*Write left elements*/
	n = 0;
	for (auto const &lr : cells.leftright)
	{
		f1 << std::setw(w) << lr.first+1;
		n++;

		if (n == 6)
		{
			f1 << endl;
			n = 0;
		}
	}
	f1 << endl;

	/*Write right elements*/
	n = 0;
	for (auto const &lr : cells.leftright)
	{
		if(lr.second < 0)
		{
			f1 << std::setw(w) << 0;
		}
		else
		{
			f1 << std::setw(w) << lr.second+1;
		}
		n++;

		if(n == 6)
		{
			f1 << endl;
			n = 0;
		}
	}
	f1 << endl	<< endl;
	f1.close();
	// exit(0);
}

#endif