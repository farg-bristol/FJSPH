/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include "Var.h"
#include <string.h>
#include <tecio/TECIO.h>
// #include "../NetCDF/netcdf"
// using namespace netCDF;
// using namespace netCDF::exceptions;
// #define NC_ERR 2

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };
/*************************************************************************/
/**************************** BINARY INPUTS ******************************/
/*************************************************************************/



/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/

void Write_Binary_Timestep(const SIM &svar, const State &pnp1, 
	const uint start, const uint end, const uint zone)
{
	uint size = end - start;

	#ifdef DEBUG
		cout << start << "  " << end << "  " << size << endl;
		cout << svar.bndPts << "  " << svar.totPts << "  " << svar.simPts << endl;
	#endif

    INTEGER4 IMax = size;
    INTEGER4 JMax = 1;
    INTEGER4 KMax = 1;
	INTEGER4 ICellMax                 = 0;
    INTEGER4 JCellMax                 = 0;
    INTEGER4 KCellMax                 = 0;
    INTEGER4 DIsDouble                = 0;
    if(sizeof(ldouble) == 8)	/*If */
		DIsDouble            = 1;

    double   SolTime                  = svar.t;
    INTEGER4 StrandID;      
    INTEGER4 ParentZn                 = 0;
    INTEGER4 IsBlock                  = 1;      /* Block */
    INTEGER4 NFConns                  = 0;
    INTEGER4 FNMode                   = 0;
    INTEGER4 TotalNumFaceNodes        = 0;
    INTEGER4 TotalNumBndryFaces       = 0;
    INTEGER4 TotalNumBndryConnections = 0;
    INTEGER4 ShrConn                  = 0;
	INTEGER4 ZoneType = 0;
	INTEGER4 I = 0;

	char* group = (char*)"Boundary";
	if(zone == 1)
	{
		group = (char*)"Boundary";
		StrandID = 0;
	}
	else if (zone == 2)
	{
		group = (char*)"Fuel";
		StrandID = 1;
	}

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
    double* x = new double[size];
    for(uint dim = 0; dim < SIMDIM; ++dim)
    {
    	for(uint ii = start; ii < end; ++ii)
  			x[ii-start] = pnp1[ii].xi(dim);

 		I   = TECDAT142(&IMax, x, &DIsDouble);
    }
    
    if(svar.outform == 1)
    {
    	double* v = new double[size];
	  	double* a = new double[size];
	  	double* rho = new double[size];
	  	double* p = new double[size];
	  	for(uint ii = start; ii < end; ++ii)
  		{
	  		v[ii-start] = pnp1[ii].v.norm();
	  		a[ii-start] = pnp1[ii].f.norm();
	  		rho[ii-start] = pnp1[ii].rho;
	  		p[ii-start] = pnp1[ii].p;
	  	}
	  	I   = TECDAT142(&IMax, v, &DIsDouble);
	    I   = TECDAT142(&IMax, a, &DIsDouble);
	    I   = TECDAT142(&IMax, rho, &DIsDouble);
	    I   = TECDAT142(&IMax, p, &DIsDouble);
    }
    else if (svar.outform == 2)
    {
    	double* a = new double[size];
	  	double* Af = new double[size];
	  	double* Sf = new double[size];
	  	double* ax = new double[size];
	  	int* b = new int[size];
	  	double* theta = new double[size];

	  	for(uint ii = start; ii < end; ++ii)
	  	{
	  		a[ii-start] = pnp1[ii].f.norm();
	  		Af[ii-start] = pnp1[ii].Af.norm();
	  		Sf[ii-start] = pnp1[ii].Sf.norm();
	  		b[ii-start] = pnp1[ii].b;
	  		theta[ii-start] = pnp1[ii].theta;
	  	}

	  	I   = TECDAT142(&IMax, a, &DIsDouble);
	    I   = TECDAT142(&IMax, Af, &DIsDouble);
	    I   = TECDAT142(&IMax, Sf, &DIsDouble);
	    for(uint dim = 0; dim < SIMDIM; ++dim)
	    {
	    	for(uint ii = start; ii < end; ++ii)
  				ax[ii-start] = pnp1[ii].Af(dim);
	  			
	 		I   = TECDAT142(&IMax, ax, &DIsDouble);
	    }
	    I   = TECDAT142(&IMax, b, &DIsDouble);
	    I   = TECDAT142(&IMax, theta, &DIsDouble);
    }
}

void Init_Binary_PLT(const SIM &svar)
{
	INTEGER4 Debug      = 0;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 fileFormat = 0; // 0 == PLT, 1 == SZPLT
    INTEGER4 I          = 0; /* Used to track return codes */

    std::string file = svar.outfolder;
    file.append("/Fuel.plt");
    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    #if SIMDIM == 2
    	std::string variables = "x y";	
		if (svar.outform == 1)
		{
			variables = "x y v a rho P";
		}
		else if (svar.outform == 2)
		{
			variables = "x y a Af Sf ax ay b theta";
		}
	#endif

	#if SIMDIM == 3
		std::string variables = "x y z";  
		if (svar.outform == 1)
		{
			variables = "x y z v a rho P";
		}
		else if (svar.outform == 2)
		{
			variables = "x y z a Af Sf ax ay az b theta";
		}
	#endif

	I = TECINI142((char*)"Simulation Particles", 
						  variables.c_str(),  
		                  file.c_str(),
		                  svar.outfolder.c_str(),      /* Scratch Directory */
		                  &fileFormat,
		                  &FileType,
		                  &Debug,
		                  &VIsDouble);

    if(I == -1)
    	exit(-1);
}

void Write_Boundary_Binary(const SIM &svar, const State &pnp1)
{
	INTEGER4 Debug      = 0;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 fileFormat = 0; // 0 == PLT, 1 == SZPLT
    INTEGER4 I          = 0; /* Used to track return codes */

    std::string file = svar.outfolder;
    file.append("/Boundary.plt");
    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    #if SIMDIM == 2
    	std::string variables = "x y";	
		if (svar.outform == 1)
		{
			variables = "x y v a rho P";
		}
		else if (svar.outform == 2)
		{
			variables = "x y a Af Sf ax ay b theta";
		}
	#endif

	#if SIMDIM == 3
		std::string variables = "x y z";  
		if (svar.outform == 1)
		{
			variables = "x y z v a rho P";
		}
		else if (svar.outform == 2)
		{
			variables = "x y z a Af Sf ax ay az b theta";
		}
	#endif

	I = TECINI142((char*)"Boundary Particles", 
						  variables.c_str(),  
		                  file.c_str(),
		                  svar.outfolder.c_str(),      /* Scratch Directory */
		                  &fileFormat,
		                  &FileType,
		                  &Debug,
		                  &VIsDouble);

    if(I == -1)
    {
    	exit(-1);
	}

  	Write_Binary_Timestep(svar, pnp1, 0, svar.bndPts, 1); /*Final input: 1 = boundary, 2 = sim*/
    
    I = TECEND142();

    if(I == -1)
    	exit(-1);
}

// void Write_Group(NcGroup &zone, State &pnp1, int start, int end, int outform)
// {
// 	int size = end - start;
// 	NcDim posDim = zone.addDim("Point",size);

// 	if(outform == 5)
// 	{	/*Write fluid data*/
// 		NcVar press = zone.addVar("Pressure", ncDouble,posDim);
// 		NcVar dens = zone.addVar("Density",ncDouble,posDim);
// 		NcVar acc = zone.addVar("Acceleration",ncDouble,posDim);
// 		NcVar vel = zone.addVar("Velocity",ncDouble,posDim);

// 		press.putAtt("units","Pa");
// 		dens.putAtt("units","kg/m^3");
// 		acc.putAtt("units","m/s^2");
// 		vel.putAtt("units","m/s");

// 		double* p = new double[size]{0.0};
// 		double* rho = new double[size]{0.0};
// 		double* f = new double[size]{0.0};
// 		double* v = new double[size]{0.0};

// 		for(int i = start; i < end; ++i)
// 		{
// 			p[i-start] = pnp1[i].p;
// 			rho[i-start] = pnp1[i].rho;
// 			f[i-start] = pnp1[i].f.norm();
// 			v[i-start] = pnp1[i].v.norm();
// 		}	

// 		press.putVar(p);
// 		dens.putVar(rho);
// 		acc.putVar(f);
// 		vel.putVar(v);
// 	}

// 	if(outform == 6)
// 	{	/*Write research data*/
// 		NcVar aero = zone.addVar("Aerodynamic Force", ncDouble,posDim);
// 		NcVar surf = zone.addVar("Surface Tension", ncDouble,posDim);
// 		NcVar Sx = zone.addVar("S_x", ncDouble,posDim);
// 		NcVar Sy = zone.addVar("S_y", ncDouble,posDim);
// 		NcVar B = zone.addVar("Boundary",ncInt,posDim);
// 		NcVar theta = zone.addVar("Theta",ncInt,posDim);

// 		aero.putAtt("units","N");
// 		surf.putAtt("units","N");
// 		Sx.putAtt("units","N");
// 		Sy.putAtt("units","N");
	
// 		double* af = new double[size]{0.0};
// 		double* sf = new double[size]{0.0};
// 		double* sx = new double[size]{0.0};
// 		double* sy = new double[size]{0.0};
// 		int* b = new int[size]{0};
// 		int* t = new int[size]{0};

// 		for(int i = start; i < end; ++i)
// 		{
// 			af[i-start] = pnp1[i].Af.norm();
// 			sf[i-start] = pnp1[i].Sf.norm();
// 			sx[i-start] = pnp1[i].Sf(0);
// 			sy[i-start] = pnp1[i].Sf(1);
// 			b[i-start] = pnp1[i].b;
// 			t[i-start] = pnp1[i].theta;
// 		}		

// 		aero.putVar(af);
// 		surf.putVar(sf);
// 		Sx.putVar(sx);
// 		Sy.putVar(sy);
// 		B.putVar(b);
// 		theta.putVar(t);

// 		if (SIMDIM == 3)
// 		{
// 			NcVar Sz = zone.addVar("S_z",ncDouble,posDim);
// 			Sz.putAtt("units","m");
// 			double* sz = new double[size]{0.0};
// 			for(int i = start; i < end; ++i)
// 				sz[i-start] = pnp1[i].Sf[2];

// 			Sz.putVar(sz);
// 		}
// 	}



// 	NcVar xVar = zone.addVar("x (m)",ncDouble,posDim);
// 	NcVar yVar = zone.addVar("y (m)",ncDouble,posDim);
	
// 	xVar.putAtt("units","m");
// 	yVar.putAtt("units","m");
// 	double* x = new double[size]{0.0};
// 	double* y = new double[size]{0.0};
	
// 	for(int i = start; i < end; ++i)
// 	{
// 		x[i-start] = pnp1[i].xi[0];
// 		y[i-start] = pnp1[i].xi[1];	
// 	}

// 	xVar.putVar(x);
// 	yVar.putVar(y);

// 	if (SIMDIM == 3)
// 	{
// 		NcVar zVar = zone.addVar("z (m)",ncDouble,posDim);
// 		zVar.putAtt("units","m");
// 		double* z = new double[size]{0.0};
// 		for(int i = start; i < end; ++i)
// 			z[i-start] = pnp1[i].xi[2];

// 		zVar.putVar(z);
// 	}
	
// }

// void Write_CDF_File(NcFile &nf, SIM &svar, State &pnp1)
// {	/*Write the timestep to a group*/
// 	NcGroup ts(nf.addGroup(std::to_string(svar.t)));

// 	/*Write the Particles group*/
// 	Write_Group(ts, pnp1, svar.bndPts, svar.totPts, svar.outform);
// }

// void Write_Boundary_CDF(SIM &svar, State &pnp1)
// {	
// 	std::string file = svar.outfolder;
// 	file.append("/Boundary.h5");
// 	NcFile bound(file, NcFile::replace);
// 	NcGroup part(bound.addGroup("Boundary"));
// 	Write_Group(bound, pnp1, 0, svar.bndPts, svar.outform);
// }

#endif