/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef BINARYIO_H
#define BINARYIO_H

#include "Var.h"
#include "Neighbours.h"
#include "Kernel.h"
#include <string.h>
#include <tecio/TECIO.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
using std::string;

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

/*************************************************************************/
/**************************** BINARY INPUTS ******************************/
/*************************************************************************/
struct Point {
	Point(){};

	Point(const uint size, const uint type)
	{
		dataType = type;
		verts = vector<StateVecD>(size);
		press = vector<ldouble>(size);
		rho = vector<ldouble>(size);
		m = vector<ldouble>(size);

		if (type == 1)
		{
			vnorm = vector<ldouble>(size);
			anorm = vector<ldouble>(size);
		}
		else if (type == 2)
		{
			vel = vector<StateVecD>(size);
			acc = vector<StateVecD>(size);
		}
		else if(type == 3)
		{
			vel = vector<StateVecD>(size);
			acc = vector<StateVecD>(size);
			cellV = vector<StateVecD>(size);
			cellP = vector<ldouble>(size);
			cellRho = vector<ldouble>(size);
		}
	}

	void operator=(const Point& pi)
	{
		verts = pi.verts; vel = pi.vel; acc = pi.acc;
		vnorm = pi.vnorm; anorm = pi.anorm; press = pi.press;
		rho = pi.rho; m = pi.m; cellV = pi.cellV; 
		cellP = pi.cellP; cellRho = pi.cellRho;
		time = pi.time;
	}

	void append(const Point& pi)
	{
		verts.insert(verts.end(),pi.verts.begin(),pi.verts.end());
		vel.insert(vel.end(),pi.vel.begin(),pi.vel.end());
		acc.insert(acc.end(),pi.acc.begin(),pi.acc.end());
		vnorm.insert(vnorm.end(),pi.vnorm.begin(),pi.vnorm.end());
		anorm.insert(anorm.end(),pi.anorm.begin(),pi.anorm.end());
		press.insert(press.end(),pi.press.begin(),pi.press.end());
		rho.insert(rho.end(),pi.rho.begin(),pi.rho.end());
		m.insert(m.end(),pi.m.begin(),pi.m.end());
		cellV.insert(cellV.end(),pi.cellV.begin(),pi.cellV.end());
		cellP.insert(cellP.end(),pi.cellP.begin(),pi.cellP.end());
		cellRho.insert(cellRho.end(),pi.cellRho.begin(),pi.cellRho.end());
		time = pi.time;
	}

	uint size()
	{
		return verts.size();
	}

	uint dataType;
	ldouble time;
	vector<StateVecD> verts, vel, acc;
	vector<ldouble> vnorm, anorm;
	vector<ldouble> press, rho, m;
	vector<StateVecD> cellV;
	vector<ldouble> cellP, cellRho;
};

void* Read_Binary_Info(SIM& svar, const string filename)
{
	string file = svar.outfolder;
	file.append(filename);
	INTEGER4 I;
	void* inputHandle = NULL;
	I = tecFileReaderOpen(file.c_str(),&inputHandle);
	if(I == -1)
	{
		cout << "Error opening szplt file. Path:" << endl;
		cout << file << endl;
		exit(-1);
	}

	INTEGER4 numVars;
	I = tecDataSetGetNumVars(inputHandle, &numVars);
	cout << "Number of variables in output file: " << numVars << endl;
	if(numVars == 2 || numVars == 3)
    {
    	cout << "Only basic data has been output. No post processing can be done." << endl;
    	exit(-1);
    }

	vector<string> varNames;
    for (INTEGER4 var = 1; var <= numVars; ++var)
    {	
        char* name = NULL;
        I = tecVarGetName(inputHandle, var, &name);
        varNames.emplace_back(name);
        tecStringFree(&name);
    }

    /*Check if the data has components or not*/
    uint dataType=0;
	if(varNames[2] == "v_x" || varNames[3] == "v_x")
	{	/*Data has components*/
		if(varNames.size()>13)
			dataType = 3;
		else if (varNames.size() > 8)
			dataType = 2;
		else/*Data has no components*/
			dataType = 1;
	}

	if(svar.outform != dataType)
	{
		cout << "Warning: Output data is different from defined in the settings" << endl;
	}
	svar.outform = dataType;

	#if SIMDIM == 2
		if(svar.numVars == 8 || svar.numVars == 12 || svar.numVars == 18)
	    {	/*Output type is 3D fluid data*/
	    	cout << "Error: Performing post processing using 2D code on 3D solution." << endl;
	    	exit(-1);
	    }
    #endif

	#if SIMDIM == 3
	    if(svar.numVars == 7 || svar.numVars == 9 || svar.numVars == 14)
	    {	/*Output type is 3D fluid data*/
	    	cout << "Error: Performing post processing using 3D code on 2D solution." << endl;
	    	exit(-1);
	    }

	#endif

    INTEGER4 numZones;
    I = tecDataSetGetNumZones(inputHandle, &numZones);
    cout << "Number of zones (timesteps) in output file: " << numZones << endl;

    if(filename == "Boundary.szplt") 
    {
    	if(svar.boutform == 1)
    	{
    		if (svar.afterSim == 1)
		    {
		    	if(static_cast<uint>(numZones) != svar.Nframe)
			    {
		cout << "Warning: number of written timesteps is different from the frame count." << endl;
			    }
		    }
		    
		    svar.Nframe = numZones;
    	}
    }
    else
    {
	    if (svar.afterSim == 1)
	    {
	    	if(static_cast<uint>(numZones) != svar.Nframe)
		    {
		cout << "Warning: number of written timesteps is different from the frame count." << endl;
		    }
	    }
	    
	    svar.Nframe = numZones;
	}

    return inputHandle;
}


Point Read_Binary_Timestep(void* inputHandle, const INTEGER4 zoneNum, const uint dataType)
{
// variables = "x y z rho P m v_x v_y v_z a_x a_y a_z Cell_ID Cell_Vx Cell_Vy Cell_Vz Cell_P Cell_Rho";
    cout << "Reading zone number: " << zoneNum << endl;

    INTEGER4 I;
    INTEGER8 iMax, jMax, kMax;
    I = tecZoneGetIJK(inputHandle, zoneNum, &iMax, &jMax, &kMax);
    if(I == -1)
    {
    	cout << "Reading zone data caused a tecplot error. Stopping." << endl;
    	exit(-1);
    }

    cout << "Number of particles in zone: " << iMax << endl;

    Point tsData(iMax,dataType);
    I = tecZoneGetSolutionTime(inputHandle, zoneNum, &tsData.time);
    vector<double> xVar(iMax);
    vector<double> yVar(iMax);
        
    INTEGER4 varCount = 1;
    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &xVar[0]);
    ++varCount;
    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &yVar[0]);
	++varCount;

    #if SIMDIM == 3
	    vector<double> zVar(iMax);
	    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &zVar[0]);
	    ++varCount;
    #endif

	#pragma omp parallel for
	for(uint ii = 0; ii < xVar.size(); ++ii)
	{
		tsData.verts[ii](0) = xVar[ii];
		tsData.verts[ii](1) = yVar[ii];
		
		#if SIMDIM == 3
			tsData.verts[ii](2) = zVar[ii];
		#endif

		// cout << tsData.verts[ii](0) << "  " << tsData.verts[ii](1) << endl;
	}

	/*Get density, pressure and mass*/
	I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.rho[0]);
	++varCount;
	I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.press[0]);
	++varCount;    
    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.m[0]);
	++varCount;  

	/*Check if the data has components or not*/
	if(dataType == 2 || dataType == 3)
	{	/*Data has components*/

	 	vector<double> xVel(iMax);
	    vector<double> yVel(iMax);
	    vector<double> xAcc(iMax);
	    vector<double> yAcc(iMax);

	    /*Get velocity components*/
		I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &xVel[0]);
	    ++varCount;
	    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &yVel[0]);
		++varCount;

		#if SIMDIM == 3
			vector<double> zVel(iMax);
		    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &zVel[0]);
		    ++varCount;
		
		#endif

		/*Get acceleration components*/
		I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &xAcc[0]);
	    ++varCount;
	    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &yAcc[0]);
		++varCount;

		#if SIMDIM == 3
		    vector<double> zAcc(iMax);
		    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &zAcc[0]);
			++varCount;	
		#endif


		/*Write data to a StateVecD vector*/
		#pragma omp parallel for
		for (uint ii = 0; ii < iMax; ++ii)
		{
			tsData.vel[ii](0) = xVel[ii];
			tsData.vel[ii](1) = yVel[ii];

			tsData.acc[ii](0) = xAcc[ii];
			tsData.acc[ii](1) = yAcc[ii];

			#if SIMDIM == 3
			tsData.vel[ii](2) = zVel[ii];
			tsData.acc[ii](2) = zAcc[ii];
			#endif
		}

		if(dataType == 3)
		{	/*Get the cell information for the points*/
			vector<double> xVel(iMax);
		    vector<double> yVel(iMax);

		    /*Get velocity components*/
			I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &xVel[0]);
		    ++varCount;
		    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &yVel[0]);
			++varCount;

			#if SIMDIM == 3
				vector<double> zVel(iMax);
			    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &zVel[0]);
			    ++varCount;
			#endif

			#pragma omp parallel for
			for(uint ii = 0; ii < xVar.size(); ++ii)
			{
				tsData.cellV[ii](0) = xVel[ii];
				tsData.cellV[ii](1) = yVel[ii];
				
				#if SIMDIM == 3
					tsData.cellV[ii](2) = zVel[ii];
				#endif
			}

			/*Get density and pressure*/
			I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.cellRho[0]);
		    ++varCount;
		    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.cellP[0]);
			++varCount;

		}
	}
	else
	{
	    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.vnorm[0]);
		++varCount;
		I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.anorm[0]);
		++varCount;
	}
	

	
    
	/*End of that timestep*/
	return tsData;
}


typedef class TECMESH {
	public:
		void FindGridSize(const SIM& svar)
		{
			if(svar.afterSim == 1)
			{	/* If post processing is being done immediately */
				/* the data is already stored*/
				minC = svar.minC;
				maxC = svar.maxC;
			}
			else
			{	/*If not, data needs attaining from frame info*/
				string file = svar.outfolder;
				file.append("frame.info");
				std::ifstream f1(file,std::ios::in);

				/*The max and min will be stored at the end of file,*/
				/*so find the size of the file*/
				if(f1.is_open())
				{
					uint lsize = 0;
					while(f1.ignore(std::numeric_limits<std::streamsize>::max(),'\n'))	
						lsize++;
					
					cout << "Found file size: " << lsize 
						 << ". Now getting max dimensions." << endl;
					f1.clear();
					GotoLine(f1, lsize-3);
					string line;
					getline(f1, line);
					// cout << line << endl;
					std::istringstream sline(line);
					StateVecD temp;
					sline >> temp(0); sline >> temp(1);
					#if SIMDIM == 3
						sline >> temp(2);
					#endif
					minC = temp;

					getline(f1, line);
					getline(f1, line);
					// cout << line << endl;
					sline = std::istringstream(line);
					sline >> temp(0); sline >> temp(1);
					#if SIMDIM == 3
						sline >> temp(2);
					#endif
					maxC = temp;
				}
				else
				{
					cout << "Couldn't open frame.info file. Please check it exists in the output folder" << endl;
					exit(-1);
				}
			}

			nX = ceil((maxC(0)-minC(0)+4*svar.postRadius)/svar.cellSize);
			step(0) = (maxC(0)-minC(0)+4*svar.postRadius)/ldouble(nX-1);
			nY = ceil((maxC(1)-minC(1)+4*svar.postRadius)/svar.cellSize);
			step(1) = (maxC(1)-minC(1)+4*svar.postRadius)/ldouble(nY-1);

			#if SIMDIM == 3
				nZ = ceil((maxC(2)-minC(2)+4*svar.postRadius)/svar.cellSize);
				step(2) = (maxC(2)-minC(2)+4*svar.postRadius)/ldouble(nZ);
				nVerts = nX*nY*nZ;
				nCells = (nX-1)*(nY-1)*(nZ-1);
				nConns = 8*nCells;
				nFaces = 6;
			#else 
				nVerts = nX*nY;
				nCells = (nX-1)*(nY-1);
				nConns = 4*nCells;
				nFaces = 4;
			#endif

			cout << "Grid Size Found. Vertices: " << nVerts << " Cells: " << nCells << endl;
			
			dataType = svar.outform;
			init(nVerts); 			
		}

		void Create_Grid(const SIM& svar)
		{

			cout << "Creating grid..." << endl;

			#if SIMDIM == 3
				
			verts = vector<StateVecD>(nVerts);
			for(INTEGER4 jj = 0; jj < nY; jj++)
				for(INTEGER4 ii = 0; ii < nX; ii++)
					for(INTEGER4 kk = 0; kk < nZ; kk++)
					{
						ldouble x = minC(0)-2*svar.postRadius + step(0)*ii;
						ldouble y = minC(1)-2*svar.postRadius + step(1)*jj;
						ldouble z = minC(2)-2*svar.postRadius + step(2)*kk;
						verts[index2(ii,jj)] = StateVecD(x,y,z);
					}

			// cout << maxC(0)+2*svar.postRadius << "  " << maxC(1)+2*svar.postRadius << endl;
			// cout << verts[nVerts-1](0) << "  " << verts[nVerts-1](1) << endl;

			xC = new double[verts.size()];
			yC = new double[verts.size()];
			zC = new double[verts.size()];

			for(uint jj = 0; jj < verts.size(); ++jj)
			{
				xC[jj] = verts[jj](0);
				yC[jj] = verts[jj](1);
				zC[jj] = verts[jj](2);
			}

			#else

			verts = vector<StateVecD>(nVerts);
			for(INTEGER4 jj = 0; jj < nY; jj++)
				for(INTEGER4 ii = 0; ii < nX; ii++)
				{
					ldouble x = minC(0)-2*svar.postRadius + step(0)*ii;
					ldouble y = minC(1)-2*svar.postRadius + step(1)*jj;
					verts[index2(ii,jj)] = StateVecD(x,y);
				}

			// cout << maxC(0)+2*svar.postRadius << "  " << maxC(1)+2*svar.postRadius << endl;
			// cout << verts[nVerts-1](0) << "  " << verts[nVerts-1](1) << endl;

			xC = new double[verts.size()];
			yC = new double[verts.size()];

			for(uint jj = 0; jj < verts.size(); ++jj)
			{
				xC[jj] = verts[jj](0);
				yC[jj] = verts[jj](1);
			}

			#endif
			
			if(verts.size() != static_cast<uint>(nVerts))
			{
				cout << "Mismatch of post process mesh dimensions." << endl;
				cout << "Mesh actual size: " << verts.size() << " Predicted: " << nVerts << endl;
				exit(-1);
			}
			// FindCellNeighbours()

			cout << "Grid made. Vertices: " << nVerts << " Cells: " << nCells << endl;
		}

		void Interp_Data(const Point& fluid, const FLUID& fvar)
		{
			INTEGER4 I;
			INTEGER4 DIsDouble                = 0;
		    if(sizeof(ldouble) == 8)
				DIsDouble            = 1;

		    const INTEGER4 strandID = 1;     
		    const INTEGER4 parentZn                 = 0;
		    const INTEGER4 isBlock                  = 1;      /* Block */
		    const INTEGER4 nFConns                  = 0;
		    const INTEGER4 fNMode                   = 0;
		    const INTEGER4 shrConn                  = 0;

		 //    #if SIMDIM == 2
			// const INTEGER4 zoneType = 0;
			// #endif
			// #if SIMDIM == 3
			// 	const INTEGER4 zoneType = 5;
			// #endif

			const INTEGER4 zoneType = 0;

			const INTEGER4 iMax = nX;
			const INTEGER4 jMax = nY;
			#if SIMDIM == 3
			 const INTEGER4 kMax = nZ;
			#else
			 const INTEGER4 kMax = 1;
			#endif

			// cout << fluid.size() << endl;
			/*Fluid from the timestep is read. build a neighbour tree.*/
			const Vec_Tree INDEX(SIMDIM,fluid.verts,50);

			/*Now that there is the neighbour tree, interpolate the properties to the grid*/
			const nanoflann::SearchParams params;
			const ldouble search_radius = fvar.sr;

			#pragma omp parallel for
			for(uint ii=0; ii < verts.size(); ++ii)
			{
				StateVecD testp = verts[ii];
				std::vector<std::pair<size_t, ldouble>> matches; /* Nearest Neighbour Search*/
				INDEX.index->radiusSearch(&testp[0], search_radius, matches, params);

				ldouble ktemp = 0.0;
				ldouble pressure = 0.0;
				ldouble mass = 0.0;
				ldouble dens = 0.0;

				if (dataType == 1)
				{
					ldouble vels = 0.0;
					ldouble accs = 0.0;
					for(auto temp:matches)
					{
						uint jj = temp.first;
						/*Find the kernel*/
						ldouble Rij = (fluid.verts[jj] - testp).norm();
						ldouble k = W2Kernel(Rij,fvar.H,1.0)/kernsum;

						ktemp += k;
						pressure += k*fluid.press[jj];
						mass += k*fluid.m[jj];
						dens += k*fluid.rho[jj];
						vels += k*fluid.vnorm[jj];
						accs += k*fluid.anorm[jj];
					}
					input(ii,ktemp,vels,accs,pressure,dens,mass);
				}

				else if (dataType == 2)
				{
					StateVecD vels = StateVecD::Zero();
					StateVecD accs = StateVecD::Zero();
					for(auto temp:matches)
					{
						uint jj = temp.first;
						/*Find the kernel*/
						ldouble Rij = (fluid.verts[jj] - testp).norm();
						ldouble k = W2Kernel(Rij,fvar.H,1.0)/kernsum;

						ktemp += k;
						pressure += k*fluid.press[jj];
						mass += k*fluid.m[jj];
						dens += k*fluid.rho[jj];
						vels += k*fluid.vel[jj];
						accs += k*fluid.acc[jj];
					}
					input(ii,ktemp,vels,accs,pressure,dens,mass);
				}
				else if (dataType == 3)
				{
					StateVecD vels = StateVecD::Zero();
					StateVecD accs = StateVecD::Zero();
					StateVecD cV = StateVecD::Zero();
					ldouble cR = 0.0;
					ldouble cP = 0.0;
					for(auto temp:matches)
					{
						uint jj = temp.first;
						/*Find the kernel*/
						ldouble Rij = (fluid.verts[jj] - testp).norm();
						ldouble k = W2Kernel(Rij,fvar.H,1.0)/kernsum;

						ktemp += k;
						vels += k*fluid.vel[jj];
						accs += k*fluid.acc[jj];
						pressure += k*fluid.press[jj];
						dens += k*fluid.rho[jj];
						mass += k*fluid.m[jj];
						cV += k*fluid.cellV[jj];
						cR +=k*fluid.cellRho[jj];
						cP +=k*fluid.cellP[jj];
					}

					input(ii,ktemp,vels,accs,pressure,dens,mass,cV,cP,cR);
				}
			}

			
			/*Data has been interpolated to the grid, so now write to the file*/
			time = fluid.time;
			string zoneName = std::to_string(time);
			I = TECZNE142(zoneName.c_str(),
			                  &zoneType,
			                  &iMax,
			                  &jMax,
			                  &kMax,
			                  0,
			                  0,
			                  0,
			                  &time,
			                  &strandID,
			                  &parentZn,
			                  &isBlock,
			                  &nFConns,
			                  &fNMode,
			                  0,              /* TotalNumFaceNodes */
			                  0,              /* NumConnectedBoundaryFaces */
			                  0,              /* TotalNumBoundaryConnections */
			                  NULL,           /* PassiveVarList */
			                  NULL, 		  /* ValueLocation = Nodal */
			                  NULL,           /* SharVarFromZone */
			                  &shrConn);

			if(I == -1)
			{
				cout << "A tecplot error occured when making a new zone" << endl;
				exit(-1);
			}

			I = TECDAT142(&nVerts, xC, &DIsDouble);
			I = TECDAT142(&nVerts, yC, &DIsDouble);

			#if SIMDIM == 3
				I = TECDAT142(&nVerts, zC, &DIsDouble);
			#endif

			I = TECDAT142(&nVerts, kvalue, &DIsDouble);
			I = TECDAT142(&nVerts, rho, &DIsDouble);
			I = TECDAT142(&nVerts, press, &DIsDouble);
			I = TECDAT142(&nVerts, m, &DIsDouble);

			
			if(dataType == 1)
			{
				I = TECDAT142(&nVerts, vnorm, &DIsDouble);
				I = TECDAT142(&nVerts, anorm, &DIsDouble);
			}
			else if(dataType == 2 || dataType == 3)
			{
				I = TECDAT142(&nVerts, vX, &DIsDouble);
				I = TECDAT142(&nVerts, vY, &DIsDouble);
				#if SIMDIM == 3
					I = TECDAT142(&nVerts, vZ, &DIsDouble);
				#endif

				I = TECDAT142(&nVerts, aX, &DIsDouble);
				I = TECDAT142(&nVerts, aY, &DIsDouble);
				#if SIMDIM == 3
					I = TECDAT142(&nVerts, aZ, &DIsDouble);
				#endif

				if (dataType == 3)
				{
					I = TECDAT142(&nVerts, cVX, &DIsDouble);
					I = TECDAT142(&nVerts, cVY, &DIsDouble);
					#if SIMDIM == 3
						I = TECDAT142(&nVerts, cVZ, &DIsDouble);
					#endif

					I = TECDAT142(&nVerts, cellRho, &DIsDouble);
					I = TECDAT142(&nVerts, cellP, &DIsDouble);
				}
			}

			INTEGER4 numZonesToRetain = 0;
		    I = TECFLUSH142(&numZonesToRetain,NULL);
		}

		void Write_Data(SIM& svar, FLUID& fvar)
		{

			/*Read in the timestep*/
			void* boundFile;
			if(svar.Bcase != 0 && svar.Bcase !=5)
				boundFile = Read_Binary_Info(svar, "Boundary.szplt");

			void* fluidFile = Read_Binary_Info(svar, "Fuel.szplt");


			#if DEBUG
				const INTEGER4 Debug = 1;
			#else 
				const INTEGER4 Debug = 0;
			#endif
			INTEGER4 VIsDouble       = 0; 
			if(sizeof(ldouble) == 8)
				VIsDouble            = 1;
		    
		    const INTEGER4 FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
		    const INTEGER4 fileFormat = 1; // 0 == PLT, 1 == SZPLT
		    INTEGER4 I          = 0; /* Used to track return codes */
		    int* valueLocation = NULL;


		    string file = svar.outfolder;
		    file.append("grid.szplt");

		    #if SIMDIM == 2
		    	std::string variables = "x y";	
				if (dataType == 1)
				{
					variables = "x y k rho P m v a";
					valueLocation = (int*)realloc(valueLocation, 7 * sizeof(int));
					nVar = 7;
				}
				else if (dataType == 2)
				{
					variables = "x y k rho P m v_x v_y a_x a_y";
					valueLocation = (int*)realloc(valueLocation, 9 * sizeof(int));
					nVar = 9;
				}
				else if (dataType == 3)
				{
		variables = "x y k rho P m v_x v_y a_x a_y Cell_Vx Cell_Vy Cell_Rho Cell_P";
					valueLocation = (int*)realloc(valueLocation, 13 * sizeof(int));
					nVar = 13;
				}
			#endif

			#if SIMDIM == 3
				std::string variables = "x y z";  
				if (dataType == 1)
				{
					variables = "x y z k rho P m v a";
					valueLocation = (int*)realloc(valueLocation, 8 * sizeof(int));
					nVar = 8;
				}
				else if (dataType == 2)
				{
					variables = "x y z k rho P m v_x v_y v_z a_x a_y a_z";
					valueLocation = (int*)realloc(valueLocation, 12 * sizeof(int));
					nVar = 12;
				}
				else if (dataType == 3)
				{
					variables = 
		"x y z k rho P m v_x v_y v_z a_x a_y a_z Cell_Vx Cell_Vy Cell_Vz Cell_Rho Cell_P";
					valueLocation = (int*)realloc(valueLocation, 17 * sizeof(int));
					nVar = 17;
				}
			#endif

			for(uint ii = 0; ii < nVar; ++ii)
				valueLocation[ii] = 1;

			cout << "Opening the grid file." << endl;
			I = TECINI142((char*)"SPH Fluid Grid",
						  variables.c_str(),
						  file.c_str(),
						  (char*)".",
						  &fileFormat,
		                  &FileType,
		                  &Debug,
		                  &VIsDouble);
			if(I == -1)
			{
				cout << "Error opening output grid file for post processing." << endl;
				exit(-1);
			}

			GetKernelSum(svar,fvar,fluidFile);
			
			/*Read the data from the files into something.*/
			if(svar.Bcase !=0 && svar.Bcase != 5 && svar.boutform == 1)
			{	/*The boundary is being written each timestep*/
				for( INTEGER4 zone = 1; zone <= static_cast<INTEGER4>(svar.Nframe); ++zone)
				{
					Point boundary = Read_Binary_Timestep(boundFile, zone, dataType);
					Point fuel = Read_Binary_Timestep(fluidFile, zone, dataType);

					Point fluid = boundary;
					fluid.append(fuel);
					
					Interp_Data(fluid, fvar);
				}
			}
			else
			{
				Point boundary;

				if(svar.Bcase != 0 && svar.Bcase !=5)
				{
					boundary = Read_Binary_Timestep(boundFile, 1, dataType);
					/*Close the boundary file*/
					I = tecFileReaderClose(&boundFile);
				}

				for( INTEGER4 zone = 1; zone <= static_cast<INTEGER4>(svar.Nframe); ++zone)
				{
					Point fuel = Read_Binary_Timestep(fluidFile, zone, dataType);

					Point fluid = boundary;
					fluid.append(fuel);
					
					Interp_Data(fluid, fvar);
				}

			}

			TECEND142();
			free(valueLocation);
			
		}
		

		void DoPostProcessing(SIM& svar, FLUID& fvar)
		{
			/*Check the folder for the grid szplt files*/
			string file = svar.outfolder;
	  		file.append("grid.szplt.szdat");
	  		struct stat info;
	  		if(stat( file.c_str(), &info ) == 0)
	  		{
		  		string cmd = "exec rm -r \"";
		  		cmd.append(svar.outfolder);
		  		cmd.append("\"*.szplt.sz*");
		  		if(system(cmd.c_str()))
		  		{
			    	cout << "System command failed to execute." << endl;
			    	exit(-1);
			    }
			}

			FindGridSize(svar);
			Create_Grid(svar);
			Write_Data(svar,fvar);
		}

		
	private:

		void GetKernelSum(const SIM& svar, const FLUID& fvar, void* fluidFile)
		{

			Point fluid = Read_Binary_Timestep(fluidFile, 1, dataType);
			/*Fluid from the timestep is read. build a neighbour tree.*/
			const Vec_Tree INDEX(SIMDIM,fluid.verts,50);

			/*Now that there is the neighbour tree, interpolate the properties to the grid*/
			const nanoflann::SearchParams params;
			const ldouble search_radius = fvar.sr;

			/*Find the maximum kernel sum to normalise by*/
			ldouble ksum = 0.0;
			#pragma omp parallel for
			for(uint ii=0; ii < verts.size(); ++ii)
			{
				StateVecD testp = verts[ii];
				std::vector<std::pair<size_t, ldouble>> matches; /* Nearest Neighbour Search*/
				INDEX.index->radiusSearch(&testp[0], search_radius, matches, params);

				ldouble ktemp = 0.0;
				for(auto temp:matches)
				{
					uint jj = temp.first;
					/*Find the kernel*/
					ldouble Rij = (fluid.verts[jj] - testp).norm();
					ldouble k = W2Kernel(Rij,fvar.H,1.0);

					ktemp += k;
				}

				if (ktemp > ksum)
					ksum = ktemp;
			}

			kernsum = ksum;
			cout << ksum << endl;
		}

		const uint index2(const uint ii, const uint jj)
		{
			return(ii + jj*nX);
		}
		const uint index3(const uint ii, const uint jj, const uint kk)
		{
			return ((kk*nY + jj)*nX + ii);
		}

		void init(const uint size)
		{
			press = new ldouble[size];
			rho = new ldouble[size];
			m = new ldouble[size];
			kvalue = new ldouble[size];

			if (dataType == 1)
			{
				vnorm = new ldouble[size]; anorm = new ldouble[size]; 
			}
			else if(dataType == 2 || dataType == 3)
			{
				vX = new ldouble[size]; vY = new ldouble[size]; 
				aX = new ldouble[size]; aY = new ldouble[size];
				#if SIMDIM == 3
				vZ = new ldouble[size]; aZ = new ldouble[size]; 
				#endif  

				if(dataType == 3)
				{	
					cVX = new ldouble[size]; cVY = new ldouble[size];  
					cellP = new ldouble[size]; cellRho = new ldouble[size]; 
					#if SIMDIM == 3
						cVZ = new ldouble[size];
					#endif
				}
			}
		}
		void input(const uint ii, const ldouble k, const ldouble& vels, 
			const ldouble& acc, const ldouble P, const ldouble dens, const ldouble mass)
		{
			vnorm[ii] = vels; anorm[ii] = acc; kvalue[ii] = k;
			press[ii] = P; rho[ii] = dens; m[ii] = mass;
		}

		void input(const uint ii, const ldouble k, const StateVecD& vels, 
		const StateVecD& acc, const ldouble P, const ldouble dens, const ldouble mass)
		{
			vX[ii] = vels(0); vY[ii] = vels(1); aX[ii] = acc(0); aY[ii] = acc(1); 
			press[ii] = P; rho[ii] = dens; m[ii] = mass; kvalue[ii] = k;
			#if SIMDIM ==3
			vZ[ii] = vels(2);  aZ[ii] = acc(2); 
			#endif
		}

		void input(const uint ii, const ldouble k, const StateVecD& vels, const StateVecD& acc, 
		const ldouble P, const ldouble dens, const ldouble mass, const StateVecD& cellV, 
		const ldouble cP, const ldouble cR)
		{
			kvalue[ii] = k;
			vX[ii] = vels(0); vY[ii] = vels(1); aX[ii] = acc(0); aY[ii] = acc(1); 
			press[ii] = P; rho[ii] = dens; m[ii] = mass; cellP[ii] = cP; cellRho[ii] = cR;
			cVX[ii] = cellV(0); cVY[ii] = cellV(1); 
			#if SIMDIM ==3
			vZ[ii] = vels(2);  aZ[ii] = acc(2); cVZ[ii] = cellV(2);
			#endif
		}

		std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
		{
		    file.seekg(file.beg);
		    for(uint ii=0; ii < num - 1; ++ii){
		        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
		    }
		    return file;
		}

		StateVecD minC, maxC, step;
		INTEGER4 nVerts, nCells, nFaces, nConns, nX, nY, nZ;
		ldouble time;
		INTEGER4 zone;
		uint dataType, nVar;
		vector<StateVecD> verts;

		ldouble kernsum;
		ldouble* xC;
		ldouble* yC;
		ldouble* zC;
		ldouble* kvalue;
		ldouble* vnorm;
		ldouble* anorm;
		ldouble* vX;
		ldouble* vY;
		ldouble* vZ;
		ldouble* aX;
		ldouble* aY;
		ldouble* aZ;
		ldouble* cVX;
		ldouble* cVY;
		ldouble* cVZ;
		ldouble* press;
		ldouble* rho;
		ldouble* m;
		ldouble* cellP;
		ldouble* cellRho;
}TECMESH;



/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/

void Write_Binary_Timestep(const SIM& svar, const State& pnp1, 
	const uint start, const uint end, const char* group, const INTEGER4 StrandID)
{
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
    INTEGER4 DIsDouble                = 0;
    if(sizeof(ldouble) == 8)
		DIsDouble            = 1;

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
    double* x = new double[size];
    
    for(uint dim = 0; dim < SIMDIM; ++dim)
    {
    	#pragma omp parallel for
    	for(uint ii = start; ii < end; ++ii)
  			x[ii-start] = pnp1[ii].xi(dim)/svar.scale;

 		I   = TECDAT142(&IMax, x, &DIsDouble);
    }

	if(svar.outform != 0)
	{
		double* rho = new double[size];
	    double* p = new double[size];
		double* m = new double[size];

		#pragma omp parallel for
	  	for(uint ii = start; ii < end; ++ii)
		{
			rho[ii-start] = pnp1[ii].rho;
	  		p[ii-start] = pnp1[ii].p;
	  		m[ii-start] = pnp1[ii].m;
		}

	    I   = TECDAT142(&IMax, rho, &DIsDouble);
		I   = TECDAT142(&IMax, p, &DIsDouble);
		I   = TECDAT142(&IMax, m, &DIsDouble);
    }

    if(svar.outform == 1)
    {
    	double* v = new double[size];
	  	double* a = new double[size];
	  	
	  	#pragma omp parallel for
	  	for(uint ii = start; ii < end; ++ii)
  		{
	  		v[ii-start] = pnp1[ii].v.norm();
	  		a[ii-start] = pnp1[ii].f.norm();
	  	}

	  	I   = TECDAT142(&IMax, v, &DIsDouble);
	    I   = TECDAT142(&IMax, a, &DIsDouble);
    }
    else if (svar.outform == 2 || svar.outform == 3)
    {
// variables = "x y z rho P m v_x v_y v_z a_x a_y a_z Cell_ID Cell_Vx Cell_Vy Cell_Vz Cell_P Cell_Rho";
    	double* ax = new double[size];
    	double* ay = new double[size];
    	double* vx = new double[size];
    	double* vy = new double[size];
	  	
	  	#pragma omp parallel for
	  	for(uint ii = start; ii < end; ++ii)
	  	{
	  		ax[ii-start] = pnp1[ii].f(0);
	  		ay[ii-start] = pnp1[ii].f(1);
	  		vx[ii-start] = pnp1[ii].v(0);
	  		vy[ii-start] = pnp1[ii].v(1);
	  	}

	  	// #pragma omp critical
	  	I   = TECDAT142(&IMax, vx, &DIsDouble);
		I   = TECDAT142(&IMax, vy, &DIsDouble);
		#if SIMDIM == 3
			double* vz = new double[size];
			#pragma omp parallel for
	  		for(uint ii = start; ii < end; ++ii)
				vz[ii-start] = pnp1[ii].f(2);
			I   = TECDAT142(&IMax, vz, &DIsDouble);
		#endif

		I   = TECDAT142(&IMax, ax, &DIsDouble);
		I   = TECDAT142(&IMax, ay, &DIsDouble);
		#if SIMDIM == 3
			double* az = new double[size];
			#pragma omp parallel for
	  		for(uint ii = start; ii < end; ++ii)
				az[ii-start] = pnp1[ii].v(2);
			I   = TECDAT142(&IMax, az, &DIsDouble);
		#endif

	    if (svar.outform == 3)
		{
			double* cID = new double[size];
		  	double* cVx = new double[size];
		  	double* cVy = new double[size];
		  	double* cP = new double[size];
		  	double* cRho = new double[size];
			#pragma omp parallel for
		  	for(uint ii = start; ii < end; ++ii)
		  	{
				cID[ii-start] = pnp1[ii].cellID;
		  		cVx[ii-start] = pnp1[ii].cellV(0);
		  		cVy[ii-start] = pnp1[ii].cellV(1);
		  		cP[ii-start] = pnp1[ii].cellP;
		  		cRho[ii-start] = pnp1[ii].cellRho;
	  		}

		    
		    I   = TECDAT142(&IMax, cVx, &DIsDouble);
		    I   = TECDAT142(&IMax, cVy, &DIsDouble);

		  	#if SIMDIM == 3
		  		double* cVz = new double[size];
		  		#pragma omp parallel for
		  		for(uint ii = start; ii < end; ++ii)
			  		cVz[ii-start] = pnp1[ii].cellV(2);
			  	I   = TECDAT142(&IMax, cVz, &DIsDouble);
		  	#endif

		    I   = TECDAT142(&IMax, cRho, &DIsDouble);
		    I   = TECDAT142(&IMax, cP, &DIsDouble);
		    I   = TECDAT142(&IMax, cID, &DIsDouble);
	    }
	    
	}
	else if (svar.outform == 4)
		{
			// variables = "x y z rho P m v a b Neighbours Aero";

			double* v = new double[size];
		  	double* a = new double[size];
			float* b = new float[size];
			float* nNb = new float[size];
			double* aF = new double[size];

			#pragma omp parallel for
		  	for(uint ii = start; ii < end; ++ii)
		  	{
		  		v[ii-start] = pnp1[ii].v.norm();
		  		a[ii-start] = pnp1[ii].f.norm();
				b[ii-start] = pnp1[ii].b;
		  		nNb[ii-start] = pnp1[ii].theta;
		  		aF[ii-start] = pnp1[ii].Af.norm();
	  		}

	  		INTEGER4 IsFloat = 0;
			I   = TECDAT142(&IMax, v, &DIsDouble);
		    I   = TECDAT142(&IMax, a, &DIsDouble);
	  		I   = TECDAT142(&IMax, b, &IsFloat);
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
	INTEGER4 VIsDouble       = 0; 
	if(sizeof(ldouble) == 8)
		VIsDouble            = 1;
    
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
    	std::string variables = "x y";	
		if (svar.outform == 1)
		{
			variables = "x y rho P m v a";
		}
		else if (svar.outform == 2)
		{
			variables = "x y rho P m v_x v_y a_x a_y";
		}
		else if (svar.outform == 3)
		{
			variables = "x y rho P m v_x v_y a_x a_y Cell_Vx Cell_Vy Cell_Rho Cell_P Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "x y rho P m v a b Neighbours Aero";
		}

	#endif

	#if SIMDIM == 3
		std::string variables = "x y z";  
		if (svar.outform == 1)
		{
			variables = "x y z rho P m v a";
		}
		else if (svar.outform == 2)
		{
			variables = "x y z rho P m v_x v_y v_z a_x a_y a_z";
		}
		else if (svar.outform == 3)
		{
			variables = 
		"x y z rho P m v_x v_y v_z a_x a_y a_z Cell_Vx Cell_Vy Cell_Vz Cell_Rho Cell_P Cell_ID";
		}
		else if (svar.outform == 4)
		{
			variables = "x y z rho P m v a b Neighbours Aero";
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
    	exit(-1);
}

#endif