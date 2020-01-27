#ifndef CDFIO_H
#define CDFIO_H

#include "Var.h"
#include "Crossing.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>

// #include <H5hut.h>
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2

using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string;

struct ZONE
{
	ZONE()
	{
		lineNo = 0;
	}

	string name, ETtype;
	uint lineNo;
	uint ctype;
	uint nF, nCverts, nFverts, nFverts2;
	uint nP, nE;
	uint pressOrcp, nvar;
	uint veltype, velstart, cpstart, densstart;
};

uint index(uint ii, uint jj, uint nPts)
{
	return(ii*nPts + jj);
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(uint ii=0; ii < num - 1; ++ii){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

/*FOR DEBUGGING*/
void Write_Zone(string input, ZONE& zn, MESH& cells)
{
	std::size_t found1 = input.find_first_of("/\\");
	std::size_t found2 = input.find_last_of("/\\");
	string path = input.substr(found1,found2-found1+1);
	
	string filename="Outputs";
	filename.append(path);
	filename.append(zn.name);
	filename.append(".plt");
	
	ofstream fout(filename,std::ios::out);
	if (!fout.is_open())
	{
		cout << "Failed to open data file for writing mesh." << endl;
		exit(-1);
	}

	fout << "TITLE = \"3D TAU Solution\"\n";
	fout << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\" \"x_velocity\" \"y_velocity\" \"z_velocity\"\n";
	fout << "ZONE T=\"" << zn.name << "\"" << endl;
	fout << " N=" << zn.nP << ", E=" << zn.nE << 
	", F=FEBLOCK, ET=" 	<< zn.ETtype << endl << endl;

	fout << std::left << std::scientific << std::setprecision(6);
	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	uint kk = 0;
		for(uint jj = cells.verts.size() - zn.nP; jj < cells.verts.size(); ++jj)
		{
			fout << cells.verts[jj][ii] << std::setw(15);
			kk++;

			if(kk == 5)
			{
				fout << "\n";
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fout << "\n";
	}

	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	uint kk = 0;
		for(uint jj = cells.pVel.size() - zn.nP; jj < cells.pVel.size(); ++jj)
		{
			fout << cells.pVel[jj][ii] << std::setw(15);
			kk++;

			if(kk == 5)
			{
				fout << "\n";
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fout << "\n";
	}

	// cout << "Velocitices written" << endl;
	fout << std::fixed;
	if (zn.ctype == 1 || zn.ctype == 2)
	{
		for(uint ii = cells.elems.size()-zn.nE; ii < cells.elems.size(); ++ii)
		{	
			for(auto elem:cells.elems[ii])
			{
				fout << elem+1 << std::setw(15);
			}
			fout << "\n";
		}
	}


	if (zn.ctype == 3)
	{
		for(uint ii = cells.elems.size()-zn.nE; ii < cells.elems.size(); ++ii)
		{	
			for(uint jj = 0; jj < 3; ++jj)
			{
				fout << cells.elems[ii][jj]+1 << " ";
			}
			fout << cells.elems[ii][2]+1 << " ";
			for(uint jj = 3; jj < 6; ++jj)
			{
				fout << cells.elems[ii][jj]+1 << " ";
			}
			fout << cells.elems[ii][5]+1;
			fout << "\n";
		}
	}

	if (zn.ctype == 4)
	{
		for(uint ii = cells.elems.size()-zn.nE; ii < cells.elems.size(); ++ii)
		{	
			for(uint jj = 0; jj < 4; ++jj)
			{
				fout << cells.elems[ii][jj]+1 << " ";
			}
			fout << cells.elems[ii][4]+1 << " ";
			fout << cells.elems[ii][4]+1 << " ";
			fout << cells.elems[ii][4]+1 << " ";
			fout << cells.elems[ii][4]+1 << " ";
			fout << "\n";
		}
	}
	fout.close();
}

template <class T>
void Average_Point_to_Cell(const vector<T>& pData, vector<T>& cData,
							const vector<vector<uint>>& elems, const T zero)
{
	vector<T> sum(elems.size(),zero);

	#pragma omp parallel for reduction(+:sum)
	for(uint ii = 0; ii < elems.size(); ++ii)
	{
		
		const uint nVerts = elems[ii].size();
		for (auto jj:elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum; 
}

vector<ldouble> CpToPressure(const vector<ldouble>& Cp, const FLUID& fvar)
{
	vector<ldouble> press(Cp.size());
	#pragma omp parallel for shared(Cp)
	for (uint ii = 0; ii < Cp.size(); ++ii)
	{
		press[ii] = Cp[ii]*fvar.gasDynamic /*+ fvar.gasPress*/;
	}
	return press;
}

/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTIONS ********************/
/*****************************************************************************/

/*To run on the mesh file*/
vector<vector<uint>> Get_Element(NcFile& fin, string variable)
{
	#ifdef DEBUG
		dbout << "Reading Element: " << variable << endl;
	#endif
	NcVar elemData = fin.getVar(variable);	
	if(!elemData.isNull())
	{
		NcDim dim = elemData.getDim(0);
		size_t nElem = dim.getSize();

		dim = elemData.getDim(1);
		size_t nPoints = dim.getSize();
		// cout << nElem << "  " << nPoints << endl;
		#ifdef DEBUG
		dbout << "Allocating array of: " << nElem << " by " << nPoints << endl;
		#endif

		/*Allocate on the heap (can be big datasets)*/		
		int* elemArray = new int[nElem*nPoints];
				
		/*Get the actual data from the file*/
		vector<size_t> startp,countp;
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(nElem);
		countp.push_back(nPoints);
		
		#ifdef DEBUG
		dbout << "Attempting to read NetCDF elements." << endl;
		#endif

		elemData.getVar(startp,countp,elemArray);

		cout << "Successfully read: " << variable << endl;
		cout << "Number of cells: " << nElem << endl;

		#ifdef DEBUG
		dbout << "Successfully read elements" << endl;
		#endif

		/*Convert it to a vector to store*/
		uint ii,jj;
		vector<vector<uint>> elemVec(nElem,vector<uint>(nPoints));
		for (ii = 0; ii < nElem; ++ii)
		{
			for(jj = 0; jj < nPoints; ++jj)
				elemVec[ii][jj] = static_cast<uint>(elemArray[index(ii,jj,nPoints)]);
		}

		#ifdef DEBUG
		dbout << "Returning vector" << endl;
		#endif
		return elemVec;
	}
	else return vector<vector<uint>>(0);
}

/*To run on the mesh file*/
vector<StateVecD> Get_Coordinates(NcFile& fin)
{
	#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
	#endif

	NcVar coordXD = fin.getVar("points_xc");	
	if(coordXD.isNull()) exit(NC_ERR);

	NcVar coordYD = fin.getVar("points_yc");	
	if(coordYD.isNull()) exit(NC_ERR);

	NcVar coordZD = fin.getVar("points_zc");	
	if(coordZD.isNull()) exit(NC_ERR);
	
	NcDim dim = coordXD.getDim(0);
	size_t nPts = dim.getSize();

	#ifdef DEBUG
		dbout << "Number of points: " << nPts << endl;
	#endif
	/*Allocate on the heap (can be big datasets)*/
	double* coordX = new double[nPts];
	double* coordY = new double[nPts];
	double* coordZ = new double[nPts];

	/*Get the actual data from the file*/
	coordXD.getVar(coordX);
	coordYD.getVar(coordY);
	coordZD.getVar(coordZ);

	/*Convert it to a vector to store*/
	vector<StateVecD> coordVec(nPts);
#if SIMDIM == 3
	for (uint ii = 0; ii < nPts; ++ii)
	{
		coordVec[ii] = StateVecD(coordX[ii],coordY[ii],coordZ[ii]);
	}
#else
	/*Need to find which coordinate to ignore*/
#ifdef DEBUG
	dbout << "Checking which dimension to ignore" << endl;
#endif
	uint counts[3] = {0};
	for(uint ii = 0; ii < nPts; ++ii)
	{
		if(coordX[ii] == 0 || coordX[ii] == -1)
		{
			counts[0]++;
		}

		if(coordY[ii] == 0 || coordY[ii] == -1)
		{
			counts[1]++;
		}

		if(coordZ[ii] == 0 || coordZ[ii] == -1)
		{
			counts[2]++;
		}
	}
	uint ignore=0;
	if (counts[0] == nPts)
		ignore = 1;
	else if (counts[1] == nPts)
		ignore = 2;
	else if (counts[2] == nPts)
		ignore = 3;
	else
		cout << "Couldn't determine which dimension is false." << endl;

#ifdef DEBUG
	dbout << "Ignored dimension: " << ignore << endl;
#endif
	cout << "Ignored dimension: " << ignore << endl;

	for (uint ii = 0; ii < nPts; ++ii)
	{
		if(ignore == 1)
			coordVec[ii] = StateVecD(coordY[ii],coordZ[ii]);
		else if(ignore == 2)
			coordVec[ii] = StateVecD(coordX[ii],coordZ[ii]);
		else if(ignore == 3)
			coordVec[ii] = StateVecD(coordX[ii],coordY[ii]);
	}
#endif
#ifdef DEBUG
		dbout << "Returning coordinates." << endl;
#endif
	return coordVec;	
}

/*To run on the solution file*/
template <class T>
int Get_Scalar_Property(NcFile& fin, string variable, vector<T>& var)
{
	#ifdef DEBUG
		dbout << "Reading variable: " << variable << endl;
	#endif
	NcVar propData = fin.getVar(variable);	
	if(!propData.isNull())
	{
		NcDim dim = propData.getDim(0);
		size_t nPts = dim.getSize();
		#ifdef DEBUG
		dbout << "Allocating array of: " << nPts << endl;
		#endif
		T* array = new T[nPts];
		
		#ifdef DEBUG
		dbout << "Attempting to read NetCDF variable:  " << variable << endl;
		#endif

		/*Get the actual data from the file*/
		propData.getVar(array);

		/*Convert it to a vector to store*/
		vector<T> propVec;
		propVec.insert(propVec.end(),&array[0],&array[nPts]);
		var = propVec;

		#ifdef DEBUG
		dbout << "returning vector" << endl;
		#endif
		
		return 0;
	}
	else return 1;
}

void Get_Cell_Faces(const vector<vector<uint>>& cell,
	const vector<vector<uint>>& facenum, std::vector<std::vector<std::vector<uint>>>& cFaces)
{
	for(uint ii = 0; ii < cell.size(); ++ii)
	{	
		uint jj = 0; 
		cFaces.emplace_back();
		for(auto faces:facenum)
		{	
			cFaces[ii].emplace_back();
			for (auto vert:faces)
			{
				cFaces[ii][jj].emplace_back(cell[ii][vert]);
			}
			++jj;
		}
	}
}

/*****************************************************************************/
/***************** READING NETCDF SOLUTION DATA FUNCTIONS ********************/
/*****************************************************************************/

void Read_SOLUTION(const string& solIn, const FLUID& fvar, MESH& cells)
{
	NcFile sol(solIn, NcFile::read);

	NcDim solPN = sol.getDim("no_of_points");
	uint solPts = static_cast<uint>(solPN.getSize());

	#ifdef DEBUG
	dbout << "Solution points: " << solPts << endl;
	#endif

	if(solPts!=cells.numPoint)
	{
		cout << "Solution file does not have the same number of vertices as the mesh." << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}

	/*Get the velocities*/
	vector<double> xvel, yvel, zvel;
	Get_Scalar_Property(sol, "x_velocity", xvel);
	#if SIMDIM == 3
	Get_Scalar_Property(sol, "y_velocity", yvel);
	#endif
	Get_Scalar_Property(sol, "z_velocity", zvel);

	/*Test for size*/
	if(xvel.size() == solPts)
	{	/*Turn the arrays into a state vector*/
		vector<StateVecD> vel(solPts);
		
		for(uint ii = 0; ii < solPts; ++ii)
		{
			#if SIMDIM == 3
			vel[ii] = StateVecD(xvel[ii],yvel[ii],zvel[ii]);
			#else
			vel[ii] = StateVecD(xvel[ii],zvel[ii]);
			#endif
		}
		cells.pVel = vel;
	}
	else 
	{
		cout << "velocities do not have the same number of vertices as the mesh." << endl;
		cout << xvel.size() << "  " << solPts << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}

	/*Get other scalar data that may or may not exist*/
	vector<double> press;
	Get_Scalar_Property(sol,"pressure",press);
	if(press.size()==0)
	{
		/*Pressure data doesn't exist, so go to cp*/
		vector<double> cp;
		Get_Scalar_Property(sol,"cp",cp);

		/*Then convert cp to pressure. can get global attributes from sol file.. Maybe...*/
		press = CpToPressure(cp, fvar);
	}
	else
	{
		for(uint ii = 0; ii < press.size(); ++ii)
			press[ii] -= fvar.gasPress;
	}

	vector<double> dens(solPts);

	for(uint ii = 0; ii < press.size(); ++ii)
	{
		dens[ii] = fvar.rho0 * pow((press[ii]/fvar.B + 1),1/fvar.gam);
	}	

	if(press.size() == 0)
	{
		cout << "The solution data has not been correctly ingested. Please check the solution file." << endl;
	}

	/*Update the point arrays with the actual stuff*/
	cells.pointP = press;
	cells.pointRho = dens;	

	
	/*Average the data to a the cell*/
	cells.SetCells();
	StateVecD zero = StateVecD::Zero();
	Average_Point_to_Cell(cells.pVel,cells.cVel, cells.elems, zero);
	Average_Point_to_Cell(cells.pointRho,cells.cellRho,cells.elems,0.0);
	Average_Point_to_Cell(cells.pointP,cells.cellP,cells.elems,0.0);
	
	/*Find cell centres*/
	Average_Point_to_Cell(cells.verts,cells.cCentre,cells.elems,zero);

	/*Check the data*/
	// for(auto value:dens)
	// {
	// 	if(value == 0)
	// 	{
	// 		cout << "Cell Rho data has not been initialised" << endl;
	// 	}
	// }
}


/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTION *********************/
/*****************************************************************************/
#if SIMDIM == 3
void Read_TAUMESH(const SIM& svar, MESH& cells, const FLUID& fvar)
{
	
	string meshIn = svar.infolder;
	string solIn = svar.infolder;

	meshIn.append(svar.meshfile);
	solIn.append(svar.solfile);
	
	#ifdef DEBUG 
		dbout << "Attempting read of NetCDF file." << endl;
		dbout << "Mesh file: " << meshIn << endl;
		dbout << "Solution file: " << solIn << endl;
	#endif
	/*Read the mesh data*/
	NcFile fin(meshIn, NcFile::read);
	cout << "Mesh file open. Reading cell data..." << endl;

	// Retrieve how many elements there are.
	NcDim elemNo = fin.getDim("no_of_elements");
	NcDim pointNo = fin.getDim("no_of_points");
	uint nElem = static_cast<uint>(elemNo.getSize());
	uint nPts = static_cast<uint>(pointNo.getSize());

	#ifdef DEBUG
		dbout << "nElem : " << nElem << " nPts: " << nPts << endl;
	#endif

	cout << "nElem : " << nElem << " nPts: " << nPts << endl;
	cells.numPoint = nPts;
	cells.numElem = nElem;

	/*Retrieve the cells*/
	vector<vector<uint>> tets = Get_Element(fin,"points_of_tetraeders");
	vector<vector<uint>> prism = Get_Element(fin,"points_of_prisms");
	vector<vector<uint>> pyra = Get_Element(fin,"points_of_pyramids");
	vector<vector<uint>> hex = Get_Element(fin,"points_of_hexaeders");
	
	#ifdef DEBUG
		if(tets.size() == 0)
			dbout << "No tetraeders in mesh file." << endl;

		if(prism.size() == 0)
			dbout << "No prisms in mesh file." << endl;

		if(pyra.size() == 0)
			dbout << "No pyramids in mesh file." << endl;

		if(hex.size() == 0)
			dbout << "No hexaeders in mesh file." << endl;
	#endif


	/*Get the coordinates of the mesh*/
	cells.verts = Get_Coordinates(fin);
	if(cells.verts.size()!= nPts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}
	

	#ifdef DEBUG
	dbout << "All element data ingested" << endl;
	#endif


	/*Put data into cell structure, and generate face based data*/
	if(tets.size()!=0)
	{
		cells.elems.insert(cells.elems.end(),tets.begin(),tets.end());
		#if SIMDIM == 3
		cout << "Buiding cell faces for points_of_tetraeders" << endl;
		std::vector<std::vector<uint>> facenum;
		vector<vector<vector<uint>>> cFaces;
		#ifdef DEBUG
		dbout << "Getting tetraeder cell faces." << endl;
		#endif
		facenum = {{0,1,2},{0,2,3},{0,3,1},{1,3,2}};
		Get_Cell_Faces(tets,facenum,cFaces);
		
		cells.cFaces.insert(cells.cFaces.end(), cFaces.begin(), cFaces.end());
		#endif
	}

	if(prism.size()!=0)
	{
		cells.elems.insert(cells.elems.end(),prism.begin(),prism.end());
		#if SIMDIM == 3
		cout << "Buiding cell faces for points_of_prisms" << endl;
		std::vector<std::vector<uint>> facenum;
		vector<vector<vector<uint>>> cFaces;
		#ifdef DEBUG
		dbout << "Getting prism cell faces." << endl;
		#endif
		facenum = {{1,4,5},{1,5,2},{0,3,5},{0,5,2},
				   {0,3,4},{0,4,1},{0,1,2},{5,4,3}};
		Get_Cell_Faces(prism,facenum,cFaces);
		
		cells.cFaces.insert(cells.cFaces.end(), cFaces.begin(), cFaces.end());
		#endif
	}

	if(pyra.size()!=0)
	{
		cells.elems.insert(cells.elems.end(),pyra.begin(),pyra.end());
		#if SIMDIM == 3
		cout << "Buiding cell faces for points_of_pyramids" << endl;
		std::vector<std::vector<uint>> facenum;
		vector<vector<vector<uint>>> cFaces;
		#ifdef DEBUG
		dbout << "Getting pyramid cell faces." << endl;
		#endif
		facenum = {{0,2,1},{0,3,2},{0,1,4},{0,4,3},{1,4,2},{2,3,4}};
		Get_Cell_Faces(pyra,facenum,cFaces);
		
		cells.cFaces.insert(cells.cFaces.end(), cFaces.begin(), cFaces.end());
		#endif
	}

	if(hex.size()!=0)
	{
		cells.elems.insert(cells.elems.end(),hex.begin(),hex.end());
		#if SIMDIM == 3
		cout << "Buiding cell faces for points_of_hexaeders" << endl;
		std::vector<std::vector<uint>> facenum;
		vector<vector<vector<uint>>> cFaces;
		#ifdef DEBUG
		dbout << "Getting hexaeder cell faces." << endl;
		#endif
		facenum = {{1,5,2},{6,2,5},{2,6,7},{7,3,2},{0,3,4},{3,7,4},
				   {0,4,1},{5,1,4},{0,1,2},{0,2,3},{7,6,4},{6,5,4}};
		Get_Cell_Faces(hex,facenum,cFaces);
		
		cells.cFaces.insert(cells.cFaces.end(), cFaces.begin(), cFaces.end());
		#endif
	}


	if(cells.elems.size()!=nElem)
	{
		cout << "Some data has been missed.\nPlease check which kinds of volumes used." << endl;
	}

	#ifdef DEBUG
	dbout << "End of interaction with mesh file and ingested data." << endl << endl;
	dbout << "Opening solultion file." << endl;
	#endif
	// /*End of reading the mesh file.*/
	// cout << "Trying solution file: " << solIn << endl;
	Read_SOLUTION(solIn,fvar,cells);
}
#endif

/*****************************************************************************/
/******************** READING NETCDF 2D DATA FUNCTION ************************/
/*****************************************************************************/
#if SIMDIM == 2


int Find_Bmap_Markers(const string& bmapIn)
{
	#ifdef DEBUG
	dbout << "Entering Find_Bmap_Markers..." << endl;
	#endif

	std::ifstream fin(bmapIn,std::ios::in);

	if(!fin.is_open())
	{
		cout << "Couldn't open the boundary map file." << endl;
		cout << "Attempted path: " << bmapIn << endl;
		exit(-1);
	}

	int marker = 0;
	/*Search for the markers that have either farfield or symmetry plane*/
	uint lineno = 0;
	string line;

	uint blockno = 0;
	uint blockstart = 1; /*A store of when this block starts to search through*/

	while(getline(fin,line))
	{

		// cout << line << endl;
		if(line.find("Type: symmetry plane")!=string::npos)
		{
			
			/*This marker is a far field, so store it.*/
			/*Go to the start of the block, and search for the marker ID*/
			GotoLine(fin,blockstart);
			while(getline(fin,line))
			{	
				// cout << "inner:\t" << line << endl;
				if(line.find("Markers:")!=string::npos)
				{
					// cout << "Found a boundary marker" << endl;
					std::stringstream sstr;

					sstr << line;

					string temp;
					int found;

					while(!sstr.eof())
					{
						sstr >> temp;
						if(std::stringstream(temp) >> found)
						{
							marker = found;
						}

						temp = "";
					}

					goto markerfound;
				}
			}
		}

		if(line.find("block end")!= string::npos)
		{	/*We're on a new block*/
			blockno++;
			blockstart = lineno+1;
		}

		lineno++;
	}

markerfound:
	cout << "Symmetry plane marker:" << marker << endl;

	#ifdef DEBUG
	dbout << "Exiting Find_Bmap_Markers..." << endl;
	#endif

	return marker;
}

vector<vector<uint>> Get_Symmetry_Plane(NcFile& fin, const int marker)
{

#ifdef DEBUG
	dbout << "Entering Get_Symmetry_Plane..."<< endl;
	dbout << "Opening variable: boundarymarker_of_surfaces" << endl;
#endif

	NcVar bmarkers = fin.getVar("boundarymarker_of_surfaces");
	if(bmarkers.isNull()) exit(NC_ERR);

	NcDim dim = bmarkers.getDim(0);
	size_t nMarkers = dim.getSize();
	
#ifdef DEBUG
		dbout << "Number of surface markers: " << nMarkers << endl;
#endif

	int* markers = new int[nMarkers];

	/*Get the actual data from the file*/
	bmarkers.getVar(markers);


	cout << "Successfully read boundary markers." << endl;

	#ifdef DEBUG
	dbout << "Successfully read boundary markers." << endl;
	#endif

	/*Now compare against the marker for the symmetry plane, */ 
	/*and add index to a vector to retrieve the faces*/
	vector<uint> symPlaneM;

	for(uint ii = 0; ii < nMarkers; ++ii)
	{
		if(markers[ii] == marker)
			symPlaneM.emplace_back(ii);
	}

	/*Get the surface indexes*/
	vector<vector<uint>> temp;
	uint nQuads = 0;
	uint nTris = 0;
	/*Start with surface quadrilaterals*/
#ifdef DEBUG
	dbout << "Attempting to read surface quadrilaterals." << endl;
#endif

	NcVar quadVar = fin.getVar("points_of_surfacequadrilaterals");
	if(!quadVar.isNull()) 
	{
		dim = quadVar.getDim(0);
		NcDim dim2 = quadVar.getDim(1);
		nQuads = dim.getSize();
		size_t nPPQuad = dim2.getSize();
			
#ifdef DEBUG
		dbout << "Number of surface quadrilaterals: " << nQuads << endl;
#endif

		int* quads = new int[nQuads*nPPQuad];

		/*Get the actual data from the file*/
		quadVar.getVar(quads);

		cout << "Successfully read surface quadrilaterals." << endl;

#ifdef DEBUG
		dbout << "Successfully read surface quadrilaterals." << endl;
#endif
		/*Create the vector of surfaces*/
		for (uint ii = 0; ii < nQuads; ++ii)
		{
			temp.emplace_back();
			for(uint jj = 0; jj < nPPQuad; ++jj)
				temp[ii].emplace_back(static_cast<uint>(quads[index(ii,jj,nPPQuad)]));
		}
	}
	else 
	{

#ifdef DEBUG
		dbout << "No surface quadrilaterals." << endl;
#endif
	}

#ifdef DEBUG
	dbout << "Attempting to read surface triangles." << endl;
#endif

	NcVar triVar = fin.getVar("points_of_surfacetriangles");
	if(!triVar.isNull())
	{
		dim = triVar.getDim(0);
		NcDim dim2 = triVar.getDim(1);
		nTris = dim.getSize();
		size_t nPPTri = dim2.getSize();
			
#ifdef DEBUG
			dbout << "Number of surface triangles: " << nTris << endl;
#endif

		int* tris = new int[nTris*nPPTri];

		/*Get the actual data from the file*/
		triVar.getVar(tris);

		cout << "Successfully read surface triangles." << endl;

#ifdef DEBUG
		dbout << "Successfully read surface triangles." << endl;
#endif

		/*Create the vector of surfaces*/
		for (uint ii = nQuads; ii < nQuads + nTris; ++ii)
		{
			temp.emplace_back();
			for(uint jj = 0; jj < nPPTri; ++jj)
				temp[ii].emplace_back(static_cast<uint>(tris[index(ii,jj,nPPTri)]));
		}
	}
	else 
	{
#ifdef DEBUG
		dbout << "No surface triangles." << endl;
#endif
	}

	if(temp.size()!= nMarkers)
	{
		cout << "Some surface faces have been missed." << endl;
		exit(-1);
	}

	vector<vector<uint>> symPlane;
	for(auto const& index:symPlaneM)
	{
		symPlane.emplace_back(temp[index]);
	}

	return symPlane;
}

void Read_TAUMESH(const SIM& svar, MESH& cells, const FLUID& fvar)
{
	string bmapIn = svar.infolder;
	string meshIn = svar.infolder;
	string solIn = svar.infolder;

	bmapIn.append(svar.bmapfile);
	meshIn.append(svar.meshfile);
	solIn.append(svar.solfile);
	
	#ifdef DEBUG 
		dbout << "Attempting read of NetCDF file." << endl;
		dbout << "Bmap file: " << bmapIn << endl;
		dbout << "Mesh file: " << meshIn << endl;	
		dbout << "Solution file: " << solIn << endl;
	#endif

	/*Find which marker we need to read for the surfaces*/
	int marker = Find_Bmap_Markers(bmapIn);


	/*Read the mesh data*/
	NcFile fin(meshIn, NcFile::read);
	cout << "Mesh file open. Reading cell data..." << endl;

	// Retrieve how many elements there are.
	NcDim elemNo = fin.getDim("no_of_elements");
	NcDim pointNo = fin.getDim("no_of_points");
	NcDim nSurfaces = fin.getDim("no_of_surfaceelements");
	uint nElem = static_cast<uint>(elemNo.getSize());
	uint nPts = static_cast<uint>(pointNo.getSize());

	#ifdef DEBUG
		dbout << "nElem : " << nElem << " nPts: " << nPts << endl;
	#endif

	cells.numPoint = nPts;
	cells.numElem = nElem;

	/*Retrieve the symmetry plane faces*/
	cells.elems = Get_Symmetry_Plane(fin, marker);


	/*Get the coordinates of the mesh*/
	cells.verts = Get_Coordinates(fin);
	if(cells.verts.size()!= nPts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}
	

	#ifdef DEBUG
	dbout << "All element data ingested" << endl;
	dbout << "End of interaction with mesh file and ingested data." << endl << endl;
	dbout << "Opening solultion file." << endl;
	#endif
	// /*End of reading the mesh file.*/
	// cout << "Trying solution file: " << solIn << endl;
	Read_SOLUTION(solIn,fvar,cells);
}
#endif


/*****************************************************************************/
/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
vector<vector<uint>> Get_Faces(NcFile& fin)
{
	#ifdef DEBUG
		dbout << "Reading face data" << endl;
	#endif

	NcVar faceData = fin.getVar("points_of_element_faces");	
	if(faceData.isNull())
	{	/*No data is available*/
		return vector<vector<uint>>(0);
	}

	NcDim dim = faceData.getDim(0);
	size_t nFace = dim.getSize();

	dim = faceData.getDim(1);
	size_t nPnts = dim.getSize();
	// cout << nElem << "  " << nPoints << endl;
	#ifdef DEBUG
	dbout << "Allocating array of: " << nFace << " by " << nPnts << endl;
	#endif

	/*Allocate on the heap (can be big datasets)*/		
	int* faceArray = new int[nFace*nPnts];
			
	/*Get the actual data from the file*/
	vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(nFace);
	countp.push_back(nPnts);
	
	#ifdef DEBUG
	dbout << "Attempting to read NetCDF faces." << endl;
	#endif

	faceData.getVar(startp,countp,faceArray);

	cout << "Successfully read face data." << endl;

	#ifdef DEBUG
	dbout << "Successfully read face data." << endl;
	#endif

	/*Convert it to a vector to store*/
	uint ii,jj;
	vector<vector<uint>> faceVec(nFace,vector<uint>(nPnts));
	for (ii = 0; ii < nFace; ++ii)
	{
		for(jj = 0; jj < nPnts; ++jj)
			faceVec[ii][jj] = static_cast<uint>(faceArray[index(ii,jj,nPnts)]);
	}

	#ifdef DEBUG
	dbout << "Returning vector" << endl;
	#endif
	return faceVec;
	
}

void Place_Faces(NcFile& fin, const vector<vector<uint>> faces, MESH& cells)
{
	#ifdef DEBUG
		dbout << "Reading element left/right and placing faces" << endl;
	#endif

	/*Read how many faces there should be per element/cell*/
	// NcVar nFaceElems = fin.getVar("faces_per_element");
	// if(nFaceElems.isNull())
	// {	/*No data is available*/
	// 	cout << "No data on how many faces per cell. Stopping" << endl;
	// 	exit(-1);
	// }

	// NcDim dim = nFaceElems.getDim(0);
	// size_t nElems = dim.getSize();

	// uint* elemFaceNo = new uint[nElems];
	
	// #ifdef DEBUG
	// dbout << "Attempting to read NetCDF element face count." << endl;
	// #endif

	// nFaceElems.getVar(elemFaceNo);

	// cout << "Successfully read element face count data." << endl;

	// #ifdef DEBUG
	// dbout << "Successfully read element face count data." << endl;
	// #endif

	NcVar leftData = fin.getVar("left_element_of_faces");	
	if(leftData.isNull())
	{	/*No data is available*/
		cout << "No cell left data. Stopping" << endl;
		exit(-1);
	}

	NcDim dim = leftData.getDim(0);
	size_t nLeft = dim.getSize();

	int* left = new int[nLeft];
	
	#ifdef DEBUG
	dbout << "Attempting to read NetCDF left elements." << endl;
	#endif

	leftData.getVar(left);

	cout << "Successfully read left element data." << endl;

	#ifdef DEBUG
	dbout << "Successfully read left element data." << endl;
	#endif

	NcVar rightData = fin.getVar("right_element_of_faces");	
	if(leftData.isNull())
	{	/*No data is available*/
		cout << "No cell right data. Stopping" << endl;
		exit(-1);
	}

	dim = rightData.getDim(0);
	size_t nRight = dim.getSize();

	int* right = new int[nRight];
	
	#ifdef DEBUG
	dbout << "Attempting to read NetCDF right elements." << endl;
	#endif

	rightData.getVar(right);

	cout << "Successfully read right element data." << endl;

	#ifdef DEBUG
	dbout << "Successfully read right element data." << endl;
	#endif

	// vector<vector<vector<uint>>> cFaces(nLeft);

	for(uint ii = 0; ii < nLeft; ++ii)
	{
		int lindex = left[ii];
		int rindex = right[ii];

		cells.cFaces[lindex].emplace_back(faces[ii]);

		if(rindex == -2)
		{	/*Farfield face*/
			cells.farfield.emplace_back(faces[ii]);
		}
		else if (rindex == -1)
		{	/*Surface face*/
			cells.surface.emplace_back(faces[ii]);
		}
		else
		{
			cells.cFaces[rindex].emplace_back(faces[ii]);
		}
	}

	/*Check if each element has the correct number of faces*/
	// for(uint ii = 0; ii < nElems; ++ii)
	// {
	// 	if(cells.cFaces[ii].size() != elemFaceNo[ii])
	// 	{
	// 		cout << "Mismatch of number of faces in an element. " <<
	// 			"Element " << ii <<  " should have " << elemFaceNo[ii] << "  but has " 
	// 			<< cells.cFaces[ii].size() << " faces." << endl;
	// 		exit(-1);
	// 	}
	// }

	#ifdef DEBUG
		dbout << "End of placing faces in elements." << endl;
	#endif

	/*Now go through the faces and see which vertices are unique, to get element data*/
	cout << "Building elements..." << endl;
#pragma omp parallel
{
	// vector<vector<uint>> local;
	#pragma omp for schedule(static) nowait
	for(uint ii = 0; ii < cells.cFaces.size(); ++ii)
	{
		for(auto const& faces:cells.cFaces[ii])
		{
			for(auto const& vert:faces)
			{
				if(std::find(cells.elems[ii].begin(),cells.elems[ii].end(),vert)
					==cells.elems[ii].end())
				{	/*Vertex doesn't exist in the elems vector yet.*/
					// #pragma omp critical
					cells.elems[ii].emplace_back(vert);
				}
			}
		}
	}
}
	cout << "Elements built." << endl;

	#ifdef DEBUG
	dbout << "All elements defined." << endl << endl;
	#endif
	

}

void Read_TAUMESH_FACE(SIM& svar, MESH& cells, const FLUID& fvar)
{
	string meshIn = svar.infolder;
	string solIn = svar.infolder;
	meshIn.append(svar.meshfile);
	solIn.append(svar.solfile);
	
	#ifdef DEBUG 
		dbout << "Attempting read of NetCDF file." << endl;
		dbout << "Mesh file: " << meshIn << endl;
		dbout << "Solution file: " << solIn << endl;
	#endif
	/*Read the mesh data*/
	NcFile fin(meshIn, NcFile::read);
	cout << "Mesh file open. Reading face data..." << endl;

	// Retrieve how many elements there are.
	NcDim elemNo = fin.getDim("no_of_elements");
	NcDim pointNo = fin.getDim("no_of_points");
	NcDim faceNo = fin.getDim("no_of_faces");
	uint nElem = static_cast<uint>(elemNo.getSize());
	uint nPnts = static_cast<uint>(pointNo.getSize());
	uint nFace = static_cast<uint>(faceNo.getSize());

	#ifdef DEBUG
		dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
	#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
	
	cells.numPoint = nPnts;
	cells.numElem = nElem;

	cells.elems = vector<vector<uint>>(nElem);
	cells.cFaces = vector<vector<vector<uint>>>(nElem);
	cells.verts = vector<StateVecD>(nPnts);

	vector<vector<uint>> faces = Get_Faces(fin);

	/*Get the coordinates of the mesh*/
	cells.verts = Get_Coordinates(fin);
	if(cells.verts.size()!= nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for(auto& vert:cells.verts)
	{
		vert/=cells.scale;
	}
	svar.Start/=cells.scale;
	// svar.Jet/=cells.scale;

	/*Get face left and right, and put the faces in the elements*/
	Place_Faces(fin, faces, cells);

	#ifdef DEBUG
	dbout << "End of interaction with mesh file and ingested data." << endl << endl;
	dbout << "Opening solultion file." << endl;
	#endif

	Read_SOLUTION(solIn,fvar,cells);
}


/*****************************************************************************/
/****************** WRITING NETCDF DATA FUNCTIONS ****************************/
/*****************************************************************************/

void Write_Group(NcGroup &zone, const State &pnp1, 
	const uint start, const uint end, const uint outform)
{
	/*Write the Particles group*/
	
	const uint size = end - start;
	NcDim posDim = zone.addDim("Point",size);

	NcVar xVar = zone.addVar("x",ncDouble,posDim);
	NcVar yVar = zone.addVar("y",ncDouble,posDim);
	NcVar zVar = zone.addVar("z",ncDouble,posDim);
	NcVar idVar = zone.addVar("id",ncInt,posDim);

	// xVar.putAtt("units","m");
	// yVar.putAtt("units","m");
	// zVar.putAtt("units","m");
	double* x = new double[size];
	double* y = new double[size];
	double* z = new double[size];
	int* id = new int[size];

	for(uint i = start; i < end; ++i)
	{
		x[i-start] = pnp1[i].xi[0];
		y[i-start] = pnp1[i].xi[1];	
		#if SIMDIM == 3
			z[i-start] = pnp1[i].xi[2];
		#else
			z[i-start] = 0.0;
		#endif

		id[i-start] = pnp1[i].partID;
	}			

	xVar.putVar(x);
	yVar.putVar(y);
	zVar.putVar(z);

	idVar.putVar(id);

	if(outform == 1)
	{	/*Write fluid data*/
		NcVar press = zone.addVar("Pressure", ncDouble,posDim);
		NcVar dens = zone.addVar("Density",ncDouble,posDim);
		NcVar acc = zone.addVar("Acceleration",ncDouble,posDim);
		NcVar vel = zone.addVar("Velocity",ncDouble,posDim);

		// press.putAtt("units","Pa");
		// dens.putAtt("units","kg/m^3");
		// acc.putAtt("units","m/s^2");
		// vel.putAtt("units","m/s");

		double* p = new double[size];
		double* rho = new double[size];
		double* f = new double[size];
		double* v = new double[size];

		for(uint i = start; i < end; ++i)
		{
			p[i-start] = pnp1[i].p;
			rho[i-start] = pnp1[i].rho;
			f[i-start] = pnp1[i].f.norm();
			v[i-start] = pnp1[i].v.norm();
		}	

		press.putVar(p);
		dens.putVar(rho);
		acc.putVar(f);
		vel.putVar(v);
	}

	if(outform == 2)
	{	/*Write research data*/
		NcVar aero = zone.addVar("Aerodynamic Force", ncDouble,posDim);
		NcVar surf = zone.addVar("Surface Tension", ncDouble,posDim);
		NcVar Sx = zone.addVar("S_x", ncDouble,posDim);
		NcVar Sy = zone.addVar("S_y", ncDouble,posDim);
		NcVar B = zone.addVar("Boundary",ncInt,posDim);
		NcVar theta = zone.addVar("Theta",ncInt,posDim);

		// aero.putAtt("units","N");
		// surf.putAtt("units","N");
		// Sx.putAtt("units","N");
		// Sy.putAtt("units","N");
	
		double* af = new double[size]{0.0};
		double* sf = new double[size]{0.0};
		double* sx = new double[size]{0.0};
		double* sy = new double[size]{0.0};
		int* b = new int[size]{0};
		int* t = new int[size]{0};

		for(uint i = start; i < end; ++i)
		{
			af[i-start] = pnp1[i].Af.norm();
			sf[i-start] = pnp1[i].Sf.norm();
			sx[i-start] = pnp1[i].Sf(0);
			sy[i-start] = pnp1[i].Sf(1);
			b[i-start] = pnp1[i].b;
			t[i-start] = pnp1[i].theta;
		}		

		aero.putVar(af);
		surf.putVar(sf);
		Sx.putVar(sx);
		Sy.putVar(sy);
		B.putVar(b);
		theta.putVar(t);

		if (SIMDIM == 3)
		{
			NcVar Sz = zone.addVar("S_z",ncDouble,posDim);
			// Sz.putAtt("units","m");
			double* sz = new double[size]{0.0};
			for(uint i = start; i < end; ++i)
				sz[i-start] = pnp1[i].Sf[2];

			Sz.putVar(sz);
		}
	}	
}

// void Write_CDF_File(NcFile &nf, SIM &svar, State &pnp1)
// {	/*Write the timestep to a group*/
// 	string time = "Step #";
// 	// time.append(std::to_string(svar.t));
// 	time.append(std::to_string(svar.stepno));
// 	NcGroup ts(nf.addGroup(time));

// 	/*Write the Particles group*/
// 	Write_Group(ts, pnp1, svar.bndPts, svar.totPts, svar.outform);
// }


void Write_CDF_File(NcFile &nd, SIM &svar, const State &pnp1)
{	/*Write the timestep to a group*/
	string timefile = svar.outfolder;
	timefile.append("/h5/fuel_");
	std::ostringstream step;
	step << std::setfill('0') << std::setw(4);
	step << svar.frame;
	timefile.append(step.str());
	timefile.append(".h5part");

	NcFile nf(timefile, NcFile::replace);
	string steps = "Step#";
	steps.append(std::to_string(svar.frame));
	NcGroup ts(nf.addGroup(steps));
	NcDim tDim = ts.addDim("time",1);
	NcVar tVar = ts.addVar("time",ncDouble,tDim);
	tVar.putVar(&svar.t);

	/*Write the Particles group*/
	Write_Group(ts, pnp1, svar.bndPts, svar.totPts, svar.outform);

}

void Write_Boundary_CDF(const SIM &svar, const State &pnp1)
{	
	std::string file = svar.outfolder;
	file.append("/h5/Boundary.h5part");
	NcFile bound(file, NcFile::replace);
	NcGroup part(bound.addGroup("Boundary"));
	Write_Group(bound, pnp1, 0, svar.bndPts, svar.outform);
}

#endif