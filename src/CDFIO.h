#ifndef CDFIO_H
#define CDFIO_H

#include "Var.h"
#include "IOFunctions.h"
#include "Geometry.h"

#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2

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
		for(uint jj = cells.cVel.size() - zn.nP; jj < cells.cVel.size(); ++jj)
		{
			fout << cells.cVel[jj][ii] << std::setw(15);
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


void Average_Point_to_Cell(vector<StateVecD> const& pData, vector<StateVecD>& cData,
							vector<vector<size_t>> const& elems)
{
	vector<StateVecD> sum(elems.size(),StateVecD::Zero());

	#pragma omp parallel for reduction(+:sum)
	for(uint ii = 0; ii < elems.size(); ++ii)
	{
		
		uint const nVerts = elems[ii].size();
		for (auto jj:elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum; 
}

void Average_Point_to_Cell(vector<real> const& pData, vector<real>& cData,
							vector<vector<size_t>> const& elems)
{
	vector<real> sum(elems.size(),0.0);

	#pragma omp parallel for reduction(+:sum)
	for(uint ii = 0; ii < elems.size(); ++ii)
	{
		
		uint const nVerts = elems[ii].size();
		for (auto jj:elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum; 
}

vector<real> CpToPressure(vector<real> const& Cp, AERO const& avar)
{
	vector<real> press(Cp.size());
	#pragma omp parallel for shared(Cp)
	for (uint ii = 0; ii < Cp.size(); ++ii)
	{
		press[ii] = Cp[ii]*avar.qRef /*+ fvar.gasPress*/;
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
uint Get_Coordinates(NcFile& fin, vector<StateVecD>& coordVec)
{
#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
#endif

#if SIMDIM == 2
	uint ignore=0;
#endif

	NcVar coordXD = fin.getVar("points_xc");
	NcVar coordYD = fin.getVar("points_yc");
	NcVar coordZD = fin.getVar("points_zc");	

	if(coordXD.isNull()) 
	{
#if SIMDIM == 3
		cout << "Missing x coordinate dimension!" << endl;
		exit(NC_ERR);
#else
		ignore = 1;
#endif
	}

	if(coordYD.isNull()) 
	{
#if SIMDIM == 3
		cout << "Missing y coordinate dimension!" << endl;
		exit(NC_ERR);
#else
		ignore = 2;
#endif
	}

	if(coordZD.isNull()) 
	{
#if SIMDIM == 3
		cout << "Missing z coordinate dimension!" << endl;
		exit(NC_ERR);
#else
		ignore = 3;
#endif
	}

#if SIMDIM == 2
#ifdef DEBUG
	dbout << "Ignored dimension: " << ignore << endl;
#endif
	cout << "Ignored dimension: " << ignore << endl;
#endif
	
	NcDim dim = coordXD.getDim(0);
	size_t nPts = dim.getSize();

	#ifdef DEBUG
		dbout << "Number of points: " << nPts << endl;
	#endif
	/*Allocate on the heap (can be big datasets)*/
	real* coordX = new real[nPts];
	real* coordY = new real[nPts];

	coordVec = vector<StateVecD>(nPts);
#if SIMDIM == 2	
	/*Get the actual data from the file*/
	if(ignore == 1)
	{
		coordYD.getVar(coordX);
		coordZD.getVar(coordY);
	}
	else if (ignore == 2)
	{
		coordXD.getVar(coordX);
		coordZD.getVar(coordY);
	}
	else if (ignore == 3)
	{
		coordXD.getVar(coordX);
		coordYD.getVar(coordY);
	}
	else 
	{
		cout << "The ignored dimension was not found." << endl;
		exit(-1);
	}

	for (uint ii = 0; ii < nPts; ++ii)
	{	/*Convert it to a vector to store*/
		coordVec[ii] = StateVecD(coordX[ii],coordY[ii]);
	}
	
#else
	real* coordZ = new real[nPts];
	coordXD.getVar(coordX);
	coordYD.getVar(coordY);	
	coordZD.getVar(coordZ);
	
	for (uint ii = 0; ii < nPts; ++ii)
	{	/*Convert it to a vector to store*/
		coordVec[ii] = StateVecD(coordX[ii],coordY[ii],coordZ[ii]);
	}
#endif

#ifdef DEBUG
	dbout << "Returning coordinates." << endl;
#endif

#if SIMDIM == 3
	return 0;	
#else
	return ignore;
#endif
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

void Get_Cell_Faces(vector<vector<uint>> const& cell,
	vector<vector<uint>> const& facenum, std::vector<std::vector<std::vector<uint>>>& cFaces)
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

void Read_SOLUTION(string const& solIn, FLUID const& fvar, AERO const& avar, 
	uint const ignored, MESH& cells, vector<uint> usedVerts, const uint& vinuse)
{
	cout << "Reading solution file..." << endl;

	NcFile sol;

	try
	{
		sol.open(solIn, NcFile::read);
	}
	catch (netCDF::exceptions::NcException& e)
	{
		cout << "A netCDF error occured whilst trying to open solution file:" << endl;
		cout << "what(): " << e.what() << endl << endl;

		cout << "Attemted file path: " << solIn << endl;
		exit(-1);
	}
	
	cout << "Solution file opened, reading contents..." << endl;

	NcDim solPN = sol.getDim("no_of_points");
	uint solPts = static_cast<uint>(solPN.getSize());

	#ifdef DEBUG
	dbout << "Solution points: " << solPts << endl;
	#endif


#if SIMDIM ==3
	if(solPts!=cells.numPoint)
	{
		cout << "Solution file does not have the same number of vertices as the mesh." << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}
#else
	if(solPts/2 != usedVerts.size() &&  vinuse == 1)
	{
		cout << "Solution file size does not match size of mesh. Please check inputs." << endl;

	}
#endif

	vector<real> realDens;
	Get_Scalar_Property(sol,"density",realDens);


	/*Get the velocities*/
	vector<real> xvel, yvel, zvel;
#if SIMDIM == 3
	Get_Scalar_Property(sol, "x_velocity", xvel);
	Get_Scalar_Property(sol, "y_velocity", yvel);
	Get_Scalar_Property(sol, "z_velocity", zvel);
#else 
	if(ignored == 1)
	{
		Get_Scalar_Property(sol, "y_velocity", xvel);
		Get_Scalar_Property(sol, "z_velocity", zvel);
	}
	else if (ignored == 2)
	{
		Get_Scalar_Property(sol, "x_velocity", xvel);
		Get_Scalar_Property(sol, "z_velocity", zvel);
	}
	else if(ignored == 3)
	{
		Get_Scalar_Property(sol, "x_velocity", xvel);
		Get_Scalar_Property(sol, "y_velocity", zvel);
	}
#endif

	vector<StateVecD> vel(solPts);
	/*Test for size*/
	if(xvel.size() == solPts)
	{	/*Turn the arrays into a state vector*/
		#pragma omp parallel for
		for(uint ii = 0; ii < solPts; ++ii)
		{
			#if SIMDIM == 3
			vel[ii] = StateVecD(xvel[ii],yvel[ii],zvel[ii]);
			#else
			vel[ii] = StateVecD(xvel[ii],zvel[ii]);
			#endif
		}


		if(usedVerts.size()!=0)
		{	/*Get the used vertices*/
			vector<StateVecD> newVel;
			for(auto index:usedVerts)
			{
				newVel.emplace_back(vel[index]);
				
			}
			vel = newVel;
		}
		else if (cells.verts.size() == solPts)
		{	// Do nothing
			// vel = vel;
		}
		else
		{
			cout << "Don't know which values to use!" << endl;
			exit(-1);
		}
	}
	else 
	{
		cout << "velocities do not have the same number of vertices as the mesh." << endl;
		cout << xvel.size() << "  " << solPts << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}

	/*Get other scalar data that may or may not exist*/
	vector<real> press;
	Get_Scalar_Property(sol,"pressure",press);
	if(press.size()==0)
	{
		/*Pressure data doesn't exist, so go to cp*/
		vector<real> cp;
		Get_Scalar_Property(sol,"cp",cp);

		/*Then convert cp to pressure. can get global attributes from sol file.. Maybe...*/
		press = CpToPressure(cp, avar);
	}
	else
	{
		#pragma omp parallel for
		for(uint ii = 0; ii < press.size(); ++ii)
			press[ii] -= avar.pRef;
	}

	
	vector<real> dens(solPts);

	#pragma omp parallel for
	for(uint ii = 0; ii < press.size(); ++ii)
	{
		dens[ii] = fvar.rho0 * pow((press[ii]/fvar.B + 1),1/fvar.gam);
	}	

	if(press.size() == 0)
	{
		cout << "The solution data has not been correctly ingested. Please check the solution file." << endl;
	}

	/*Update the point arrays with the actual stuff*/
	if(usedVerts.size()!=0)
	{	/*Get the used vertices*/
		vector<real> newP;
		vector<real> newR;
		for(auto index:usedVerts)
		{
			newP.emplace_back(press[index]);
			newR.emplace_back(dens[index]);
			
		}
		press = newP;
		dens = newR;	
	}
	else if (cells.verts.size() == solPts)
	{	//Do nothing
		// cells.pointP = press;
		// cells.pointRho = dens;	
	}
	else
	{
		cout << "Don't know which values to use!" << endl;
		exit(-1);
	}

	/*Average the data to a the cell*/
	cells.SetCells();
	

	cout << "Averaging points to cell centres..." << endl;
	Average_Point_to_Cell(vel,cells.cVel, cells.elems);
	// Average_Point_to_Cell(dens,cells.SPHRho,cells.elems);
	Average_Point_to_Cell(press,cells.cP,cells.elems);
	Average_Point_to_Cell(realDens,cells.cRho,cells.elems);
	
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
/*************** READING NETCDF EDGE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/

#if SIMDIM == 2
vector<vector<size_t>> Get_Edges(NcFile& fin)
{
#ifdef DEBUG
		dbout << "Reading edge data" << endl;
#endif

	NcVar edgeData = fin.getVar("points_of_element_edges");	
	if(edgeData.isNull())
	{	/*No data is available*/
		cout << "No edge data. Stopping." << endl;
		exit(-1);
	}

	NcDim dim = edgeData.getDim(0);
	size_t nEdge = dim.getSize();

	dim = edgeData.getDim(1);
	size_t nPnts = dim.getSize();
	// cout << nElem << "  " << nPoints << endl;
#ifdef DEBUG
	dbout << "Allocating array of: " << nEdge << " by " << nPnts << endl;
#endif

	/*Allocate on the heap (can be big datasets)*/		
	int* edgeArray = new int[nEdge*nPnts];
			
	/*Get the actual data from the file*/
	vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(nEdge);
	countp.push_back(nPnts);
	
#ifdef DEBUG
	dbout << "Attempting to read NetCDF points_of_element_edges." << endl;
#endif

	edgeData.getVar(startp,countp,edgeArray);

	cout << "Successfully read points_of_element_edges data." << endl;

#ifdef DEBUG
	dbout << "Successfully read points_of_element_edges data." << endl;
#endif

	/*Convert it to a vector to store*/
	size_t ii,jj;
	vector<vector<size_t>> edgeVec(nEdge,vector<size_t>(nPnts));
	for (ii = 0; ii < nEdge; ++ii)
	{
		for(jj = 0; jj < nPnts; ++jj)
			edgeVec[ii][jj] = static_cast<size_t>(edgeArray[index(ii,jj,nPnts)]);
	}

#ifdef DEBUG
	dbout << "Returning vector" << endl;
#endif
	return edgeVec;
	
}

void Place_Edges(NcFile& fin, MESH& cells)
{
	#ifdef DEBUG
		dbout << "Reading element left/right and placing faces" << endl;
	#endif

	NcVar leftData = fin.getVar("left_element_of_edges");	
	if(leftData.isNull())
	{	/*No data is available*/
		cout << "No cell left data. Stopping" << endl;
		exit(-1);
	}

	NcDim dim = leftData.getDim(0);
	size_t nLeft = dim.getSize();

	int* left = new int[nLeft];
	
	#ifdef DEBUG
	dbout << "Attempting to read NetCDF left_element_of_edges." << endl;
	#endif

	leftData.getVar(left);

	cout << "Successfully read left_element_of_edges data." << endl;

	#ifdef DEBUG
	dbout << "Successfully read left_element_of_edges data." << endl;
	#endif

	NcVar rightData = fin.getVar("right_element_of_edges");	
	if(leftData.isNull())
	{	/*No data is available*/
		cout << "No cell right data. Stopping" << endl;
		exit(-1);
	}

	dim = rightData.getDim(0);
	size_t nRight = dim.getSize();

	int* right = new int[nRight];
	
	#ifdef DEBUG
	dbout << "Attempting to read NetCDF right_element_of_edges." << endl;
	#endif

	rightData.getVar(right);

	cout << "Successfully read right_element_of_edges data." << endl;

	#ifdef DEBUG
	dbout << "Successfully read right_element_of_edges data." << endl;
	#endif
	vector<std::pair<int,int>> leftright(nLeft);

	#pragma omp parallel shared(leftright)
	{ 
		/*Create local of */
		

		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < nLeft; ++ii)
		{
			int lindex = left[ii];
			int rindex = right[ii];
			leftright[ii] = std::pair<int,int>(lindex,rindex);
			#pragma omp critical
			{
				cells.cFaces[lindex].emplace_back(ii);
				if(rindex >=0)
					cells.cFaces[rindex].emplace_back(ii);
			}
		}
	}

	cells.leftright = leftright;

	#ifdef DEBUG
		dbout << "End of placing edges in elements." << endl;
	#endif

	/*Now go through the faces and see which vertices are unique, to get element data*/
	cout << "Building elements..." << endl;
#pragma omp parallel
{
	// vector<vector<uint>> local;
	#pragma omp for schedule(static) nowait
	for(uint ii = 0; ii < cells.cFaces.size(); ++ii)
	{
		for(auto const& index:cells.cFaces[ii])
		{
			vector<size_t> const face = cells.faces[index];
			for(auto const& vert:face)
			{
				if(std::find(cells.elems[ii].begin(),cells.elems[ii].end(),vert)
					==cells.elems[ii].end())
				{	/*Vertex doesn't exist in the elems vector yet.*/
					#pragma omp critical
					{
						cells.elems[ii].emplace_back(vert);
					}
				}
			}
		}
	}


	/*Find cell centres*/
	#pragma omp single
	{
		cout << "Finding cell centres..." << endl;
	}

	Average_Point_to_Cell(cells.verts,cells.cCentre,cells.elems);

	/*Find cell centres*/
	#pragma omp single
	{
		cout << "Finding cell volumes..." << endl;
	}

	// Find cell volumes
	#pragma omp for schedule(static) nowait
	for(size_t ii = 0; ii < cells.elems.size(); ++ii)
	{
		cells.cVol[ii] = Cell_Volume(cells.verts,cells.faces,cells.elems[ii],cells.cFaces[ii],cells.cCentre[ii]);
	}
}
	cout << "Elements built." << endl;

	#ifdef DEBUG
	dbout << "All elements defined." << endl << endl;
	#endif
	

}

void Read_TAUMESH_EDGE(SIM& svar, MESH& cells, FLUID const& fvar, AERO const& avar)
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
	NcDim faceNo = fin.getDim("no_of_edges");
	uint nElem = static_cast<uint>(elemNo.getSize());
	uint nPnts = static_cast<uint>(pointNo.getSize());
	uint nEdge = static_cast<uint>(faceNo.getSize());

	#ifdef DEBUG
		dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << endl;
	#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << endl;
	
	cells.numPoint = nPnts;
	cells.numElem = nElem;

	cells.elems = vector<vector<size_t>>(nElem);
	cells.cFaces = vector<vector<size_t>>(nElem);
	cells.verts = vector<StateVecD>(nPnts);
	cells.cVol = vector<real>(nElem);
	// cells.leftright = vector<std::pair<int,int>>(nFace);

	/*Get the faces of the mesh*/
	cells.faces = Get_Edges(fin);

	/*Get the vertices that are in use to take the values from the solution file*/
	int* usedVerts = new int[nPnts];
	uint vinuse = 0;
#ifdef DEBUG 
	dbout << "Attempting to read NetCDF vertices_in_use data." << endl;
#endif
	NcVar uVar = fin.getVar("vertices_in_use");
	if(uVar.isNull())
	{
		cout << "Warning: no information on vertices used in solution file." << endl;
		vinuse = 0;
	}
	else
	{
		uVar.getVar(usedVerts);
		vinuse = 1;
	}

#ifdef DEBUG 
	dbout << "Successfully read vertices_in_use data." << endl;
#endif

	vector<uint> uVerts(nPnts);
	for(size_t ii = 0; ii < nPnts; ++ii)
		uVerts[ii] = usedVerts[ii];


	/*Get the coordinates of the mesh*/
	uint ignored = Get_Coordinates(fin,cells.verts);
	if(cells.verts.size()!= nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for(auto& vert:cells.verts)
	{
		vert*=cells.scale;
	}
	svar.Start*=cells.scale;
	// svar.Jet/=cells.scale;

	/*Get face left and right, and put the faces in the elements*/
	Place_Edges(fin, cells);

	#ifdef DEBUG
	dbout << "End of interaction with mesh file and ingested data." << endl << endl;
	dbout << "Opening solultion file." << endl;
	#endif

	Read_SOLUTION(solIn,fvar,avar,ignored,cells,uVerts, vinuse);

	for(size_t ii = 0; ii < cells.cMass.size(); ii++)
	{
		cells.cMass[ii] = cells.cRho[ii]*cells.cVol[ii]; 
	}
}

#endif


/*****************************************************************************/
/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
vector<vector<size_t>> Get_Faces(NcFile& fin)
{
#ifdef DEBUG
		dbout << "Reading face data" << endl;
#endif

	NcVar faceData = fin.getVar("points_of_element_faces");	
	if(faceData.isNull())
	{	/*No data is available*/
		cout << "No face data. Stopping." << endl;
		exit(-1);
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
	size_t ii,jj;
	vector<vector<size_t>> faceVec(nFace,vector<size_t>(nPnts));
	for (ii = 0; ii < nFace; ++ii)
	{
		for(jj = 0; jj < nPnts; ++jj)
			faceVec[ii][jj] = static_cast<size_t>(faceArray[index(ii,jj,nPnts)]);
	}

#ifdef DEBUG
	dbout << "Returning vector" << endl;
#endif
	return faceVec;
	
}

void Place_Faces(NcFile& fin, MESH& cells)
{
	#ifdef DEBUG
		dbout << "Reading element left/right and placing faces" << endl;
	#endif

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
	vector<std::pair<int,int>> leftright(nLeft);
	vector<size_t> bIndex;
	#pragma omp parallel shared(leftright)
	{ 
		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < nLeft; ++ii)
		{
			int lindex = left[ii];
			int rindex = right[ii];
			leftright[ii] = std::pair<int,int>(lindex,rindex);
			#pragma omp critical
			{
				cells.cFaces[lindex].emplace_back(ii);
				if(rindex >=0)
					cells.cFaces[rindex].emplace_back(ii);
				else if(rindex == -1)
				{	/*Create a list of vertices on the boundary*/
					for(auto const& vindex:cells.faces[ii])
					{
						bIndex.emplace_back(vindex);
					}
				}
			}
		}
	}

	#ifdef DEBUG
		dbout << "End of placing faces in elements." << endl;
		dbout << "Building boundary indexes." << endl;
	#endif


	/*Sort and erase duplicates*/
	std::sort(bIndex.begin(),bIndex.end());
	bIndex.erase(std::unique(bIndex.begin(),bIndex.end()),bIndex.end());
	
	vector<StateVecD> bVerts(bIndex.size());
	cells.bIndex = bIndex;
	cells.leftright = leftright;	

	/*Now go through the faces and see which vertices are unique, to get element data*/
	cout << "Building elements..." << endl;
	#pragma omp parallel
	{
		
		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < bIndex.size(); ++ii)
		{
			bVerts[ii] = cells.verts[bIndex[ii]];
		}

		// Create list of unique vertex indices for each element
		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < cells.cFaces.size(); ++ii)
		{
			for(auto const& index:cells.cFaces[ii])
			{
				vector<size_t> const face = cells.faces[index];
				for(auto const& vert:face)
				{
					if(std::find(cells.elems[ii].begin(),cells.elems[ii].end(),vert)
						==cells.elems[ii].end())
					{	/*Vertex doesn't exist in the elems vector yet.*/
						#pragma omp critical
						{
							cells.elems[ii].emplace_back(vert);
						}
					}
				}
			}
		}

		/*Find cell centres*/
		#pragma omp single
		{
			cout << "Finding cell centres..." << endl;
		}

		Average_Point_to_Cell(cells.verts,cells.cCentre,cells.elems);

		/*Find cell centres*/
		#pragma omp single
		{
			cout << "Finding cell volumes..." << endl;
		}

		// Find cell volumes
		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < cells.elems.size(); ++ii)
		{
			cells.cVol[ii] = Cell_Volume(cells.verts,cells.faces,cells.elems[ii],cells.cFaces[ii],cells.cCentre[ii]);
		}
	}

	cells.bVerts = bVerts;
	cout << "Elements built." << endl;

	#ifdef DEBUG
	dbout << "All elements defined." << endl << endl;
	#endif


}

void Read_TAUMESH_FACE(SIM& svar, MESH& cells, FLUID const& fvar, AERO const& avar)
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
	NcDim faceNo = fin.getDim("no_of_faces");
	NcDim pointNo = fin.getDim("no_of_points");
	
	uint nElem = static_cast<uint>(elemNo.getSize());
	uint nPnts = static_cast<uint>(pointNo.getSize());
	uint nFace = static_cast<uint>(faceNo.getSize());

	#ifdef DEBUG
		dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
	#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
	
	cells.numPoint = nPnts;
	cells.numElem = nElem;

	cells.elems = vector<vector<size_t>>(nElem);
	cells.cFaces = vector<vector<size_t>>(nElem);
	cells.cVol = vector<real>(nElem);

	cells.verts = vector<StateVecD>(nPnts);
	// cells.leftright = vector<std::pair<int,int>>(nFace);

	/*Get the faces of the mesh*/
	cells.faces = Get_Faces(fin);

	/*Get the coordinates of the mesh*/
	uint ignored = Get_Coordinates(fin, cells.verts);
	if(cells.verts.size()!= nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for(auto& vert:cells.verts)
	{
		vert*=cells.scale;
	}
	svar.Start*=cells.scale;
	// svar.Jet/=cells.scale;
	

	/*Get face left and right, and put the faces in the elements*/
	Place_Faces(fin, cells);

	#ifdef DEBUG
	dbout << "End of interaction with mesh file and ingested data." << endl << endl;
	dbout << "Opening solultion file." << endl;
	#endif
	vector<uint> empty;
	uint empt;
	Read_SOLUTION(solIn,fvar,avar,ignored,cells,empty,empt);

	for(size_t ii = 0; ii < cells.cMass.size(); ii++)
	{
		cells.cMass[ii] = cells.cRho[ii]*cells.cVol[ii]; 
	}
}


/*****************************************************************************/
/****************** WRITING NETCDF DATA FUNCTIONS ****************************/
/*****************************************************************************/

void Write_Group(NcGroup& zone, State const& pnp1, 
	uint const start, uint const end, uint const outform)
{
	/*Write the Particles group*/
	
	uint const size = end - start;
	NcDim posDim = zone.addDim("Point",size);

	NcVar xVar = zone.addVar("x",ncFloat,posDim);
	NcVar yVar = zone.addVar("y",ncFloat,posDim);
	NcVar zVar = zone.addVar("z",ncFloat,posDim);
	NcVar idVar = zone.addVar("id",ncInt,posDim);

	// xVar.putAtt("units","m");
	// yVar.putAtt("units","m");
	// zVar.putAtt("units","m");
	real* x = new real[size];
	real* y = new real[size];
	real* z = new real[size];
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
		NcVar press = zone.addVar("Pressure", ncFloat,posDim);
		NcVar dens = zone.addVar("Density",ncFloat,posDim);
		NcVar acc = zone.addVar("Acceleration",ncFloat,posDim);
		NcVar vel = zone.addVar("Velocity",ncFloat,posDim);

		// press.putAtt("units","Pa");
		// dens.putAtt("units","kg/m^3");
		// acc.putAtt("units","m/s^2");
		// vel.putAtt("units","m/s");

		real* p = new real[size];
		real* rho = new real[size];
		real* f = new real[size];
		real* v = new real[size];

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
		NcVar aero = zone.addVar("Aerodynamic Force", ncFloat,posDim);
		NcVar surf = zone.addVar("Surface Tension", ncFloat,posDim);
		NcVar Sx = zone.addVar("S_x", ncFloat,posDim);
		NcVar Sy = zone.addVar("S_y", ncFloat,posDim);
		NcVar B = zone.addVar("Boundary",ncInt,posDim);
		NcVar theta = zone.addVar("Theta",ncInt,posDim);

		// aero.putAtt("units","N");
		// surf.putAtt("units","N");
		// Sx.putAtt("units","N");
		// Sy.putAtt("units","N");
	
		real* af = new real[size]{0.0};
		real* sf = new real[size]{0.0};
		real* sx = new real[size]{0.0};
		real* sy = new real[size]{0.0};
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
			NcVar Sz = zone.addVar("S_z",ncFloat,posDim);
			// Sz.putAtt("units","m");
			real* sz = new real[size]{0.0};
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


void Write_CDF_File(NcFile& nd, SIM& svar, State const& pnp1)
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
	NcVar tVar = ts.addVar("time",ncFloat,tDim);
	tVar.putVar(&svar.t);

	/*Write the Particles group*/
	Write_Group(ts, pnp1, svar.bndPts, svar.totPts, svar.outform);

}

void Write_Boundary_CDF(SIM const& svar, State const& pnp1)
{	
	std::string file = svar.outfolder;
	file.append("/h5/Boundary.h5part");
	NcFile bound(file, NcFile::replace);
	NcGroup part(bound.addGroup("Boundary"));
	Write_Group(bound, pnp1, 0, svar.bndPts, svar.outform);
}

#endif