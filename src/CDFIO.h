#ifndef CDFIO_H
#define CDFIO_H

#include "Var.h"
#include "IOFunctions.h"
#include "Geometry.h"

#include <netcdf>
// using namespace netCDF;
// using namespace netCDF::exceptions;
#define NC_ERR 2
#define ERR(e)                                 \
	{                                          \
		printf("Error: %s\n", nc_strerror(e)); \
		exit(-1);                              \
	}

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
void Write_Zone(string input, ZONE &zn, MESH &cells)
{
	std::size_t found1 = input.find_first_of("/\\");
	std::size_t found2 = input.find_last_of("/\\");
	string path = input.substr(found1, found2 - found1 + 1);

	string filename = "Outputs";
	filename.append(path);
	filename.append(zn.name);
	filename.append(".plt");

	ofstream fout(filename, std::ios::out);
	if (!fout.is_open())
	{
		cout << "Failed to open data file for writing mesh." << endl;
		exit(-1);
	}

	fout << "TITLE = \"3D TAU Solution\"\n";
	fout << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\" \"x_velocity\" \"y_velocity\" \"z_velocity\"\n";
	fout << "ZONE T=\"" << zn.name << "\"" << endl;
	fout << " N=" << zn.nP << ", E=" << zn.nE << ", F=FEBLOCK, ET=" << zn.ETtype << endl
		 << endl;

	fout << std::left << std::scientific << std::setprecision(6);
	for (uint ii = 0; ii < SIMDIM; ++ii)
	{
		uint kk = 0;
		for (uint jj = cells.verts.size() - zn.nP; jj < cells.verts.size(); ++jj)
		{
			fout << cells.verts[jj][ii] << std::setw(15);
			kk++;

			if (kk == 5)
			{
				fout << "\n";
				kk = 0;
			}
		}

		if (kk % 5 != 0)
			fout << "\n";
	}

	for (uint ii = 0; ii < SIMDIM; ++ii)
	{
		uint kk = 0;
		for (uint jj = cells.cVel.size() - zn.nP; jj < cells.cVel.size(); ++jj)
		{
			fout << cells.cVel[jj][ii] << std::setw(15);
			kk++;

			if (kk == 5)
			{
				fout << "\n";
				kk = 0;
			}
		}

		if (kk % 5 != 0)
			fout << "\n";
	}

	// cout << "Velocitices written" << endl;
	fout << std::fixed;
	if (zn.ctype == 1 || zn.ctype == 2)
	{
		for (uint ii = cells.elems.size() - zn.nE; ii < cells.elems.size(); ++ii)
		{
			for (auto elem : cells.elems[ii])
			{
				fout << elem + 1 << std::setw(15);
			}
			fout << "\n";
		}
	}

	if (zn.ctype == 3)
	{
		for (uint ii = cells.elems.size() - zn.nE; ii < cells.elems.size(); ++ii)
		{
			for (uint jj = 0; jj < 3; ++jj)
			{
				fout << cells.elems[ii][jj] + 1 << " ";
			}
			fout << cells.elems[ii][2] + 1 << " ";
			for (uint jj = 3; jj < 6; ++jj)
			{
				fout << cells.elems[ii][jj] + 1 << " ";
			}
			fout << cells.elems[ii][5] + 1;
			fout << "\n";
		}
	}

	if (zn.ctype == 4)
	{
		for (uint ii = cells.elems.size() - zn.nE; ii < cells.elems.size(); ++ii)
		{
			for (uint jj = 0; jj < 4; ++jj)
			{
				fout << cells.elems[ii][jj] + 1 << " ";
			}
			fout << cells.elems[ii][4] + 1 << " ";
			fout << cells.elems[ii][4] + 1 << " ";
			fout << cells.elems[ii][4] + 1 << " ";
			fout << cells.elems[ii][4] + 1 << " ";
			fout << "\n";
		}
	}
	fout.close();
}

void Average_Point_to_Cell(vector<StateVecD> const &pData, vector<StateVecD> &cData,
						   vector<vector<size_t>> const &elems)
{
	vector<StateVecD> sum(elems.size(), StateVecD::Zero());

#pragma omp parallel for reduction(+ \
								   : sum)
	for (uint ii = 0; ii < elems.size(); ++ii)
	{

		uint const nVerts = elems[ii].size();
		for (auto jj : elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum;
}

void Average_Point_to_Cell(vector<real> const &pData, vector<real> &cData,
						   vector<vector<size_t>> const &elems)
{
	vector<real> sum(elems.size(), 0.0);

#pragma omp parallel for reduction(+ \
								   : sum)
	for (uint ii = 0; ii < elems.size(); ++ii)
	{

		uint const nVerts = elems[ii].size();
		for (auto jj : elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum;
}

vector<real> CpToPressure(vector<real> const &Cp, AERO const &avar)
{
	vector<real> press(Cp.size());
#pragma omp parallel for shared(Cp)
	for (uint ii = 0; ii < Cp.size(); ++ii)
	{
		press[ii] = Cp[ii] * avar.qInf /*+ fvar.gasPress*/;
	}
	return press;
}

/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
/*To run on the solution file*/
vector<real> Get_Scalar_Property_real(int &fin, string const &variable, size_t const &nPts)
{
#ifdef DEBUG
	dbout << "Reading variable: " << variable << endl;
#endif
	int retval;
	int varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Failed to get variable id for: " << variable << endl;
		ERR(retval);
	}

#ifdef DEBUG
	dbout << "Allocating array of: " << nPts << endl;
#endif

	double *array = new double[nPts];

#ifdef DEBUG
	dbout << "Attempting to read NetCDF variable:  " << variable << endl;
#endif

	if ((retval = nc_get_var_double(fin, varID, &array[0])))
	{
		cout << "Failed to get variable id for: " << variable << endl;
		ERR(retval);
	}

	/*Convert it to a vector to store*/
	vector<double> propVec;
	propVec.insert(propVec.end(), &array[0], &array[nPts]);
	vector<real> var = propVec;

#ifdef DEBUG
	dbout << "returning vector" << endl;
#endif
	return var;
}

vector<int> Get_Scalar_Property_int(int &fin, string variable, int nPts)
{
#ifdef DEBUG
	dbout << "Reading variable: " << variable << endl;
#endif
	int retval;
	int varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Failed to get variable id for: " << variable << endl;
		ERR(retval);
	}

#ifdef DEBUG
	dbout << "Allocating array of: " << nPts << endl;
#endif

	int *array = new int[nPts];

#ifdef DEBUG
	dbout << "Attempting to read NetCDF variable:  " << variable << endl;
#endif

	if ((retval = nc_get_var_int(fin, varID, &array[0])))
	{
		cout << "Failed to get variable data for: " << variable << endl;
		ERR(retval);
	}

	/*Convert it to a vector to store*/
	vector<int> propVec;
	propVec.insert(propVec.end(), &array[0], &array[nPts]);

#ifdef DEBUG
	dbout << "returning vector" << endl;
#endif
	return propVec;
}

/*To run on the mesh file*/
vector<vector<size_t>> Get_Element(int &fin, string const &variable, size_t const &nElem, size_t const &nPpEd)
{
#ifdef DEBUG
	dbout << "Reading \"" << variable << "\"" << endl;
#endif
	int retval;
	int varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Failed to get the variable ID of \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

// cout << nElem << "  " << nPoints << endl;
#ifdef DEBUG
	dbout << "Allocating array of: " << nElem << " by " << nPpEd << endl;
#endif

	vector<vector<size_t>> elemVec(nElem, vector<size_t>(nPpEd));
	
	try
	{
		/*Allocate on the heap (can be big datasets)*/
		int *elemArray = new int [nElem*nPpEd];
	
		size_t start[] = {0,0};
		size_t end[] = {nElem,nPpEd};

		/*Get the actual data from the file*/
		if ((retval = nc_get_vara_int(fin, varID, start,end, &elemArray[0])))
		{
			cout << "Failed to get the variable data of \"" << variable << "\"" << endl;
			ERR(retval);
			exit(-1);
		}

	#ifdef DEBUG
		dbout << "Attempting to read NetCDF elements." << endl;
	#endif

		cout << "Successfully read \"" << variable << "\"" << endl;
		cout << "Number of cells: " << nElem << endl;

	#ifdef DEBUG
		dbout << "Successfully read \"" << variable << "\"" << endl;
	#endif

		/*Convert it to a vector to store*/
		for (size_t ii = 0; ii < nElem; ++ii)
			for (size_t jj = 0; jj < nPpEd; ++jj)
				elemVec[ii][jj] = static_cast<size_t>(elemArray[index(ii,jj,nPpEd)]);

	}
	catch (std::bad_alloc &ba)
	{

		std::cerr << "Bad alloc caught. Failed to allocate \"" << variable << "\"" << endl;
		exit(-1);
	}

#ifdef DEBUG
	dbout << "Returning vector" << endl;
#endif
	return elemVec;
}

/*To run on the mesh file*/
uint Get_Coordinates(int &fin, size_t const &nPnts, vector<StateVecD> &coordVec)
{
	#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
	#endif
		
	int retval;
	int xcID, ycID, zcID;

	#if SIMDIM == 3

		if ((retval = nc_inq_varid(fin, "points_xc", &xcID)))
		{
			cout << "Failed to get the variable ID of: "
				<< "points_xc" << endl;
			cout << "Stopping" << endl;
			exit(-1);
		}

		if ((retval = nc_inq_varid(fin, "points_yc", &ycID)))
		{
			cout << "Failed to get the variable ID of: "
				<< "points_yc" << endl;
			cout << "Stopping" << endl;
			exit(-1);
		}

		if ((retval = nc_inq_varid(fin, "points_zc", &zcID)))
		{
			cout << "Failed to get the variable ID of: "
				<< "points_zc" << endl;
			cout << "Stopping" << endl;
			exit(-1);
		}

	#else
		uint ignore = 0;
		uint nFail = 0;
		if ((retval = nc_inq_varid(fin, "points_xc", &xcID)))
		{
			cout << "Failed to get the variable ID of: "
				<< "points_xc" << endl;
			cout << "Assuming this is the ignored dimension" << endl;
			ignore = 1;
			nFail++;
		}

		if ((retval = nc_inq_varid(fin, "points_yc", &ycID)))
		{
			cout << "Failed to get the variable ID of: "
				<< "points_yc" << endl;
			cout << "Assuming this is the ignored dimension" << endl;
			ignore = 2;
			nFail++;
		}

		if ((retval = nc_inq_varid(fin, "points_zc", &zcID)))
		{
			cout << "Failed to get the variable ID of: "
				<< "points_zc" << endl;
			cout << "Assuming this is the ignored dimension" << endl;
			ignore = 3;
			nFail++;
		}

		if(nFail > 1)
		{
			cout << "More than one dimension was not aquired, meaning something went wrong.\n\tStopping" << endl;
			exit(-1);
		}

	#endif 
	

	#if SIMDIM == 2
		#ifdef DEBUG
			dbout << "Ignored dimension: " << ignore << endl;
		#endif
		cout << "Ignored dimension: " << ignore << endl;
	#endif

	#ifdef DEBUG
		dbout << "Number of points: " << nPnts << endl;
	#endif
	/*Allocate on the heap (can be big datasets)*/
	real *coordX = new real[nPnts];
	real *coordY = new real[nPnts];

	coordVec = vector<StateVecD>(nPnts);
	#if SIMDIM == 2
		/*Get the actual data from the file*/
		if (ignore == 1)
		{
			if ((retval = nc_get_var_double(fin, ycID, &coordX[0])))
				ERR(retval);
			if ((retval = nc_get_var_double(fin, zcID, &coordY[0])))
				ERR(retval);
		}
		else if (ignore == 2)
		{
			if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
				ERR(retval);
			if ((retval = nc_get_var_double(fin, zcID, &coordY[0])))
				ERR(retval);
		}
		else if (ignore == 3)
		{
			if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
				ERR(retval);
			if ((retval = nc_get_var_double(fin, ycID, &coordY[0])))
				ERR(retval);
		}
		else
		{
			cout << "The ignored dimension was not found." << endl;
			exit(-1);
		}

		for (uint ii = 0; ii < nPnts; ++ii)
		{ /*Convert it to a vector to store*/
			coordVec[ii] = StateVecD(coordX[ii], coordY[ii]);
		}

	#else
		real *coordZ = new real[nPnts];
		if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
			ERR(retval);
		if ((retval = nc_get_var_double(fin, ycID, &coordY[0])))
			ERR(retval);
		if ((retval = nc_get_var_double(fin, zcID, &coordZ[0])))
			ERR(retval);

		for (uint ii = 0; ii < nPnts; ++ii)
		{ /*Convert it to a vector to store*/
			coordVec[ii] = StateVecD(coordX[ii], coordY[ii], coordZ[ii]);
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

void Get_Cell_Faces(vector<vector<uint>> const &cell,
					vector<vector<uint>> const &facenum, std::vector<std::vector<std::vector<uint>>> &cFaces)
{
	for (uint ii = 0; ii < cell.size(); ++ii)
	{
		uint jj = 0;
		cFaces.emplace_back();
		for (auto faces : facenum)
		{
			cFaces[ii].emplace_back();
			for (auto vert : faces)
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

void Read_SOLUTION(string const &solIn, FLUID const &fvar, AERO const &avar,
				   uint const ignored, MESH &cells, vector<uint> usedVerts)
{
	cout << "Reading solution file..." << endl;

	int retval;
	int solID;

	if ((retval = nc_open(solIn.c_str(), NC_NOWRITE, &solID)))
	{
		cout << "A netCDF error occured whilst trying to open the solution file:" << endl;
		cout << "\t" << solIn << endl
			 << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Solution file opened, reading contents..." << endl;

	int ptDimID;
	size_t solPts;

	if ((retval = nc_inq_dimid(solID, "no_of_points", &ptDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(solID, ptDimID, &solPts)))
	{
		ERR(retval);
		exit(-1);
	}

#ifdef DEBUG
	dbout << "Solution points: " << solPts << endl;
#endif

#if SIMDIM == 3
	if (solPts != cells.numPoint)
	{
		cout << "Solution file does not have the same number of vertices as the mesh." << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}
#else
	if (solPts / 2 != usedVerts.size())
	{
		cout << "Solution file size does not match size of mesh. Please check inputs." << endl;
	}
#endif

	vector<real> realDens = Get_Scalar_Property_real(solID, "density", solPts);

	/*Get the velocities*/
	vector<real> xvel, yvel, zvel;
#if SIMDIM == 3
	xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
	yvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
	zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);
#else
	if (ignored == 1)
	{
		xvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
		zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);
	}
	else if (ignored == 2)
	{
		xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
		zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);
	}
	else if (ignored == 3)
	{
		xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
		zvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
	}
#endif

	vector<StateVecD> vel(solPts);
	/*Test for size*/
	if (xvel.size() == solPts)
	{ /*Turn the arrays into a state vector*/
#pragma omp parallel for
		for (uint ii = 0; ii < solPts; ++ii)
		{
#if SIMDIM == 3
			vel[ii] = StateVecD(xvel[ii], yvel[ii], zvel[ii]);
#else
			vel[ii] = StateVecD(xvel[ii], zvel[ii]);
#endif
		}

		if (usedVerts.size() != 0)
		{ /*Get the used vertices*/
			vector<StateVecD> newVel;
			for (auto index : usedVerts)
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
	vector<real> press = Get_Scalar_Property_real(solID, "pressure", solPts);
	vector<real> dens = Get_Scalar_Property_real(solID, "density",solPts);

	// vector<real> dens(solPts);

	#pragma omp parallel for
	for (uint ii = 0; ii < press.size(); ++ii)
	{
		press[ii] -= avar.pRef; /* Want value to be gauge pressure */
		// dens[ii] = fvar.rho0 * pow((press[ii] / fvar.B + 1), 1 / fvar.gam);
	}

	if (press.size() == 0)
	{
		cout << "The solution data has not been correctly ingested. Please check the solution file." << endl;
	}

	/*Update the point arrays with the actual stuff*/
	if (usedVerts.size() != 0)
	{ /*Get the used vertices*/
		vector<real> newP;
		vector<real> newR;
		newP.reserve(usedVerts.size());
		newR.reserve(usedVerts.size());

		for (auto index : usedVerts)
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
	Average_Point_to_Cell(vel, cells.cVel, cells.elems);
	// Average_Point_to_Cell(dens,cells.SPHRho,cells.elems);
	Average_Point_to_Cell(press, cells.cP, cells.elems);
	Average_Point_to_Cell(dens, cells.cRho, cells.elems);

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

void Place_Edges(int &fin, size_t const &nElem, size_t const &nPnts, size_t const &nEdge, MESH &cells)
{
	#ifdef DEBUG
		dbout << "Reading element left/right and placing faces" << endl;
	#endif

	vector<int> left = Get_Scalar_Property_int(fin,"left_element_of_edges",nEdge);
	vector<int> right = Get_Scalar_Property_int(fin,"right_element_of_edges",nEdge);

	vector<std::pair<int, int>> leftright(nEdge);

	#pragma omp parallel shared(leftright)
	{
		/*Create local of */

		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < nEdge; ++ii)
		{
			int lindex = left[ii];
			int rindex = right[ii];
			leftright[ii] = std::pair<int, int>(lindex, rindex);
			#pragma omp critical
			{
				cells.cFaces[lindex].emplace_back(ii);
				if (rindex >= 0)
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
		for (uint ii = 0; ii < cells.cFaces.size(); ++ii)
		{
			for (auto const &index : cells.cFaces[ii])
			{
				vector<size_t> const face = cells.faces[index];
				for (auto const &vert : face)
				{
					if (std::find(cells.elems[ii].begin(), cells.elems[ii].end(), vert) == cells.elems[ii].end())
					{ /*Vertex doesn't exist in the elems vector yet.*/
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

		Average_Point_to_Cell(cells.verts, cells.cCentre, cells.elems);

		/*Find cell centres*/
		#pragma omp single
		{
			cout << "Finding cell volumes..." << endl;
		}

		// Find cell volumes
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < cells.elems.size(); ++ii)
		{
			cells.cVol[ii] = Cell_Volume(cells.verts, cells.faces, cells.elems[ii], cells.cFaces[ii], cells.cCentre[ii]);
		}
	}
	cout << "Elements built." << endl;

	#ifdef DEBUG
		dbout << "All elements defined." << endl
			<< endl;
	#endif
}

void Read_TAUMESH_EDGE(SIM &svar, MESH &cells, FLUID const &fvar, AERO const &avar)
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
	int retval;
	int meshID;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{
		cout << "A netCDF error occured whilst trying to open the mesh file:" << endl;
		cout << "\t" << meshIn << endl
			 << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading face data..." << endl;

	int ptDimID, elemDimID, edgeDimID, nPpEDimID;
	size_t nPnts, nElem, nEdge, nPpEd;

	

	// Retrieve how many elements there are.
	if ((retval = nc_inq_dimid(meshID, "no_of_elements", &elemDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, elemDimID, &nElem)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve edge number
	if ((retval = nc_inq_dimid(meshID, "no_of_edges", &edgeDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, edgeDimID, &nEdge)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve points per edge
	if ((retval = nc_inq_dimid(meshID, "points_per_edge", &nPpEDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, nPpEDimID, &nPpEd)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimid(meshID, "no_of_points", &ptDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, ptDimID, &nPnts)))
	{
		ERR(retval);
		exit(-1);
	}

	#ifdef DEBUG
		dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << endl;
	#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << endl;

	cells.numPoint = nPnts;
	cells.numElem = nElem;
	cells.numFace = nEdge;

	cells.elems = vector<vector<size_t>>(nElem);
	cells.cFaces = vector<vector<size_t>>(nElem);
	cells.verts = vector<StateVecD>(nPnts);
	cells.cVol = vector<real>(nElem);
	// cells.leftright = vector<std::pair<int,int>>(nFace);

	/*Get the faces of the mesh*/
	
	cells.faces = Get_Element(meshID, "points_of_element_edges", nEdge, nPpEd);

	/*Get the vertices that are in use to take the values from the solution file*/
	vector<int> usedVerts = Get_Scalar_Property_int(meshID, "vertices_in_use", nPnts);

	#ifdef DEBUG
		dbout << "Successfully read vertices_in_use data." << endl;
	#endif

	vector<uint> uVerts(nPnts);
	for (size_t ii = 0; ii < nPnts; ++ii)
		uVerts[ii] = static_cast<uint>(usedVerts[ii]);

	/*Get the coordinates of the mesh*/
	uint ignored = Get_Coordinates(meshID, nPnts, cells.verts);
	if (cells.verts.size() != nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for (auto &vert : cells.verts)
	{
		vert *= cells.scale;
	}
	svar.Start *= cells.scale;
	// svar.Jet/=cells.scale;

	/*Get face left and right, and put the faces in the elements*/
	Place_Edges(meshID, nElem, nPnts, nEdge, cells);

	#ifdef DEBUG
		dbout << "End of interaction with mesh file and ingested data." << endl
			<< endl;
		dbout << "Opening solultion file." << endl;
	#endif

	Read_SOLUTION(solIn, fvar, avar, ignored, cells, uVerts);

	for (size_t ii = 0; ii < cells.cMass.size(); ii++)
	{
		cells.cMass[ii] = cells.cRho[ii] * cells.cVol[ii];
	}
}

#endif

/*****************************************************************************/
/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/

void Place_Faces(int &fin, size_t const &nFace, MESH &cells)
{
	#ifdef DEBUG
		dbout << "Reading element left/right and placing faces" << endl;
	#endif

	vector<int> left = Get_Scalar_Property_int(fin,"left_element_of_faces",nFace);
	vector<int> right = Get_Scalar_Property_int(fin,"right_element_of_faces",nFace);

	vector<std::pair<int, int>> leftright(nFace);
	#pragma omp parallel shared(leftright)
	{
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < nFace; ++ii)
		{
			int lindex = left[ii];
			int rindex = right[ii];
			leftright[ii] = std::pair<int, int>(lindex, rindex);
			#pragma omp critical
			{
				cells.cFaces[lindex].emplace_back(ii);
				cells.cFaces[rindex].emplace_back(ii);
			}
		}
	}

	#ifdef DEBUG
		dbout << "End of placing faces in elements." << endl;
		dbout << "Building boundary indexes." << endl;
	#endif

	cells.leftright = leftright;

	/*Now go through the faces and see which vertices are unique, to get element data*/
	cout << "Building elements..." << endl;
	#pragma omp parallel
	{
		// Create list of unique vertex indices for each element
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < cells.cFaces.size(); ++ii)
		{
			for (auto const &index : cells.cFaces[ii])
			{
				vector<size_t> const face = cells.faces[index];
				for (auto const& vert : face)
				{
					if (std::find(cells.elems[ii].begin(), cells.elems[ii].end(), vert) == cells.elems[ii].end())
					{ /*Vertex doesn't exist in the elems vector yet.*/
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

		Average_Point_to_Cell(cells.verts, cells.cCentre, cells.elems);

		/*Find cell centres*/
		#pragma omp single
		{
			cout << "Finding cell volumes..." << endl;
		}

		// Find cell volumes
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < cells.elems.size(); ++ii)
		{
			cells.cVol[ii] = Cell_Volume(cells.verts, cells.faces, cells.elems[ii], cells.cFaces[ii], cells.cCentre[ii]);
		}
	}

	cout << "Elements built." << endl;

	#ifdef DEBUG
		dbout << "All elements defined." << endl << endl;
	#endif
}

void Read_TAUMESH_FACE(SIM &svar, MESH &cells, FLUID const &fvar, AERO const &avar)
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
	int retval;
	int meshID;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{
		cout << "A netCDF error occured whilst trying to open the mesh file:" << endl;
		cout << "\t" << meshIn << endl
			 << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading face data..." << endl;

	int ptDimID, elemDimID, faceDimID, nPpFDimID;
	size_t nPnts, nElem, nFace, nPpFc;

	
	// Retrieve how many elements there are.
	if ((retval = nc_inq_dimid(meshID, "no_of_elements", &elemDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, elemDimID, &nElem)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve face number
	if ((retval = nc_inq_dimid(meshID, "no_of_faces", &faceDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, faceDimID, &nFace)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve points per face 
	if ((retval = nc_inq_dimid(meshID, "points_per_face", &nPpFDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, nPpFDimID, &nPpFc)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimid(meshID, "no_of_points", &ptDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, ptDimID, &nPnts)))
	{
		ERR(retval);
		exit(-1);
	}

#ifdef DEBUG
	dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;

	cells.numPoint = nPnts;
	cells.numElem = nElem;
	cells.numFace = nFace;

	cells.elems = vector<vector<size_t>>(nElem);
	cells.cFaces = vector<vector<size_t>>(nElem);
	cells.cVol = vector<real>(nElem);

	cells.verts = vector<StateVecD>(nPnts);
	// cells.leftright = vector<std::pair<int,int>>(nFace);

	/*Get the faces of the mesh*/
	cells.faces = Get_Element(meshID, "points_of_element_faces", nFace, nPpFc);

	/*Get the coordinates of the mesh*/
	uint ignored = Get_Coordinates(meshID, nPnts, cells.verts);
	if (cells.verts.size() != nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for (auto &vert : cells.verts)
	{
		vert *= cells.scale;
	}
	svar.Start *= cells.scale;
	// svar.Jet *= cells.scale;

	/*Get face left and right, and put the faces in the elements*/
	Place_Faces(meshID, nFace, cells);

#ifdef DEBUG
	dbout << "End of interaction with mesh file and ingested data." << endl
		  << endl;
	dbout << "Opening solultion file." << endl;
#endif
	vector<uint> empty;
	Read_SOLUTION(solIn, fvar, avar, ignored, cells, empty);

	for (size_t ii = 0; ii < cells.cMass.size(); ii++)
	{
		cells.cMass[ii] = cells.cRho[ii] * cells.cVol[ii];
	}
}

#endif