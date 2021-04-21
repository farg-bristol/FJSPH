/*Cell based to Face based data converter*/

#include <vector>
#include <fstream>
#include <iomanip>
#include <netcdf>
#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"

// using namespace netCDF;
// using namespace netCDF::exceptions;

#define NC_ERR 2
#define ERR(e)                                 \
	{                                          \
		printf("Error: %s\n", nc_strerror(e)); \
		exit(-1);                              \
	}


#ifdef DEBUG
	/*Open debug file to write to*/
	std::ofstream dbout("Cell2Face.log",std::ios::out);
#endif

// #ifndef NTHREADS
// #define NTHREADS 4
// #endif

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 3
#endif

/* Define data type. */
typedef double real;
typedef unsigned int uint;

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real,SIMDIM,1> StateVecD;
typedef Eigen::Matrix<int,SIMDIM,1> StateVecI;

using std::vector;
using std::cout;
using std::endl;
using std::string; 
using std::setw;
typedef struct CELL
{
	/*Standard contructor*/
	CELL(const uint nElem, const uint nPts)
	{
		numElem = nElem; 
		numPoint = nPts;
		verts.reserve(numPoint);
		elems.reserve(numElem);
	}
	
	/*Zone info*/
	uint numElem, numPoint;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Cell based data*/
	vector<vector<uint>> elems;

	/*Surface faces*/
	vector<vector<uint>> sfaces;
}CELL;

typedef class FACE
{
	public:
	FACE(const CELL& cdata): numElem(cdata.numElem), numPoint(cdata.numPoint)
	{
		verts = cdata.verts;
		// numElem = cdata.numElem;
		// numPoint = cdata.numPoint;
		numFaces = 0; nFar = 0; nWall = 0;
		// nFacesPElem = new int[cdata.numElem];
	}	

	FACE() : numElem(0), numPoint(0)
	{
		numFaces = 0; nFar = 0; nWall = 0;
	};

	void insert(const FACE& flocal)
	{
		faces.insert(faces.end(),flocal.faces.begin(),flocal.faces.end());
		celllr.insert(celllr.end(),flocal.celllr.begin(),flocal.celllr.end());
		numFaces += flocal.numFaces;
		nFar += flocal.nFar;
		nWall += flocal.nWall;
	}

	uint numFaces, nFar, nWall;
	vector<StateVecD> verts;
	vector<vector<uint>> faces; /*Face indexes*/
	vector<std::pair<int,int>> celllr; /*Cell left and right of the face*/

	// int* nFacesPElem;

	const uint numElem, numPoint;
}FACE;

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

vector<int> Find_Bmap_Markers(string const& bmapIn)
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

	vector<int> markers;
	/*Search for the markers that have either farfield or symmetry plane*/
	uint lineno = 0;
	string line;

	uint blockno = 0;
	uint blockstart = 1; /*A store of when this block starts to search through*/

	while(getline(fin,line))
	{

		// cout << line << endl;
		if(line.find("Type: laminar wall")!=string::npos ||
			line.find("Type: euler wall")!=string::npos || 
			line.find("Type: viscous wall")!=string::npos ||
			line.find("Type: sharp edge")!=string::npos)
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
							markers.emplace_back(found);
						}

						temp = "";
					}

					/*Go back to where we were*/
					blockstart = lineno+2;
					GotoLine(fin,lineno+2);
					break;
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

	cout << "Wall markers:" << endl;

	for(auto mark:markers)
	{
		cout << mark << "  " ;
	}
	cout << endl;

	#ifdef DEBUG
	dbout << "Exiting Find_Bmap_Markers..." << endl;
	#endif

	return markers;
}

void Get_Dimension(int& fin, string variable, int& dimID, size_t& dimVal)
{
	int retval;
	if ((retval = nc_inq_dimid(fin, variable.c_str(), &dimID)))
	{
		cout << "Error: failed getting the ID of dimension \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(fin, dimID, &dimVal)))
	{
		cout << "Error: failed getting the length of dimension \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}
}

vector<vector<uint>> Get_Surface(int& fin, vector<int> const& markers)
{
	#ifdef DEBUG
	dbout << "Entering Get_Surface..." << endl;
	#endif
	int retval;

	vector<vector<uint>> faceVec;

	int surfTDim, surfQDim, pPSTDim, pPSQDim;
	size_t nsurfT, nsurfQ, npPST, npPSQ;

	/* Get surface dimensions */
	if ((retval = nc_inq_dimid(fin, "no_of_surfacetriangles", &surfTDim)))
	{
		cout << "No surfacetriangle data" << endl;
	}
	else
	{
		Get_Dimension(fin, "no_of_surfacetriangles", surfTDim, nsurfT);
		Get_Dimension(fin, "points_per_surfacetriangle", pPSTDim, npPST);

		#ifdef DEBUG
			dbout << "Number of triangles: " << nsurfT << endl;
		#endif

		int *tris = new int[nsurfT * npPST];

		/*Get the actual data from the file*/
		size_t startp[] = {0,0};
		size_t countp[] = {nsurfT, npPST};
	

		#ifdef DEBUG
			dbout << "Attempting to read NetCDF surface triangles." << endl;
		#endif

		int varID;
		if ((retval = nc_inq_varid(fin, "points_of_surfacetriangles", &varID)))
		{
			cout << "Error: failed to get ID for variable \"points_of_surfacetriangles\"" << endl;
			ERR(retval);
			exit(-1);
		}
			

		if ((retval = nc_get_vara_int(fin, varID, startp, countp, &tris[0])))
		{
			cout << "Error: failed to get data for variable  \"points_of_surfacetriangles\"" << endl;
			ERR(retval);
			exit(-1);
		}

		#ifdef DEBUG
			dbout << "Putting surface faces into a vector." << endl;
		#endif

		vector<vector<uint>> localVec(nsurfT, vector<uint>(npPST));
		for (uint ii = 0; ii < nsurfT; ++ii)
		{
			for (uint jj = 0; jj < npPST; ++jj)
				localVec[ii][jj] = (static_cast<uint>(tris[index(ii, jj, npPST)]));
		}

		faceVec.insert(faceVec.end(), localVec.begin(), localVec.end());
	}


	if ((retval = nc_inq_dimid(fin, "no_of_surfacequadrilaterals", &surfTDim)))
	{
		cout << "No surfacequadrilateral data" << endl;
	}
	else
	{
		Get_Dimension(fin, "no_of_surfacequadrilaterals", surfQDim, nsurfQ);
		Get_Dimension(fin, "points_per_surfacequadrilateral", pPSQDim, npPSQ);
		
		#ifdef DEBUG
			dbout << "Number of quadrilaterals: " << nsurfQ << endl;
		#endif

		/*Allocate on the heap (can be big datasets)*/
		
		int* quads = new int[nsurfQ*npPSQ];

		size_t startp[] = {0, 0};
		size_t countp[] = {nsurfQ, npPSQ};

		#ifdef DEBUG
			dbout << "Attempting to read NetCDF surface quadrilaterals." << endl;
		#endif

		int varID;
		if ((retval = nc_inq_varid(fin, "points_of_surfacequadrilaterals", &varID)))
		{
			cout << "Error: failed to get ID for variable \"points_of_surfacequadrilaterals\"" << endl;
			ERR(retval);
			exit(-1);
		}

		if ((retval = nc_get_vara_int(fin, varID, startp, countp, &quads[0])))
		{
			cout << "Error: failed to get data for variable  \"points_of_surfacequadrilaterals\"" << endl;
			ERR(retval);
			exit(-1);
		}

		#ifdef DEBUG
			dbout << "Putting surface faces into a vector." << endl;
		#endif
		/*Convert it to a vector to store*/
		vector<vector<uint>> localVec(nsurfQ,vector<uint>(npPSQ));

		for (uint ii = 0; ii < nsurfQ; ++ii)
		{
			for(uint jj = 0; jj < npPSQ; ++jj)
				localVec[ii][jj] = (static_cast<uint>(quads[index(ii,jj,npPSQ)]));
		}

		faceVec.insert(faceVec.end(),localVec.begin(),localVec.end());
	}

	/*Get the boundarymarkers*/
	int boundMDim;
	size_t nbound;
	Get_Dimension(fin, "no_of_surfaceelements", boundMDim, nbound);

	uint nMarkers = static_cast<uint>(nbound);
	if(faceVec.size() != nMarkers)
	{
		cout << "Mismatch of number of surface elements defined and number ingested." << endl;
		cout << "Number of surface elements: " << nMarkers << 
		"  Number in vector: " << faceVec.size() << endl;
	}

	int* faceMarkers = new int[nMarkers];

	int varID;
	if ((retval = nc_inq_varid(fin,"boundarymarker_of_surfaces", &varID)))
	{
		cout << "Error: failed to get ID for variable \"boundarymarker_of_surfaces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_get_var_int(fin, varID, &faceMarkers[0])))
	{
		cout << "Error: failed to get data for variable  \"boundarymarker_of_surfaces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	vector<vector<uint>> farVec;

	/*Want to find which surfaces are the ones I want to look for*/
	for(size_t ii = 0; ii < nMarkers; ++ii)
	{
		if(std::find(markers.begin(),markers.end(),faceMarkers[ii])!= markers.end())
		{	/*The face is an inner boundary*/
			/*Pre sort to save time in the loop*/
			vector<uint> v = faceVec[ii];
			std::sort(v.begin(),v.end());
			farVec.emplace_back(v);
		}
	}

	#ifdef DEBUG
	dbout << "Exiting Get_Surface..." << endl;
	#endif

	return farVec;
}

double* Get_Real_Scalar(int& fin, string variable, int const& nPnts)
{
	int retval, varID;
	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Error: failed to get ID for variable \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	double* var = new double[nPnts];

	if ((retval = nc_get_var_double(fin, varID, &var[0])))
	{
		cout << "Error: failed to get data for variable  \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	return var;
}

/*To run on the mesh file*/
vector<vector<uint>> Get_Element(int& fin, string variable, size_t& nElem, size_t& nPnts)
{
	#ifdef DEBUG
		dbout << "Reading Element: " << variable << endl;
	#endif
	int retval, varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Error: failed to get ID for variable \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	// cout << nElem << "  " << nPoints << endl;
	#ifdef DEBUG
	dbout << "Allocating array of: " << nElem << " by " << nPnts << endl;
	#endif

	/*Allocate on the heap (can be big datasets)*/		
	int* elemArray = new int[nElem*nPnts];
			
	/*Get the actual data from the file*/
	size_t startp[] = {0,0};
	size_t countp[] = {nElem,nPnts};

	#ifdef DEBUG
	dbout << "Attempting to read NetCDF elements." << endl;
	#endif

	if ((retval = nc_get_vara_int(fin, varID, startp, countp, &elemArray[0])))
	{
		cout << "Error: failed to get data for variable  \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Successfully read: " << variable << endl;
	cout << "Number of elements: " << nElem << endl;

	#ifdef DEBUG
	dbout << "Successfully read elements" << endl;
	#endif

	/*Convert it to a vector to store*/
	uint ii,jj;
	vector<vector<uint>> elemVec(nElem,vector<uint>(nPnts));
	for (ii = 0; ii < nElem; ++ii)
	{
		for(jj = 0; jj < nPnts; ++jj)
			elemVec[ii][jj] = static_cast<uint>(elemArray[index(ii,jj,nPnts)]);
	}

	#ifdef DEBUG
	dbout << "Returning vector" << endl;
	#endif
	return elemVec;
}

/*To run on the mesh file*/
vector<StateVecD> Get_Coordinates(int& fin, size_t const& nPnts)
{
	#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
	#endif

	/*Get the coordinate data*/
	double* coordX = Get_Real_Scalar(fin, "points_xc", nPnts);
	double* coordY = Get_Real_Scalar(fin, "points_yc", nPnts);
	double* coordZ = Get_Real_Scalar(fin, "points_zc", nPnts);

	/*Convert it to a vector to store*/
	vector<StateVecD> coordVec(nPnts);
	for (uint ii = 0; ii < nPnts; ++ii)
	{
		#if SIMDIM == 3
		coordVec[ii] = StateVecD(coordX[ii],coordY[ii],coordZ[ii]);
		#else
		coordVec[ii] = StateVecD(coordX[ii],coordY[ii]);
		#endif
	}
	#ifdef DEBUG
		dbout << "Returning coordinates." << endl;
	#endif
	return coordVec;	
}

void Add_Face(vector<uint> const& face, std::pair<int,int> const& leftright,
	FACE& fdata)
{
	/*Break the face down into two triangles*/
	if(face.size() == 4)
	{	/*Break the face down into two triangles*/
		fdata.celllr.emplace_back(leftright);
		fdata.celllr.emplace_back(leftright);
		vector<uint> face1 = {face[0],face[1],face[2]};
		vector<uint> face2 = {face[0],face[2],face[3]};
		fdata.faces.emplace_back(face1);
		fdata.faces.emplace_back(face2);
		fdata.numFaces+=2;
	}
	else
	{	/*Face is already a triangle. No work to be done*/
		fdata.celllr.emplace_back(leftright);
		fdata.faces.emplace_back(face);
		fdata.numFaces++;
	}
}

void CheckFaces(const vector<vector<uint>>& vertincells,
	 const vector<vector<uint>>& lfaces, const uint lindex, const uint lgeom, 
	const CELL& cdata, vector<vector<uint>>& colour, FACE& fdata)
{
	std::pair<int,int> leftright;
	uint lfaceindex = 0;
	for(auto const& face:lfaces)
	{	/*Define that the cell of the top-level cell search is on the 'left'*/

		/*Create an ordered list of the left side of the face*/
		vector<uint> lface = face;
		std::sort(lface.begin(),lface.end());

		/* Search through each vertex of the face, and see if this face exists*/
		for(auto const& vert:face)
		{	/*Search through the cells where the vertex is mentioned*/
			for (auto const& rindex:vertincells[vert])
			{	
				/*If the cell is the current cell, ignore*/
				if(rindex==lindex)
					continue;	

				/*Define that the cell of the inner search is on the 'right'*/
				vector<uint> rcell = cdata.elems[rindex];
				vector<vector<uint>> rfaces;

				if(rcell.size() == 4)
				{	/*Tetraeder*/
					/*If left cell is a hex, they can't share a face*/
					if(lgeom == 3)
						continue;

					rfaces = {{rcell[0],rcell[1],rcell[2]},
							  {rcell[1],rcell[0],rcell[3]},
					   		  {rcell[2],rcell[3],rcell[0]},
					   		  {rcell[3],rcell[1],rcell[2]}};
				}
				else if(rcell.size() == 5)
				{	/*Pyra*/
					rfaces = {{rcell[0],rcell[3],rcell[2],rcell[1]},
					   		  {rcell[1],rcell[4],rcell[0]},
					   		  {rcell[2],rcell[4],rcell[1]},
					   		  {rcell[3],rcell[4],rcell[2]},
					   		  {rcell[4],rcell[3],rcell[0]}};						
				}
				else if(rcell.size() == 6)
				{	/*Prism*/
					rfaces = {{rcell[1],rcell[2],rcell[5],rcell[4]},
					   		  {rcell[4],rcell[3],rcell[0],rcell[1]},
					   		  {rcell[3],rcell[5],rcell[2],rcell[0]},
					   		  {rcell[2],rcell[0],rcell[1]},
					   		  {rcell[5],rcell[3],rcell[4]}};
				}
				else if(rcell.size() == 8)
				{	/*Hexaeder*/

					/*If left cell is a tet, they can't share a face*/
					if(lgeom == 0)
						continue;

					rfaces = {{rcell[0],rcell[1],rcell[2],rcell[3]},
							  {rcell[4],rcell[0],rcell[3],rcell[7]},
							  {rcell[1],rcell[0],rcell[4],rcell[5]},
							  {rcell[2],rcell[1],rcell[5],rcell[6]},
							  {rcell[3],rcell[2],rcell[6],rcell[7]},
							  {rcell[7],rcell[6],rcell[5],rcell[4]}};				 
				}

				uint rfaceindex = 0;
				for(auto& rface:rfaces)
				{
					std::sort(rface.begin(),rface.end());
					if(lface == rface)
					{
						/*Then the face is a match*/
						if(colour[lindex][lfaceindex] == 0 && colour[rindex][rfaceindex] == 0)
						{
							leftright = std::pair<int,int>(lindex,rindex);
							colour[lindex][lfaceindex] = 1;
							colour[rindex][rfaceindex] = 1;
							Add_Face(face,leftright,fdata);
							goto matchfound;
						}

					}
					rfaceindex++;
				}	
			}
		}
		
		/*If a match has not been found, the face must be a boundary face*/	
		if(colour[lindex][lfaceindex] == 0)
		{			
			for (auto sface:cdata.sfaces)
			{	/*Search through the surface faces to identify */
				/*if the face is an internal face or not*/
				if(sface == lface)
				{	/*Face is an internal face*/
					if (lface.size() == 4)
						fdata.nFar+=2;
					else
						fdata.nFar++;
					
					leftright = std::pair<int,int>(lindex,-1);
					colour[lindex][lfaceindex] = 1;
					Add_Face(face,leftright,fdata);
					goto matchfound;
				}
			}

			/*If still uncoloured, then face is an external boundary*/
			if (lface.size() == 4)
				fdata.nWall += 2;
			else
				fdata.nWall++;

			leftright = std::pair<int,int>(lindex,-2);
			colour[lindex][lfaceindex] = 1;
			Add_Face(face,leftright,fdata);
			
		}

matchfound:			
		lfaceindex++;
	}
}

void BuildFaces(const CELL& cdata, FACE& fdata)
{
	#ifdef DEBUG
	dbout << "Entering BuildFaces..." << endl;
	#endif
	uint nElem = cdata.numElem;
	uint nPts = cdata.numPoint;
	uint unmfaces=0;

	// for(auto const& vert:cdata.elems[33])
	// 	cout << vert << "  ";

	// cout << endl;

	/*Vector of which cells the vertex is referenced in*/
	vector<vector<uint>> vertincells(nPts);

	for(uint ii = 0; ii < nElem; ++ii)
	{
		for(auto const& vert:cdata.elems[ii])
		{	/*Add count for every mention of the vertex, and push back the cell id*/
			vertincells[vert].emplace_back(ii);
		}
	}

	// uint numFaces = 0;
	// cout << vertmentions[0] << endl;
	// for(auto const& vert:vertincells[0])
	// 	cout << vert << "  ";
	
	// cout << endl;

	/*Create a colour vector to identify if a face has been identified*/
	vector<vector<uint>> colour(cdata.elems.size());
	for(uint ii = 0; ii < nElem; ++ii)
	{
		if(cdata.elems[ii].size() == 8)
		{
			colour[ii] = vector<uint>(6,0);
		}
		else if(cdata.elems[ii].size() == 6)
		{
			colour[ii] = vector<uint>(5,0);
		}
		else if(cdata.elems[ii].size() == 5)
		{
			colour[ii] = vector<uint>(5,0);
		}
		else if(cdata.elems[ii].size() == 4)
		{
			colour[ii] = vector<uint>(4,0);
		}
	}

	#ifdef DEBUG
	dbout << "Starting main loop..." << endl;
	#endif

	// uint cellSum=0;
	#pragma omp parallel shared(vertincells)
	{
		/*Create local copy of the face data*/
		FACE flocal;
		// uint cellCount = 0;
		// uint reported_Count = 0;
		// uint localCount = 0;

		#pragma omp for schedule(static) nowait 
		for(uint lindex = 0; lindex < nElem; ++lindex)
		{	
			vector<uint> lcell = cdata.elems[lindex];
			vector<vector<uint>> lfaces;
			uint lgeom=0;
			/*Check size of the inner vector, as that will show which geometry it is*/
			if(lcell.size() == 4)
			{	/*Tetraeder*/
				lgeom = 0;
				// fdata.nFacesPElem[lindex] = 4;
				lfaces = {{lcell[0],lcell[1],lcell[2]},
						  {lcell[1],lcell[0],lcell[3]},
				   		  {lcell[2],lcell[3],lcell[0]},
				   		  {lcell[3],lcell[1],lcell[2]}};
			}
			else if(lcell.size() == 5)
			{	/*Pyra*/
				lgeom = 1;
				// fdata.nFacesPElem[lindex] = 6;
				lfaces = {{lcell[0],lcell[3],lcell[2],lcell[1]},
				   		  {lcell[1],lcell[4],lcell[0]},
				   		  {lcell[2],lcell[4],lcell[1]},
				   		  {lcell[3],lcell[4],lcell[2]},
				   		  {lcell[4],lcell[3],lcell[0]}};						
			}
			else if(lcell.size() == 6)
			{	/*Prism*/
				lgeom = 2;
				// fdata.nFacesPElem[lindex] = 8;
				lfaces = {{lcell[1],lcell[2],lcell[5],lcell[4]},
				   		  {lcell[4],lcell[3],lcell[0],lcell[1]},
				   		  {lcell[3],lcell[5],lcell[2],lcell[0]},
				   		  {lcell[2],lcell[0],lcell[1]},
				   		  {lcell[5],lcell[3],lcell[4]}};
			}
			
			else if(lcell.size() == 8)
			{	/*Hexaeder*/
				lgeom = 3;
				// fdata.nFacesPElem[lindex] = 12;
				lfaces = {{lcell[0],lcell[1],lcell[2],lcell[3]},
						  {lcell[4],lcell[0],lcell[3],lcell[7]},
						  {lcell[1],lcell[0],lcell[4],lcell[5]},
						  {lcell[2],lcell[1],lcell[5],lcell[6]},
						  {lcell[3],lcell[2],lcell[6],lcell[7]},
						  {lcell[7],lcell[6],lcell[5],lcell[4]}};
			}

			CheckFaces(vertincells, lfaces, lindex, lgeom, cdata, colour, flocal);
			
// 			if (cellCount >= 25000)
// 		    {
// 				#pragma omp atomic
// 				cellSum += 25000;
// 				cellCount = 0;
// #ifdef DEBUG
// 				#pragma omp critical
// 				cout << "Thread " << omp_get_thread_num() << " processed cells: " << localCount << " (" << 
// 				100.0 * float(localCount)/(float(cdata.numElem))  << "%)" <<  endl; 
// #endif
// 		    }
// 		    else
// 		    {
// 				++cellCount;
// 		    }

// 		    localCount++;
// 		    // size_t tid = 0;
// 			size_t tid = omp_get_thread_num();
// 			if(tid == 0)
// 			{
// 				if(cellSum - reported_Count >= 100000)
// 				{
// 					cout << "Processed cells: " << cellSum << " (" << 
// 						100.0 * float(cellSum)/(float(cdata.numElem))  << "%)" <<  endl; 

// 					reported_Count = cellSum;

// 				}
// 			}

		}

		cout << "Thread " << omp_get_thread_num() << ": cell processing complete!" << endl;

		#pragma omp for schedule(static) ordered
		for(int ii=0; ii<omp_get_num_threads(); ii++)
		{
			#pragma omp ordered
			fdata.insert(flocal);
		}

		
		#pragma omp for schedule(static) nowait reduction(+:unmfaces)
		for (auto cell:colour)
			for (auto face:cell)
			{
				if(face!=1)
				{
					unmfaces++;
				}
			}
	}

	if(unmfaces != 0)
	cout << "There are " << unmfaces << " unmatched faces." << endl;

	if (fdata.faces.size() != fdata.numFaces)
	{
		cout << "numFaces is not being measured correctly." << endl;
		cout << "faces vector size: "<< fdata.faces.size() << "  numFaces: " << fdata.numFaces << endl;
	}

	if(fdata.celllr.size() != fdata.numFaces)
	{
		cout << "Not all faces have left and right cells identified." << endl;
		cout << "Celllr size: " << fdata.celllr.size();
		cout << "  N_Faces: " << fdata.numFaces << endl;
	}

	#ifdef DEBUG
	dbout << "Building faces complete. Number of faces: " << fdata.numFaces << endl;
	dbout << "Average number of faces per element:" << float(fdata.numFaces)/float(fdata.numElem) << endl;
	dbout << "Exiting BuildFaces..." << endl;
	#endif

	cout << "Building faces complete. Number of faces: " << fdata.numFaces << endl;
	cout << "Average number of faces per element:" << float(fdata.numFaces)/float(fdata.numElem) << endl;
	cout << "Number of wall faces: " << fdata.nWall << " Far field faces: " << fdata.nFar << endl;
}

void Write_Face_Data(string const& meshIn, FACE const& fdata)
{
	string meshOut = meshIn;
	meshOut.append(".faces");

	#ifdef DEBUG
	dbout << "Entering Write_Face_Data..." << endl;
	dbout << "Output file: " << meshOut << endl;
	#endif

	#ifdef DEBUG 
	cout << "Attempting write output file." << endl;
	cout << "File: " << meshOut << endl;
	#endif
	
	/* Create the file. */
	int retval;
	int meshID;
	if ((retval = nc_create(meshOut.c_str(), NC_CLOBBER, &meshID)))
	{	
		cout << "Error: Failed to open file \"" << meshOut << "\" when outputting." << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_put_att_text(meshID, NC_GLOBAL, "converted_mesh_filename", meshIn.length(), meshIn.c_str())))
	{
		cout << "Error: Failed to attach filename attribute" << endl;
		ERR(retval);
		exit(-1);
	}

	/* Define the dimensions. */
	int nElemID, nFaceID, nPpFcID, nWallID, nFarID, nPntsID;

	if ((retval = nc_def_dim(meshID, "no_of_elements", fdata.numElem, &nElemID)))
	{
		cout << "Error: Failed to define \"no_of_elements\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_faces", fdata.numFaces, &nFaceID)))
	{
		cout << "Error: Failed to define \"no_of_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "points_per_face", 3, &nPpFcID)))
	{
		cout << "Error: Failed to define \"points_per_face\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_wall_faces", fdata.nWall, &nWallID)))
	{
		cout << "Error: Failed to define \"no_of_wall_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_farfield_faces", fdata.nFar, &nFarID)))
	{
		cout << "Error: Failed to define \"no_of_farfield_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_points", fdata.numPoint, &nPntsID)))
	{
		cout << "Error: Failed to define \"no_of_points\"" << endl;
		ERR(retval);
		exit(-1);
	}

	/* Define the variables */
	int faceVarID, leftVarID, rightVarID, ptsxID, ptsyID, ptszID;
	
	/*Define the faces*/
	int dimIDs[] = {nFaceID,nPpFcID};
	if ((retval = nc_def_var(meshID, "points_of_element_faces", NC_INT, 2,
							 dimIDs, &faceVarID)))
	{
		cout << "Error: Failed to define \"points_of_element_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "left_element_of_faces", NC_INT, 1,
								&nFaceID, &leftVarID)))
	{
		cout << "Error: Failed to define \"left_element_of_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "right_element_of_faces", NC_INT, 1,
							 &nFaceID, &rightVarID)))
	{
		cout << "Error: Failed to define \"right_element_of_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}


	/*Define the points*/
	if ((retval = nc_def_var(meshID, "points_xc", NC_DOUBLE, 1, &nPntsID, &ptsxID)))
	{
		cout << "Error: Failed to define \"points_xc\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "points_yc", NC_DOUBLE, 1, &nPntsID, &ptsyID)))
	{
		cout << "Error: Failed to define \"points_yc\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "points_zc", NC_DOUBLE, 1, &nPntsID, &ptszID)))
	{
		cout << "Error: Failed to define \"points_zc\"" << endl;
		ERR(retval);
		exit(-1);
	}

	/* End define mode. */
	if ((retval = nc_enddef(meshID)))
		ERR(retval);

	/*Create the C array for the faces*/
	int* faces = new int[fdata.numFaces*3];
	for(uint ii = 0; ii < fdata.numFaces; ++ii)
		for(uint jj = 0; jj < 3; ++jj)
		{
			faces[index(ii,jj,3)] = static_cast<int>(fdata.faces[ii][jj]);
		}

	/*Put faces into the file*/
	size_t start[] = {0,0};
	size_t end[] = {fdata.numFaces,3};
	if ((retval = nc_put_vara_int(meshID, faceVarID, start, end, &faces[0])))
	{
		cout << "Failed to write face data" << endl;
		ERR(retval);
		exit(-1);
	}
	// elemFaces.putVar(faces);

	// /*State how many faces there are per element*/
	// int* nFacesPElem = new int[fdata.numElem];
	// for(uint ii = 0; ii < fdata.numElem; ++ii)
	// 	nFacesPElem[ii] = 0; 

	// for(uint ii = 0; ii < fdata.numFaces; ++ii)
	// {	/*Count the mentions of a cell*/
	// 	nFacesPElem[fdata.celllr[ii].first]++;
	// }

	// nElemFaces.putVar(fdata.nFacesPElem);

	/*Put face left and right into the file*/
	int* left = new int[fdata.numFaces];
	int* right = new int[fdata.numFaces];

	for(uint ii = 0; ii < fdata.numFaces; ++ii)
	{
		left[ii] = fdata.celllr[ii].first;
		right[ii] = fdata.celllr[ii].second;
	}


	if ((retval = nc_put_var_int(meshID, leftVarID, &left[0])))
	{
		cout << "Failed to write left element data" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_put_var_int(meshID, rightVarID, &right[0])))
	{
		cout << "Failed to write right element data" << endl;
		ERR(retval);
		exit(-1);
	}

	/*Create the C arrays for the vertices*/
	double* x = new double[fdata.numPoint];
	double* y = new double[fdata.numPoint];
	double* z = new double[fdata.numPoint];

	for(uint ii = 0; ii < fdata.numPoint; ++ii)
	{
		x[ii] = fdata.verts[ii](0);
		y[ii] = fdata.verts[ii](1);
		z[ii] = fdata.verts[ii](2);
	}

	/*Put them in the file*/
	if ((retval = nc_put_var_double(meshID, ptsxID, &x[0])))
	{
		cout << "Failed to write x coordinate." << endl;
		ERR(retval);
		exit(-1);
	}
		
	if ((retval = nc_put_var_double(meshID, ptsyID, &y[0])))
	{
		cout << "Failed to write y coordinate." << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_put_var_double(meshID, ptszID, &z[0])))
	{
		cout << "Failed to write z coordinate." << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_close(meshID)))
	{
		ERR(retval);
		exit(-1);
	}
}



void Write_ASCII_Face_Data(const FACE& fdata)
{
	#ifdef DEBUG
	dbout << "Entering Write_ASCII_Face_Data..." << endl;
	cout << "Attempting write output file." << endl;
	cout << "File: " << "Test.dat" << endl;
	#endif
	std::ofstream fout("Test.dat",std::ios::out);
	if(!fout.is_open())
	{
		cout << "Couldn't open the output file." << endl;
		exit(-1);
	}

	uint TotalNumFaceNodes = 0;
	for(uint ii = 0; ii < fdata.faces.size(); ++ii)
	{
		TotalNumFaceNodes += fdata.faces[ii].size();
	}

	fout << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
	fout << "ZONE T=\"FEPOLYHEDRON Test\"" << endl;
	fout << "ZONETYPE=FEPOLYHEDRON" << endl;
	fout << "NODES=" << fdata.numPoint << " ELEMENTS=" << fdata.numElem << " FACES=" << fdata.numFaces << endl;
	fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
	fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

	uint w = 15;
	uint preci = 6;
	fout << std::left << std::scientific << std::setprecision(preci);
	// fout << fdata.numElem << " " << fdata.numFaces << "  " <<  fdata.numPoint << endl;
	
	/*Write vertices in block format (Each dimension in turn)*/
	uint newl = 0;
	fout << std::setw(1);
	for(uint DIM = 0; DIM < SIMDIM; ++DIM)
	{
		for(uint ii = 0; ii < fdata.verts.size(); ++ii)
		{
			fout << std::setw(w) << fdata.verts[ii](DIM);
			newl++;

			if(newl>4)
			{
				fout << endl;
				fout << " ";
				newl=0;
			}
		}
	}
	fout << endl;
	

	fout << std::left << std::fixed;
	w = 9;
	/*Inform of how many vertices in each face*/
	fout << "#node count per face" << endl;
	newl = 0;
	for (uint ii = 0; ii < fdata.faces.size(); ++ii)
	{
		fout << std::setw(w) << fdata.faces[ii].size();
		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}
	fout << endl;
	/*Write the face data*/
	fout << "#face nodes" << endl;
	for (uint ii = 0; ii < fdata.faces.size(); ++ii)
	{
		for(auto const& vertex:fdata.faces[ii])
		{	/*Write face vertex indexes*/
			fout << std::setw(w) << vertex+1;
			if (vertex > fdata.numPoint)
			{
				cout << "Trying to write a vertex outside of the number of points." << endl;
			}
		}
		fout << endl;
	}

	/*Write face left and right*/
	newl = 0;
	fout << "#left elements" << endl;
	for (uint ii = 0; ii < fdata.celllr.size(); ++ii)
	{
		fout << std::setw(w) << fdata.celllr[ii].first+1 ;
		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}
	fout << endl;

	fout << "#right elements" << endl;
	newl = 0;
	for (uint ii = 0; ii < fdata.celllr.size(); ++ii)
	{
		if(fdata.celllr[ii].second < 0)
			fout<< std::setw(w) << 0 ;
		else
			fout << std::setw(w) << fdata.celllr[ii].second+1;


		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}

	fout.close();

	#ifdef DEBUG
	dbout << "Exiting Write_Face_Data..." << endl;
	#endif
}

void Write_Cell_Data(const CELL& cdata)
{
	#ifdef DEBUG
	dbout << "Entering Write_Cell_Data..." << endl;
	#endif

	cout << "Writing cell based data." << endl;

	std::ofstream fout("Cell.dat",std::ios::out);
	if (!fout.is_open())
	{
		cout << "Failed to open data file for writing mesh." << endl;
		exit(-1);
	}

	fout << "TITLE = \"3D TAU Solution\"\n";
	fout << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\"\n";
	fout << "ZONE T=\"Cell Data\"" << endl;
	fout << "N=" << cdata.numPoint << ", E=" << cdata.numElem << 
	", F=FEBLOCK, ET=BRICK"  << endl << endl;

	/*Write vertices*/
	fout << std::left << std::scientific << std::setprecision(6);
	fout << std::setw(1);
	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	
		uint kk = 0;
		for(uint jj = 0; jj < cdata.verts.size(); ++jj)
		{
			fout << std::setw(15) << cdata.verts[jj][ii];
			kk++;

			if(kk == 5)
			{
				fout << endl;
				fout << std::setw(1);
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fout << "\n";
	}

	/*Write element indexing*/
	fout << std::fixed;
	for(uint ii = 0; ii < cdata.elems.size(); ++ii)
	{	
		for(auto elem:cdata.elems[ii])
		{
			fout << std::setw(6) << elem+1;
		}
		fout << "\n";
	}

	fout.close();

	#ifdef DEBUG
	dbout << "Exiting Write_Cell_Data..." << endl;
	#endif

}

void Write_Griduns(const FACE& fdata)
{
	std::ofstream fout("griduns",std::ios::out);

	if(!fout.is_open())
	{
		cout << "Couldn't open griduns." << endl;
		exit(-1);
	}

	uint w = 12;
	fout << std::left << std::fixed;
	fout << setw(w) << fdata.numElem << setw(w) << 
		fdata.numFaces << setw(w) << fdata.numPoint << endl;

	for(uint ii = 0; ii < fdata.numFaces; ++ii)
	{
		for(auto const& vert:fdata.faces[ii])
			fout << setw(w) << vert;

		fout << setw(w) << fdata.celllr[ii].first << setw(w) << fdata.celllr[ii].second << endl;
	}

	fout << std::scientific << std::setprecision(6);
	w = 15;
	for(uint ii = 0; ii < fdata.numPoint; ++ii)
	{
		fout << setw(w) << ii;
		for(uint jj = 0; jj < SIMDIM; ++jj)
		{
			fout << setw(w) << fdata.verts[ii](jj);
		}
		fout << endl;
	}

	fout.close();
}

int main (int argc, char** argv)
{
	// omp_set_num_threads(NTHREADS);

	/*Idea: Take TAU cell based mesh, and convert to a face based data in NetCDF or TECIO*/
	string meshIn = argv[1];
	string bmapIn = argv[2];

	#ifdef DEBUG 
		cout << "Attempting read of NetCDF file." << endl;
		cout << "Mesh file: " << meshIn << endl;
		// cout << "Solution file: " << solIn << endl;
	#endif

	int meshID;
	int retval;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{	
		cout << "Failed to open mesh file \"" << meshIn << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading cell data..." << endl;

	int ptDimID, elemDimID;
	size_t nPnts, nElem;

		// Retrieve how many elements there are.
	Get_Dimension(meshID, "no_of_elements", elemDimID, nElem);
	Get_Dimension(meshID, "no_of_points", ptDimID, nPnts);

	CELL cdata(nElem, nPnts);

	cout << "nElem : " << nElem << " nPts: " << nPnts << endl;

	/*Retrieve the cdata*/
	int nTetDimID, nPpTDimID;
	size_t nTets, nPpTe;
	vector<vector<uint>> tets; 
	if ((retval = nc_inq_dimid(meshID, "no_of_tetraeders", &nTetDimID)))
	{
		cout << "No tetraeders in mesh file." << endl;
	}
	else
	{
		Get_Dimension(meshID, "no_of_tetraeders", nTetDimID, nTets);
		Get_Dimension(meshID, "points_per_tetraeder", nPpTDimID, nPpTe);
		tets = Get_Element(meshID, "points_of_tetraeders",nTets,nPpTe);
	}

	/* Retrieve prism data */
	int nPriDimID, nPpPrDimID;
	size_t nPris, nPpPr;
	vector<vector<uint>> prism;
	if ((retval = nc_inq_dimid(meshID, "no_of_prisms", &nPriDimID)))
	{
		cout << "No prisms in mesh file." << endl;
	}
	else
	{
		Get_Dimension(meshID, "no_of_prisms", nPriDimID, nPris);
		Get_Dimension(meshID, "points_per_prism", nPpPrDimID, nPpPr);
		prism = Get_Element(meshID, "points_of_prisms", nPris, nPpPr);
	}

	/* Get pyramid data */
	int nPyrDimID, nPpPyDimID;
	size_t nPyrs, nPpPy;
	vector<vector<uint>> pyra;
	if ((retval = nc_inq_dimid(meshID, "no_of_pyramids", &nPyrDimID)))
	{
		cout << "No pyramids in mesh file." << endl;
	}
	else
	{
		Get_Dimension(meshID, "no_of_pyramids", nPyrDimID, nPyrs);
		Get_Dimension(meshID, "points_per_pyramid", nPpPyDimID, nPpPy);
		pyra = Get_Element(meshID, "points_of_pyramids", nPyrs, nPpPy);
	}

	/* Get hexahedral data */
	int nHexDimID, nPpHeDimID;
	size_t nHexs, nPpHe;
	vector<vector<uint>> hexa;
	if ((retval = nc_inq_dimid(meshID, "no_of_hexaeders", &nHexDimID)))
	{
		cout << "No hexaeders in mesh file." << endl;
	}
	else
	{
		Get_Dimension(meshID, "no_of_hexaeders", nHexDimID, nHexs);
		Get_Dimension(meshID, "points_per_hexaeder", nPpHeDimID, nPpHe);
		hexa = Get_Element(meshID, "points_of_hexaeders", nHexs, nPpHe);
	}
	
	/*Get the coordinates of the mesh*/
	cdata.verts = Get_Coordinates(meshID, nPnts);

	if (cdata.verts.size() != nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Get surface faces*/
	vector<int> markers = Find_Bmap_Markers(bmapIn);
	cdata.sfaces = Get_Surface(meshID, markers);	

	/*Put data into cell structure, and generate face based data*/
	if(tets.size()!=0)
		cdata.elems.insert(cdata.elems.end(),tets.begin(),tets.end());

	if(prism.size()!=0)
		cdata.elems.insert(cdata.elems.end(),prism.begin(),prism.end());

	if(pyra.size()!=0)
		cdata.elems.insert(cdata.elems.end(),pyra.begin(),pyra.end());


	if(hexa.size()!=0)
		cdata.elems.insert(cdata.elems.end(),hexa.begin(),hexa.end());


	if(cdata.elems.size()!=nElem)
		cout << "Some data has been missed.\nPlease check which kinds of volumes used." << endl;
	else
		cout << "All element data ingested" << endl;


	// Write_Cell_Data(cdata);

	/*Now build the face based data*/
	cout << "Building the face-based data..." << endl;
	FACE fdata(cdata);
	BuildFaces(cdata,fdata);

	Write_Face_Data(meshIn, fdata);
	// Write_Griduns(fdata);
	return 0;
}