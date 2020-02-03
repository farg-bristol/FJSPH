/*Cell based to Face based data converter*/

#include <vector>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2
#include <fstream>
#include <iomanip>

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

vector<int> Find_Bmap_Markers(const string& bmapIn)
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
		if(line.find("Type: euler wall")!=string::npos || 
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

vector<vector<uint>> Get_Surface(NcFile& fin, const vector<int>& markers)
{
	#ifdef DEBUG
	dbout << "Entering Get_Surface..." << endl;
	#endif
	vector<vector<uint>> faceVec;

	NcVar sTria = fin.getVar("points_of_surfacetriangles");	
	if(sTria.isNull()) 
	{
		cout << "No triangle surfaces." << endl;
		#ifdef DEBUG
		dbout  << "No triangle surfaces." << endl;
		#endif
	}
	else
	{
		NcDim nTriD = sTria.getDim(0);
		uint nTri = static_cast<uint>(nTriD.getSize());

		#ifdef DEBUG
			dbout << "Number of triangles: " << nTri << endl;
		#endif

		int* tris = new int[nTri*3];

		/*Get the actual data from the file*/
		vector<size_t> startp,countp;
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(nTri);
		countp.push_back(3);

		#ifdef DEBUG
		dbout << "Attempting to read NetCDF surface triangles." << endl;
		#endif

		sTria.getVar(startp,countp,tris);

		#ifdef DEBUG
		dbout << "Putting surface faces into a vector." << endl;
		#endif

		vector<vector<uint>> localVec(nTri,vector<uint>());
		for (uint ii = 0; ii < nTri; ++ii)
		{
			for(uint jj = 0; jj < 3; ++jj)
				localVec[ii].emplace_back(static_cast<uint>(tris[index(ii,jj,3)]));
		}

		faceVec.insert(faceVec.end(),localVec.begin(),localVec.end());
	}

	NcVar sQuad = fin.getVar("points_of_surfacequadrilaterals");	
	if(sQuad.isNull()) 
	{
		cout << "No quadrilateral surfaces." << endl;
		#ifdef DEBUG
		dbout  << "No quadrilateral surfaces." << endl;
		#endif


	}
	else
	{
		NcDim nQuadD = sQuad.getDim(0);
		uint nQuad = static_cast<uint>(nQuadD.getSize());
		
		#ifdef DEBUG
			dbout << "Number of quadrilaterals: " << nQuad << endl;
		#endif

		/*Allocate on the heap (can be big datasets)*/
		
		int* quads = new int[nQuad*4];

		vector<size_t> startp,countp;
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(nQuad);
		countp.push_back(4);

		#ifdef DEBUG
		dbout << "Attempting to read NetCDF surface quadrilaterals." << endl;
		#endif

		sQuad.getVar(startp,countp,quads);

		#ifdef DEBUG
		dbout << "Putting surface faces into a vector." << endl;
		#endif
		/*Convert it to a vector to store*/
		vector<vector<uint>> localVec(nQuad,vector<uint>());

		for (uint ii = 0; ii < nQuad; ++ii)
		{
			for(uint jj = 0; jj < 4; ++jj)
				localVec[ii].emplace_back(static_cast<uint>(quads[index(ii,jj,4)]));
		}

		faceVec.insert(faceVec.end(),localVec.begin(),localVec.end());
	}

	/*Get the boundarymarkers*/
	NcVar surfaceMarkers = fin.getVar("boundarymarker_of_surfaces");
	if(surfaceMarkers.isNull())
	{
		cout << "No data available on surface markers..." << endl;
	}
	
	NcDim nMarkersD = surfaceMarkers.getDim(0);
	uint nMarkers = static_cast<uint>(nMarkersD.getSize());
	if(faceVec.size() != nMarkers)
	{
		cout << "Mismatch of number of surface elements defined and number ingested." << endl;
		cout << "Number of surface elements: " << nMarkers << 
		"  Number in vector: " << faceVec.size() << endl;
	}

	int* faceMarkers = new int[nMarkers];

	surfaceMarkers.getVar(faceMarkers);

	vector<vector<uint>> farVec;

	/*Want to find which surfaces are the ones I want to look for*/
	for(uint ii = 0; ii < nMarkers; ++ii)
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
		cout << "Number of elements: " << nElem << endl;

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
	uint nPts = static_cast<uint>(dim.getSize());
	
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
	for (uint ii = 0; ii < nPts; ++ii)
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

void Write_Face_Data(const string& meshIn, const FACE& fdata)
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

	NcFile fout(meshOut, NcFile::replace);

	/*Dimensions needed*/
	NcDim nElems = fout.addDim("no_of_elements",fdata.numElem);
	NcDim nFaces = fout.addDim("no_of_faces",fdata.numFaces);
	NcDim ppFace = fout.addDim("points_per_face",3);
	NcDim nWall = fout.addDim("no_of_wall_faces",fdata.nWall);
	NcDim nFar = fout.addDim("no_of_farfield_faces",fdata.nFar);
	NcDim nPoint = fout.addDim("no_of_points",fdata.numPoint);
	
	/*Define the faces*/
	vector<NcDim> faceVar;
	faceVar.emplace_back(nFaces);
	faceVar.emplace_back(ppFace);
	NcVar elemFaces = fout.addVar("points_of_element_faces",ncInt,faceVar);
	// NcVar nElemFaces = fout.addVar("faces_per_element",ncInt,nElems);
	NcVar leftElems = fout.addVar("left_element_of_faces",ncInt,nFaces);
	NcVar rightElems = fout.addVar("right_element_of_faces",ncInt,nFaces);

	/*Define the points*/
	NcVar vertsX = fout.addVar("points_xc",ncDouble,nPoint);
	NcVar vertsY = fout.addVar("points_yc",ncDouble,nPoint);
	NcVar vertsZ = fout.addVar("points_zc",ncDouble,nPoint);

	/*Create the C array for the faces*/
	int* faces = new int[fdata.numFaces*3];
	for(uint ii = 0; ii < fdata.numFaces; ++ii)
		for(uint jj = 0; jj < 3; ++jj)
		{
			faces[index(ii,jj,3)] = static_cast<int>(fdata.faces[ii][jj]);
		}

	/*Put faces into the file*/
	elemFaces.putVar(faces);

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

	leftElems.putVar(left);
	rightElems.putVar(right);

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
	vertsX.putVar(x);
	vertsY.putVar(y);
	vertsZ.putVar(z);

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

	NcFile fin(meshIn, NcFile::read);
	cout << "Mesh file open. Reading cell data..." << endl;	

	NcDim elemNo = fin.getDim("no_of_elements");
	NcDim pointNo = fin.getDim("no_of_points");
	uint nElem = static_cast<uint>(elemNo.getSize());
	uint nPts = static_cast<uint>(pointNo.getSize());

	CELL cdata(nElem,nPts);

	cout << "nElem : " << nElem << " nPts: " << nPts << endl;
	

	/*Retrieve the cdata*/
	vector<vector<uint>> tets = Get_Element(fin,"points_of_tetraeders");
	vector<vector<uint>> prism = Get_Element(fin,"points_of_prisms");
	vector<vector<uint>> pyra = Get_Element(fin,"points_of_pyramids");
	vector<vector<uint>> hex = Get_Element(fin,"points_of_hexaeders");
	
	
	if(tets.size() == 0)
		cout << "No tetraeders in mesh file." << endl;

	if(prism.size() == 0)
		cout << "No prisms in mesh file." << endl;

	if(pyra.size() == 0)
		cout << "No pyramids in mesh file." << endl;

	if(hex.size() == 0)
		cout << "No hexaeders in mesh file." << endl;
	

	/*Get the coordinates of the mesh*/
	cdata.verts = Get_Coordinates(fin);

	if(cdata.verts.size()!= nPts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Get surface faces*/
	vector<int> markers = Find_Bmap_Markers(bmapIn);
	cdata.sfaces = Get_Surface(fin, markers);	

	/*Put data into cell structure, and generate face based data*/
	if(tets.size()!=0)
		cdata.elems.insert(cdata.elems.end(),tets.begin(),tets.end());

	if(prism.size()!=0)
		cdata.elems.insert(cdata.elems.end(),prism.begin(),prism.end());

	if(pyra.size()!=0)
		cdata.elems.insert(cdata.elems.end(),pyra.begin(),pyra.end());


	if(hex.size()!=0)
		cdata.elems.insert(cdata.elems.end(),hex.begin(),hex.end());


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