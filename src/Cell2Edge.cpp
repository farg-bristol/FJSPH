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
	std::ofstream dbout("Cell2Edge.log",std::ios::out);
#endif

// #ifndef NTHREADS
// #define NTHREADS 4
// #endif

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 2
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
typedef struct FACE
{
	/*Standard contructor*/
	FACE(const uint nElem, const uint nPts)
	{
		numElem = nElem;
		numPoint = nPts;
		verts.reserve(numPoint);
	}
	
	/*Zone info*/
	uint numPoint, numElem;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Surface faces*/
	vector<vector<uint>> faces;
	vector<vector<uint>> internal;

	/*Face markers*/
	// vector<int> markers;
}FACE;

typedef class EDGE
{
	public:
	EDGE(const FACE& fdata): numPoint(fdata.numPoint)
	{
		verts = fdata.verts;
		numElem = fdata.numElem;
		// numPoint = fdata.numPoint;
		numEdges = 0; nFar = 0; nWall = 0;
		// nFacesPElem = new int[fdata.numElem];
	}	

	EDGE() : numPoint(0)
	{
		numEdges = 0; nFar = 0; nWall = 0;
	};

	void insert(const EDGE& elocal)
	{
		edges.insert(edges.end(),elocal.edges.begin(),elocal.edges.end());
		wall.insert(wall.end(),elocal.wall.begin(),elocal.wall.end());
		far.insert(far.end(),elocal.far.begin(),elocal.far.end());
		celllr.insert(celllr.end(),elocal.celllr.begin(),elocal.celllr.end());
		numEdges += elocal.numEdges;
		nFar += elocal.nFar;
		nWall += elocal.nWall;
	}

	uint numElem, numEdges, nFar, nWall;
	vector<StateVecD> verts;
	vector<vector<uint>> edges; /*edge indexes*/
	vector<std::pair<int,int>> celllr; /*Cell left and right of the face*/
	vector<vector<uint>> wall, far;
	vector<uint> usedVerts;
	// int* nFacesPElem;

	uint numPoint;
}EDGE;

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

vector<int> Find_Bmap_Markers(const string& bmapIn, int& symPlane1, int& symPlane2)
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
	uint foundSymP = 0;

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
					std::stringstream sstr(line);
					string temp;
					int found;

					while(!sstr.eof())
					{
						sstr >> temp;
						if(std::stringstream(temp) >> found)
						{
							if(foundSymP == 0)
							{
								symPlane1 = found;
								foundSymP++;
							}
							else
							{
								symPlane2 = found;
							}
							
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
		else if(line.find("Type: euler wall")!=string::npos || 
			line.find("Type: viscous wall")!=string::npos ||
			line.find("Type: sharp edge")!=string::npos)
		{
			
			/*This marker is an internal edge so store it.*/
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

	cout << "Symmetry plane markers: " << symPlane1 << "  " << symPlane2 << endl;

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
		if(coordX[ii] == 0 || coordX[ii] == -1 || coordX[ii] == 1)
		{
			counts[0]++;
		}

		if(coordY[ii] == 0 || coordY[ii] == -1 || coordY[ii] == 1)
		{
			counts[1]++;
		}

		if(coordZ[ii] == 0 || coordZ[ii] == -1 || coordZ[ii] == 1)
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

void Get_Data(NcFile& fin, FACE& fdata, vector<int> markers, int symPlane)
{
	#ifdef DEBUG
	dbout << "Entering Get_Data..." << endl;
	#endif

	/*Faces*/
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
		cout << "Mismatch of number of surface markers and faces ingested." << endl;
		cout << "Number of surface markers: " << nMarkers << 
		"  Number in vector: " << faceVec.size() << endl;
	}

	int* faceMarkers = new int[nMarkers];

	surfaceMarkers.getVar(faceMarkers);
	cout << "symPlane: " << symPlane << endl;

	/*Create the faces of the symmetry plane*/
	for(uint ii = 0; ii < nMarkers; ++ii)
	{
		if(faceMarkers[ii] == symPlane)
		{	/*The face is part of the symmetry plane*/
			
			fdata.faces.emplace_back(faceVec[ii]);
		}
		else if(std::find(markers.begin(),markers.end(),faceMarkers[ii])!= markers.end())
		{
			/*Its an internal face*/
			fdata.internal.emplace_back(faceVec[ii]);
		}
	}

	cout << "Cell faces: " << fdata.faces.size();
	cout << "   Internal faces: " << fdata.internal.size() << endl;

#ifdef DEBUG
	dbout << "Cell faces: " << fdata.faces.size();
	dbout << "   Internal faces: " << fdata.internal.size() << endl;
#endif

	/*Get the coordinates*/
	fdata.verts = Get_Coordinates(fin);

	#ifdef DEBUG
	dbout << "Exiting Get_Data..." << endl;
	#endif
}


void Add_Edge(vector<uint> const& edge, std::pair<int,int> const& leftright,
	EDGE& edata)
{
	edata.celllr.emplace_back(leftright);
	edata.edges.emplace_back(edge);
	edata.numEdges++;
}


void BuildEdges(const FACE& fdata, EDGE& edata)
{
	#ifdef DEBUG
	dbout << "Entering BuildFaces..." << endl;
	#endif
	uint nFace = fdata.faces.size();
	uint nPts = fdata.numPoint;
	uint unmfaces=0;
	uint nEdges = 0;


	/*Vector of which cells the vertex is referenced in*/
	vector<vector<uint>> vertincells(nPts);

	for(uint ii = 0; ii < nFace; ++ii)
	{
		for(auto const& vert:fdata.faces[ii])
		{	/*Add count for every mention of the vertex, and push back the cell id*/
			vertincells[vert].emplace_back(ii);
		}
	}

	/*Create a colour vector to identify if a face has been identified*/

	vector<vector<uint>> colour(nFace);
	for(uint ii = 0; ii < nFace; ++ii)
	{
		nEdges+=fdata.faces[ii].size();
		colour[ii] = vector<uint>(fdata.faces[ii].size(),0);
	}

	#ifdef DEBUG
	dbout << "Starting main loop..." << endl;
	dbout << "Number of faces: " << fdata.faces.size() << endl;
	#endif

	// uint cellSum=0;
	#pragma omp parallel shared(vertincells)
	{
		/*Create local copy of the face data*/
		EDGE elocal;
		// uint cellCount = 0;
		// uint reported_Count = 0;
		// uint localCount = 0;
		
		#pragma omp for schedule(static) nowait 
		for(uint lindex = 0; lindex < nFace; ++lindex)
		{	
			const vector<uint> lcell = fdata.faces[lindex];
	
			uint numl = lcell.size();
			uint l0 = lcell[numl-1];
			uint lfaceindex = 0;
			
			for(auto const& l1:lcell)
			{	/*Define that the cell of the top-level cell search is on the 'left'*/
				
				for (auto const& rindex:vertincells[l1])
				{	
					/*If the cell is the current cell, ignore*/
					if(rindex==lindex)
						continue;	

					/*Define that the cell of the inner search is on the 'right'*/
					const vector<uint> rcell = fdata.faces[rindex];

					uint numr = rcell.size();
					uint r0 = rcell[numr-1];
					uint rfaceindex = 0;

					for(auto const& r1:rcell)
					{
						if((r0 == l0 && r1 == l1) || (r1 == l0 && r0 == l1))
						{
							/*Then the edge is a match*/
							if(colour[lindex][lfaceindex] == 0 && colour[rindex][rfaceindex] == 0)
							{
								std::pair<int,int> leftright(lindex,rindex);
								#pragma omp atomic
								colour[lindex][lfaceindex]++;
								#pragma omp atomic
								colour[rindex][rfaceindex]++;

								vector<uint> edge = {l0,l1};
								Add_Edge(edge,leftright,elocal);
								goto matchfound;
							}
						}
						r0 = r1;
						rfaceindex++;
					}	
				}
				
				
				/*If a match has not been found, the face must be a boundary face*/	
				if(colour[lindex][lfaceindex] == 0)
				{			
					for (auto sface:fdata.internal)
					{	/*Search through the surface faces to identify */
						/*if the face is an internal face or not*/
						int a = (std::find(sface.begin(),sface.end(),l0) != sface.end());
						int b = (std::find(sface.begin(),sface.end(),l1) != sface.end());

						if(a && b)
						{	/*edge is an internal face*/
							
							std::pair<int,int> leftright(lindex,-1);
							#pragma omp atomic
								colour[lindex][lfaceindex]++;
							vector<uint> edge = {l0,l1};
							Add_Edge(edge,leftright,elocal);
							elocal.wall.emplace_back(edge);
							elocal.nWall++;
							goto matchfound;
						}
					}

					/*If still uncoloured, then face is an external boundary*/
					std::pair<int,int> leftright(lindex,-2);
					#pragma omp atomic
					colour[lindex][lfaceindex]++;
					vector<uint> edge = {l0,l1};
					Add_Edge(edge,leftright,elocal);
					elocal.far.emplace_back(edge);
					elocal.nFar++;
				}

		matchfound:	
				l0 = l1;		
				lfaceindex++;
			}
			
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

		

		#pragma omp for schedule(static) ordered
		for(int ii=0; ii<omp_get_num_threads(); ii++)
		{
			#pragma omp ordered
			edata.insert(elocal);
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

	if (edata.edges.size() != edata.numEdges)
	{
		cout << "Not all edges have been found." << endl;
		cout << "Edges vector size: "<< edata.edges.size() << "  nEdges: " << edata.numEdges << endl;
	}

	if(edata.celllr.size() != edata.numEdges)
	{
		cout << "Not all edges have left and right cells identified." << endl;
		cout << "Celllr size: " << edata.celllr.size();
		cout << "  nEdges: " << edata.numEdges << endl;
	}

	#ifdef DEBUG
	dbout << "Building edges complete. Number of edges: " << edata.numEdges << endl;
	dbout << "Average number of edges per element: " << float(edata.numEdges)/float(nFace) << endl;
	dbout << "Number of wall edges: " << edata.nWall << " Far field edges: " << edata.nFar << endl;
	dbout << "Exiting BuildFaces..." << endl;
	#endif

	cout << "Building edges complete. Number of edges: " << edata.numEdges << endl;
	cout << "Average number of edges per element: " << float(edata.numEdges)/float(nFace) << endl;
	cout << "Number of wall edges: " << edata.nWall << " Far field edges: " << edata.nFar << endl;
	// cout << edata.wall.size() << "  " << edata.far.size() << endl;
}

/*Recast data to ignore unused vertices*/
void Recast_Data(EDGE& edata)
{
#ifdef DEBUG
	dbout << "Entering recast data..." << endl;
#endif
	vector<uint> vertInUse;
	for(auto const& edge:edata.edges)
	{
		vertInUse.emplace_back(edge[0]);
		vertInUse.emplace_back(edge[1]);
	}

	/*Sort and erase duplicates*/
	std::sort(vertInUse.begin(),vertInUse.end());
	vertInUse.erase(std::unique(vertInUse.begin(),vertInUse.end()),vertInUse.end());

#ifdef DEBUG
	dbout << "vertInUse size: " << vertInUse.size() << endl;
#endif
	cout << "vertInUse size: " << vertInUse.size() << endl;
	/*Now the edges need recasting to the index of the vector*/
	vector<StateVecD> newVerts(vertInUse.size());
	#pragma omp parallel for schedule(static)
	for(uint ii = 0; ii < vertInUse.size(); ++ii)
	{
		for(auto& edge:edata.edges)
		{
			if(edge[0] == vertInUse[ii])
				edge[0] = ii;
			
			if(edge[1] == vertInUse[ii])
				edge[1] = ii;
		}
		newVerts[ii] = edata.verts[ii];
	}
	edata.usedVerts = vertInUse;
	edata.verts = newVerts;
	edata.numPoint = newVerts.size();

#ifdef DEBUG
	dbout << "Exiting recast data..." << endl;
#endif
}

void Write_Edge_Data(const string& meshIn, const EDGE& edata)
{
	string meshOut = meshIn;
	meshOut.append(".edges");

	#ifdef DEBUG
	dbout << "Entering Write_Edge_Data..." << endl;
	dbout << "Output file: " << meshOut << endl;
	#endif

	#ifdef DEBUG 
	cout << "Attempting write output file." << endl;
	cout << "File: " << meshOut << endl;
	#endif

	NcFile fout(meshOut, NcFile::replace);

	/*Dimensions needed*/
	NcDim nElems = fout.addDim("no_of_elements",edata.numElem);
	NcDim nEdges = fout.addDim("no_of_edges",edata.numEdges);
	NcDim ppFace = fout.addDim("points_per_edge",2);
	NcDim nWall = fout.addDim("no_of_wall_edges",edata.nWall);
	NcDim nFar = fout.addDim("no_of_farfield_edges",edata.nFar);
	NcDim nPoint = fout.addDim("no_of_points",edata.numPoint);
	
	/*Define the faces*/
	vector<NcDim> faceVar;
	faceVar.emplace_back(nEdges);
	faceVar.emplace_back(ppFace);
	NcVar elemEdges = fout.addVar("points_of_element_edges",ncInt,faceVar);
	NcVar leftElems = fout.addVar("left_element_of_edges",ncInt,nEdges);
	NcVar rightElems = fout.addVar("right_element_of_edges",ncInt,nEdges);
	NcVar usedVerts = fout.addVar("vertices_in_use",ncInt,nPoint);

	/*Define the points*/
	NcVar vertsX = fout.addVar("points_xc",ncDouble,nPoint);
	NcVar vertsZ = fout.addVar("points_zc",ncDouble,nPoint);

	/*Create the C array for the faces*/
	int* edges = new int[edata.numEdges*2];
	for(uint ii = 0; ii < edata.numEdges; ++ii)
		for(uint jj = 0; jj < 2; ++jj)
		{
			edges[index(ii,jj,2)] = static_cast<int>(edata.edges[ii][jj]);
		}

	/*Put edges into the file*/
	elemEdges.putVar(edges);

	/*Put face left and right into the file*/
	int* left = new int[edata.numEdges];
	int* right = new int[edata.numEdges];

	for(uint ii = 0; ii < edata.numEdges; ++ii)
	{
		left[ii] = edata.celllr[ii].first;
		right[ii] = edata.celllr[ii].second;
	}

	leftElems.putVar(left);
	rightElems.putVar(right);

	/*Write which vertices are in use, for reading the solution file*/
	int* uVert = new int[edata.numPoint];

	for(size_t ii = 0; ii < edata.numPoint; ++ii)
	{
		uVert[ii] = edata.usedVerts[ii]; 
	}

	usedVerts.putVar(uVert);

	/*Create the C arrays for the vertices*/
	double* x = new double[edata.numPoint];
	double* z = new double[edata.numPoint];

	for(uint ii = 0; ii < edata.numPoint; ++ii)
	{
		x[ii] = edata.verts[ii](0);
		z[ii] = edata.verts[ii](1);
	}

	/*Put them in the file*/
	vertsX.putVar(x);
	vertsZ.putVar(z);

}

void Write_Edge_Tecplot(const EDGE& edata)
{
	#ifdef DEBUG
	dbout << "Entering Write_Edge_Tecplot..." << endl;
	cout << "Attempting write output file." << endl;
	cout << "File: " << "griduns.dat" << endl;
	#endif
	std::ofstream fout("griduns.dat",std::ios::out);
	if(!fout.is_open())
	{
		cout << "Couldn't open the output file." << endl;
		exit(-1);
	}


	fout << "VARIABLES= \"X\" \"Z\" " << endl;
	fout << "ZONE T=\"FEPOLYGON Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "NODES=" << edata.numPoint << " ELEMENTS=" << edata.numElem << " FACES=" << edata.numEdges << endl;
	fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

	uint w = 15;
	uint preci = 6;
	fout << std::left << std::scientific << std::setprecision(preci);
	
	/*Write vertices in block format (Each dimension in turn)*/
	uint newl = 0;
	fout << std::setw(1);
	for(uint DIM = 0; DIM < SIMDIM; ++DIM)
	{
		for(uint ii = 0; ii < edata.verts.size(); ++ii)
		{
			fout << std::setw(w) << edata.verts[ii](DIM);
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
	// /*Inform of how many vertices in each face*/
	// fout << "#node count per face" << endl;
	// newl = 0;
	// for (uint ii = 0; ii < edata.edges.size(); ++ii)
	// {
	// 	fout << std::setw(w) << edata.edges[ii].size();
	// 	newl++;

	// 	if(newl>4)
	// 	{
	// 		fout << endl;
	// 		newl=0;
	// 	}
	// }
	// fout << endl;

	/*Write the face data*/
	fout << "#face nodes" << endl;
	for (uint ii = 0; ii < edata.edges.size(); ++ii)
	{
		for(auto const& vertex:edata.edges[ii])
		{	/*Write face vertex indexes*/
			fout << std::setw(w) << vertex+1;
			if (vertex > edata.numPoint)
			{
				cout << "Trying to write a vertex outside of the number of points." << endl;
			}
		}
		fout << endl;
	}

	/*Write face left and right*/
	newl = 0;
	fout << "#left elements" << endl;
	for (uint ii = 0; ii < edata.celllr.size(); ++ii)
	{
		fout << std::setw(w) << edata.celllr[ii].first+1 ;
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
	for (uint ii = 0; ii < edata.celllr.size(); ++ii)
	{
		if(edata.celllr[ii].second < 0)
			fout<< std::setw(w) << 0 ;
		else
			fout << std::setw(w) << edata.celllr[ii].second+1;


		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}

	fout.close();

	#ifdef DEBUG
	dbout << "Exiting Write_Edge_Tecplot..." << endl;
	#endif
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
	NcDim faceNo = fin.getDim("no_of_surfaceelements");
	uint nFace = static_cast<uint>(elemNo.getSize());
	uint nPts = static_cast<uint>(pointNo.getSize());
	uint nFcs = static_cast<uint>(faceNo.getSize());

	FACE fdata(nFace,nPts);

	cout << "nFace : " << nFace << " nPts: " << nPts << " nFcs: " << nFcs << endl;
	
	/*Get surface faces*/
	int symPlane1 = 0;
	int symPlane2 = 0;
	vector<int> markers = Find_Bmap_Markers(bmapIn,symPlane1, symPlane2);
	Get_Data(fin,fdata,markers,symPlane1);

	cout << "Face size: " << fdata.faces.size() << endl;

	EDGE edata(fdata);

	cout << "Building edge based data..." << endl;
	BuildEdges(fdata,edata);

	Recast_Data(edata);


	Write_Edge_Tecplot(edata);
	Write_Edge_Data(meshIn,edata);
	return 0;
}