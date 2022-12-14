#include "Convert.h"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <TECIO.h>

/*Define Simulation Dimension*/
#define SIMDIM 2

using namespace netCDF;

/****** Eigen vector definitions ************/
typedef std::array<real,2> StateVecD;
typedef std::array<int,2> StateVecI;

typedef struct FACE
{
	/*Standard contructor*/
	FACE(uint const nElem, uint const nPts) : nPnts(nPts), nElem(nElem)
	{
		verts.reserve(nPts);
	}
	
	/*Zone info*/
	uint nPnts, nElem;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Surface faces*/
	vector<vector<size_t>> faces;
	vector<vector<size_t>> sfaces;

	/*Face markers*/
	vector<int> smarkers;
}FACE;

typedef class EDGE
{
	public:
	EDGE(const FACE& fdata) : nPnts(fdata.nPnts), nEdge(0), nSurf(0)
	{
		verts = fdata.verts;
		nElem = fdata.nElem;
	}	

	EDGE() : nPnts(0), nEdge(0), nSurf(0) {}

	void insert(const EDGE& elocal)
	{
		edges.insert(edges.end(),elocal.edges.begin(),elocal.edges.end());
		celllr.insert(celllr.end(),elocal.celllr.begin(),elocal.celllr.end());
		smarkers.insert(smarkers.end(),elocal.smarkers.begin(),elocal.smarkers.end());
		nEdge += elocal.nEdge;
		nSurf += elocal.nSurf;
		// nFar += elocal.nFar;
		// nWall += elocal.nWall;
	}

	size_t nPnts, nElem, nEdge, nSurf;
	vector<StateVecD> verts;
	vector<std::pair<size_t,size_t>> edges; /*edge indexes*/
	vector<std::pair<int,int>> celllr; /*Cell left and right of the face*/
	vector<int> smarkers;
	vector<uint> usedVerts;
	// int* nFacesPElem;

}EDGE;


std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(" \n\r\t\f\v");
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(" \n\r\t\f\v");
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string Get_Parameter_Value(string const& line)
{
    size_t pos = line.find(":");

    string value = line.substr(pos + 1);
    return ltrim(rtrim(value));
}

template<typename T>
void Get_Number(string const& line, string const& param, T &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        iss >> value;
    }
}

int Distance(vector<StateVecD> const &verts, std::pair<size_t, size_t> const &edge, size_t const &point)
{
    real  ty, tx;
    StateVecD vtx0, vtx1;

    tx = verts[point][0];
    ty = verts[point][1];

    vtx0 = verts[edge.first];
    vtx1 = verts[edge.second];

    real dist = ((vtx1[0]-vtx0[0])*(ty-vtx0[1])-(tx-vtx0[0])*(vtx1[1]-vtx0[0]))
    		/sqrt(pow((vtx1[0]-vtx0[0]),2) + pow((vtx1[1]-vtx0[1]),2));

    if(dist > 0)
    	return 1;
    else
    	return 0;
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

		size_t end = line.find_first_of('#');
		if(end != std::string::npos)
			line = line.substr(0,end+1);

		// cout << line << endl;
		if(line.find("Type")!=string::npos && line.find("symmetry plane")!=string::npos) 
		{
			string line2;
			/*This marker is a far field, so store it.*/
			/*Go to the start of the block, and search for the marker ID*/
			GotoLine(fin,blockstart);
			while(getline(fin,line2))
			{	
				size_t end = line2.find_first_of('#');
				if(end != std::string::npos)
					line2 = line2.substr(0,end+1);
				// cout << "inner:\t" << line << endl;

				if(line2.find("Markers")!=string::npos)
				{
					// cout << "Found a boundary marker" << endl;
					string temp = Get_Parameter_Value(line2);
					std::istringstream sstr(temp);
					int found;
					sstr >> found;

					if(foundSymP == 0)
					{
						symPlane1 = found;
						foundSymP++;
					}
					else
					{
						symPlane2 = found;
					}

					/*Go back to where we were*/
					blockstart = lineno+2;
					GotoLine(fin,lineno+2);
					break;
				}
			}
		}
		else if(line.find("Type")!=string::npos )
		{
			string line2;
			/*This marker is not a symmetry plane so store it.*/
			/*Go to the start of the block, and search for the marker ID*/
			GotoLine(fin,blockstart);
			while(getline(fin,line2))
			{	
				size_t end = line2.find_first_of('#');
				if(end != std::string::npos)
					line2 = line2.substr(0,end+1);
				// cout << "inner:\t" << line << endl;
				if(line2.find("Markers")!=string::npos)
				{
					string temp = Get_Parameter_Value(line2);
					std::istringstream sstr(temp);
					int found;
					sstr >> found;
					markers.emplace_back(found);
					
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

	cout << endl;

	#ifdef DEBUG
	dbout << "Exiting Find_Bmap_Markers..." << endl;
	#endif

	return markers;
}

/*To run on the mesh file*/
vector<std::array<real,2>> Get_Coordinates(int& fin, size_t const& nPnts, vector<uint> const& ptIndex)
{
	#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
	#endif

	vector<std::array<real,3>> coords = Get_Coordinate_Vector(fin, nPnts);
	/*Need to find which coordinate to ignore*/
	#ifdef DEBUG
		dbout << "Checking which dimension to ignore" << endl;
	#endif
	double tolerance = 1e-6;

	uint counts1[3] = {0};
	uint counts2[3] = {0};
	uint counts3[3] = {0};
	std::array<real,3> point1(coords[0]);

	#pragma omp parallel for
	for(auto const ii : ptIndex)
	{
		/* x plane */
		if (abs(coords[ii][0]) < tolerance) 
		{
			#pragma omp atomic
			counts1[0]++;
		}

		if(abs(abs(coords[ii][0])-1.0) < tolerance)
		{
			#pragma omp atomic
			counts2[0]++;
		}

		if(abs(coords[ii][0]-point1[0]) < tolerance)
		{
			#pragma omp atomic
			counts3[0]++;
		}

		/*  y plane */
		if (abs(coords[ii][1]) < tolerance)
		{
			#pragma omp atomic
			counts1[1]++;
		}

		if(abs(abs(coords[ii][1])-1.0) < tolerance)
		{
			#pragma omp atomic
			counts2[1]++;
		}

		if(abs(coords[ii][1]-point1[1]) < tolerance)
		{
			#pragma omp atomic
			counts3[1]++;
		}		

		/* z plane */
		if (abs(coords[ii][2]) < tolerance)
		{
			#pragma omp atomic
			counts1[2]++;
		}

		if(abs(abs(coords[ii][2])-1.0) < tolerance)
		{
			#pragma omp atomic
			counts2[2]++;
		}

		if(abs(coords[ii][2]-point1[2]) < tolerance)
		{
			#pragma omp atomic
			counts3[2]++;
		}

	}

	uint ignore1 = 0;
	uint ignore2 = 0;
	uint ignore3 = 0;

	if (counts1[0] == nPnts/2)
		ignore1 = 1;
	else if (counts1[1] == nPnts/2)
		ignore1 = 2;
	else if (counts1[2] == nPnts/2)
		ignore1 = 3;

	if(counts2[0] == nPnts/2)
		ignore2 = 1;
	else if (counts2[1] == nPnts/2)
		ignore2 = 2;
	else if (counts2[2] == nPnts/2)
		ignore2 = 3;

	if(counts3[0] == nPnts/2)
		ignore3 = 1;
	else if (counts3[1] == nPnts/2)
		ignore3 = 2;
	else if (counts3[2] == nPnts/2)
		ignore3 = 3;


	// cout << counts[0] << "  " << counts[1] << "  " << counts[2] << endl;
	if(ignore1 == 0 && ignore2 == 0 && ignore3 == 0)
	{
		cout << "Couldn't determine which dimension is false." << endl;
		exit(-1);
	}

	/* Default to plane = 0, and if not then use subsequent  */
	uint ignore = 0;

	if(ignore1 != 0)
		ignore = ignore1;
	else if (ignore2 != 0)
		ignore = ignore2;
	else if (ignore3 != 0)
		ignore = ignore3;	

	#ifdef DEBUG
		dbout << "Ignored dimension: " << ignore << endl;
	#endif
	cout << "Ignored dimension: " << ignore << endl;


	vector<StateVecD> coordVec(nPnts/2);
	for (auto const& ii : ptIndex)
	{
		if (ignore == 1)
			coordVec[ii] = StateVecD{coords[ii][1], coords[ii][2]};
		else if (ignore == 2)
			coordVec[ii] = StateVecD{coords[ii][0], coords[ii][2]};
		else if (ignore == 3)
			coordVec[ii] = StateVecD{coords[ii][0], coords[ii][1]};
	}


	#ifdef DEBUG
		dbout << "Returning coordinates." << endl;
	#endif
	return coordVec;	
}

void Get_Data(int& fin, size_t const& nPnts, FACE& fdata, vector<int> markers, 
				int const& symPlane1, int const& symPlane2)
{
	#ifdef DEBUG
	dbout << "Entering Get_Data..." << endl;
	#endif
	int retval;

	/*Faces*/
	vector<vector<size_t>> faceVec;

	int nSTDimID, nPpSTDimID;
	size_t nSTri, nPpST;
	if ((retval = nc_inq_dimid(fin, "no_of_surfacetriangles", &nSTDimID)))
	{
		cout << "No surface triangles in mesh file." << endl;
		#ifdef DEBUG
		dbout << "No surface triangles in mesh file." << endl;
		#endif
	}
	else
	{
		Get_Dimension(fin, "no_of_surfacetriangles", nSTDimID, nSTri);
		Get_Dimension(fin, "points_per_surfacetriangle", nPpSTDimID, nPpST);
		vector<vector<size_t>> localVec = Get_Element(fin, "points_of_surfacetriangles", nSTri, nPpST);

		#ifdef DEBUG
		dbout << "Number of triangles: " << nSTri << endl;
		#endif

		faceVec.insert(faceVec.end(),localVec.begin(),localVec.end());
	}

	int nSQDimID, nPpSQDimID;
	size_t nSQua, nPpSQ;
	if ((retval = nc_inq_dimid(fin, "no_of_surfacequadrilaterals", &nSTDimID)))
	{
		cout << "No surface quadrilaterals in mesh file." << endl;
		#ifdef DEBUG
		dbout << "No surface quadrilaterals in mesh file." << endl;
		#endif
	}
	else
	{
		Get_Dimension(fin, "no_of_surfacequadrilaterals", nSQDimID, nSQua);
		Get_Dimension(fin, "points_per_surfacequadrilateral", nPpSQDimID, nPpSQ);
		vector<vector<size_t>> localVec = Get_Element(fin, "points_of_surfacequadrilaterals", nSQua, nPpSQ);

		#ifdef DEBUG
		dbout << "Number of quadrilaterals: " << nSTri << endl;
		#endif

		faceVec.insert(faceVec.end(), localVec.begin(), localVec.end());
	}

	/*Get the boundarymarkers*/
	int boundMDim;
	size_t nMarkers;
	Get_Dimension(fin, "no_of_surfaceelements", boundMDim, nMarkers);

	if(faceVec.size() != nMarkers)
	{
		cout << "Mismatch of number of surface markers and faces ingested." << endl;
		cout << "Number of surface markers: " << nMarkers << 
		"  Number in vector: " << faceVec.size() << endl;
	}

	int *faceMarkers = new int[nMarkers];

	faceMarkers = Get_Int_Scalar(fin, "boundarymarker_of_surfaces", nMarkers);

	// cout << "symPlane: " << symPlane << endl;

	/*Create the faces of the symmetry plane*/
	vector<vector<size_t>> symFace1, symFace2; 
	for(uint ii = 0; ii < nMarkers; ++ii)
	{
		if(faceMarkers[ii] == symPlane1)
		{	/*The face is part of the symmetry plane*/
			symFace1.emplace_back(faceVec[ii]);
		}
		else if (faceMarkers[ii] == symPlane2)
		{
			symFace2.emplace_back(faceVec[ii]);
		}
		else if(std::find(markers.begin(),markers.end(),faceMarkers[ii])!= markers.end())
		{	/*Its a boundary face*/
			fdata.smarkers.emplace_back(faceMarkers[ii]);
			fdata.sfaces.emplace_back(faceVec[ii]);
		}
	}

	cout << "Symmetry plane 1 faces: " << symFace1.size() << endl;
	cout << "Symmetry plane 2 faces: " << symFace1.size() << endl;
	cout << "        Boundary faces: " << fdata.smarkers.size() << endl;

	#ifdef DEBUG
		dbout << "Symmetry plane 1 faces: " << symFace1.size() << endl;
		dbout << "Symmetry plane 2 faces: " << symFace1.size() << endl;
		dbout << "        Boundary faces: " << fdata.smarkers.size() << endl;
	#endif

	vector<uint> index1;
	for(auto const& face:symFace1)
	{
		index1.insert(index1.end(),face.begin(),face.end());
	}

	vector<uint> index2;
	for (auto const& face : symFace2)
	{
		index2.insert(index2.end(), face.begin(), face.end());
	}

	/* Find the max value in each */
	auto max1 = std::max_element(index1.begin(), index1.end());
	auto max2 = std::max_element(index2.begin(), index2.end());

	cout << "Index 1: " << *max1+1 << "  Index 2: " << *max2+1 << endl;
	vector<uint> ptIndex;
	if (*max1+1 == nPnts/2)
	{
		cout << "Using the first symmetry plane" << endl;
		std::sort(index1.begin(), index1.end());
		index1.erase(std::unique(index1.begin(), index1.end()), index1.end());
		fdata.faces = symFace1;
		ptIndex = index1;
	}
	else // if (*max2+1 == nPnts/2)
	{
		cout << "Using the second symmetry plane" << endl;
		std::sort(index2.begin(), index2.end());
		index2.erase(std::unique(index2.begin(), index2.end()), index2.end());
		fdata.faces = symFace2;
		ptIndex = index2;
	}

	/*Get the coordinates*/
	fdata.verts = Get_Coordinates(fin,nPnts,ptIndex);

	#ifdef DEBUG
	dbout << "Exiting Get_Data..." << endl;
	#endif
}

void Add_Edge(std::pair<uint,uint> const& edge, std::pair<int,int> const& leftright,
	EDGE& edata)
{
	edata.celllr.emplace_back(leftright);
	edata.edges.emplace_back(edge);
	edata.nEdge++;
}


void BuildEdges(const FACE& fdata, EDGE& edata)
{
	#ifdef DEBUG
	dbout << "Entering BuildFaces..." << endl;
	#endif
	uint nFace = fdata.faces.size();
	uint nPts = fdata.nPnts;
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
			const vector<size_t> lcell = fdata.faces[lindex];
	
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
					const vector<size_t> rcell = fdata.faces[rindex];

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

								std::pair<uint,uint> edge(l0,l1);
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

					for (size_t ii = 0; ii < fdata.smarkers.size(); ++ii)
					{	/*Search through the surface faces to identify */
						/*if the face is an boundary face or not*/
						vector<size_t> const& sface = fdata.sfaces[ii];

						int a = (std::find(sface.begin(),sface.end(),l0) != sface.end());
						int b = (std::find(sface.begin(),sface.end(),l1) != sface.end());

						if(a && b)
						{	/*edge is an boundary face*/
							
							std::pair<int,int> leftright(lindex,-1);
							#pragma omp atomic
								colour[lindex][lfaceindex]++;
							std::pair<uint,uint> edge(l0,l1);
							Add_Edge(edge,leftright,elocal);
							elocal.smarkers.emplace_back(fdata.smarkers[ii]);
							elocal.nSurf++;
							goto matchfound;
						}
					}

					/*If still uncoloured, then face is an external boundary*/
					// std::pair<int,int> leftright(lindex,-2);
					// #pragma omp atomic
					// colour[lindex][lfaceindex]++;
					// std::pair<uint,uint> edge(l0,l1);
					// Add_Edge(edge,leftright,elocal);
					// elocal.far.emplace_back(edge);
					// elocal.nSurf++;
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

	if (edata.edges.size() != edata.nEdge)
	{
		cout << "Not all edges have been found." << endl;
		cout << "Edges vector size: "<< edata.edges.size() << "  nEdges: " << edata.nEdge << endl;
	}

	if(edata.celllr.size() != edata.nEdge)
	{
		cout << "Not all edges have left and right cells identified." << endl;
		cout << "Celllr size: " << edata.celllr.size();
		cout << "  nEdges: " << edata.nEdge << endl;
	}

	#ifdef DEBUG
	dbout << "Building edges complete. Number of edges: " << edata.nEdge << endl;
	dbout << "Average number of edges per element: " << float(edata.nEdge)/float(nFace) << endl;
	dbout << "Number of surface edges: " << edata.nSurf << endl;
	dbout << "Exiting BuildFaces..." << endl;
	#endif

	cout << "Building edges complete. Number of edges: " << edata.nEdge << endl;
	cout << "Average number of edges per element: " << float(edata.nEdge)/float(nFace) << endl;
	cout << "Number of surface edges: " << edata.nSurf << endl;
	// cout << edata.wall.size() << "  " << edata.far.size() << endl;
}

int search(vector<uint> array, uint valueToFind)
{
	int pos = 0;
	int length = array.size();
	int limit = std::min(length, 1);

	while (limit < length && array[limit] < valueToFind)
	{
		pos = limit + 1;
		limit = std::min(length, limit * 2 + 1);
	}
	while (pos < limit)
	{
		int testpos = pos + ((limit - pos) >> 1);

		if (array[testpos] < valueToFind)
			pos = testpos + 1;
		else
			limit = testpos;
	}
	return (pos < length && array[pos] == valueToFind ? pos : -1);
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
		vertInUse.emplace_back(edge.first);
		vertInUse.emplace_back(edge.second);
	}

	/*Sort and erase duplicates*/
	std::sort(vertInUse.begin(),vertInUse.end());
	vertInUse.erase(std::unique(vertInUse.begin(),vertInUse.end()),vertInUse.end());

	size_t const newsize = vertInUse.size();
#ifdef DEBUG
	dbout << "vertInUse size: " << newsize << endl;
#endif
	cout << "vertInUse size: " << newsize << endl;
	if(newsize != edata.nPnts/2)
	{
		cout << "Used points do not equal exactly half the mesh points.\n Something has gone wrong." << endl;
		exit(-1); 
	}
	/*Now the edges need recasting to the index of the vector*/
	vector<StateVecD> newVerts(newsize);
	// vector<std::pair<size_t,size_t>> edges;
	#pragma omp parallel shared(newVerts, edata, vertInUse)
	{
		// #pragma omp for schedule(static) nowait
		// for(auto& edge:edata.edges)
		// {
		// 	if(edge.first != vertInUse[edge.first] || edge.second != vertInUse[edge.second])
		// 		cout << "A value that is not in vinuse is used" << endl;
		// } 
	// 		int ifirst = search(vertInUse,edge.first);
	// 		// cout << "Found first index " << ifirst << endl;

	// 		auto lower1 = std::lower_bound(vertInUse.begin(),vertInUse.end(),edge.first);
	// 		if(lower1 != vertInUse.end())
	// 		{
	// 			auto idx = std::distance(vertInUse.begin(),lower1);
	// 			// cout << "Found first index " << idx << endl;
	// 			// edge.first = idx;
	// 			if(edge.first!= ifirst)
	// 			cout << "Point: " << edge.first << "  index1: " << ifirst << "  lower: " << idx << endl;   
	// 		}
	// 		// auto lower2 = std::lower_bound(vertInUse.begin(), vertInUse.end(), edge.second);
	// 		// if (lower2 != vertInUse.end())
	// 		// {
	// 		// 	auto idx = std::distance(vertInUse.begin(), lower2);
	// 		// 	cout << "Found second index " << idx << endl;
	// 		// 	edge.second = idx;
	// 		// }
			


	// 		// int isecond = search(vertInUse,edge.second);
	// 		// // cout << "Found second index " << isecond << endl;
	// 		// edge.first = ifirst;
	// 		// edge.second = isecond;
		// 	if(edge.first > newsize || edge.second > newsize)
		// 		cout << "Vertex used that isn't in the used vertices list" << endl;
		// }

		#pragma omp for 
		for(size_t ii = 0; ii < newsize; ii++)
			newVerts[ii] = edata.verts[vertInUse[ii]];

	}
	edata.usedVerts = vertInUse;
	edata.verts = newVerts;
	edata.nPnts = newsize;

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
	int nElemID, nEdgeID, nPpEcID, nPntsID, nMarkID;

	if ((retval = nc_def_dim(meshID, "no_of_elements", edata.nElem, &nElemID)))
	{
		cout << "Error: Failed to define \"no_of_elements\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_edges", edata.nEdge, &nEdgeID)))
	{
		cout << "Error: Failed to define \"no_of_edges\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "points_per_edge", 2, &nPpEcID)))
	{
		cout << "Error: Failed to define \"points_per_edge\"" << endl;
		ERR(retval);
		exit(-1);
	}

	// if ((retval = nc_def_dim(meshID, "no_of_wall_edges", edata.nWall, &nWallID)))
	// {
	// 	cout << "Error: Failed to define \"no_of_wall_edges\"" << endl;
	// 	ERR(retval);
	// 	exit(-1);
	// }

	// if ((retval = nc_def_dim(meshID, "no_of_farfield_edges", edata.nFar, &nFarID)))
	// {
	// 	cout << "Error: Failed to define \"no_of_farfield_edges\"" << endl;
	// 	ERR(retval);
	// 	exit(-1);
	// }

	if ((retval = nc_def_dim(meshID, "no_of_surfaceelements", edata.nSurf, &nMarkID)))
	{
		cout << "Error: Failed to define \"no_of_surfaceelements\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_points", edata.nPnts, &nPntsID)))
	{
		cout << "Error: Failed to define \"no_of_points\"" << endl;
		ERR(retval);
		exit(-1);
	}

	/* Define the variables */
	int edgeVarID, leftVarID, rightVarID, usedPID, ptsxID, ptszID, markID;

	/*Define the faces*/
	int dimIDs[] = {nEdgeID, nPpEcID};
	if ((retval = nc_def_var(meshID, "points_of_element_edges", NC_INT, 2,
							 dimIDs, &edgeVarID)))
	{
		cout << "Error: Failed to define \"points_of_element_edges\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "left_element_of_edges", NC_INT, 1,
							 &nEdgeID, &leftVarID)))
	{
		cout << "Error: Failed to define \"left_element_of_edges\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "right_element_of_edges", NC_INT, 1,
							 &nEdgeID, &rightVarID)))
	{
		cout << "Error: Failed to define \"right_element_of_edges\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "boundarymarker_of_surfaces", NC_INT, 1,
							 &nMarkID, &markID)))
	{
		cout << "Error: Failed to define \"boundarymarker_of_surfaces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_var(meshID, "vertices_in_use", NC_INT, 1,
							 &nPntsID, &usedPID)))
	{
		cout << "Error: Failed to define \"vertices_in_use\"" << endl;
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
	int* edges = new int[edata.nEdge*2];
	for(uint ii = 0; ii < edata.nEdge; ++ii)
	{
		edges[index(ii,0,2)] = static_cast<int>(edata.edges[ii].first);
		edges[index(ii,1,2)] = static_cast<int>(edata.edges[ii].second);
	}

	/*Put edges into the file*/
	size_t start[] = {0, 0};
	size_t end[] = {edata.nEdge, 2};
	if ((retval = nc_put_vara_int(meshID, edgeVarID, start, end, &edges[0])))
	{
		cout << "Failed to write edge data" << endl;
		ERR(retval);
		exit(-1);
	}

	/*Put face left and right into the file*/
	int* left = new int[edata.nEdge];
	int* right = new int[edata.nEdge];

	for(uint ii = 0; ii < edata.nEdge; ++ii)
	{
		left[ii] = edata.celllr[ii].first;
		right[ii] = edata.celllr[ii].second;
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

	

	if ((retval = nc_put_var_int(meshID, markID, &edata.smarkers[0])))
	{
		cout << "Failed to write surface marker data" << endl;
		ERR(retval);
		exit(-1);
	}

	/*Write which vertices are in use, for reading the solution file*/
	int* uVert = new int[edata.nPnts];

	for(size_t ii = 0; ii < edata.nPnts; ++ii)
	{
		uVert[ii] = edata.usedVerts[ii]; 
	}

	if ((retval = nc_put_var_int(meshID, usedPID, &uVert[0])))
	{
		cout << "Failed to write used vertices data" << endl;
		ERR(retval);
		exit(-1);
	}

	/*Create the C arrays for the vertices*/
	double* x = new double[edata.nPnts];
	double* z = new double[edata.nPnts];

	for(uint ii = 0; ii < edata.nPnts; ++ii)
	{
		x[ii] = edata.verts[ii][0];
		z[ii] = edata.verts[ii][1];
	}

	/*Put them in the file*/
	if ((retval = nc_put_var_double(meshID, ptsxID, &x[0])))
	{
		cout << "Failed to write x coordinate." << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_put_var_double(meshID, ptszID, &z[0])))
	{
		cout << "Failed to write z coordinate." << endl;
		ERR(retval);
		exit(-1);
	}
}


void Write_Real_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
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
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_Int_Vector(void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
						vector<int32_t> const& varVec, string const& varName)
{
	if(tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &varVec[0]))
	{
		cout << "Failed to write \"" << varName << "\". frame: " << 
			outputZone << ". varCount: " << varCount << endl;
		exit(-1);
	}
	varCount++;
}

void Write_Edge_Tecplot_Binary(const EDGE& edata)
{
	int32_t debug = 0;
	#ifdef DEBUG
	debug = 0;
	#endif
	string file = "griduns.plt";
	string varNames = "X,Z";
	string title = "2D Edge based grid conversion output";
	string group = "Mesh data";
	int32_t fileType = 0;
	int32_t fileFormat = 0;
	int32_t fileIsDouble = 1;

	if(TECINI142(title.c_str(),varNames.c_str(),file.c_str(),".",
				&fileFormat,&fileType,&debug,&fileIsDouble))
	{
		cout << "Failed to open " << file << endl;
		exit(-1);
	}

	int32_t zoneType = 6; /* FE Polygon */	
	int32_t numNodes = edata.nPnts;
	int32_t numElems = edata.nElem;
	int64_t numFaces = edata.nEdge;

	/* Need to find how many face nodes there are */
	int64_t numFaceNodes = 2*edata.edges.size();
    double  solTime = 0.0;
    int32_t strandID = 0;     // Static Zone
    int32_t unused = 0;       // ParentZone is no longer used
    int32_t numBoundaryFaces = 0;
    int32_t numBoundaryConnections = 0;
    int32_t shareConnectivity = 0;

    if(TECPOLYZNE142(
        "Polygonal Quad Zone",
        &zoneType,
        &numNodes,
        &numElems,
        &numFaces,
        &numFaceNodes,
        &solTime,
        &strandID,
        &unused,
        &numBoundaryFaces,
        &numBoundaryConnections,
        NULL,
        NULL,  // All nodal variables
        NULL,
        &shareConnectivity))
    {
        printf("Polyquads: error calling TECPOLYZNE\n");
        exit(-1);
    }

	real* x  = new real[edata.nPnts];
	int32_t nPnts = edata.nPnts;
	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(size_t ii = 0; ii < edata.nPnts; ++ii)
			x[ii] = edata.verts[ii][dim];

		string name = "position coordinate " + std::to_string(dim);

		TECDAT142(&nPnts, x, &fileIsDouble);
	}

	delete[] x;

	int32_t* faceNodes = new int32_t[2*numFaces];
	int32_t* left = new int32_t[numFaces];
	int32_t* right = new int32_t[numFaces];

	/*Inform of how many vertices in each face*/
	for (size_t ii = 0; ii < edata.nEdge; ++ii)
	{
		faceNodes[ii*2] = edata.edges[ii].first+1; /*Write face vertex indexes*/
		faceNodes[ii*2+1] = edata.edges[ii].second+1;
		left[ii] = edata.celllr[ii].first+1;
		right[ii] = edata.celllr[ii].second > -1 ? edata.celllr[ii].second+1 : 0;
	}
	int32_t faceOffset = edata.nEdge;
	if(TECPOLYFACE142(
                &faceOffset,
                NULL,
                &faceNodes[0],
                &left[0],
                &right[0]))
	{
		printf("Error calling TECPOLYFACE\n");
        exit(-1);
	}

	delete [] faceNodes;
	delete [] left;
	delete [] right;

	if(TECEND142())
	{
		printf("Polyquads: error calling TECEND\n");
        exit(-1);
	}

}

void Write_Edge_szplt(const EDGE& edata)
{
	#if FOD == 1
		int32_t realType = 2;
	#else
		int32_t realType = 1;
	#endif

	void* fileHandle;
	int32_t outputZone;

	string file = "griduns.plt";
	vector<int32_t> varTypes = {realType,realType};
	string varNames = "X,Z";
	string title = "2D Edge based grid conversion output";
	string group = "Mesh data";
	int32_t fileType = 0;
	int32_t fileFormat = FILEFORMAT_SZL;

	if(tecFileWriterOpen(file.c_str(),title.c_str(),varNames.c_str(),fileType,fileFormat,1,NULL,&fileHandle))
	{
		cout << "Failed to open " << file << endl;
		exit(-1);
	}
	#ifdef DEBUG
	if(tecFileSetDiagnosticsLevel(fileHandle, 1))
	{
		std::cerr << "Failed to set debug option for output file: " << file << endl;
		exit(-1);
	}
	#endif

	vector<int32_t> shareVarFromZone(varTypes.size(),0);
	vector<int32_t> valueLocation(varTypes.size(),1);
	vector<int32_t> passiveVarList(varTypes.size(),0);
 
	int32_t zoneType = 6; /* FE Polygon */	
	int64_t totalNumFaceNodes = 2*edata.edges.size();
	if(tecZoneCreatePoly(fileHandle,group.c_str(),zoneType,edata.nPnts,edata.nEdge,edata.nElem,
		totalNumFaceNodes,&varTypes[0],&shareVarFromZone[0],&valueLocation[0],
		&passiveVarList[0],0,0,0,&outputZone))
	{
		std::cerr << "Failed to create polygonal zone." << endl;
		exit(-1);
	}

	vector<real> x(edata.nPnts);
	int32_t var = 1;

	for(uint dim = 0; dim < SIMDIM; ++dim)
	{
		#pragma omp parallel for
		for(size_t ii = 0; ii < edata.nPnts; ++ii)
			x[ii] = edata.verts[ii][dim];

		string name = "position coordinate " + std::to_string(dim);
		Write_Real_Vector(fileHandle, outputZone, var, edata.nPnts, x, name);
	}

	vector<int32_t> faceCounts(edata.edges.size());
	vector<int32_t> left(edata.edges.size());
	vector<int32_t> right(edata.edges.size());
	vector<int32_t> faceNodes(totalNumFaceNodes);

	/*Inform of how many vertices in each face*/
	for (size_t ii = 0; ii < edata.edges.size(); ++ii)
	{
		faceCounts[ii] =  2;
		faceNodes[ii*2] = edata.edges[ii].first+1; /*Write face vertex indexes*/
		faceNodes[ii*2+1] = edata.edges[ii].second+1;
		left[ii] = edata.celllr[ii].first+1;
		right[ii] = edata.celllr[ii].second > -1 ? edata.celllr[ii].second+1 : 0;
	}

	tecZoneWritePolyFaces32(fileHandle,outputZone,0,edata.nEdge,&faceCounts[0],
					&faceNodes[0],&left[0],&right[0],1);

	tecFileWriterClose(&fileHandle);
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
	fout << "NODES=" << edata.nPnts << " ELEMENTS=" << edata.nElem << " FACES=" << edata.nEdge << endl;
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
			fout << std::setw(w) << edata.verts[ii][DIM];
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
			/*Write face vertex indexes*/
			fout << std::setw(w) << edata.edges[ii].first+1;
			if (edata.edges[ii].first > edata.nPnts)
			{
				cout << "Trying to write a vertex outside of the number of points." << endl;
			}
			fout << std::setw(w) << edata.edges[ii].second+1;
			if (edata.edges[ii].first > edata.nPnts)
			{
				cout << "Trying to write a vertex outside of the number of points." << endl;
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

void Write_griduns(const EDGE& edata)
{
	/*File format:*/
	/*<nCell> <nEdge> <nVert>*/
	/*<edge1Vert1> <edge1Vert2> <edge1LeftCell> <edge1RightCell> */
	/*<edge2Vert1> <edge2Vert2> <edge2LeftCell> <edge2RightCell> */
	/*...*/
	/*<vert1X> <vert1Y>*/
	/*<vert2X> <vert2Y>*/
	/*...*/

	#ifdef DEBUG
	dbout << "Entering Write_griduns..." << endl;
	cout << "Attempting write output file." << endl;
	cout << "File: " << "griduns" << endl;
	#endif

	/*Establish edge orientation*/
	/*Find an edge where cell left is common between two neighboring edges*/
	int leftright = 0;
	size_t jj = 0;
	while(edata.celllr[jj].first != edata.celllr[jj+1].first)
	{
		jj++;
	}

	/*Share the same cell, so have three vertices to use to test*/
	/*make sure the three points are unique*/
	if(edata.edges[jj].first != edata.edges[jj+1].first)
	{
		if(edata.edges[jj].second != edata.edges[jj+1].first)
		{
			/*Have 3 unique points*/
			leftright = Distance(edata.verts,edata.edges[jj],edata.edges[jj+1].first);
		}
		
	}
	else
	{
		if(edata.edges[jj].second != edata.edges[jj+1].second)
		{
			leftright = Distance(edata.verts,edata.edges[jj],edata.edges[jj+1].second);
		}
	}

	/*If leftright = 0, then first cell is on the right*/


	std::ofstream fout("griduns",std::ios::out);
	if(!fout.is_open())
	{
		cout << "Couldn't open the output file." << endl;
		exit(-1);
	}


	fout << std::left << std::fixed;
	uint w = 9;

	fout << std::setw(w) << edata.nElem << std::setw(w) << edata.nEdge
		 << std::setw(w) << edata.nPnts << endl;

	if(leftright == 1)
	{

		for(size_t ii = 0; ii < edata.edges.size(); ++ii)
		{	/*Add 1 (for fortran bull) if not a boundary edge*/
			if(edata.celllr[ii].second == -1)
			{
				fout << std::setw(w) << edata.edges[ii].second+1 << std::setw(w) << edata.edges[ii].first+1;
				fout << std::setw(w) << edata.celllr[ii].second << std::setw(w) << edata.celllr[ii].first+1;
			}
			else if (edata.celllr[ii].second == -2)
			{
				fout << std::setw(w) << edata.edges[ii].first+1 << std::setw(w) << edata.edges[ii].second+1;
				fout << std::setw(w) << edata.celllr[ii].first+1 << std::setw(w) << edata.celllr[ii].second;
			}
			else
			{
				fout << std::setw(w) << edata.edges[ii].first+1 << std::setw(w) << edata.edges[ii].second+1;
				fout << std::setw(w) << edata.celllr[ii].first+1 << std::setw(w) << edata.celllr[ii].second+1;
			}
		
			 
			fout << endl;
		}
	}
	else
	{
		for(size_t ii = 0; ii < edata.edges.size(); ++ii)
		{	/*Add 1 (for fortran bull) if not a boundary edge*/
			if(edata.celllr[ii].second == -1)
			{
				fout << std::setw(w) << edata.edges[ii].first+1 << std::setw(w) << edata.edges[ii].second+1;
				fout << std::setw(w) << edata.celllr[ii].second << std::setw(w) << edata.celllr[ii].first+1;
			}
			else if (edata.celllr[ii].second == -2)
			{
				fout << std::setw(w) << edata.edges[ii].second+1 << std::setw(w) << edata.edges[ii].first+1;
				fout << std::setw(w) << edata.celllr[ii].first+1 << std::setw(w) << edata.celllr[ii].second;
			}
			else
			{
				fout << std::setw(w) << edata.edges[ii].first+1 << std::setw(w) << edata.edges[ii].second+1;
				fout << std::setw(w) << edata.celllr[ii].second+1 << std::setw(w) << edata.celllr[ii].first+1;
			}
		
			 
			fout << endl;
		}
	}

	w = 17;
	uint preci = 8;
	fout << std::left << std::scientific << std::setprecision(preci);

	for(size_t ii = 0 ; ii < edata.verts.size(); ++ii)
	{
		fout << std::setw(8) << ii+1 << std::setw(w) << edata.verts[ii][0] << std::setw(w) << edata.verts[ii][1] << endl;
	}

	fout.close();

	#ifdef DEBUG
	dbout << "Exiting Write_griduns..." << endl;
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

	int meshID;
	int retval;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{
		cout << "Failed to open mesh file \"" << meshIn << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading cell data..." << endl;

	int ptDimID, elemDimID, faceDimID;
	size_t nPnts, nElem, nFace;

	// Retrieve how many elements there are.
	Get_Dimension(meshID, "no_of_elements", elemDimID, nElem);
	Get_Dimension(meshID, "no_of_points", ptDimID, nPnts);
	Get_Dimension(meshID, "no_of_surfaceelements", faceDimID, nFace);


	FACE fdata(nElem,nPnts);

	cout << "nElem: " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
	
	/*Get surface faces*/
	int symPlane1 = 0;
	int symPlane2 = 0;
	vector<int> markers = Find_Bmap_Markers(bmapIn,symPlane1, symPlane2);
	Get_Data(meshID,nPnts,fdata,markers,symPlane1,symPlane2);

	cout << "Face size: " << fdata.faces.size() << endl;

	EDGE edata(fdata);

	cout << "Building edge based data..." << endl;
	BuildEdges(fdata,edata);

	Recast_Data(edata);

	cout << "Writing data..." << endl;
	Write_Edge_Data(meshIn,edata);

	// Write_Edge_Tecplot(edata);
	Write_Edge_Tecplot_Binary(edata);
	// Write_griduns(edata);

	return 0;
}