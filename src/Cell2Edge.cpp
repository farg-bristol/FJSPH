#include "Convert.h"

#ifdef DEBUG
	/*Open debug file to write to*/
	std::ofstream dbout("Cell2Edge.log",std::ios::out);
#endif

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 2
#endif

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real,SIMDIM,1> StateVecD;
typedef Eigen::Matrix<int,SIMDIM,1> StateVecI;

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

	size_t numElem, numEdges, nFar, nWall;
	vector<StateVecD> verts;
	vector<std::pair<size_t,size_t>> edges; /*edge indexes*/
	vector<std::pair<int,int>> celllr; /*Cell left and right of the face*/
	vector<std::pair<size_t,size_t>> wall, far;
	vector<uint> usedVerts;
	// int* nFacesPElem;

	uint numPoint;
}EDGE;

int Distance(vector<StateVecD> const &verts, std::pair<size_t, size_t> const &edge, size_t const &point)
{
    real  ty, tx;
    StateVecD vtx0, vtx1;

    tx = verts[point](0);
    ty = verts[point](1);

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
		else if(line.find("Type: laminar wall")!=string::npos ||
			line.find("Type: euler wall")!=string::npos || 
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
vector<StateVecD> Get_Coordinates(int& fin, size_t const& nPnts, vector<uint> const& ptIndex)
{
	#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
	#endif

	vector<Eigen::Vector3d> coords = Get_Coordinate_Vector(fin, nPnts);
	/*Need to find which coordinate to ignore*/
	#ifdef DEBUG
		dbout << "Checking which dimension to ignore" << endl;
	#endif
	double tolerance = 1e-6;

	uint counts[3] = {0};

	#pragma omp parallel for
	for(auto const ii : ptIndex)
	{
		if (abs(coords[ii](0)) < tolerance || (abs(coords[ii](0)) < 1 + tolerance && abs(coords[ii](0)) > 1 - tolerance))
		{
			#pragma omp atomic
			counts[0]++;
		}

		if (abs(coords[ii](1)) < tolerance || (abs(coords[ii](1)) < 1 + tolerance && abs(coords[ii](1)) > 1 - tolerance))
		{
			#pragma omp atomic
			counts[1]++;
		}

		if (abs(coords[ii](2)) < tolerance || (abs(coords[ii](2)) < 1 + tolerance && abs(coords[ii](2)) > 1 - tolerance))
		{
			#pragma omp atomic
			counts[2]++;
		}

	}

	uint ignore=0;

	if (counts[0] == nPnts/2)
		ignore = 1;
	else if (counts[1] == nPnts/2)
		ignore = 2;
	else if (counts[2] == nPnts/2)
		ignore = 3;

	// cout << counts[0] << "  " << counts[1] << "  " << counts[2] << endl;
	if(ignore == 0)
	{
		cout << "Couldn't determine which dimension is false." << endl;
		exit(-1);
	}

	#ifdef DEBUG
		dbout << "Ignored dimension: " << ignore << endl;
	#endif
	cout << "Ignored dimension: " << ignore << endl;


	vector<StateVecD> coordVec(nPnts/2);
	for (auto const& ii : ptIndex)
	{
		if (ignore == 1)
			coordVec[ii] = StateVecD(coords[ii](1), coords[ii](2));
		else if (ignore == 2)
			coordVec[ii] = StateVecD(coords[ii](0), coords[ii](2));
		else if (ignore == 3)
			coordVec[ii] = StateVecD(coords[ii](0), coords[ii](1));
	}


	#ifdef DEBUG
		dbout << "Returning coordinates." << endl;
	#endif
	return coordVec;	
}

void Get_Data(int& fin, size_t const& nPnts, FACE& fdata, vector<int> markers, int const& symPlane1, int const& symPlane2)
{
	#ifdef DEBUG
	dbout << "Entering Get_Data..." << endl;
	#endif
	int retval;

	/*Faces*/
	vector<vector<uint>> faceVec;

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
		vector<vector<uint>> localVec = Get_Element(fin, "points_of_surfacetriangles", nSTri, nPpST);

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
		vector<vector<uint>> localVec = Get_Element(fin, "points_of_surfacequadrilaterals", nSQua, nPpSQ);

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
	vector<vector<uint>> symFace1, symFace2; 
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
		{	/*Its an internal face*/
			fdata.internal.emplace_back(faceVec[ii]);
		}
	}

	// cout << "Symmetry plane 1 faces: " << symFace1.size() << endl;
	// cout << "Symmetry plane 2 faces: " << symFace1.size() << endl;
	// cout << "   Internal faces: " << fdata.internal.size() << endl;

	#ifdef DEBUG
		dbout << "Symmetry plane 1 faces: " << symFace1.size() << endl;
		dbout << "Symmetry plane 2 faces: " << symFace1.size() << endl;
		dbout << "   Internal faces: " << fdata.internal.size() << endl;
	#endif

	vector<uint> index1;
	for(auto const face:symFace1)
	{
		index1.insert(index1.end(),face.begin(),face.end());
	}

	vector<uint> index2;
	for (auto const face : symFace2)
	{
		index2.insert(index2.end(), face.begin(), face.end());
	}

	/* Find the max value in each */
	auto max1 = std::max_element(index1.begin(), index1.end());
	auto max2 = std::max_element(index2.begin(), index2.end());

	cout << *max1 << "  " << *max2 << endl;
	vector<uint> ptIndex;
	if (*max1+1 == nPnts/2)
	{
		cout << "Using the first symmetry plane" << endl;
		std::sort(index1.begin(), index1.end());
		index1.erase(std::unique(index1.begin(), index1.end()), index1.end());
		fdata.faces = symFace1;
		ptIndex = index1;
	}

	if (*max2+1 == nPnts/2)
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
							std::pair<uint,uint> edge(l0,l1);
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
					std::pair<uint,uint> edge(l0,l1);
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
	/*Now the edges need recasting to the index of the vector*/
	vector<StateVecD> newVerts(newsize);
	// vector<std::pair<size_t,size_t>> edges;
	#pragma omp parallel shared(newVerts, edata, vertInUse)
	{
		// #pragma omp for schedule(static) nowait
		for(auto& edge:edata.edges)
		{
			if(edge.first != vertInUse[edge.first] || edge.second != vertInUse[edge.second])
				cout << "A value that is not in vinuse is used" << endl;
		} 
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
	edata.numPoint = newsize;

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
	int nElemID, nEdgeID, nPpEcID, nWallID, nFarID, nPntsID;

	if ((retval = nc_def_dim(meshID, "no_of_elements", edata.numElem, &nElemID)))
	{
		cout << "Error: Failed to define \"no_of_elements\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_edges", edata.numEdges, &nEdgeID)))
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

	if ((retval = nc_def_dim(meshID, "no_of_wall_edges", edata.nWall, &nWallID)))
	{
		cout << "Error: Failed to define \"no_of_wall_edges\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_farfield_edges", edata.nFar, &nFarID)))
	{
		cout << "Error: Failed to define \"no_of_farfield_edges\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_points", edata.numPoint, &nPntsID)))
	{
		cout << "Error: Failed to define \"no_of_points\"" << endl;
		ERR(retval);
		exit(-1);
	}

	/* Define the variables */
	int edgeVarID, leftVarID, rightVarID, usedPID, ptsxID, ptszID;

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
	int* edges = new int[edata.numEdges*2];
	for(uint ii = 0; ii < edata.numEdges; ++ii)
	{
		edges[index(ii,0,2)] = static_cast<int>(edata.edges[ii].first);
		edges[index(ii,1,2)] = static_cast<int>(edata.edges[ii].second);
	}

	/*Put edges into the file*/
	size_t start[] = {0, 0};
	size_t end[] = {edata.numEdges, 2};
	if ((retval = nc_put_vara_int(meshID, edgeVarID, start, end, &edges[0])))
	{
		cout << "Failed to write edge data" << endl;
		ERR(retval);
		exit(-1);
	}

	/*Put face left and right into the file*/
	int* left = new int[edata.numEdges];
	int* right = new int[edata.numEdges];

	for(uint ii = 0; ii < edata.numEdges; ++ii)
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

	/*Write which vertices are in use, for reading the solution file*/
	int* uVert = new int[edata.numPoint];

	for(size_t ii = 0; ii < edata.numPoint; ++ii)
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
	double* x = new double[edata.numPoint];
	double* z = new double[edata.numPoint];

	for(uint ii = 0; ii < edata.numPoint; ++ii)
	{
		x[ii] = edata.verts[ii](0);
		z[ii] = edata.verts[ii](1);
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
			/*Write face vertex indexes*/
			fout << std::setw(w) << edata.edges[ii].first+1;
			if (edata.edges[ii].first > edata.numPoint)
			{
				cout << "Trying to write a vertex outside of the number of points." << endl;
			}
			fout << std::setw(w) << edata.edges[ii].second+1;
			if (edata.edges[ii].first > edata.numPoint)
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

	fout << std::setw(w) << edata.numElem << std::setw(w) << edata.numEdges << std::setw(w) << edata.numPoint << endl;

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
		fout << std::setw(8) << ii+1 << std::setw(w) << edata.verts[ii](0) << std::setw(w) << edata.verts[ii](1) << endl;
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
	// Write_griduns(edata);
	return 0;
}