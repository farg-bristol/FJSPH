/*Cell based to Face based data converter*/

#include <algorithm>
#include <omp.h>
#include "Convert.h"

// #include "Third_Party/Eigen/Core"
// #include "Third_Party/Eigen/StdVector"

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

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real,3,1> StateVecD;
typedef Eigen::Matrix<int,3,1> StateVecI;

const std::string WHITESPACE = " \n\r\t\f\v";

typedef struct CELL
{

	/*Standard contructor*/
	CELL(const size_t nElem_, const size_t nPts_)
	{
		nElem = nElem_; 
		nPnts = nPts_;
		verts.reserve(nPts_);
		elems.reserve(nElem_);
	}
	
	/*Zone info*/
	size_t nElem, nPnts;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Cell based data*/
	vector<vector<size_t>> elems;

	/*Surface faces*/
	size_t nSurf, nSurfQ, nSurfT;
	vector<vector<size_t>> stfaces, sqfaces;
	vector<int> smarkers;
}CELL;

typedef struct FACE
{
	// EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	FACE(const CELL &cdata) : nElem(cdata.nElem), nPnts(cdata.nPnts)
	{
		verts = cdata.verts;
		// numElem = cdata.numElem;
		// numPoint = cdata.numPoint;
		nFaces = 0;
		ntFaces = 0;
		nqFaces = 0;
		nSurf = 0;
		// nFacesPElem = new int[cdata.numElem];
	}	

	FACE() : nElem(0), nPnts(0)
	{
		nSurf = 0;
		nFaces = 0;
		ntFaces = 0;
		nqFaces = 0;
	};

	void insert(const FACE& flocal)
	{
		trig_faces.insert(trig_faces.end(),flocal.trig_faces.begin(),flocal.trig_faces.end());
		quad_faces.insert(quad_faces.end(),flocal.quad_faces.begin(),flocal.quad_faces.end());
		qcelllr.insert(qcelllr.end(),flocal.qcelllr.begin(),flocal.qcelllr.end());
		tcelllr.insert(tcelllr.end(),flocal.tcelllr.begin(),flocal.tcelllr.end());
		
		surf_trigs.insert(surf_trigs.end(),flocal.surf_trigs.begin(),flocal.surf_trigs.end());
		surf_quads.insert(surf_quads.end(),flocal.surf_quads.begin(),flocal.surf_quads.end());
		strigslr.insert(strigslr.end(),flocal.strigslr.begin(),flocal.strigslr.end());
		squadslr.insert(squadslr.end(),flocal.squadslr.begin(),flocal.squadslr.end());

		tsmarkers.insert(tsmarkers.end(),flocal.tsmarkers.begin(),flocal.tsmarkers.end());
		qsmarkers.insert(qsmarkers.end(),flocal.qsmarkers.begin(),flocal.qsmarkers.end());

		nFaces += flocal.nFaces;
		ntFaces += flocal.ntFaces;
		nqFaces += flocal.nqFaces;
		nSurf += flocal.nSurf;
	}

	size_t nFaces, ntFaces, nqFaces;
	size_t nSurf;
	vector<StateVecD> verts;
	vector<vector<size_t>> trig_faces; /*Face indexes*/
	vector<vector<size_t>> quad_faces; /*Face indexes*/
	vector<std::pair<int,int>> qcelllr, tcelllr; /*Cell left and right of the face*/


	/* Surface data */
	vector<vector<size_t>> surf_trigs, surf_quads;
	vector<std::pair<int,int>> strigslr, squadslr;
	vector<int> tsmarkers, qsmarkers;


	const size_t nElem, nPnts;
}FACE;

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string Get_Parameter_Value(string const& line)
{
    size_t pos = line.find(":");
    size_t end = line.find("#",pos+1); /* Check if a comment exists on the line */

    if (end != string::npos)
    {
        string value = line.substr(pos + 1, (end-pos+2) );
        return ltrim(rtrim(value));
    }

    string value = line.substr(pos + 1);
    return ltrim(rtrim(value));
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

	/* Find out how many blocks, i.e. how many boundaries to expect. */
	uint nBlocks = 0;
	string line;
	while (getline(fin, line))
	{
		if (line.find("block end") != string::npos)
		{ /*We're on a new block*/
			nBlocks++;
		}
	}
	fin.seekg(0,fin.beg);

	vector<int> markers(nBlocks);

	size_t blockno = 0;

	while(getline(fin,line))
	{

		if (line.find("Markers") != string::npos)
		{
			string value = Get_Parameter_Value(line);
			std::istringstream sstr(value);
			sstr >> markers[blockno];
		}

		if(line.find("block end")!= string::npos)
		{	/*We're on a new block*/
			blockno++;
		}
	}

	#ifdef DEBUG
	dbout << "Exiting Find_Bmap_Markers..." << endl;
	#endif

	return markers;
}


void Get_Surface(int& fin, CELL& cdata)
{
	#ifdef DEBUG
	dbout << "Entering Get_Surface..." << endl;
	#endif
	int retval;

	int surfTDim, surfQDim, pPSTDim, pPSQDim;
	size_t nsurfT, nsurfQ, nPpST, nPpSQ;

	/* Get surface dimensions */
	if ((retval = nc_inq_dimid(fin, "no_of_surfacetriangles", &surfTDim)))
	{
		cout << "No surfacetriangle data" << endl;
		cdata.nSurfT = 0;
	}
	else
	{
		Get_Dimension(fin, "no_of_surfacetriangles", surfTDim, nsurfT);
		Get_Dimension(fin, "points_per_surfacetriangle", pPSTDim, nPpST);

		#ifdef DEBUG
			dbout << "Number of triangles: " << nsurfT << endl;
		#endif

		vector<vector<size_t>> localVec = Get_Element(fin, "points_of_surfacetriangles", nsurfT, nPpST);

		cdata.stfaces.insert(cdata.stfaces.end(), localVec.begin(), localVec.end());
		cdata.nSurfT = nsurfT;
	}


	if ((retval = nc_inq_dimid(fin, "no_of_surfacequadrilaterals", &surfTDim)))
	{
		cout << "No surfacequadrilateral data" << endl;
		cdata.nSurfQ = 0;
	}
	else
	{
		Get_Dimension(fin, "no_of_surfacequadrilaterals", surfQDim, nsurfQ);
		Get_Dimension(fin, "points_per_surfacequadrilateral", pPSQDim, nPpSQ);
		
		#ifdef DEBUG
			dbout << "Number of quadrilaterals: " << nsurfQ << endl;
		#endif

		vector<vector<size_t>> localVec = Get_Element(fin, "points_of_surfacequadrilaterals", nsurfQ, nPpSQ);

		cdata.sqfaces.insert(cdata.sqfaces.end(),localVec.begin(),localVec.end());
		cdata.nSurfQ = nsurfQ;
	}

	/*Get the boundarymarkers*/
	int boundMDim;
	size_t nMarkers;
	Get_Dimension(fin, "no_of_surfaceelements", boundMDim, nMarkers);

	cdata.nSurf = cdata.nSurfT + cdata.nSurfQ;

	if(cdata.nSurf != nMarkers)
	{
		cout << "Mismatch of number of surface elements defined and number ingested." << endl;
		cout << "Number of surface elements: " << nMarkers << 
		"  Number in vectors: " << cdata.nSurf << endl;
	}

	int* faceMarkers = new int[nMarkers];

	faceMarkers = Get_Int_Scalar(fin, "boundarymarker_of_surfaces", nMarkers);

	cdata.smarkers = vector<int>(nMarkers);
	for(size_t ii = 0; ii < nMarkers; ++ii)
	{
		cdata.smarkers[ii] = faceMarkers[ii];
	}
	
	// vector<vector<size_t>> farVec;

	// /*Want to find which surfaces are the ones I want to look for*/
	// for(size_t ii = 0; ii < nMarkers; ++ii)
	// {
	// 	if(std::find(markers.begin(),markers.end(),faceMarkers[ii])!= markers.end())
	// 	{	/*The face is an inner boundary*/
	// 		/*Pre sort to save time in the loop*/
	// 		vector<size_t> v = faceVec[ii];
	// 		std::sort(v.begin(),v.end());
	// 		farVec.emplace_back(v);
	// 	}
	// }

	#ifdef DEBUG
	dbout << "Exiting Get_Surface..." << endl;
	#endif

	// return farVec;
}

void Add_Face(vector<size_t> const& face, std::pair<int,int> const& leftright,
	FACE& fdata)
{
	/*Break the face down into two triangles*/
	if(face.size() == 4)
	{	/*Break the face down into two triangles*/
		// 	fdata.celllr.emplace_back(leftright);
		// 	fdata.celllr.emplace_back(leftright);
		// 	vector<size_t> face1 = {face[0],face[1],face[2]};
		// 	vector<size_t> face2 = {face[0],face[2],face[3]};
		// 	fdata.faces.emplace_back(face1);
		// 	fdata.faces.emplace_back(face2);
		// 	fdata.numFaces+=2;
	
		fdata.qcelllr.emplace_back(leftright);
		fdata.quad_faces.emplace_back(face);
		fdata.nqFaces++;
	}
	else
	{	/*Face is already a triangle. No work to be done*/
		fdata.tcelllr.emplace_back(leftright);
		fdata.trig_faces.emplace_back(face);
		fdata.ntFaces++;
	}

	fdata.nFaces++;
}

void Add_Surface_Face(vector<size_t> const& face, std::pair<int,int> const& leftright,
	FACE& fdata)
{
	/*Break the face down into two triangles*/
	if(face.size() == 4)
	{	/*Break the face down into two triangles*/
		// 	fdata.celllr.emplace_back(leftright);
		// 	fdata.celllr.emplace_back(leftright);
		// 	vector<size_t> face1 = {face[0],face[1],face[2]};
		// 	vector<size_t> face2 = {face[0],face[2],face[3]};
		// 	fdata.faces.emplace_back(face1);
		// 	fdata.faces.emplace_back(face2);
		// 	fdata.numFaces+=2;
	
		fdata.squadslr.emplace_back(leftright);
		fdata.surf_quads.emplace_back(face);
		fdata.nqFaces++;
	}
	else
	{	/*Face is already a triangle. No work to be done*/
		fdata.strigslr.emplace_back(leftright);
		fdata.surf_trigs.emplace_back(face);
		fdata.ntFaces++;
	}

	fdata.nFaces++;
}


void CheckFaces(const vector<vector<size_t>>& vertincells,
	 const vector<vector<size_t>>& lfaces, const size_t lindex, const size_t lgeom, 
	const CELL& cdata, vector<vector<size_t>>& colour, FACE& fdata)
{
	std::pair<int,int> leftright;
	size_t lfaceindex = 0;
	for(auto const& face:lfaces)
	{	/*Define that the cell of the top-level cell search is on the 'left'*/

		/*Create an ordered list of the left side of the face*/
		vector<size_t> lface = face;
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
				vector<size_t> rcell = cdata.elems[rindex];
				vector<vector<size_t>> rfaces;

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

				size_t rfaceindex = 0;
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
			if(lface.size() == 3)
			{		
				for (size_t ii = 0; ii < cdata.nSurfT; ++ii)
				{	/*Search through the surface faces to identify */
					/*if the face is an internal face or not*/
					vector<size_t> sface = cdata.stfaces[ii];
					std::sort(sface.begin(),sface.end());
					if(sface == lface)
					{	/*Face is an internal face*/
						leftright = std::pair<int,int>(lindex,-1);
						fdata.tsmarkers.emplace_back(cdata.smarkers[ii]);
						colour[lindex][lfaceindex] = 1;
						Add_Surface_Face(face,leftright,fdata);
						fdata.nSurf++;
						goto matchfound;
					}
				}
			}
			else
			{		
				for (size_t ii = 0; ii < cdata.nSurfQ; ++ii)
				{	/*Search through the surface faces to identify */
					/*if the face is an internal face or not*/
					vector<size_t> sface = cdata.sqfaces[ii];
					std::sort(sface.begin(),sface.end());
					if(sface == lface)
					{	
						leftright = std::pair<int,int>(lindex,-1);
						fdata.qsmarkers.emplace_back(cdata.smarkers[ii+cdata.nSurfT]);
						colour[lindex][lfaceindex] = 1;
						Add_Surface_Face(face,leftright,fdata);
						fdata.nSurf++;
						goto matchfound;
					}
				}
			}
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
	size_t nElem = cdata.nElem;
	size_t nPts = cdata.nPnts;
	size_t unmfaces=0;

	// for(auto const& vert:cdata.elems[33])
	// 	cout << vert << "  ";

	// cout << endl;

	/*Vector of which cells the vertex is referenced in*/
	vector<vector<size_t>> vertincells(nPts);

	for(size_t ii = 0; ii < nElem; ++ii)
	{
		for(auto const& vert:cdata.elems[ii])
		{	/*Add count for every mention of the vertex, and push back the cell id*/
			vertincells[vert].emplace_back(ii);
		}
	}

	// size_t numFaces = 0;
	// cout << vertmentions[0] << endl;
	// for(auto const& vert:vertincells[0])
	// 	cout << vert << "  ";
	
	// cout << endl;

	/*Create a colour vector to identify if a face has been identified*/
	vector<vector<size_t>> colour(cdata.elems.size());
	for(size_t ii = 0; ii < nElem; ++ii)
	{
		if(cdata.elems[ii].size() == 8)
		{
			colour[ii] = vector<size_t>(6,0);
		}
		else if(cdata.elems[ii].size() == 6)
		{
			colour[ii] = vector<size_t>(5,0);
		}
		else if(cdata.elems[ii].size() == 5)
		{
			colour[ii] = vector<size_t>(5,0);
		}
		else if(cdata.elems[ii].size() == 4)
		{
			colour[ii] = vector<size_t>(4,0);
		}
	}

	#ifdef DEBUG
	dbout << "Starting main loop..." << endl;
	#endif

	// size_t cellSum=0;
	#pragma omp parallel shared(vertincells)
	{
		/*Create local copy of the face data*/
		FACE flocal;
		// size_t cellCount = 0;
		// size_t reported_Count = 0;
		// size_t localCount = 0;

		#pragma omp for schedule(static) nowait 
		for(size_t lindex = 0; lindex < nElem; ++lindex)
		{	
			vector<size_t> lcell = cdata.elems[lindex];
			vector<vector<size_t>> lfaces;
			size_t lgeom=0;
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

	/* Put the surface faces into the main structures */
	fdata.trig_faces.insert(fdata.trig_faces.end(),fdata.surf_trigs.begin(),fdata.surf_trigs.end());
	fdata.quad_faces.insert(fdata.quad_faces.end(),fdata.surf_quads.begin(),fdata.surf_quads.end());
	fdata.tcelllr.insert(fdata.tcelllr.end(),fdata.strigslr.begin(),fdata.strigslr.end());
	fdata.qcelllr.insert(fdata.qcelllr.end(),fdata.squadslr.begin(),fdata.squadslr.end());

	if (fdata.trig_faces.size() != fdata.ntFaces)
	{
		cout << "numFaces is not being measured correctly." << endl;
		cout << "faces vector size: "<< fdata.trig_faces.size() << "  numFaces: " << fdata.ntFaces << endl;
	}

	if (fdata.quad_faces.size() != fdata.nqFaces)
	{
		cout << "number of quad faces is not being measured correctly." << endl;
		cout << "faces vector size: "<< fdata.quad_faces.size() << "  numFaces: " << fdata.nqFaces << endl;
	}

	if(fdata.tcelllr.size() != fdata.ntFaces)
	{
		cout << "Not all triangle faces have left and right cells identified." << endl;
		cout << "Celllr size: " << fdata.tcelllr.size();
		cout << "  N_Faces: " << fdata.ntFaces << endl;
	}

	if(fdata.qcelllr.size() != fdata.nqFaces)
	{
		cout << "Not all quadrilateral faces have left and right cells identified." << endl;
		cout << "Celllr size: " << fdata.qcelllr.size();
		cout << "  N_Faces: " << fdata.nqFaces << endl;
	}

	#ifdef DEBUG
	dbout << "Building faces complete. Number of faces: " << fdata.nFaces << endl;
	dbout << "Average number of faces per element: " << float(fdata.nFaces)/float(fdata.nElem) << endl;
	dbout << "Exiting BuildFaces..." << endl;
	#endif

	cout << "Building faces complete. Number of faces: " << fdata.nFaces << endl;
	cout << "Average number of faces per element: " << float(fdata.nFaces)/float(fdata.nElem) << endl;
	cout << "Number of surface faces: " << fdata.nSurf << endl;
}


void Write_Face_Data(string const& meshIn, FACE const& fdata)
{
	string meshOut = meshIn;
	meshOut.append(".faces");

	#ifdef DEBUG
	dbout << "Entering Write_Face_Data..." << endl;
	dbout << "Output file: " << meshOut << endl;
	#endif

	cout << "Attempting write output file." << endl;
	cout << "File: " << meshOut << endl;
	
	/* Create the file. */
	int retval;
	int meshID;
	if ((retval = nc_create(meshOut.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &meshID)))
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
	int nElemID, nFaceID, nTFaceID, nPpTFcID, nQFaceID, nPpQFcID, nMarkID, nPntsID;

	if ((retval = nc_def_dim(meshID, "no_of_elements", fdata.nElem, &nElemID)))
	{
		cout << "Error: Failed to define \"no_of_elements\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_faces", fdata.nFaces, &nFaceID)))
	{
		cout << "Error: Failed to define \"no_of_faces\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if(fdata.ntFaces != 0)
	{
		if ((retval = nc_def_dim(meshID, "no_of_triangles", fdata.ntFaces, &nTFaceID)))
		{
			cout << "Error: Failed to define \"no_of_triangles\"" << endl;
			ERR(retval);
			exit(-1);
		}

		if ((retval = nc_def_dim(meshID, "points_per_triangle", 3, &nPpTFcID)))
		{
			cout << "Error: Failed to define \"points_per_triangle\"" << endl;
			ERR(retval);
			exit(-1);
		}	
	}

	if(fdata.nqFaces != 0)
	{
		if ((retval = nc_def_dim(meshID, "no_of_quadrilaterals", fdata.nqFaces, &nQFaceID)))
		{
			cout << "Error: Failed to define \"no_of_quadrilaterals\"" << endl;
			ERR(retval);
			exit(-1);
		}

		if ((retval = nc_def_dim(meshID, "points_per_quadrilateral", 4, &nPpQFcID)))
		{
			cout << "Error: Failed to define \"points_per_quadrilateral\"" << endl;
			ERR(retval);
			exit(-1);
		}	
	}

	// if ((retval = nc_def_dim(meshID, "no_of_surfacequadrilaterals", fdata.nWall, &nWallID)))
	// {
	// 	cout << "Error: Failed to define \"no_of_surfacequadrilaterals\"" << endl;
	// 	ERR(retval);
	// 	exit(-1);
	// }

	// if ((retval = nc_def_dim(meshID, "no_of_surfacetriangles", fdata.nWall, &nWallID)))
	// {
	// 	cout << "Error: Failed to define \"no_of_surfacetriangles\"" << endl;
	// 	ERR(retval);
	// 	exit(-1);
	// }
	
	if ((retval = nc_def_dim(meshID, "no_of_surfaceelements", fdata.nSurf, &nMarkID)))
	{
		cout << "Error: Failed to define \"no_of_surfaceelements\"" << endl;
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_def_dim(meshID, "no_of_points", fdata.nPnts, &nPntsID)))
	{
		cout << "Error: Failed to define \"no_of_points\"" << endl;
		ERR(retval);
		exit(-1);
	}

	/* Define the variables */
	int faceTVarID, faceQVarID, leftVarID, rightVarID, markID, ptsxID, ptsyID, ptszID;
	
	/*Define the faces*/
	int dimIDs[] = {nTFaceID,nPpTFcID};

	if(fdata.ntFaces != 0)
	{
		if ((retval = nc_def_var(meshID, "points_of_triangles", NC_INT, 2,
								dimIDs, &faceTVarID)))
		{
			cout << "Error: Failed to define \"points_of_triangles\"" << endl;
			ERR(retval);
			exit(-1);
		}
	}

	if(fdata.nqFaces != 0)
	{
		dimIDs[0] = nQFaceID;
		dimIDs[1] = nPpQFcID;
		if ((retval = nc_def_var(meshID, "points_of_quadrilaterals", NC_INT, 2,
								dimIDs, &faceQVarID)))
		{
			cout << "Error: Failed to define \"points_of_quadrilaterals\"" << endl;
			ERR(retval);
			exit(-1);
		}	
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

	if ((retval = nc_def_var(meshID, "boundarymarker_of_surfaces", NC_INT, 1,
							 &nMarkID, &markID)))
	{
		cout << "Error: Failed to define \"boundarymarker_of_surfaces\"" << endl;
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
	int* faces = new int[fdata.ntFaces*3];
	size_t start[] = {0,0};
	size_t end[] = {fdata.ntFaces,3};

	if(fdata.ntFaces != 0)
	{
		for(size_t ii = 0; ii < fdata.ntFaces; ++ii)
			for(size_t jj = 0; jj < 3; ++jj)
			{
				faces[index(ii,jj,3)] = static_cast<int>(fdata.trig_faces[ii][jj]);
			}

		/*Put faces into the file*/
		
		if ((retval = nc_put_vara_int(meshID, faceTVarID, start, end, &faces[0])))
		{
			cout << "Failed to write triangle face data" << endl;
			ERR(retval);
			exit(-1);
		}	

	}

	if(fdata.nqFaces != 0)
	{
		/*Create the C array for the faces*/
		faces = new int[fdata.nqFaces*4];
		for(size_t ii = 0; ii < fdata.nqFaces; ++ii)
			for(size_t jj = 0; jj < 4; ++jj)
			{
				faces[index(ii,jj,4)] = static_cast<int>(fdata.quad_faces[ii][jj]);
			}

		/*Put faces into the file*/
		end[0] = fdata.nqFaces;
		end[1] = 4;
		if ((retval = nc_put_vara_int(meshID, faceQVarID, start, end, &faces[0])))
		{
			cout << "Failed to write quadrilateral face data" << endl;
			ERR(retval);
			exit(-1);
		}
	}

	// elemFaces.putVar(faces);

	// /*State how many faces there are per element*/
	// int* nFacesPElem = new int[fdata.numElem];
	// for(size_t ii = 0; ii < fdata.numElem; ++ii)
	// 	nFacesPElem[ii] = 0; 

	// for(size_t ii = 0; ii < fdata.numFaces; ++ii)
	// {	/*Count the mentions of a cell*/
	// 	nFacesPElem[fdata.celllr[ii].first]++;
	// }

	// nElemFaces.putVar(fdata.nFacesPElem);

	/*Put face left and right into the file*/
	int* left = new int[fdata.nFaces];
	int* right = new int[fdata.nFaces];

	for(size_t ii = 0; ii < fdata.ntFaces; ++ii)
	{
		left[ii] = fdata.tcelllr[ii].first;
		right[ii] = fdata.tcelllr[ii].second;
	}

	for(size_t ii = 0; ii < fdata.nqFaces; ++ii)
	{
		left[ii+fdata.ntFaces] = fdata.qcelllr[ii].first;
		right[ii+fdata.ntFaces] = fdata.qcelllr[ii].second;
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

	/* Write the boundary markers */
	int* markers = new int[fdata.nSurf];
	size_t ntSurf = fdata.tsmarkers.size();
	for(size_t ii = 0; ii< ntSurf; ++ii )
	{
		markers[ii] = fdata.tsmarkers[ii];
	}

	size_t nqSurf = fdata.qsmarkers.size();
	for(size_t ii = 0; ii< nqSurf; ++ii )
	{
		markers[ii+ntSurf] = fdata.qsmarkers[ii];
	}

	if ((retval = nc_put_var_int(meshID, markID, &markers[0])))
	{
		cout << "Failed to write surface marker data" << endl;
		ERR(retval);
		exit(-1);
	}

	/*Create the C arrays for the vertices*/
	double* x = new double[fdata.nPnts];
	double* y = new double[fdata.nPnts];
	double* z = new double[fdata.nPnts];

	for(size_t ii = 0; ii < fdata.nPnts; ++ii)
	{
		x[ii] = fdata.verts[ii][0];
		y[ii] = fdata.verts[ii][1];
		z[ii] = fdata.verts[ii][2];
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

	size_t TotalNumFaceNodes = 0;
	for(size_t ii = 0; ii < fdata.trig_faces.size(); ++ii)
	{
		TotalNumFaceNodes += fdata.trig_faces[ii].size();
	}

	for(size_t ii = 0; ii < fdata.quad_faces.size(); ++ii)
	{
		TotalNumFaceNodes += fdata.quad_faces[ii].size();
	}

	fout << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
	fout << "ZONE T=\"FEPOLYHEDRON Test\"" << endl;
	fout << "ZONETYPE=FEPOLYHEDRON" << endl;
	fout << "NODES=" << fdata.nPnts << " ELEMENTS=" << fdata.nElem << " FACES=" << fdata.nFaces << endl;
	fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
	fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

	size_t w = 15;
	size_t preci = 6;
	fout << std::left << std::scientific << std::setprecision(preci);
	// fout << fdata.numElem << " " << fdata.numFaces << "  " <<  fdata.numPoint << endl;
	
	/*Write vertices in block format (Each dimension in turn)*/
	size_t newl = 0;
	fout << std::setw(1);
	for(size_t DIM = 0; DIM < 3; ++DIM)
	{
		for(size_t ii = 0; ii < fdata.verts.size(); ++ii)
		{
			fout << std::setw(w) << fdata.verts[ii][DIM];
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
	for (size_t ii = 0; ii < fdata.trig_faces.size(); ++ii)
	{
		fout << std::setw(w) << fdata.trig_faces[ii].size();
		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}

	for (size_t ii = 0; ii < fdata.quad_faces.size(); ++ii)
	{
		fout << std::setw(w) << fdata.quad_faces[ii].size();
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
	for (size_t ii = 0; ii < fdata.trig_faces.size(); ++ii)
	{
		for(auto const& vertex:fdata.trig_faces[ii])
		{	/*Write face vertex indexes*/
			fout << std::setw(w) << vertex+1;
			// if (vertex > fdata.nPnts)
			// {
			// 	cout << "Trying to write a vertex outside of the number of points." << endl;
			// }
		}
		fout << endl;
	}

	for (size_t ii = 0; ii < fdata.quad_faces.size(); ++ii)
	{
		for(auto const& vertex:fdata.quad_faces[ii])
		{	/*Write face vertex indexes*/
			fout << std::setw(w) << vertex+1;
			// if (vertex > fdata.nPnts)
			// {
			// 	cout << "Trying to write a vertex outside of the number of points." << endl;
			// }
		}
		fout << endl;
	}

	/*Write face left and right*/
	newl = 0;
	fout << "#left elements" << endl;
	for (size_t ii = 0; ii < fdata.tcelllr.size(); ++ii)
	{
		fout << std::setw(w) << fdata.tcelllr[ii].first+1 ;
		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}

	for (size_t ii = 0; ii < fdata.qcelllr.size(); ++ii)
	{
		fout << std::setw(w) << fdata.qcelllr[ii].first+1 ;
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
	for (size_t ii = 0; ii < fdata.tcelllr.size(); ++ii)
	{
		if(fdata.tcelllr[ii].second < 0)
			fout<< std::setw(w) << 0 ;
		else
			fout << std::setw(w) << fdata.tcelllr[ii].second+1;


		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}

	for (size_t ii = 0; ii < fdata.qcelllr.size(); ++ii)
	{
		if(fdata.qcelllr[ii].second < 0)
			fout<< std::setw(w) << 0 ;
		else
			fout << std::setw(w) << fdata.qcelllr[ii].second+1;


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
	fout << "N=" << cdata.nPnts << ", E=" << cdata.nElem << 
	", F=FEBLOCK, ET=BRICK"  << endl << endl;

	/*Write vertices*/
	fout << std::left << std::scientific << std::setprecision(6);
	fout << std::setw(1);
	for(size_t ii = 0; ii < 3; ++ii)
	{	
		size_t kk = 0;
		for(size_t jj = 0; jj < cdata.verts.size(); ++jj)
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
	for(size_t ii = 0; ii < cdata.elems.size(); ++ii)
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

	size_t w = 12;
	fout << std::left << std::fixed;
	fout << setw(w) << fdata.nElem << setw(w) << 
		fdata.nFaces << setw(w) << fdata.nPnts << endl;

	for(size_t ii = 0; ii < fdata.ntFaces; ++ii)
	{
		for(auto const& vert:fdata.trig_faces[ii])
			fout << setw(w) << vert;

		fout << setw(w) << fdata.tcelllr[ii].first << setw(w) << fdata.tcelllr[ii].second << endl;
	}

	for(size_t ii = 0; ii < fdata.nqFaces; ++ii)
	{
		for(auto const& vert:fdata.quad_faces[ii])
			fout << setw(w) << vert;

		fout << setw(w) << fdata.qcelllr[ii].first << setw(w) << fdata.qcelllr[ii].second << endl;
	}

	fout << std::scientific << std::setprecision(6);
	w = 15;
	for(size_t ii = 0; ii < fdata.nPnts; ++ii)
	{
		fout << setw(w) << ii;
		for(size_t jj = 0; jj < 3; ++jj)
		{
			fout << setw(w) << fdata.verts[ii][jj];
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
	vector<vector<size_t>> tets; 
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
	vector<vector<size_t>> prism;
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
	vector<vector<size_t>> pyra;
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
	vector<vector<size_t>> hexa;
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
	cdata.verts = Get_Coordinate_Vector(meshID, nPnts);

	if (cdata.verts.size() != nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Get surface faces*/
	// vector<int> markers = Find_Bmap_Markers(bmapIn);
	Get_Surface(meshID, cdata);	

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
	// Write_ASCII_Face_Data(fdata);
	// Write_Griduns(fdata);
	return 0;
}