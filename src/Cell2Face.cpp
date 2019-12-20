/*Cell based to Face based data converter*/

#include <vector>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2

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
	uint numPoint, numElem;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Cell based data*/
	vector<vector<uint>> elems;
}CELL;

typedef struct FACE
{
	FACE(const vector<StateVecD>& cverts)
	{
		verts = cverts;
		numPoint = cverts.size();
	}	

	uint numElem, numPoint, nfaces;
	vector<StateVecD> verts;
	vector<vector<uint>> faces; /*Face indexes*/
	vector<std::pair<int,int>> celllr; /*Cell left and right of the face*/

	real* x;
	real* y;
	real* z;
	int* left;
	int* right;
	int** faceindex;
}FACE;

uint index(uint ii, uint jj, uint nPts)
{
	return(ii*nPts + jj);
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
		cout << "Number of cdata: " << nElem << endl;

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

void CheckFaces(const vector<uint>& vertmentions, const vector<vector<uint>>& vertincells,
	vector<vector<uint>>& colour, vector<vector<uint>>& lfaces, uint lindex, 
	const CELL& cdata, FACE& fdata)
{
	uint lfaceindex = 0;
	for(auto const& face:lfaces)
	{	
		uint match = 0;
/*Search through the cells where the first vertex is mentioned, and see if this face exists*/
		for (uint jj = 0; jj < vertmentions[face[0]]; ++jj)
		{	/*If the cell is the current cell, ignore*/
			uint rindex = vertincells[face[0]][jj];
			if(rindex!=lindex)
			{	
				vector<uint> rcell = cdata.elems[rindex];

				/*Create an ordered list of the left side of the face*/
				vector<uint> lface = face;
				std::sort(lface.begin(),lface.end());

				uint rfaceindex = 0;

				if(rcell.size() == 4)
				{	/*Check if cell is a tet, since it cant share a face*/
					break;
				}

				if(rcell.size() == 6)
				{	/*Prism*/
					vector<vector<uint>> rfaces = {{rcell[1],rcell[4],rcell[5],rcell[2]},
										   		   {rcell[0],rcell[1],rcell[4],rcell[3]},
										   		   {rcell[0],rcell[3],rcell[5],rcell[2]}};

					for(auto rface:rfaces)
					{	
						std::sort(rface.begin(),rface.end());
						if(lface == rface)
						{
							match = 1;
							break;
						}
						rfaceindex++;
					}
				}

				if(rcell.size() == 5)
				{	/*Pyra*/
					vector<uint> rface = {rcell[0],rcell[3],rcell[2],rcell[1]};
										  
					std::sort(rface.begin(),rface.end());
					if(lface == rface)
					{
						match = 1;
					}							
				}

				if(rcell.size() == 8)
				{
					vector<vector<uint>> rfaces ={{rcell[0],rcell[1],rcell[2],rcell[3]},
												  {rcell[4],rcell[0],rcell[3],rcell[7]},
												  {rcell[1],rcell[0],rcell[4],rcell[5]},
												  {rcell[2],rcell[1],rcell[5],rcell[6]},
												  {rcell[3],rcell[2],rcell[6],rcell[7]},
												  {rcell[7],rcell[6],rcell[5],rcell[4]}};

					for(auto rface:rfaces)
					{
						std::sort(rface.begin(),rface.end());
						if(lface == rface)
						{
							match = 1;
							break;
						}
						rfaceindex++;
					}						 
				}


				if(match == 1)
				{	/*Then the face is a match*/
					fdata.nfaces++;
					fdata.celllr.emplace_back(std::pair<int,int>(lindex,rindex));
					fdata.celllr.emplace_back(std::pair<int,int>(lindex,rindex));
					/*Break the face down into two triangles*/
					vector<uint> face1 = {face[0],face[1],face[2]};
					vector<uint> face2 = {face[0],face[2],face[3]};
					fdata.faces.emplace_back(face1);
					fdata.faces.emplace_back(face2);
					colour[lindex][lfaceindex] = 1;
					colour[rindex][rfaceindex] = 1;
				}
			}
		}
		

		if(match==0)
		{
			if(colour[lindex][lfaceindex] == 0)
			{
				fdata.nfaces++;
				/*Break the face down into two triangles*/
				vector<uint> face1 = {face[0],face[1],face[2]};
				vector<uint> face2 = {face[0],face[2],face[3]};
				fdata.faces.emplace_back(face1);
				fdata.faces.emplace_back(face2);
				colour[lindex][lfaceindex] = 1;
				/*Hexes will always be next to the surface geometry...*/
				fdata.celllr.emplace_back(std::pair<int,int>(lindex,-1));
				fdata.celllr.emplace_back(std::pair<int,int>(lindex,-1));
			}
		}

		lfaceindex++;
	}
}

vector<vector<uint>> BuildFaces(const CELL& cdata, FACE& fdata)
{
	uint nElem = cdata.numElem;
	uint nPts = cdata.numPoint;
	/*Build the faces, and check there are no duplicate faces*/
	/*Have a list of cells that have unidentified faces, that could have this face*/

	/*Vector of how many times a vertex has been mentioned (ie how many cells is it in)*/
	vector<uint> vertmentions(nPts,0); 

	/*Vector of which cells the vertex is referenced in*/
	vector<vector<uint>> vertincells(nPts);

	for(uint ii = 0; ii < nElem; ++ii)
	{
		auto const& cell = cdata.elems[ii];
		for(auto const& vert:cell)
		{	/*Add count for every mention of the vertex, and push back the cell id*/
			vertmentions[vert] +=1;
			vertincells[vert].emplace_back(ii);
		}
	}

	uint nfaces = 0;

	vector<vector<uint>> colour(cdata.elems.size());
	for(uint ii = 0; ii < nElem; ++ii)
		colour[ii].emplace_back(vector<uint>(cdata.elems[ii].size(),0));


	for(uint lindex = 0; lindex < nElem; ++lindex)
	{	
		vector<uint> lcell = cdata.elems[lindex];

		/*Check size of the inner vector, as that will show which geometry it is*/
		if(lcell.size() == 8)
		{	/*Hex*/
			/*Create faces*/
			
			vector<vector<uint>> faces = {{lcell[0],lcell[1],lcell[2],lcell[3]},
										  {lcell[4],lcell[0],lcell[3],lcell[7]},
										  {lcell[1],lcell[0],lcell[4],lcell[5]},
										  {lcell[2],lcell[1],lcell[5],lcell[6]},
										  {lcell[3],lcell[2],lcell[6],lcell[7]},
										  {lcell[7],lcell[6],lcell[5],lcell[4]}};

			CheckFaces(vertmentions, vertincells, colour, faces, lindex, cdata, fdata);
		}	/*End of hex check*/

		if(lcell.size() == 6)
		{	/*Prism*/
			vector<vector<uint>> faces = {{lcell[1],lcell[4],lcell[5],lcell[2]},
								   		  {lcell[0],lcell[1],lcell[4],lcell[3]},
								   		  {lcell[0],lcell[3],lcell[5],lcell[2]},
								   		  {lcell[0],lcell[1],lcell[2]},
								   		  {lcell[5],lcell[4],lcell[3]}};

			CheckFaces(vertmentions, vertincells, colour, faces, lindex, cdata, fdata);
		}
		if(lcell.size() == 5)
		{	/*Pyra*/
			vector<vector<uint>> faces = {{lcell[0],lcell[3],lcell[2],lcell[1]},
								   		  {lcell[0],lcell[1],lcell[4]},
								   		  {lcell[2],lcell[4],lcell[1]},
								   		  {lcell[3],lcell[4],lcell[2]},
								   		  {lcell[0],lcell[4],lcell[3]}};

			CheckFaces(vertmentions, vertincells, colour, faces, lindex, cdata, fdata);
		}
		if(lcell.size() == 4)
		{	/*Tet*/
			vector<vector<uint>> faces = {{lcell[0],lcell[1],lcell[2]},
								   		  {lcell[0],lcell[2],lcell[3]},
								   		  {lcell[0],lcell[3],lcell[1]},
								   		  {lcell[3],lcell[2],lcell[1]}};

			CheckFaces(vertmentions, vertincells, colour, faces, lindex, cdata, fdata);
		}
	}

	return vector<vector<uint>>(0);
}

int main (int argc, char** argv)
{
	/*Idea: Take TAU cell based mesh, and convert to a face based data in NetCDF or TECIO*/
	string meshIn = argv[1];

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

	uint nfaces = 0;
	for(auto const& cell: cdata.elems)
	{
		if(cell.size() == 8)
			nfaces+=12; /*Hex*/
		if(cell.size() == 6)
			nfaces+=8;	/*Prism*/
		if(cell.size() == 5)
			nfaces+=6;	/*Pyra*/
		if(cell.size() == 4)
			nfaces+=4;	/*Tet*/
	}

	/*Now build the face based data*/
	FACE fdata(cdata.verts);


	return 0;
}