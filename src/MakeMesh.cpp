// Make mesh

#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include <omp.h>
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <dirent.h>

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

#ifdef DEBUG
	/*Open debug file to write to*/	
	std::ofstream dbout("MakeMesh.log",std::ios::out);
#endif

using std::vector;
using std::pair;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string; 
using std::setw;

typedef unsigned int uint;
typedef double real;

uint DIM = 2;


/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real,Eigen::Dynamic,1,0,3,1> StateVecD;
typedef Eigen::Matrix<size_t,Eigen::Dynamic,1,0,3,1> StateVecI;
typedef Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic,0,3,3> RotMat;

struct MESH
{
	void make(size_t nV, size_t nE, size_t nF)
	{
		nVerts = nV;
		nElems = nE;
		nFaces = nF;
		nWall = 0;
		nFar = 0;

		x = new real[nVerts];
		y = new real[nVerts];

		if(DIM == 2)
		{
			faces = new int[nFaces*2];
		}
		else 
		{
			z = new real[nVerts];
			faces = new int[nFaces*3];
		}	
		
		left = new int[nFaces];
		right = new int[nFaces];

		vx = new real[nVerts];
		vy = new real[nVerts];
		vz = new real[nVerts];
		p = new real[nVerts];	
		rho = new real[nVerts];
	}

	size_t nVerts, nElems, nFaces, nWall, nFar;

	real* x;
	real* y;
	real* z;

	int* faces;
	int* left;
	int* right;

	real* vx;
	real* vy;
	real* vz;

	real* p;
	real* rho;
};

struct SIM
{
	string infolder,meshname,solname;
	StateVecD start, end;
	StateVecI ijk;
	StateVecD dx;
	StateVecD vel;
	real p, rho;
};

size_t index(size_t ii, size_t jj, size_t nI)
{
	return(ii + jj*nI);
}

size_t index(size_t ii, size_t jj, size_t kk, size_t nI, size_t nJ)
{
	return  (kk * nJ + jj) * nI + ii;
}

real getDouble(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	std::stringstream sstr;
	sstr << line;
	real d;
	if(sstr >> d)
	{
#ifdef DEBUG
		dbout << name << ": " << d << endl;
#endif	
		return d; 
	}
	else
	{
		cout << "Line does not contain a value. Please check your input." << endl;
		cout << "Expecting: " << name << endl;
	}
		
	return 0;
}

std::string getString(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	size_t ptr = line.find_first_of(' ');
	string result = line.substr(0,ptr);
#ifdef DEBUG
	dbout << name << ": " << result << endl;
#endif	
	return result; 
}

template<typename T>
Eigen::Matrix<T,Eigen::Dynamic,1> getVector(ifstream& In, string const& name, uint& lineno)
{
	string line;
	getline(In,line);
	lineno++;
	std::stringstream sline(line);
	
	vector<T> vec;
	real value;
	string temp;
	for(uint ii = 0; ii < 3; ii++)
	{
		sline >> temp;
		if(std::stringstream(temp)>> value)
		{
			vec.emplace_back(value);
		}
		temp = "";
	}

	if(vec.size() == 0)
	{
		cout << "Error: Dimension is zero. Check the input order." << endl;
		cout << "Line " << lineno << ": " << line << endl; 
		exit(-1);
	}

	DIM = vec.size();

	Eigen::Matrix<T,Eigen::Dynamic,1> x(DIM);

	for(size_t ii = 0; ii < DIM; ii++)
	{
		x(ii) = vec[ii];
	}

#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1);
	if(DIM == 3)
		dbout << "  " << x(2);
	dbout << endl;
#endif	

	return x;
}


void GetInput(int argc, char **argv, SIM& svar)
{
	StateVecD angle;
	if (argc > 2) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only one input argument accepted,\n";
		cout << "1: Input file directoy.\n";
		cout << "Other inputs will be ignored." << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tERROR: No inputs provided. Stopping... \n";
    	exit(-1);    	
    }


	/*Get parameters if it has been provided*/
	// cout << argv[1] << endl;
	string file = argv[1];
	    	
	if(file.back() != '/')
    	file.append("/");

    svar.infolder = file;
	file.append("Preprocessing");
#ifdef DEBUG
	dbout << "Reading settings file. Path:" << endl << file << endl;
#endif
	std::ifstream in(file);
  	if(!in.is_open()) 
  	{	
  		cerr << "Error opening the settings file." << endl;
	    exit(-1);
  	}

  	
  	// First define start coordinate.
  	// Number in i, j, k (if 3D)
  	// spacing for each cell (uniformly)
  	// Or define the end coordinate?

  	uint lineno = 0;
  	svar.start = getVector<real>(in,"Start Vector",lineno);
  	const int dim = DIM;
  	svar.end = getVector<real>(in,"End Vector",lineno);
  	svar.ijk = getVector<size_t>(in,"IJK Values",lineno);
  	svar.vel = getVector<real>(in,"Mesh velocity",lineno);
  	svar.p = getDouble(in,lineno,"Mesh Pressure");
  	svar.rho = getDouble(in,lineno,"Mesh Density");
  	svar.meshname = getString(in,lineno,"Mesh filename");
  	svar.solname = getString(in,lineno,"Solution filename");

  	if (svar.end.size()!= dim || svar.ijk.size() != dim || svar.vel.size() != dim)
  	{
  		cout << "Mismatch of dimensions in file. Please check your input parameters." << endl;
  		exit(-1);
  	}

  	cout << "Dimension: " << DIM << endl;
  	// Calculate the steps:
  	StateVecD dx(dim);

  	dx(0) = (svar.end(0)-svar.start(0))/(1.0*(svar.ijk(0)-1));
  	dx(1) = (svar.end(1)-svar.start(1))/(1.0*(svar.ijk(1)-1));
  	if(dim == 3)
	  	dx(2) = (svar.end(2)-svar.start(2))/(1.0*(svar.ijk(2)-1));

	cout << "Cell dimensions: " << dx(0) << "  " << dx(1);
	if(dim == 3)
		cout << "  " << dx(2);
	cout << endl;

	real cVol = dx(0) * dx(1);
	if(dim == 3)
		cVol *= dx(2);

	cout << "Cell volume: " << abs(cVol) << endl;

	svar.dx = dx;
} 

void Create_Mesh(SIM const& svar, MESH& cells)
{
	size_t const& nI = svar.ijk(0);
	size_t const& nJ = svar.ijk(1);
	size_t nK = 1;
	if(DIM == 3)
	{
		nK = svar.ijk(2);
	}

	StateVecD const& dx = svar.dx;

	size_t nVerts = nI*nJ;
	if (DIM == 3)
		nVerts*=nK;

	size_t nElems = (nI-1)*(nJ-1);
	if (DIM == 3)
		nElems*=(nK-1);
	
	size_t nFaces = 0;
	if(DIM == 2)
	{
		nFaces = 2*nI*nJ -nI -nJ;
	}
	else if (DIM == 3)
	{
		nFaces = nK*(nI-1)*(nJ-1) + nI*(nJ-1)*(nK-1)+ nJ*(nI-1)*(nK-1);
		nFaces *= 2;
	}
	
	cells.make(nVerts,nElems,nFaces);

	if(DIM == 3)
	{
		cout << "Building vertices..." << endl;
		for(size_t kk = 0; kk < nK; kk++)
		{
			for(size_t jj = 0; jj < nJ; jj++)
			{
				for(size_t ii = 0; ii < nI; ii++)
				{
					cells.x[index(ii,jj,kk,nI,nJ)] = svar.start(0)+ii*dx(0);
					cells.y[index(ii,jj,kk,nI,nJ)] = svar.start(1)+jj*dx(1);
					cells.z[index(ii,jj,kk,nI,nJ)] = svar.start(2)+kk*dx(2);

					cells.vx[index(ii,jj,kk,nI,nJ)] = svar.vel(0);
					cells.vy[index(ii,jj,kk,nI,nJ)] = svar.vel(1);
					cells.vz[index(ii,jj,kk,nI,nJ)] = svar.vel(2);
					cells.p[index(ii,jj,kk,nI,nJ)] = svar.p;
					cells.rho[index(ii,jj,kk,nI,nJ)] = svar.rho;
				}
			}
		} 

		// Create the faces
		cout << "Building faces..." << endl;
		size_t faceno = 0;

		// Do the front and back faces
		for(size_t jj = 0; jj < nJ-1; jj++)
		{
			for(size_t ii = 0; ii < nI-1; ii++)
			{
				cells.faces[faceno*3+0] = index(ii,jj,0,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,jj,0,nI,nJ);
				cells.faces[faceno*3+2] = index(ii+1,jj+1,0,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,jj,0,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;

				cells.faces[faceno*3+0] = index(ii,jj,0,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,jj+1,0,nI,nJ);
				cells.faces[faceno*3+2] = index(ii,jj+1,0,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,jj,0,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;


				cells.faces[faceno*3+0] = index(ii,jj,nK-1,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,jj,nK-1,nI,nJ);
				cells.faces[faceno*3+2] = index(ii+1,jj+1,nK-1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,jj,nK-2,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;

				cells.faces[faceno*3+0] = index(ii,jj,nK-1,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,jj+1,nK-1,nI,nJ);
				cells.faces[faceno*3+2] = index(ii,jj+1,nK-1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,jj,nK-2,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;
			}
		}

		// Do the bottom and top faces
		for(size_t kk = 0; kk < nK-1; kk++)
		{
			for(size_t ii = 0; ii < nI-1; ii++)
			{
				cells.faces[faceno*3+0] = index(ii,0,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,0,kk,nI,nJ);
				cells.faces[faceno*3+2] = index(ii+1,0,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,0,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;

				cells.faces[faceno*3+0] = index(ii,0,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,0,kk+1,nI,nJ);
				cells.faces[faceno*3+2] = index(ii,0,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,0,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;


				cells.faces[faceno*3+0] = index(ii,nJ-1,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,nJ-1,kk,nI,nJ);
				cells.faces[faceno*3+2] = index(ii+1,nJ-1,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,nJ-2,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;

				cells.faces[faceno*3+0] = index(ii,nJ-1,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(ii+1,nJ-1,kk+1,nI,nJ);
				cells.faces[faceno*3+2] = index(ii,nJ-1,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(ii,nJ-2,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;
			}
		}


		// Do the left and right faces
		for(size_t jj = 0; jj < nJ-1; jj++)
		{
			for(size_t kk = 0; kk < nK-1; kk++)
			{
				cells.faces[faceno*3+0] = index(0,jj,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(0,jj+1,kk,nI,nJ);
				cells.faces[faceno*3+2] = index(0,jj+1,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(0,jj,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;

				cells.faces[faceno*3+0] = index(0,jj,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(0,jj+1,kk+1,nI,nJ);
				cells.faces[faceno*3+2] = index(0,jj,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(0,jj,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;


				cells.faces[faceno*3+0] = index(nI-1,jj,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(nI-1,jj+1,kk,nI,nJ);
				cells.faces[faceno*3+2] = index(nI-1,jj+1,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(nI-2,jj,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;

				cells.faces[faceno*3+0] = index(nI-1,jj,kk,nI,nJ);
				cells.faces[faceno*3+1] = index(nI-1,jj+1,kk+1,nI,nJ);
				cells.faces[faceno*3+2] = index(nI-1,jj,kk+1,nI,nJ);
				cells.left[faceno] = static_cast<int>(index(nI-2,jj,kk,nI-1,nJ-1));
				cells.right[faceno] = -2;
				cells.nFar++;
				faceno++;
			}
		}

		

		// Do faces in i
		for(size_t ii = 1; ii < nI-1; ii++)
		{
			for(size_t jj = 0; jj < nJ-1; jj++)
			{
				for(size_t kk = 0; kk < nK-1; kk++)
				{
					cells.faces[faceno*3+0] = index(ii,jj,kk,nI,nJ);
					cells.faces[faceno*3+1] = index(ii,jj+1,kk,nI,nJ);
					cells.faces[faceno*3+2] = index(ii,jj+1,kk+1,nI,nJ);
					cells.left[faceno] = static_cast<int>(index(ii-1,jj,kk,nI-1,nJ-1));
					cells.right[faceno] = static_cast<int>(index(ii,jj,kk,nI-1,nJ-1));
					faceno++;

					cells.faces[faceno*3+0] = index(ii,jj,kk,nI,nJ);
					cells.faces[faceno*3+1] = index(ii,jj+1,kk+1,nI,nJ);
					cells.faces[faceno*3+2] = index(ii,jj,kk+1,nI,nJ);
					cells.left[faceno] = static_cast<int>(index(ii-1,jj,kk,nI-1,nJ-1));
					cells.right[faceno] = static_cast<int>(index(ii,jj,kk,nI-1,nJ-1));
					faceno++;
				}
			}
		}

		// Do faces in j
		for(size_t jj = 1; jj < nJ-1; jj++)
		{
			for(size_t ii = 0; ii < nI-1; ii++)
			{
				for(size_t kk = 0; kk < nK-1; kk++)
				{
					cells.faces[faceno*3+0] = index(ii,jj,kk,nI,nJ);
					cells.faces[faceno*3+1] = index(ii+1,jj,kk,nI,nJ);
					cells.faces[faceno*3+2] = index(ii+1,jj,kk+1,nI,nJ);
					cells.left[faceno] = static_cast<int>(index(ii,jj-1,kk,nI-1,nJ-1));
					cells.right[faceno] = static_cast<int>(index(ii,jj,kk,nI-1,nJ-1));
					faceno++;

					cells.faces[faceno*3+0] = index(ii,jj,kk,nI,nJ);
					cells.faces[faceno*3+1] = index(ii+1,jj,kk+1,nI,nJ);
					cells.faces[faceno*3+2] = index(ii,jj,kk+1,nI,nJ);
					cells.left[faceno] = static_cast<int>(index(ii,jj-1,kk,nI-1,nJ-1));
					cells.right[faceno] = static_cast<int>(index(ii,jj,kk,nI-1,nJ-1));
					faceno++;
				}
			}
		}

		// Do faces in k
		for(size_t kk = 1; kk < nK-1; kk++)
		{
			for(size_t ii = 0; ii < nI-1; ii++)
			{
				for(size_t jj = 0; jj < nJ-1; jj++)
				{
					cells.faces[faceno*3+0] = index(ii,jj,kk,nI,nJ);
					cells.faces[faceno*3+1] = index(ii+1,jj,kk,nI,nJ);
					cells.faces[faceno*3+2] = index(ii+1,jj+1,kk,nI,nJ);
					cells.left[faceno] = static_cast<int>(index(ii,jj,kk-1,nI-1,nJ-1));
					cells.right[faceno] = static_cast<int>(index(ii,jj,kk,nI-1,nJ-1));
					faceno++;

					cells.faces[faceno*3+0] = index(ii,jj,kk,nI,nJ);
					cells.faces[faceno*3+1] = index(ii+1,jj+1,kk,nI,nJ);
					cells.faces[faceno*3+2] = index(ii,jj+1,kk,nI,nJ);
					cells.left[faceno] = static_cast<int>(index(ii,jj,kk-1,nI-1,nJ-1));
					cells.right[faceno] = static_cast<int>(index(ii,jj,kk,nI-1,nJ-1));
					faceno++;
				}
			}
		}
	}
	else if (DIM == 2)
	{
		cout << "Building vertices..." << endl;
		for(size_t jj = 0; jj < nJ; jj++)
		{
			for(size_t ii = 0; ii < nI; ii++)
			{
				cells.x[index(ii,jj,nI)] = svar.start(0)+ii*dx(0);
				cells.y[index(ii,jj,nI)] = svar.start(1)+jj*dx(1);

				cells.vx[index(ii,jj,nI)] = svar.vel(0);
				cells.vy[index(ii,jj,nI)] = svar.vel(1);
				cells.p[index(ii,jj,nI)] = svar.p;
				cells.rho[index(ii,jj,nI)] = svar.rho;
			}
		}

		// Create the faces
		
		size_t faceno = 0;
		cout << "Building faces..." << endl;
		// Do the top and bottom edges
		for(size_t ii = 0; ii < nI-1; ii++)
		{
			cells.faces[faceno*2+0] = index(ii,0,nI);
			cells.faces[faceno*2+1] = index(ii+1,0,nI);
			cells.left[faceno] = static_cast<int>(index(ii,0,nI-1));
			cells.right[faceno] = -2;
			cells.nFar++;
			faceno++;

			
			
			cells.faces[faceno*2+0] = index(ii+1,nJ-1,nI);
			cells.faces[faceno*2+1] = index(ii,nJ-1,nI);
			cells.left[faceno] = static_cast<int>(index(ii,nJ-2,nI-1));
			cells.right[faceno] = -2;
			cells.nFar++;
			faceno++;				
		}

		// Do the left & right edges
		for(size_t jj = 0; jj < nJ-1; jj++)
		{
			cells.faces[faceno*2+0] = index(0,jj+1,nI);
			cells.faces[faceno*2+1] = index(0,jj,nI);
			cells.left[faceno] = static_cast<int>(index(0,jj,nI-1));
			cells.right[faceno] = -2;
			cells.nFar++;
			faceno++;

			cells.faces[faceno*2+0] = index(nI-1,jj,nI);
			cells.faces[faceno*2+1] = index(nI-1,jj+1,nI);
			cells.left[faceno] = static_cast<int>(index(nI-2,jj,nI-1));
			cells.right[faceno] = -2;
			cells.nFar++;
			faceno++;
						
		}

		// Do the horizontal faces
		for(size_t jj = 1; jj < nJ-1; jj++)
			for(size_t ii = 0; ii < nI-1; ii++)
			{
				cells.faces[faceno*2+0] = index(ii,jj,nI);
				cells.faces[faceno*2+1] = index(ii+1,jj,nI);
				cells.left[faceno] = static_cast<int>(index(ii,jj,nI-1));
				cells.right[faceno] = static_cast<int>(index(ii,jj-1,nI-1));
				faceno++;
			}

		// Do the vertical faces
		for(size_t ii = 1; ii < nI-1; ii++)
			for(size_t jj = 0; jj < nJ-1; jj++)
			{
				cells.faces[faceno*2+0] = index(ii,jj,nI);
				cells.faces[faceno*2+1] = index(ii,jj+1,nI);
				cells.left[faceno] = static_cast<int>(index(ii-1,jj,nI-1));
				cells.right[faceno] = static_cast<int>(index(ii,jj,nI-1));
				faceno++;
			}
	}
}


void Write_mesh(SIM const& svar, MESH const& cells)
{
	cout << "Writing mesh..." << endl;
	string meshOut = svar.infolder;
	meshOut.append(svar.meshname);

	if(DIM == 2)
	{
		meshOut.append(".edges");
	}
	else
	{
		meshOut.append(".faces");
	}

	#ifdef DEBUG
	dbout << "Entering Write_Face_Data..." << endl;
	dbout << "Output file: " << meshOut << endl;
	#endif

	NcFile fout(meshOut, NcFile::replace);

	/*Dimensions needed*/
	NcDim nElems = fout.addDim("no_of_elements",cells.nElems);
	NcDim nFaces, ppFace, nWall, nFar;

	if(DIM == 2)
	{
		nFaces = fout.addDim("no_of_edges",cells.nFaces);
		ppFace = fout.addDim("points_per_edge",DIM);
		nWall = fout.addDim("no_of_wall_edges",cells.nWall);
		nFar = fout.addDim("no_of_farfield_edges",cells.nFar);
	}
	else
	{
		nFaces = fout.addDim("no_of_faces",cells.nFaces);
		ppFace = fout.addDim("points_per_face",DIM);
		nWall = fout.addDim("no_of_wall_faces",cells.nWall);
		nFar = fout.addDim("no_of_farfield_faces",cells.nFar);
	}
	NcDim nPoint = fout.addDim("no_of_points",cells.nVerts);

	/*Define the faces*/
	vector<NcDim> faceVar;
	faceVar.emplace_back(nFaces);
	faceVar.emplace_back(ppFace);
	NcVar elemFaces, leftElems, rightElems;
	if (DIM == 2)
	{
		elemFaces = fout.addVar("points_of_element_edges",ncInt,faceVar);
		leftElems = fout.addVar("left_element_of_edges",ncInt,nFaces);
		rightElems = fout.addVar("right_element_of_edges",ncInt,nFaces);
	}
	else
	{
		elemFaces = fout.addVar("points_of_element_faces",ncInt,faceVar);
		leftElems = fout.addVar("left_element_of_faces",ncInt,nFaces);
		rightElems = fout.addVar("right_element_of_faces",ncInt,nFaces);
		// nElemFaces = fout.addVar("faces_per_element",ncInt,nElems);
	}
	/*Define the points*/
	NcVar vertsX = fout.addVar("points_xc",ncDouble,nPoint);
	NcVar vertsY = fout.addVar("points_yc",ncDouble,nPoint);
	

	/*Put faces into the file*/
	elemFaces.putVar(cells.faces);
	leftElems.putVar(cells.left);
	rightElems.putVar(cells.right);

	vertsX.putVar(cells.x);
	vertsY.putVar(cells.y);
	if(DIM == 3)
	{
		NcVar vertsZ = fout.addVar("points_zc",ncDouble,nPoint);
		vertsZ.putVar(cells.z);
	}

	cout << "Mesh written." << endl;
}

void Write_Solution(SIM const& svar, MESH const& cells)
{
	cout << "Writing solution..." << endl;
	string solOut = svar.infolder;
	solOut.append(svar.solname);


	NcFile sout(solOut, NcFile::replace);

	NcDim nVert = sout.addDim("no_of_points",cells.nVerts);

	NcVar dens = sout.addVar("density",ncDouble,nVert);
	NcVar velX = sout.addVar("x_velocity",ncDouble,nVert);
	NcVar velY = sout.addVar("y_velocity",ncDouble,nVert);
	NcVar press = sout.addVar("pressure",ncDouble,nVert);

	dens.putVar(cells.rho);
	velX.putVar(cells.vx);
	velY.putVar(cells.vy);
	if(DIM == 3)
	{
		NcVar velZ = sout.addVar("z_velocity",ncDouble,nVert);
		velZ.putVar(cells.vz);
	}
	press.putVar(cells.p);

	cout << "Solution written." << endl;
}

void Write_Tecplot(SIM const& svar, MESH const& cells)
{
	string meshOut = svar.infolder;
	meshOut.append(svar.meshname);
	meshOut.append(".dat");

	#ifdef DEBUG
	dbout << "Entering Write_ASCII_Face_Data..." << endl;
	cout << "Attempting write output file." << endl;
	cout << "File: " << meshOut << endl;
	#endif
	std::ofstream fout(meshOut,std::ios::out);
	if(!fout.is_open())
	{
		cout << "Couldn't open the output file." << endl;
		exit(-1);
	}

	uint TotalNumFaceNodes = 0;
	if(DIM == 3)
		TotalNumFaceNodes = cells.nFaces*3;



	if(DIM == 3)
	{
		fout << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
		fout << "ZONE T=\"FEPOLYHEDRON Test\"" << endl;
		fout << "ZONETYPE=FEPOLYHEDRON" << endl;
		fout << "NODES=" << cells.nVerts << " ELEMENTS=" << cells.nElems << " FACES=" << cells.nFaces << endl;
		fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
		fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;
	}
	else
	{
		fout << "VARIABLES= \"X\" \"Z\" " << endl;
		fout << "ZONE T=\"FEPOLYGON Test\"" << endl;
		fout << "ZONETYPE=FEPOLYGON" << endl;
		fout << "NODES=" << cells.nVerts << " ELEMENTS=" << cells.nElems << " FACES=" << cells.nFaces << endl;
		// fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
		fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;
	}

	uint w = 15;
	uint preci = 6;
	fout << std::left << std::scientific << std::setprecision(preci);

	uint newl = 0;
	fout << std::setw(1);
	for(uint ii = 0; ii < cells.nVerts; ++ii)
	{
		fout << std::setw(w) << cells.x[ii];
		newl++;

		if(newl>4)
		{
			fout << endl;
			fout << " ";
			newl=0;
		}
	}
	fout << endl;

	for(uint ii = 0; ii < cells.nVerts; ++ii)
	{
		fout << std::setw(w) << cells.y[ii];
		newl++;

		if(newl>4)
		{
			fout << endl;
			fout << " ";
			newl=0;
		}
	}
	fout << endl;

	if(DIM == 3)
	{
		for(uint ii = 0; ii < cells.nVerts; ++ii)
		{
			fout << std::setw(w) << cells.z[ii];
			newl++;

			if(newl>4)
			{
				fout << endl;
				fout << " ";
				newl=0;
			}
			
		}
		fout << endl;
	}

	fout << std::left << std::fixed;
	w = 9;

	size_t nPpF = 2;

	if(DIM == 3)
	{
		nPpF = 3;
	
		fout << "#node count per face" << endl;
		newl = 0;
		for (uint ii = 0; ii < cells.nFaces; ++ii)
		{
			fout << std::setw(w) << nPpF;
			newl++;

			if(newl>4)
			{
				fout << endl;
				newl=0;
			}
		}
		fout << endl;
	}

	/*Write the face data*/
	fout << "#face nodes" << endl;
	size_t facecount = 0;
	for (uint ii = 0; ii < cells.nFaces; ++ii)
	{
		for(size_t jj = 0; jj < nPpF; jj++)
		{	/*Write face vertex indexes*/
			fout << std::setw(w) << cells.faces[facecount*nPpF+jj]+1;
		}
		facecount++;
		fout << endl;
	}

	/*Write face left and right*/
	newl = 0;
	fout << "#left elements" << endl;
	for (uint ii = 0; ii < cells.nFaces; ++ii)
	{
		fout << std::setw(w) << cells.left[ii]+1 ;
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
	for (uint ii = 0; ii < cells.nFaces; ++ii)
	{
		if(cells.right[ii] < 0)
			fout<< std::setw(w) << 0 ;
		else
			fout << std::setw(w) << cells.right[ii]+1;


		newl++;

		if(newl>4)
		{
			fout << endl;
			newl=0;
		}
	}

	fout.close();

	#ifdef DEBUG
	dbout << "Exiting Write_Tecplot..." << endl;
	#endif

}

int main(int argc, char* argv[])
{
	SIM svar;
	MESH cells;

	GetInput(argc,argv,svar);

	Create_Mesh(svar,cells);

	// Write_Tecplot(svar,cells);

	Write_mesh(svar,cells);

	Write_Solution(svar,cells);


	return 0;
}