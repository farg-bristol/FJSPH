/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef POSTVAR_H
#define POSTVAR_H

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
#include "NanoFLANN/nanoflann.hpp"
#include "NanoFLANN/utils.h"
#include "NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"

#ifdef DEBUG
	/*Open debug file to write to*/	
	std::ofstream dbout("PostProcessing.log",std::ios::out);
#endif

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using std::vector;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string; 
using std::setw;

/* Define data type. */
typedef unsigned int uint;
typedef double real;

uint DIM = 2;

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real,Eigen::Dynamic,1,0,3,1> StateVecD;
typedef Eigen::Matrix<int,Eigen::Dynamic,1,0,3,1> StateVecI;
typedef Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic,0,3,3> RotMat;

#pragma omp declare reduction(+: std::vector<StateVecD> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
        	[](StateVecD lhs, StateVecD rhs){return lhs + rhs;})) \
                    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(+: std::vector<real> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
        	[](real lhs, real rhs){return lhs + rhs;})) \
                    initializer(omp_priv = omp_orig)                  

#pragma omp declare reduction(+:StateVecD : omp_out=omp_out+omp_in)\
                    initializer(omp_priv = omp_orig) 

/*Define particle type indexes*/
#define BOUND 0
#define PISTON 1
#define BACK 2            
#define START 3
#define PIPE 4
#define FREE 5
#define GHOST 6


/*Simulation parameters*/
typedef struct SIM {
	size_t simPts,bndPts,totPts;	    /*Simulation particles, Boundary particles, total particles*/
	uint Nframe; 			        /*Max number of frames to output*/
	uint frame;						/*Current frame number*/
	int Bcase;  		/*What boundary shape to take*/
	uint outtype;                   /*ASCII or binary output*/
	uint outform, boutform;         /*Output type. Fluid properties or Research.*/
	uint framecount;

	std::string infolder, outfolder;
	real scale;			

	/*Post Processessing settings*/
	uint sliceOrSet;
	StateVecD planeOrigin, planeNorm; /*Plane conditions if in 3D*/

	int maxCells;				/*Number of cells along longest axis*/
	real cellSize;					/*Size of each cell to make*/
	StateVecD maxC, minC;			/*Max and min coords to make grid*/
	uint numVars;

	real Hfac;
	real H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
} SIM;

// SPH Point data structure
struct Point {
	Point(){};

	Point(const uint size, const uint type)
	{
		dataType = type;
		verts = vector<StateVecD>(size);
		rho = vector<real>(size);
		m = vector<real>(size);

		if (type == 1 || type == 4)
		{
			vnorm = vector<real>(size);
			anorm = vector<real>(size);
			if(type == 4)
			{
				Af = vector<real>(size);
			}
		}
		else if (type == 2 || type == 5)
		{
			vel = vector<StateVecD>(size);
			acc = vector<StateVecD>(size);
		}
		else if(type == 3)
		{
			vel = vector<StateVecD>(size);
			acc = vector<StateVecD>(size);
			cellV = vector<StateVecD>(size);
			cellP = vector<real>(size);
			cellRho = vector<real>(size);
		}
	}

	void operator=(const Point& pi)
	{
		verts = pi.verts; vel = pi.vel; acc = pi.acc;
		vnorm = pi.vnorm; anorm = pi.anorm;
		rho = pi.rho; m = pi.m; Af = pi.Af;
		cellV = pi.cellV; cellP = pi.cellP; cellRho = pi.cellRho;
		time = pi.time;
	}

	void append(const Point& pi)
	{
		verts.insert(verts.end(),pi.verts.begin(),pi.verts.end());
		vel.insert(vel.end(),pi.vel.begin(),pi.vel.end());
		acc.insert(acc.end(),pi.acc.begin(),pi.acc.end());
		vnorm.insert(vnorm.end(),pi.vnorm.begin(),pi.vnorm.end());
		anorm.insert(anorm.end(),pi.anorm.begin(),pi.anorm.end());
		rho.insert(rho.end(),pi.rho.begin(),pi.rho.end());
		m.insert(m.end(),pi.m.begin(),pi.m.end());
		cellV.insert(cellV.end(),pi.cellV.begin(),pi.cellV.end());
		cellP.insert(cellP.end(),pi.cellP.begin(),pi.cellP.end());
		cellRho.insert(cellRho.end(),pi.cellRho.begin(),pi.cellRho.end());
		time = pi.time;
	}

	uint size()
	{
		return verts.size();
	}

	uint dataType;
	double time;
	vector<StateVecD> verts, vel, acc;
	vector<real> vnorm, anorm;
	vector<real> rho, m, Af;
	vector<StateVecD> cellV;
	vector<real> cellP, cellRho;
};

/* Neighbour search tree containers */
typedef std::vector<std::vector<size_t>> outl;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,real,-1,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

#endif /* POSTVAR_H */
