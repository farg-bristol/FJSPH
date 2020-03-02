/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VAR_H
#define VAR_H

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 2
#endif

// #ifndef NTHREADS
// #define NTHREADS 4
// #endif



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
	std::ofstream dbout("WCSPH.log",std::ios::out);
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
#ifndef FOD
#define FOD 1 /*0 = float, 1 = double*/
#endif

typedef unsigned int uint;

#if FOD == 1
typedef double real;
#else
typedef float real;
#endif

/*Get machine bit precision for Simulation of Simplicity*/
#ifndef MEPSILON
#define MEPSILON std::numeric_limits<real>::epsilon() /*For  float, power is -24*/
#endif

#ifndef MERROR
#define MERROR (7*MEPSILON + 56*MEPSILON*MEPSILON)
#endif

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real,SIMDIM,1> StateVecD;
typedef Eigen::Matrix<int,SIMDIM,1> StateVecI;
typedef Eigen::Matrix<real,SIMDIM,SIMDIM> RotMat;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<real, SIMDIM+1,1> DensVecD;
typedef Eigen::Matrix<real, SIMDIM+1, SIMDIM+1> DensMatD;

#if SIMDIM == 3
	#include "VLM.h"
#endif

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
	StateVecI xyPART; 				/*Starting sim particles in x and y box*/
	size_t simPts,bndPts,totPts;	    /*Simulation particles, Boundary particles, total particles*/
	size_t finPts;					/*How many points there will be at simulation end*/
	size_t psnPts;					/*Piston Points*/
	size_t nrefresh; 					/*Crossflow refresh particle number*/
	size_t delNum;                  /*Number of deleted particles*/
	size_t intNum;                  /*Number of internal particles*/
	uint nmax;                      /*Max number of particles*/
	uint outframe;	                /*Terminal output frame interval*/
	uint addcount;					/*Current Number of add-particle calls*/
	real Pstep,Bstep, dx, clear;	/*Initial spacings for particles and boundary*/
	uint nfull;						/*Full neighbour list amount*/
	StateVecD Box;					/*Box dimensions*/
	StateVecD Start; 				/*Sim box bottom left coordinate*/
	RotMat Rotate;				    /*Starting rotation matrix*/ 
	RotMat Transp;
	Eigen::Vector2d Jet;			/*Jet properties*/
	uint subits;                    /*Max number of sub-iterations*/
	uint Nframe; 			        /*Max number of frames to output*/
	uint frame;						/*Current frame number*/
	real dt, t, framet;			/*Timestep, Simulation + frame times*/
	real beta,gamma;				/*Newmark-Beta Parameters*/
	real maxmu;                   /*Maximum viscosity component (CFL)*/
	int Bcase, Bclosed, ghost;		/*What boundary shape to take*/
	uint outtype;                   /*ASCII or binary output*/
	uint outform, boutform, gout;   /*Output type. Fluid properties or Research.*/
	uint framecount;
	vector<size_t> back;            /*Particles at the back of the pipe*/

	std::string infolder, outfolder;
	std::string meshfile, bmapfile, solfile;
	real scale;			
	#if SIMDIM == 3
		VLM vortex;
	#endif
	StateVecD Force;					/*Total Force*/

	uint restart;

	/*Post Processessing settings*/
	uint afterSim;
	uint sliceOrSet;
	StateVecD planeOrigin, planeNorm; /*Plane conditions if in 3D*/
	real cellSize;				/*Size of each cell to make*/
	StateVecD maxC, minC;			/*Max and min coords to make grid*/
	uint numVars, wrongDim;
	real postRadius;
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	real Hfac;
	real H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
	real rho0, rhoJ; 			/*Resting Fluid density*/
	real pPress, gasPress;		/*Starting pressure in pipe*/
	real gasDynamic, gasVel;    /*Reference gas velocity*/
	real simM, bndM;			/*Particle and boundary masses*/
	real gasM;					/*A gas particle mass*/
	real correc;				/*Smoothing Kernel Correction*/
	real alpha,Cs,mu;		    /*}*/
	real sig;					/* Fluid properties*/
	real gam, B; 				/*}*/
	real mug;					/* Gas Properties*/
	real rhog;
	real T;						/*Temperature*/
	real Rgas;                  /*Specific gas constant*/
	real contangb;				/*Boundary contact angle*/
	real resVel;				/*Reservoir velocity*/
	//real front, height, height0;		/*Dam Break validation parameters*/
}FLUID;

/*Aerodynamic Properties*/
typedef class AERO
{
	public:
		AERO()
		{
			Cf = 1.0/3.0;
			Ck = 8;
			Cd = 5;
			Cb = 0.5;
		}
		real L;							/*Gissler Parameters*/
		real td;							/* }*/
		real omega;						/* }*/
		real tmax;						/* }*/
		real ycoef;						/* }*/
		real woccl;						/* }*/
		real Cf, Ck, Cd, Cb, Cdef;		/* }*/
		real nfull;						/* }*/

		int acase;	                        /*Aerodynamic force case*/
		StateVecD vJet, vInf;               /*Jet + Freestream velocity*/
		real Acorrect;					/*Correction factor for aero force*/
		real a;                          /*Case 3 tuning parameters*/
		real b;                          /*}*/
		real h1;                         /*}*/
		real h2;                         /*}*/
}AERO;

typedef struct MESH
{
	/*Standard contructor*/
	MESH(){}

	// void Create_Tree()
	// {

	// }

	void SetCells()
	{
		uint nC = elems.size();
		cVel = std::vector<StateVecD>(nC);		
		cellCp = std::vector<real>(nC);
		cellP = std::vector<real>(nC);
		cellRho = std::vector<real>(nC);
	}
	
	/*Zone info*/
	std::string zone;
	uint numPoint, numElem;
	real scale;

	/*Point based data*/
	std::vector<StateVecD> verts;
	std::vector<StateVecD> pVel;
	std::vector<real> pointCp;
	std::vector<real> pointP;
	std::vector<real> pointRho;

	/*Face based data*/
	vector<vector<size_t>> faces;
	vector<std::pair<int,int>> leftright;

	/*Cell based data*/
	vector<vector<size_t>> elems;
	vector<vector<StateVecD>> cVerts;
	vector<StateVecD> cCentre;
	vector<vector<size_t>> cFaces;
	// vector<vector<size_t>> cNeighb;

	/*Boundary data*/
	vector<size_t> bIndex;
	vector<StateVecD> bVerts;

	/*Solution vectors*/
	vector<StateVecD> cVel;
	vector<real> cellCp;
	vector<real> cellP;
	vector<real> cellRho;
}MESH;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(const StateVecD& X, const StateVecD& Vi, const real Rhoi, const real Mi, 
			const real press, const int bound, const uint pID)
		{
			xi = X;	v = Vi; 
			rho = Rhoi; m = Mi; b = bound;
			partID = pID;
			f = StateVecD::Zero();
			Sf = StateVecD::Zero(); 
			Af = StateVecD::Zero();
			normal = StateVecD::Zero();
			vPert = StateVecD::Zero();
			p = press;
			Rrho = 0.0;
			theta = 0.0;
			cellV = StateVecD::Zero();
			cellID = 0;
			faceID = 0;
			cellP = 0.0;
			cellRho = 0.0;
			internal = 0;

		}

		/*To add particles dynamically for boundary layer*/
		Particle(const StateVecD& X,  const Particle& pj, const size_t pID)
		{
			xi = X;	v = pj.v; 
			rho = pj.rho; m = pj.m; b = pj.b;
			partID = pID;
			f = StateVecD::Zero();
			Sf = StateVecD::Zero(); 
			Af = StateVecD::Zero();
			normal = StateVecD::Zero();
			p = pj.p;
			Rrho = 0.0;
			theta = 0.0;
			cellV = StateVecD::Zero();
			cellID = 0;
			cellP = 0.0;
			cellRho = 0.0;
		}

		Particle(){};

		int size() const
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		real operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		StateVecD xi, v, f, Sf, Af, normal, vPert;
		real Rrho, rho, p, m, theta;
		uint b; //What state is a particle. Boundary, forced particle or unforced
		StateVecD cellV;
		size_t partID, cellID, faceID;
		real cellP, cellRho;
		uint internal; 
		StateVecD bNorm;
		real y;
}Particle;

typedef class Part {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Part(const Particle& pi)
		{
			xi = pi.xi; v = pi.v; Sf = pi.Sf;
			normal = pi.normal;	vPert = pi.vPert;
			rho = pi.rho; p = pi.p;	m = pi.m; b = pi.b;
			cellV = pi.cellV;
			partID = pi.partID;
			cellID = pi.cellID;
			cellP = pi.cellP;
			cellRho = pi.cellRho;
			internal = pi.internal;
			bNorm = pi.bNorm;
			y = pi.y;
		}

		Part(const StateVecD& xin, const StateVecD& vin, const real pin, 
				const real rhoin, const real min, const uint bin, const uint pID)
		{
			xi = xin;
			v = vin;
			p = pin;
			rho = rhoin;
			m = min;
			b = bin;
			partID = pID;
			cellV = StateVecD::Zero();
			cellP = 0.0;
			cellRho = 0.0;
		}

		int size() const
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		real operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		void operator=(const Particle& pi)
		{
			xi = pi.xi; v = pi.v; Sf = pi.Sf; vPert = pi.vPert;
			rho = pi.rho;	p = pi.p;	m = pi.m;	b = pi.b;
			cellV = pi.cellV;
			partID = pi.partID;
			cellID = pi.cellID;
			cellP = pi.cellP;
			cellRho = pi.cellRho;
			internal = pi.internal;
			bNorm = pi.bNorm;
			y = pi.y;
		}

		StateVecD xi, v, Sf, normal, vPert;
		real rho, p, m;
		uint b; //What state is a particle. Boundary, forced particle or unforced, or air
		StateVecD cellV;
		size_t partID, cellID/*, faceID*/;
		real cellP, cellRho;
		uint internal;
		StateVecD bNorm;
		real y;
}Part;

Particle PartToParticle(Part& pj)
{
	Particle pi;
	pi.partID = pj.partID;
	pi.xi = pj.xi; pi.v = pj.v; 
	pi.rho = pj.rho; 
	pi.p = pj.p;
	pi.m = pj.m;
	pi.b = pj.b;
	pi.Sf = pj.Sf;
	pi.cellV = pj.cellV;
	pi.cellID = pj.cellID;
	pi.cellP = pj.cellP;
	pi.cellRho = pj.cellRho;
	return pi;
}

typedef std::vector<Particle> State;

/* Neighbour search tree containers */
typedef std::vector<std::vector<size_t>> outl;
typedef KDTreeVectorOfVectorsAdaptor<State,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

typedef struct KDTREE
{
	KDTREE(const State& pnp1, const MESH& cells): NP1(SIMDIM,pnp1,20), 
	CELL(SIMDIM,cells.cCentre,20), BOUNDARY(SIMDIM,cells.bVerts,20) {}
	Sim_Tree NP1;
	Vec_Tree CELL;
	Vec_Tree BOUNDARY;	
}KDTREE;

#endif /* VAR_H */
