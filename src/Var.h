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

#if SIMDIM == 3
	#include "VLM.h"
#endif

#include <limits>
#include <vector>
#include <fstream>
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
using std::endl;

/* Define data type. */
typedef double ldouble;
typedef unsigned int uint;

/*Get machine bit precision for Simulation of Simplicity*/
#ifndef MEPSILON
#define MEPSILON std::numeric_limits<ldouble>::epsilon() /*For  float, power is -24*/
#endif

#ifndef MERROR
#define MERROR (7*MEPSILON + 56*MEPSILON*MEPSILON)
#endif

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<ldouble,SIMDIM,1> StateVecD;
typedef Eigen::Matrix<int,SIMDIM,1> StateVecI;
typedef Eigen::Matrix<ldouble,SIMDIM,SIMDIM> RotMat;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<ldouble, SIMDIM+1,1> DensVecD;
typedef Eigen::Matrix<ldouble, SIMDIM+1, SIMDIM+1> DensMatD;

/*Define particle type indexes*/
#define BOUND 0
#define START 1
#define PIPE 2
#define FREE 3
#define GHOST 4
#define PISTON 5

/*Simulation parameters*/
typedef struct SIM {
	StateVecI xyPART; 				/*Starting sim particles in x and y box*/
	size_t simPts,bndPts,totPts;	    /*Simulation particles, Boundary particles, total particles*/
	size_t finPts;					/*How many points there will be at simulation end*/
	size_t psnPts;					/*Piston Points*/
	size_t nrefresh; 					/*Crossflow refresh particle number*/
	uint nmax;                      /*Max number of particles*/
	uint outframe;	                /*Terminal output frame interval*/
	uint addcount;					/*Current Number of add-particle calls*/
	double Pstep,Bstep, dx, clear;	/*Initial spacings for particles and boundary*/
	uint nfull;						/*Full neighbour list amount*/
	StateVecD Box;					/*Box dimensions*/
	StateVecD Start; 				/*Sim box bottom left coordinate*/
	RotMat Rotate;				    /*Starting rotation matrix*/ 
	RotMat Transp;
	Eigen::Vector2d Jet;			/*Jet properties*/
	uint subits;                    /*Max number of sub-iterations*/
	uint Nframe; 			        /*Max number of frames to output*/
	uint frame;						/*Current frame number*/
	double dt, t, framet;			/*Timestep, Simulation + frame times*/
	double beta,gamma;				/*Newmark-Beta Parameters*/
	double maxmu;                   /*Maximum viscosity component (CFL)*/
	int Bcase, Bclosed, ghost;		/*What boundary shape to take*/
	uint outtype;                   /*ASCII or binary output*/
	uint outform, boutform, gout;   /*Output type. Fluid properties or Research.*/
	uint framecount;
	vector<size_t> back;            /*Particles at the back of the pipe*/

	std::string infolder, outfolder;
	std::string meshfile, bmapfile, solfile;
	ldouble scale;			
	#if SIMDIM == 3
		VLM vortex;
	#endif
	StateVecD Force;					/*Total Force*/

	/*Post Processessing settings*/
	uint afterSim;
	uint sliceOrSet;
	StateVecD planeOrigin, planeNorm; /*Plane conditions if in 3D*/
	ldouble cellSize;				/*Size of each cell to make*/
	StateVecD maxC, minC;			/*Max and min coords to make grid*/
	uint numVars, wrongDim;
	ldouble postRadius;
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	ldouble H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
	ldouble rho0, rhoJ; 					/*Resting Fluid density*/
	ldouble pPress, gasPress;		/*Starting pressure in pipe*/
	ldouble gasDynamic, gasVel;     /*Reference gas velocity*/
	ldouble Simmass, Boundmass;		/*Particle and boundary masses*/
	ldouble correc;					/*Smoothing Kernel Correction*/
	ldouble alpha,Cs,mu;		    /*}*/
	ldouble sig;					/* Fluid properties*/
	ldouble gam, B; 				/*}*/
	ldouble mug;					/* Gas Properties*/
	ldouble rhog;
	ldouble T;						/*Temperature*/
	ldouble Rgas;                   /*Specific gas constant*/
	ldouble contangb;				/*Boundary contact angle*/
	ldouble resVel;					/*Reservoir velocity*/
	//double front, height, height0;		/*Dam Break validation parameters*/
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
		ldouble L;							/*Gissler Parameters*/
		ldouble td;							/* }*/
		ldouble omega;						/* }*/
		ldouble tmax;						/* }*/
		ldouble ycoef;						/* }*/
		ldouble woccl;						/* }*/
		ldouble Cf, Ck, Cd, Cb, Cdef;		/* }*/
		ldouble nfull;						/* }*/

		int acase;	                        /*Aerodynamic force case*/
		StateVecD vJet, vInf;               /*Jet + Freestream velocity*/
		ldouble Acorrect;					/*Correction factor for aero force*/
		ldouble a;                          /*Case 3 tuning parameters*/
		ldouble b;                          /*}*/
		ldouble h1;                         /*}*/
		ldouble h2;                         /*}*/
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
		cellCp = std::vector<double>(nC);
		cellP = std::vector<double>(nC);
		cellRho = std::vector<double>(nC);
	}
	
	/*Zone info*/
	std::string zone;
	uint numPoint, numElem;
	ldouble scale;

	/*Point based data*/
	std::vector<StateVecD> verts;
	std::vector<StateVecD> pVel;
	std::vector<double> pointCp;
	std::vector<double> pointP;
	std::vector<double> pointRho;

	/*Face based data*/
	vector<vector<size_t>> faces;
	vector<std::pair<int,int>> leftright;

	/*Cell based data*/
	vector<vector<size_t>> elems;
	vector<vector<StateVecD>> cVerts;
	vector<StateVecD> cCentre;
	vector<vector<size_t>> cFaces;
	vector<vector<size_t>> cNeighb;

	/*Solution vectors*/
	vector<StateVecD> cVel;
	vector<double> cellCp;
	vector<double> cellP;
	vector<double> cellRho;
}MESH;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(const StateVecD& X, const StateVecD& Vi, const ldouble Rhoi, const ldouble Mi, 
			const ldouble press, const int bound, const uint pID)
		{
			xi = X;	v = Vi; 
			rho = Rhoi; m = Mi; b = bound;
			partID = pID;
			f = StateVecD::Zero();
			Sf = StateVecD::Zero(); 
			Af = StateVecD::Zero();
			normal = StateVecD::Zero();
			p = press;
			Rrho = 0.0;
			theta = 0.0;
			cellV = StateVecD::Zero();
			cellID = 0;
			cellP = 0.0;
			cellRho = 0.0;
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

		ldouble operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		// void erase_list(void)
		// {
		// 	list.erase(list.begin(),list.end());
		// }
		// std::vector<uint> list; /*Neighbour list*/

		StateVecD xi, v, Sf, Af, normal, f;
		ldouble Rrho, rho, p, m, theta;
		uint b; //What state is a particle. Boundary, forced particle or unforced
		StateVecD cellV;
		size_t partID, cellID;
		ldouble cellP, cellRho;
}Particle;

typedef class Part {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Part(const Particle& pi)
		{
			partID = pi.partID;
			xi = pi.xi; v = pi.v; 
			rho = pi.rho; 
			p = pi.p;
			m = pi.m;
			b = pi.b;
			Sf = pi.Sf;
			normal = pi.normal;
			cellV = pi.cellV;
			cellID = pi.cellID;
			cellP = pi.cellP;
			cellRho = pi.cellRho;
		}

		Part(const StateVecD xin, const StateVecD vin, const ldouble pin, 
				const ldouble rhoin, const ldouble min, const uint bin, const uint pID)
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

		ldouble operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		void operator=(const Particle pi)
		{
			partID = pi.partID;
			xi = pi.xi; v = pi.v; 
			rho = pi.rho; 
			p = pi.p;
			m = pi.m;
			b = pi.b;
			Sf = pi.Sf;
			cellV = pi.cellV;
			cellID = pi.cellID;
			cellP = pi.cellP;
			cellRho = pi.cellRho;
		}



		// void erase_list(void)
		// {
		// 	list.erase(list.begin(),list.end());
		// }
		// std::vector<uint> list; /*Neighbour list*/

		StateVecD xi, v, Sf, normal;
		ldouble rho, p, m;
		uint b; //What state is a particle. Boundary, forced particle or unforced, or air
		StateVecD cellV;
		uint partID, cellID;
		ldouble cellP, cellRho;
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
typedef KDTreeVectorOfVectorsAdaptor<State,ldouble,SIMDIM,nanoflann::metric_L2_Simple,size_t> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,ldouble,SIMDIM,nanoflann::metric_L2_Simple,size_t> Vec_Tree;
#endif /* VAR_H */
