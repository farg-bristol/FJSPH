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
#include <string.h>
#include <sstream>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include "NanoFLANN/nanoflann.hpp"
#include "NanoFLANN/utils.h"
#include "NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"
#include <omp.h>

#ifdef DEBUG
	/*Open debug file to write to*/	
	std::ofstream dbout("WCSPH.log",std::ios::out);
#endif

// std::ofstream pertLog("cellPert.log",std::ios::out);

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

#if FOD == 1
typedef double real;
#else
typedef float real;
#endif

typedef unsigned int uint;

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
typedef Eigen::Matrix<real,SIMDIM,SIMDIM> StateMatD;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<real, SIMDIM+1,1> StateP1VecD;
typedef Eigen::Matrix<real, SIMDIM+1, SIMDIM+1> StateP1MatD;

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
typedef struct PState{
	PState()
	{
		BOUND_ = 0;
	    PISTON_ = 1;
		BACK_ = 2; 
		START_ = 3;
		PIPE_ = 4;
		FREE_ = 5;
		GHOST_ = 6;
	}

	size_t BOUND_ ,
    PISTON_,
	BACK_,            
	START_,
	PIPE_,
	FREE_,
	GHOST_;
} PState;

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

PState PartState;

/*Simulation parameters*/
typedef struct SIM {
	SIM()
	{
#if SIMDIM == 3
		nfull = 257;
#else
		nfull = 48;
#endif
	}

	StateVecI xyPART; 				/*Starting sim particles in x and y box*/
	size_t simPts,bndPts,totPts;	/*Simulation particles, Boundary particles, total particles*/
	size_t finPts;					/*How many points there will be at simulation end*/
	size_t psnPts;					/*Piston Points*/
	size_t nrefresh; 				/*Crossflow refresh particle number*/
	size_t delNum;                  /*Number of deleted particles*/
	size_t intNum;                  /*Number of internal particles*/
	uint nmax;                      /*Max number of particles*/
	uint outframe;	                /*Terminal output frame interval*/
	uint outtype;                   /*ASCII or binary output*/
	uint outform, boutform, gout;   /*Output type. Fluid properties or Research.*/
	uint addcount;					/*Current Number of add-particle calls*/
	real Pstep,Bstep, dx, clear;	/*Initial spacings for particles and boundary*/
	real diam;                      /*Droplet diameter*/
	uint nrad;                      /*Points along the radius of the jet/droplet*/
	uint nfull;						/*Full neighbour list amount*/
	StateVecD Box;					/*Box dimensions*/
	StateVecD Start; 				/*Sim box bottom left coordinate*/
	StateVecD Angle;				/*Rotations in degrees*/
	StateMatD Rotate;				/*Starting rotation matrix*/ 
	StateMatD Transp;               /*Transpose of rotation matrix*/
	Eigen::Vector2d Jet;			/*Jet properties*/
	uint subits;                    /*Max number of sub-iterations*/
	uint Nframe; 			        /*Max number of frames to output*/
	uint frame;						/*Current frame number*/
	real dt, framet;			    /*Timestep, frame times*/
	double t;                       /*Simulation time*/
	real beta,gamma;				/*Newmark-Beta Parameters*/
	real maxmu;                     /*Maximum viscosity component (CFL)*/
	int Bcase, Bclosed, ghost;		/*What initial shape to take*/
	int Asource;                     /*Source of aerodynamic solution*/
	
	uint framecount;                /*How many frames have been output*/
	vector<size_t> back;            /*Particles at the back of the pipe*/
	uint iter;                      /*Current iteration number to renormalise*/

	std::string infolder, outfolder, outdir;
	std::string meshfile, bmapfile, solfile;
	void* boundFile; /*TECIO file handles*/
	void* fuelFile;
	void* ghostFile;
	vector<int32_t> varTypes;

	real scale;			
	#if SIMDIM == 3
		VLM vortex;
	#endif
	StateVecD Force, AForce;					/*Total Force*/
	real mass;
	real tMom, aMom;


	uint restart;

	/*Post Processessing settings*/
	uint afterSim;
	uint sliceOrSet;
	StateVecD planeOrigin, planeNorm; /*Plane conditions if in 3D*/
	real cellSize;				      /*Size of each cell to make*/
	StateVecD maxC, minC;			  /*Max and min coords to make grid*/
	uint numVars, wrongDim;
	real postRadius;
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	real H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
	real Wdx;                   /*Kernel value at the initial particle spacing distance.*/
	real rho0, rhoJ; 			/*Resting Fluid density*/
	real pPress;		/*Starting pressure in pipe*/
	
	real simM, bndM;			/*Particle and boundary masses*/
	real correc;				/*Smoothing Kernel Correction*/
	real alpha,Cs,mu,nu;		/*}*/
	real sig;					/* Fluid properties*/
	real gam, B; 				/*}*/
	real artMu;					/*Artificial viscosity*/
	// real mug;
	// real rhog;
	
	real contangb;				/*Boundary contact angle*/
	real resVel;				/*Reservoir velocity*/
	//real front, height, height0;		/*Dam Break validation parameters*/

	/*delta SPH terms*/
	real delta; /*delta-SPH contribution*/
	real dCont; /*delta-SPH continuity constant term, since density and speed of sound are constant*/
	real dMom;  /*delta-SPH momentum constant term (dCont * rho0)*/

	real maxU; /*Maximum particle velocity shift*/
}FLUID;

/*Aerodynamic Properties*/
typedef struct AERO
{
	AERO()
	{
		Cf = 1.0/3.0;
		Ck = 8;
		Cd = 5;
		Cb = 0.5;
#if SIMDIM == 3		
		nfull = 1.713333e+02;
#else
		nfull = 28;
#endif
	}

		/*Gissler Parameters*/
	real L;
	real td;
	real omega;
	real tmax;
	real ycoef;
	real Cf, Ck, Cd, Cb, Cdef;
	real nfull;

	real pVol;                     /*Volume of a particle*/
	real aPlate;                   /*Area of a plate*/

		/* Gas Properties*/
	real qInf, vRef, pRef;         /*Reference gas values*/
	real rhog, mug;
	real gasM;					   /*A gas particle mass*/
	real T;						   /*Temperature*/
	real Rgas;                     /*Specific gas constant*/

	int acase;	                   /*Aerodynamic force case*/
	StateVecD vJet, vInf;          /*Jet + Freestream velocity*/
	real vJetMag;
	
	real Acorrect;				   /*Correction factor for aero force*/
	real a;                        /*Case 3 tuning parameters*/
	real b;                        /*}*/
	real h1;                       /*}*/
	real h2;                       /*}*/
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
		cCp = std::vector<real>(nC);
		cP = std::vector<real>(nC);
		cRho = std::vector<real>(nC);
		cMass = vector<real>(nC);
		fNum = vector<size_t>(nC,0);
		fMass = vector<real>(nC,0);
		cPertn = vector<StateVecD>(nC,StateVecD::Zero());
		cPertnp1 = vector<StateVecD>(nC,StateVecD::Zero());
		vFn = vector<StateVecD>(nC,StateVecD::Zero());
		vFnp1 = vector<StateVecD>(nC,StateVecD::Zero());
	}
	

	size_t size()
	{
		return elems.size();
	}

	/*Zone info*/
	std::string zone;
	uint numPoint, numElem;
	real scale;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Face based data*/
	vector<vector<size_t>> faces;
	vector<std::pair<int,int>> leftright;

	/*Cell based data*/
	vector<vector<size_t>> elems;
	vector<StateVecD> cCentre;
	vector<vector<size_t>> cFaces;

	/*Boundary data*/
	vector<size_t> bIndex;
	vector<StateVecD> bVerts;

	/*Solution vectors*/
	vector<StateVecD> cVel;
	vector<real> cCp;
	vector<real> cP;
	vector<real> cMass;
	vector<real> cRho;
	
	// Cell information for the momentum balance
	vector<StateVecD> cPertn;
	vector<StateVecD> cPertnp1;
	vector<real> cVol;


	vector<size_t> fNum; // number of fuel particles in cell
	vector<real> fMass;

	vector<StateVecD> vFn;
	vector<StateVecD> vFnp1;

}MESH;

/*Container for delta-plus SPH calculations*/
typedef class DELTAP {
	public: 
		DELTAP(size_t const& size)
		{
			L = vector<StateMatD>(size,StateMatD::Zero());
			gradRho = vector<StateVecD>(size,StateVecD::Zero());
			norm = vector<StateVecD>(size,StateVecD::Zero());
			avgV = vector<StateVecD>(size,StateVecD::Zero());
			
			lam = vector<real>(size,0.0);
			lam_nb = vector<real>(size,0.0);
			kernsum = vector<real>(size,0.0);
		}

		DELTAP(){}

		void realloc(size_t const& size)
		{
			if(L.size() != 0)
				clear();

			alloc(size);
		}

		void update(vector<StateMatD> const& L_, vector<StateVecD> const& gradRho_, 
			vector<StateVecD> const& norm_, vector<StateVecD> const& avgV_,
			vector<real> const& lam_, vector<real> const& lam_nb_, 
			vector<real> const& kernsum_)
		{
			L = L_; gradRho = gradRho_; norm = norm_; 
			avgV = avgV_; lam = lam_; lam_nb = lam_nb_;
			kernsum = kernsum_;
		}

		void clear()
		{
			L.clear(); gradRho.clear(); norm.clear(); 
			avgV.clear(); lam.clear(); lam_nb.clear();
			kernsum.clear();

		}

		vector<StateMatD> L;
		vector<StateVecD> gradRho;
		vector<StateVecD> norm;
		vector<StateVecD> avgV;

		vector<real> lam, lam_nb;
		vector<real> kernsum;


	private:

		void alloc(size_t const& size)
		{
			L = vector<StateMatD>(size);
			gradRho = vector<StateVecD>(size);
			norm = vector<StateVecD>(size);
			avgV = vector<StateVecD>(size);
			lam = vector<real>(size);
			lam_nb = vector<real>(size);
			kernsum = vector<real>(size);
		}


}DELTAP;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(StateVecD const& X, StateVecD const& Vi, real const Rhoi, real const Mi, 
			real const press, int const bound, uint const pID)
		{
			partID = pID; cellID = 0; faceID = 0;
			b = bound; surf = 0;

			xi = X;	v = Vi; f = StateVecD::Zero(); Af = StateVecD::Zero();
			Rrho = 0.0; rho = Rhoi; p = press; m = Mi; curve = 0.0;
			theta = 0.0; nNeigb = 0.0; s = 0.0; woccl = 0.0; pDist = 0.0;
						
			cellV = StateVecD::Zero();
			cellP = 0.0;
			internal = 0;

			Sf = StateVecD::Zero();
			normal = StateVecD::Zero();
			vPert = StateVecD::Zero(); 
		}

		/*To add particles dynamically for boundary layer*/
		Particle(StateVecD const& X, Particle const& pj, size_t const pID, int const bound)
		{
			partID = pID; cellID = 0; faceID = 0;
			b = bound; surf = 0;

			xi = X;	v = pj.v; f = StateVecD::Zero(); Af = StateVecD::Zero();
			Rrho = 0.0; rho = pj.rho; p = pj.p; m = pj.m; curve = 0.0;
			theta = 0.0; nNeigb = 0.0; s = 0.0; woccl = 0.0; pDist = 0.0;

			cellV = StateVecD::Zero();
			cellP = 0.0;
			internal = 0;

			Sf = StateVecD::Zero();
			normal = StateVecD::Zero();
			vPert = StateVecD::Zero();
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

		size_t partID, cellID, faceID;
		uint b; //What state is a particle. See PartState above for possible options
		uint surf; /*Is a particle a surface? 1 = yes, 0 = no*/
		StateVecD xi, v, f, Af;
		real Rrho, rho, p, m, curve, theta, nNeigb, s, woccl, pDist;
		StateVecD cellV;
		real cellP;
		uint internal; 

		/*Mesh surface repulsion */
		StateVecD bNorm;
		real y;
		
		StateVecD Sf, normal, vPert;		
}Particle;

typedef class Part {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Part(Particle const& pi)
		{
			xi = pi.xi; v = pi.v; Sf = pi.Sf;
			normal = pi.normal;	vPert = pi.vPert;
			rho = pi.rho; p = pi.p;	m = pi.m; curve = pi.curve;
			woccl = pi.woccl; pDist = pi.pDist;
			b = pi.b; surf = pi.surf;
			cellV = pi.cellV;
			partID = pi.partID;
			cellID = pi.cellID;
			cellP = pi.cellP;
			internal = pi.internal;
			bNorm = pi.bNorm;
			y = pi.y;
			nDist = 0.0;
		}

		Part(Particle const &pi, real const& dist)
		{
			xi = pi.xi;
			v = pi.v;
			Sf = pi.Sf;
			normal = pi.normal;
			vPert = pi.vPert;
			rho = pi.rho;
			p = pi.p;
			m = pi.m;
			curve = pi.curve;
			woccl = pi.woccl;
			pDist = pi.pDist;
			nDist = dist; // Neighbour distance (squared)
			b = pi.b;
			surf = pi.surf;
			cellV = pi.cellV;
			partID = pi.partID;
			cellID = pi.cellID;
			cellP = pi.cellP;
			internal = pi.internal;
			bNorm = pi.bNorm;
			y = pi.y;
		}

		Part(StateVecD const& xin, StateVecD const& vin, real const pin, 
				real const rhoin, real const min, uint const bin, uint const pID)
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
			surf = 0;
			internal = 0;

		}

		int size() const
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		real operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		void operator=(Particle const& pi)
		{
			xi = pi.xi; v = pi.v; Sf = pi.Sf; 
			normal = pi.normal; vPert = pi.vPert;
			rho = pi.rho;	p = pi.p; m = pi.m; curve = pi.curve;
			woccl = pi.woccl; pDist = pi.pDist; 
			b = pi.b; surf = pi.surf;
			cellV = pi.cellV;
			partID = pi.partID;
			cellID = pi.cellID;
			cellP = pi.cellP;
			internal = pi.internal;
			bNorm = pi.bNorm;
			y = pi.y;
		}

		StateVecD xi, v, Sf, normal, vPert;
		real rho, p, m, curve, woccl, pDist, nDist;
		uint b, surf; //What state is a particle.
		StateVecD cellV;
		size_t partID, cellID/*, faceID*/;
		real cellP;
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
	return pi;
}

typedef std::vector<Particle> State;

/* Neighbour search tree containers */
typedef std::vector<std::vector<std::pair<size_t,real>>> outl;
typedef std::vector<std::vector<size_t>> celll;
typedef KDTreeVectorOfVectorsAdaptor<State,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

typedef struct KDTREE
{
	KDTREE(State const& pnp1, MESH const& cells): NP1(SIMDIM,pnp1,20), 
	CELL(SIMDIM,cells.cCentre,20), BOUNDARY(SIMDIM,cells.bVerts,20) {}

	Sim_Tree NP1;
	Vec_Tree CELL;
	Vec_Tree BOUNDARY;	

}KDTREE;

#endif /* VAR_H */
