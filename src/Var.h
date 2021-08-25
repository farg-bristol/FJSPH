/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
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
#include <omp.h>

// Third party includes
#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"
#include "Third_Party/NanoFLANN/nanoflann.hpp"
#include "Third_Party/NanoFLANN/utils.h"
#include "Third_Party/NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"

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
		BUFFER_ = 2;
		BACK_ = 3; 
		PIPE_ = 4;
		FREE_ = 5;
		GHOST_ = 6;
	}

	size_t BOUND_ ,
    PISTON_,
	BUFFER_, 
	BACK_,
	PIPE_,
	FREE_,
	GHOST_;
} PState;

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

PState PartState;

real random(int const& interval)
{
	return real(rand() % interval - interval / 2) * MERROR;
}

typedef struct SIM {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	SIM()
	{
		/* Input files */
		CDForFOAM = 0;
		isBinary = 0; labelSize = 32; scalarSize = 64;
		offset_axis = 0; angle_alpha = 0;

		out_encoding = 1; /* default to binary */
		outform = 5; gout = 0; /* Default to restartable data */

		restart = 0;
		scale = 1.0;

		partID = 0;
		simPts = 0; bndPts = 0; totPts = 0;
		gstPts = 0;  psnPts = 0; delNum = 0; 
		intNum = 0;	nrefresh = 0; addcount = 0;
		finPts = 9999999; /* Default to basically no particle limit */
		
		Pstep = -1.0; dx = -1.0; Bstep = 0.7;

		bound_type = "(none)"; /* Default to no boundary */
		Scase = 0; Bcase = 0;
		Bclosed = 0; ghost = 0;

		/* Universal geometry parameters */
		sim_start.setConstant(-1);
		bound_start.setConstant(-1);
		Angle = StateVecD::Zero();
		Rotate = StateMatD::Zero();
		Transp = StateMatD::Zero();

		/* Jet geometry parameters */
		jet_diam = -1;
		jet_depth = -1;
		nrad = 0;
		
		/* Droplet geometry parameters */
		diam = -1.0;
		
		/* Rectangle */
		xyPART.setConstant(-1);
		sim_box.setConstant(-1);
		bound_box.setConstant(-1);
	
		/* Integration parameters */
		subits = 20;
		Nframe = -1;
		frame = 0;
		t = 0.0;
		dt = 2e-10;
		framet = -1;
		dt_max = 1;
		dt_min = 0.0;
		beta = 0.25; gamma = 0.5;
		maxmu = 0;

		/*Gravity Vector*/
		#if SIMDIM == 3
			grav = StateVecD(0.0,0.0,-9.81);
		#else
			grav = StateVecD(0.0,-9.81);
		#endif

		framecount = 0;
		Force = StateVecD::Zero(); AForce = StateVecD::Zero();					/*Total Force*/
		mass = 0;
		tMom = 0; aMom = 0;
	
		/* SPHPart tracking settings */
		eqOrder = 2;
		max_x_sph = 9999999;
		max_x = 9999999;
		nSuccess = 0; nFailed = 0;
		IPT_diam = 0.0; IPT_area = 0.0;
		cellsout = 0; streakout = 1; partout = 0;
	}

	/* File input parameters */
	uint CDForFOAM;
	std::string infile, output_prefix, outdir;
	std::string restart_prefix;

	/* OpenFOAM files */
	std::string foamdir, foamsol;
	int isBinary, labelSize, scalarSize;
	
	/* TAU files */
	std::string taumesh, taubmap, tausol;
	int offset_axis;
	real angle_alpha;
	vector<int> markers;
	vector<string> bnames;
	vector<int> bwrite; /* If surface wants data writing */

	/* TECIO File outputs */
	void* boundFile; /*TECIO file handles*/
	void* fuelFile;
	void* ghostFile;
	vector<int32_t> varTypes;

	ofstream surfacefile; /* Surface impact file */

	/* Output type */
	uint out_encoding;              /*ASCII or binary output*/
	uint outform, gout;   			/*Output type.*/

	uint restart;					/* If starting from existing solution */
	real scale;						/* Simulation scale */

	/* SPHPart counts */
	size_t partID;                   /* Track the particle ID separately */
	size_t totPts;					 /* Total count */
	size_t simPts;					 /* Fluid count */
	size_t bndPts;			 		 /* Boundary count */
	size_t gstPts;                   /* Ghost count */
	size_t finPts;					 /* End count */
	size_t psnPts;					 /* Piston count */
	size_t delNum;                   /* Number of deleted particles */
	size_t intNum;                   /* Number of internal particles */
	size_t nrefresh; 				 /* last add call particle number */
	uint addcount;					 /* Current Number of add-particle calls */

	/* Particle size value */
	real Pstep, Bstep, dx;			 /* Initial spacings for particles and boundary */

	/* Geometry parameters */
	string start_type, bound_type;
	int Scase, Bcase;				 /* What initial shape to take */
	StateVecD sim_start, bound_start;/* Starting coordinates (For box, bottom left) */
	int Bclosed, ghost;		  		 /* If the jet is closed or not */
	int Asource;                     /* Source of aerodynamic solution */

	/* Jet geometry parameters */
	real jet_diam;			        /* Jet diameter */
	real jet_depth;					/* Jet depth down - recommended > diam */
	StateVecD Angle;				/* Rotations in degrees */
	StateMatD Rotate;				/* Starting rotation matrix */ 
	StateMatD Transp;               /* Transpose of rotation matrix */
	real diam;                      /* Droplet diameter */
	uint nrad;                      /* Points along the radius of the jet/droplet */
	/* Box geometry parameters */
	StateVecI xyPART; 				/* Starting sim particles in x and y box */
	StateVecD sim_box, bound_box;	/* Box dimensions */
	
	/* Integration parameters */	
	uint subits;                    /* Max number of sub-iterations */
	uint Nframe; 			        /* Max number of frames to output */
	uint frame;						/* Current frame number */
	double t;                       /* Simulation time */
	real dt, framet;			    /* Timestep, frame times */
	real dt_max;					/* Maximum timestep */
	real dt_min;					/* Minimum timestep */
	real beta, gamma;				/* Newmark-Beta Parameters */
	real maxmu;                     /* Maximum viscosity component (CFL) */
	StateVecD grav;                 /* Gravity vector (assumed to be on z-axis) */

	uint framecount;                /* How many frames have been output */
	vector<size_t> back;            /* Particles at the back of the pipe */
	vector<vector<size_t>> buffer;  /* ID of particles inside the buffer zone */

	#if SIMDIM == 3
		VLM vortex;
	#endif
	StateVecD Force, AForce;					/*Total Force*/
	real mass;
	real tMom, aMom;
	
	/* Particle tracking settings */
	int eqOrder;
	real max_x_sph, max_x;
	size_t nSuccess, nFailed;
	real IPT_diam, IPT_area;
	uint cellsout, streakout, partout; /* Whether to output for particle tracking */
	ofstream cellfile, streakfile, partfile; /* File handles for particle tracking */

	/*Post Processessing settings*/
	uint afterSim;
	uint sliceOrSet;
	StateVecD planeOrigin, planeNorm; /* Plane conditions if in 3D */
	real cellSize;				      /* Size of each cell to make */
	StateVecD maxC, minC;			  /* Max and min coords to make grid */
	uint numVars, wrongDim;
	real postRadius;
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	/* Add default values */
	FLUID()
	{
		H = -1;
		HSQ = -1;
		sr = -1;
		Wdx = -1;

		rho0 = 1000;
		rhoJ = 1000;
		pPress = 0;

		alpha = 0.1;
		Cs = 300;
		mu = 8.94e-4;
		nu = mu/rho0;
		sig = 0.0708;
		gam = 7; B = -1;
		contangb = 150;
		delta = 0.1;
	}
	real H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
	real Wdx;                   /*Kernel value at the initial particle spacing distance.*/
	real rho0, rhoJ; 			/*Resting Fluid density*/
	real pPress;		/*Starting pressure in pipe*/
	
	real simM, bndM;			/*SPHPart and boundary masses*/
	real correc;				/*Smoothing Kernel Correction*/
	real alpha,Cs,mu,nu;		/*}*/
	real sig;					/* Fluid properties*/
	real gam, B; 				/*}*/
	
	real contangb;				/*Boundary contact angle*/
	real resVel;				/*Reservoir velocity*/
	//real front, height, height0;		/*Dam Break validation parameters*/

	/*delta SPH terms*/
	real delta; /*delta-SPH contribution*/
	real dCont; /*delta-SPH continuity constant term, since density and speed of sound are constant*/
	real dMom;  /*delta-SPH momentum constant term (dCont * rho0)*/

	// real maxU; /*Maximum particle velocity shift*/ (No longer needed)
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

		qInf = 0;
		vRef = 0;
		pRef = 101353;
		rhog = 1.29251;
		mug = 1.716e-5;
		T = 298;
		Rgas = 287;

		acase = 0;
		vStart = StateVecD::Zero();
		vInf = StateVecD::Zero();
		vJetMag = -1;
	}

	/*Gissler Parameters*/
	real L;
	real td;
	real omega;
	real tmax;
	real ycoef;
	real Cf, Ck, Cd, Cb, Cdef;

	real pVol;                     /*Volume of a particle*/
	real aPlate;                   /*Area of a plate*/

	/* Gas Properties*/
	real qInf, vRef, pRef;         /*Reference gas values*/
	real MRef;
	real rhog, mug;
	real gasM;					   /*A gas particle mass*/
	real T;						   /*Temperature*/
	real Rgas;                     /*Specific gas constant*/

	string aero_case;
	int acase;	                   /*Aerodynamic force case*/
	StateVecD vStart, vInf;          /*Jet & Freestream velocity*/
	real vJetMag;

	// real a;                        /*Case 3 tuning parameters*/
	// real b;                        /*}*/
	// real h1;                       /*}*/
	// real h2;                       /*}*/
}AERO;

typedef struct MESH
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*Standard contructor*/
	MESH(){}

	size_t size()
	{
		return elems.size();
	}

	/*Zone info*/
	std::string zone;
	size_t nPnts, nElem, nFace, nSurf;
	real scale;

	real maxlength;

	/*Point based data*/
	vector<StateVecD> verts;

	/*Face based data*/
	vector<vector<size_t>> faces;
	vector<std::pair<int,int>> leftright;
	vector<std::pair<size_t,int>> smarkers;

	/*Cell based data*/
	vector<vector<size_t>> elems;
	vector<StateVecD> cCentre;
	vector<vector<size_t>> cFaces;

	/*Solution vectors*/
	vector<StateVecD> cVel;
	vector<real> cCp;
	vector<real> cP;
	vector<real> cMass;
	vector<real> cRho;
	// vector<real> cRho_SPH;
	
	// Cell information for the momentum balance
	vector<StateVecD> cPertn;
	vector<StateVecD> cPertnp1;
	vector<real> cVol;


	vector<size_t> fNum; // number of fuel particles in cell
	vector<real> fMass;

	vector<StateVecD> vFn;
	vector<StateVecD> vFnp1;

}MESH;

/* Structure to track impacts and data for a surface marker of the mesh */
struct SURF
{
	SURF()
	{ 
		marker_count = 0;
		marker = 0;
		output = 0;
	}

	SURF(SURF const& in)
	{
		name = in.name;
		marker = in.marker;
		output = in.output;
		marker_count = 0;

		faceIDs = in.faceIDs;
		face_count = vector<uint>(faceIDs.size(),0);
	}

	/* Marker data */
	string name;
	int marker;
	int output;
	
	unsigned marker_count;
	real marker_beta; /* Collision efficiency (same as collection eff) */
	real marker_area;
	vector<size_t> marker_pIDs; /* ID for particles that hit the surface */
	vector<size_t> impacted_face; /* ID of the face the particle hit */
	vector<StateVecD> end_pos;
	vector<StateVecD> start_pos;

	/* Face data */
	vector<size_t> faceIDs;
	vector<uint> face_count;
	vector<real> face_beta;
	vector<real> face_area;
};


/*Container for delta-plus SPH calculations*/
typedef class DELTAP {
	public: 
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		DELTAP(size_t const& size)
		{
			L = vector<StateMatD>(size,StateMatD::Zero());
			gradRho = vector<StateVecD>(size,StateVecD::Zero());
			norm = vector<StateVecD>(size,StateVecD::Zero());
			avgV = vector<StateVecD>(size,StateVecD::Zero());
			
			lam = vector<real>(size,0.0);
			lam_nb = vector<real>(size,0.0);
			lam_ng = vector<real>(size,0.0);
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
			vector<real> const& lam_, vector<real> const& lam_nb_, vector<real> const& lam_ng_,
			vector<real> const& kernsum_)
		{
			L = L_; gradRho = gradRho_; norm = norm_; 
			avgV = avgV_; lam = lam_; lam_nb = lam_nb_; lam_ng = lam_ng_;
			kernsum = kernsum_;
		}

		void clear()
		{
			L.clear(); gradRho.clear(); norm.clear(); 
			avgV.clear(); lam.clear(); lam_nb.clear(); lam_ng.clear();
			kernsum.clear();

		}

		void erase(size_t const& start, size_t const& end)
		{
			L.erase(L.begin()+start, L.begin()+end);
			gradRho.erase(gradRho.begin()+start, gradRho.begin()+end);
			norm.erase(norm.begin()+start, norm.begin()+end);
			avgV.erase(avgV.begin()+start, avgV.begin()+end);
			lam.erase(lam.begin()+start, lam.begin()+end);
			lam_nb.erase(lam_nb.begin()+start, lam_nb.begin()+end);
			lam_ng.erase(lam_ng.begin()+start, lam_ng.begin()+end);
			kernsum.erase(kernsum.begin()+start, kernsum.begin()+end);
		}

		vector<StateMatD> L;
		vector<StateVecD> gradRho;
		vector<StateVecD> norm;
		vector<StateVecD> avgV;

		vector<real> lam;		/* Eigenvalues with all particles considered */
		vector<real> lam_nb; 	/* Eigenvalues without boundary particles (and ghost particles) */
		vector<real> lam_ng; 	/* Eigenvalues without ghost particles */
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

/*SPHPart data class*/
typedef struct SPHPart {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	SPHPart(StateVecD const& X, StateVecD const& Vi, real const Rhoi, real const Mi, 
		real const press, int const bound, uint const pID)
	{
		partID = pID; cellID = 0; faceID = 0;
		b = bound; surf = 0;

		xi = X;	v = Vi; acc = StateVecD::Zero(); Af = StateVecD::Zero();
		Rrho = 0.0; rho = Rhoi; p = press; m = Mi; 
		s = 0.0; woccl = 0.0; pDist = 0.0;
					
		cellV = StateVecD::Zero();
		cellP = 0.0; cellRho = 0.0;
		internal = 0;

		vPert = StateVecD::Zero(); 
	}

	/*To add particles dynamically for boundary layer*/
	SPHPart(StateVecD const& X, SPHPart const& pj, int const bound, size_t const pID)
	{
		partID = pID; cellID = 0; faceID = 0;
		b = bound; surf = 0;

		xi = X;	v = pj.v; acc = StateVecD::Zero(); Af = StateVecD::Zero();
		Rrho = 0.0; rho = pj.rho; p = pj.p; m = pj.m;
		s = 0.0; woccl = 0.0; pDist = 0.0;

		cellV = StateVecD::Zero();
		cellP = 0.0; cellRho = 0.0;
		internal = 0;

		vPert = StateVecD::Zero();
	}

	SPHPart(){};

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
	StateVecD xi, v, acc, Af;
	real Rrho, rho, p, m, s, woccl, pDist;
	StateVecD cellV;
	real cellP, cellRho;
	uint internal; 

	/*Mesh surface repulsion */
	StateVecD bNorm;
	real y;
	
	StateVecD vPert;		
}SPHPart;

typedef struct IPTPart
{
	public:
	IPTPart()
	{
		partID = 0;
		going = 1;
		nIters = 0;
		nNotFound = 0;
		failed = 0;
		t = 0; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = StateVecD::Zero();
		faceRho = 0.0;
	
		cellID = -3;
		cellV = StateVecD::Zero();
		cellRho = 0.0; 
			
		acc = 0.0;
		v = StateVecD::Zero();
		xi = StateVecD::Zero();
		
		mass = 0.0;
		d = 0.0;
		A = 0.0;
	}

	IPTPart(StateVecD const xi_, real const mass_, real const d_)
	{
		partID = 0;
		going = 1;
		nIters = 0;
		nNotFound = 0;
		failed = 0;
		t = 0; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = StateVecD::Zero();
		faceRho = 0.0;
	
		cellID = -3;
		cellV = StateVecD::Zero();
		cellRho = 0.0; 
			
		acc = 0.0;
		v = StateVecD::Zero();
		xi = xi_;
		
		mass = mass_;
		d = d_;
		A = M_PI * d_*d_/4.0; 
	}

	IPTPart(StateVecD const xi_, IPTPart const& pi_, size_t const& pID_)
	{
		partID = pID_;
		going = 1;
		nIters = 0;
		nNotFound = 0;
		failed = 0;
		t = 0; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = StateVecD::Zero();
		faceRho = 0.0;
	
		cellID = pi_.cellID;
		cellV = pi_.cellV; 
		cellRho = pi_.cellRho; 
			
		acc = pi_.acc;
		v = pi_.v;
		xi = xi_;
		
		mass = pi_.mass;
		d = pi_.d;
		A = pi_.A; 
	}

	IPTPart(SPHPart const& pi, real const& time, real const& diam, real const& area)
	{
		partID = pi.partID;
		going = 1;
		nIters = 0;
		nNotFound = 0;
		failed = 0;
		t = time; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = StateVecD::Zero();
		faceRho = 0.0;

		/* Make the current values the same as those of the SPH particle */
		cellID = pi.cellID;
		cellV = pi.cellV; 
		cellRho = pi.cellRho; 

		acc = 0.0;
		v = pi.v;
		xi = pi.xi;
		
		mass = pi.m;
		/* Derive diameter and area from the mass and resting density */
		d = diam;
		A = area;
	}

	size_t partID;
	uint going; /* Is particle still being integrated */
	uint nIters; /* How many integration steps it's gone through */
	size_t nNotFound;
	size_t failed;

	/* Timestep properties */
	real t, dt;

	/* Face properties */
	uint faceID;
	StateVecD faceV;
	real faceRho;

	/* Containing cell properties */
	uint cellID;
	StateVecD cellV;
	real cellRho;

	real acc; /* Implicit acceleration */
	StateVecD v, xi; /* State variables */

	real mass; /* SPHPart mass */
	real d, A; /* SPHPart diameter and area */
}IPTPart;


typedef std::vector<SPHPart> SPHState;
typedef std::vector<IPTPart> IPTState;
typedef std::vector<SURF> SURFS;
/* Neighbour search tree containers */
typedef std::vector<std::vector<std::pair<size_t,real>>> OUTL;
typedef std::vector<std::vector<size_t>> celll;
typedef KDTreeVectorOfVectorsAdaptor<SPHState,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

typedef struct KDTREE
{
	KDTREE(SPHState const& pnp1, MESH const& cells): NP1(SIMDIM,pnp1,20), 
	CELL(SIMDIM,cells.cCentre,20)/* , BOUNDARY(SIMDIM,cells.bVerts,20) */ {}

	Sim_Tree NP1;
	Vec_Tree CELL;
	// Vec_Tree BOUNDARY;	

}KDTREE;

#endif /* VAR_H */
