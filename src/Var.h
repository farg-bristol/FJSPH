/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VAR_H
#define VAR_H

#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <omp.h>

// Third party includes
#include <Eigen/Core>
#include <Eigen/StdVector>
#include "Third_Party/NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"

using std::vector;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string; 
using std::setw;

#ifdef DEBUG
    /*Open debug file to write to*/
    extern FILE* dbout;
#endif

#ifdef DAMBREAK
std::ofstream dambreak("Dam_Data.log",std::ios::out);
#endif
// std::ofstream pertLog("cellPert.log",std::ios::out);

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 2
#endif

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#define M_PI_4 M_PI/4.0
#endif

/* Define data type. */
#ifndef FOD
#define FOD 1 /*0 = float, 1 = double*/
#endif

#if FOD == 1
// #define real double
typedef double real;
int32_t const static realType = 2;
#else
// #define real float
typedef float real;
int32_t const static realType = 1;
#endif

/* Define the default surface tension model */
#if !defined(HESF) && !defined(CSF) && !defined(PAIRWISE) && !defined(NOST) 
#define PAIRWISE
#endif

#if !defined(ISOEOS) && !defined(COLEEOS)
#define COLEEOS
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

const std::string WHITESPACE = " \n\r\t\f\v";

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

#pragma omp declare reduction(min: StateVecD : omp_out = \
            omp_out.norm() < omp_in.norm() ? omp_in : omp_out)\
                    initializer(omp_priv = omp_orig) 

#pragma omp declare reduction(max: StateVecD : omp_out = \
            omp_out.norm() > omp_in.norm() ? omp_in : omp_out)\
                    initializer(omp_priv = omp_orig) 


/*Define particle type indexes*/
enum partType{BOUND=0,PISTON,BUFFER,BACK,PIPE,GHOST,FREE,OUTLET,LOST};

/* Define shapes for boundary and fluid initialisation */
enum shapeType{NONE=0,BOX,CYLINDER,SPHERE,JET,CONVJET};

// Shapes for new intitialisation
enum shape_type {fineLine = 0, linePlane, squareCube, circleSphere, cylinder, arcSection, coordDef,
                inletZone, hollow, solid};

enum solve_type {newmark_beta = 0, runge_kutta, DBC, pressure_G, ghost};

enum aero_force {none = 0, Gissler, Induced_Pressure, SkinFric};

enum aero_source {constVel = 0, meshInfl, VLMInfl};

real const static default_val = 9999999.0;

inline bool check_vector(StateVecD const& v)
{
    return (v[0] != default_val && v[1] != default_val 
            #if SIMDIM == 3
            && v[2] != default_val 
            #endif
            );
}

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

// PState PartState;

inline real random(int const& interval)
{
    return real(rand() % interval - interval / 2) * MERROR;
}

struct SIM 
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SIM()
    {   
        // If OMP_NUM_THREADS exists, use it to set the threads, 
        // but if not, use get num threads.
        if (std::getenv("OMP_NUM_THREADS"))
            numThreads = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
        else
            numThreads = omp_get_num_threads();
        

        /* Input files */
        CDForFOAM = 0;
        isBinary = 0; buoyantSim = 0; incomp = 0;
        labelSize = 32; scalarSize = 64;
        #if SIMDIM == 3
        offset_axis = 0;
        #else
        offset_axis = 2; /* Default to XZ plane in 2D */
        #endif 
        
        angle_alpha = 0;

        boundFile = NULL;
        fluidFile = NULL;

        write_tecio = 1;
        write_h5part = 0;
        out_encoding = 1; /* default to binary */
        gout = 0; /* Default to restartable data */
        single_file = 1;
        restart_header_written = 0;
        restart = 0;
        scale = 1.0;

        nbound = 0; nfluid = 0;
        partID = 0;
        simPts = 0; bndPts = 0; totPts = 0;
        gstPts = 0; delNum = 0; 
        intNum = 0;	nrefresh = 0; addcount = 0;
        finPts = 99999999; /* Default to basically no particle limit */
    
        Pstep = -1.0; Bstep = 0.7; dx = -1.0;

        bound_type = "(none)"; /* Default to no boundary */
        Scase = 0; Bcase = 0; ghost = 0;
        Asource = constVel;
        init_hydro_pressure = 0;
        use_global_gas_law = 1;
        hydro_height = -1;

        /* Universal geometry parameters */
        offset_vec = StateVecD::Zero();
        
        /* Droplet geometry parameters */
        nrad = 0;
        diam = -1.0;
        
        /* Integration parameters */
        solver_type = newmark_beta;
        subits = 20;
        Nframe = -1;
        frame = 0;
        bound_solver = 0;
        nStable = 0;
        nStable_Limit = 10;
        nUnstable = 0;
        nUnstable_Limit = 3;
        subits_factor = 0.333;
        minRes = -7.0;
        cfl = 1.0;
        cfl_step = 0.05;
        cfl_max = 2.0;
        cfl_min = 0.1;
        t = 0.0;
        tframem1 = 0.0;
        dt = 2e-10;
        framet = -1;
        dt_max = 1;
        dt_min = 0.0;
        beta = 0.25; gamma = 0.5;
        maxmu = 0;
        maxshift = 9999999;

        /*Gravity Vector*/
        #if SIMDIM == 3
            grav = StateVecD(0.0,0.0,-9.81);
        #else
            grav = StateVecD(0.0,-9.81);
        #endif

        framecount = 0;
        restart_tol = 0.001;
        Force = StateVecD::Zero(); AForce = StateVecD::Zero();					/*Total Force*/
        mass = 0;
        tMom = 0; aMom = 0;
    
        /* SPHPart tracking settings */
        using_ipt = 0;
        eqOrder = 2;
        max_x_sph = 9999999;
        max_x = 9999999;
        nSuccess = 0; nFailed = 0;
        IPT_diam = 0.0; IPT_area = 0.0;
        relax = 0.6; nrelax = 5;
        cellsout = 0; streakout = 1; partout = 0;

        dropDragSweep = 0;
        speedTest = 0;
        nRuns = 1;
    }

    uint numThreads;

    /* File input parameters */
    uint CDForFOAM;
    std::string infile, output_prefix, outdir;
    std::string fluidfile, boundfile;
    std::string restart_prefix;
    std::string vlm_file;

    /* OpenFOAM files */
    std::string foamdir, foamsol;
    int isBinary, buoyantSim, incomp, labelSize, scalarSize;
    
    /* TAU files */
    std::string taumesh, taubmap, tausol;
    int offset_axis;
    real angle_alpha;
    vector<int> markers;
    vector<string> bnames;
    vector<int> bwrite; /* If surface wants data writing */

    /* TECIO File outputs */
    void* boundFile; /*TECIO file handles*/
    void* fluidFile;
    void* ghostFile;
    vector<int32_t> varTypes;
    std::string output_names;   	/**< Names of variables to output */
    std::vector<uint> outvar; 	/**< Variables to output */
    std::vector<int32_t> var_types; /**< Data types for output */
    std::string var_names;      /**< Variable names for tecplot output */

    /* H5Part file outputs */
    int64_t ffile, bfile;
	uint write_tecio;			/**< Write Tecplot file (this or write_h5part needs to be true) */
	uint write_h5part;			/**< Write H5part file (this or write_tecio needs to be true) */

    ofstream surfacefile; /* Surface impact file */
    void* surfaceHandle;

    /* Output type */
    uint out_encoding;              /* ASCII or binary output*/
    uint gout;                      /* Output type.*/
    uint single_file;               /* Use single file or not for output */
    uint restart_header_written;    /* Whether the restart header has been written yet */
    uint restart;					/* If starting from existing solution */
    real scale;						/* Simulation scale */

    /* SPHPart counts */
    size_t nbound;					 /* Number of boundary blocks */ 
    size_t nfluid;					 /* Number of fluid blocks */
    size_t partID;                   /* Track the particle ID separately */
    size_t totPts;					 /* Total count */
    size_t simPts;					 /* Fluid count */
    size_t bndPts;			 		 /* Boundary count */
    size_t gstPts;                   /* Ghost count */
    size_t finPts;					 /* End count */
    size_t delNum;                   /* Number of deleted particles */
    size_t intNum;                   /* Number of internal particles */
    size_t nrefresh; 				 /* last add call particle number */
    uint addcount;					 /* Current Number of add-particle calls */

    /* Particle size value */
    real Pstep, Bstep, dx;			 /* Initial spacings for particles and boundary */

    /* Geometry parameters */
    string start_type, bound_type;
    int Scase, Bcase;				  /* What initial shape to take */
    StateVecD offset_vec;             /* Global offset coordinate */
    int ghost;		  		          /* If the jet is closed or not */
    int Asource;                      /* Source of aerodynamic solution */
    int init_hydro_pressure;          /* Initialise fluid with hydrostatic pressure? */
    int use_global_gas_law;			  /* Whether to use block specific gas laws (Currently not active) */
    real hydro_height;                /* Hydrostatic height to initialise using */

    /* Jet geometry parameters */
    real diam;                      /* Droplet diameter */
    uint nrad;                      /* Points along the radius of the jet/droplet */

    /* Integration parameters */	
    string solver_name;             /* Name of solver */
    uint solver_type;               /* Use Runge-Kutta or Newmark-Beta */
    uint subits;                    /* Max number of sub-iterations */
    uint Nframe; 			        /* Max number of frames to output */
    uint frame;						/* Current frame number */
    uint bound_solver;              /* Use boundary pressure or ghost particles */
    uint nStable;                   /* Count for number of overly stable timesteps to alter CFL */
    uint nStable_Limit;             /* Limit before changing CFL */
    uint nUnstable;                 /* Count for number of unstable timestep to alter CFL */
    uint nUnstable_Limit;           /* Limit before changing CFL */
    real subits_factor;             /* Factor * max subits, under which is considered overly stable CFL */
    double cfl;                     /* CFL criterion number */
    double cfl_step;                /* CFL step to perform if unstable */
    double cfl_max;                 /* Maximum CFL for the simulation */
    double cfl_min;                 /* Minimum CFL for the simulation */
    double t;                       /* Simulation time */
    double tframem1;				/* Last written frame time */
    real minRes;                    /* Minimum solver residual */
    real dt, framet;			    /* Timestep, frame times */
    real dt_max;					/* Maximum timestep */
    real dt_min;					/* Minimum timestep */
    real beta, gamma;				/* Newmark-Beta Parameters */
    real maxmu;                     /* Maximum viscosity component (CFL) */
    real maxshift;					/* Maximum shifting velocity */
    StateVecD grav;                 /* Gravity vector (assumed to be on z-axis) */

    uint framecount;                /* How many frames have been output */
    real restart_tol;               /* Tolerance on buffer particles fitting*/

    StateVecD Force, AForce;		/*Total Force*/
    real mass;
    real tMom, aMom;
    
    /* Particle tracking settings */
    int using_ipt;
    int eqOrder;
    real max_x_sph, max_x; /* Transition distance for sph (naive assumption) and termination for IPT */
    size_t nSuccess, nFailed;
    real IPT_diam, IPT_area;
    real relax, nrelax; // Relaxation parameters
    uint cellsout, streakout, partout; /* Whether to output for particle tracking */
    FILE* cellfile; 
    FILE* streakfile;
    FILE* partfile; /* File handles for particle tracking */
    void* cellHandle; /*TECIO file handles*/
    void* streakHandle;
    void* partHandle;

    /* Droplet sweep parameters */
    int dropDragSweep;
    vector<size_t> nacross;
    vector<real> diameters;
    vector<real> velocities;
    vector<real> Reynolds;

    /* Speed test options */
    int speedTest;
    int nRuns;

    /*Post Processessing settings*/
    uint afterSim;
    uint sliceOrSet;
    StateVecD planeOrigin, planeNorm; /* Plane conditions if in 3D */
    real cellSize;				      /* Size of each cell to make */
    StateVecD maxC, minC;			  /* Max and min coords to make grid */
    uint numVars, wrongDim;
    real postRadius;
};

/*Fluid and smoothing parameters*/
struct FLUID 
{
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /* Add default values */
    FLUID()
    {
        H = -1.0;
        HSQ = -1.0;
        sr = -1.0;
        Hfac = 2.0;
        Wdx = -1.0;

        rho0 = 1000.0;
        rhoJ = 1000.0;
        pPress = 0.0;
        backP = 0.0;
        rhoMax = 1500;	/* Allow a large variation by default */
        rhoMin = 500;	/* making it essentially unbounded other */
        rhoVar = 50.0;	/* than truly unphysical pressures */
        rhoMaxIter = 1.0;

        simM = -1.0;
        bndM = -1.0;
        correc = -1.0;
        alpha = 0.1;
        Cs = 300.0;
        mu = 8.94e-4;
        nu = mu/rho0;
        sig = 0.0708;
        gam = 7.0; B = -1.0;
        contangb = 150.0;
        resVel = 0.0;
        delta = 0.1;
        dCont = -1.0;
        dMom = -1.0;
    }

    real H, HSQ, sr; 			/* Support Radius, SR squared, Search radius*/
    real Hfac;
    real Wdx;                   /* Kernel value at the initial particle spacing distance.*/
    real rho0, rhoJ; 			/* Resting Fluid density*/
    real pPress;				/* Starting pressure in pipe*/
    real backP;					/* Background pressure */
    real rhoMax, rhoMin;		/* Minimum and maximum pressure */
    real rhoVar;				/* Maximum density variation allowed (in %) */
    real rhoMaxIter;            /* Maximum value to reduce the timestep */

    real simM, bndM;			/* SPHPart and boundary masses*/
    real correc;				/* Smoothing Kernel Correction*/
    real alpha,Cs,mu,nu;		/*}*/
    real sig;					/* Fluid properties*/
    real gam, B; 				/*}*/
    
    real contangb;				/* Boundary contact angle*/
    real resVel;				/* Reservoir velocity*/
    //real front, height, height0;		/* Dam Break validation parameters*/

    /*delta SPH terms*/
    real delta; /*delta-SPH contribution*/
    real dCont; /*delta-SPH continuity constant term, since density and speed of sound are constant*/
    real dMom;  /*delta-SPH momentum constant term (dCont * rho0)*/

    // real maxU; /*Maximum particle velocity shift*/ (No longer needed)
};

/*Aerodynamic Properties*/
struct AERO
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    AERO()
    {
        vInf = StateVecD::Zero();

        Cf = 1.0/3.0;
        Ck = 8.0;
        Cd = 5.0;
        Cb = 0.5;

        vRef = 0.0;
        pRef = 101353.0;
        rhog = 1.29251;
        mug = 1.716e-5;
        T = 298.0;
        Rgas = 287.0;
        gamma = 1.403;
        sos = sqrt(T*Rgas*gamma);
        isossqr = 1.0/(sos*sos);
        MRef = -1.0;

        vJetMag = -1.0;
        cutoff = 0.75;
        interp_fac = 2.0;
        iinterp_fac = 0.5;
        #if SIMDIM == 2
        nfull = 30;
        #else
        nfull = 200;
        #endif
        infull = 1.0/nfull;
        incount = 1.0/(iinterp_fac * nfull);
        acase = 0;
        use_lam = 1;
        useDef = 0;
        use_dx = 0;
    }

    void GetYcoef(const FLUID& fvar, const real diam)
    {
        #if SIMDIM == 3
            L = diam * std::cbrt(3.0/(4.0*M_PI));
            aSphere = M_PI*L*L;
            aPlate = diam*diam;
        #else
            L = diam/sqrt(M_PI);
            aSphere = 2*L;
            aPlate = diam;
        #endif

        td = (2.0*fvar.rho0*pow(L,SIMDIM-1))/(Cd*fvar.mu);

        omega = sqrt((Ck*fvar.sig)/(fvar.rho0*pow(L,SIMDIM))-1.0/pow(td,2.0));

        tmax = -2.0 *(atan(sqrt(pow(td*omega,2.0)+1)
                        +td*omega) - M_PI)/omega;

        Cdef = 1.0 - exp(-tmax/td)*(cos(omega*tmax)+
            1/(omega*td)*sin(omega*tmax));
        ycoef = 0.5*Cdef*(Cf/(Ck*Cb))*(rhog*L)/fvar.sig;

        //cout << ycoef << "  " << Cdef << "  " << tmax << "  " << endl;
    }

    StateVecD vInf;        /* Freestream velocity*/
    
    /*Gissler Parameters*/
    real L;
    real td;
    real omega;
    real tmax;
    real ycoef;
    real Cf, Ck, Cd, Cb, Cdef;

    real pVol;                     /* Volume of a particle */
    real aSphere;				   /* Area of a sphere/cylinder */
    real aPlate;                   /* Area of a plate*/

    /* Gas Properties*/
    real vRef, pRef;               /* Reference gas values*/
    real rhog, mug;
    real gasM;					   /* A gas particle mass*/
    real T;						   /* Temperature*/
    real Rgas;                     /* Specific gas constant*/
    real gamma;                    /* Ratio of specific heats */
    real sos;                      /* Speed of sound */
    real isossqr;				   /* Inverse speed of sound squared */
    real MRef;    
    
    real vJetMag;
    real cutoff;				   /* Lambda value cutoff */
    real interp_fac;               /* Factor for the interpolation value */
    real iinterp_fac;              /* Inverse factor for the interpolation */
    real nfull;                    /* 2/3ds full neighbourhood count */
    real infull;                   /* inverse of nfull */
    real incount;                  /* inverse nfull * factor */

    string aero_case;
    int acase;	                   /* Aerodynamic force case*/
    int use_lam;                    
    int useDef;
    int use_dx;
};

struct MESH
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /*Standard contructor*/
    MESH():nPnts(0), nElem(0), nFace(0), nSurf(0), scale(1.0), maxlength(0.0),minlength(9999999.9){}

    inline size_t size()
    {
        return nElem;
    }

    inline void alloc(size_t nPnts_, size_t nElem_, size_t nFace_, size_t nSurf_)
    {
        nPnts = nPnts_; nElem = nElem_; nFace = nFace_; nSurf = nSurf_;
        verts = vector<StateVecD>(nPnts);
        faces.reserve(nFace*2);
        leftright.reserve(nFace);

        // elems = vector<vector<size_t>>(nElem);
        cCentre = vector<StateVecD>(nElem);
        cFaces = vector<vector<size_t>>(nElem);

        cVel = vector<StateVecD>(nElem);
        cP = vector<real>(nElem);
        cRho = vector<real>(nElem);
    }

    /*Zone info*/
    std::string zone;
    size_t nPnts, nElem, nFace, nSurf;
    real scale;

    real maxlength;
    real minlength;

    /*Point based data*/
    vector<StateVecD> verts;

    /*Face based data*/
    vector<vector<size_t>> faces;
    vector<std::pair<int,int>> leftright;
    vector<std::pair<size_t,int>> smarkers;

    /*Cell based data*/
    // vector<vector<size_t>> elems;
    vector<StateVecD> cCentre;
    vector<vector<size_t>> cFaces;

    /*Solution vectors*/
    vector<StateVecD> cVel;
    vector<real> cCp;
    vector<real> cP;
    vector<real> cRho;
    // vector<real> cRho_SPH;
    
    // Cell information for the momentum balance
    vector<StateVecD> cPertn;
    vector<StateVecD> cPertnp1;
    // vector<real> cVol;
    // vector<real> cMass;


    vector<size_t> fNum; // number of fluid particles in cell
    vector<real> fMass;

    vector<StateVecD> vFn;
    vector<StateVecD> vFnp1;

};

/* Structure to track impacts and data for a surface marker of the mesh */
struct SURF
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SURF()
    { 
        marker_count = 0;
        marker_beta = 0.0;
        marker_area = 0.0;
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

/*SPHPart data class*/
struct SPHPart 
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SPHPart(StateVecD const& X, StateVecD const& Vi, real const Rhoi, real const Mi, 
        real const press, int const bound, uint const pID)
    {
        partID = pID; faceID = 0; cellID = -1;
        b = bound; surf = 0; surfzone = 0; nFailed = 0;

        xi = X;	v = Vi; acc = StateVecD::Zero(); Af = StateVecD::Zero();
        aVisc = StateVecD::Zero(); norm = StateVecD::Zero();
        Rrho = 0.0; rho = Rhoi; p = press; m = Mi; 
        curve = 0.0; s = 0.0; woccl = 0.0; pDist = 0.0; deltaD = 0.0;
                    
        cellV = StateVecD::Zero();
        cellP = 0.0; cellRho = 0.0;
        internal = 0;

        L = StateMatD::Zero(); gradRho = StateVecD::Zero(); norm = StateVecD::Zero();
        colourG = 0.0; colour = 0.0; lam = 0.0; lam_nb = 0.0; lam_ng = 0.0; kernsum = 0.0;
        bNorm = StateVecD::Zero(); y = 0.0;
        vPert = StateVecD::Zero(); 
    }

    /*To add particles dynamically for fictitious particles*/
    SPHPart(StateVecD const& X, SPHPart const& pj, int const bound, size_t const pID)
    {
        partID = pID; faceID = 0; cellID = -1;
        b = bound; surf = 0; surfzone = 0; nFailed = 0;

        xi = X;	v = pj.v; acc = StateVecD::Zero(); Af = StateVecD::Zero();
        aVisc = StateVecD::Zero(); norm = StateVecD::Zero();
        Rrho = 0.0; rho = pj.rho; p = pj.p; m = pj.m;
        curve = 0.0; s = 0.0; woccl = 0.0; pDist = 0.0; deltaD = 0.0;

        cellV = StateVecD::Zero();
        cellP = pj.cellP; cellRho = pj.cellRho;
        internal = 0;

        L = StateMatD::Zero(); gradRho = StateVecD::Zero(); norm = StateVecD::Zero();
        colourG = 0.0; colour = 0.0; lam = 0.0; lam_nb = 0.0; lam_ng = 0.0; kernsum = 0.0;
        bNorm = StateVecD::Zero(); y = 0.0;
        vPert = StateVecD::Zero();
    }

    SPHPart():partID(0),faceID(0),cellID(-1),b(0),surf(0),surfzone(0),nFailed(0),xi(StateVecD::Zero()),
         v(StateVecD::Zero()), acc(StateVecD::Zero()), Af(StateVecD::Zero()), aVisc(StateVecD::Zero()),
         Rrho(0.0),rho(0.0),p(0.0),m(0.0),curve(0.0),s(0.0),woccl(0.0),pDist(0.0),deltaD(0.0),
         cellV(StateVecD::Zero()),cellP(0.0),cellRho(0.0),internal(0),
         L(StateMatD::Zero()), gradRho(StateVecD::Zero()), norm(StateVecD::Zero()),  
         colourG(0.0), colour(0.0), lam(0.0), lam_nb(0.0), lam_ng(0.0), kernsum(0.0),
         bNorm(StateVecD::Zero()),y(0.0), vPert(StateVecD::Zero()) {};


    inline int size() const
    {	/*For neighbour search, return size of xi vector*/
        return(xi.size());
    }

    inline real operator[](int a) const
    {	/*For neighbour search, return index of xi vector*/
        return(xi[a]);
    }

    size_t partID, faceID;
    long int cellID; // Make it signed so that it can be used for checks.
    uint b; //What state is a particle. See PartState above for possible options
    uint surf; /*Is a particle a surface? 1 = yes, 0 = no*/
    uint surfzone; /* Is particle in the surface area? */
    uint nFailed; /* How many times has the containment query failed */

    StateVecD xi, v, acc, Af, aVisc;
    real Rrho, rho, p, m, curve, s, woccl, pDist, deltaD;
    StateVecD cellV;
    real cellP, cellRho;
    uint internal; 

    StateMatD L;
    StateVecD gradRho;
    StateVecD norm;

    real colourG;
    real colour;    /* Kernel sum with volume considered. */
    real lam;		/* Eigenvalues with all particles considered */
    real lam_nb; 	/* Eigenvalues without boundary particles (and ghost particles) */
    real lam_ng; 	/* Eigenvalues without ghost particles */
    real kernsum;	/* Summation of the kernel */

    /*Mesh surface repulsion */
    StateVecD bNorm;
    real y;
    
    StateVecD vPert;		
};

struct IPTPart
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    IPTPart()
    {
        partID = 0;
        going = 1;
        nIters = 0;
        failed = 0;
        t = 0; dt = 0;

        /* Set initial IDs to a nonsense value, so that they don't interfere */
        faceID = -1; 
        faceV = StateVecD::Zero();
        faceRho = 0.0;
    
        cellID = -3;
        cellV = StateVecD::Zero();
        cellRho = 0.0; 
            
        acc = 0.0;
        relax = 1.0;
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
        failed = 0;
        t = 0; dt = 0;

        /* Set initial IDs to a nonsense value, so that they don't interfere */
        faceID = -1; 
        faceV = StateVecD::Zero();
        faceRho = 0.0;
    
        cellID = -3;
        cellV = StateVecD::Zero();
        cellRho = 0.0; 
            
        acc = 0.0;
        relax = 1.0;
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
        failed = 0;
        t = 0; dt = 0;

        /* Set initial IDs to a nonsense value, so that they don't interfere */
        faceID = -1; 
        faceV = StateVecD::Zero();
        faceRho = 0.0;
    
        cellID = pi_.cellID;
        cellV = pi_.cellV; 
        cellRho = pi_.cellRho; 
            
        acc = pi_.acc;
        relax = pi_.relax;
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
        failed = 0;
        t = time; dt = 0;

        /* Set initial IDs to a nonsense value, so that they don't interfere */
        faceID = -1; 
        faceV = StateVecD::Zero();
        faceRho = 0.0;

        /* Make the current values the same as those of the SPH particle */
        cellID = pi.cellID;
        cellV = pi.cellV; 
        cellRho = pi.cellRho; 

        acc = 0.0;
        relax = 1.0;
        v = pi.v;
        xi = pi.xi;
        
        mass = pi.m;
        /* Derive diameter and area from the mass and resting density */
        d = diam;
        A = area;
    }

    void reset(SPHPart const& pi, real const& time)
    {
        partID = pi.partID;
        going = 1;
        nIters = 0;
        failed = 0;
        t = time; dt = 0;

        /* Set initial IDs to a nonsense value, so that they don't interfere */
        faceID = -1; 
        faceV = StateVecD::Zero();
        faceRho = 0.0;

        /* Make the current values the same as those of the SPH particle */
        cellID = pi.cellID;
        cellV = pi.cellV; 
        cellRho = pi.cellRho; 

        acc = 0.0;
        relax = 1.0;
        v = pi.v;
        xi = pi.xi;
    }

    size_t partID;
    uint going; /* Is particle still being integrated */
    uint nIters; /* How many integration steps it's gone through */
    size_t failed;

    /* Timestep properties */
    real t, dt;

    /* Face properties */
    uint faceID;
    StateVecD faceV;
    real faceRho;

    /* Containing cell properties */
    long int cellID;
    StateVecD cellV;
    real cellRho;

    real acc; /* Implicit acceleration */
    real relax;
    StateVecD v, xi; /* State variables */

    real mass; /* SPHPart mass */
    real d, A; /* SPHPart diameter and area */
};

// Structure to define the inlet and outlet conditions for either fluid or boundaries
struct bound_block
{
    bound_block(size_t const& a) // Initialise starting index only first
    {
        index.first = a;
        index.second = a;
        delete_norm = StateVecD::Constant(default_val);
        delconst = default_val;
        insert_norm = StateVecD::Constant(default_val);
        insconst = default_val;
        aero_norm = StateVecD::Constant(default_val);
        aeroconst = default_val;
        nTimes = 0;
        fixed_vel_or_dynamic = 0;
        hcpl = 0;
        bound_solver = 0;
        no_slip = 0;
        block_type = 0;

    }

    bound_block(size_t const& a, size_t const& b)
    {
        index.first = a;
        index.second = b;
        delete_norm = StateVecD::Constant(default_val);
        delconst = default_val;
        insert_norm = StateVecD::Constant(default_val);
        insconst = default_val;
        aero_norm = StateVecD::Constant(default_val);
        aeroconst = default_val;
        nTimes = 0;
        fixed_vel_or_dynamic = 0;
        hcpl = 0;
        bound_solver = 0;
        no_slip = 0;
        block_type = 0;
    }

    string name;
    std::pair<size_t,size_t> index;
    
    StateVecD delete_norm; // Deletion plane
    real delconst;
    StateVecD insert_norm; // Insertion plane
    real insconst;
    StateVecD aero_norm;
    real aeroconst;

    vector<size_t> back;            /* Particles at the back of the pipe */
    vector<vector<size_t>> buffer;  /* ID of particles inside the buffer zone */

    vector<real> times;
    vector<StateVecD> vels;

    size_t nTimes;
    int fixed_vel_or_dynamic; //Used for inlets, but also for boundaries to identify if moving or static
    int hcpl;
    int no_slip;
    int bound_solver;
    int block_type;
};

typedef std::vector<bound_block> LIMITS;
typedef std::vector<SPHPart> SPHState;
typedef std::vector<IPTPart> IPTState;
typedef std::vector<SURF> SURFS;

/* Neighbour search tree containers */
typedef nanoflann::ResultItem<size_t,real> neighbour_index;
typedef std::vector<std::vector<neighbour_index>> OUTL;
typedef std::vector<std::vector<size_t>> celll;
typedef KDTreeVectorOfVectorsAdaptor<SPHState,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,real,SIMDIM,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

// Default to no approximate neighbours, and unsorted list.
nanoflann::SearchParameters const flann_params(0, false); 

// struct KDTREE
// {
// 	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

// 	KDTREE(SPHState const& pnp1, MESH const& cells): NP1(SIMDIM,pnp1,20), 
// 	CELL(SIMDIM,cells.cCentre,20)/* , BOUNDARY(SIMDIM,cells.bVerts,20) */ {}

// 	Sim_Tree NP1;
// 	Vec_Tree CELL;
// 	// Vec_Tree BOUNDARY;	

// };

#endif /* VAR_H */
