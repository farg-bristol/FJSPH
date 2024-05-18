/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef VAR_H
#define VAR_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>

// Third party includes
#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"
#include "Third_Party/NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::string;
using std::vector;

#ifdef DEBUG
/*Open debug file to write to*/
extern FILE* dbout;
#endif

#ifdef DAMBREAK
std::ofstream dambreak("Dam_Data.log", std::ios::out);
#endif
// std::ofstream pertLog("cellPert.log",std::ios::out);

/* Define Simulation Dimension */
#ifndef SIMDIM
#define SIMDIM 2
#endif

/* Define the default surface tension model */
#if !defined(HESF) && !defined(CSF) && !defined(PAIRWISE) && !defined(NOST)
#define PAIRWISE
#endif

#if !defined(ISOEOS) && !defined(COLEEOS)
#define COLEEOS
#endif

// Define pi
#ifndef M_PI
#define M_PI (4.0 * atan(1.0))
#define M_PI_4 M_PI / 4.0
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

// Tecplot data types
int32_t const static int32Type = 3;
int32_t const static int16Type = 4;
int32_t const static uint8Type = 5;

typedef unsigned int uint;

/*Get machine bit precision for Simulation of Simplicity*/
#ifndef MEPSILON
#define MEPSILON std::numeric_limits<real>::epsilon() /*For  float, power is -24*/
#endif

#ifndef MERROR
#define MERROR (7 * MEPSILON + 56 * MEPSILON * MEPSILON)
#endif

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<real, SIMDIM, 1> StateVecD;
typedef Eigen::Matrix<int, SIMDIM, 1> StateVecI;
typedef Eigen::Matrix<real, SIMDIM, SIMDIM> StateMatD;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<real, SIMDIM + 1, 1> StateP1VecD;
typedef Eigen::Matrix<real, SIMDIM + 1, SIMDIM + 1> StateP1MatD;

const string WHITESPACE = " \n\r\t\f\v";

#pragma omp declare reduction(+ : std::vector<StateVecD> : std::transform(                              \
        omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(),                                \
            [](StateVecD lhs, StateVecD rhs){return lhs + rhs;}                                         \
)) initializer(omp_priv = omp_orig)

#pragma omp declare reduction(+ : std::vector<real> : std::transform(                                   \
        omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(),                                \
            [](real lhs, real rhs){return lhs + rhs;}                                                   \
)) initializer(omp_priv = omp_orig)

#pragma omp declare reduction(+ : StateVecD : omp_out = omp_out + omp_in)                               \
    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(                                                                          \
        min:StateVecD : omp_out = omp_out.norm() < omp_in.norm() ? omp_in : omp_out                     \
) initializer(omp_priv = omp_orig)

#pragma omp declare reduction(                                                                          \
        max:StateVecD : omp_out = omp_out.norm() > omp_in.norm() ? omp_in : omp_out                     \
) initializer(omp_priv = omp_orig)

/*Define particle type indexes*/
enum partType
{
    BOUND = 0,
    PISTON,
    BUFFER,
    BACK,
    PIPE,
    FREE,
    OUTLET,
    LOST
};

/* Define shapes for boundary and fluid initialisation */
enum shapeType
{
    NONE = 0,
    BOX,
    CYLINDER,
    SPHERE,
    JET,
    CONVJET
};

// Shapes for new intitialisation
enum shape_type
{
    linePlane = 0,
    squareCube,
    circleSphere,
    cylinder,
    arcSection,
    coordDef,
    inletZone,
    hollow,
    solid
};

enum integrate_type
{
    newmark_beta = 0,
    runge_kutta
};

enum solve_type
{
    DBC = 0,
    pressure_G,
    ghost
};

enum inlet_vel_type
{
    fixedVel,
    dynamicVel
};

enum particle_order
{
    grid = 0,
    hcp
};

enum aero_force
{
    NoAero = 0,
    Gissler,
    InducedPressure,
    SkinFric
};

enum aero_source
{
    constVel = 0,
    meshInfl,
    VLMInfl
};

enum mesh_source
{
    TAU_CDF = 0,
    OpenFOAM
};
enum plane_axis
{
    no_offset = 0,
    x_axis,
    y_axis,
    z_axis
};

real const static default_val = 9999999.0;

inline bool check_vector(StateVecD const& v)
{
    return (
        v[0] != default_val && v[1] != default_val
#if SIMDIM == 3
        && v[2] != default_val
#endif
    );
}

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

// PState PartState;

inline real random(int const& interval) { return real(rand() % interval - interval / 2) * MERROR; }

// Structure for mapping
struct OutputVariable
{
    OutputVariable(string name, size_t type, bool write, bool vec)
        : output_name(name), data_type(type), write(write), is_vector(vec){};
    string output_name;
    size_t data_type;
    bool write, is_vector;
};

typedef std::map<string, OutputVariable> OutputMap;

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

#if SIMDIM == 3
        offset_axis = 0;
#else
        offset_axis = 2; /* Default to XZ plane in 2D */
#endif

/*Gravity Vector*/
#if SIMDIM == 3
        grav = StateVecD(0.0, 0.0, -9.81);
#else
        grav = StateVecD(0.0, -9.81);
#endif
    }

    uint numThreads; /* Number of threads that the simulation will use. */

    /* File input parameters */
    string input_file = "";       /* Input settings file. */
    string output_prefix = "";    /* Output file prefix before _boundary or _fluid. */
    string input_fluid_file = ""; /* Fluid definition settings file. */
    string input_bound_file = ""; /* Boundary definition settings file. */
    string restart_prefix = "";   /* Restart file prefix before _boundary or _fluid. */
    string vlm_file = "";         /* Vortex Lattice Method (VLM) settings file. Can be input_file. */
    uint mesh_source = TAU_CDF;   /* Using TAU or OpenFOAM for mesh input */

    /* OpenFOAM files */
    string foam_dir = "";          /* OpenFOAM project root directory. */
    string foam_sol = "";          /* OpenFOAM time folder to use. */
    bool foam_is_binary = false;   /* OpenFOAM files are stored in ASCII or binary format. */
    bool foam_buoyant_sim = false; /* OpenFOAM solution is using a buoyant sim or not. */
    bool foam_is_incomp = false;   /* OpenFOAM solution is foam_is_incompressible or not. */
    int foam_label_size = 32;      /* OpenFOAM label size for binary files. */
    int foam_scalar_size = 64;     /* OpenFOAM scalar size for binary files. */

    /* TAU files */
    string tau_mesh = "";      /* TAU mesh file name. */
    string tau_bmap = "";      /* TAU boundary map definition file name. */
    string tau_sol = "";       /* TAU solution file name. */
    int offset_axis;           /* TAU offset axis for two-dimensional  */
    real angle_alpha = 0.0;    /* TAU angle of attack. */
    vector<int> tau_markers;   /* TAU boundary surface tau_markers. */
    vector<string> tau_bnames; /* TAU boundary surface names  */
    vector<int> tau_bwrite;    /* TAU if surface wants data writing */

    /* TECIO File outputs */
    void* output_bound_file = NULL; /* TECIO SPH boundary file handle. */
    void* output_fluid_file = NULL; /* TECIO SPH fluid file handle. */
    std::vector<int32_t> var_types; /* TECIO data types for variables being output */
    string var_names;               /* TECIO variable name string output */
    string output_names = "";       /* Names of variables to output */
    OutputMap output_variables;     /* Internal map of variables storing which to be output */

    /* H5Part file outputs */
    int64_t h5part_bound_file = -1; /* HDF5 h5part SPH boundary file handle */
    int64_t h5part_fluid_file = -1; /* HDF5 h5part SPH fluid file handle */
    bool write_tecio = true;        /* Write Tecplot file (this or write_h5part needs to be true) */
    bool write_h5part = false;      /* Write H5part file (this or write_tecio needs to be true) */

    ofstream ascii_surface_file;       /* ACSII surface impact file */
    void* tecio_surface_handle = NULL; /* TECIO Surface impact file */

    /* Output type */
    uint out_encoding = 1;           /* ASCII or binary output*/
    uint single_file = 1;            /* Write timesteps to a single file or a file for each. */
    uint restart_header_written = 0; /* Whether the restart header has been written yet */
    uint restart = 0;                /* If starting from existing solution */
    real scale = 1.0;                /* Simulation scale */

    /* SPHPart counts */
    size_t n_bound_blocks = 0;   /* Number of boundary blocks */
    size_t n_fluid_blocks = 0;   /* Number of fluid blocks */
    size_t part_id = 0;          /* Track the particle ID separately */
    size_t total_points = 0;     /* Total count */
    size_t fluid_points = 0;     /* Fluid count */
    size_t bound_points = 0;     /* Boundary count */
    size_t max_points = 9999999; /* End count */
    size_t delete_count = 0;     /* Number of deleted particles */
    size_t internal_count = 0;   /* Number of internal particles */

    /* Particle size value */
    real particle_step = -1.0;    /* Particle step for creating particle blocks. Must be specified. */
    real bound_step_factor = 1.0; /* Boundary spacing factor of dx */
    real dx = -1.0;               /* Initial spacings for particles and boundary. Must be specified. */

    /* Geometry parameters */
    StateVecD offset_vec = StateVecD::Zero(); /* Global offset coordinate */
    int Asource = constVel;                   /* Source of aerodynamic solution */
    int init_hydro_pressure = 0;              /* Initialise fluid with hydrostatic pressure? */
    int use_global_gas_law = 1; /* Whether to use block specific gas laws (Currently not active) */
    real hydro_height = -1;     /* Hydrostatic height to initialise using */

    /* Integration parameters */
    string solver_name = "";         /* Name of solver */
    uint solver_type = newmark_beta; /* Use Runge-Kutta or Newmark-Beta. */
    uint max_subits = 20;            /* Max number of sub-iterations. */
    uint max_frames = -1;            /* Max number of frames to output. Must be specified. */
    uint current_frame = 0;          /* Current frame number. */
    uint bound_solver = 0;           /* Use boundary pressure or ghost particles. */
    uint n_stable = 0;               /* Count for number of overly stable timesteps to alter CFL */
    uint n_stable_limit = 10;        /* Limit before changing CFL */
    uint n_unstable = 0;             /* Count for number of unstable timestep to alter CFL */
    uint n_unstable_limit = 3;       /* Limit before changing CFL */
    real subits_factor = 0.333;   /* Factor * max subits, under which is considered overly stable CFL */
    double cfl = 1.0;             /* CFL criterion number. */
    double cfl_step = 0.05;       /* CFL step to perform if unstable. */
    double cfl_max = 2.0;         /* Maximum CFL for the simulation. */
    double cfl_min = 0.1;         /* Minimum CFL for the simulation. */
    double current_time = 0.0;    /* Simulation time. */
    double last_frame_time = 0.0; /* Last written frame time. */
    real frame_time_interval;     /* Frame time interval. */
    real min_residual = -7.0;     /* Minimum solver residual. */
    real delta_t = 2e-10;         /* Integration timestep. */
    real delta_t_max = 1.0;       /* Maximum timestep. */
    real delta_t_min = 0.0;       /* Minimum timestep. */
    real nb_beta = 0.25;          /* Newmark-Beta beta parameter. */
    real nb_gamma = 0.5;          /* Newmark-Beta gamma parameter. */
    real max_shift_vel = 9999999; /* Maximum shifting velocity. */
    StateVecD grav;               /* Gravity vector (assumed to be on z-axis). */

    real restart_tol = 0.001; /* Tolerance on buffer particles fitting. */

    /* Particle tracking settings */
    int using_ipt = 0;        /* Transition to IPT after reaching a point. */
    int ipt_eq_order = 2;     /* Order of equation to use for IPT. */
    real max_x_sph = 9999999; /* Transition x-coordinate to convert to IPT. */
    real max_x = 9999999; /* Transition distance for sph (naive assumption) and termination for IPT */
    size_t ipt_n_success = 0;   /* Number of successful IPT particles. */
    size_t ipt_n_failed = 0;    /* Number of failed IPT particles. */
    real ipt_diam = 0.0;        /* Diameter of IPT particle. */
    real ipt_area = 0.0;        /* Cross-sectional area of IPT particle. */
    real relax = 0.6;           /* Relaxation factor for IPT iteration. */
    real n_relax = 5;           /* Number of relaxation steps to try. */
    uint cells_out = 0;         /* IPT output cells traversed by particle. */
    uint streak_out = 1;        /* IPT output streamlines. */
    uint part_out = 0;          /* Whether to output for particle tracking. */
    FILE* cell_file = NULL;     /* ACSII file handle for IPT cells intersected. */
    FILE* streak_file = NULL;   /* ASCII file handle for IPT streamlines. */
    FILE* part_file = NULL;     /* ASCII file handle for IPT particle scatter data. */
    void* cell_handle = NULL;   /* TECIO file handle for IPT cells intersected. */
    void* streak_handle = NULL; /* TECIO file handle for IPT streamlines. */
    void* part_handle = NULL;   /* TECIO file handle for IPT particle scatter data. */
};

/*Fluid and smoothing parameters*/
struct FLUID
{
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    real H = -1.0;          /* Support radius */
    real H_sq = -1.0;       /* Support radius squared */
    real sr = -1.0;         /* Search radius (support radius * factor)*/
    real H_fac = 2.0;       /* Search radius factor */
    real W_dx = -1.0;       /* Kernel value at the initial particle spacing distance.*/
    real rho_rest = 1000.0; /* Resting fluid density */
    real rho_pipe = 1000.0; /* Pipe starting density */
    real press_pipe = 0.0;  /* Starting pressure in pipe*/
    real press_back = 0.0;  /* Background pressure */

    /* Allow a large variation by default */
    /* making it essentially unbounded other */
    /* than truly unphysical pressures */
    real rho_max = 1500.0;   /* Maximum density value */
    real rho_min = 500.0;    /* Minimum density value */
    real rho_var = 50.0;     /* Maximum density variation allowed (in %) */
    real rho_max_iter = 1.0; /* Maximum value to reduce the timestep */

    real sim_mass = -1.0;     /* SPH fluid particle mass */
    real bnd_mass = -1.0;     /* SPH boundary mass */
    real W_correc = -1.0;     /* Smoothing Kernel Correction*/
    real visc_alpha = 0.1;    /* Artificial viscosity factor */
    real speed_sound = 300.0; /* Speed of sound */
    real mu = 8.94e-4;        /* Dynamic viscosity */
    real nu = mu / rho_rest;  /* Kinematic viscosity */
    real sig = 0.0708;        /* Surface tension coefficient */
    real gam = 7.0;           /* Tait equation gamma power*/
    real B = -1.0;            /* Tait equation speed of sound factor */

    real contangb = 150.0; /* Boundary contact angle */
    real res_vel = 0.0;    /* Reservoir velocity */
    // real front, height, height0; /* Dam Break validation parameters*/

    /* delta SPH terms */
    real dsph_delta = 0.1; /* delta-SPH contribution */
    real dsph_cont = -1.0; /* delta-SPH continuity constant term */
    real dsph_mom = -1.0;  /* delta-SPH momentum constant term (dsph_cont * rho0) */
};

/*Aerodynamic Properties*/
struct AERO
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    void GetYcoef(const FLUID& fvar, const real diam)
    {
#if SIMDIM == 3
        L = diam * std::cbrt(3.0 / (4.0 * M_PI));
        A_sphere = M_PI * L * L;
        A_plate = diam * diam;
#else
        L = diam / sqrt(M_PI);
        A_sphere = 2 * L;
        A_plate = diam;
#endif

        td = (2.0 * fvar.rho_rest * pow(L, SIMDIM - 1)) / (Cd * fvar.mu);

        omega = sqrt((Ck * fvar.sig) / (fvar.rho_rest * pow(L, SIMDIM)) - 1.0 / pow(td, 2.0));

        tmax = -2.0 * (atan(sqrt(pow(td * omega, 2.0) + 1) + td * omega) - M_PI) / omega;

        Cdef = 1.0 - exp(-tmax / td) * (cos(omega * tmax) + 1 / (omega * td) * sin(omega * tmax));
        ycoef = 0.5 * Cdef * (Cf / (Ck * Cb)) * (rho_g * L) / fvar.sig;

        // cout << ycoef << "  " << Cdef << "  " << tmax << "  " << endl;
    }

    StateVecD v_inf = StateVecD::Zero(); /* Freestream velocity */

    /*Gissler Parameters*/
    real L = -1.0; /* Aerodynamic length */
    real td = -1.0;
    real omega = -1.0;
    real tmax = -1.0;
    real ycoef = -1.0;
    real Cf = 1.0 / 3.0;
    real Ck = 8.0;
    real Cd = 5.0;
    real Cb = 0.5;
    real Cdef = -1.0;

    real A_sphere = -1.0; /* Area of a sphere/cylinder */
    real A_plate = -1.0;  /* Area of a plate*/

    /* Gas Properties*/
    real v_ref = 0.0;                      /* Reference velocity */
    real p_ref = 101353.0;                 /* Reference pressure */
    real rho_g = 1.29251;                  /* Gas density */
    real mu_g = 1.716e-5;                  /* Gas dynamic viscosity */
    real mass_g = -1.0;                    /* Gas particle mass */
    real temp_g = 298.0;                   /* Gas temperature */
    real R_g = 287.0;                      /* Specific gas constant */
    real gamma = 1.403;                    /* Ratio of specific heats */
    real sos = sqrt(temp_g * R_g * gamma); /* Speed of sound */
    real i_sos_sq = 1.0 / (sos * sos);     /* Inverse speed of sound squared */
    real M_ref = -1.0;                     /* Reference Mach */

    real lam_cutoff = 0.75;                         /* Lambda value cutoff */
    real interp_fac = 2.0;                          /* Factor for the interpolation value */
    real i_interp_fac = 0.5;                        /* Inverse factor for the interpolation */
    real n_full = SIMDIM == 2 ? 30 : 200;           /* 2/3ds full neighbourhood count */
    real i_n_full = 1.0 / n_full;                   /* inverse of n_full */
    real i_n_count = 1.0 / (i_interp_fac * n_full); /* inverse nfull * factor */

    string aero_case;    /* String input to define aero force type */
    int acase = NoAero;  /* Aerodynamic force case */
    int use_lam = 1;     /* Use lambda value for aero interpolation */
    int use_TAB_def = 0; /* Use TAB deformation estimate  */
    int use_dx = 0;      /* Use particle spacing or support radius for aero force. */
};

struct MESH
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /*Standard contructor*/
    MESH() : nPnts(0), nElem(0), nFace(0), nSurf(0), scale(1.0), maxlength(0.0), minlength(9999999.9) {}

    inline size_t size() { return nElem; }

    inline void alloc(size_t nPnts_, size_t nElem_, size_t nFace_, size_t nSurf_)
    {
        nPnts = nPnts_;
        nElem = nElem_;
        nFace = nFace_;
        nSurf = nSurf_;
        verts = vector<StateVecD>(nPnts);
        faces.reserve(nFace * 2);
        leftright.reserve(nFace);

        // elems = vector<vector<size_t>>(nElem);
        cCentre = vector<StateVecD>(nElem);
        cFaces = vector<vector<size_t>>(nElem);

        cVel = vector<StateVecD>(nElem);
        cP = vector<real>(nElem);
        cRho = vector<real>(nElem);
    }

    /*Zone info*/
    size_t nPnts = 0; /* Number of points */
    size_t nElem = 0; /* Number of elements */
    size_t nFace = 0; /* Number of faces (internal and boundary) */
    size_t nSurf = 0; /* Number of boundary faces */
    real scale;       /* Coordinate scale */

    real maxlength; /* Max cell edge length */
    real minlength; /* Min cell edge length */

    /*Point based data*/
    vector<StateVecD> verts; /* Point coordinates */

    /*Face based data*/
    vector<vector<size_t>> faces;            /* Face indexes */
    vector<std::pair<int, int>> leftright;   /* Face cell left/right indexes */
    vector<std::pair<size_t, int>> smarkers; /* Surface indicators; cell, and bmap ID */

    /*Cell based data*/
    // vector<vector<size_t>> elems;
    vector<StateVecD> cCentre;     /* Cell centriod coordinates */
    vector<vector<size_t>> cFaces; /* Cell face IDs  */

    /*Solution vectors*/
    vector<StateVecD> cVel; /* Cell velocities */
    vector<real> cCp;       /* Cell coefficient of pressure */
    vector<real> cP;        /* Cell pressure */
    vector<real> cRho;      /* Cell density */
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

        face_IDs = in.face_IDs;
        face_count = vector<uint>(face_IDs.size(), 0);
    }

    /* Marker data */
    string name;
    int marker;
    int output;

    unsigned marker_count;
    real marker_beta; /* Collision efficiency (same as collection eff) */
    real marker_area;
    vector<size_t> marker_pIDs;   /* ID for particles that hit the surface */
    vector<size_t> impacted_face; /* ID of the face the particle hit */
    vector<StateVecD> end_pos;
    vector<StateVecD> start_pos;

    /* Face data */
    vector<size_t> face_IDs;
    vector<uint> face_count;
    vector<real> face_beta;
    vector<real> face_area;
};

/*SPHPart data class*/
struct SPHPart
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SPHPart(
        StateVecD const& X, StateVecD const& Vi, real const Rhoi, real const Mi, real const press,
        int const bound, uint const p_id
    )
    {
        part_id = p_id;
        faceID = 0;
        cellID = -1;
        b = bound;
        surf = 0;
        surfzone = 0;
        ipt_n_failed = 0;

        xi = X;
        v = Vi;
        acc = StateVecD::Zero();
        Af = StateVecD::Zero();
        aVisc = StateVecD::Zero();
        norm = StateVecD::Zero();
        Rrho = 0.0;
        rho = Rhoi;
        p = press;
        m = Mi;
        curve = 0.0;
        norm_curve = 0.0;
        woccl = 0.0;
        pDist = 0.0;
        deltaD = 0.0;

        cellV = StateVecD::Zero();
        cellP = 0.0;
        cellRho = 0.0;
        internal = 0;

        L = StateMatD::Zero();
        gradRho = StateVecD::Zero();
        norm = StateVecD::Zero();
        colourG = 0.0;
        colour = 0.0;
        lam = 0.0;
        lam_nb = 0.0;
        kernsum = 0.0;
        bNorm = StateVecD::Zero();
        y = 0.0;
        vPert = StateVecD::Zero();
    }

    /*To add particles dynamically for fictitious particles*/
    SPHPart(StateVecD const& X, SPHPart const& pj, int const bound, size_t const p_id)
    {
        part_id = p_id;
        faceID = 0;
        cellID = -1;
        b = bound;
        surf = 0;
        surfzone = 0;
        ipt_n_failed = 0;

        xi = X;
        v = pj.v;
        acc = StateVecD::Zero();
        Af = StateVecD::Zero();
        aVisc = StateVecD::Zero();
        norm = StateVecD::Zero();
        Rrho = 0.0;
        rho = pj.rho;
        p = pj.p;
        m = pj.m;
        curve = 0.0;
        norm_curve = 0.0;
        woccl = 0.0;
        pDist = 0.0;
        deltaD = 0.0;

        cellV = StateVecD::Zero();
        cellP = pj.cellP;
        cellRho = pj.cellRho;
        internal = 0;

        L = StateMatD::Zero();
        gradRho = StateVecD::Zero();
        norm = StateVecD::Zero();
        colourG = 0.0;
        colour = 0.0;
        lam = 0.0;
        lam_nb = 0.0;
        kernsum = 0.0;
        bNorm = StateVecD::Zero();
        y = 0.0;
        vPert = StateVecD::Zero();
    }

    SPHPart(){};

    inline int size() const { /* Return size of xi vector*/ return xi.size(); }

    inline real operator[](int a) const { /* Return index of xi vector*/ return xi[a]; }

    size_t part_id = 0;    /* Particle ID */
    size_t faceID = 0;     /* Mesh face ID */
    long int cellID = -1;  /* Cell ID. -1 is not in a cell. */
    uint b = 0;            /* Status flag. See PartState above for possible options */
    uint surf = 0;         /* Is a particle a surface? 1 = yes, 0 = no */
    uint surfzone = 0;     /* Is particle in the surface area? */
    uint ipt_n_failed = 0; /* How many times has the containment query failed */

    StateVecD xi = StateVecD::Zero();    /* Particle position */
    StateVecD v = StateVecD::Zero();     /* Particle velocity */
    StateVecD acc = StateVecD::Zero();   /* Particle acceleration */
    StateVecD Af = StateVecD::Zero();    /* Particle aerodynamic force. */
    StateVecD aVisc = StateVecD::Zero(); /* Artificial viscosity contribution */
    real Rrho = 0.0;                     /* Density gradient */
    real rho = 0.0;                      /* Density */
    real p = 0.0;                        /* Pressure */
    real m = 0.0;                        /* Mass */
    real curve = 0.0;                    /* Surface curvature */
    real norm_curve = 0.0;               /* Normalised curvature */
    real woccl = 0.0;                    /* Occlusion factor */
    real pDist = 0.0;                    /* Particle distance */
    real deltaD = 0.0;                   /* delta-SPH contribution */
    StateVecD cellV = StateVecD::Zero(); /* Cell velocity */
    real cellP = 0.0;                    /* Cell pressure */
    real cellRho = 0.0;                  /* Cell density */
    uint internal = 0;                   /* Is particle across an internal mesh boundary? */

    StateMatD L = StateMatD::Zero();       /* L gradient renormalisation matrix */
    StateVecD gradRho = StateVecD::Zero(); /* Density Gradient Jacobian */
    StateVecD norm = StateVecD::Zero();    /* Surface normal */

    real colourG = 0.0; /* Colour gradient. */
    real colour = 0.0;  /* Kernel sum with volume considered. */
    real lam = 0.0;     /* Eigenvalues with all particles considered */
    real lam_nb = 0.0;  /* Eigenvalues without boundary particles */
    real kernsum = 0.0; /* Summation of the kernel */

    /* Mesh surface repulsion */
    StateVecD bNorm = StateVecD::Zero(); /* Boundary normal value. */
    real y = 0.0;                        /* Distance from boundary */

    StateVecD vPert = StateVecD::Zero(); /* ALE velocity perturbation */
};

struct IPTPart
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    IPTPart()
    {
        part_id = 0;
        going = 1;
        nIters = 0;
        failed = 0;
        t = 0;
        dt = 0;

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
        part_id = 0;
        going = 1;
        nIters = 0;
        failed = 0;
        t = 0;
        dt = 0;

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
        A = M_PI * d_ * d_ / 4.0;
    }

    IPTPart(StateVecD const xi_, IPTPart const& pi_, size_t const& pID_)
    {
        part_id = pID_;
        going = 1;
        nIters = 0;
        failed = 0;
        t = 0;
        dt = 0;

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
        part_id = pi.part_id;
        going = 1;
        nIters = 0;
        failed = 0;
        t = time;
        dt = 0;

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
        part_id = pi.part_id;
        going = 1;
        nIters = 0;
        failed = 0;
        t = time;
        dt = 0;

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

    size_t part_id;
    uint going;  /* Is particle still being integrated */
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
        particle_order = 0;
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
        particle_order = 0;
        bound_solver = 0;
        no_slip = 0;
        block_type = 0;
    }

    string name;                     /* Block name */
    std::pair<size_t, size_t> index; /* boundary left and right */

    StateVecD delete_norm; /* Deletion plane */
    real delconst;         /* Deletion constant to adjust plane */
    StateVecD insert_norm; /* Insertion plane */
    real insconst;         /* Insertion constant to adjust plane */
    StateVecD aero_norm;   /* Aerodynamic beginning plane */
    real aeroconst;        /* Aerodynamic constant to adjust plane */

    vector<size_t> back;           /* Particles at the back of the pipe */
    vector<vector<size_t>> buffer; /* ID of particles inside the buffer zone */

    vector<real> times;     /* Times for boundary motion */
    vector<StateVecD> vels; /* Velocities for boundary motion */

    size_t nTimes;            /* Number of times for boundary motion */
    int fixed_vel_or_dynamic; /* Used for inlets, but also for boundaries to identify if moving or static
                               */
    int particle_order;       /* Use regular grid or HCP packing */
    int no_slip;              /* Is wall no slip */
    int bound_solver;         /* Solver to use for boundaries */
    int block_type;           /* Block type from definition files. */
};

typedef std::vector<bound_block> LIMITS;
typedef std::vector<SPHPart> SPHState;
typedef std::vector<IPTPart> IPTState;
typedef std::vector<SURF> SURFS;

/* Neighbour search tree containers */
typedef nanoflann::ResultItem<size_t, real> neighbour_index;
typedef std::vector<std::vector<neighbour_index>> OUTL;
typedef std::vector<std::vector<size_t>> celll;
typedef KDTreeVectorOfVectorsAdaptor<SPHState, real, SIMDIM, nanoflann::metric_L2_Simple, size_t>
    Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<
    std::vector<StateVecD>, real, SIMDIM, nanoflann::metric_L2_Simple, size_t>
    Vec_Tree;

// Default to no approximate neighbours, and unsorted list.
nanoflann::SearchParameters const flann_params(0, false);

#endif /* VAR_H */
