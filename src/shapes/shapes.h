/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef SHAPES_H
#define SHAPES_H

#include "../Var.h"
#include <random>

struct shape_block
{ /* Basically everything needs to be defined. Provide defaults to check against */
    /* Boolean parameters */
    bool write_data = false;
    bool no_slip = false;
    bool particle_order = false; /* Use grid of HCP ordering */

    /* String parameters */
    std::string name;               /* Name to write in files for particle block */
    std::string shape;              /* Predominant shape of the block. */
    std::string subshape;           // For inlets (square or circle)
    std::string filename;           /* Name of the file used to read the block */
    std::string position_filename;  /* File for defined motion of this boundary*/
    std::string solver_name;        /* Name of the boundary solver to use */
    std::string particle_order_str; /* Grid or HCP string */
    std::string inlet_bc_type;      /* Inlet velocity definition */

    /* Integer parameters */
    int block_index = 0;           /* Index of the block in the file */
    int bound_type = -1;           /* Integer value for what boundary type it is */
    int sub_bound_type = -1;       /* Integer value for sub-type */
    int fixed_vel_or_dynamic = 0;  /* Integer value for inlet velocity type */
    int bound_solver = pressure_G; /* Integer value for boundary solver type */

    int ni = 0, nj = 0, nk = 0; /* counters in i,j,k */

    size_t npts = 0;   /* Number of points */
    size_t ntimes = 0; /* Number of timestamps */

    /* Float parameters */
    real insconst = default_val;
    real delconst = default_val;
    real aeroconst = default_val;

    real dx = -1;                 /* Particle spacing */
    real radius = -1;             /* Radius */
    real arc_start = -1;          /* Arc start */
    real arc_end = -1;            /* Arc end */
    real arclength = default_val; /* Arc length in degrees */
    real thickness = -1;          /* How thick is the wall? */
    real length = -1;             /* Archway length for 3D */
    real sstraight = 0;           /* Straight to the start point */
    real estraight = 0;           /* Straight from the end point */

    // Starting physical properties
    real vmag = 0;          /* Velocity magnitude in jet direction */
    real press = 0;         /* Starting pressure */
    real dens = 1000;       /* Starting density (derived from pressure or specified?) */
    real mass = -1;         /* Starting mass (derived from spacing and density) */
    real renorm_vol = -1;   /* Volume to renormalise mass using */
    real nu = -1;           /* Kinematic viscosity */
    real rho0 = 1000;       /* Resting density */
    real gamma = 7;         /* Cole gamma value */
    real speedOfSound = -1; /* Speed of sound */
    real backgroundP = 0;   /* Background pressure */

    /* Vector parameters */
    StateVecD stretch = StateVecD::Constant(1.0); /* Stretching coefficient to test tension */

    StateVecD insert_norm = StateVecD::UnitX();
    StateVecD delete_norm = StateVecD::UnitX();
    StateVecD aero_norm = StateVecD::UnitX();

    // Inlet variables
    StateVecD normal = StateVecD::UnitX();    /* Normal vector to rotate */
    StateVecD angles = StateVecD::Zero();     /* Angles for rotation */
    StateMatD rotmat = StateMatD::Identity(); /* Rotation matrix */

    StateVecD start = StateVecD::Constant(default_val); /* Start coordinate */
    StateVecD end = StateVecD::Constant(default_val);   /* End coordinate */
    StateVecD right = StateVecD::Constant(default_val); /* Right coordinate */

    StateVecD mid = StateVecD::Constant(default_val);    /* Arc midpoint */
    StateVecD centre = StateVecD::Constant(default_val); /* Arc/circle centrepoint */

    StateVecD static_vel = StateVecD::Zero(); /* Static velocity for boundaries */
    StateVecD vel = StateVecD::Zero();        /* Starting velocity */

    /* Array parameters */
    std::vector<size_t> back;
    std::vector<std::vector<size_t>> buffer;

    std::vector<uint> bc;       /* Boundary condition of point */
    std::vector<int> intersect; /* For fluid intersection with fibres */

    std::vector<real> times;     /* Timestamps */
    std::vector<StateVecD> pos;  /* Positions at each time*/
    std::vector<StateVecD> vels; /* Velocities at times-1 */

    std::vector<StateVecD> coords; /* Coordinates */

    void check_input(SIM const& svar, FLUID const& fvar, real& globalspacing, int& fault);
    void generate_points(real const& globalspacing);
};

struct Shapes
{
    std::vector<shape_block> block;
    size_t totPts = 0;  /* Total boundary points */
    size_t nblocks = 0; /* Number of boundaries */

    inline void emplace_back()
    {
        block.emplace_back();
        nblocks++;
    }

    inline shape_block& back() { return block.back(); }
};

Shapes
read_shapes_JSON(std::string const& filename, SIM const& svar, FLUID const& fvar, real& globalspacing);

#endif
