#ifndef SHAPES_H
#define SHAPES_H

#include "../Var.h"
#include <random>


struct shape_block
{   /* Basically everything needs to be defined. Provide defaults to check against */
    shape_block() : bound_type(-1), sub_bound_type(-1), hcpl(0), fixed_vel_or_dynamic(0), 
        bound_solver(1), no_slip(0), npts(0),
        insert_norm(StateVecD::UnitX()), insconst(default_val),
        delete_norm(StateVecD::UnitX()), delconst(default_val),
        aero_norm(StateVecD::UnitX()), aeroconst(default_val),  
        ntimes(0), dx(-1), write_data(0),
        stretch(StateVecD::Constant(1.0)), ni(-1), nj(-1), nk(-1), 
        normal(StateVecD::UnitX()),
        angles(StateVecD::Zero()), rotmat(StateMatD::Identity()), 
        start(StateVecD::Constant(default_val)), end(StateVecD::Constant(default_val)), 
        right(StateVecD::Constant(default_val)), mid(StateVecD::Constant(default_val)), 
        centre(StateVecD::Constant(default_val)), radius(-1), arc_start(-1), arc_end(-1), arclength(default_val),
        thickness(-1), length(-1), sstraight(0), estraight(0), vel(StateVecD::Constant(0.0)), vmag(0.0), press(0.0), 
        dens(1000), mass(-1), renorm_vol(-1), nu(-1), rho0(1000), gamma(7), speedOfSound(-1), backgroundP(0) {}

    std::string name;
    std::string shape;
    std::string subshape; // For inlets (square or circle)
    std::string filename;
    std::string position_filename;  /* File for defined motion of this boundary*/
    std::string solver_name;

	int bound_type; /* What boundary type is it? */
    int sub_bound_type;
    int hcpl; /* Use grid of HCP ordering */
    int fixed_vel_or_dynamic;
    int bound_solver;
    int no_slip;
    
	size_t npts;			/* Number of points */
	std::vector<StateVecD> coords; /* Coordinates */
    std::vector<uint> bc; /* Boundary condition of point */
    std::vector<int> intersect; /* For fluid intersection with fibres */
    StateVecD insert_norm;
    real insconst;
    StateVecD delete_norm;
    real delconst;
    StateVecD aero_norm;
    real aeroconst;

    std::vector<size_t> back;
    std::vector<std::vector<size_t>> buffer;

    size_t ntimes;   /* Number of timestamps */
    std::vector<real> times; /* Timestamps */
    std::vector<StateVecD> pos; /* Positions at each time*/
    std::vector<StateVecD> vels; /* Velocities at times-1 */

	real dx;		/* Particle spacing */
    int write_data;

    StateVecD stretch;                 /* Stretching coefficient to test tension */

    int ni, nj, nk;       /* counters in i,j,k*/
    
    // Inlet variables
    StateVecD normal;   /* Normal vector to rotate */
    StateVecD angles;   /* Angles for rotation */
    StateMatD rotmat;   /* Rotation matrix */

	StateVecD start;	/* Start coordinate */
	StateVecD end;		/* End coordinate */
    StateVecD right;    /* Right coordinate */

	StateVecD mid;		/* Arc midpoint */
	StateVecD centre;   /* Arc/circle centrepoint */
	real radius;	    /* Radius */
    real arc_start;     /* Arc start */
    real arc_end;       /* Arc end */
    real arclength;     /* Arc length in degrees */
    real thickness;     /* How thick is the wall? */
    real length;        // Archway length for 3D
    real sstraight;     // Straight to the start point
    real estraight;     // Straight from the end point
	// int nthick;      /* Number of particles along the thickness */
    
    // Starting physical properties
    StateVecD static_vel; /* Static velocity for boundaries */
    StateVecD vel;      /* Starting velocity */
    real vmag;          /* Velocity magnitude in jet direction */
    real press;         /* Starting pressure */
    real dens;          /* Starting density (derived from pressure or specified?) */
    real mass;          /* Starting mass (derived from spacing and density) */
    real renorm_vol;    /* Volume to renormalise mass using */
    real nu;            /* Kinematic viscosity */
    /* Gas law properties */
    real rho0;          /* Resting density */
    real gamma;         /* Cole gamma value */
    real speedOfSound;  /* Speed of sound */
    real backgroundP;   /* Background pressure */
};

struct Shapes
{
	Shapes():totPts(0),nblocks(0){};

	std::vector<shape_block> block;
	size_t totPts; /* Total boundary points */
	size_t nblocks;	/* Number of boundaries */
};



#endif
