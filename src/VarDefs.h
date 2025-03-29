/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef VARDEFS_H
#define VARDEFS_H

#include <string.h>

#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"

/* Define Simulation Dimension */
#ifndef SIMDIM
#define SIMDIM 2
#endif

/* Define the default surface tension model */
#if !defined(HESF) && !defined(CSF) && !defined(PAIRWISE) && !defined(NOST)
#define PAIRWISE
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

const std::string WHITESPACE = " \n\r\t\f\v";

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

enum bound_solve_type
{
    DBC = 0,
    pressure_G,
    ghost
};

enum pressure_solver
{
    COLE_EOS = 0,
    ISO_EOS
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

enum output_encoding
{
    ascii = 0,
    binary
};

real const static default_val = 9999999.0;
int const static c_no_face = -1;
int const static c_no_cell = -3;

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#endif