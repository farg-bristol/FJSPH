/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VAR_H
#define VAR_H

#include <vector>
#include "../Eigen/Core"
#include "../Eigen/StdVector"
#include "../NanoFLANN/nanoflann.hpp"
#include "../NanoFLANN/utils.h"
#include "../NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"


// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

/* Define Simulation dimension, and data type. */
/* Want to have long double at some point,     */
/* but neighbour search won't have it...       */
constexpr int simDim = 2;
typedef double ldouble;

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<ldouble,simDim,1> StateVecD;
typedef Eigen::Matrix<int,simDim,1> StateVecI;
//typedef Eigen::Vector2d StateVecD;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<ldouble, simDim+1,1> DensVecD;
typedef Eigen::Matrix<ldouble, simDim+1, simDim+1> DensMatD;

/*Simulation parameters*/
typedef struct SIM {
	StateVecI xyPART; 						/*Starting sim particles in x and y box*/
	unsigned int simPts,bndPts,totPts;	    /*Simulation particles, Boundary particles, total particles*/
	unsigned int nrefresh; 					/*Crossflow refresh particle number*/
	unsigned int nmax;                      /*Max add-particle calls*/
	unsigned int outframe;	                /*Terminal output frame interval*/
	unsigned int addcount;					/*Current Number of add-particle calls*/
	unsigned int aircount;					/*Ghost particle count*/
	double Pstep,Bstep;						/*Initial spacings for particles and boundary*/
	StateVecD Box;							/*Box dimensions*/
	StateVecD Start; 						/*Sim box bottom left coordinate*/
	unsigned int subits;                    /*Max number of sub-iterations*/
	unsigned int Nframe; 			        /*Max number of frames to output*/
	double dt, t, framet;					/*Timestep, Simulation + frame times*/
	double beta,gamma;						/*Newmark-Beta Parameters*/
	double maxmu;                           /*Maximum viscosity component (CFL)*/
	int Bcase, Bclosed;						/*What boundary shape to take*/
	int outform;                            /*Output type. Fluid properties or Research.*/
	int frameout;                           /**/
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	ldouble H, HSQ, sr; 					/*Support Radius, SR squared, Search radius*/
	ldouble rho0; 						/*Resting Fluid density*/
	ldouble Simmass, Boundmass;			/*Particle and boundary masses*/
	ldouble correc;						/*Smoothing Kernel Correction*/
	ldouble alpha,eps,Cs,mu;			/*}*/
	ldouble sig;							/*} Fluid properties*/
	ldouble gam, B; 						/*}*/
	ldouble contangb;					/*Boundary contact angle*/

	//double front, height, height0;		/*Dam Break validation parameters*/

}FLUID;

/*Crossflow parameters*/
typedef struct CROSS
{
	int acase;	                        /*Aerodynamic force case*/
	StateVecD vJet, vInf;               /*Crossflow Parameters: Jet + Freestream velocity*/
	ldouble rhog;
	ldouble mug;
	ldouble Acorrect;					/*Correction factor for aero force*/
	ldouble a;                          /*Tuning parameters*/
	ldouble b;                          /*Tuning parameters*/
	ldouble h1;                          /*Tuning parameters*/
	ldouble h2;                          /*Tuning parameters*/
	StateVecD vortexPos;
}CROSS;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(StateVecD X, StateVecD Vi, StateVecD Fi,
			ldouble Rhoi, ldouble Mi, int bound)
		{
			xi = X;	v = Vi; f = Fi;
			V(0.0,0.0);	Sf(0.0,0.0); Af(0.0,0.0);
			rho = Rhoi;
			Rrho = 0;
			m = Mi;
			b = bound;
			theta = 0.0;
		}


		int size() const
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		double operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		StateVecD xi, v, V, f, Sf, Af;
		ldouble rho, p, Rrho, m, theta;
		int b; //What state is a particle. Boundary, forced particle or unforced
}Particle;
typedef std::vector<Particle> State;

/* Neighbour search tree containers */
typedef std::vector<std::vector<size_t>> outl;
typedef KDTreeVectorOfVectorsAdaptor<State, ldouble> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,ldouble> Temp_Tree;


#endif /* VAR_H */
