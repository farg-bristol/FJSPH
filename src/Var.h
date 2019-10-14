/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VAR_H
#define VAR_H

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 2
#endif

#ifndef NTHREADS
#define NTHREADS 4
#endif

#if SIMDIM == 3
	#include "VLM.h"
#endif

#include <vector>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include "NanoFLANN/nanoflann.hpp"
#include "NanoFLANN/utils.h"
#include "NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using std::vector;

/* Define data type. */
/* Want to have long double at some point,     */
/* but neighbour search won't have it...       */
typedef double ldouble;
typedef unsigned int uint;

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<ldouble,SIMDIM,1> StateVecD;
typedef Eigen::Matrix<int,SIMDIM,1> StateVecI;
typedef Eigen::Matrix<ldouble,SIMDIM,SIMDIM> RotMat;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<ldouble, SIMDIM+1,1> DensVecD;
typedef Eigen::Matrix<ldouble, SIMDIM+1, SIMDIM+1> DensMatD;

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

/*Simulation parameters*/
typedef struct SIM {
	StateVecI xyPART; 				/*Starting sim particles in x and y box*/
	uint simPts,bndPts,totPts;	    /*Simulation particles, Boundary particles, total particles*/
	uint finPts;					/*How many points there will be at simulation end*/
	uint psnPts;					/*Piston Points*/
	uint nrefresh; 					/*Crossflow refresh particle number*/
	uint nmax;                      /*Max add-particle calls*/
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
	uint outform, boutform;         /*Output type. Fluid properties or Research.*/
	uint frameout;                  /**/
	uint framecount;

	std::string infolder, outfolder;
	std::string meshfile;			
	#if SIMDIM == 3
		VLM vortex;
	#endif
	StateVecD Force;					/*Total Force*/
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	ldouble H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
	ldouble rho0; 					/*Resting Fluid density*/
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

typedef struct MESH
{
	void reserve(const uint nV, const uint nC, const uint nCverts, const uint nfaces, const uint nFverts)
	{
		verts = std::vector<StateVecD>(nV);
		pVel = std::vector<StateVecD>(nV);
		pointCp = std::vector<double>(nV);
		pointP = std::vector<double>(nV);
		pointRho = std::vector<double>(nV);

		elems = std::vector<std::vector<uint>>(nC,std::vector<uint>(nCverts));
		#if SIMDIM == 2
			cVerts = std::vector<std::vector<StateVecD>>(nC,std::vector<StateVecD>(nCverts));
			if(nfaces != 0 || nFverts !=0)
			{
				std::cout << "Some 3D data has been initialised.\n"
				<< "Please check the simulation dimensions" << std::endl;
			}
		#else /*Don't ask...*/
			cFaces = std::vector<std::vector<std::vector<StateVecD>>>
				(nC,std::vector<std::vector<StateVecD>>(nfaces,std::vector<StateVecD>(nFverts)));
		#endif
		cVel = std::vector<StateVecD>(nC);		
		cellCp = std::vector<double>(nC);
		cellP = std::vector<double>(nC);
		cellRho = std::vector<double>(nC);
		// cNeighb = std::vector<std::vector<uint>>(nC,std::vector<uint>());
		cNeighb.reserve(nC);

		numPoint = nV;
		numElem = nC;
		nFaces = nfaces;
	}

	std::string zone;
	uint numPoint, numElem, nFaces;
	/*Point based data*/
	std::vector<StateVecD> verts;
	std::vector<StateVecD> pVel;
	std::vector<double> pointCp;
	std::vector<double> pointP;
	std::vector<double> pointRho;

	/*Cell based data*/
	std::vector<std::vector<uint>> elems;
	#if SIMDIM == 2
		std::vector<std::vector<StateVecD>> cVerts;
	#endif
	/*Yep... a quad layered vector... I don't like it any more than you*/
	#if SIMDIM == 3
		std::vector<std::vector<std::vector<StateVecD>>> cFaces; 
	#endif
	std::vector<StateVecD> cVel;
	std::vector<std::vector<uint>> cNeighb;
	std::vector<double> cellCp;
	std::vector<double> cellP;
	std::vector<double> cellRho;
}MESH;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(const StateVecD X, const StateVecD& Vi, const ldouble Rhoi, const ldouble Mi, 
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

		StateVecD xi, v, f, Sf, Af, normal;
		ldouble rho, p, Rrho, m, theta;
		uint b; //What state is a particle. Boundary, forced particle or unforced
		StateVecD cellV;
		uint partID, cellID;
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

typedef std::vector<Particle> State;

/* Neighbour search tree containers */
typedef std::vector<std::vector<uint>> outl;
typedef KDTreeVectorOfVectorsAdaptor<State, ldouble> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,ldouble> Temp_Tree;

#endif /* VAR_H */