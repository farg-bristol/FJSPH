/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VAR_H
#define VAR_H

#include <vector>
#include "../Eigen/Core"
#include "../Eigen/StdVector"
#include "../Eigen/LU"
#include "../Eigen/Geometry"
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
constexpr int simDim = 3;
typedef double ldouble;
typedef unsigned int uint;
constexpr uint nthreads = 12;

/****** Eigen vector definitions ************/
typedef Eigen::Matrix<ldouble,simDim,1> StateVecD;
typedef Eigen::Matrix<int,simDim,1> StateVecI;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<ldouble, simDim+1,1> DensVecD;
typedef Eigen::Matrix<ldouble, simDim+1, simDim+1> DensMatD;


/*A structure for the coordinates to define a panel.*/
typedef struct Panel
{ /*A and B are 1/4 chord bounds for the vortex.
	C is the control Point location*/
	Panel(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C,
			Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d p4)
	: A(A), B(B), C(C), p1(p1), p2(p2), p3(p3), p4(p4) {}

	Eigen::Vector3d A, B, C;
	Eigen::Vector3d p1, p2, p3, p4;
}Panel;

typedef class VLM
{
	public:
		VLM()
		{
			/*Define x and y end coordinates of the wing*/
			coords[0] = 5;
			coords[1] = 2;

			/*Split it up into this many panels (Will be doubled on the other side)*/
			panels(0) = 10;
			panels(1) = 5;
			npanels = 2*panels[0]*panels[1];
			// nverts = 2*(panels[0]+1)*(panels[1]+1);

			/*Initialise matrices*/
			aInf = Eigen::MatrixXd(npanels,npanels);
			gamma = Eigen::VectorXd(npanels);
			RHS = Eigen::VectorXd(npanels);

			panelData.reserve(npanels);
			// panelxyz.reserve(nverts);

			/*Define Angle of attack*/
			AoA = 12 * M_PI/180;

			/*Want freestream to be aligned with jet axis*/
			Freestream[0]= 0;
			Freestream[1]= 1;
			Freestream[2]= 0;

			Freestream = Freestream.normalized();

			MakeMatrix();

		}

		void GetGamma(Eigen::Vector3d inf)
		{
			Freestream = inf;
			/*Find the inflence matrix, and invert it to find gamma*/

			/*Find the influence matrix aInf*/
 			for(int i=0; i < npanels; ++i)
 			{
 				for(int j=0; j < npanels; ++j)
 				{
 					/*aInf[i,j] = the influence of vortex j on control point i*/
 					Eigen::Vector3d A, B, C;
 					A = panelData[i].A;
 					B = panelData[i].B;
 					C = panelData[j].C;

					aInf(j,i) = FindInfluence(A,B,C).dot(norm);
					
				}

			RHS(i) = -Freestream.dot(norm);
			
			}
			/*Find Gamma */
 			gamma = aInf.fullPivLu().solve(RHS);

 			// for (int i= 0; i < gamma.rows(); ++i)
 			// {
 			// 	std::cout << gamma[i] << std::endl;
 			// }
 		}

		Eigen::Vector3d getVelocity(Eigen::Vector3d pos)
		{	/*Find velocity for a particle at its position*/
			Eigen::Vector3d vel = Eigen::Vector3d::Zero();

			for(int i = 0; i<npanels; ++i)
				vel += gamma(i)*FindInfluence(panelData[i].A,panelData[i].B,pos);
			
			return vel+Freestream;
		}

		std::vector<Panel> panelData;

	protected:
		void MakeMatrix(void)
		{
			/*Define the steps for each dimension*/
			double dz = coords[0]/double(panels[0]);
			double dy = cos(AoA)*coords[1]/double(panels[1]);
			double dx = sin(AoA)*coords[1]/double(panels[1]);

			double y0 = cos(AoA)*coords[1];
			double x0 = sin(AoA)*coords[1];

			/*Find 1/4 panel points (assuming symmetry) */
			/*Find 3/4 Control Points*/
			for(int i = 0; i < 2*panels[0]; ++i)
			{
				for (int j=0; j < panels[1]; ++j)
				{
					/*Panel Verticies (for visual)*/
					Eigen::Vector3d p1(x0-dx*(double(j)),dy*(double(j))-y0,double(i)*dz);
					Eigen::Vector3d p2(x0-dx*(double(j+1)),dy*(double(j+1))-y0,double(i)*dz);
					Eigen::Vector3d p3(x0-dx*(double(j+1)),dy*(double(j+1))-y0,double(i+1)*dz);
					Eigen::Vector3d p4(x0-dx*(double(j)),dy*(double(j))-y0,double(i+1)*dz);


					/*1/4 chord points*/
					Eigen::Vector3d A(x0-dx*(0.25+double(j)),dy*(0.25+double(j))-y0,double(i)*dz);
					Eigen::Vector3d B(x0-dx*(0.25+double(j)),dy*(0.25+double(j))-y0,double(i+1)*dz);
					/*3/4 control point*/
					Eigen::Vector3d C(x0-dx*(0.75+double(j)),dy*(0.75+double(j))-y0,dz*(0.5+double(i)));
					panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
				}
			}
	
			
			/*Find normal*/
			Eigen::Vector3d a = panelData[0].C-panelData[0].A;
			Eigen::Vector3d b = panelData[0].C-panelData[0].B;

			norm = a.cross(b);
			norm = norm.normalized();

			/*End initialisation*/
		}

		Eigen::Vector3d FindInfluence(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
 		{	
			Eigen::Vector3d r0, r1, r2, inf;

			/*Bounded Vortex*/
			Eigen::Vector3d coefAB;

			r0 = B-A;
			r1 = C-A;
			r2 = C-B;

			coefAB = (1/(4*M_PI))*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(r0.dot(r1.normalized())-r0.dot(r2.normalized()));


			/*Horseshoe vortex from point A*/
			Eigen::Vector3d coefA;
		
			inf = A + Freestream;
			r2 = C - A;
			r1 = C - inf;
			r0 = A - inf;

			coefA = (1/(4*M_PI))*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(1-r0.dot(r2)/(r0.norm()*r2.norm()));

			/*Horseshoe vortex from point B*/
			Eigen::Vector3d coefB;

			/*Vector B to infinity*/
			inf = B + Freestream;
			r1 = C - B;
			r2 = C - inf;
			r0 = inf - B;

			coefB = (1/(4*M_PI))*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(r0.dot(r1)/(r1.norm()*r0.norm())+1.0);

			return (coefAB+coefB+coefA);			
 		}

		int npanels, nverts;
		double AoA;
		
		Eigen::MatrixXd aInf;

		Eigen::VectorXd gamma;
		Eigen::VectorXd RHS;

		Eigen::Vector2d coords;
		Eigen::Vector2i panels;

		Eigen::Vector3d norm;
		Eigen::Vector3d Freestream;

}VLM;

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

		ldouble tmax;
		ldouble ycoef;
		ldouble woccl;
		ldouble L;
		ldouble Cf, Ck, Cd, Cb, Cdef;
		ldouble td;
		ldouble omega;
		ldouble nfull;
}AERO;

/*Simulation parameters*/
typedef struct SIM {
	StateVecI xyPART; 				/*Starting sim particles in x and y box*/
	uint simPts,bndPts,totPts;	    /*Simulation particles, Boundary particles, total particles*/
	uint nrefresh; 					/*Crossflow refresh particle number*/
	uint nmax;                      /*Max add-particle calls*/
	uint outframe;	                /*Terminal output frame interval*/
	uint addcount;					/*Current Number of add-particle calls*/
	double Pstep,Bstep;				/*Initial spacings for particles and boundary*/
	StateVecD Box;					/*Box dimensions*/
	StateVecD Start; 				/*Sim box bottom left coordinate*/
	uint subits;                    /*Max number of sub-iterations*/
	uint Nframe; 			        /*Max number of frames to output*/
	double dt, t, framet;			/*Timestep, Simulation + frame times*/
	double beta,gamma;				/*Newmark-Beta Parameters*/
	double maxmu;                   /*Maximum viscosity component (CFL)*/
	int Bcase, Bclosed;				/*What boundary shape to take*/
	int outform;                    /*Output type. Fluid properties or Research.*/
	int frameout;                   /**/
	VLM vortex;
} SIM;

/*Fluid and smoothing parameters*/
typedef struct FLUID {
	ldouble H, HSQ, sr; 			/*Support Radius, SR squared, Search radius*/
	ldouble rho0; 					/*Resting Fluid density*/
	ldouble Simmass, Boundmass;		/*Particle and boundary masses*/
	ldouble correc;					/*Smoothing Kernel Correction*/
	ldouble alpha,Cs,mu;		/*}*/
	ldouble sig;					/* Fluid properties*/
	ldouble gam, B; 				/*}*/
	ldouble mug;					/* Gas Properties*/
	ldouble rhog;
	ldouble contangb;				/*Boundary contact angle*/
	ldouble volume;					/*Particle volume*/
	AERO avar;
	//double front, height, height0;		/*Dam Break validation parameters*/

}FLUID;

/*Crossflow parameters*/
typedef struct CROSS
{
	int acase;	                        /*Aerodynamic force case*/
	StateVecD vJet, vInf;               /*Crossflow Parameters: Jet + Freestream velocity*/
	ldouble Acorrect;					/*Correction factor for aero force*/
	ldouble a;                          /*Tuning parameters*/
	ldouble b;                          /*Tuning parameters*/
	ldouble h1;                         /*Tuning parameters*/
	ldouble h2;                         /*Tuning parameters*/
}CROSS;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(StateVecD X, StateVecD Vi, StateVecD Fi,
			ldouble Rhoi, ldouble Mi, int bound)
		{
			xi = X;	v = Vi; f = Fi;
			rho = Rhoi; m = Mi; b = bound;
			Sf = StateVecD::Zero(); 
			Af = StateVecD::Zero();
			normal = StateVecD::Zero();
			Rrho = 0.0;
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

		// void erase_list(void)
		// {
		// 	list.erase(list.begin(),list.end());
		// }
		// std::vector<uint> list; /*Neighbour list*/

		StateVecD xi, v, f, Sf, Af, normal;
		ldouble rho, p, Rrho, m, theta;
		uint b; //What state is a particle. Boundary, forced particle or unforced
}Particle;

typedef class Part {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Part(Particle pi)
		{
			xi = pi.xi; v = pi.v; 
			rho = pi.rho; 
			p = pi.p;
			m = pi.m;
			b = pi.b;
			Sf = pi.Sf;
			normal = pi.normal;
		}

		int size() const
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		double operator[](int a) const
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}

		void operator=(Particle pi)
		{
			xi = pi.xi; v = pi.v; 
			rho = pi.rho; 
			p = pi.p;
			m = pi.m;
			b = pi.b;
			Sf = pi.Sf;
		}

		// void erase_list(void)
		// {
		// 	list.erase(list.begin(),list.end());
		// }
		// std::vector<uint> list; /*Neighbour list*/

		StateVecD xi, v, Sf, normal;
		ldouble rho, p, m;
		uint b; //What state is a particle. Boundary, forced particle or unforced
}Part;

typedef std::vector<Particle> State;

/* Neighbour search tree containers */
typedef std::vector<std::vector<uint>> outl;
typedef KDTreeVectorOfVectorsAdaptor<State, ldouble> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<StateVecD>,ldouble> Temp_Tree;

#endif /* VAR_H */
