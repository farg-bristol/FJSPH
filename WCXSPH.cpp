/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/***			        + XSPH Correction (Also described in Monaghan)                    ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)                                  ***/
/*** Density Reinitialisation: Colagrossi, A. and Landrini, M. (2003): MLS                ***/
/*** Smoothing Kernel: Wendland's C2 ***/
/*** Integrator: Newmark-Beta ****/
/*** Variable Timestep Criteria: CFL + Monaghan, J.J. (1989) conditions ***/

#include <omp.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <string.h>
#include <sstream>
#include <chrono>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include "Eigen/LU"
#include "NanoFLANN/nanoflann.hpp"
#include "NanoFLANN/utils.h"
#include "NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"
#include "PDS/poisson_disk_sampling.h"


#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace nanoflann;

double maxmu = 0.0;			/*CFL Parameter*/

typedef struct SIM {
	Vector2i xyPART; 						/*Starting sim particles in x and y box*/
	unsigned int SimPts,bound_parts,npts;	/*Number of particles*/
	unsigned int nrefresh, nmax, outframe;	/*Crossflow Particles*/
	unsigned int addcount, aircount;		/*Number of add particle calls*/
	double Pstep,Bstep;						/*Initial spacings*/
	Vector2d Box;							/*Box dimensions*/
	Vector2d Start; 				/*Starting sim box dimensions*/
	unsigned int subits,Nframe; 			/*Timestep values*/
	double dt, t, framet;					/*Simulation + frame times*/
	double tNorm;
	double beta,gamma;						/*Newmark-Beta Parameters*/
	int Bcase, Bclosed;						/*What boundary shape to take*/
	double H,HSQ, sr; 						/*Support Radius + Search radius*/
	int acase;								/*Aerodynamic force case*/
} SIM;

typedef struct FLUID {
	double rho0; 						/*Resting density*/
	double Simmass, Boundmass;			/*Particle and boundary masses*/
	double alpha,eps,Cs, mu;			/*}*/
	double sig;							/*} Fluid properties*/
	double gam, B; 						/*}*/
	double contangb;					/*Boundary contact angle*/
	double front, height, height0;		/*Dam Break validation parameters*/
	Vector2d vJet, vInf;
}FLUID;

SIM svar;
FLUID fvar;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(Vector2d X, Vector2d Vi, Vector2d Fi, 
			double Rhoi, double Mi, bool bound, bool l)	
		{
			xi = X;	v = Vi; f = Fi;
			V(0.0,0.0);	Sf(0.0,0.0);
			rho = Rhoi;
			Rrho = 0;
			m = Mi;
			b = bound;
			left = l;
		}

		int size() const 
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		double operator[](int a) const 
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}
	
		Vector2d xi, v, V, f, Sf;
		double rho, p, Rrho, m;
		bool b, left;
}Particle;
typedef std::vector<Particle> State;
///****** Initialise the particles memory *********/
State pn;	/*Particles at n*/
State pnp1; 	/*Particles at n+1*/

/*Define Neighbour search types*/
typedef std::vector<std::vector<size_t>> outl;
outl outlist;
typedef KDTreeVectorOfVectorsAdaptor<State, double> KD_Tree;
typedef KDTreeVectorOfVectorsAdaptor<vector<Vector2d>, double> Temp_Tree;
nanoflann::SearchParams params;

struct EVecTraits
{
	typedef double ValueType;

	static constexpr auto kSize = 2;

	static ValueType Get(const Eigen::Vector2d& v, const std::size_t i)
	{
		return *(&v[0] + i);
	}

	static void Set(Eigen::Vector2d* const v, const std::size_t i, const ValueType val)
	{
		*(&v[0][0] + i) = val;
	} 
};

namespace thinks {
	namespace poisson_disk_sampling {
		template<>
		struct VecTraits<Eigen::Vector2d>
		{
			typedef double ValueType;

			static constexpr auto kSize = 2;

			static ValueType Get(const Eigen::Vector2d& v, const std::size_t i)
			{
				return *(&v[0] + i);
			}

			static void Set(Eigen::Vector2d* const v, const std::size_t i, const ValueType val)
			{
				*(&v[0][0] + i) = val;
			} 
		};
	}
}

void write_header() 
{
	cout << "******************************************************************" << endl << endl;
	cout << "                              WCXSPH                              " << endl << endl;
	cout << "        Weakly Compressible Smoothed Particle Hydrodynamics       " << endl;
	cout << "                       with XSPH correction                       " << endl;
	cout << "                      for tide-breaking case                      " << endl << endl;
	cout << "                         James O. MacLeod                         " << endl;
	cout << "                    University of Bristol, U.K.                   " << endl << endl;
	cout << "******************************************************************" << endl << endl;
}

int getInt(ifstream& In)
{
	string line;
	getline(In,line);
	int i = stoi(line);
	return i;
}

double getDouble(ifstream& In)
{
	string line;
	getline(In,line);
	double d = stod(line);
	return d; 
}

std::string getString(ifstream& In)
{
	string line;
	getline(In,line);
	return line; 
}

void DefaultInput(void) 
{
	//Timestep Parameters
		svar.framet = 0.1;		/*Frame timestep*/
		svar.Nframe = 2500;		/*Number of output frames*/
		svar.outframe = 50;		/*Terminal output frame interval*/
		svar.subits = 10;			/*Newmark-Beta iterations*/
		svar.nmax = 2000;			/*Max No particles (if dynamically allocating)*/	
		svar.beta = 0.25;			/*Newmark-Beta parameters*/
		svar.gamma = 0.5;	
		
		//Simulation window parameters
		svar.xyPART(40,40); /*Number of particles in (x,y) directions*/
		svar.Start(0.2,0.2); /*Simulation particles start + end coords*/
		svar.Bcase = 1; 		/*Boundary case*/
		svar.Box(3,2); 		/*Boundary dimensions*/
		svar.Pstep = 0.01;	/*Initial particle spacing*/
		svar.Bstep = 0.6; 	/*Boundary factor of particle spacing (dx = Pstep*Bstep)*/
		
		// SPH Parameters
		svar.H =  3.0*svar.Pstep; /*Support Radius*/
		
		//Fluid Properties
		fvar.alpha = 0.025;			/*Artificial viscosity factor*/
  		fvar.eps = 0.25;			/*XSPH Correction factor*/
		fvar.contangb = 150.0;		/*Boundary contact angle*/
		fvar.rho0 = 1000.0;  		/*Rest density*/
		fvar.Cs = 50; 	 			/*Speed of sound*/
		fvar.mu = 0.0001002;		/*Viscosity*/
		fvar.sig = 0.0728;			/*Surface tension*/
		fvar.vJet(1) = 10.0;

		/*Universal parameters based on input values*/
	  	svar.dt = 0.00002; 		/*Initial timestep*/
	  	svar.t = 0.0;				/*Total simulation time*/
	  	svar.HSQ = svar.H*svar.H; 
		svar.sr = 4*svar.HSQ; 	/*KDtree search radius*/
		svar.Bclosed = 0; 		/*Boundary begins open*/
	  	svar.SimPts = svar.xyPART(0)*svar.xyPART(1); /*total sim particles*/
	  	svar.aircount = 0;

	  	Vector2d Finish(svar.Start(0)+1.0*svar.Pstep*svar.xyPART(0),
						svar.Start(1)+1.0*svar.Pstep*svar.xyPART(1));
		
		fvar.Simmass = fvar.rho0* /*Mass from spacing and density*/
	  		((Finish(0)-svar.Start(0))*(Finish(1)-svar.Start(1)))/(1.0*svar.SimPts);
		fvar.Boundmass = fvar.Simmass*svar.Bcase;
		
		fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
		fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
		fvar.height = Finish(1);

}

void GetInput(char* InFile)
{

	std::ifstream in(InFile);
  	if(in.is_open()) 
  	{	/*Simulation parameters*/
  		cout << "Input file opened. Reading settings..." << endl;
  		svar.framet = getDouble(in);
  		svar.Nframe = getInt(in);
  		svar.outframe = getInt(in);
  		svar.subits = getInt(in);
  		svar.nmax = getInt(in);	
  		svar.beta = getDouble(in);
  		svar.gamma = getDouble(in);
  		svar.xyPART(0) = getInt(in); 
  		svar.xyPART(1) = getInt(in);
  		svar.Start(0) = getDouble(in);
  		svar.Start(1) = getDouble(in);
  		svar.Bcase = getInt(in);
  		svar.Box(0) = getDouble(in);
  		svar.Box(1) = getDouble(in);
  		svar.Pstep = getDouble(in);
  		svar.Bstep = getDouble(in);
  		double Hfac = getDouble(in); /*End of state read*/
  		svar.H= Hfac*svar.Pstep;

		/*Fluid parameters read*/
  		fvar.alpha = getDouble(in);
  		fvar.eps = getDouble(in);
  		fvar.contangb = getDouble(in);
  		fvar.rho0 = getDouble(in);
  		fvar.Cs = getDouble(in);
  		fvar.mu = getDouble(in);
  		fvar.sig = getDouble(in);
  		fvar.vJet(0) = 0.0; fvar.vJet(1) = getDouble(in); 
  		fvar.vInf(0) = getDouble(in); fvar.vInf(1) = 0.0;
  		svar.acase = getInt(in);

		in.close();

  	}
  	else {
	    cerr << "Error opening the input file." << endl;
	    exit(-1);
  	}
  	
  	/*Universal parameters based on input values*/
  	svar.dt = 0.00002; 		/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	svar.HSQ = svar.H*svar.H; 
	svar.sr = 4*svar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
  	svar.SimPts = svar.xyPART(0)*svar.xyPART(1); /*total sim particles*/
  	svar.aircount = 0;

  	Vector2d Finish(svar.Start(0)+1.0*svar.Pstep*svar.xyPART(0),
					svar.Start(1)+1.0*svar.Pstep*svar.xyPART(1));
	
	fvar.Simmass = fvar.rho0* /*Mass from spacing and density*/
  		((Finish(0)-svar.Start(0))*(Finish(1)-svar.Start(1)))/(1.0*svar.SimPts);
	fvar.Boundmass = fvar.Simmass*svar.Bcase;
	
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
	fvar.height = Finish(1);
}

///******Wendland's C2 Quintic Kernel*******
double W2Kernel(double dist,double H) 
{
	double q = dist/H;
	double correc = (7/(4*M_PI*H*H));
	return (pow(1-0.5*q,4))*(2*q+1)*correc;
}

/*Gradient*/
Vector2d W2GradK(Vector2d Rij, double dist,double H)
{
	double q = dist/H;
	double correc = (7/(4*M_PI*H*H));
	return 5.0*(Rij/(H*H))*pow(1-0.5*q,3)*correc;
}

/*2nd Gradient*/
double W2Grad2(Vector2d Rij, double dist,double H) 
{
	double q = dist/H;
	double correc = (7/(4*M_PI*H*H));
	return Rij.dot(Rij)*(5.0*correc/(H*H))*(2*q-1)*pow(1-0.5*q,2);
}

///**************** Update neighbour list **************
void FindNeighbours(KD_Tree &mat_index)
{
	outlist.erase(outlist.begin(),outlist.end());
	double search_radius = svar.sr;
	/*Find neighbour list*/
	for(size_t i=0; i<pnp1.size(); ++i)
	{		
		std::vector<std::pair<size_t,double>> matches; /* Nearest Neighbour Search*/
		
		mat_index.index->radiusSearch(&pnp1[i].xi[0], search_radius, matches, params);

		std::vector<size_t> temp;
		for (auto &j:matches) 
		{
			temp.emplace_back(j.first); 	
		}
		outlist.emplace_back(temp);		
	}		
}

void AddPoints(void)
{	
	// cout << "Adding points..." << endl;
	Vector2d v = fvar.vJet;  /*Jet velocity*/
	Vector2d f = Vector2d::Zero(); 
	double rho=fvar.rho0; 
	double jetS = svar.Start(0)+2*svar.H;
	double jetE = svar.Start(0)+svar.Start(1) - 2*svar.H;
	svar.nrefresh = 0;
	

	if (svar.acase == 5)
	{
		/*Create the simulation particles*/
		for( double x = jetS; x<jetE - (jetE-jetS)/2; x+=svar.Pstep) 
		{ /*Do the left set of points*/
			Vector2d xi(x,-svar.Box[1]);
			int pos = svar.npts+svar.nrefresh;		
			pn.insert(pn.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,false,true));
			pnp1.insert(pnp1.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,false,true));
			++svar.SimPts;
			++svar.nrefresh;
		}

		double start2 = pn.back().xi[0]+svar.Pstep;
		for( double x = start2; x<=jetE; x+=svar.Pstep) 
		{	/*Do the right set of points*/
			Vector2d xi(x,-svar.Box[1]);
			int pos = svar.npts+svar.nrefresh;		
			pn.insert(pn.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,false,false));
			pnp1.insert(pnp1.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,false,false));
			++svar.SimPts;
			++svar.nrefresh;
		}
	}
	else
	{
		/*Create the simulation particles*/
		for( double x = jetS; x<jetE - (jetE-jetS)/2; x+=svar.Pstep) 
		{ /*Do the left set of points*/
			Vector2d xi(x,-svar.Box[1]);
			pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,true));
			pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,true));
			++svar.SimPts;
			++svar.nrefresh;
		}

		double start2 = pn.back().xi[0]+svar.Pstep;
		for( double x = start2; x<=jetE; x+=svar.Pstep) 
		{	/*Do the right set of points*/
			Vector2d xi(x,-svar.Box[1]);		
			pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,false));
			pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,false));
			++svar.SimPts;
			++svar.nrefresh;
		}
	}

	svar.npts += svar.nrefresh;
	++svar.addcount;
	// cout << "New points: " << svar.nrefresh << "  npts: " << 
		// svar.npts << " SimPts: "<< svar.SimPts <<  endl;
}

void CloseBoundary(void) 
{
	cout << "Closing boundary..." << endl;
	Vector2d v = Vector2d::Zero();
	Vector2d f = Vector2d::Zero();
	double rho=fvar.rho0;

	double holeS = svar.Start(0); /*Distance from Box start to hole*/
	double holeD = svar.Start(1); /*Diameter of hole (or width)*/
	double stepb = (svar.Pstep*svar.Bstep);
	int Nb = ceil((holeD)/stepb);
	stepb = holeD/(1.0*Nb);
	State temp;

	for(double x = holeS; x<holeS+holeD; x+=stepb)
	{
		Vector2d xi(x,0.0);
		temp.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
	}
	/*Insert particles at the end of boundary particles*/
	pn.insert(pn.begin()+svar.bound_parts,temp.begin(),temp.end());
	pnp1.insert(pnp1.begin()+svar.bound_parts,temp.begin(),temp.end());
	svar.bound_parts+=temp.size(); /*Adjust counter values*/
	svar.npts +=temp.size();
}

///**************** RESID calculation **************
void Forces(KD_Tree &mat_index)
{
	maxmu=0; 					/* CFL Parameter */
	double alpha = fvar.alpha; 	/* Artificial Viscosity Parameter*/
	double eps = fvar.eps; 		/* XSPH Influence Parameter*/
	double numpartdens = 0.0;
	/*Surface tension factor*/
	const static double lam = (6.0/81.0*pow((2.0*svar.H),4.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0)); 

	/********* LOOP 1 - all points: Calculate numpartdens ************/
	for (size_t i=0; i< svar.npts; ++i) 
	{	
		Vector2d pi = pnp1[i].xi;
		for (size_t j=0; j<outlist[i].size(); ++j) 
		{ /* Surface Tension calcs */
			Vector2d pj = pnp1[outlist[i][j]].xi;
			double r = (pj-pi).norm();
			numpartdens += W2Kernel(r,svar.H);
		}	
	}
	numpartdens=numpartdens/(1.0*svar.npts);

	// #pragma omp parallel
	{
		// #pragma omp for
		
	/******** LOOP 2 - Boundary points: Calculate density and pressure. **********/
		for (size_t i=0; i< svar.bound_parts; ++i) 
		{	/*Find variation in density for the boundary (but don't bother with force)*/
			double Rrhocontr = 0.0;
			Particle pi = pnp1[i];

			for(size_t j=0; j<outlist[i].size(); ++j)
			{
				Particle pj = pnp1[outlist[i][j]];
				Vector2d Rij = pj.xi-pi.xi;
				Vector2d Vij = pj.v-pi.v;
				double r = Rij.norm();
				Vector2d Grad = W2GradK(Rij, r,svar.H);
				Rrhocontr -= pj.m*(Vij.dot(Grad));  
			}
			pnp1[i].Rrho = Rrhocontr; /*drho/dt*/
		}

		
	
	/******* LOOP 3 - only for ghost particle case: Find surface points. *********/
		if(svar.Bcase == 3 && svar.acase == 5)
		{
			/*Delete previous air particles*/
			//for (index p=std::next(pnp1.begin(),svar.bound_parts+svar.SimPts); p!=pnp1.end(); ++p)
			while (pnp1.size()!=svar.npts)
			{
				pnp1.pop_back();
				pn.pop_back();
			}

			std::vector<Vector2d> temp; /*Temporary storage for air particles*/
			for (size_t i=svar.bound_parts; i< svar.npts; ++i) 
			{	/*Find the surface of fluid particles.*/
				Particle pi = pnp1[i];
				Vector2d SurfC= Vector2d::Zero(); 
				pi.left = false;

				for (size_t j=0; j < outlist[i].size(); ++j) 
				{
					Particle pj = pnp1[outlist[i][j]];
					/*Check if the position is the same, and skip the particle if yes*/
					if(pi.xi == pj.xi)
						continue;

					Vector2d Rij = pj.xi-pi.xi;
					double r = Rij.norm();

					/*Surface Tension as described by Nair & Poeschel (2017)*/
					double fac=1.0;
					if(pj.b==true) 
			            fac=(1+0.5*cos(M_PI*(fvar.contangb/180)));
			        //cout << lam(svar.H) << endl;
					double sij = 0.5*pow(numpartdens,-2.0)*(fvar.sig/lam)*fac;
					SurfC -= (Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*svar.H))/pj.m;
				}

				if (SurfC.norm()>0.05)
				{
					/*Create particles... */
					namespace pds = thinks::poisson_disk_sampling;
					double radius = svar.Pstep;
					std::array<double,2> xmin = {pn[i].xi[0]-2.0*svar.H, pn[i].xi[1]-2.0*svar.H};
					std::array<double,2> xmax = {pn[i].xi[0]+2.0*svar.H,pn[i].xi[1]+2.0*svar.H};
					// cout << "Centre Coord: " << pn[i].xi[0] << " " << pn[i].xi[1] << endl;
					// cout << "Min Coords: " << xmin[0] << "  " << xmin[1] << endl;
					// cout << "Max Coords: " << xmax[0] << "  " << xmax[1] << endl;

					std::vector<Vector2d> samples = 
					pds::PoissonDiskSampling<double,2,Vector2d,EVecTraits>(radius,xmin,xmax);
					
					// cout << samples.size() << endl;
					for (auto j:samples)
					{
						Vector2d xi(j[0],j[1]); 
						temp.emplace_back(xi);
					}
				}

			}

			if(temp.size()!=0)
			{
				Temp_Tree temp_index(2,temp,10);
				temp_index.index->buildIndex();
				double search_radius = svar.Pstep*svar.Pstep;
				//std::vector<size_t> delete_list;
				//cout << temp.size() << endl;

				for (auto i=temp.begin(); i!=temp.end(); )
				{	/*Check for duplicate particles and delete them when too close.*/
					Vector2d xi = *i;
					std::vector<std::pair<size_t,double>> matches;
					temp_index.index->radiusSearch(&xi[0],search_radius,matches,params);
					//cout << matches.size() << endl;
					if (matches.size()!=1)
					{

						for (auto j:matches)
						{
							if (j.second == 0.0)
								continue;

							temp[j.first] = temp.back();
							temp.pop_back();
						}
						temp_index.index->buildIndex();
					}
					else
						++i;
				}

				for (size_t i=svar.bound_parts; i < svar.npts; ++i)
				{	/*Check for particles inside the fluid*/
					std::vector<std::pair<size_t,double>> matches;
					temp_index.index->radiusSearch(&pn[i].xi[0],search_radius,matches,params);
					//cout << matches.size()<< endl;
					if (matches.size()!=0)
					{
						for (auto j:matches)
						{
							temp[j.first] = temp.back();
							temp.pop_back();
						}
						temp_index.index->buildIndex();
					}
				}

				// search_radius = 2*svar.Pstep*svar.Pstep;
				// for (size_t i=0; i< svar.bound_parts; ++i)
				// {	/*Check for particles next to the boudnary*/
				// 	std::vector<std::pair<size_t,double>> matches;
				// 	temp_index.index->radiusSearch(&pn[i].xi[0],search_radius,matches,params);
				// 	//cout << matches.size()<< endl;
				// 	if (matches.size()!=0)
				// 	{
				// 		for (auto j:matches)
				// 		{
				// 			temp[j.first] = temp.back();
				// 			temp.pop_back();
				// 		}
				// 		temp_index.index->buildIndex();
				// 	}
				// }

				/*Place particles in the simulation vector*/
				svar.aircount = temp.size();

				double rho = 1.225;
				Vector2d f = Vector2d::Zero();
				double airmass = fvar.Simmass*rho/fvar.rho0;
				for (auto i:temp)
				{
					pnp1.emplace_back(Particle(i,fvar.vInf,f,rho,airmass,false,false));
					pn.emplace_back(Particle(i,fvar.vInf,f,rho,airmass,false,false));
				}

				mat_index.index->buildIndex();
				FindNeighbours(mat_index);
			}
		}
		
		// cout << "Npts: " << svar.npts << " Air Count: " << svar.aircount << endl;
		// cout << "Bound Parts: " << svar.bound_parts << " Sim Points: " << svar.SimPts << endl; 
		// cout << "Pn size: " << pn.size() << " PnP1 Size: " << pnp1.size() << endl;
		// cout << "Outlist Size: " << outlist.size() << endl;

		for (size_t i=svar.bound_parts; i< svar.npts + svar.aircount; ++i) 
		{	/*Do force calculation for fluid particles.*/
			Particle pi = pnp1[i];
			pi.V = pi.v;
			if (svar.Bcase == 3 && svar.acase == 3)
				pi.left = true;

			double Rrhocontr = 0.0;
			Vector2d contrib= Vector2d::Zero();
			Vector2d visc = Vector2d::Zero();
			Vector2d SurfC= Vector2d::Zero(); 

			vector<double> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/

			for (size_t j=0; j < outlist[i].size(); ++j) 
			{	/*Find force and density variation for particles*/
				Particle pj = pnp1[outlist[i][j]];

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.xi == pj.xi)
					continue;

				Vector2d Rij = pj.xi-pi.xi;
				Vector2d Vij = pj.v-pi.v;
				double r = Rij.norm();
				double Kern = W2Kernel(r,svar.H);
				Vector2d Grad = W2GradK(Rij, r,svar.H);
				
				/*Pressure and artificial viscosity calc - Monaghan 1994 p.400*/
				double rhoij = 0.5*(pi.rho+pj.rho);
				double cbar= 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));
				double vdotr = Vij.dot(Rij);
				double muij= svar.H*vdotr/(r*r+0.01*svar.HSQ);
				mu.emplace_back(muij);
				double pifac = alpha*cbar*muij/rhoij;
		
				if (vdotr > 0) pifac = 0;
				contrib += pj.m*Grad*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));

				/*Laminar Viscosity (Morris)*/
				visc -= Vij*(pj.m*fvar.mu)/(pi.rho*pj.rho)
					*(1.0/(r*r+0.01*svar.HSQ))*Rij.dot(Grad);

				/*Surface Tension as described by Nair & Poeschel (2017)*/
				double fac=1.0;
				if(pj.b==true) 
		            fac=(1+0.5*cos(M_PI*(fvar.contangb/180)));
		        //cout << lam(svar.H) << endl;
				double sij = 0.5*pow(numpartdens,-2.0)*(fvar.sig/lam)*fac;
				SurfC -= (Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*svar.H))/pj.m;
				
				//SurfC += (pj.m/pj.rho)*Rij*Kern; 

				/* XSPH Influence*/
				pi.V+=eps*(pj.m/rhoij)*Kern*Vij; 

				/*drho/dt*/
				Rrhocontr -= pj.m*(Vij.dot(Grad));	

				if (svar.Bcase == 3 && svar.acase == 3)
				{
					double num = -Rij.dot(fvar.vInf);
					double denom = Rij.norm()*fvar.vInf.norm();
					if (num/denom > 0.9)
						pi.left = false;
				}

			}
			
			/*Crossflow force*/
			Vector2d Fd= Vector2d::Zero();
			if (svar.Bcase == 3 && pi.xi[1] > svar.Pstep)
			{
				switch(svar.acase)
				{
					default:
						break;
					
					case 0: /*No aero force*/
						break;

					case 1:
					{  /*All upstream particles at start*/
						if( pi.left == true)
						{
							Vector2d Vdiff = fvar.vInf - pi.V;
							double Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
							double Cd;
							
							if (Re < 3500)
							 	Cd = 0.01*(1+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/Re);
							else 
								Cd = 0.01*(1+0.197*pow(Re,0.63)+2.6e-4*pow(Re,1.38))*(24.0/Re);

							// cout << "Reynolds: " << Re << " Cd: " << Cd << endl;

							Fd = Vdiff.normalized()*Cd*(2*svar.Pstep)*1.225*Vdiff.squaredNorm();
							Fd[1] = 0.0;
						}
						
						break;
					}

					case 2:
					{	/* Surface particles*/
						if (SurfC.norm() > 0.05)
						{
							Vector2d Vdiff = fvar.vInf - pi.V;
							double Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
							double Cd;
							
							if (Re < 3500)
							 	Cd = 0.01*(1+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/Re);
							else 
								Cd = 0.01*(1+0.197*pow(Re,0.63)+2.6e-4*pow(Re,1.38))*(24.0/Re);

							// cout << "Reynolds: " << Re << " Cd: " << Cd << endl;

							Fd = Vdiff.normalized()*Cd*(2*svar.Pstep)*1.225*Vdiff.squaredNorm();
							Fd[1] = 0.0;
						}
						break;
					}

					case 3:
					{	/* Left Surface particles*/
						if(pi.left == true && SurfC.norm() > 0.05)
						{
							
							Vector2d Vdiff = fvar.vInf - pi.V;
							double Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
							double Cd;
							
							if (Re < 3500)
							 	Cd = 0.1*(1+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/Re);
							else 
								Cd = 0.1*(1+0.197*pow(Re,0.63)+2.6e-4*pow(Re,1.38))*(24.0/Re);

							// cout << "Reynolds: " << Re << " Cd: " << Cd << endl;

							Fd = Vdiff.normalized()*Cd*(2*svar.Pstep)*1.225*Vdiff.squaredNorm();
							Fd[1] = 0.0;
						
						}
						break;
					}

					case 4:
					{	/* Surface particles proportional to ST*/
						/*Work in progress...*/
						break;
					}

				}
				
			}
			

			pi.Rrho = Rrhocontr; /*drho/dt*/
			pi.f= contrib - SurfC*fvar.sig/pi.m + Fd/pi.m;

			pi.Sf = Fd/pi.m;
			pi.f(1) -= 9.81; /*Add gravity*/

			pnp1[i]=pi; //Update the actual structure

			//CFL f_cv Calc
			double it = *max_element(mu.begin(),mu.end());
			if (it > maxmu)
				maxmu=it;
		} /*End of sim parts*/
	}
}

///*Density Reinitialisation using Least Moving Squares as in A. Colagrossi (2003)*
void DensityReinit() 
{
	Vector3d one(1.0,0.0,0.0);
	for(size_t i=0; i< svar.npts; ++i)
	{
		Matrix3d A= Matrix3d::Zero();
		//Find matrix A.
		Particle pi = pnp1[i];
		for (size_t j=0; j< outlist[i].size(); ++j) 
		{
			Particle pj = pnp1[outlist[i][j]];
			Vector2d Rij = pi.xi-pj.xi;
			Matrix3d Abar;	
			Abar << 1      , Rij(0)        , Rij(1)        ,
				    Rij(0) , Rij(0)*Rij(0) , Rij(1)*Rij(0) ,
				    Rij(1) , Rij(1)*Rij(0) , Rij(1)*Rij(1) ;

			A+= W2Kernel(Rij.norm(),svar.H)*Abar*pj.m/pj.rho;
		}
				
		Vector3d Beta;
		//Check if A is invertible
		FullPivLU<Matrix3d> lu(A);
		if (lu.isInvertible())
			Beta = lu.inverse()*one;
		else
			Beta = (1/A(0,0))*one;

		//Find corrected kernel
		double rho = 0.0;
		for (size_t j=0; j< outlist[i].size(); ++j) 
		{
			Vector2d Rij = pi.xi-pnp1[outlist[i][j]].xi;
			rho += pnp1[outlist[i][j]].m*W2Kernel(Rij.norm(),svar.H)*
			(Beta(0)+Beta(1)*Rij(0)+Beta(2)*Rij(1));
		}

		pnp1[i].rho = rho;
	}
}

///**************** Integration loop **************
double Newmark_Beta(KD_Tree &mat_index)
{
	vector<Vector2d> xih;
	xih.reserve(svar.npts);
	double errsum = 1.0;
	double logbase = 0.0;
	unsigned int k = 0;
	while (log10(sqrt(errsum/(1.0*svar.npts))) - logbase > -7.0)
	{	
		// cout << "K: " << k << endl;
		Forces(mat_index); /*Guess force at time n+1*/

		/*Previous State for error calc*/
		for (size_t  i=0; i< svar.npts; ++i)
			xih.emplace_back(pnp1[i].xi);

		/*Update the state at time n+1*/
		if (svar.Bcase == 4)
		{
				for(size_t i = 0; i<svar.bound_parts; ++i)
				{
					if(pnp1[i].b == true)
						pnp1[i].xi= pn[i].xi +svar.dt*pnp1[i].v;
				}	
		}
		/*Update the state at time n+1*/
		for (size_t i=0; i <svar.bound_parts; ++i) 
		{	/*Boundary Particles*/
			pnp1[i].rho = pn[i].rho+0.5*svar.dt*(pn[i].Rrho+pnp1[i].Rrho);
			pnp1[i].p = fvar.B*(pow(pnp1[i].rho/fvar.rho0,fvar.gam)-1);
		}
		for (size_t i=svar.bound_parts; i < svar.npts ; ++i )
		{	/*Fluid particles*/
			pnp1[i].v = pn[i].v+svar.dt*((1-svar.gamma)*pn[i].f+svar.gamma*pnp1[i].f);
			pnp1[i].rho = pn[i].rho+svar.dt*((1-svar.gamma)*pn[i].Rrho+svar.gamma*pnp1[i].Rrho);
			pnp1[i].xi = pn[i].xi+svar.dt*pn[i].V+0.5*(svar.dt*svar.dt)*(1-2*svar.beta)*pn[i].f
							+(svar.dt*svar.dt*svar.beta)*pnp1[i].f;
			pnp1[i].p = fvar.B*(pow(pnp1[i].rho/fvar.rho0,fvar.gam)-1);
		}
		
		mat_index.index->buildIndex();
		FindNeighbours(mat_index);

		errsum = 0.0;
		for (size_t i=0; i < pnp1.size(); ++i)
		{
			Vector2d r = pnp1[i].xi-xih[i];
			errsum += r.squaredNorm();
		}

		if(k == 0)
			logbase=log10(sqrt(errsum/(1.0*svar.npts)));

		if(k > svar.subits)
			break;

		++k;
		// cout << log10(sqrt(errsum/(1.0*svar.npts))) - logbase << endl;
	}

	
	/*Find maximum safe timestep*/
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});
	double maxf = maxfi->f.norm();
	double dtf = sqrt(svar.H/maxf);
	double dtcv = svar.H/(fvar.Cs+maxmu);
	svar.dt = 0.5*min(dtf,dtcv);

	
	/*Check if more particles need to be created*/
	if(svar.Bcase == 3)
	{
		switch(svar.Bclosed)
		{
			case 0:
			{
				int refresh = 1;
				for (size_t i = svar.npts-svar.nrefresh; i<svar.npts; ++i)
				{	/*Check that the starting area is clear first...*/
					if(pn[i].xi[1]<svar.Pstep-svar.Box[1])
						refresh = 0;
				}

				if(refresh == 1)
				{	/*...If it is, then check if we've exceeded the max points...*/	
					if (svar.addcount < svar.nmax) 
					{	/*...If we havent, then add points. */
						AddPoints();
					}
					else 
					{	/*...If we have, then check we're sufficiently 
						clear to close the boundary*/
						for (size_t i = svar.npts-svar.nrefresh; i<svar.npts; ++i)
						{
							if(pn[i].xi(1)<svar.H*2)
								refresh = 0;
						}
						if (refresh ==1)
						{
							CloseBoundary();
							svar.Bclosed = 1;
						}
					}
					KD_Tree mat_index(2,pnp1,10);
					mat_index.index->buildIndex();
					FindNeighbours(mat_index);
				}
				break;
			}
			case 1:
				break;
		}
	}
	
	
	//Update the state at time n
	pn = pnp1;

	return log10(sqrt(errsum/(1.0*svar.npts)))-logbase;
}


void InitSPH()
{
	switch (svar.Bcase)
	{
		default:
			cout << "Initialising simulation with " << svar.SimPts << " particles" << endl;
			break;

		case 3:
			cout << "Initialising simulation..." << endl;
			break;
	}
	
	
	//Structure input initialiser
	//Particle(Vector2d x, Vector2d v, Vector2d vh, Vector2d f, float rho, float Rrho, bool bound) :
	Vector2d v = Vector2d::Zero();  
	Vector2d f = Vector2d::Zero(); 
	double rho=fvar.rho0;  
	 
	/*create the boundary particles*/ 	 
	double stepx = svar.Pstep*svar.Bstep;
	double stepy = svar.Pstep*svar.Bstep;
	
	int Ny = ceil(svar.Box(1)/stepy);
	stepy = svar.Box(1)/Ny;
	
	int Nx = ceil(svar.Box(0)/stepx);
	stepx = svar.Box(0)/Nx;
	switch (svar.Bcase) 
	{
		case 0:
		{ /*No boundary*/
			break;
		}
		case 1: /*Rectangle*/
		{
			for(int i = 0; i <= Ny ; ++i) {
				Vector2d xi(-svar.Start(0),i*stepy-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}
			// for(int i = 1; i <Nx ; ++i) {
			// 	Vector2d xi(i*stepx,Box(1));
			// 	particles.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,Bound));	
			// }
			// Vector2d x(stepx-svar.Start(0),(Ny+0.5)*stepy);
			// pn.emplace_back(Particle(x,v,f,rho,fvar.Boundmass,true));
			// x(0) = svar.Box(0) -stepx;
			// pn.emplace_back(Particle(x,v,f,rho,fvar.Boundmass,true));

			// for(int i= Ny; i>0; --i) {
			// 	Vector2d xi(svar.Box(0)-svar.Start(0),i*stepy-svar.Start(1));
			// 	pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true));	
			// }
			for(int i = Nx; i > 0; --i) {
				Vector2d xi(i*stepx-svar.Start(0),-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}
			break;
		}
		case 2: /*Bowl*/
		{
			
			double r= svar.Box(0);
			double dtheta = atan((svar.Pstep*svar.Bstep)/r);
			for (double theta = 0.0; theta < M_PI; theta+=dtheta)
			{
				Vector2d xi(-r*cos(theta),r*(1-sin(theta)));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}
			break;
		}
		case 3:
		{	/*Jet in Crossflow*/
			double holeS = svar.Start(0); /*Distance from Box start to hole*/
			double holeD = svar.Start(1); /*Diameter of hole (or width)*/
			double stepb = (svar.Pstep*svar.Bstep);
			int Nb = ceil((holeS)/stepb);
			stepb = holeS/(1.0*Nb);
			
			for (double x=0.0; x<holeS; x+= stepb)
			{
				Vector2d xi(x,0.0);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}

			for (double y = -stepb; y >= -svar.Box[1]-stepb; y-=stepb)			{
				Vector2d xi(pn.back().xi[0],y);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}

			for (double y = pn.back().xi[1]; y < 0.0 ; y+=stepb)
			{
				Vector2d xi(holeS+holeD,y);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}


			for(double x = holeS+holeD; x<=svar.Box(0); x+=stepb)
			{
				Vector2d xi(x,0.0);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true,false));
			}
			break;
		}
		default: 
		{
			cerr << "Boundary case is not within the design. 0 <= Bcase <= 4." << endl;
			exit(-1);
		}
	}
	
	svar.bound_parts = pn.size();
	
	/*Create the simulation particles*/

	svar.addcount = 0;
	switch(svar.Bcase)
	{
		case 3: 
		{	/*Crossflow case*/
			svar.SimPts = 0;
			svar.npts = pn.size();
			/*Update n+1 before adding sim particles*/
			for (auto p: pn)
				pnp1.emplace_back(p);

			AddPoints();
			break;
		}

		default:
		{
			int which = 2 ;
			if (which ==1) 
			{	/*Hexahedral start*/
				svar.xyPART(0) = 75;
				svar.xyPART(1) = 150;
				svar.SimPts = 75*150;
				stepx = sqrt(2)*svar.Pstep;
				stepy = stepx/2;

				for(int j=0; j< svar.xyPART(1); ++j)
				{
					for( int i=0; i< svar.xyPART(0); ++i) 
					{
						double indent = 0.0;
						if (i != svar.xyPART(0) && (j+1)%2 == 0) indent = stepy;

						double x = /*svar.Start(0)+*/i*stepx+indent;
						double y = /*svar.Start(1)+*/j*stepy;
						Vector2d xi(x,y);		
						pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,false));
					}
				}
				fvar.height0 = /*stepy*(svar.xyPART(1)-1)*/
						(svar.xyPART(1)-1)*svar.Pstep/sqrt(2);
			}
			else if (which == 2)
			{	/*Square lattice start*/
				for( int i=0; i< svar.xyPART(0); ++i) 
				{
					for(int j=0; j< svar.xyPART(1); ++j)
					{				
							Vector2d xi(i*svar.Pstep,j*svar.Pstep);		
							pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,false));
					}
				}
				fvar.height0 = /*stepy*(svar.xyPART(1)-1)*/(svar.xyPART(1)-1)*svar.Pstep;
			}

			for (auto p: pn)
				pnp1.emplace_back(p);

			break;
		}
	}
		
	// svar.SimPts+=10*10;

	svar.npts = pn.size();
	
	if(svar.npts!=svar.bound_parts+svar.SimPts)
	{
		cerr<< "Mismatch of particle count." << endl;
		cerr<< "Particle array size doesn't match defined values." << endl;
		exit(-1);
	}
	cout << "Total Particles: " << svar.npts << endl;

	

}

void write_settings()
{
	std::ofstream fp("Test_Settings.txt", std::ios::out);

  if(fp.is_open()) {
    //fp << "VERSION: " << VERSION_TAG << std::endl << std::endl; //Write version
    fp << "SIMULATION PARAMETERS:" << std::endl; //State the system parameters.
    fp << "\tNumber of frames (" << svar.Nframe <<")" << std::endl;
    fp << "\tParticle Spacing ("<< svar.Pstep << ")" << std::endl;
    fp << "\tParticle Mass ("<< fvar.Simmass << ")" << std::endl;
    fp << "\tReference density ("<< fvar.rho0 << ")" << std::endl;
    fp << "\tSupport Radius ("<< svar.H << ")" << std::endl;
    fp << "\tGravitational strength ("<< 9.81 << ")" << std::endl;
    fp << "\tNumber of boundary points (" << svar.bound_parts << ")" << std::endl;
    fp << "\tNumber of simulation points (" << svar.SimPts << ")" << std::endl;
    fp << "\tIntegrator type (Newmark_Beta)" << std::endl;

    fp.close();
  }
  else {
    cerr << "Error opening the settings output file." << endl;
    exit(-1);
  }
}

void write_frame_data(std::ofstream& fp)
{	
	if (svar.Bcase >0)
	{
			fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bound_parts << ", F=POINT" <<
		    ", STRANDID=1, SOLUTIONTIME=" << svar.t << std::endl;
		  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bound_parts); ++b)
			{
		        fp << b->xi(0) << " " << b->xi(1) << " ";
		        fp << b->v.norm() << " ";
		        fp << b->f.norm() << " ";
		        fp << b->rho << " "  << b->p  << " " << b->Sf.norm() << std::endl;
		  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.SimPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << std::endl;
    unsigned int i=0;
  	for (auto p=std::next(pnp1.begin(),svar.bound_parts); p!=std::next(pnp1.begin(),svar.bound_parts+svar.SimPts); ++p)
	{
		/*if (p->xi!=p->xi || p->v!=p->v || p->f!=p->f) {
			cerr << endl << "Simulation is broken. A value is nan." << endl;
			cerr << "Broken line..." << endl;
			cerr << p->xi(0) << " " << p->xi(1) << " ";
	        cerr << p->v.norm() << " ";
	        cerr << p->f.norm() << " ";
	        cerr << p->rho << " " << p->p << std::endl; 
	        fp.close();
			exit(-1);
		}*/
        fp << p->xi(0) << " " << p->xi(1) << " ";
        fp << p->v.norm() << " ";
        fp << p->f.norm() << " ";
        fp << p->rho << " "  << p->p 
        << " " << p->Sf.norm() << std::endl; 
        ++i;
  	}

  	if (svar.Bcase ==3 && svar.acase == 5 && svar.aircount !=0)
  	{
		fp <<  "ZONE T=\"Ghost Air\"" <<", I=" << svar.aircount << ", F=POINT" <<
	    ", STRANDID=3, SOLUTIONTIME=" << svar.t  << std::endl;
	  	for (auto p=std::next(pnp1.begin(),svar.npts); p!=pnp1.end(); ++p)
		{
	        fp << p->xi(0) << " " << p->xi(1) << " ";
	        fp << p->v.norm() << " ";
	        fp << p->f.norm() << " ";
	        fp << p->rho << " "  << p->p 
	        << " " << p->Sf.norm() << std::endl; 
	        ++i;
	  	}
  	}
}

int main(int argc, char *argv[]) 
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
    double duration;
    double error = 0;

    write_header();

    /******* Define the global simulation parameters ******/
    if (argc > 3) 
	{	/*Check number of input arguments*/
		cout << "WARNING: only two input arguments accepted,\n";
		cout << "1: Input file   2: Output file.\n";
		cout << "Other inputs will be ignored." << endl;

	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "WARNING: No inputs provided.\n";
    	cout << "Program will assume a default set of parameters.";
    	cout << "Output file is \'Test.plt\'" << endl;
    	DefaultInput();
    }
    else if (argc > 1)
    {	/*Get parameters if it has been provided*/
    	GetInput(argv[1]);
    }

    /*Check for output file name*/
	std::ofstream f1;
	if(argc == 3)
	{	/*Open that file if it has*/
		f1.open(argv[2], std::ios::out);
	} 
	else
	{
		cout << "WARNING: output file not provided.\nWill write to Test.plt" << endl;
		f1.open("Test.plt", std::ios::out);
	}
    
	InitSPH();
		
	///********* Tree algorithm stuff ************/
	KD_Tree mat_index(2,pnp1,10);
	mat_index.index->buildIndex();
	outlist.reserve(svar.npts);
	FindNeighbours(mat_index);

	///*** Perform an iteration to populate the vectors *****/
	Forces(mat_index); 

	write_settings();

	///*************** Open simulation files ***************/
	std::ofstream f2("frame.info", std::ios::out);
	
	if (f1.is_open() && f2.is_open())
	{
		f1 << std::fixed << setw(10);
		f2 << std::fixed << setw(10);
		
		cout << std::fixed << setprecision(4);
		/* Write file header defining variable names */
		f1 << "TITLE = \"WCXSPH Output\"" << std::endl;
		f1 << "VARIABLES = \"x (m)\", \"y (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\", \"Aero Force\"" << std::endl;
		write_frame_data(f1);
		
		/*Timing calculation + error sum output*/
		t2 = high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;
		cout << "Frame: " << 0 << "  Sim Time: " << svar.t << "  Compute Time: " 
		<< duration <<"  Error: " << error << endl;
		f2 << "Frame: " << 0 << "  S Time: " << svar.t << "  C Time: " 
			<< duration << "  Error: " << error << 
			" Its: " << 0 << endl; 

		///************************* MAIN LOOP ********************/
		
		for (unsigned int frame = 1; frame<= svar.Nframe; ++frame) 
		{	
			int stepits=0;	
			double stept=0.0;		  
			while (stept<svar.framet) 
			{
			    error = Newmark_Beta(mat_index);
			    svar.t+=svar.dt;
			    stept+=svar.dt;
			    ++stepits;
			}
			
			
			t2= high_resolution_clock::now();
			duration = duration_cast<microseconds>(t2-t1).count()/1e6;

			/*Write each frame info to file (Useful to debug for example)*/
			f2 << "Frame: " << frame << "  S-Time: " << svar.t << "  C-Time: " 
			<< duration <<"  Error: " << error << 
			" Its: " << stepits << endl;  
			if(svar.outframe !=0)
			{
				if (frame % svar.outframe == 0 )
				{	/*Output to console every 20 or so steps*/
				  	cout << "Frame: " << frame << "  Sim Time: " << svar.t-svar.dt << "  Compute Time: " 
				  	<< duration <<"  Error: " << error << endl;
				}
			}
			

			DensityReinit();
			write_frame_data(f1);
		}
		f1.close();
		f2.close();
	}
	else
	{
		cerr << "Error opening frame output file." << endl;
		exit(-1);
	}

	cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.t << " seconds" << endl;
	return 0;
}