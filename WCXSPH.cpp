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
#include <cfloat>
#include <vector>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iterator>
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
	Vector2d Start; 						/*Starting sim box dimensions*/
	unsigned int subits,Nframe; 			/*Timestep values*/
	double dt, t, framet;					/*Simulation + frame times*/
	double tNorm, dtfac;						/*Timestep factor*/
	double beta,gamma;						/*Newmark-Beta Parameters*/
	int Bcase, Bclosed;						/*What boundary shape to take*/
	double H,HSQ, sr; 						/*Support Radius + Search radius*/
	int acase, outform;						/*Aerodynamic force case*/
} SIM;

typedef struct FLUID {
	double rho0; 						/*Resting density*/
	double Simmass, Boundmass;			/*Particle and boundary masses*/
	double alpha,eps,Cs, mu;			/*}*/
	double sig;							/*} Fluid properties*/
	double gam, B; 						/*}*/
	double contangb;					/*Boundary contact angle*/
	//double front, height, height0;		/*Dam Break validation parameters*/
	Vector2d vJet, vInf;
	double Acorrect;					/*Correction factor for aero force*/
}FLUID;

SIM svar;
FLUID fvar;

/*Particle data class*/
typedef class Particle {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Particle(Vector2d X, Vector2d Vi, Vector2d Fi, 
			double Rhoi, double Mi, int bound)	
		{
			xi = X;	v = Vi; f = Fi;
			V(0.0,0.0);	Sf(0.0,0.0); Af(0.0,0.0);
			rho = Rhoi;
			Rrho = 0;
			m = Mi;
			b = bound;
			
		}

		int size() const 
		{	/*For neighbour search, return size of xi vector*/
			return(xi.size());
		}

		double operator[](int a) const 
		{	/*For neighbour search, return index of xi vector*/
			return(xi[a]);
		}
	
		Vector2d xi, v, V, f, Sf, Af;
		double rho, p, Rrho, m;
		int b; //What state is a particle. Boundary, forced particle or unforced
}Particle;
typedef std::vector<Particle> State;
///****** Initialise the particles memory *********/
State pn;	/*Particles at n*/
State pnp1; 	/*Particles at n+1*/

/*Define Neighbour search types*/
typedef std::vector<std::vector<size_t>> outl;
outl outlist;
typedef KDTreeVectorOfVectorsAdaptor<State, double> Sim_Tree;
typedef KDTreeVectorOfVectorsAdaptor<vector<Vector2d>, double> Temp_Tree;
nanoflann::SearchParams params;

struct EVecTraits
{
	typedef double ValueType;

	static constexpr auto kSize = 2;

	static ValueType Get(const Eigen::Vector2d& v, const std::size_t i)
	{return *(&v[0] + i);}

	static void Set(Eigen::Vector2d* const v, const std::size_t i, const ValueType val)
	{*(&v[0][0] + i) = val;} 
};

namespace thinks {
	namespace poisson_disk_sampling {
		template<>
		struct VecTraits<Eigen::Vector2d>
		{
			typedef double ValueType;

			static constexpr auto kSize = 2;

			static ValueType Get(const Eigen::Vector2d& v, const std::size_t i)
			{return *(&v[0] + i);}

			static void Set(Eigen::Vector2d* const v, const std::size_t i, const ValueType val)
			{*(&v[0][0] + i) = val;} 
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

Eigen::Vector2i getIVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::stringstream sline(line);

	Eigen::Vector2i x;
	sline >> x[0]; sline >> x[1]; 
	return x;
}

Eigen::Vector2d getDVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::istringstream sline(line);
	
	Eigen::Vector2d x;
	sline >> x[0]; sline >> x[1];  
	return x;
}

void DefaultInput(void) 
{
	//Timestep Parameters
	svar.framet = 0.1;		/*Frame timestep*/
	svar.Nframe = 2500;		/*Number of output frames*/
	svar.outframe = 50;		/*Terminal output frame interval*/
	svar.outform = 1;		/*Output format*/
	svar.subits = 10;		/*Newmark-Beta iterations*/
	svar.dtfac = 0.3;
	svar.nmax = 2000;		/*Max No particles (if dynamically allocating)*/	
	
	//Simulation window parameters
	svar.xyPART(40,40); 	/*Number of particles in (x,y) directions*/
	svar.Start(0.2,0.2); 	/*Simulation particles start + end coords*/
	svar.Bcase = 1; 		/*Boundary case*/
	svar.Box(3,2); 			/*Boundary dimensions*/
	svar.Pstep = 0.01;		/*Initial particle spacing*/
	svar.Bstep = 0.6; 		/*Boundary factor of particle spacing (dx = Pstep*Bstep)*/
	
	// SPH Parameters
	svar.H =  3.0*svar.Pstep; /*Support Radius*/
	
}

void GetInput(int argc, char **argv)
{
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
    	// cout << argv[1] << endl;
    	std::ifstream in(argv[1]);
	  	if(in.is_open()) 
	  	{	/*Simulation parameters*/
	  		cout << "Input file opened. Reading settings..." << endl;
	  		svar.framet = getDouble(in);
	  		svar.Nframe = getInt(in);
	  		svar.outframe = getInt(in);
	  		svar.outform = getInt(in);
	  		svar.subits = getInt(in);
	  		svar.dtfac = getDouble(in);
	  		svar.nmax = getInt(in);	
	  		svar.xyPART = getIVector(in);
	  		svar.Start = getDVector(in);
	  		svar.Box = getDVector(in);
	  		svar.Pstep = getDouble(in);
	  		svar.Bstep = getDouble(in);
	  		svar.Bcase = getInt(in);
	  		svar.acase = getInt(in);
	  		fvar.vJet = getDVector(in); 
	  		fvar.vInf = getDVector(in);
	  		fvar.Acorrect = getDouble(in);

			in.close();
	  	}
	  	else {
		    cerr << "Error opening the input file." << endl;
		    exit(-1);
	  	}
	}

	/*Get fluid properties from fluid.dat*/
	std::ifstream fluid("fluid.dat");
	if (fluid.is_open())
	{	/*Fluid parameters read*/
		Vector2d nb = getDVector(fluid);
		svar.beta = nb[0];	svar.gamma = nb[1];
		double Hfac = getDouble(fluid); /*End of state read*/
	  	svar.H= Hfac*svar.Pstep;
	  	fvar.alpha = getDouble(fluid);
  		fvar.eps = getDouble(fluid);
  		fvar.contangb = getDouble(fluid);
  		fvar.rho0 = getDouble(fluid);
  		fvar.Cs = getDouble(fluid);
  		fvar.mu = getDouble(fluid);
  		fvar.sig = getDouble(fluid);

  		fluid.close();
	}
	else 
	{
		cerr << "fluid.dat not found. Assuming standard set of parameters (water)." << endl;
		svar.beta = 0.25;	svar.gamma = 0.5;
		svar.H= 2.0*svar.Pstep;
	  	fvar.alpha = 0.1;
  		fvar.eps = 0.05;
  		fvar.contangb = 150.0;
  		fvar.rho0 = 1000.0;
  		fvar.Cs = 100.0;
  		fvar.mu = 1.0;
  		fvar.sig = 0.0728;
	}

  	/*Universal parameters based on input values*/
  	svar.dt = 2E-07; 		/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	svar.HSQ = svar.H*svar.H; 
	svar.sr = 4*svar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
  	svar.SimPts = svar.xyPART[0]*svar.xyPART[1]; /*total sim particles*/
  	svar.aircount = 0;

  	/*Mass from spacing and density*/
	fvar.Simmass = double(fvar.rho0*svar.Pstep*svar.Pstep); 
	fvar.Boundmass = fvar.Simmass*svar.Bcase;
	
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
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
void FindNeighbours(Sim_Tree &NP1_INDEX)
{
	outlist.erase(outlist.begin(),outlist.end());
	double search_radius = svar.sr;
	/*Find neighbour list*/
	for(size_t i=0; i<pnp1.size(); ++i)
	{		
		std::vector<std::pair<size_t,double>> matches; /* Nearest Neighbour Search*/
		
		NP1_INDEX.index->radiusSearch(&pnp1[i].xi[0], search_radius, matches, params);

		std::vector<size_t> temp;
		for (auto &j:matches) 
		{
			temp.emplace_back(j.first); 	
		}
		outlist.emplace_back(temp);		
	}		
}

void AddPoints(double y)
{	
	// cout << "Adding points..." << endl;
	Vector2d v = fvar.vJet;  /*Jet velocity*/
	Vector2d f = Vector2d::Zero(); 
	double rho=fvar.rho0; 
	double jetS = svar.Start(0)+2*svar.H;
	double jetE = svar.Start(0)+svar.Start(1);
	svar.nrefresh = 0;
	

	if (svar.acase == 5)
	{
		/*Create the simulation particles*/
		for( double x = jetS; x<jetE - (jetE-jetS)/2; x+=svar.Pstep) 
		{ /*Do the left set of points*/
			Vector2d xi(x,y);
			int pos = svar.npts+svar.nrefresh;		
			pn.insert(pn.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,2));
			pnp1.insert(pnp1.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,2));
			++svar.SimPts;
			++svar.nrefresh;
		}

		double start2 = pn.back().xi[0]+svar.Pstep;
		for( double x = start2; x<=jetE; x+=svar.Pstep) 
		{	/*Do the right set of points*/
			Vector2d xi(x,y);
			int pos = svar.npts+svar.nrefresh;		
			pn.insert(pn.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,1));
			pnp1.insert(pnp1.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,1));
			++svar.SimPts;
			++svar.nrefresh;
		}
	}
	else
	{
		/*Create the simulation particles*/
		for( double x = jetS; x<jetE - (jetE-jetS)/2; x+=svar.Pstep) 
		{ /*Do the left set of points*/
			Vector2d xi(x,y);
			pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
			pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
			++svar.SimPts;
			++svar.nrefresh;
		}

		double start2 = pn.back().xi[0]+svar.Pstep;
		for( double x = start2; x<=jetE; x+=svar.Pstep) 
		{	/*Do the right set of points*/
			Vector2d xi(x,y);		
			pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,1));
			pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,1));
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
		temp.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
	}
	/*Insert particles at the end of boundary particles*/
	pn.insert(pn.begin()+svar.bound_parts,temp.begin(),temp.end());
	pnp1.insert(pnp1.begin()+svar.bound_parts,temp.begin(),temp.end());
	svar.bound_parts+=temp.size(); /*Adjust counter values*/
	svar.npts +=temp.size();
}



void Ghost_Particles(Sim_Tree &NP1_INDEX, double lam, double numpartdens)
{
/***** Find particles outside of the influence of the Liquid *******/
	
	/*Default to an air particle that is outside the influence of the simulation*/
	for (size_t i=svar.npts; i< svar.npts+svar.aircount; ++i)
		pnp1[i].b = 4; 

	for (size_t i=svar.bound_parts; i< svar.npts; ++i)
	{
		for (size_t j=0; j < outlist[i].size(); ++j) 
		{	/*If it's inside the list of neighbours of the liquid, keep it.*/
			if (pnp1[outlist[i][j]].b == 4) pnp1[outlist[i][j]].b = 3;
		}
	}

	for (size_t i=svar.npts; i< svar.npts+svar.aircount; ++i)
	{	/*If it's still outside the influence of the liquid, then delete it.*/
		if(pnp1[i].b == 4)
		{
			pnp1[i] = pnp1.back();
			pnp1.pop_back();
		}
	}

/****** Find particles with reduced density to create particles around ********/


/********** Create more air particles *********/
	std::vector<Vector2d> temp; /*Temporary storage for air particles*/
	for (size_t i=svar.bound_parts; i< svar.npts; ++i) 
	{	/*Find the surface of fluid particles.*/
		Particle pi = pnp1[i];
		Vector2d SurfC= Vector2d::Zero(); 
		pi.b = 1;

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
			if(pj.b==0) 
	            fac=(1+0.5*cos(M_PI*(fvar.contangb/180)));
	        //cout << lam(svar.H) << endl;
			double sij = 0.5*pow(numpartdens,-2.0)*(fvar.sig/lam)*fac;
			SurfC -= (Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*svar.H));
		}

		if (SurfC.norm()>0.05)
		{
			/*Create particles... */
			namespace pds = thinks::poisson_disk_sampling;
			double radius = svar.Pstep;
			std::array<double,2> xmin = {pn[i].xi[0]-2.0*svar.H, pn[i].xi[1]-2.0*svar.H};
			std::array<double,2> xmax = {pn[i].xi[0]+2.0*svar.H,pn[i].xi[1]+2.0*svar.H};

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
			pnp1.emplace_back(Particle(i,fvar.vInf,f,rho,airmass,1));
			pn.emplace_back(Particle(i,fvar.vInf,f,rho,airmass,1));
		}

		NP1_INDEX.index->buildIndex();
		FindNeighbours(NP1_INDEX);
	}
}

Vector2d AeroForce(Vector2d &Vdiff)
{
	Vector2d Fd = Vector2d::Zero();
	double Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
	double Cd = 0.0;

	if (Re < 3500)
	 	Cd = (1.0+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/(Re+0.000001));
	else 
		Cd = (1+0.197*pow(Re,0.63)+2.6e-4*pow(Re,1.38))*(24.0/(Re+0.0001));

	// cout << "Reynolds: " << Re << " Cd: " << Cd << endl;
	//cout << fvar.Acorrect << endl;
	Fd = fvar.Acorrect*(2*svar.Pstep)*Cd*1.225*Vdiff.normalized()*Vdiff.squaredNorm()*svar.Pstep;
	//Fd[1] = 0.0;
	// cout << "Reynolds: " << Re  << " Cd: " << Cd << " Fd: " << Fd[0] << " " << Fd[1] << endl ;
	return Fd;
}

///**************** RESID calculation **************
void Forces(Sim_Tree &NP1_INDEX)
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
		{	
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
			Ghost_Particles(NP1_INDEX, lam, numpartdens);
		
/*		cout << "Npts: " << svar.npts << " Air Count: " << svar.aircount << endl;
		cout << "Bound Parts: " << svar.bound_parts << " Sim Points: " << svar.SimPts << endl; 
		cout << "Pn size: " << pn.size() << " PnP1 Size: " << pnp1.size() << endl;
		cout << "Outlist Size: " << outlist.size() << endl;*/


/******* LOOP 4 - All simulation points: Calculate forces on the fluid. *********/
		for (size_t i=svar.bound_parts; i< svar.npts + svar.aircount; ++i) 
		{	
			Particle pi = pnp1[i];
			pi.V = pi.v;
			if (svar.Bcase == 3 && svar.acase == 3)
				pi.b = 2;

			double Rrhocontr = 0.0;
			Vector2d contrib= Vector2d::Zero();
			Vector2d visc = Vector2d::Zero();
			Vector2d SurfC= Vector2d::Zero(); 
			Vector2d surfNorm = Vector2d::Zero();

			vector<double> mu;  /*Vector to find largest mu value for CFL stability*/
			mu.emplace_back(0);	/*Avoid dereference of empty vector*/

			for (size_t j=0; j < outlist[i].size(); ++j) 
			{	/* Neighbour list loop. */
				Particle pj = pnp1[outlist[i][j]];

				/*Check if the position is the same, and skip the particle if yes*/
				if(pi.xi == pj.xi)
					continue;

				Vector2d Rij = pj.xi-pi.xi;
				Vector2d Vij = pj.v-pi.v;
				double r = Rij.norm();
				double Kern = W2Kernel(r,svar.H);
				Vector2d Grad = W2GradK(Rij, r,svar.H);
				
				/*Pressure and artificial viscosity - Monaghan (1994) p.400*/
				double rhoij = 0.5*(pi.rho+pj.rho);
				double cbar= 0.5*(sqrt((fvar.B*fvar.gam)/pi.rho)+sqrt((fvar.B*fvar.gam)/pj.rho));
				double vdotr = Vij.dot(Rij);
				double muij= svar.H*vdotr/(r*r+0.01*svar.HSQ);
				mu.emplace_back(muij);
				double pifac = alpha*cbar*muij/rhoij;
		
				if (vdotr > 0) pifac = 0;
				contrib += pj.m*Grad*(pifac - pi.p*pow(pi.rho,-2)-pj.p*pow(pj.rho,-2));

				/*Laminar Viscosity - Morris (2003)*/
				visc -= Vij*(pj.m*fvar.mu)/(pi.rho*pj.rho)
					*(1.0/(r*r+0.01*svar.HSQ))*Rij.dot(Grad);

				/*Surface Tension - Nair & Poeschel (2017)*/
				double fac=1.0;
				if(pj.b==true) fac=(1+0.5*cos(M_PI*(fvar.contangb/180)));
				double sij = 0.5*pow(numpartdens,-2.0)*(fvar.sig/lam)*fac;
				SurfC -= (Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*svar.H));

				surfNorm += (Rij/r)*sij*sin((3.0*M_PI*r)/(4.0*svar.H));
				/* XSPH Influence*/
				pi.V+=eps*(pj.m/rhoij)*Kern*Vij; 

				/*drho/dt*/
				Rrhocontr -= pj.m*(Vij.dot(Grad));	

				if (svar.Bcase == 3 && svar.acase == 3)
				{
					double num = -Rij.dot(fvar.vInf);
					double denom = Rij.norm()*fvar.vInf.norm();
					if (num/denom > 0.98)
						pi.b = 1;
				}

			}
			
			/*Crossflow force*/
			Vector2d Fd= Vector2d::Zero();
			if (svar.Bcase == 3 && pi.xi[1] > svar.Pstep)
			{
				switch(svar.acase)
				{
					case 0: /*No aero force*/
						break;

					case 1:	{  /*All upstream particles at start*/
						if( pi.b == 2)
						{
							Vector2d Vdiff = fvar.vInf - pi.V;
							Fd = AeroForce(Vdiff);
						}
						
						break;
					}
					case 2:	{	/* Surface particles */
						if (SurfC.norm()*pow(svar.Pstep*100,3.7622)/pi.m > 2.5)
						{  /*				^ Need to tune this parameter... */
							Vector2d Vdiff = fvar.vInf - pi.V;
							Fd = AeroForce(Vdiff);
						}
						break;
					}
					case 3: {	/* All upstream particles */
						if(pi.b == 2 /*&& SurfC.norm()*pow(svar.Pstep*100,3.7622)/pi.m > 2.5*/)
						{
							Vector2d Vdiff = fvar.vInf - pi.V;
							Fd = AeroForce(Vdiff);
						}
						break;
					}
					case 4:	{	/* Surface particles proportional to ST*/
						/*Work in progress...*/
						break;
					}
				}
			} 
			else if (svar.Bcase == 3 && pi.xi[1] < svar.Pstep)
			{
				Vector2d Vdiff = fvar.vJet - pi.V;
				double Re = Vdiff.norm()*2*svar.Pstep/fvar.mu;
				double Cd = 0.1*(1+0.197*pow(Re,0.63)+2.6*pow(Re,1.38))*(24.0/(Re+0.0001));

				Fd = Vdiff.normalized()*Cd*(2*svar.Pstep)*fvar.rho0*Vdiff.squaredNorm();
				Fd[1] += 9.81*pi.m;
				//cout << Re << endl;
				//cout << Fd[0] << "  " << Fd[1] << endl;
			}
			

			pi.Rrho = Rrhocontr; /*drho/dt*/
			pi.f= contrib - SurfC*fvar.sig/pi.m + Fd/pi.m;

			pi.Sf = SurfC*fvar.sig/pi.m;
			pi.Af = Fd/pi.m;
			pi.f[1] -= 9.81; /*Add gravity*/

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
double Newmark_Beta(Sim_Tree &NP1_INDEX)
{
	vector<Vector2d> xih;
	xih.reserve(svar.npts);
	double errsum = 1.0;
	double logbase = 0.0;
	unsigned int k = 0;
	while (log10(sqrt(errsum/(double(svar.npts)))) - logbase > -7.0)
	{	
		// cout << "K: " << k << endl;
		Forces(NP1_INDEX); /*Guess force at time n+1*/

		/*Previous State for error calc*/
		xih.erase(xih.begin(),xih.end());
		for (size_t  i=0; i< svar.npts; ++i)
			xih.emplace_back(pnp1[i].xi);

		/*Update the state at time n+1*/
		for (size_t i=0; i <svar.bound_parts; ++i) 
		{	/****** BOUNDARY PARTICLES ***********/
			pnp1[i].rho = pn[i].rho+0.5*svar.dt*(pn[i].Rrho+pnp1[i].Rrho);
			pnp1[i].p = fvar.B*(pow(pnp1[i].rho/fvar.rho0,fvar.gam)-1);
		}

		for (size_t i=svar.bound_parts; i < svar.npts ; ++i )
		{	/****** FLUID PARTICLES ***********/
			pnp1[i].v = pn[i].v+svar.dt*((1-svar.gamma)*pn[i].f+svar.gamma*pnp1[i].f);
			pnp1[i].rho = pn[i].rho+svar.dt*((1-svar.gamma)*pn[i].Rrho+svar.gamma*pnp1[i].Rrho);
			pnp1[i].xi = pn[i].xi+svar.dt*pn[i].V+0.5*(svar.dt*svar.dt)*(1-2*svar.beta)*pn[i].f
							+(svar.dt*svar.dt*svar.beta)*pnp1[i].f;
			pnp1[i].p = fvar.B*(pow(pnp1[i].rho/fvar.rho0,fvar.gam)-1);
		}
		
		/****** UPDATE TREE ***********/
		NP1_INDEX.index->buildIndex();
		FindNeighbours(NP1_INDEX);

		/****** FIND ERROR ***********/
		errsum = 0.0;
		for (size_t i=svar.bound_parts; i < svar.npts; ++i)
		{
			Vector2d r = pnp1[i].xi-xih[i];
			errsum += r.squaredNorm();
		}

		if(k == 0)
			logbase=log10(sqrt(errsum/(double(svar.npts))));

		if(k > svar.subits)
			break;

		// cout << k << "  " << log10(sqrt(errsum/(double(svar.npts)))) - logbase << endl;
		++k;
	}

	
	/*Find maximum safe timestep*/
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});
	double maxf = maxfi->f.norm();
	double dtf = sqrt(svar.H/maxf);
	double dtcv = svar.H/(fvar.Cs+maxmu);
	 
	
	if (std::isinf(maxf))
	{
		cerr << "Forces are quasi-infinite. Stopping..." << endl;
		exit(-1);
	}

/***********************************************************************************/
/**************CODE IS NOW VERY SENSITIVE TO DT. MODIFY WITH CAUTION****************/
/***********************************************************************************/
	svar.dt = svar.dtfac*min(dtf,dtcv); 
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
// cout << "Timestep Params: " << maxf << " " << fvar.Cs + maxmu << " " << dtf << " " << dtcv << endl;
// cout << "New Timestep: " << svar.dt << endl;

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
						AddPoints(pnp1[svar.npts-1].xi[1]-svar.Pstep);
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
					Sim_Tree NP1_INDEX(2,pnp1,10);
					NP1_INDEX.index->buildIndex();
					FindNeighbours(NP1_INDEX);
				}
				break;
			}
			case 1:
				break;
		}
	}
	
	
	/****** UPDATE TIME N ***********/
	pn = pnp1;

	return log10(sqrt(errsum/(double(svar.npts))))-logbase;
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
	 
/************** Create the boundary particles  *****************/ 	 
	double stepx = svar.Pstep*svar.Bstep;
	double stepy = svar.Pstep*svar.Bstep;
	
	int Ny = ceil(svar.Box(1)/stepy);
	stepy = svar.Box(1)/Ny;	/*Find exact step to meet dimensions*/
	
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
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}
/*			
			for(int i = 1; i <Nx ; ++i) {
				Vector2d xi(i*stepx,Box(1));
				particles.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,Bound));	
			}
			Vector2d x(stepx-svar.Start(0),(Ny+0.5)*stepy);
			pn.emplace_back(Particle(x,v,f,rho,fvar.Boundmass,true));
			x(0) = svar.Box(0) -stepx;
			pn.emplace_back(Particle(x,v,f,rho,fvar.Boundmass,true));

			for(int i= Ny; i>0; --i) {
				Vector2d xi(svar.Box(0)-svar.Start(0),i*stepy-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true));	
			}*/
			for(int i = Nx; i > 0; --i) {
				Vector2d xi(i*stepx-svar.Start(0),-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
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
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}
			break;
		}
		case 3:
		{	/*Jet in Crossflow*/
			double holeS = svar.Start(0); /*Distance from Box start to hole*/
			double holeD = svar.Start(1)+4*svar.Pstep; /*Diameter of hole (or width)*/
			double stepb = (svar.Pstep*svar.Bstep);
			int Nb = ceil((holeS)/stepb);
			stepb = holeS/(1.0*Nb);
			
			for (double x=0.0; x<holeS; x+= stepb)
			{
				Vector2d xi(x,0.0);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}

			/*Create a bit of the pipe downward.*/
			for (double y = -stepb; y >= -svar.Box[1]-stepb; y-=stepb)			{
				Vector2d xi(pn.back().xi[0],y);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}

			for (double y = pn.back().xi[1]; y < 0.0 ; y+=stepb)
			{
				Vector2d xi(holeS+holeD,y);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}


			for(double x = holeS+holeD; x<=svar.Box(0); x+=stepb)
			{
				Vector2d xi(x,0.0);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
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
	
/***********  Create the simulation particles  **************/

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

			for (double y = 0.0; y > -svar.Box[1]; y-=svar.Pstep)
				AddPoints(y);

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
						pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
					}
				}
				// fvar.height0 = (svar.xyPART(1)-1)*svar.Pstep/sqrt(2);
			}
			else if (which == 2)
			{	/*Square lattice start*/
				for( int i=0; i< svar.xyPART(0); ++i) 
				{
					for(int j=0; j< svar.xyPART(1); ++j)
					{				
							Vector2d xi(i*svar.Pstep,j*svar.Pstep);		
							pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
					}
				}
				// fvar.height0 = (svar.xyPART(1)-1)*svar.Pstep;
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

void write_research_data(std::ofstream& fp)
{	
	if (svar.Bcase >0)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bound_parts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << std::endl;
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bound_parts); ++b)
		{
	        fp << b->xi(0) << " " << b->xi(1) << " ";
	        fp << b->f.norm() << " ";
	        fp << b->Af[0] << " " << b->Af[1] << " "; 
	        fp << b->Sf[0] << " " << b->Sf[1] << endl; 
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.SimPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << std::endl;
  	for (auto p=std::next(pnp1.begin(),svar.bound_parts); p!=std::next(pnp1.begin(),svar.bound_parts+svar.SimPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " ";
        fp << p->f.norm() << " ";
        fp << p->Af.norm() <<  " "; 
        fp << p->Sf[0] << " " << p->Sf[1] << " ";
        fp << p->b << endl; 
  	}

  	if (svar.Bcase ==3 && svar.acase == 5 && svar.aircount !=0)
  	{
		fp <<  "ZONE T=\"Ghost Air\"" <<", I=" << svar.aircount << ", F=POINT" <<
	    ", STRANDID=3, SOLUTIONTIME=" << svar.t  << std::endl;
	  	for (auto p=std::next(pnp1.begin(),svar.npts); p!=pnp1.end(); ++p)
		{
	        fp << p->xi(0) << " " << p->xi(1) << " ";
	        fp << p->f.norm() << " ";
	        fp << p->Af[0] << " " << p->Af[1] << " "; 
	        fp << p->Sf[0] << " " << p->Sf[1] << "\n";  
	  	}
  	}
}

void write_fluid_data(std::ofstream& fp)
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
		        fp << b->rho << " "  << b->p  << " " << b->Sf << std::endl;
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
        << " " << p->Sf << std::endl; 
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
	        << " " << p->Sf << std::endl; 
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
    cout << std::scientific << setw(6);

    write_header();

    /******* Define the global simulation parameters ******/

    GetInput(argc,argv);
    

    
	std::ofstream f1;
	/*Check for output file name*/
	if(argc == 3)
	{	/*Open that file if it has*/
		f1.open(argv[2], std::ios::out);
	} 
	else
	{	/*Otherwise, open a standard file name*/
		cout << "WARNING: output file not provided.\nWill write to Test.plt" << endl;
		f1.open("Test.plt", std::ios::out);
	}
    
	InitSPH();
		
	///********* Tree algorithm stuff ************/
	Sim_Tree NP1_INDEX(2,pnp1,10);
	NP1_INDEX.index->buildIndex();
	outlist.reserve(svar.npts);
	FindNeighbours(NP1_INDEX);

	///*** Perform an iteration to populate the vectors *****/
	Forces(NP1_INDEX); 

	write_settings();

	///*************** Open simulation files ***************/
	std::ofstream f2("frame.info", std::ios::out);
	std::ofstream f3("Crossflow.txt", std::ios::out);
	if (f1.is_open() && f2.is_open() && f3.is_open())
	{
		f1 << std::scientific << setprecision(5);
		f2 << std::scientific<< setw(10);
		f3 << std::scientific << setw(10);

		
		/* Write file header defining variable names */
		f1 << "TITLE = \"WCXSPH Output\"" << std::endl;
		
		switch (svar.outform)
		{	
			case 1:
				f1 << "VARIABLES = \"x (m)\", \"y (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\", \"Aero Force\"" << std::endl;
				write_fluid_data(f1);
				break;
			case 2:
				f1 << "VARIABLES = \"x (m)\", \"y (m)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"A<sub>f</sub>\", \"S<sub>fx</sub>\", \"S<sub>fy</sub>\", \"B\""<< std::endl;
				write_research_data(f1);
				break;
		}
		
		
		/*Timing calculation + error sum output*/
		t2 = high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;
		cout << "Frame: " << 0 << "  Sim Time: " << svar.t << "  Compute Time: " 
		<< duration <<"  Error: " << error << endl;
		f2 << "Frame:   Pts:    S-Time:    C-Time    Error:   Its:" << endl;
		f2 << 0 << "        " << svar.npts << "    " << svar.t << "    " << duration 
			<< "    " << error << "  " << 0 << endl; 

		///************************* MAIN LOOP ********************/
		
		for (unsigned int frame = 1; frame<= svar.Nframe; ++frame) 
		{	
			int stepits=0;	
			double stept=0.0;		  
			while (stept<svar.framet) 
			{
			    error = Newmark_Beta(NP1_INDEX);
			    svar.t+=svar.dt;
			    stept+=svar.dt;
			    ++stepits;
			    //cout << svar.t << "  " << svar.dt << endl;
			}
			
			
			t2= high_resolution_clock::now();
			duration = duration_cast<microseconds>(t2-t1).count()/1e6;

			/*Write each frame info to file (Useful to debug for example)*/
			f2 << frame << "         " << svar.npts << "  " << svar.t << "  " << duration 
				<< "  " << error << "  " << stepits << endl;  
			if(svar.outframe !=0)
			{
				if (frame % svar.outframe == 0 )
				{	/*Output to console every 20 or so steps*/
				  	cout << "Frame: " << frame << "  Sim Time: " << svar.t-svar.dt << "  Compute Time: " 
				  	<< duration <<"  Error: " << error << endl;
				}
			}
			

			DensityReinit();
			switch (svar.outform)
			{	
				case 1:
					write_fluid_data(f1);
					break;
				case 2:
					write_research_data(f1);
					break;
			}
		
		}
		f1.close();
		f2.close();
		f3.close();
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