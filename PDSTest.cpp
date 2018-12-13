#include <omp.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
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

using namespace std;
using namespace Eigen;

typedef struct SIM {
	Vector2i xyPART; 						/*Starting sim particles in x and y box*/
	unsigned int SimPts,bound_parts,npts;	/*Number of particles*/
	unsigned int nrefresh,nmax;	 			/*Crossflow Particles*/
	double Pstep,Bstep;						/*Initial spacings*/
	Vector2d Box;							/*Box dimensions*/
	Vector2d Start; 				/*Starting sim box dimensions*/
	unsigned int subits,Nframe; 			/*Timestep values*/
	double dt, t, framet;					/*Simulation + frame times*/
	double tNorm;
	double beta,gamma;						/*Newmark-Beta Parameters*/
	int Bcase, Bclosed;						/*What boundary shape to take*/
	double H,HSQ, sr; 						/*Support Radius + Search radius*/
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

typedef std::vector<std::vector<size_t>> outl;
outl outlist;
typedef KDTreeVectorOfVectorsAdaptor<vector<Vector2d>, double> KD_Tree;
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




void Initialise()
{
	svar.framet = 0.1;		/*Frame timestep*/
	svar.Nframe = 2500;		/*Number of output frames*/
	svar.subits = 10;			/*Newmark-Beta iterations*/
	svar.nmax = 2000;			/*Max No particles (if dynamically allocating)*/	
	svar.beta = 0.25;			/*Newmark-Beta parameters*/
	svar.gamma = 0.5;	
	
	//Simulation window parameters
	svar.xyPART= {6,6}; /*Number of particles in (x,y) directions*/
	svar.Start = {0.2,0.2}; /*Simulation particles start + end coords*/
	svar.Bcase = 1; 		/*Boundary case*/
	svar.Box = {3,2}; 		/*Boundary dimensions*/
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
	fvar.vJet[1] = 10.0;	

	/*Universal parameters based on input values*/
  	svar.dt = 0.00002; 		/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	svar.HSQ = svar.H*svar.H; 
	svar.sr = 4*svar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
  	svar.SimPts = svar.xyPART(0)*svar.xyPART(1); /*total sim particles*/

  	Vector2d Finish(svar.Start(0)+1.0*svar.Pstep*svar.xyPART(0),
					svar.Start(1)+1.0*svar.Pstep*svar.xyPART(1));
	
	fvar.Simmass = fvar.rho0* /*Mass from spacing and density*/
  		((Finish(0)-svar.Start(0))*(Finish(1)-svar.Start(1)))/(1.0*svar.SimPts);
	fvar.Boundmass = fvar.Simmass*svar.Bcase;
	
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
}


void InitSPH(void)
{
	//Structure input initialiser
	//Particle(Vector2d x, Vector2d v, Vector2d vh, Vector2d f, float rho, float Rrho, bool bound) :
	Vector2d v = Vector2d::Zero();  
	Vector2d f = Vector2d::Zero(); 
	double rho=fvar.rho0;  
	
	/*Square lattice start*/
	cout << "Particle coords: " << endl;
	for( int i=0; i< svar.xyPART[0]; ++i) 
	{
		for(int j=0; j< svar.xyPART[1]; ++j)
		{				
				Vector2d xi(i*svar.Pstep,j*svar.Pstep);		
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,false,false));
				cout << xi[0] << " " << xi[1] << endl; 
		}
	}
	svar.bound_parts = pn.size();
}

void Foo(void)
{
	namespace pds = thinks::poisson_disk_sampling;
	Vector2d v = Vector2d::Zero();  
	Vector2d f = Vector2d::Zero(); 
	double rho=fvar.rho0; 
	double radius = svar.Pstep;
	svar.SimPts=0;
	std::vector<Vector2d> temp;

	for (size_t i = 0; i < svar.bound_parts; ++i )
	{
		if ((pn[i].xi[0]== 0.0 || pn[i].xi[0]== 1.0*(svar.xyPART[0]-1)*svar.Pstep) ||
			(pn[i].xi[1]== 0.0 || pn[i].xi[1]== 1.0*(svar.xyPART[1]-1)*svar.Pstep))
		{
			std::array<double,2> xmin = {pn[i].xi[0]-2.0*svar.H, pn[i].xi[1]-2.0*svar.H};
			std::array<double,2> xmax = {pn[i].xi[0]+2.0*svar.H,pn[i].xi[1]+2.0*svar.H};
			cout << "Centre Coord: " << pn[i].xi[0] << " " << pn[i].xi[1] << endl;
			cout << "Min Coords: " << xmin[0] << "  " << xmin[1] << endl;
			cout << "Max Coords: " << xmax[0] << "  " << xmax[1] << endl;

			std::vector<Vector2d> samples = 
			pds::PoissonDiskSampling<double,2,Vector2d,EVecTraits>(radius,xmin,xmax);
			
			cout << samples.size() << endl;
			for (auto j:samples)
			{
				Vector2d xi(j[0],j[1]); 
				temp.emplace_back(xi);
			}
		}
	}

	KD_Tree mat_index(2,temp,10);
	mat_index.index->buildIndex();
	double search_radius = svar.Pstep*svar.Pstep;
	//std::vector<size_t> delete_list;
	cout << temp.size() << endl;

	for (auto i=temp.begin(); i!=temp.end(); )
	{
		Vector2d xi = *i;
		std::vector<std::pair<size_t,double>> matches;
		mat_index.index->radiusSearch(&xi[0],search_radius,matches,params);
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
			mat_index.index->buildIndex();
		}
		else
			++i;

	}
	cout << temp.size() << endl;

	for (size_t i = 0; i < svar.bound_parts; ++i )
	{
		std::vector<std::pair<size_t,double>> matches;
		mat_index.index->radiusSearch(&pn[i].xi[0],search_radius,matches,params);
		//cout << matches.size()<< endl;
		if (matches.size()!=0)
		{
			for (auto j:matches)
			{
				temp[j.first] = temp.back();
				temp.pop_back();
			}
			mat_index.index->buildIndex();
		}
	}
	svar.SimPts = temp.size();
	// std::sort(delete_list.begin(),delete_list.end());
	// auto last = std::unique(delete_list.begin(),delete_list.end());
	// delete_list.erase(last,delete_list.end());

	
		

	//FindNeighbours(mat_index);

	for (auto i:temp)
		pn.emplace_back(Particle(i,v,f,rho,fvar.Simmass,false,false)); 
	

}

void write_frame_data(std::ofstream& fp)
{	
	if (svar.Bcase >0)
	{
			fp <<  "ZONE T=\"Particle Data\"" << ", I=" << svar.bound_parts << ", F=POINT" << std::endl;
		  	for (auto b=pn.begin(); b!=std::next(pn.begin(),svar.bound_parts); ++b)
			{
		        fp << b->xi(0) << " " << b->xi(1) << " ";
		        fp << b->v.norm() << " ";
		        fp << b->f.norm() << " ";
		        fp << b->rho << " "  << b->p  << " " << b->Sf.norm() << std::endl;
		  	}
	}
    
    fp <<  "ZONE T=\"Ghost Air\"" <<", I=" << svar.SimPts << ", F=POINT" << std::endl;
    unsigned int i=0;
  	for (auto p=std::next(pn.begin(),svar.bound_parts); p!=pn.end(); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " ";
        fp << p->v.norm() << " ";
        fp << p->f.norm() << " ";
        fp << p->rho << " "  << p->p 
        << " " << p->Sf.norm() << std::endl; 
        ++i;
  	}
}


int main()
{
	Initialise();

	InitSPH();

	Foo();

	std::ofstream f1("PDS_Test.plt", std::ios::out);

	if (f1.is_open())
	{
		f1 << std::fixed << setw(10);
		/* Write file header defining variable names */
		f1 << "TITLE = \"WCXSPH Output\"" << std::endl;
		f1 << "VARIABLES = \"x (m)\", \"y (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\", \"SurfC\"" << std::endl;


		write_frame_data(f1);
	}
	else
	{
		cerr << "Error opening frame output file." << endl;
		exit(-1);
	}

	cout << "Simulation complete!" << endl;

	return 0;
}