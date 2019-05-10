/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include "../Eigen/Core"
#include "../Eigen/StdVector"
#include "../Eigen/LU"

#include "Var.h"

#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using namespace std; 


void write_header() 
{
	cout << "******************************************************************" << endl << endl;
	cout << "                              WCXSPH                              " << endl << endl;
	cout << "        Weakly Compressible Smoothed Particle Hydrodynamics       " << endl;
	cout << "                       with XSPH correction                       " << endl;
	cout << "                    for Jet in Crossflow case                     " << endl << endl;
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

StateVecI getIVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::stringstream sline(line);

	StateVecI x;
	sline >> x[0]; sline >> x[1]; 
	return x;
}

StateVecD getDVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::istringstream sline(line);
	
	StateVecD x;
	sline >> x[0]; sline >> x[1];  
	return x;
}

void DefaultInput(SIM &svar) 
{
	//Timestep Parameters
	svar.framet = 0.1;		/*Frame timestep*/
	svar.Nframe = 2500;		/*Number of output frames*/
	svar.outframe = 50;		/*Terminal output frame interval*/
	svar.outform = 1;		/*Output format*/
	svar.frameout = 1;		/*Output frame info*/
	svar.subits = 10;		/*Newmark-Beta iterations*/
	svar.nmax = 2000;		/*Max No particles (if dynamically allocating)*/	
	
	//Simulation window parameters
	svar.xyPART[0] = 40; svar.xyPART[1]= 40; 	/*Number of particles in (x,y) directions*/
	svar.Start[0] = 0.2; svar.Start[1]= 0.2; 	/*Simulation particles start + end coords*/
	svar.Box[0] = 3; svar.Box[1]= 2; 			/*Boundary dimensions*/
	svar.Pstep = 0.01;		/*Initial particle spacing*/
	svar.Bstep = 0.6; 		/*Boundary factor of particle spacing (dx = Pstep*Bstep)*/
	svar.Bcase = 1; 		/*Boundary case - Rectangle */
	
}

void GetInput(int argc, char **argv, SIM &svar, FLUID &fvar, CROSS &cvar)
{
	if (argc > 3) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only two input arguments accepted,\n";
		cout << "1: Input file   2: Output file.\n";
		cout << "Other inputs will be ignored." << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tWARNING: No inputs provided.\n";
    	cout << "Program will assume a default set of parameters.";
    	cout << "Output file is \'Test.plt\'" << endl;
    	DefaultInput(svar);
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
	  		svar.frameout = getInt(in);
	  		svar.subits = getInt(in);
	  		svar.nmax = getInt(in);	
	  		svar.xyPART = getIVector(in);
	  		svar.Start = getDVector(in);
	  		svar.Box = getDVector(in);
	  		svar.Pstep = getDouble(in);
	  		svar.Bstep = getDouble(in);
	  		svar.Bcase = getInt(in);
	  		if(svar.Bcase == 3)
	  		{
	  			cvar.acase = getInt(in);
		  		cvar.vJet = getDVector(in); 
		  		cvar.vInf = getDVector(in);
		  		cvar.Acorrect = getDouble(in);
		  		if(cvar.acase >= 3 )
		  		{
		  			cvar.a = getDouble(in);
		  			cvar.h1 = getDouble(in);
		  			cvar.b = getDouble(in);
		  			cvar.h2 = getDouble(in);
		  			if(cvar.acase == 5)
		  			{
		  				cvar.vortexPos = getDVector(in);
		  			}
		  		}
	  		}
	  		
	  		
			in.close();
	  	}
	  	else {
		    cerr << "Error opening the input file." << endl;
		    exit(-1);
	  	}
	}

	/*Get fluid properties from fluid.dat*/
	std::ifstream fluid("Fluid.dat");
	if (fluid.is_open())
	{	/*Fluid parameters read*/
		StateVecD nb = getDVector(fluid);
		svar.beta = nb[0];	svar.gamma = nb[1];
		double Hfac = getDouble(fluid); /*End of state read*/
	  	fvar.H= Hfac*svar.Pstep;
	  	fvar.alpha = getDouble(fluid);
  		fvar.eps = getDouble(fluid);
  		fvar.contangb = getDouble(fluid);
  		fvar.rho0 = getDouble(fluid);
  		cvar.rhog = getDouble(fluid);
  		fvar.Cs = getDouble(fluid);
  		fvar.mu = getDouble(fluid);
  		cvar.mug = getDouble(fluid);
  		fvar.sig = getDouble(fluid);

  		fluid.close();
	}
	else 
	{
		cerr << "fluid.dat not found. Assuming standard set of parameters (water)." << endl;
		svar.beta = 0.25;	svar.gamma = 0.5;
		fvar.H= 2.0*svar.Pstep;
	  	fvar.alpha = 0.1;
  		fvar.eps = 0.05;
  		fvar.contangb = 150.0;
  		fvar.rho0 = 1000.0;
  		fvar.Cs = 100.0;
  		fvar.mu = 1.0;
  		fvar.sig = 0.0728;
  		cvar.rhog = 1.225;
  		cvar.mug = 18.5E-06;
	}

  	/*Universal parameters based on input values*/
  	svar.addcount = 0;
  	svar.dt = 2E-07; 		/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	fvar.HSQ = fvar.H*fvar.H; 
  	fvar.correc = (7/(4*M_PI*fvar.H*fvar.H));
	fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
  	svar.simPts = svar.xyPART[0]*svar.xyPART[1]; /*total sim particles*/
  	svar.aircount = 0;

  	/*Mass from spacing and density*/
	fvar.Simmass = double(fvar.rho0*svar.Pstep*svar.Pstep); 
	fvar.Boundmass = fvar.Simmass*svar.Bcase;
	
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
}

/******************* OUTPUTS *********************/
void write_settings(SIM &svar, FLUID &fvar)
{
	std::ofstream fp("Test_Settings.txt", std::ios::out);

  if(fp.is_open()) {
    //fp << "VERSION: " << VERSION_TAG << std::endl << std::endl; //Write version
    fp << "SIMULATION PARAMETERS:" << std::endl; //State the system parameters.
    fp << "\tNumber of frames (" << svar.Nframe <<")" << std::endl;
    fp << "\tParticle Spacing ("<< svar.Pstep << ")" << std::endl;
    fp << "\tParticle Mass ("<< fvar.Simmass << ")" << std::endl;
    fp << "\tReference density ("<< fvar.rho0 << ")" << std::endl;
    fp << "\tSupport Radius ("<< fvar.H << ")" << std::endl;
    fp << "\tGravitational strength ("<< 9.81 << ")" << std::endl;
    fp << "\tNumber of boundary points (" << svar.bndPts << ")" << std::endl;
    fp << "\tNumber of simulation points (" << svar.simPts << ")" << std::endl;
    fp << "\tIntegrator type (Newmark_Beta)" << std::endl;

    fp.close();
  }
  else {
    cerr << "Error opening the settings output file." << endl;
    exit(-1);
  }
}

void write_research_data(std::ofstream& fp, SIM &svar, State &pnp1)
{	
	if (svar.Bcase >0)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << std::endl;
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
		{
	        fp << b->xi(0) << " " << b->xi(1) << " ";
	        fp << b->f.norm() << " ";
	        fp << b->Af.norm() << " " << b->Sf.norm() << " "; 
	        fp << b->Sf[0] << " " << b->Sf[1] << " ";
	        fp << b->b << " " << b->theta << endl; 
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << std::endl;
  	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " ";
        fp << p->f.norm() << " ";
        fp << p->Af[0] << " " << p->Sf.norm() << " "; 
        fp << p->Sf[0] << " " << p->Sf[1] << " ";
        fp << p->b << " " << p->theta  << endl; 
  	}
}

void write_fluid_data(std::ofstream& fp, SIM &svar, State &pnp1)
{	
	if (svar.Bcase >0)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << std::endl;
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
		{
	        fp << b->xi[0] << " " << b->xi[1] << " ";
	        fp << b->v.norm() << " ";
	        fp << b->f.norm() << " ";
	        fp << b->rho << " "  << b->p << std::endl;
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << std::endl;
  	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " ";
        fp << p->v.norm() << " ";
        fp << p->f.norm() << " ";
        fp << p->rho << " "  << p->p  << std::endl;  
  	}
}

void write_basic_data(std::ofstream& fp, SIM &svar, State &pnp1)
{	
	if (svar.Bcase >0)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << std::endl;
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
		{
	        fp << b->xi[0] << " " << b->xi[1] <<  std::endl;
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << std::endl;
  	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << std::endl;  
  	}

}

void write_file_header(std::ofstream& fp, SIM &svar, State &pnp1)
{
	switch (svar.outform)
		{	
			case 1:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\"" << std::endl;
				write_fluid_data(fp, svar, pnp1);
				break;
			case 2:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"A<sub>f</sub>\", \"S<sub>f</sub>\", \"S<sub>fx</sub>\", \"S<sub>fy</sub>\", \"B\", \"Theta\""<< std::endl;
				write_research_data(fp, svar, pnp1);
				break;
			case 3:
				fp << "VARIABLES = \"x (m)\", \"y (m)\""<< std::endl;
				write_basic_data(fp, svar, pnp1);
				break;
		}
}

#endif