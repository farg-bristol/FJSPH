/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <dirent.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include "../Eigen/Core"
#include "../Eigen/StdVector"
#include "../Eigen/LU"
// #include "TECIO.h"

// #include <vtkXMLUnstructuredGrid.h>
// #include <vtkXMLUnstructuredGridWriter.h>
// #include <vtkSmartPointer.h>
// #include <vtkPoints.h>
// #include <vtkMarchingCubes.h>

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

#include "Var.h"

#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using namespace std; 

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

void write_header() 
{
	cout << "******************************************************************" << endl << endl;
	cout << "                              WCSPH                               " << endl << endl;
	cout << "        Weakly Compressible Smoothed Particle Hydrodynamics       " << endl;
	cout << "                      for Fuel Jettison case                      " << endl << endl;
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

int MakeOutputDir(int argc, char *argv[], SIM &svar, std::ofstream& f1)
{
	char cCurrentPath[FILENAME_MAX];
	struct stat info;
	  	
  	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	return errno;

	/*Open an output directory in the name of the input file, under Outputs*/
  	string pathname = cCurrentPath;
  	string input = argv[1];
  	uint pos1 = input.find("/"); 	uint pos2 = input.find(".");
  	string file = input.substr(pos1,pos2-pos1);
  	// cout << file << endl;
  	pathname.append("/Outputs");
  	pathname.append(file);
  	svar.outfolder = pathname;
	int sreturn=0;

  	if( stat( pathname.c_str(), &info ) != 0 )
  	{	/*Output directory doesn't exist, so create it*/
  		string cmd = "mkdir ";
	  	cmd.append(pathname);
	    sreturn = system(cmd.c_str());
	    if(sreturn == -1)
	    {
	    	cout << "System command failed to execute." << endl;
	    	exit(-1);
	    }
  	}
	else if( info.st_mode & S_IFDIR ) 
	{	/*If it exists, check if there are contents inside directory */
		DIR *dir;
		struct dirent *ent;
		if ((dir = opendir (pathname.c_str())) != NULL) 
		{
			ent = readdir (dir);
			ent = readdir (dir);
			if((ent = readdir (dir)) != NULL) 
			{ /* Delete contents inside the directory */
			    cout << "Deleting files in output directory..." << endl;
				
			    string cmd = "exec rm ";
			    string fuel = pathname;
				fuel.append("/Fuel.plt");
				
				string vlmpanel = pathname;
				vlmpanel.append("/VLM_Panels.plt");
				
				string vlmvortices = pathname;
				vlmvortices.append("/VLM_Vortices.plt");
				

				if(file_exists(fuel))
				{
					fuel.insert(0,cmd);
					sreturn = system(fuel.c_str());
				}

				if(file_exists(vlmpanel))
				{
					vlmpanel.insert(0,cmd);
					sreturn = system(vlmpanel.c_str());
				}

				if(file_exists(vlmvortices))
				{
					vlmvortices.insert(0,cmd);
					sreturn = system(vlmvortices.c_str());
				}
				
				
				if(sreturn == -1)
			    {
			    	cout << "System command failed to execute." << endl;
			    	exit(-1);
			    }
		  }
		  closedir (dir);
		}
	}
	else
	{
	    cerr << "Can't access or create output directory. Stopping." << endl;
	    exit(-1);
	}

	/*Check for output file name*/
	if(argc == 2)
	{
		string path = pathname;
		path.append("/Fuel.plt");
		f1.open(path, std::ios::out);
	}
	else if(argc == 3)
	{	/*Open that file if it has*/
		f1.open(argv[2], std::ios::out);
	}
	else
	{	/*Otherwise, open a standard file name*/
		cout << "\tWARNING: output file not provided.\nWill write to Test.plt" << endl;
		f1.open("Test.plt", std::ios::out);
	}

	return 0;
}

StateVecI getIVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::istringstream sline(line);
	// cout << sline.str() << endl;
	StateVecI x;
	sline >> x[0]; sline >> x[1]; sline >> x(2);

	// cout << sline.str() << endl;

	if (!sline)
	{	
		cout << "2D input provided. Please provide a 3D file." << endl;
		cout << "Incorrect line:" << endl;
		cout << sline.str() << endl;
		exit(-1);
	}

	return x;
}

StateVecD getDVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::istringstream sline(line);
	
	StateVecD x;
	sline >> x[0]; sline >> x[1]; sline >> x[2]; 

	if (!sline)
	{	
		cout << "2D input provided. Please provide a 3D file." << endl;
		cout << "Incorrect line:" << endl;
		cout << sline.str() << endl;
		exit(-1);
	}
		
	return x;
}

/*Function for a 2D Vector (e.g. Newmark Beta parameters)*/
StateVecD getvector(ifstream& In)
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

void GetAero(FLUID &fvar, ldouble rad)
{
	fvar.avar.L = rad * std::cbrt(3.0/(4.0*M_PI));
	fvar.avar.td = (2.0*fvar.rho0*fvar.avar.L*fvar.avar.L)/(fvar.avar.Cd*fvar.mu);
	fvar.avar.omega = sqrt((fvar.avar.Ck*fvar.sig)/(fvar.rho0*pow(fvar.avar.L,3.0))-1.0/pow(fvar.avar.td,2.0));

	fvar.avar.tmax = -2.0 *(atan(sqrt(pow(fvar.avar.td*fvar.avar.omega,2.0)+1)
					+fvar.avar.td*fvar.avar.omega) - M_PI)/fvar.avar.omega;

	fvar.avar.Cdef = 1.0 - exp(-fvar.avar.tmax/fvar.avar.td)*(cos(fvar.avar.omega*fvar.avar.tmax)+
		1/(fvar.avar.omega*fvar.avar.td)*sin(fvar.avar.omega*fvar.avar.tmax));
	fvar.avar.ycoef = 0.5*fvar.avar.Cdef*(fvar.avar.Cf/(fvar.avar.Ck*fvar.avar.Cb))*(fvar.rhog*fvar.avar.L)/fvar.sig;

	//cout << fvar.avar.ycoef << "  " << fvar.avar.Cdef << "  " << fvar.avar.tmax << "  " << endl;
}

RotMat GetRotationMat(StateVecD &angles)
{
	RotMat rotx, roty, rotz;
	rotx << 1.0, 0.0            , 0.0           ,
		    0.0, cos(angles(2)) , sin(angles(2)),
			0.0, -sin(angles(2)), cos(angles(2));

	roty << cos(angles(1)) , 0.0 , -sin(angles(1)),
		    0.0            , 1.0 , 0.0            ,
			sin(angles(1)) , 0.0 , cos(angles(1));

	rotz << cos(angles(0)) , sin(angles(0)) , 0.0 ,
		    -sin(angles(0)), cos(angles(0)) , 0.0 ,
			0.0            , 0.0            , 1.0 ;

	return rotx*roty*rotz;

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
    	cout << "Program will assume a default set of parameters.\n";
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
	  		svar.Pstep = getDouble(in);
	  		svar.Bstep = getDouble(in);
	  		svar.Bcase = getInt(in);
	  		cvar.acase = getInt(in);
	  		svar.Start = getDVector(in);
	  		if(svar.Bcase < 3)
	  		{
	  			svar.Box= getDVector(in);
	  		}
	  		if(svar.Bcase >= 3)
	  		{	
	  			StateVecD angles = getDVector(in);
	  			angles = angles *M_PI/180;
	  			svar.Rotate = GetRotationMat(angles);
		  		svar.Jet = getvec(in); /*Defined in VLM.h. Reused here*/
		  		fvar.pPress = getDouble(in);
		  		cvar.vJet = StateVecD::Zero(); cvar.vInf = StateVecD::Zero();
		  		cvar.vJet(1) = getDouble(in);  cvar.vInf(1) = getDouble(in);
		  		cvar.vJet = svar.Rotate*cvar.vJet;
		  		cvar.Acorrect = getDouble(in);
		  		if(cvar.acase >= 2)
		  		{
		  			cvar.a = getDouble(in);
		  			cvar.h1 = getDouble(in);
		  			cvar.b = getDouble(in);
		  			cvar.h2 = getDouble(in);
		  		}
		  		
		  		if(cvar.acase > 5)
		  		{
		  			cout << "Aerodynamic case is not in design. Stopping..." << endl;
		  			exit(-1);
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
		StateVecD nb = getvector(fluid);
		svar.beta = nb[0];	svar.gamma = nb[1];
		double Hfac = getDouble(fluid); /*End of state read*/
	  	fvar.H= Hfac*svar.Pstep;
	  	fvar.alpha = getDouble(fluid);
  		fvar.contangb = getDouble(fluid);
  		fvar.rho0 = getDouble(fluid);
  		fvar.rhog = getDouble(fluid);
  		fvar.Cs = getDouble(fluid);
  		fvar.mu = getDouble(fluid);
  		fvar.mug = getDouble(fluid);
  		fvar.sig = getDouble(fluid);

  		fluid.close();
	}
	else 
	{
		cerr << "fluid.dat not found. Assuming standard set of parameters (water)." << endl;
		svar.beta = 0.25;	svar.gamma = 0.5;
		fvar.H= 2.0*svar.Pstep;
	  	fvar.alpha = 0.1;
  		fvar.contangb = 150.0;
  		fvar.rho0 = 1000.0;
  		fvar.Cs = 100.0;
  		fvar.mu = 1.0;
  		fvar.sig = 0.0728;
  		fvar.rhog = 1.225;
  		fvar.mug = 18.5E-06;
	}

  	/*Universal parameters based on input values*/
  	svar.addcount = 0;
  	svar.dt = 2E-010; 		/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	fvar.HSQ = fvar.H*fvar.H; 
  	
  	if(simDim == 2)
  		fvar.correc = (7/(4*M_PI*fvar.H*fvar.H));
  	else if(simDim == 3)
  		fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
  	else
  	{
  		cout << "Simulation Dimension mismatch. Stopping" << endl;
  		exit(-1);
  	}

	fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
  	svar.simPts = svar.xyPART[0]*svar.xyPART[1]*svar.xyPART[2]; /*total sim particles*/

  	/*Mass from spacing and density*/
	fvar.Simmass = fvar.rho0*pow(svar.Pstep,simDim); 
	fvar.Boundmass = fvar.Simmass;
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
	ldouble rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
	svar.dx = pow(fvar.Simmass/rho, 1.0/double(simDim));
	cout << rho << "  " << svar.dx << endl;
	svar.vortex.GetGamma(cvar.vInf);
	GetAero(fvar, fvar.H);
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
	if (svar.Bcase >0 && svar.Bcase != 5)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << "\n";
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
		{
	        fp << b->xi(0) << " " << b->xi(1) << " " << b->xi(2) << " ";
	        fp << b->f.norm() << " ";
	        fp << b->Af.norm() << " " << b->Sf.norm() << " "; 
	        fp << b->Af[0] << " " << b->Af[1] << " " << b->Af[2] << " ";
	        fp << b->b << " " << b->theta << "\n"; 
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
  	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " " << p->xi(2) << " ";
        fp << p->f.norm() << " ";
        fp << p->Af[0] << " " << p->Sf.norm() << " "; 
        fp << p->Af[0] << " " << p->Af[1] << " " << p->Af(2) << " ";
        fp << p->b << " " << p->theta  << "\n"; 
  	}
  	fp << std::flush;
}

void write_fluid_data(std::ofstream& fp, SIM &svar, State &pnp1)
{	
	if (svar.Bcase >0 && svar.Bcase != 5)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << "\n";
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
		{
	        fp << b->xi[0] << " " << b->xi[1] << " " << b->xi[2] << " ";
	        fp << b->v.norm() << " ";
	        fp << b->f.norm() << " ";
	        fp << b->rho << " "  << b->p << "\n";
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
  	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " " << p->xi[2] << " ";
        fp << p->v.norm() << " ";
        fp << p->f.norm() << " ";
        fp << p->rho << " "  << p->p  << "\n";  
  	}  	
  	fp << std::flush;
}

void write_basic_data(std::ofstream& fp, SIM &svar, State &pnp1)
{	
	if (svar.Bcase >0 && svar.Bcase != 5)
	{
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" <<
	    ", STRANDID=1, SOLUTIONTIME=" << svar.t << "\n";
	  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
		{
	        fp << b->xi[0] << " " << b->xi[1] << " " << b->xi[2] << std::endl;
	  	}
	}
    
    fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
  	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
	{
        fp << p->xi(0) << " " << p->xi(1) << " " << p->xi[2] << "\n";  
  	}
  	fp << std::flush;
}


// void write_VTK_data(char* file, SIM &svar, State &pnp1)
// {	/*File will be the name of the input file*/
// 	string string;
// 	string = file;
// 	string.append("");

// 	std::ofstream fp("VLM_Panels.plt", std::ios::out);

// }

void write_file_header(std::ofstream& fp, SIM &svar, State &pnp1)
{
	switch (svar.outform)
		{	
			case 1:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\"" << "\n";
				write_fluid_data(fp, svar, pnp1);
				break;
			case 2:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"A<sub>f</sub>\", \"S<sub>f</sub>\", \"S<sub>fx</sub>\", \"S<sub>fy</sub>\", \"S<sub>fz</sub>\", \"B\", \"Theta\"" << "\n";
				write_research_data(fp, svar, pnp1);	
				break;
			case 3:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\""<< "\n";
				write_basic_data(fp, svar, pnp1);
				break;
		}
}

#endif