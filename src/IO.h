/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <dirent.h>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include "Eigen/LU"
#include "Var.h"
#include "BinaryIO.h"
#include "CDFIO.h"
// #include "TauIO.h"

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string; 
using std::setw;

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/


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

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

void check_folder(string pathname)
{	
#ifdef DEBUG
	dbout << "Checking existence of folder:" << pathname << endl;
#endif

	struct stat info;
	if( stat( pathname.c_str(), &info ) != 0 )
  	{	/*Output directory doesn't exist, so create it*/
		pathname.insert(0,"\"");
		pathname.append("\"");
  		string cmd = "mkdir ";
	  	cmd.append(pathname);
#ifdef DEBUG
		dbout << "Folder doesn't exist. Trying to create folder. Command:" << endl;
		dbout << cmd << endl;
#endif
	    if(system(cmd.c_str()))
	    {
	    	cout << "System command failed to execute." << endl;
	    	cout << "Command: " << cmd << endl;
	    	exit(-1);
	    }
  	}
	else if( info.st_mode & S_IFDIR ) 
	{	/*If it exists, Check that directory can be accessed*/
		DIR *dir;
		if ((dir = opendir (pathname.c_str())) != NULL) 
		    closedir (dir);
	}
	else
	{
	    cerr << "Can't access or create output directory. Stopping." << endl;
	    exit(-1);
	}
}

int MakeOutputDir(int argc, char *argv[], SIM &svar)
{
	char cCurrentPath[FILENAME_MAX];
	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	return errno;

	/*Open an output directory in the name of the input file, under Outputs*/
	string pathname = cCurrentPath;
  	string input = argv[1];
  	uint pos1 = input.find("/"); 	uint pos2 = input.find(".");
  	string file = input.substr(pos1+1,pos2-pos1-1);
  	pathname.append("/Outputs/");

  
  	/*Check for output file name*/
	if(argc == 2)
	{	
		pathname.append(file);
		string path = pathname;
		check_folder(pathname);
	}
	else if(argc == 3)
	{	/*Open that file if it has*/
		pathname.append(argv[2]);
		string path = pathname;
		check_folder(pathname);
	}
	else
	{	/*Otherwise, open a standard file name*/
		cout << "\tWARNING: no files provided. Stopping..." << endl;
		exit(-1);
	}

	/*Check if there is a slash at the end.*/

  	if (pathname.back() != '/')
  	{
  		pathname.append("/");
  	}

  	svar.outfolder = pathname;

  	/*Check output folder for any prexisting files*/
  	if(svar.outtype == 0)
  	{
  		string file = pathname;
  		file.append("Fuel.szplt.szdat");
#ifdef DEBUG
  		dbout << "Checking for existence of previous szplt files." << endl;
  		dbout << "Path: " << file << endl;
#endif
  		struct stat info;
  		if(stat( file.c_str(), &info ) == 0)
  		{
	  		string cmd = "exec rm -r \"";
	  		cmd.append(pathname);
	  		cmd.append("\"*.szplt.sz*");
#ifdef DEBUG
  		dbout << "Files found. Attempting to remove." << endl;
  		dbout << "Command: " << cmd << endl;
#endif
	  		if(system(cmd.c_str()))
	  		{
		    	cout << "System command failed to execute." << endl;
		    	cout << "Command: " << cmd << endl;
		    	exit(-1);
		    }
		}
  	}
  	else if(svar.outtype == 2)
  	{	/*Create h5 folder*/
  		pathname.append("h5/");
  		check_folder(pathname);
  		string file = pathname;
  		file.append("fuel_0.00e+00.h5part");
  		struct stat info;
  		if(stat( file.c_str(), &info ) == 0)
  		{
	  		string cmd = "exec rm -r \"";
	  		cmd.append(pathname);
	  		cmd.append("\"*.h5part");
	  		if(system(cmd.c_str()))
	  		{
		    	cout << "System command failed to execute." << endl;
		    	cout << "Command: " << cmd << endl;
		    	exit(-1);
		    }
		}
  	}

	return 0;
}

int getInt(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	std::stringstream sstr;
	sstr << line;
	int i;
	if(sstr >> i)
	{
#ifdef DEBUG
		dbout << name << ": " << i << endl;
#endif	
		return i; 
	}
	else
	{
		cout << "Line does not contain a value. Please check your input." << endl;
		cout << "Expecting: " << name << endl;
	}
		
	return 0;
}

real getDouble(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	std::stringstream sstr;
	sstr << line;
	real d;
	if(sstr >> d)
	{
#ifdef DEBUG
		dbout << name << ": " << d << endl;
#endif	
		return d; 
	}
	else
	{
		cout << "Line does not contain a value. Please check your input." << endl;
		cout << "Expecting: " << name << endl;
	}
		
	return 0;
}

std::string getString(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	size_t ptr = line.find_first_of(' ');
	string result = line.substr(0,ptr);
#ifdef DEBUG
	dbout << name << ": " << result << endl;
#endif	
	return result; 
}

StateVecI getIVector(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);

	lineno++;
	std::istringstream sline(line);
	// cout << sline.str() << endl;
	StateVecI x;
	sline >> x(0); sline >> x(1); 

	#if SIMDIM == 2
		int temp;
		if(sline >> temp)
		{
			cout << "\tWARNING: 3D Input provided for a 2D Simulation." << endl;
			cout << "\t         The third dimension shall be ignored." << endl;
			cout << "\tLine " << lineno << ": " << endl;
			cout << "\t" << sline.str() << endl;
		}
#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1) << endl;
#endif	
#endif
#if SIMDIM == 3
	sline >> x(2);
	if (!sline)
	{	
		cout << "2D input provided. Please provide a 3D file." << endl;
		cout << "Incorrect line " << lineno << ": " << endl;
		cout << sline.str() << endl;
		exit(-1);
	}
#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1) << "  " << x(2) << endl;
#endif	
#endif

	return x;
}

StateVecD getDVector(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	std::istringstream sline(line);
	
	StateVecD x;
	sline >> x(0); sline >> x(1);

#if SIMDIM == 2
	real temp;
	if(sline >> temp)
	{
		cout << "\tWARNING: 3D Input provided for a 2D Simulation." << endl;
		cout << "\t         The third dimension shall be ignored." << endl;
		cout << "\tLine " << lineno << ": " << endl;
		cout << "\t" << sline.str() << endl;
	}
#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1) << endl;
#endif	
#endif
#if (SIMDIM == 3)
	sline >> x(2);
	if (!sline)
	{	
		cout << "2D input provided. Please provide a 3D file." << endl;
		cout << "Incorrect line " << lineno << ": " << endl;
		cout << sline.str() << endl;
		exit(-1);
	}
#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1) << "  " << x(2) << endl;
#endif	
#endif	

	return x;
}

/*Function for a 2D Vector (e.g. Newmark Beta parameters)*/
Eigen::Vector2d getvector(ifstream& In, uint& lineno, const string& name)
{
	string line;
	getline(In,line);
	lineno++;
	std::istringstream sline(line);
	
	Eigen::Vector2d x;
	sline >> x[0]; sline >> x[1]; 
#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1) << endl;
#endif	
	return x;
}

void GetAero(AERO& avar, const FLUID& fvar, const real rad)
{
	#if SIMDIM == 3
		avar.L = rad * std::cbrt(3.0/(4.0*M_PI));
	#endif
	#if SIMDIM == 2
		avar.L = 2*rad;
	#endif
	avar.td = (2.0*fvar.rho0*pow(avar.L,SIMDIM-1))/(avar.Cd*fvar.mu);

	avar.omega = sqrt((avar.Ck*fvar.sig)/(fvar.rho0*pow(avar.L,SIMDIM))-1.0/pow(avar.td,2.0));

	avar.tmax = -2.0 *(atan(sqrt(pow(avar.td*avar.omega,2.0)+1)
					+avar.td*avar.omega) - M_PI)/avar.omega;

	avar.Cdef = 1.0 - exp(-avar.tmax/avar.td)*(cos(avar.omega*avar.tmax)+
		1/(avar.omega*avar.td)*sin(avar.omega*avar.tmax));
	avar.ycoef = 0.5*avar.Cdef*(avar.Cf/(avar.Ck*avar.Cb))*(fvar.rhog*avar.L)/fvar.sig;

	//cout << avar.ycoef << "  " << avar.Cdef << "  " << avar.tmax << "  " << endl;
}

RotMat GetRotationMat(StateVecD& angles)
{
	if (SIMDIM == 3)
	{
		RotMat rotx, roty, rotz;
		rotx << 1.0, 0.0            , 0.0           ,
			    0.0, cos(angles(0)) , sin(angles(0)),
				0.0, -sin(angles(0)), cos(angles(0));

		roty << cos(angles(1)) , 0.0 , -sin(angles(1)),
			    0.0            , 1.0 , 0.0            ,
				sin(angles(1)) , 0.0 , cos(angles(1));

		rotz << cos(angles(2)) , sin(angles(2)) , 0.0 ,
			    -sin(angles(2)), cos(angles(2)) , 0.0 ,
				0.0            , 0.0            , 1.0 ;

		return rotx*roty*rotz;
	}
	else if (SIMDIM == 2)
	{
		RotMat rot;
		rot << cos(angles(0)), -sin(angles(0)),
		       sin(angles(0)),  cos(angles(0));

		return rot;
	}

	return RotMat::Zero();

}

void GetInput(int argc, char **argv, SIM& svar, FLUID& fvar, AERO& avar)
{
	uint justPost = 0;
	real Hfac = 0;
	StateVecD angle;
	if (argc > 3) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only two input arguments accepted,\n";
		cout << "1: Input file directoy   2: Output file directory.\n";
		cout << "Other inputs will be ignored." << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tERROR: No inputs provided. Stopping... \n";
    	exit(-1);    	
    }
    else if (argc > 1)
    {	/*Get parameters if it has been provided*/
    	// cout << argv[1] << endl;
    	string file = argv[1];
    	
    	
    	if(file == "-post")
    	{
    		cout << "Post processing option selected." << endl;
    		justPost = 1;
    		file = argv[2];
    	}


    	if(file.back() != '/')
	    	file.append("/");

	    svar.infolder = file;
		file.append("Settings.dat");
#ifdef DEBUG
		dbout << "Reading settings file. Path:" << endl << file << endl;
#endif
    	std::ifstream in(file);
	  	if(in.is_open()) 
	  	{	/*Simulation parameters*/
	  		cout << "Input file opened. Reading settings..." << endl;
		  	uint lineno = 0;
	  		svar.framet = getDouble(in, lineno, "Frame time");
	  		svar.Nframe = getInt(in, lineno, "Number of frames");
	  		svar.outframe = getInt(in, lineno, "Output frame info");
	  		svar.outtype = getInt(in, lineno, "Output data type");
	  		svar.outform = getInt(in, lineno, "Output contents");
	  		svar.boutform = getInt(in, lineno, "Boundary time output");
	  		svar.gout = getInt(in, lineno, "Output ghost particles to file");
	  		svar.subits = getInt(in, lineno, "Max sub iterations");
	  		svar.nmax = getInt(in, lineno, "Max particle add rounds");	
	  		svar.xyPART = getIVector(in, lineno, "Particles in each coordinate");
	  		svar.Pstep = getDouble(in, lineno, "Particle initial spacing");
	  		svar.Bstep = getDouble(in, lineno, "Boundary spacing factor");
	  		svar.Bcase = getInt(in, lineno, "Simulation boundary case");
	  		avar.acase = getInt(in, lineno, "Simulation Aerodynamic case");
	  		svar.ghost = getInt(in, lineno, "Ghost particles?");
	  		svar.Start = getDVector(in, lineno, "Starting position");
	  		if(svar.Bcase < 2)
	  		{
	  			svar.Box= getDVector(in, lineno, "Box dimensions");
	  			(void)getDouble(in, lineno, "Skip");
	  			fvar.pPress = getDouble(in, lineno, "Pipe pressure");
	  			if(svar.Box(0) < 0 || svar.Box(1) < 0 
	  			#if SIMDIM == 3
	  				|| svar.Box(2) < 0
  				#endif
	  			)
	  			{
	  				cout << "Box dimensions are negative. Please check the input and try again." << endl;
	  				cout << "Line " << lineno << endl;
	  				exit(-1);
	  			}
	  		}
	  		else if(svar.Bcase > 1 && svar.Bcase < 8)
	  		{	
	  			StateVecD angles = getDVector(in, lineno, "Starting angle");
	  			angle = angles;
	  			angles = angles *M_PI/180;
	  			svar.Rotate = GetRotationMat(angles);
	  			svar.Transp = svar.Rotate.transpose();
		  		svar.Jet = getvector(in, lineno, "Jet dimensions");
		  		fvar.pPress = getDouble(in, lineno, "Pipe pressure");
		  		avar.vJet = StateVecD::Zero(); avar.vInf = StateVecD::Zero();
		  		avar.vJet(1) = getDouble(in, lineno, "Jet velocity");  
		  		avar.vJet = svar.Rotate*avar.vJet;
		  		avar.vInf = getDVector(in, lineno, "Freestream velocity");
		  		if(avar.acase == 2 || avar.acase == 3)
		  		{
		  			avar.a = getDouble(in, lineno, "a");
		  			avar.h1 = getDouble(in, lineno, "h1");
		  			avar.b = getDouble(in, lineno, "b");
		  			avar.h2 = getDouble(in, lineno, "h2");
		  		}
		  		if(avar.acase > 5)
		  		{
		  			cout << "Aerodynamic case is not in design. Stopping..." << endl;
		  			exit(-1);
		  		}
	  		}
	  		else
	  		{
	  			cout << "Boundary case not within design. Stopping." << endl;
	  			exit(-1);
	  		}
	  		/*Get post processing options*/
	  		svar.afterSim = getInt(in, lineno, "Post processing");
	  		svar.cellSize = getDouble(in, lineno, "Post processing mesh size");
	  		svar.postRadius = getDouble(in, lineno, "Post processing support radius");

#ifdef DEBUG
			dbout << "Closing settings file"  << endl;
#endif
			in.close();
	  	}
	  	else {
		    cerr << "Error opening the settings file." << endl;
		    exit(-1);
	  	}
	}

	/*Get fluid properties from fluid.dat*/
	svar.scale = 1.0;
	string file = svar.infolder;
	file.append("Fluid.dat");
#ifdef DEBUG
	dbout << "Reading fluid file. Path:" << endl << file << endl;
#endif
	std::ifstream fluid(file);
	if (fluid.is_open())
	{	/*Fluid parameters read*/
		uint lineno = 0;
		Eigen::Vector2d nb = getvector(fluid, lineno, "Newmark-Beta terms");
		svar.beta = nb[0];	svar.gamma = nb[1];
		Hfac = getDouble(fluid, lineno, "Smoothing length factor"); /*End of state read*/

	  	fvar.alpha = getDouble(fluid, lineno, "Artificial visc");
  		fvar.contangb = getDouble(fluid, lineno, "contact angle");
  		fvar.rho0 = getDouble(fluid, lineno, "Fluid density rho0");
  		fvar.rhog = getDouble(fluid, lineno, "Air density rhog");
  		fvar.Cs = getDouble(fluid, lineno, "Speed of sound");
  		fvar.mu = getDouble(fluid, lineno, "Fluid viscosity");
  		fvar.mug = getDouble(fluid, lineno, "Air viscosity");
  		fvar.sig = getDouble(fluid, lineno, "Surface Tension");
  		if(svar.Bcase == 6 || svar.Bcase == 4)
  		{
	  		fvar.gasVel = getDouble(fluid, lineno, "Gas ref Vel");
	  		fvar.gasPress = getDouble(fluid, lineno, "Get ref Press");
	  		fvar.T = getDouble(fluid, lineno, "Gas ref Temp");

	  		if(svar.Bcase == 6)
	  		{
		  		svar.meshfile = getString(fluid, lineno, "Mesh input file");
		  		svar.solfile = getString(fluid,lineno, "Mesh solution file");
		  		svar.scale = getDouble(fluid,lineno, "Mesh scale");
	  		}
  		}
#ifdef DEBUG
		dbout << "Closing fluid file"  << endl;
#endif
  		fluid.close();

  		if(svar.Bcase == 4)
  		{
  			fvar.gasVel=avar.vInf(1);
  		}
	}
	else 
	{
		cerr << "Fluid.dat not found. Assuming standard set of parameters (water)." << endl;
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
	  	

	  	
		fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
		fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/

		/*Pipe Pressure calc*/
		real rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);

		/*Defining pstep then finding dx*/
		// fvar.Simmass = fvar.rho0*pow(svar.Pstep,SIMDIM);
		// svar.dx = pow(fvar.Simmass/rho, 1.0/real(SIMDIM));
		

		// Defining dx to fit the pipe, then find rest spacing
		uint ndiam = ceil(abs(svar.Jet(0)/svar.Pstep));
		svar.dx = 0.999*svar.Jet(0)/real(ndiam);	

	 	svar.Pstep = svar.dx * pow(rho/fvar.rho0,1.0/3.0);

	 	/*Mass from spacing and density*/
		fvar.Simmass = fvar.rho0*pow(svar.Pstep,SIMDIM); 
		fvar.Boundmass = fvar.Simmass;
		fvar.gasDynamic = 0.5*fvar.rhog*fvar.gasVel;

		svar.addcount = 0;
	  	svar.dt = 2E-010; 			/*Initial timestep*/
	  	svar.t = 0.0;				/*Total simulation time*/
	  	fvar.H = Hfac*svar.Pstep;
	  	fvar.HSQ = fvar.H*fvar.H; 

		fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/
		svar.Bclosed = 0; 		/*Boundary begins open*/
		svar.psnPts = 0; 		/*Start with no pitson points*/
	  	svar.delNum = 0;

#if SIMDIM == 2
				fvar.correc = (7/(4*M_PI*fvar.H*fvar.H));
				svar.simPts = svar.xyPART[0]*svar.xyPART[1];
#endif
#if SIMDIM == 3
				fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
				svar.simPts = svar.xyPART[0]*svar.xyPART[1]*svar.xyPART[2]; /*total sim particles*/
#endif

#ifdef DEBUG
		dbout << "Tait Gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
		dbout << "Pipe rho: " << rho << endl;
		dbout << "Number of fluid particles along diameter: " << ndiam << endl;
		dbout << "Pipe step (dx): " << svar.dx << endl;
		dbout << "Particle mass: " << fvar.Simmass << endl;
		dbout << "Freestream initial spacing: " << svar.Pstep << endl;
		dbout << "Support Radius: " << fvar.H << endl;
		dbout << "Gas Dynamic pressure: " << fvar.gasDynamic << endl << endl;
#endif

		cout << "****** SIMULATION SETTINGS *******" << endl;
		cout << "Tait gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
		cout << "Newmark-Beta parameters: " << svar.beta << ", " << svar.gamma << endl;
		cout << "Boundary case: " << svar.Bcase << endl;
		cout << "Aerodynamic case:" << avar.acase << endl;
		cout << "Speed of sound: " << fvar.Cs << endl;
		cout << endl;
		cout << "****** PIPE SETTINGS *******" << endl;
		cout << "Pipe density: " << rho << endl;
		cout << "Pipe step (dx): " << svar.dx << endl;
		cout << "Pipe diameter: " << svar.Jet(0) << endl;
		cout << "Number of fluid particles along diameter: " << ndiam << endl;
		cout << "Pipe start position: " << svar.Start(0) << "  " << svar.Start(1);
		#if SIMDIM == 3
		cout << "  " << svar.Start(2);
		#endif
		cout << endl;
		cout << "Pipe start rotation: " << angle(0) << "  " << angle(1);
		#if SIMDIM == 3
		cout << "  " << angle(2);
		#endif
		cout << endl;
		cout << "Jet velocity: " << avar.vJet << endl;

		cout << endl;		
		cout << "****** RESTING FUEL SETTINGS *******" << endl;
		cout << "Particle spacing: " << svar.Pstep << endl;
		cout << "Support radius: " << fvar.H << endl;
		cout << "Particle mass: " << fvar.Simmass << endl;
		cout << "Liquid viscosity: " << fvar.mu << endl;
		cout << "Resting density: " << fvar.rho0 << endl;
		
		cout << endl;
		cout << "****** FREESTREAM SETTINGS ******" << endl;
		if(svar.Bcase != 6)
			cout << "Gas Velocity: " << avar.vInf << endl;

		if(svar.Bcase == 6)
		{
			cout << "Reference velocity: " << fvar.gasVel << endl;
			cout << "Reference pressure: " << fvar.gasPress << endl;
			cout << "Reference temperature: " << fvar.T << endl;
		}
		cout << "Gas density: " << fvar.rhog << endl;
		cout << "Gas viscosity: " << fvar.mug << endl;
		cout << "Gas dynamic pressure: " << fvar.gasDynamic << endl << endl;

		#if SIMDIM == 3
			if(svar.Bcase == 4)
			{
				svar.vortex.Init(svar.infolder);	
				svar.vortex.GetGamma(avar.vInf);	
			}
		#endif

		GetAero(avar, fvar, fvar.H);
	
	if (justPost==1)
	{
		cout << "Starting post processing from output of the same case." << endl;
	
		svar.afterSim = 0;
		/*Open an output directory in the name of the input file, under Outputs*/
		char cCurrentPath[FILENAME_MAX];
		if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
			exit(-1);

		/*Check the folder for the grid szplt files*/

		string pathname = cCurrentPath;
	  	string input = argv[2];
	  	uint pos1 = input.find("/"); 	uint pos2 = input.find(".");
	  	string file = input.substr(pos1+1,pos2-pos1-1);
	  	pathname.append("/Outputs/");
	  	pathname.append(file);
	  	svar.outfolder = pathname;

	  	fvar.H = svar.postRadius;
	  	fvar.HSQ = fvar.H*fvar.H; 
		fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/

	  	// cout << "Found path: " << pathname << endl;
		TECMESH postgrid;

		postgrid.DoPostProcessing(svar,fvar);

		cout << "Post Processing complete!" << endl;
		exit(0);
	}
}

/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/
void Write_ASCII_Timestep(std::ofstream& fp, SIM &svar, const State &pnp1, 
	const uint bwrite, const uint start, const uint end, const string& name)
{
	if(bwrite == 1)
	 	fp <<  "ZONE T=\"" << name << "\"";
	else
		fp <<  "ZONE T=\"" << name << "\"";

    fp <<", I=" << end - start << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
    fp << std::left << std::scientific << std::setprecision(6);
    const static uint width = 15;
    
    if(svar.outform == 0)
    {
		for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
		{
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale; 
			fp << "\n";  
	  	}

	  	// if (airP.size() > 0 )
	  	// {
		  // 	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
		  //   ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
		  // 	for(auto p:airP)
		  // 	{
		  // 		for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.xi(i) << " "; 
				// fp << "\n"; 
		  // 	}
	  	// }
	}	
    if(svar.outform == 1)
    {
		for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
		{	
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;

			fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;
	        fp << setw(width) << p->v.norm();
	        fp << setw(width) << p->f.norm() << "\n";
	  	}

	  // 	if (airP.size() > 0 )
	  // 	{
		 //  	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
		 //    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";

		 //  	for(auto p:airP)
		 //  	{
		 //  		for(uint i = 0; i < SIMDIM; ++i)
		 //        	fp << p.xi(i) << " "; 

		 //        fp << p.v.norm() << " ";
		 //        fp << p.f.norm() << " ";
		 //        fp << p.rho << " "  << p.p  << "\n";
		 //  	}
	 	// }
	  	fp << std::flush; 	
    }
    else if (svar.outform == 2)
	{
		for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
		{
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;
	        
	        fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->v(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
				fp << setw(width) << p->f(i); 

			fp << "\n";
	  	}  

	  	// if (airP.size() > 0 )
	  	// {
		  // 	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
		  //   ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
		    
		  // 	for(auto p:airP)
		  // 	{
		  // 		for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.xi(i) << " "; 

		  //     	fp << p.f.norm() << " " << p.Af.norm() << " " << p.Sf.norm() << " ";
		  //       for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.cellV(i) << " "; 

		  //       fp << p.b << " " << p.theta  << "\n"; 
		  // 	}
	  	// }	
	  	fp << std::flush;
    }
    else if (svar.outform == 3)
    {
    	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
		{
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;
	        
	        fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->v(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->f(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
				fp << setw(width) << p->cellV(i);

			fp << setw(width) << p->cellRho << setw(width) << p->cellP;
	        fp << setw(width) << p->cellID << "\n"; 
	  	}
    }
    else if (svar.outform == 4)
    {
    	for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
		{	
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;

			fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;
	        fp << setw(width) << p->v.norm();
	        fp << setw(width) << p->f.norm();
	        fp << setw(width) << p->b << setw(width) << p->theta;
	        fp << setw(width) << p->Af.norm() << "\n";
	  	}

    }
}

void Write_ASCII_header(std::ofstream& fp, SIM &svar)
{
	string variables;

#if SIMDIM == 2
	variables = "\"X\" \"Z\"";  
	if (svar.outform == 1)
	{
		variables = "\"X\" \"Z\" \"rho\" \"P\" \"m\" \"v\" \"a\"";
	}
	else if (svar.outform == 2)
	{
		variables = "\"X\" \"Z\" \"rho\" \"P\" \"m\" \"v_x\" \"v_z\" \"a_x\" \"a_z\"";
	}
	else if (svar.outform == 3)
	{
		variables = 
"\"X\" \"Z\" \"rho\" \"P\" \"m\" \"v_x\" \"v_z\" \"a_x\" \"a_z\" \"Cell_Vx\" \"Cell_Vz\" \"Cell_Rho\" \"Cell_P\" \"Cell_ID\"";
	}
	else if (svar.outform == 4)
	{
		variables = "\"X\" \"Z\" \"rho\" \"P\" \"m\" \"v\" \"a\" \"b\" \"Neighbours\" \"Aero\"";
	}

#endif

#if SIMDIM == 3
	variables = "\"X\" \"Y\" \"Z\"";  
	if (svar.outform == 1)
	{
		variables = "\"X\" \"Y\" \"Z\" \"rho\" \"P\" \"m\" \"v\" \"a\"";
	}
	else if (svar.outform == 2)
	{
		variables = "\"X\" \"Y\" \"Z\" \"rho\" \"P\" \"m\" \"v_x\" \"v_y\" \"v_z\" \"a_x\" \"a_y\" \"a_z\"";
	}
	else if (svar.outform == 3)
	{
		variables = 
"\"X\" \"Y\" \"Z\" \"rho\" \"P\" \"m\" \"v_x\" \"v_y\" \"v_z\" \"a_x\" \"a_y\" \"a_z\" \"Cell_Vx\" \"Cell_Vy\" \"Cell_Vz\" \"Cell_Rho\" \"Cell_P\" \"Cell_ID\"";
	}
	else if (svar.outform == 4)
	{
		variables = "\"X\" \"Y\" \"Z\" \"rho\" \"P\" \"m\" \"v\" \"a\" \"b\" \"Neighbours\" \"Aero\"";
	}

#endif

	fp << variables << "\n";
}

void Write_Boundary_ASCII(std::ofstream& fp, SIM &svar, State &pnp1)
{	
	
		Write_ASCII_header(fp,svar);
		fp <<  "ZONE T=\"Boundary Data\"" << ", I=" << svar.bndPts << ", F=POINT" << "\n";
		if (svar.outform == 0)
		{
		  	for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
			{
		        for(uint i = 0; i < SIMDIM; ++i)
		        	fp << b->xi(i) << " "; 
		        fp << "\n";
		  	}
		}
		if (svar.outform == 1)
		{	/*Fluid Data*/
		    for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
			{	
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << b->xi(i) << " "; 
		        
		        fp << b->v.norm() << " ";
		        fp << b->f.norm() << " ";
		        fp << b->rho << " "  << b->p << "\n";
		  	}
		}
		if (svar.outform == 2)
		{	/*Research Data*/    
	        for (auto b=pnp1.begin(); b!=std::next(pnp1.begin(),svar.bndPts); ++b)
			{
		        for(uint i = 0; i < SIMDIM; ++i)
    				fp << b->xi(i) << " "; 
		        
		        fp << b->f.norm() << " " << b->Af.norm() << " " << b->Sf.norm() << " "; 
		        for(uint i = 0; i < SIMDIM; ++i)
		        	fp << b->Af(i) << " ";

		        fp << b->b << " " << b->theta << "\n"; 
		  	}		  	
		}
}

#endif