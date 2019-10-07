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
#include "Eigen/Core"
#include "Eigen/StdVector"
#include "Eigen/LU"
#include "Var.h"
#include "BinaryIO.h"

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
	struct stat info;
	int sreturn=0;
	if( stat( pathname.c_str(), &info ) != 0 )
  	{	/*Output directory doesn't exist, so create it*/
		pathname.insert(0,"\"");
		pathname.append("\"");
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

  	svar.outfolder = pathname;
	return 0;
}

int getInt(ifstream& In, uint& lineno)
{
	string line;
	getline(In,line);
	lineno++;
	int i = stoi(line);
	return i;
}

double getDouble(ifstream& In, uint& lineno)
{
	string line;
	getline(In,line);
	lineno++;
	double d = stod(line);
	return d; 
}

std::string getString(ifstream& In, uint& lineno)
{
	string line;
	getline(In,line,' ');
	lineno++;
	return line; 
}

StateVecI getIVector(ifstream& In, uint& lineno)
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
	#endif
	return x;
}

StateVecD getDVector(ifstream& In, uint& lineno)
{
	string line;
	getline(In,line);
	lineno++;
	std::istringstream sline(line);
	
	StateVecD x;
	sline >> x(0); sline >> x(1);

	#if SIMDIM == 2
		double temp;
		if(sline >> temp)
		{
			cout << "\tWARNING: 3D Input provided for a 2D Simulation." << endl;
			cout << "\t         The third dimension shall be ignored." << endl;
			cout << "\tLine " << lineno << ": " << endl;
			cout << "\t" << sline.str() << endl;
		}
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
	#endif	

	return x;
}

/*Function for a 2D Vector (e.g. Newmark Beta parameters)*/
Eigen::Vector2d getvector(ifstream& In, uint& lineno)
{
	string line;
	getline(In,line);
	lineno++;
	std::istringstream sline(line);
	
	Eigen::Vector2d x;
	sline >> x[0]; sline >> x[1]; 
		
	return x;
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

void GetInput(int argc, char **argv, SIM &svar, FLUID &fvar, CROSS &cvar)
{
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
    	if(file.back() != '/')
	    	file.append("/");

	    svar.infolder = file;
		file.append("Settings.dat");
    	std::ifstream in(file);
	  	if(in.is_open()) 
	  	{	/*Simulation parameters*/
	  		cout << "Input file opened. Reading settings..." << endl;
		  	uint lineno = 0;
	  		svar.framet = getDouble(in, lineno);
	  		svar.Nframe = getInt(in, lineno);
	  		svar.outframe = getInt(in, lineno);
	  		svar.outtype = getInt(in, lineno);
	  		svar.outform = getInt(in, lineno);
	  		svar.frameout = getInt(in, lineno);
	  		svar.subits = getInt(in, lineno);
	  		svar.nmax = getInt(in, lineno);	
	  		svar.xyPART = getIVector(in, lineno);
	  		svar.Pstep = getDouble(in, lineno);
	  		svar.Bstep = getDouble(in, lineno);
	  		svar.Bcase = getInt(in, lineno);
	  		cvar.acase = getInt(in, lineno);
	  		svar.ghost = getInt(in, lineno);
	  		svar.Start = getDVector(in, lineno);
	  		if(svar.Bcase < 2)
	  		{
	  			svar.Box= getDVector(in, lineno);
	  			(void)getDouble(in, lineno);
	  			fvar.pPress = getDouble(in, lineno);
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
	  		else if(svar.Bcase > 1 && svar.Bcase < 7)
	  		{	
	  			StateVecD angles = getDVector(in, lineno);
	  			angles = angles *M_PI/180;
	  			svar.Rotate = GetRotationMat(angles);
	  			svar.Transp = svar.Rotate.transpose();
		  		svar.Jet = getvector(in, lineno); /*Defined in VLM.h. Reused here*/
		  		fvar.pPress = getDouble(in, lineno);
		  		cvar.vJet = StateVecD::Zero(); cvar.vInf = StateVecD::Zero();
		  		cvar.vJet(1) = getDouble(in, lineno);  
		  		cvar.vJet = svar.Rotate*cvar.vJet;
		  		cvar.vInf = getDVector(in, lineno);
		  		if(cvar.acase >= 2)
		  		{
		  			cvar.a = getDouble(in, lineno);
		  			cvar.h1 = getDouble(in, lineno);
		  			cvar.b = getDouble(in, lineno);
		  			cvar.h2 = getDouble(in, lineno);
		  		}
		  		
		  		if(cvar.acase > 5)
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
			in.close();
	  	}
	  	else {
		    cerr << "Error opening the settings file." << endl;
		    exit(-1);
	  	}
	}

	/*Get fluid properties from fluid.dat*/
	string file = svar.infolder;
	file.append("Fluid.dat");
	std::ifstream fluid(file);
	if (fluid.is_open())
	{	/*Fluid parameters read*/
		uint lineno = 0;
		Eigen::Vector2d nb = getvector(fluid, lineno);
		svar.beta = nb[0];	svar.gamma = nb[1];
		double Hfac = getDouble(fluid, lineno); /*End of state read*/
	  	fvar.H= Hfac*svar.Pstep;
	  	fvar.alpha = getDouble(fluid, lineno);
  		fvar.contangb = getDouble(fluid, lineno);
  		fvar.rho0 = getDouble(fluid, lineno);
  		fvar.rhog = getDouble(fluid, lineno);
  		fvar.Cs = getDouble(fluid, lineno);
  		fvar.mu = getDouble(fluid, lineno);
  		fvar.mug = getDouble(fluid, lineno);
  		fvar.sig = getDouble(fluid, lineno);
  		fvar.gasVel = getDouble(fluid, lineno);
  		fvar.gasPress = getDouble(fluid, lineno);
  		fvar.T = getDouble(fluid, lineno);
  		svar.meshfile = getString(fluid, lineno);
  		fluid.close();
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
  	svar.addcount = 0;
  	svar.dt = 2E-010; 		/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	fvar.HSQ = fvar.H*fvar.H; 

	fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
	switch(SIMDIM)
	{
		case 2:
			fvar.correc = (7/(4*M_PI*fvar.H*fvar.H));
			svar.simPts = svar.xyPART[0]*svar.xyPART[1];
			break;
		case 3:
			fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
			svar.simPts = svar.xyPART[0]*svar.xyPART[1]*svar.xyPART[2]; /*total sim particles*/
			break;
		default:
			cout << "Simulation Dimension mismatch. Stopping" << endl;
	  		exit(-1);
	  		break;
	}
  	
  	/*Mass from spacing and density*/
	fvar.Simmass = fvar.rho0*pow(svar.Pstep,SIMDIM); 
	fvar.Boundmass = fvar.Simmass;
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/
	/*Pipe Pressure calc*/
	ldouble rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
	svar.dx = pow(fvar.Simmass/rho, 1.0/double(SIMDIM));
	// cout << rho << "  " << svar.dx << endl;
	fvar.gasDynamic = 0.5*fvar.rhog*fvar.gasVel;


	#if SIMDIM == 3
		if(svar.Bcase == 4)
		{
			svar.vortex.Init(svar.infolder);	
			svar.vortex.GetGamma(cvar.vInf);	
		}
	#endif

	GetAero(fvar, fvar.H);
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(uint ii=0; ii < num - 1; ++ii){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

void Skip_Variable(ifstream& fin, const int np)
{
	for(int ii = 0; ii < ceil(float(np)/5.0); ++ii)
	{
		fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	}
}

void Get_Vector(ifstream& fin, const uint np, std::vector<StateVecD>& var, const uint yskip)
{
	string line;
	for(uint dim =0; dim < SIMDIM; dim++)
	{	
		uint kk = 0;
		/*If in 2D, skip the y-dimension, and read the z-dim.*/
		if(yskip == 1)
		{
			if (dim == 1)
			{
				for(uint ii = 0; ii < ceil(float(np)/5.0); ++ii)
				{
					fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
					// getline(fin,line);
				}
			}
		}
		for(int ii = 0; ii < ceil(float(np)/5.0); ++ii)
		{
			getline(fin, line);
			std::istringstream sline(line);

			for (uint jj = 0; jj < 5; ++jj)
			{	
				double temp;
				sline >> temp;

				if (!sline)
					break;
				else
				{
					var[kk](dim) = temp;
					++kk;
				}
			}
		}
		if (kk!= var.size())
		{
			cout << "Mismatch of array size.\n" << 
			" Not all of the array has been populated." << endl;
			cout << "populated: " << kk << " Array size: " << var.size() << endl;
			cout << var[kk] << endl;
		} 
	}
}

template <class T>
void Get_Scalar_Data(ifstream& fin, const uint np, T& var)
{
	string line;

	uint kk = 0;
	for(uint ii = 0; ii < ceil(float(np)/5.0); ++ii)
	{
		getline(fin, line);
		std::istringstream sline(line);

		for (uint jj = 0; jj < 5; ++jj)
		{	
			double temp;
			sline >> temp;

			if (!sline)
				break;
			else
			{
				var[kk] = temp;
				++kk;
			}
		}
	}

	if (kk!= var.size())
	{
		cout << "Mismatch of array size.\n" << 
		" Not all of the array has been populated." << endl;
		cout << "populated: " << kk << " Array size: " << var.size() << endl;
		cout << var[kk] << endl;
	}
}

void Get_Cells(ifstream& fin, const uint nE, const uint nCverts, std::vector<std::vector<uint>>& cell)
{
	string line;
	// getline(fin, line);
	// cout << line << endl;
	for(uint ii = 0; ii < nE; ++ii)
	{
		getline(fin, line);
		std::istringstream sline(line);

		for (uint jj = 0; jj < nCverts; ++jj)
		{	
			uint temp;
			sline >> temp;
			cell[ii][jj] = (temp-1);								
		}
	}
}

std::vector<ldouble> CpToPressure(const std::vector<ldouble>& Cp, const FLUID& fvar)
{
	std::vector<ldouble> press(Cp.size());
	#pragma omp parallel for shared(Cp)
	for (uint ii = 0; ii < Cp.size(); ++ii)
	{
		press[ii] = Cp[ii]*fvar.gasDynamic + fvar.gasPress;
	}
	return press;
}

void NormalisePressure(MESH &cells, const FLUID& fvar)
{
	/*Check for density and pressure information*/
	/*Create the data based on the other.*/
	if(cells.pointP[0]!=0)
	{
		cells.pointP = CpToPressure(cells.pointCp,fvar);
	}

	if(cells.pointRho[0]!=0)
	{
		for(uint ii = 0; ii < cells.pointRho.size(); ++ii)
		{
			cells.pointRho[ii] = cells.pointP[ii]/(fvar.Rgas*fvar.T);
		}
	}

	/*Normalise pressure to be in terms of the Tait equation*/
	/*I have no idea if this is conservative...*/
	for(uint ii = 0; ii < cells.pointP.size(); ++ii)
	{
		cells.pointP[ii] -= fvar.gasPress;

		// cells.pointRho[ii] = fvar.rhog*pow((cells.pointP[ii]/fvar.B +1),1/fvar.gam);
	}
}


template <class T>
void Average_Point_to_Cell(std::vector<T>& pData, std::vector<T>& cData,
							const std::vector<std::vector<uint>>& elems, const T zero)
{
	uint nVerts = elems[0].size();
	for(uint ii = 0; ii < elems.size(); ++ii)
	{
		T sum = zero;
		for (auto jj:elems[ii])
		{
			sum += pData[jj];
		}
		cData[ii] = sum/nVerts;
	}
}

void Read_TAUMESH(string input, MESH& cells, FLUID& fvar)
{
	// input.append(svar.meshfile);
	std::ifstream fin(input, std::ios::in);

	if(!fin.is_open())
	{
		cout << "Couldn't open mesh file. Stopping." << endl;
		cout << "Path attempted: " << input << endl;
		exit(-1);
	}
	else 
	{
		cout << "Mesh file open, reading data..." << endl;
	}
	std::string line;
	getline(fin,line);
	getline(fin,line);

	uint veltype = 0;
	if(line.find("\"x_velocity\"")!=string::npos)
	{
		if(line.find("\"y_velocity\"")!=string::npos)
		{
			if(line.find("\"z_velocity\"")!=string::npos)
			{
				cout << "All velocity components found!" << endl;
			}
			else
			{
				cout << "velocity components \"u\" and \"v\" found, but no \"w\"..." << endl;
				cout << "I can't work with this. Stopping..." << endl;
				exit(-1);
			}
		}
		else if (line.find("\"z_velocity\"")!=string::npos)
		{
			cout << "velocity components \"u\" and \"w\" found, but no \"v\"..." << endl;
			#if SIMDIM == 2
				cout << "I can work with this. Continuing..." << endl;
				veltype = 1;
			#else 
				cerr << "SIMDIM is 3, and \"v\" velocity component missing. Stopping." << endl;
				exit(-1);
			#endif
		}
		else
		{
			cout << "Only velocity component \"u\" found, but no \"v\" and \"w\"..." << endl;
			cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
		}
	}
	else if (line.find("\"y_velocity\"")!=string::npos)
	{
		if(line.find("\"z_velocity\"")!=string::npos)
		{
			cout << "velocity components \"v\" and \"w\" found, but no \"u\"..." << endl;
			cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
		}
		else
		{
			cout << "Only velocity component \"v\" found, but no \"u\" and \"w\"..." << endl;
			cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
		}
	}
	else if (line.find("\"z_velocity\"")!=string::npos)
	{
		cout << "Only velocity component \"w\" found, but no \"u\" and \"v\"..." << endl;
		cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
	}
	else
	{
		cout << "Warning: No velocity components provided.\n" ;
		cout << "I can't work with this. Stopping..." << endl;
		exit(-1);
	}

	uint nvar = std::count(line.begin(),line.end(),'\"');
	nvar /= 2;
	std::size_t ptr = line.find("\"x_velocity\"");
	uint velstart = std::count(line.begin(),line.begin()+ptr, '\"');
	velstart/=2;

	/*Check to see if there is pressure data available*/
	ptr = line.find("\"pressure\"");
	uint pressOrcp = 1;
	uint cpstart = 0;
	if (ptr != string::npos)
	{
		cout << "Pressure data directly available!" << endl;
		pressOrcp = 0;
		cpstart = std::count(line.begin(),line.begin()+ptr, '\"');
		cpstart/=2;
	}
	else
	{
		ptr = line.find("\"cp\"");	
		if (ptr != string::npos)
		{
			cpstart = std::count(line.begin(),line.begin()+ptr, '\"');
			cpstart/=2; 
		}
		else 
		{
			cout << "Couldn't find any pressure data" << endl;
		}
	}

	/*Check to see if there is density data available*/
	ptr = line.find("\"density\"");
	uint densstart = 0;
	if (ptr != string::npos)
	{
		cout << "Density data directly available!" << endl;
		
		densstart = std::count(line.begin(),line.begin()+ptr, '\"');
		densstart/=2;
	}
	else
	{

	}



	/*Next bit depends on dimensions. If 3D, read the hexa data.*/
	/*If 2D, skip this and read the symmetry plane data.*/

	getline(fin, line); /*Get Zone line data. Tells which type of volume*/
	uint nF, nCverts, nFverts;
	
	if(line.find("hexa")!=string::npos)
	{
		/*N verts = 8, N faces = 6, N edges = 12*/
		#if SIMDIM == 2
			nCverts = 4;
			nF = 0;
			nFverts = 0;
		#else
			nF = 6;
			nCverts = 8;
			nFverts = 4;
		#endif
	}
	else if(line.find("tetra")!=string::npos)
	{
		/*N verts = 4, N faces = 4, N edges = 6*/
		#if SIMDIM == 2
			nCverts = 3;
			nF = 0;
			nFverts = 0;
		#else
			nF = 4;
			nCverts = 4;
			nFverts = 3;
		#endif
	}
	else
	{
		cout << "Couldn't determine which type of volume used." << endl;
		exit(-1);
	}
	
	
	getline(fin, line); /*Get numbers*/
	uint nP = 0;
	uint nE = 0;
	std::size_t ptr2;
	ptr = line.find("N=");
	if(ptr!=string::npos)
	{
		ptr2 = line.find_first_not_of("0123456789",ptr+2);
		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
		nP = stoi(temp);
	}
	ptr = line.find("E=");
	if(ptr!=string::npos)
	{
		ptr2 = line.find_first_not_of("0123456789",ptr+2);
		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
		nE = stoi(temp);
	}

	#if SIMDIM == 2
		/*Skip the 3D data and read a symmetry plane*/
		uint lineNo = ceil(float(nP)/5.0)*nvar + nE+7;
		GotoLine(fin,lineNo);
		getline(fin,line);
		// cout << lineNo << endl;
		// cout << line << endl;

		ptr = line.find("N=");
		if(ptr!=string::npos)
		{
			ptr2 = line.find_first_not_of("0123456789",ptr+2);
			string temp = line.substr(ptr+2,ptr2-(ptr+2));
			
			nP = stoi(temp);
		}
		else
		{
			cout << "Number of vertices not found. Stopping." << endl;
			cout << "Line: " << lineNo << " Contents: " << line << endl;
			exit(-1);
		}

		ptr = line.find("E=");
		if(ptr!=string::npos)
		{
			ptr2 = line.find_first_not_of("0123456789",ptr+2);
			string temp = line.substr(ptr+2,ptr2-(ptr+2));
			
			nE = stoi(temp);
		}
		else
		{
			cout << "Number of elements not found. Stopping." << endl;
			cout << "Line: " << lineNo << " Contents: " << line << endl;
			exit(-1); 
		}
		// cout << nP << "  " << nE << endl;
	#endif  

	// cout << cells.numPoint << "  " << cells.numElem << endl;
	cells.reserve(nP,nE,nCverts,nF,nFverts);
	
	getline(fin,line);

	/*************** START OF VERTICES DATA *******************/
	uint varcount = 0;
	/*Get the position vectors*/
	#if SIMDIM == 2	/*Skip y component*/
		Get_Vector(fin, nP, cells.verts, 1);
	#else /*get y for 3D sim*/
		Get_Vector(fin, nP, cells.verts, 0);
	#endif

	varcount +=3;
	/*Skip variables aside from the velocity vectors*/
	for (uint ii = 3; ii < velstart ; ++ii)
	{
		if(varcount == cpstart)
		{	/*If Cp data is encountered, then read it in.*/
			if(pressOrcp == 1)
				Get_Scalar_Data(fin, nP, cells.pointCp);
			else
				Get_Scalar_Data(fin, nP, cells.pointP);
		}
		else if(varcount == densstart)
		{
			Get_Scalar_Data(fin, nP, cells.pointRho);
		}
		else
		{
			Skip_Variable(fin,nP);
		}
		varcount++;
	}
	

	#if SIMDIM == 2
		if (veltype == 1) /*Don't skip since there isnt a y vel*/
		{
			Get_Vector(fin, nP, cells.pVel, 0);
			varcount += 2;
		} 
		else /*Skip y velocity component*/
		{
			Get_Vector(fin, nP, cells.pVel, 1);
			varcount +=3;
		}
	#else /*Get the 3D velocity*/
		if(veltype == 0)
		{
			Get_Vector(fin, nP, cells.pVel, 0);
			varcount +=3;
		} 
	#endif

	uint velend;
	#if SIMDIM == 2
		if (veltype == 1)
			velend = 2;
		else
			velend = 3;
	#else
		velend = SIMDIM;
	#endif 

	/*Skip remaining variables to get to the cell data*/
	for (uint ii = 0; ii < nvar - (velstart+velend); ++ii)
	{
		if(varcount == cpstart)
		{	/*If Cp data is encountered, then read it in.*/
			if(pressOrcp == 1)
				Get_Scalar_Data(fin, nP, cells.pointCp);
			else
				Get_Scalar_Data(fin, nP, cells.pointP);
		}
		else if(varcount == densstart)
		{
			Get_Scalar_Data(fin, nP, cells.pointRho);
		}
		else
		{
			Skip_Variable(fin,nP);
		}
		varcount++;
	}

	if(varcount != nvar)
	{
		cout << "Some point data has been missed. \nCell data won't be read correctly. Stopping." << endl;
		exit(-1);
	}

	cout << "Vertex data complete. Reading cells..." << endl;

	/***************** END OF VERTICES DATA *******************/

	// cout << cells.verts.size() << "  " << cells.pointMach.size() <<  endl;
	// cout << cells.pointMach.back() << endl;
	// getline(fin,line);
	// cout << line << endl;

	/*BEGINNING OF CELL CONNECTIVITY*/
	Get_Cells(fin,nE,nCverts,cells.elems);

	cout << "Cell data complete. Closing file..." << endl;
 
	// getline(fin,line);
	// cout << line << endl;
	fin.close(); 
	
	#if SIMDIM == 2

	// cout << cells.verts.size() << endl; 
	for(uint ii = 0; ii < nE; ++ii)
	{
		for (uint jj = 0; jj < nCverts; ++jj)
		{	
			if (ii > cells.cVerts.size())
			{
				cout << "Loop attempted to access out of bounds." << endl;
				exit(-1);
			}
			if(cells.elems[ii][jj]>cells.verts.size())
			{
				cout << "Value in element list exceeds vertex list size." << endl;
				cout << ii << "  " << jj << "  " << cells.elems[ii][jj] << endl;
				exit(-1);
			}
			cells.cVerts[ii][jj] = cells.verts[cells.elems[ii][jj]];								
		}
	}

	// uint cellID = 14603;
	// for(uint jj = 0; jj < nCverts; ++jj)
	// {
	// 	cout << cells.elems[cellID][jj] + 1 << ":  " << cells.cVerts[cellID][jj][0] <<
	// 	        "  " << cells.cVerts[cellID][jj][1] << endl;
	// }
	// cout << endl;

	#endif
	
	// for(uint ii = 0; ii < 3; ++ii)
	// {	
	// 	cout << endl << "Cell: " << ii << endl;
	// 	for (auto vert: cells.elems[ii] )
	// 	{
	// 		cout << cells.verts[vert][0] << "  " << cells.verts[vert][1] 
	// 			<< "  " << cells.verts[vert][2] << endl;
	// 	}
	// }
	

	#if SIMDIM == 3
		/*Face vertices to have all faces counter-clockwise*/
		std::vector<std::vector<int>> facenum = {{1,5,6,2},{2,6,7,3},{0,3,7,4},
												{0,4,5,1},{0,1,2,3},{7,6,5,4}};
		
		// cout << nE << "  " << nF << "  " << nCverts << endl;
		// cout << cells.verts.size() << "  " << cells.elems.size() << endl;
		// cout << cells.cFaces.size() << "  " << cells.cFaces[0].size() << "  ";
		// cout << cells.cFaces[0][0].size() << endl;
		
		/*Build Face data for containement queries*/
		for(uint ii = 0; ii < nE; ++ii)
		{	
			uint jj = 0; 
			for(auto faces:facenum)
			{	
				uint kk = 0;
				for (auto vert:faces)
				{
					// cout << vert << "  ";
					if (ii > cells.cFaces.size())
					{
						cout << "Loop attempted to access out of bounds." << endl;
						exit(-1);
					}
					if(cells.elems[ii][vert]>cells.verts.size())
					{
						cout << "Value in element list exceeds vertex list size." << endl;
						cout << ii << "  " << vert << "  " << cells.elems[ii][vert] << endl;
						exit(-1);
					}
					cells.cFaces[ii][jj][kk] = cells.verts[cells.elems[ii][vert]];
					++kk;
				}
				// cout << endl;
				++jj;
			}
		}

		// uint ii = 21707;	
		// cout << endl << "Cell: " << ii << endl;
		// uint jj = 0;
		// for (std::vector<StateVecD> const& faces: cells.cFaces[ii] )
		// {	
		// 	cout << "Face: " << jj << endl;
		// 	uint kk = 0;
		// 	for (StateVecD const& verts:faces)
		// 	{
		// 		cout << facenum[jj][kk] << ":  " << verts[0] << "  " << verts[1] 
		// 			 << "  " << verts[2] << endl;
		// 		++kk;
		// 	}
		// 	++jj;
	
		// }
		
	#endif


	/*Build Cell connectivity*/
	/*Find the cell neighbours so that they can be checked 
		when a particle isn't in the cell any more*/

	/* check if a cell has 2 vertices the same in it as the checked cell
		If yes, then it's a neighbour. */
	cout << "Building cell neighbours..." << endl;
	
	
	#pragma omp parallel 
	{
		std::vector<std::vector<uint>> cNeighb = std::vector<std::vector<uint>>(nE,std::vector<uint>());
	    #pragma omp for schedule(static) nowait
		for (uint ii = 0; ii < nE; ++ii)
		{
			for (uint jj = 0; jj < nE; ++jj)
			{
				if (jj == ii)
					continue;

				uint count = 0;
				for (uint kk = 0; kk < cells.elems[ii].size(); ++kk)
				{
					if(std::find(cells.elems[jj].begin(),cells.elems[jj].end(),cells.elems[ii][kk])!=cells.elems[jj].end())
						count++;
				}

				uint thresh;
				#if SIMDIM == 2
					thresh = 1;
				#else
					thresh = 2;
				#endif

				if(count >=thresh)
					cNeighb[ii].push_back(jj);
			}
		}
		
		#pragma omp for schedule(static) ordered
    	for(int i=0; i<NTHREADS; i++)
    	{
    		#pragma omp ordered
    		cells.cNeighb.insert(cells.cNeighb.end(), cNeighb.begin(), cNeighb.end());
    	}
	       
	}

	cout << "Averaging point data to the cell..." << endl;


	StateVecD zero = StateVecD::Zero();
	/*Average data from the points to find the cell based data*/
	NormalisePressure(cells,fvar);

	Average_Point_to_Cell(cells.pVel,cells.cVel, cells.elems, zero);
	Average_Point_to_Cell(cells.pointRho,cells.cellRho,cells.elems,0.0);
	if(pressOrcp == 1)
	{
		Average_Point_to_Cell(cells.pointCp,cells.cellCp,cells.elems,0.0);
		cells.cellP = CpToPressure(cells.cellCp,fvar);
	}
	else
		Average_Point_to_Cell(cells.pointP,cells.cellP,cells.elems,0.0);
	
	// Average_Point_to_Cell(cells.pointMach,cells.cellMach, cells.elems, 0.0);

}

#if SIMDIM == 2
void Read_Radial(string input, MESH &cells)
{
	input.append("O_Mesh.plt");
	std::ifstream fin(input, std::ios::in);

	if(!fin.is_open())
	{
		cout << "Couldn't open O_Mesh.plt. Stopping." << endl;
		cout << "Path attempted: " << input << endl;
		exit(-1);
	}
	else 
	{
		cout << "Mesh file open, reading data..." << endl;
	}
	std::string line;
	getline(fin,line);
	getline(fin,line);

	/*If 2D, skip this and read the symmetry plane data.*/

	getline(fin, line); /*Get Zone line data. Tells which type of volume*/
	uint nF, nCverts, nFverts;

	nCverts = 4;
	nF = 0;
	nFverts = 0;

	getline(fin, line); /*Get numbers*/
	uint nP = 0;
	uint nE = 0;
	std::size_t ptr2;
	std::size_t ptr = line.find("N=");
	if(ptr!=string::npos)
	{
		ptr2 = line.find_first_not_of("0123456789",ptr+2);
		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
		nP = stoi(temp);
	}
	ptr = line.find("E=");
	if(ptr!=string::npos)
	{
		ptr2 = line.find_first_not_of("0123456789",ptr+2);
		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
		nE = stoi(temp);
	}

	// cout << cells.numPoint << "  " << cells.numElem << endl;
	cells.reserve(nP,nE,nCverts,nF,nFverts);
	
	getline(fin,line);

	Get_Vector(fin, nP, cells.verts, 0);

	Get_Vector(fin, nE, cells.cVel, 0); 

	/*BEGINNING OF CELL CONNECTIVITY*/
	Get_Cells(fin,nE,nCverts,cells.elems);

	// getline(fin,line);
	// cout << line << endl;
	fin.close(); 

	for(uint ii = 0; ii < nE; ++ii)
	{
		for (uint jj = 0; jj < nCverts; ++jj)
		{	
			if (ii > cells.cVerts.size())
			{
				cout << "Loop attempted to access out of bounds." << endl;
				exit(-1);
			}
			if(cells.elems[ii][jj]>cells.verts.size())
			{
				cout << "Value in element list exceeds vertex list size." << endl;
				cout << ii << "  " << jj << "  " << cells.elems[ii][jj] << endl;
				exit(-1);
			}
			cells.cVerts[ii][jj] = cells.verts[cells.elems[ii][jj]];								
		}
	}

	cout << "Building cell neighbours..." << endl;
	
	#pragma omp parallel 
	{
		std::vector<std::vector<uint>> cNeighb = std::vector<std::vector<uint>>(nE,std::vector<uint>());
	    #pragma omp for schedule(static) nowait
		for (uint ii = 0; ii < nE; ++ii)
		{
			for (uint jj = 0; jj < nE; ++jj)
			{
				if (jj == ii)
					continue;

				uint count = 0;
				for (uint kk = 0; kk < cells.elems[ii].size(); ++kk)
				{
					if(std::find(cells.elems[jj].begin(),cells.elems[jj].end(),cells.elems[ii][kk])!=cells.elems[jj].end())
						count++;
				}

				uint thresh;
				#if SIMDIM == 2
					thresh = 1;
				#else
					thresh = 2;
				#endif

				if(count >=thresh)
					cNeighb[ii].push_back(jj);
			}
		}
		
		#pragma omp for schedule(static) ordered
    	for(int i=0; i<NTHREADS; i++)
    	{
    		#pragma omp ordered
    		cells.cNeighb.insert(cells.cNeighb.end(), cNeighb.begin(), cNeighb.end());
    	}
	       
	}
}
#endif

/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/

void Write_settings(SIM &svar, FLUID &fvar)
{
	string sett = svar.outfolder;
	sett.append("/Test_Settings.txt");
	std::ofstream fp(sett, std::ios::out);

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

void Write_Mesh_Data(SIM &svar, MESH &cells)
{
	string mesh = svar.outfolder;
	mesh.append("/Mesh.plt");
	std::ofstream fm(mesh, std::ios::out);

	if(!fm.is_open())
	{ 
		cout << "Failed to open mesh output file" << endl;
		exit(-1);
	}

	#if SIMDIM == 2
		fm << "TITLE = \"2D TAU Solution\"\n";
		fm << "VARIABLES = \"x (m)\" \"z (m)\" \"x_velocity\" \"z_velocity\"\n";
		fm << "ZONE T=\"2D Solution Plane\"\n";
		fm << "VARLOCATION=([1-2]=NODAL,[3-4]=CELLCENTERED)\n";
		fm << "N=" << cells.numPoint << ", E=" << cells.numElem << ", F=FEBLOCK ET=Quadrilateral\n\n";
	#endif

	#if SIMDIM == 3
		fm << "TITLE = \"3D TAU Solution\"\n";
		fm << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\" \"x_velocity\" \"y_velocity\" \"z_velocity\"\n";
		fm << "ZONE T=\"3D Solution Plane\"\n";
		fm << "VARLOCATION=([1-3]=NODAL,[4-6]=CELLCENTERED)\n";
		fm << "N=" << cells.numPoint << ", E=" << cells.numElem << ", F=FEBLOCK ET=Brick\n\n";
	#endif

	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	uint kk = 0;
		for(uint jj = 0; jj < cells.numPoint; ++jj)
		{
			fm << cells.verts[jj][ii] << " ";
			kk++;

			if(kk == 5)
			{
				fm << "\n";
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fm << "\n";
	}

	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	uint kk = 0;
		for(uint jj = 0; jj < cells.numElem; ++jj)
		{
			fm << cells.cVel[jj][ii] << " ";
			kk++;

			if(kk == 5)
			{
				fm << "\n";
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fm << "\n";
	}

	for(uint ii = 0; ii <= cells.numElem; ++ii)
	{	
		for(auto elem:cells.elems[ii])
		{
			fm << elem+1 << " ";
		}
		fm << "\n";
	}


}

void Write_ASCII_Timestep(std::ofstream& fp, SIM &svar, State &pnp1/*, State &airP*/)
{
	 fp <<  "ZONE T=\"Particle Data\"" <<", I=" << svar.simPts << ", F=POINT" <<
    ", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
    switch(svar.outform)
    {
    	case 0:
    	{
    		for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
			{
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->xi(i) << " "; 
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
		  	fp << std::flush;
		  	break;
    	}
    	case 1:
    	{
    		for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
			{	
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->xi(i) << " ";

		        fp << p->v.norm() << " ";
		        fp << p->f.norm() << " ";
		        fp << p->rho << " "  << p->p  << "\n";
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
		  	break;
    	}
    	case 2:
    	{
    		for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
			{
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->xi(i) << " ";
		        
		        fp << p->f.norm() << " " << p->Af.norm() << " " << p->Sf.norm() << " ";
		        for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->cellV(i) << " "; 

		        fp << p->b << " " << p->theta  << "\n"; 
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
		  	break;
    	}
    }
}

void Write_ASCII_header(std::ofstream& fp, SIM &svar)
{
	switch (SIMDIM)
	{
	case 3:
		switch (svar.outform)
		{	
			case 0:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\""<< "\n";
				break;
			case 1:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\"" << "\n";
				break;
			case 2:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"A<sub>f</sub>\", \"S<sub>f</sub>\", \"S<sub>fx</sub>\", \"S<sub>fy</sub>\", \"S<sub>fz</sub>\", \"B\", \"Theta\"" << "\n";	
				break;
		}
		break;

	case 2:
		switch (svar.outform)
		{	
			case 0:
				fp << "VARIABLES = \"x (m)\", \"z (m)\""<< "\n";
				break;
			case 1:
				fp << "VARIABLES = \"x (m)\", \"z (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\"" << "\n";
				break;
			case 2:
				fp << "VARIABLES = \"x (m)\", \"z (m)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"A<sub>f</sub>\", \"S<sub>f</sub>\", \"S<sub>fx</sub>\", \"S<sub>fz</sub>\", \"B\", \"Theta\"" << "\n";
				break;
		}
		break;
	}
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