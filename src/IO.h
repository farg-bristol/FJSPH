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

using namespace std; 

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
			cout << "\tLine:" << endl;
			cout << "\t" << sline.str() << endl;
		}
	#endif
	#if SIMDIM == 3
		sline >> x(2);
		if (!sline)
		{	
			cout << "2D input provided. Please provide a 3D file." << endl;
			cout << "Incorrect line:" << endl;
			cout << sline.str() << endl;
			exit(-1);
		}
	#endif
	return x;
}

StateVecD getDVector(ifstream& In)
{
	string line;
	getline(In,line);
	std::istringstream sline(line);
	
	StateVecD x;
	sline >> x(0); sline >> x(1);

	#if SIMDIM == 2
		double temp;
		if(sline >> temp)
		{
			cout << "\tWARNING: 3D Input provided for a 2D Simulation." << endl;
			cout << "\t         The third dimension shall be ignored." << endl;
			cout << "\tLine:" << endl;
			cout << "\t" << sline.str() << endl;
		}
	#endif
	#if (SIMDIM == 3)
		sline >> x(2);
		if (!sline)
		{	
			cout << "2D input provided. Please provide a 3D file." << endl;
			cout << "Incorrect line:" << endl;
			cout << sline.str() << endl;
			exit(-1);
		}
	#endif	

	return x;
}

/*Function for a 2D Vector (e.g. Newmark Beta parameters)*/
Eigen::Vector2d getvector(ifstream& In)
{
	string line;
	getline(In,line);
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
		cout << "1: Input file   2: Output file directory.\n";
		cout << "Other inputs will be ignored." << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tWARNING: No inputs provided. Stopping... \n";
    	exit(-1);    	
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
	  		svar.outtype = getInt(in);
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
		  		svar.Jet = getvector(in); /*Defined in VLM.h. Reused here*/
		  		fvar.pPress = getDouble(in);
		  		cvar.vJet = StateVecD::Zero(); cvar.vInf = StateVecD::Zero();
		  		cvar.vJet(1) = getDouble(in);  
		  		cvar.vJet = svar.Rotate*cvar.vJet;
		  		cvar.vInf = getDVector(in);
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
		Eigen::Vector2d nb = getvector(fluid);
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
	ldouble rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
	svar.dx = pow(fvar.Simmass/rho, 1.0/double(SIMDIM));
	// cout << rho << "  " << svar.dx << endl;
	#if SIMDIM == 3
		svar.vortex.GetGamma(cvar.vInf);
	#endif
	GetAero(fvar, fvar.H);
}

void Get_Vector(ifstream &fin, int np, std::vector<StateVecD> &var, int yvalue)
{
	string line;
	for(uint dim =0; dim < SIMDIM; dim++)
	{	
		uint i = 0;
		if(SIMDIM == 2)
		{	/*If in 2D, skip the y-dimension, and read the z-dim.*/
			if(yvalue == 0)
			{
				if (dim == 1)
				{
					for(int ii = 0; ii <= floor(np/5); ++ii)
					{
						getline(fin, line);
					}
				}
			}
		}
		for(int ii = 0; ii <= floor(np/5); ++ii)
		{
			getline(fin, line);
			istringstream sline(line);

			for (uint jj = 0; jj < 5; ++jj)
			{	
				double temp;
				sline >> temp;

				if (!sline)
					break;
				else
				{
					var[i](dim) = temp;
					++i;
				}
			}
		}
	}
}

template <class T>
void Get_Scalar_Data(ifstream &fin, const int np, T &var)
{
	string line;

	uint i = 0;
	for(int ii = 0; ii <= floor(np/5); ++ii)
	{
		getline(fin, line);
		istringstream sline(line);

		for (uint jj = 0; jj < 5; ++jj)
		{	
			double temp;
			sline >> temp;

			if (!sline)
				break;
			else
			{
				var[i] = temp;
				++i;
			}
		}
	}
}

template <class T>
void Average_Point_to_Cell(std::vector<T> &pData, std::vector<T> &cData,
							const std::vector<std::vector<int>> &elems, const T zero)
{
	for(uint ii = 0; ii < elems.size(); ++ii)
	{
		T sum = zero;
		for (auto jj:elems[ii])
		{
			sum += pData[jj];
		}
		cData[ii] = sum/8.0;
	}
}

void Read_TAUMESH(MESH &cells)
{
	std::ifstream fin("sol.pval.1000.plt", std::ios::in);

	if(fin.is_open())
	{
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
					cout << "Continuing!" << endl;
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
				cout << "I can work with this. Continuing..." << endl;
				veltype = 1;
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

		std::size_t ptr = line.find("\"cp\"");
		std::size_t ptr2 = line.find("\"Mach_number\"");

		if(ptr != string::npos)
		{
			cout << "Cp data present" << endl;

		}
		if (ptr2 != string::npos)
		{
			cout << "Mach data present" << endl;
		}

		int mfirst =0;
		if(ptr != string::npos && ptr2 != string::npos)
		{
			if (ptr < ptr2)
			{
				mfirst = 0;
			}
			else
			{
				mfirst = 1;
			}
		}

		getline(fin, line); /*Get Zone line data. Tells which type of volume*/

		if(line.find("hexa")!=string::npos)
		{
			/*Do something. N verts = 8, N faces = 6, N edges = 12*/
		}

		if(line.find("tetra")!=string::npos)
		{
			/*Do something else. */
		}

		getline(fin, line); /*Get numbers*/

		ptr = line.find("N=");
		if(ptr!=string::npos)
		{
			ptr2 = line.find_first_not_of("0123456789",ptr+2);
			string temp = line.substr(ptr+2,ptr2-(ptr+2));
			
			cells.numPoint = stoi(temp);
		}
		ptr = line.find("E=");
		if(ptr!=string::npos)
		{
			ptr2 = line.find_first_not_of("0123456789",ptr+2);
			string temp = line.substr(ptr+2,ptr2-(ptr+2));
			
			cells.numElem = stoi(temp);
		}

		// cout << cells.numPoint << "  " << cells.numElem << endl;
		cells.reserve(cells.numPoint,cells.numElem);

		getline(fin,line);

		/*************** START OF VERTICES DATA *******************/
		/*Get the position vectors*/
		Get_Vector(fin, cells.numPoint, cells.verts, 0);

		/*Get velocities*/
		Get_Vector(fin, cells.numPoint, cells.pointVel, veltype);
		
		
		if(mfirst == 1)
		{	
			Get_Scalar_Data(fin, cells.numPoint, cells.pointMach);
			Get_Scalar_Data(fin, cells.numPoint, cells.pointCp);
		}
		else
		{	
			Get_Scalar_Data(fin, cells.numPoint, cells.pointCp);
			Get_Scalar_Data(fin, cells.numPoint, cells.pointMach);
		}
		/***************** END OF VERTICES DATA *******************/

		// cout << cells.verts.size() << "  " << cells.pointMach.size() <<  endl;
		// cout << cells.pointMach.back() << endl;
		

		/*BEGINNING OF CELL CONNECTIVITY*/
		for(uint ii = 0; ii <= cells.numElem; ++ii)
		{
			getline(fin, line);
			istringstream sline(line);

			for (uint jj = 0; jj < 8; ++jj)
			{	
				int temp;
				sline >> temp;
				cells.elems[ii][jj] = (temp-1);								
			}
		}

		fin.close();
	}
	else
	{
		cout << "Couldn't open sol.pval.1000.plt. Stopping" << endl;
		exit(-1);
	}

	// for(uint ii = 0; ii <= cells.numElem; ++ii)
	// {
	// 	for (uint jj = 0; jj < 8; ++jj)
	// 	{	
	// 		cells.elemverts[ii][jj] = cells.verts[cells.elems[ii][jj]];								
	// 	}
	// }

	// for(uint ii = 0; ii < 3; ++ii)
	// {	
	// 	cout << endl << "Cell: " << ii << endl;
	// 	for (auto vert: cells.elems[ii] )
	// 	{
	// 		cout << cells.verts[vert][0] << "  " << cells.verts[vert][1] 
	// 			<< "  " << cells.verts[vert][2] << endl;
	// 	}
	// }
	
	/*Face vertices to have all faces counter-clockwise*/
	std::vector<std::vector<int>> facenum = {{1,5,6,2},{2,6,7,3},{0,3,7,4},
											{0,4,5,1},{0,1,2,3},{7,6,5,4}};

	/*Build Face data for containement queries*/
	for(uint ii = 0; ii <= cells.numElem; ++ii)
	{	
		uint jj = 0; 
		for(auto faces:facenum)
		{	
			uint kk = 0;
			for (auto vert:faces)
			{
				cells.elemfaces[ii][jj][kk] = cells.verts[cells.elems[ii][vert]];
				++kk;
			}
			++jj;
		}
	}


	/*Build Cell connectivity*/
	/*Find the cell neighbours so that they can be checked 
		when a particle isn't in the cell any more*/

	/* check if a cell has 2 vertices the same in it as the checked cell
		If yes, then it's a neighbour. */
	cout << "Building cell neighbours..." << endl;

	for (uint ii = 0; ii < cells.numElem; ++ii)
	{
		cells.elemneighb.emplace_back();
		for (uint jj = 0; jj < cells.numElem; ++jj)
		{
			if (jj == ii)
				continue;

			int count = 0;
			for (uint kk = 0; kk < cells.elems[ii].size(); ++kk)
			{
				if(std::find(cells.elems[jj].begin(),cells.elems[jj].end(),cells.elems[ii][kk])!=cells.elems[jj].end())
					count++;
			}

			if(count >=2)
				cells.elemneighb[ii].emplace_back(jj);
		}
	}

	cout << "Averaging point data to the cell..." << endl;
	StateVecD zero = StateVecD::Zero();
	/*Average data from the points to find the cell based data*/
	Average_Point_to_Cell(cells.pointVel,cells.cellVel, cells.elems, zero);
	Average_Point_to_Cell(cells.pointCp,cells.cellCp, cells.elems, 0.0);
	Average_Point_to_Cell(cells.pointMach,cells.cellMach, cells.elems, 0.0);

}


/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/

void Write_settings(SIM &svar, FLUID &fvar)
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

void Write_ASCII_Timestep(std::ofstream& fp, SIM &svar, State &pnp1)
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
		  	fp << std::flush;
		  	break;
    	}
    	case 2:
    	{
    		for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
			{
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->xi(i) << " ";
		        
		        fp << p->f.norm() << " " << p->Af(0) << " " << p->Sf.norm() << " ";
		        for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->Af(i) << " "; 

		        fp << p->b << " " << p->theta  << "\n"; 
		  	}  	
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
				fp << "VARIABLES = \"x (m)\", \"y (m)\""<< "\n";
				break;
			case 1:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"v (m/s)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"<greek>r</greek> (kg/m<sup>-3</sup>)\", \"P (Pa)\"" << "\n";
				break;
			case 2:
				fp << "VARIABLES = \"x (m)\", \"y (m)\", \"a (m/s<sup>-1</sup>)\", " << 
			"\"A<sub>f</sub>\", \"S<sub>f</sub>\", \"S<sub>fx</sub>\", \"S<sub>fy</sub>\", \"B\", \"Theta\"" << "\n";
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