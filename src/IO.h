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
#include "TauIO.h"

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
	if( stat( pathname.c_str(), &info ) != 0 )
  	{	/*Output directory doesn't exist, so create it*/
		pathname.insert(0,"\"");
		pathname.append("\"");
  		string cmd = "mkdir ";
	  	cmd.append(pathname);
	    if(system(cmd.c_str()))
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

  	/*Create h5 folder*/
  	if(svar.outtype == 0)
  	{
  		string file = pathname;
  		file.append("fuel.szplt.szdat");
  		struct stat info;
  		if(stat( file.c_str(), &info ) == 0)
  		{
	  		string cmd = "exec rm -r \"";
	  		cmd.append(pathname);
	  		cmd.append("\"*.szplt.sz*");
	  		if(system(cmd.c_str()))
	  		{
		    	cout << "System command failed to execute." << endl;
		    	exit(-1);
		    }
		}
  	}
  	else if(svar.outtype == 2)
  	{
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
		    	exit(-1);
		    }
		}
  	}

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
	getline(In,line);
	lineno++;
	size_t ptr = line.find_first_of(' ');
	string result = line.substr(0,ptr);
	return result; 
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

void GetAero(AERO &avar, const FLUID& fvar, const ldouble rad)
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

void GetInput(int argc, char **argv, SIM &svar, FLUID &fvar, AERO& avar)
{
	uint justPost = 0;
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
	  		svar.boutform = getInt(in, lineno);
	  		svar.subits = getInt(in, lineno);
	  		svar.nmax = getInt(in, lineno);	
	  		svar.xyPART = getIVector(in, lineno);
	  		svar.Pstep = getDouble(in, lineno);
	  		svar.Bstep = getDouble(in, lineno);
	  		svar.Bcase = getInt(in, lineno);
	  		avar.acase = getInt(in, lineno);
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
	  		else if(svar.Bcase > 1 && svar.Bcase < 8)
	  		{	
	  			StateVecD angles = getDVector(in, lineno);
	  			angles = angles *M_PI/180;
	  			svar.Rotate = GetRotationMat(angles);
	  			svar.Transp = svar.Rotate.transpose();
		  		svar.Jet = getvector(in, lineno); /*Defined in VLM.h. Reused here*/
		  		fvar.pPress = getDouble(in, lineno);
		  		avar.vJet = StateVecD::Zero(); avar.vInf = StateVecD::Zero();
		  		avar.vJet(1) = getDouble(in, lineno);  
		  		avar.vJet = svar.Rotate*avar.vJet;
		  		avar.vInf = getDVector(in, lineno);
		  		
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
	  		svar.afterSim = getInt(in, lineno);
	  		svar.cellSize = getDouble(in, lineno);
	  		svar.postRadius = getDouble(in, lineno);

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
  		if(svar.Bcase == 6 || svar.Bcase == 4)
  		{
	  		fvar.gasVel = getDouble(fluid, lineno);
	  		fvar.gasPress = getDouble(fluid, lineno);
	  		fvar.T = getDouble(fluid, lineno);

	  		if(svar.Bcase == 6)
	  		{
		  		svar.meshfile = getString(fluid, lineno);
		  		svar.solfile = getString(fluid,lineno);
	  		}
  		}
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
	  	svar.addcount = 0;
	  	svar.dt = 2E-010; 		/*Initial timestep*/
	  	svar.t = 0.0;				/*Total simulation time*/
	  	fvar.HSQ = fvar.H*fvar.H; 

		fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/
		svar.Bclosed = 0; 		/*Boundary begins open*/
		svar.psnPts = 0; 		/*Start with no pison points*/
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


// #if SIMDIM == 2
// void Read_Radial(string input, MESH &cells)
// {
// 	input.append("O_Mesh.plt");
// 	std::ifstream fin(input, std::ios::in);
// 	ZONE zone;

// 	if(!fin.is_open())
// 	{
// 		cout << "Couldn't open O_Mesh.plt. Stopping." << endl;
// 		cout << "Path attempted: " << input << endl;
// 		exit(-1);
// 	}
// 	else 
// 	{
// 		cout << "Mesh file open, reading data..." << endl;
// 	}
// 	std::string line;
// 	getline(fin,line);
// 	getline(fin,line);

// 	/*If 2D, skip this and read the symmetry plane data.*/

// 	getline(fin, line); /*Get Zone line data. Tells which type of volume*/
// 	uint nF, nCverts, nFverts;

// 	nCverts = 4;
// 	nF = 0;
// 	nFverts = 0;

// 	getline(fin, line); /*Get numbers*/
	
// 	std::size_t ptr2;
// 	std::size_t ptr = line.find("N=");
// 	if(ptr!=string::npos)
// 	{
// 		ptr2 = line.find_first_not_of("0123456789",ptr+2);
// 		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
// 		zone.nP = stoi(temp);
// 	}
// 	ptr = line.find("E=");
// 	if(ptr!=string::npos)
// 	{
// 		ptr2 = line.find_first_not_of("0123456789",ptr+2);
// 		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
// 		zone.nE = stoi(temp);
// 	}

// 	// cout << cells.numPoint << "  " << cells.numElem << endl;
// 	// cells.reserve(zone.nP,zone.nE,nCverts,nF,nFverts);
	
// 	getline(fin,line);

// 	cells.verts = Get_Vector(fin, zone, 0);

// 	cells.cVel = Get_Vector(fin, zone, 0); 

// 	Get_2DCells(fin,cells.verts,zone,cells.elems, cells.cVerts);

// 	// getline(fin,line);
// 	// cout << line << endl;
// 	fin.close(); 

// 	cout << "Building cell neighbours..." << endl;
	
// 	#pragma omp parallel 
// 	{
// 		std::vector<std::vector<uint>> cNeighb = std::vector<std::vector<uint>>(zone.nE,std::vector<uint>());
// 	    #pragma omp for schedule(static) nowait
// 		for (uint ii = 0; ii < zone.nE; ++ii)
// 		{
// 			for (uint jj = 0; jj < zone.nE; ++jj)
// 			{
// 				if (jj == ii)
// 					continue;

// 				uint count = 0;
// 				for (uint kk = 0; kk < cells.elems[ii].size(); ++kk)
// 				{
// 					if(std::find(cells.elems[jj].begin(),cells.elems[jj].end(),cells.elems[ii][kk])!=cells.elems[jj].end())
// 						count++;
// 				}

// 				uint thresh;
// 				#if SIMDIM == 2
// 					thresh = 1;
// 				#else
// 					thresh = 2;
// 				#endif

// 				if(count >=thresh)
// 					cNeighb[ii].push_back(jj);
// 			}
// 		}
		
// 		#pragma omp for schedule(static) ordered
//     	for(int i=0; i<NTHREADS; i++)
//     	{
//     		#pragma omp ordered
//     		cells.cNeighb.insert(cells.cNeighb.end(), cNeighb.begin(), cNeighb.end());
//     	}
	       
// 	}
// }
// #endif

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

void Write_ASCII_Timestep(std::ofstream& fp, SIM &svar, const State &pnp1, 
	const uint bwrite, const uint start, const uint end, const State &airP )
{
	if(bwrite == 1)
	 	fp <<  "ZONE T=\"Boundary Data\"";
	else
		fp <<  "ZONE T=\"Particle Data\"";

    fp <<", I=" << end - start << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";

    switch(svar.outform)
    {
    	case 0:
    	{
    		for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
			{
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->xi(i) << " "; 
				fp << "\n";  
		  	}

		  	if (airP.size() > 0 )
		  	{
			  	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
			    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
			  	for(auto p:airP)
			  	{
			  		for(uint i = 0; i < SIMDIM; ++i)
			        	fp << p.xi(i) << " "; 
					fp << "\n"; 
			  	}
		  	}
		  	fp << std::flush;
		  	break;
    	}
    	case 1:
    	{
    		for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
			{	
				for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->xi(i) << " ";

		        fp << p->v.norm() << " ";
		        fp << p->f.norm() << " ";
		        fp << p->rho << " "  << p->p  << "\n";
		  	}

		  	if (airP.size() > 0 )
		  	{
			  	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
			    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";

			  	for(auto p:airP)
			  	{
			  		for(uint i = 0; i < SIMDIM; ++i)
			        	fp << p.xi(i) << " "; 

			        fp << p.v.norm() << " ";
			        fp << p.f.norm() << " ";
			        fp << p.rho << " "  << p.p  << "\n";
			  	}
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
		        
		        fp << p->f.norm() << " " << p->Af.norm() << " " << p->cellP << " ";
		        for(uint i = 0; i < SIMDIM; ++i)
		        	fp << p->cellV(i) << " "; 

		        fp << p->b << " " << p->theta  << "\n"; 
		  	}  

		  	if (airP.size() > 0 )
		  	{
			  	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
			    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
			    
			  	for(auto p:airP)
			  	{
			  		for(uint i = 0; i < SIMDIM; ++i)
			        	fp << p.xi(i) << " "; 

			      	fp << p.f.norm() << " " << p.Af.norm() << " " << p.Sf.norm() << " ";
			        for(uint i = 0; i < SIMDIM; ++i)
			        	fp << p.cellV(i) << " "; 

			        fp << p.b << " " << p.theta  << "\n"; 
			  	}
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