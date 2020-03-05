#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H
	

#include "Var.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <dirent.h>

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

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
  	string input = svar.infolder;
  	pathname.append("/");
  	pathname.append(input);
  	pathname.append(svar.outfolder);
  	pathname.append("/");
  
  	/*Check for output file name*/		
	check_folder(pathname);
	
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
			if(svar.restart == 0)
			{	// Delete if not restarting
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
			else
			{
				// Check if there is the constructed szplt
				file = pathname;
				file.append("Fuel.szplt");
				if(stat( file.c_str(), &info ) != 0)
		  		{
		  			// File needs creating.
		  			string cmd = "exec szcombine \"";
			  		cmd.append(file);
			  		cmd.append("\"");

#ifdef DEBUG
			  		dbout << "Attempting to combine szplt." << endl;
			  		dbout << "Command: " << cmd << endl;
#endif
			  		if(system(cmd.c_str()))
			  		{
				    	cout << "System command failed to execute." << endl;
				    	cout << "Command: " << cmd << endl;
				    	exit(-1);
				    }
		  		}

		  		file = pathname;
				file.append("Boundary.szplt");
				if(stat( file.c_str(), &info ) != 0)
		  		{
		  			// File needs creating.
		  			string cmd = "exec szcombine \"";
			  		cmd.append(file);
			  		cmd.append("\"");

#ifdef DEBUG
			  		dbout << "Attempting to combine szplt." << endl;
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
		}
		else if(svar.restart == 1)
		{
			cout << "No previous simulation file. Cannot restart." << endl;
			exit(-1);
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

uint index(uint ii, uint jj, uint nPts)
{
	return(ii*nPts + jj);
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(uint ii=0; ii < num - 1; ++ii){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


void Write_Input(SIM& svar, FLUID& fvar, AERO& avar, StateVecD& angle)
{
	// Open the file.
	string file = svar.outfolder;
	file.append("Settings.dat");
	uint width = 50;
	std::ofstream sett(file);

	sett << svar.framet << setw(width) << "#Frame time interval" << endl;
	sett << svar.Nframe << setw(width) << "#Number of frames" << endl;
	sett << svar.outframe << setw(width) << "#Output frame info" << endl;
	sett << svar.outtype << setw(width) << "#Output data type" << endl;
	sett << svar.outform << setw(width) << "#Output contents" << endl;
	sett << svar.boutform << setw(width) << "#Boundary time output" << endl;
	sett << svar.gout << setw(width) << "#Ghost particle output" << endl;
	sett << svar.subits << setw(width) << "#Maximum sub iterations" << endl;
	sett << svar.nmax << setw(width) << "#Maximum number of particles" << endl;
	sett << svar.afterSim << setw(width) << "#Do post processing?" << endl;
	sett << svar.cellSize << setw(width) << "#Post processing mesh size" << endl;
	sett << svar.postRadius << setw(width) << "#Post processing support radius" << endl;
	sett << svar.Pstep << setw(width) << "#Particle initial spacing" << endl;
	sett << svar.Bstep << setw(width) << "#Boundary spacing factor" << endl;
	sett << svar.Bcase << setw(width) << "#Boundary case" << endl;
	sett << avar.acase << setw(width) << "#Aerodynamic case" << endl;
	sett << svar.ghost << setw(width) << "#Ghost particles?" << endl;
	sett << svar.Start(0) << " " << svar.Start(1);
#if SIMDIM == 3
	sett << " " << svar.Start(2);
#endif
	sett << setw(width) << "#Fluid start position" << endl;
	if(svar.Bcase < 2)
	{
		sett << svar.xyPART(0) << " " << svar.xyPART(1);
#if SIMDIM == 3
		sett << " " << svar.xyPART(2);
#endif
		sett << setw(width) << "#Particle counts" << endl;
		sett << svar.Box(0) << " " << svar.Box(1);
#if SIMDIM == 3
		sett << " " << svar.Box(2);
#endif
		sett << setw(width) << "#Box dimensions" << endl << endl;
		sett << fvar.pPress << setw(width) << "#Pipe Pressure" << endl;
	}
	else
	{
		sett << angle(0) << " " << angle(1);
#if SIMDIM == 3
		sett << " " << angle(2);
#endif
		sett << "  #Fluid start rotation" << endl;
		sett << svar.Jet(0) << " " << svar.Jet(1) << setw(width) << "#Jet dimensions" << endl;
		sett << fvar.pPress << setw(width) << "#Ghost particles?" << endl;
		sett << avar.vJet(1) << setw(width) << "#Jet velocity" << endl;
		sett << avar.vInf(0) << " " << avar.vInf(1);
#if SIMDIM == 3
		sett << " " << avar.vInf(2);
#endif
		sett << "  #Freestream velocity" << endl;

		if(avar.acase == 2 || avar.acase == 3)
  		{
  			sett << avar.a << setw(width) << "#a" << endl;
  			sett << avar.h1 << setw(width) << "h1" << endl;
  			sett << avar.b << setw(width) << "#b" << endl;
  			sett << avar.h2 << setw(width) << "#h2" << endl;
  		}
	}

	/*End of settings write*/
	sett.close();

	/*Write fluid file now*/
	file = svar.outfolder;
	file.append("Fluid.dat");

	std::ofstream fluid(file);

	fluid << svar.beta << " " << svar.gamma << setw(width) << "#Newmark-Beta terms" << endl;
	fluid << fvar.Hfac << setw(width) << "#Support radius factor" << endl;
	fluid << fvar.alpha << setw(width) << "#Artificial viscosity" << endl;
	fluid << fvar.contangb << setw(width) << "#Contact angle" << endl;
	fluid << fvar.rho0 << setw(width) << "#Fluid density" << endl;
	fluid << fvar.rhog << setw(width) << "#Gas density" << endl;
	fluid << fvar.Cs << setw(width) << "#Speed of sound" << endl;
	fluid << fvar.mu << setw(width) << "#Fluid viscosity" << endl;
	fluid << fvar.mug << setw(width) << "#Gas viscosity" << endl;
	fluid << fvar.sig << setw(width) << "#Surface Tension" << endl;
	fluid << fvar.gasVel << setw(width) << "#Gas reference velocity" << endl;
	fluid << fvar.gasPress << setw(width) << "#Gas reference pressure" << endl;
	fluid << fvar.T << setw(width) << "#Gas reference temperature" << endl;
	fluid << svar.meshfile << "   #Mesh face file" << endl;
	fluid << svar.solfile << "   #Mesh Solution file" << endl;
	fluid << svar.scale << setw(width) << "#Mesh scale" << endl;

	fluid.close();
}

#endif