#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H
	

#include "Var.h"
#include "BinaryIO.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <dirent.h>

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
	#define stat _stat
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

int Combine_SZPLT(string& file)
{
	string cmd = "exec szcombine \"";
	cmd.append(file);
	cmd.append("\"");

#ifdef DEBUG
	dbout << "Attempting to combine szplt." << endl;
	dbout << "Command: " << cmd << endl;
#endif

	cout << "Combining szplt: " << file << endl;
	if(system(cmd.c_str()))
	{
    	cout << "Could not combine szplt file." << endl;
    	cout << "Command: " << cmd << endl;
    	return -1;
	}
	return 0;
}

/*Open an output directory in the name of the input file, under Outputs*/
int MakeOutputDir(int argc, char *argv[], SIM& svar)
{
  	/*Check output folder for any prexisting files*/
  	if(svar.outtype == 0)
  	{
  		string file = svar.outfolder;
  		file.append("Fuel.szplt.szdat");
#ifdef DEBUG
  		dbout << "Checking for existence of previous szplt files." << endl;
  		dbout << "Path: " << file << endl;
#endif
  		struct stat info;
  		if(stat( file.c_str(), &info ) == 0)
  		{
  			// cout << "Split files exist..." << endl;

			if(svar.restart == 0)
			{	// Delete if not restarting
		  		string cmd = "exec rm -r \"";
		  		cmd.append(svar.outfolder);
		  		cmd.append("\"*.szplt.sz*");

		  		// cout << cmd << endl;
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
				file = svar.outfolder;				
				file.append("Fuel.szplt");
				if(stat( file.c_str(), &info ) != 0)
		  		{	/*Fuel.szplt does not exist - needs creating*/
					if(Combine_SZPLT(file) == -1)
						exit(-1);
					
					if(svar.Bcase != 0 && svar.Bcase != 3 && svar.Bcase !=4)
					{
			  			file = svar.outfolder;
						file.append("Boundary.szplt");
						if(stat( file.c_str(), &info ) != 0)
				  		{ /*Boundary.szplt also doesn't exist - needs creating*/
							
					  			string szplt = file;
								szplt.append(".szdat");
								if(stat( szplt.c_str(), &info ) == 0)
								{
									if(Combine_SZPLT(file) == -1)
										exit(-1);
								}
								else
								{	//Boundaries should exist but don't.
									cout << "Boundary files cannot be found. Stopping" << endl;
									exit(-1);
								}
						}
			  		}
		  		}
		  		else
		  		{	// Both the combined and split szplt exist - check which is more recent.
		  			// cout << "Both file types exist" << endl << endl;
		  			auto sztime = info.st_mtime;
		  			file.append(".szdat");
					if(stat(file.c_str(), &info) == 0)
					{
						auto szdattime = info.st_mtime;

						cout << sztime << "   " << szdattime << endl;
 
						if (szdattime > (sztime+10))
						{
							// combine the files
							file = svar.outfolder;
							file.append("Fuel.szplt");
							if(Combine_SZPLT(file) == -1)
								exit(-1);

							if(svar.Bcase != 0 && svar.Bcase != 3 && svar.Bcase !=4)
							{	//Check if files exist
								file = svar.outfolder;
								file.append("Boundary.szplt.szdat");
								if(stat( file.c_str(), &info ) == 0)
						  		{
									file = svar.outfolder;
									file.append("Boundary.szplt");
									if(Combine_SZPLT(file) == -1)
										exit(-1);
								}
								else
								{
									cout << "Boundary files cannot be found. Stopping" << endl;
									exit(-1);
								}
							}

						}
					}

					// Do the same check for the boundary file
					file = svar.outfolder;
					file.append("Boundary.szplt");
					if(stat(file.c_str(),&info) == 0)
					{   //File exists, so get the modify time
						auto bsztime = info.st_mtime;

						file.append(".szdat");

						if(stat( file.c_str(), &info ) == 0)
				  		{   // data file also exists, so see which is more recent.
				  			auto bszdattime = info.st_mtime;

				  			if( bszdattime > bsztime)
				  			{
								file = svar.outfolder;
								file.append("Boundary.szplt");
								if(Combine_SZPLT(file) == -1)
									exit(-1);
							}
						}
					}
					else
					{
						file.append(".szdat");
						if(stat( file.c_str(), &info ) == 0)
				  		{   // data file exists
				  			file = svar.outfolder;
							file.append("Boundary.szplt");
							if(Combine_SZPLT(file) == -1)
								exit(-1);
				  		}
					}
		  		}	  		

			}
		}
		else if(svar.restart == 1)
		{
			string file = svar.outfolder;
	  		file.append("Fuel.szplt");

	  		struct stat info;
	  		if(stat( file.c_str(), &info ) == 0)
	  		{
	  			cout << "Fuel solution file exist..." << endl;
	  		}
	  		else
	  		{
				cout << "No previous fuel simulation file. Cannot restart." << endl;
				exit(-1);
			}

			file = svar.outfolder;
	  		file.append("Boundary.szplt");

	  		if(svar.Bcase != 0 && svar.Bcase != 3 && svar.Bcase != 4)
	  		{
		  		if(stat( file.c_str(), &info ) == 0)
		  		{
		  			cout << "Boundanry solution files exist..." << endl;
		  		}
		  		else
		  		{
					cout << "No previous boundary simulation file. Cannot restart." << endl;
					exit(-1);
				}
			}
		}
  	}

	return 0;
}

int getInt(ifstream& In, uint& lineno, string const& name)
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

real getDouble(ifstream& In, uint& lineno, string const& name)
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

std::string getString(ifstream& In, uint& lineno, string const& name)
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

StateVecI getIVector(ifstream& In, uint& lineno, string const& name)
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
			cout << "\t" << sline.str() << endl << endl;
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
		cout << sline.str() << endl << endl;
		exit(-1);
	}
#ifdef DEBUG
	dbout << name << ": " << x(0) << "  " << x(1) << "  " << x(2) << endl;
#endif	
#endif

	return x;
}

StateVecD getDVector(ifstream& In, uint& lineno, string const& name)
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
Eigen::Vector2d getvector(ifstream& In, uint& lineno, string const& name)
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

void CheckContents(void* const& inputHandle, SIM& svar, int32_t& numZones, double& time)
{
	int32_t numVars, I;
	I = tecDataSetGetNumVars(inputHandle, &numVars);

	std::ostringstream outputStream;
    for (int32_t var = 1; var <= numVars; ++var)
    {
        char* name = NULL;
        I = tecVarGetName(inputHandle, var, &name);
        outputStream << name;
        if (var < numVars)
            outputStream << ',';
        tecStringFree(&name);
    }

	if (I == -1)
	{
		cout << "Couldn't obtain number of variables. Stopping" << endl;
		exit(-1);
	}

	cout << "Number of variables in output file: " << numVars << endl;

	uint dataType = 0;
#if SIMDIM == 2
	if( outputStream.str() == "X,Z")
		dataType = 0;

	else if (outputStream.str() == "X,Z,rho,Rrho,m,v,a")
		dataType = 1;

	else if (outputStream.str() == "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b")
		dataType = 2;	

	else if (outputStream.str() == "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_Vx,Cell_Vz,Cell_P,Cell_ID" )
		dataType = 3;

	else if (outputStream.str() == "X,Z,rho,Rrho,m,v,a,Neighbours,Aero")
		dataType = 4;

	else if (outputStream.str() == "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_ID")
		dataType = 5;

#else
	if( outputStream.str() == "X,Y,Z")
		dataType = 0;

	else if (outputStream.str() == "X,Y,Z,rho,Rrho,m,v,a")
		dataType = 1;

	else if (outputStream.str() == "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b")
		dataType = 2;	
	else if (outputStream.str() == "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_Vx,Cell_Vy,Cell_Vz,Cell_P,Cell_ID")
		dataType = 3;

	else if (outputStream.str() == "X,Y,Z,rho,Rrho,m,v,a,Neighbours,Aero")
		dataType = 4;

	else if (outputStream.str() == "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_ID")
		dataType = 5;
	
#endif
	cout << "Variables:  " << outputStream.str() << "  Data type: " << dataType << endl;

	if(dataType == 0 || dataType == 1 || dataType == 4)
	{
		cout << "Cannot restart. File does not contain the correct information." << endl;
		exit(-1);
	}

    if(dataType!= svar.outform)
    {
    	cout << "Data is different from that in the restfings file." << endl;
    	cout << "Changing to the file data type: " << dataType << endl;
    	svar.outform = dataType;
    	if(svar.outform == 2)
		{
			if(svar.Asource == 1 || svar.Asource == 2)
			{
				cout << "No cell information. Cannot restart" << endl;
				exit(-1);
			}
		}
    }

    I = tecDataSetGetNumZones(inputHandle, &numZones);
    cout << "Number of zones in file: " << numZones << endl;

    I = tecZoneGetSolutionTime(inputHandle, numZones, &time);
    cout << "Latest zone time: " << time << endl << endl;
}

void GetYcoef(AERO& avar, const FLUID& fvar, const real diam)
{
	#if SIMDIM == 3
		avar.L = diam * std::cbrt(3.0/(4.0*M_PI));
	#endif
	#if SIMDIM == 2
		avar.L = diam/sqrt(M_PI);
	#endif

	// avar.L = diam / 2.0; pow(diam,1.25)/2.0;

	
	avar.td = (2.0*fvar.rho0*pow(avar.L,SIMDIM-1))/(avar.Cd*fvar.mu);

	avar.omega = sqrt((avar.Ck*fvar.sig)/(fvar.rho0*pow(avar.L,SIMDIM))-1.0/pow(avar.td,2.0));

	avar.tmax = -2.0 *(atan(sqrt(pow(avar.td*avar.omega,2.0)+1)
					+avar.td*avar.omega) - M_PI)/avar.omega;

	avar.Cdef = 1.0 - exp(-avar.tmax/avar.td)*(cos(avar.omega*avar.tmax)+
		1/(avar.omega*avar.td)*sin(avar.omega*avar.tmax));
	avar.ycoef = 0.5*avar.Cdef*(avar.Cf/(avar.Ck*avar.Cb))*(avar.rhog*avar.L)/fvar.sig;

	//cout << avar.ycoef << "  " << avar.Cdef << "  " << avar.tmax << "  " << endl;
}


#endif