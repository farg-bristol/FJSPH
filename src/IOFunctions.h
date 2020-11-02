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

void Combine_SZPLT(string& file)
{
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

int MakeOutputDir(int argc, char *argv[], SIM& svar)
{
	char cCurrentPath[FILENAME_MAX];
	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
		return errno;

	/*Open an output directory in the name of the input file, under Outputs*/
	string pathname = cCurrentPath;
  	pathname.append("/");
  	pathname.append(svar.infolder);
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
  			// cout << "Split files exist..." << endl;

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
		  		{	/*Fuel.szplt does not exist - needs creating*/
					Combine_SZPLT(file);
					
					if(svar.Bcase != 0 && svar.Bcase !=4)
					{
			  			file = pathname;
						file.append("Boundary.szplt");
						if(stat( file.c_str(), &info ) != 0)
				  		{ /*Boundary.szplt also doesn't exist - needs creating*/
							
					  			string szplt = file;
								szplt.append(".szdat");
								if(stat( szplt.c_str(), &info ) == 0)
								{
									Combine_SZPLT(file);
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

						// cout << sztime << "   " << szdattime << endl;
 
						if (szdattime > sztime)
						{
							// combine the files
							file = pathname;
							file.append("Fuel.szplt");
							Combine_SZPLT(file);

							if(svar.Bcase != 0 && svar.Bcase !=4)
							{	//Check if files exist
								file = pathname;
								file.append("Boundary.szplt.szdat");
								if(stat( file.c_str(), &info ) == 0)
						  		{
									file = pathname;
									file.append("Boundary.szplt");
									Combine_SZPLT(file);
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
					file = pathname;
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
								file = pathname;
								file.append("Boundary.szplt");
								Combine_SZPLT(file);
							}
						}
					}
					else
					{
						file.append(".szdat");
						if(stat( file.c_str(), &info ) == 0)
				  		{   // data file exists
				  			file = pathname;
							file.append("Boundary.szplt");
							Combine_SZPLT(file);
				  		}
					}
		  		}	  		

			}
		}
		else if(svar.restart == 1)
		{
			string file = pathname;
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

			file = pathname;
	  		file.append("Boundary.szplt");

	  		if(svar.Bcase != 0 && svar.Bcase != 4)
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

void CheckContents(void* const& inputHandle, SIM& svar)
{
	int32_t numVars, numZones, I;
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

	if( outputStream.str() == "X,Z" || outputStream.str() == "X,Y,Z")
	{
		dataType = 0;
	}
	else if (outputStream.str() == "X,Z,rho,Rrho,m,v,a" ||
			outputStream.str() == "X,Y,Z,rho,Rrho,m,v,a")
	{
		dataType = 1;
	}
	else if (outputStream.str() == "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b")
	{
		dataType = 2;	
	}
	else if (outputStream.str() == "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_Vx,Cell_Vz,Cell_P,Cell_ID" ||
		outputStream.str() == "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_Vx,Cell_Vy,Cell_Vz,Cell_P,Cell_ID")
	{
		dataType = 3;
	}
	else if (outputStream.str() == "X,Z,rho,Rrho,m,v,a,Neighbours,Aero" ||
		outputStream.str() == "X,Y,Z,rho,Rrho,m,v,a,Neighbours,Aero")
	{
		dataType = 4;
	}
	else if (outputStream.str() == "X,Z,rho,Rrho,m,v_x,v_z,a_x,a_z,b,Cell_ID" ||
			outputStream.str() == "X,Y,Z,rho,Rrho,m,v_x,v_y,v_z,a_x,a_y,a_z,b,Cell_ID")
	{
		dataType = 5;
	}

	cout << outputStream.str() << "  " << dataType << endl;

	if(dataType == 0 || dataType == 1 || dataType == 4)
	{
		cout << "Cannot restart. File does not contain the correct information." << endl;
		exit(-1);
	}

    if(dataType!= svar.outform)
    {
    	cout << "Data is different from that in the settings file." << endl;
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
    if(numZones != static_cast<INTEGER4>(svar.frame+1))
    {
    	cout << "Mismatch of number of frames and number of zones in output file." << endl;
    	cout << "Reading the last zone in file." << endl;
    	svar.frame = numZones - 1;
    }

}

/*Make a guess on how big the array is going to be (doesn't need to be totally exact)*/
int ParticleCount(SIM &svar)
{
	int partCount = 0;
	real step = svar.Pstep*svar.Bstep;
	int Ny = ceil(svar.Box(1)/step);
	int Nx = ceil(svar.Box(0)/step);

	#if(SIMDIM == 3)
		int Nz = ceil(svar.Box(2)/step);
	#endif

		if(svar.Bcase == 0)
			partCount += svar.simPts; /*Simulation pn*/
		else if (svar.Bcase == 1)
		{
			#if (SIMDIM == 3)
				partCount = Nx*Nz + 2*Nx*Ny + 2*Ny*Nz; /*Boundary particles*/
			#else
				partCount = 2*Ny + Nx; /*Boundary particles*/
			#endif
			partCount += svar.simPts; /*Simulation pn*/
		}
		else if (svar.Bcase == 2)
		{
			#if SIMDIM == 3
			{
				real holeD = svar.Jet(0)+4*svar.Pstep; /*Diameter of hole (or width)*/

	            /*Find the points on the side of the pipe (Assume a circle)*/
	            float dtheta = atan((step)/(0.5*holeD));
	            Ny = ceil(svar.Jet(1)/step);
	            int holeWall = ceil(2*M_PI/dtheta)*Ny;

				/*Simulation Points*/
				int simCount = 0;
				real jetR = 0.5*(svar.Jet(0));
				
				/*Do the centerline of points*/
				for (real z = -jetR; z <= jetR; z+= svar.dx)
					simCount++;

				for (real x = svar.dx; x < jetR ; x+=svar.dx)
				{ /*Do the either side of the centerline*/
					for (real z = -jetR; z <= jetR; z+= svar.dx)
					{	/*If the point is inside the hole diameter, add it*/
						if(((x*x)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
			    			simCount += 2;
					}
				}

				/*Need to add the pn already present*/
				int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

				partCount = holeWall + simPts;
			}
			#else
			{
	            /*Find the points on the side of the pipe*/
	            Ny = ceil((svar.Jet(1)*3)/step);
	            int holeWall = 2*Ny;

				/*Simulation Points*/
				int simCount = 0;
				real jetR = 2*(0.5*(svar.Jet(0)));
				for (real x = -jetR; x <= jetR; x+= svar.dx)
					simCount++;

				/*Need to add the pn already present*/
				int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

				partCount = holeWall + simPts;
			}
			#endif
		}
		else if(svar.Bcase == 3)
		{	
			#if SIMDIM == 3
			{
				real holeD = svar.Jet(0)+4*svar.Pstep; /*Diameter of hole (or width)*/

	            /*Find the points on the side of the pipe (Assume a circle)*/
	            float dtheta = atan((step)/(0.5*holeD));
	            Ny = ceil(svar.Jet(1)/step);
	            int holeWall = ceil(2*M_PI/dtheta)*Ny;

				/*Simulation Points*/
				int simCount = 0;
				real jetR = 0.5*(svar.Jet(0));
				
				/*Do the centerline of points*/
				for (real z = -jetR; z <= jetR; z+= svar.dx)
					simCount++;

				for (real x = svar.dx; x < jetR ; x+=svar.dx)
				{ /*Do the either side of the centerline*/
					for (real z = -jetR; z <= jetR; z+= svar.dx)
					{	/*If the point is inside the hole diameter, add it*/
						if(((x*x)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
			    			simCount += 2;
					}
				}

				/*Need to add the pn already present*/
				int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

				partCount = holeWall + simPts;
			}
			#else
			{
	            /*Find the points on the side of the pipe*/
	            Ny = ceil(svar.Jet(1)/step);
	            int holeWall = 2*Ny;

				/*Simulation Points*/
				int simCount = 0;
				real jetR = 0.5*(svar.Jet(0));
				for (real x = -jetR; x <= jetR; x+= svar.dx)
					simCount++;

				/*Need to add the pn already present*/
				int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

				partCount = holeWall + simPts;
			}
			#endif
		}
		else if (svar.Bcase == 4)
		{	
			#if SIMDIM == 3
				uint simCount = 0;
				real radius = 0.5*svar.Start(0);

				for (real y = 0; y <= radius; y+=svar.dx)
				{	
					real xradius = sqrt(radius*radius - y*y);

					simCount++;

					for (real z = svar.dx; z <= xradius; z+= svar.dx)
					{ /*Do the centerline of points*/
						simCount +=2;
					}

					for (real x = svar.dx; x <= xradius ; x+=svar.dx)
					{ /*Do the either side of the centerline*/
						simCount += 2;

						for (real z = svar.dx; z <= xradius; z+= svar.dx)
						{
							if(((x*x) + (z*z) + (y*y)) <= (radius*radius) )
				    		{   /*If the point is inside the hole diameter, add it*/
								simCount += 4;
							}
						}	
					}
				}
				
				for (real y = -svar.dx; y >= -radius; y-=svar.dx)
				{	
					real xradius = sqrt(radius*radius - y*y);

					simCount++;

					for (real z = svar.dx; z <= xradius; z+= svar.dx)
					{ /*Do the centerline of points*/
						simCount +=2;
					}

					for (real x = svar.dx; x <= xradius ; x+=svar.dx)
					{ /*Do the either side of the centerline*/
						simCount += 2;

						for (real z = svar.dx; z <= xradius; z+= svar.dx)
						{
							if(((x*x) + (z*z) + (y*y)) <= (radius*radius) )
				    		{   /*If the point is inside the hole diameter, add it*/
								simCount += 4;
							}
						}	
					}
				}

				partCount = simCount;
			#else
				uint simCount = 0;
				real radius = 0.5*svar.Start(0);
				for (real y = 0; y <= radius; y+=svar.dx)
				{	
					++simCount;
					for (real x = svar.dx; x <= radius ; x+=svar.dx)
					{ /*Do the either side of the centerline*/
							if(((x*x) + (y*y)) <= (radius*radius))
				    		{   /*If the point is inside the hole diameter, add it*/
								++simCount;
								++simCount;
							}	
					}
				}

				for (real y = -svar.dx; y >= -radius; y-=svar.dx)
				{	
					++simCount;
					for (real x = svar.dx; x <= radius ; x+=svar.dx)
					{ /*Do the either side of the centerline*/
							if(((x*x) + (y*y)) <= (radius*radius))
				    		{   /*If the point is inside the hole diameter, add it*/
								++simCount;
								++simCount;
							}	
					}
				}

				partCount = simCount;
			#endif
		}
		else if (svar.Bcase == 5) 
		{	/*Piston driven flow*/

			#if SIMDIM == 2
				/*Create the reservoir tank*/
				real tankW = svar.Start(0);
				real tankD = svar.Start(1);
				real stepb = (svar.Pstep*svar.Bstep);

				uint pisCnt = ceil((tankW + 8*svar.dx-4*svar.Pstep)/stepb);
				
				/*Create the reservoir tank*/
				uint wall = ceil((tankD+4*svar.dx+6*svar.Pstep)/stepb);
				
				/*Create the tapering section*/
				real theta = atan(svar.Jet(1)/(0.5*tankW-0.5*svar.Jet(0)));
				real stepy = stepb*sin(theta);
				uint taper = ceil((svar.Jet(1))/stepy);

				/*Create the exit bit.*/
				uint exit = ceil(svar.Jet(1)/stepb);

				/*Simulation Particles*/
				uint simCount = ceil(tankW/svar.dx)*ceil(tankD/svar.dx);	

				partCount = pisCnt + 2*(wall+taper+exit) + simCount;
			#endif
		}	

	return partCount;
}

void GetYcoef(AERO& avar, const FLUID& fvar, const real diam)
{
	// #if SIMDIM == 3
	// 	avar.L = rad * std::cbrt(3.0/(4.0*M_PI));
	// #endif
	// #if SIMDIM == 2
	// 	avar.L = rad/sqrt(M_PI);
	// #endif

	avar.L = diam / 2.0;
	
	avar.td = (2.0*fvar.rho0*pow(avar.L,SIMDIM-1))/(avar.Cd*fvar.mu);

	avar.omega = sqrt((avar.Ck*fvar.sig)/(fvar.rho0*pow(avar.L,SIMDIM))-1.0/pow(avar.td,2.0));

	avar.tmax = -2.0 *(atan(sqrt(pow(avar.td*avar.omega,2.0)+1)
					+avar.td*avar.omega) - M_PI)/avar.omega;

	avar.Cdef = 1.0 - exp(-avar.tmax/avar.td)*(cos(avar.omega*avar.tmax)+
		1/(avar.omega*avar.td)*sin(avar.omega*avar.tmax));
	avar.ycoef = 0.5*avar.Cdef*(avar.Cf/(avar.Ck*avar.Cb))*(avar.rhog*avar.L)/fvar.sig;

	//cout << avar.ycoef << "  " << avar.Cdef << "  " << avar.tmax << "  " << endl;
}

void Write_Input(SIM const& svar, FLUID const& fvar, AERO const& avar)
{
	// Open the file.
	string file = svar.outfolder;
	file.append("Settings");
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
	sett << svar.cellSize << setw(width) << "#Post processing mesh size" << endl;
	sett << svar.postRadius << setw(width) << "#Post processing support radius" << endl;
	sett << svar.Pstep << setw(width) << "#Particle initial spacing" << endl;
	sett << svar.Bstep << setw(width) << "#Boundary spacing factor" << endl;
	sett << svar.Bcase << setw(width) << "#Boundary case" << endl;
	sett << svar.Asource << setw(width) << "#Aerodynamic source" << endl; 
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
		sett << svar.Angle(0) << " " << svar.Angle(1);
#if SIMDIM == 3
		sett << " " << svar.Angle(2);
#endif
		sett << "  #Fluid start rotation" << endl;
		sett << svar.Jet(0) << " " << svar.Jet(1) << setw(width) << "#Jet dimensions" << endl;
		sett << fvar.pPress << setw(width) << "#Pipe Pressure" << endl;
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
	file.append("Fluid");

	std::ofstream fluid(file);

	fluid << svar.beta << " " << svar.gamma << setw(width) << "#Newmark-Beta terms" << endl;
	fluid << fvar.Hfac << setw(width) << "#Support radius factor" << endl;
	fluid << fvar.alpha << setw(width) << "#Artificial viscosity" << endl;
	fluid << fvar.contangb << setw(width) << "#Contact angle" << endl;
	fluid << fvar.rho0 << setw(width) << "#Fluid density" << endl;
	fluid << avar.rhog << setw(width) << "#Gas density" << endl;
	fluid << fvar.Cs << setw(width) << "#Speed of sound" << endl;
	fluid << fvar.mu << setw(width) << "#Fluid viscosity" << endl;
	fluid << avar.mug << setw(width) << "#Gas viscosity" << endl;
	fluid << fvar.sig << setw(width) << "#Surface Tension" << endl;
	fluid << svar.outdir << setw(width) << "#Output Folder" << endl;
	fluid << svar.meshfile << "   #Mesh edge/face file" << endl;
	fluid << svar.solfile << "   #Mesh Solution file" << endl;
	fluid << svar.scale << setw(width) << "#Mesh scale" << endl;
	fluid << avar.vRef << setw(width) << "#Gas reference velocity" << endl;
	fluid << avar.pRef << setw(width) << "#Gas reference pressure" << endl;
	fluid << avar.T << setw(width) << "#Gas reference temperature" << endl;
	
	fluid.close();
}


#endif