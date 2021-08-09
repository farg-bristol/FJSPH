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

const std::string WHITESPACE = " \n\r\t\f\v";

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
void Check_If_Restart_Possible(SIM const& svar)
{
  	/*Check output folder for any prexisting files*/
  	if(svar.out_encoding == 1)
  	{
  		string file = svar.restart_prefix;
  		file.append("_fuel.szplt.szdat");
		#ifdef DEBUG
			dbout << "Checking for existence of previous szplt files." << endl;
			dbout << "Path: " << file << endl;
		#endif
  		struct stat info;
  		if(stat( file.c_str(), &info ) == 0)
  		{
  			// cout << "Split files exist..." << endl;

			// Check if there is the constructed szplt
			file = svar.restart_prefix;				
			file.append("_fuel.szplt");
			if(stat( file.c_str(), &info ) != 0)
			{	/*Fuel.szplt does not exist - needs creating*/
				if(Combine_SZPLT(file) == -1)
					exit(-1);
			}
			else
			{	// Both the combined and split szplt exist - check which is more recent.
				// cout << "Both file types exist" << endl << endl;
				auto sztime = info.st_mtime;
				file.append(".szdat");
				if(stat(file.c_str(), &info) == 0)
				{
					auto szdattime = info.st_mtime;
					if (szdattime > (sztime+10))
					{	// combine the files
						file = svar.output_prefix;
						file.append("_fuel.szplt");
						if(Combine_SZPLT(file) == -1)
							exit(-1);
					}
				}
			}
		}
		else
		{
			cout << "The fuel szplt files needed to restart could not be found. Stopping." << endl;
			exit(-1);
		}

		if(svar.Bcase != 0)
		{
			file = svar.restart_prefix;
			file.append("_boundary.szplt.szdat");

			if(stat( file.c_str(), &info ) == 0)
			{
				// cout << "Split files exist..." << endl;
				// Check if there is the constructed szplt
				file = svar.restart_prefix;				
				file.append("_boundary.szplt");
				if(stat( file.c_str(), &info ) != 0)
				{	/*Fuel.szplt does not exist - needs creating*/
					if(Combine_SZPLT(file) == -1)
						exit(-1);
				}
				else
				{	// Both the combined and split szplt exist - check which is more recent.
					// cout << "Both file types exist" << endl << endl;
					auto sztime = info.st_mtime;
					file.append(".szdat");
					if(stat(file.c_str(), &info) == 0)
					{
						auto szdattime = info.st_mtime;
						if (szdattime > (sztime+10))
						{	// combine the files
							file = svar.restart_prefix;
							file.append("_boundary.szplt");
							if(Combine_SZPLT(file) == -1)
								exit(-1);
						}
					}
				}
			}
			else
			{
				cout << "The boundary szplt files needed to restart could not be found. Stopping." << endl;
				exit(-1);
			}			
		}
  	}
}

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string Get_Parameter_Value(string const& line)
{
    size_t pos = line.find(":");
    size_t end = line.find("#",pos+1); /* Check if a comment exists on the line */

    if (end != string::npos)
    {
        string value = line.substr(pos + 1, (end-pos+2) );
        return ltrim(rtrim(value));
    }

    string value = line.substr(pos + 1);
    return ltrim(rtrim(value));
}

void Get_String(string const& line, string const& param, string &value)
{
    if(line.find(param) != string::npos)
    {
        value = Get_Parameter_Value(line);
    }
}

template<typename T>
void Get_Number(string const& line, string const& param, T &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        iss >> value;
    }
}

void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<real,3,1>/* vec<real,3> */ &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        
        real a, b, c;
        string temp2;
        
        std::getline(iss,temp2,',');
        std::istringstream iss2(temp2);
        iss2 >> a;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> b;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> c;
        
        value = /* vec<real,3> */ Eigen::Matrix<real,3,1>(a,b,c);
    }
}

void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<real,2,1>/* vec<real,2> */ &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        
        real a, b;
        string temp2;
        
        std::getline(iss,temp2,',');
        std::istringstream iss2(temp2);
        iss2 >> a;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> b;
        
        value = /* vec<real,2> */ Eigen::Matrix<real,2,1>(a,b);
    }
}

void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<int,3,1>/* vec<int,3> */ &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        
        int a, b, c;
        string temp2;
        
        std::getline(iss,temp2,',');
        std::istringstream iss2(temp2);
        iss2 >> a;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> b;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> c;
        
        value = /* vec<int,3> */ Eigen::Matrix<int,3,1>(a,b,c);
    }
}

void Get_Vector(string const& line, string const& param, 
				Eigen::Matrix<int,2,1>/* vec<int,2> */ &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        
        int a, b;
        string temp2;
        
        std::getline(iss,temp2,',');
        std::istringstream iss2(temp2);
        iss2 >> a;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> b;
        
        value = Eigen::Matrix<int,2,1>(a,b) /* vec<int,2>(a,b) */;
    }
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

void Restart(SIM& svar, FLUID const& fvar, MESH const& cells, State& pn, State& pnp1)
{	
	// Read the values from the solution folder, then check. 

	// Check that a restart can be performed. 
	if(svar.outform == 0 || svar.outform == 1 || svar.outform == 4)
	{
		cout << "Output type cannot be restarted from. Make sure that you have the correct output selected." << endl;
		exit(-1);
	}

	if(svar.Asource != 0)
	{
		if(svar.outform != 3 || svar.outform != 5)
		{
			cout << "No cell information. Cannot restart" << endl;
			exit(-1);
		}
	}

	#ifdef DEBUG
		dbout << "Reading frame info file for restart information." << endl;
	#endif

	// Re-read the settings and fluid file
	// cout << svar.outfolder << endl;

	// Now get the data from the files. Start with the boundary
	if(svar.out_encoding == 1)
	{
		State boundary, fuel, rest;	

		// Read the fuel
		void* fuelHandle = NULL;
		void* boundHandle = NULL;

		int32_t fuelFrames, boundFrames;
		double fuelTime, boundTime;

		string fuelf = svar.output_prefix;
		fuelf.append("_fuel.szplt");

		if(tecFileReaderOpen(fuelf.c_str(),&fuelHandle))
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << fuelf << endl;
			exit(-1);
		}

		// Check how many frames are in the fuel file.
		cout << "Checking Fuel file..." << endl;
		CheckContents(fuelHandle,svar,fuelFrames,fuelTime);

		if (svar.Bcase !=0)
		{
			string boundf = svar.output_prefix;
			boundf.append("_boundary.szplt");

			if(tecFileReaderOpen(boundf.c_str(),&boundHandle))
			{
				cout << "Error opening szplt file. Path:" << endl;
				cout << boundf << endl;
				exit(-1);
			}

			cout << "Checking Boundary file..." << endl;
			CheckContents(boundHandle,svar,boundFrames,boundTime);

			if(fuelFrames!= boundFrames)
			{
				cout << "Caution! Number of frames is not consistent between fuel and boundary files." << endl;
			}

			if(fuelTime != boundTime)
			{
				cout << "Caution! Frame times are not consistent between fuel and boundary files." << endl;

				if(fuelTime > boundTime)
				{
					double time = 0.0;
					for(int32_t frame = fuelFrames-1; frame > 1; frame--)
					{
						if(tecZoneGetSolutionTime(fuelHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from fuel file." << endl;
							continue;
						}
						
						if(time == boundTime)
						{
							cout << "Found the correct frame" << endl;
							fuelFrames = frame;
							fuelTime = time;
							break;
						}
					}
				}
				else
				{
					double time = 0.0;
					for(int32_t frame = boundFrames-1; frame >= 1; frame--)
					{
						if(tecZoneGetSolutionTime(boundHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from boundary file." << endl;
							continue;
						}
						

						if(time == fuelTime)
						{
							cout << "Found the correct frame" << endl;
							boundFrames = frame;
							boundTime = time;
							break;
						}
					}
				}

				if(fuelTime != boundTime)
				{
					cout << "Could not find a consistent time in each file. Stopping." << endl;
					exit(-1);
				}
			}

			cout << "Attempting to read the boundary..." << endl;
			Read_Binary_Timestep(boundHandle,svar,boundFrames,boundary);
		}

	    svar.frame = fuelFrames-1;

	    // Read the actual data.
		cout  << "Attempting to read the fuel..." << endl;
		Read_Binary_Timestep(fuelHandle,svar,fuelFrames,fuel);


		pn = boundary;
		svar.bndPts = boundary.size();
		svar.simPts = fuel.size();
		pn.insert(pn.end(),fuel.begin(),fuel.end());	
		svar.totPts = pn.size();
		
		if(svar.simPts + svar.bndPts != svar.totPts)
		{
			cout << "Mismatch of array sizes. Total array is not the sum of sim and boundary arrays" << endl;
			exit(-1);
		}

		// Check the frame timings to ensure dt does not go negative.
		double frametime = svar.framet * real(fuelFrames-1);

		if(svar.t > frametime)
		{
			// There's a problem with the frame times.
			cout << "Time is further ahead than what is supposed by the frame timings." << endl;
			cout << "Need to adjust frame time to match the current time" << endl;

			real framet = svar.t/real(fuelFrames-1.0);
			cout << "Old frame time: " << svar.framet << "  New frame time: " << framet << endl << endl;
			svar.framet = framet;			
		}

		// if(svar.totPts != totPts || svar.bndPts != bndPts || svar.simPts != simPts)
		// {
		// 	cout << "Mismatch of the particle numbers in the frame file and data file." << endl;
		// }
	}
	else
	{
		// TODO: ASCII Restart.
		// Particle numbers can be found from frame file.
		// Find EOF, then walk back from there how many particles.

		cout << "ASCII file restart not yet implemented." << endl;
		exit(-1);
	}

	// Go through the particles giving them the properties of the cell
	vector<size_t> buffer;
	#pragma omp parallel for 
	for(size_t ii = 0; ii < svar.totPts; ++ii)
	{
		pn[ii].partID = ii;
		pn[ii].p =  fvar.B*(pow(pn[ii].rho/fvar.rho0,fvar.gam)-1);

		if(pn[ii].b == PartState.BACK_)
		{
			#pragma omp critical
			svar.back.emplace_back(ii);
		}
		else if(pn[ii].b == PartState.BUFFER_)
		{
			#pragma omp critical
			buffer.emplace_back(ii);
		}
		
		// Initialise the rest of the values to 0
		pn[ii].s = 0.0;
		pn[ii].woccl = 0.0;
		pn[ii].pDist = 0.0;
		pn[ii].internal = 0;

		pn[ii].vPert = StateVecD::Zero();
	}	

	/* Put the particles into the buffer vector */
	svar.buffer = vector<vector<size_t>>(svar.back.size(),vector<size_t>(4));
	real eps = 0.01*svar.dx; /* Tolerance value */
	for(size_t ii = 0; ii < svar.back.size(); ++ii)
	{
		StateVecD const test = svar.Transp*(pn[svar.back[ii]].xi-svar.sim_start);
		
		for(size_t index = 0; index < buffer.size(); ++index)
		{
			/* Check which particle in the back vector it corresponds to by checking which  */
			/* particle it lies behind */
			StateVecD const xi = svar.Transp*(pn[buffer[index]].xi-svar.sim_start);
			// cout << test[0] << "  " << test[1] << "  " << xi[0] << "  " << xi[1] << endl;

			if(xi[0] < test[0] + eps && xi[0] > test[0] - eps)
			{	/* X coordinate is within bounds, so should lie behind this point */

				/* Start wit the furthest away, and go closer */
				if (xi[1] < test[1] - 4.0*svar.dx + eps)
				{
					svar.buffer[ii][3] = buffer[index];
				}
				else if (xi[1] < test[1] - 3.0*svar.dx + eps)
				{
					svar.buffer[ii][2] = buffer[index];
				}
				else if (xi[1] < test[1] - 2.0*svar.dx + eps)
				{
					svar.buffer[ii][1] = buffer[index];
				}
				else if (xi[1] < test[1] - svar.dx + eps)
				{
					svar.buffer[ii][0] = buffer[index];
				}
				else
				{
					cout << "Couldn't identify where to place the buffer particle" << endl;
					exit(-1);
				}
			}
		}
	}	
	
	pnp1 = pn;
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