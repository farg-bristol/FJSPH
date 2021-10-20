
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








#endif