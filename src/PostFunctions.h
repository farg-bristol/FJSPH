#ifndef POSTFUNCTIONS_H
#define POSTFUNCTIONS_H

	
#include "PostVar.h"

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
	cout << "*******************************************************************" << endl << endl;
	cout << "                       WCSPH Post Processing                       " << endl << endl;
	cout << "        Weakly Compressible Smoothed Particle Hydrodynamics        " << endl << endl;
	cout << "                         James O. MacLeod                          " << endl;
	cout << "                    University of Bristol, U.K.                    " << endl << endl;
	cout << "*******************************************************************" << endl << endl;
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
  		file.append("Fuel.szplt");
#ifdef DEBUG
  		dbout << "Checking for existence of previous szplt files..." << endl;
  		dbout << "Path: " << file << endl;
#endif
  		struct stat info;
  		if(stat( file.c_str(), &info ) == 0)
  		{
#ifdef DEBUG
  			dbout << "Exists..." << endl;
#endif
		}
		else
		{
  			// File needs creating.
#ifdef DEBUG
  			dbout << "Does not exist..." << endl;
#endif
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

// StateVecI getIVector(ifstream& In, uint& lineno, const string& name)
// {
// 	string line;
// 	getline(In,line);
// 	lineno++;
// 	std::stringstream sline(line);

// 	vector<int> vec;
// 	int value;
// 	string temp;
// 	for(uint ii = 0; ii < 3; ii++)
// 	{
// 		sline >> temp;
// 		if(std::stringstream(temp)>> value)
// 		{
// 			vec.emplace_back(value);
// 		}
// 		temp = "";
// 	}

// 	if(vec.size() == 0)
// 	{
// 		cout << "Error: Dimension is zero. Check the input order." << endl;
// 		cout << "Line " << lineno << ": " << line << endl;  
// 		exit(-1);
// 	}

// 	DIM = vec.size();

// 	StateVecI x(DIM);

// 	for(size_t ii = 0; ii < DIM; ii++)
// 	{
// 		x(ii) = vec[ii];
// 	}

// #ifdef DEBUG
// 	dbout << name << ": " << x(0) << "  " << x(1);
// 	if(DIM == 3)
// 		dbout << "  " << x(2);
// 	dbout << endl;
// #endif	

// 	return x;
// }

// StateVecD getDVector(ifstream& In, uint& lineno, const string& name)
// {
// 	string line;
// 	getline(In,line);
// 	lineno++;
// 	std::stringstream sline(line);
	
// 	vector<real> vec;
// 	real value;
// 	string temp;
// 	for(uint ii = 0; ii < 3; ii++)
// 	{
// 		sline >> temp;
// 		if(std::stringstream(temp)>> value)
// 		{
// 			vec.emplace_back(value);
// 		}
// 		temp = "";
// 	}

// 	if(vec.size() == 0)
// 	{
// 		cout << "Error: Dimension is zero. Check the input order." << endl;
// 		cout << "Line " << lineno << ": " << line << endl; 
// 		exit(-1);
// 	}

// 	DIM = vec.size();

// 	StateVecD x(DIM);

// 	for(size_t ii = 0; ii < DIM; ii++)
// 	{
// 		x(ii) = vec[ii];
// 	}

// #ifdef DEBUG
// 	dbout << name << ": " << x(0) << "  " << x(1);
// 	if(DIM == 3)
// 		dbout << "  " << x(2);
// 	dbout << endl;
// #endif	

// 	return x;
// }

/*Function for a 2D Vector (e.g. Newmark Beta parameters)*/
// Eigen::Vector2d getvector(ifstream& In, uint& lineno, const string& name)
// {
// 	string line;
// 	getline(In,line);
// 	lineno++;
// 	std::istringstream sline(line);
	
// 	Eigen::Vector2d x;
// 	sline >> x[0]; sline >> x[1]; 
// #ifdef DEBUG
// 	dbout << name << ": " << x(0) << "  " << x(1) << endl;
// #endif	
// 	return x;
// }

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

///******Wendland's C2 Quintic Kernel*******
inline real W2Kernel(const real dist, const real H, const real correc) 
{
	const real q = dist/H;
	return (pow(1-0.5*q,4))*(2*q+1)*correc;
}

/*Gradient*/
inline StateVecD W2GradK(const StateVecD& Rij, const real dist, const real H, const real correc)
{
	const real q = dist/H;
	return 5.0*(Rij/(H*H))*pow(1-0.5*q,3)*correc;
}

/*2nd Gradient*/
inline real W2Grad2(const StateVecD& Rij, const real dist, const real H, const real correc) 
{
	const real q = dist/H;
	return Rij.dot(Rij)*(5.0*correc/(H*H))*(2*q-1)*pow(1-0.5*q,2);
}

#endif