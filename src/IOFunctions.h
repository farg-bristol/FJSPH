
#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H
	

#include "Var.h"
#include <filesystem>

inline void write_header() 
{
	cout << "******************************************************************" << endl << endl;
	cout << "                              WCSPH                               " << endl << endl;
	cout << "        Weakly Compressible Smoothed Particle Hydrodynamics       " << endl;
	cout << "                      for Fuel Jettison case                      " << endl << endl;
	cout << "                         James O. MacLeod                         " << endl;
	cout << "                    University of Bristol, U.K.                   " << endl << endl;
	cout << "******************************************************************" << endl << endl;
}


inline std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
inline std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

inline string Get_Parameter_Value(string const& line)
{
    size_t pos = line.find(":");
    size_t end = line.find_first_of("#",pos+1); /* Check if a comment exists on the line */

    if (end != string::npos)
    {
        string value = line.substr(pos + 1, (end-pos+2) );
        return ltrim(rtrim(value));
    }

    string value = line.substr(pos + 1);
    return ltrim(rtrim(value));
}

inline void Get_String(string const& line, string const& param, string &value)
{
    size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
			value = Get_Parameter_Value(line);
	}
}

inline void Get_Number(string const& line, string const& param, int &value)
{
	size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			iss >> value;
		}
	}
}

inline void Get_Number(string const& line, string const& param, real &value)
{
	size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			iss >> value;
		}
	}
}

inline void Get_Number(string const& line, string const& param, size_t &value)
{
	size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			iss >> value;
		}
	}
}

inline void Get_Number(string const& line, string const& param, uint &value)
{
	size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			iss >> value;
		}
	}
}

inline void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<real,2,1> &value)
{
    size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
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
			
			value = Eigen::Matrix<real,2,1>(a,b);
		}
	}
	else
		return;

}

inline void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<int,2,1> &value)
{
    size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
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
			
			value = Eigen::Matrix<int,2,1>(a,b);
		}
	}
}

inline void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<real,3,1> &value)
{
    size_t pos = line.find(":");
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
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
			
			value = Eigen::Matrix<real,3,1>(a,b,c);
		}
	}
}

inline void Get_Vector(string const& line, string const& param, 
			Eigen::Matrix<int,3,1> &value)
{
    size_t pos = line.find(":");
	string substr;
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));	
		if(substr == param)
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
			
			value = Eigen::Matrix<int,3,1>(a,b,c);
		}
	}
}

inline void Get_Array(string const& line, string const& param, 
				vector<real>& value)
{
    size_t pos = line.find(":");
	string substr;
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			
			real a;
			string temp2;
			
			while(std::getline(iss,temp2,','))
			{
				std::istringstream iss2(temp2);
				iss2 >> a;
				value.emplace_back(a);
			} 
		}
	}
}

inline void Get_Array(string const& line, string const& param, 
				vector<int>& value)
{
    size_t pos = line.find(":");
	string substr;
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			
			int a;
			string temp2;
			
			while(std::getline(iss,temp2,','))
			{
				std::istringstream iss2(temp2);
				iss2 >> a;
				value.emplace_back(a);
			} 
		}
	}
}

inline void Get_Array(string const& line, string const& param, 
				vector<size_t>& value)
{
    size_t pos = line.find(":");
	string substr;
	if(pos != string::npos)
	{
		string substr = ltrim(line.substr(0,pos));
		if(substr == param)
		{
			string temp = Get_Parameter_Value(line);
			std::istringstream iss(temp);
			
			size_t a;
			string temp2;
			
			while(std::getline(iss,temp2,','))
			{
				std::istringstream iss2(temp2);
				iss2 >> a;
				value.emplace_back(a);
			} 
		}
	}
}

inline size_t index(size_t ii, size_t jj, size_t nii)
{
	return(ii + jj*nii);
}

inline size_t index(size_t ii, size_t jj, size_t kk, size_t ni, size_t nj)
{
	return  (kk * nj + jj) * ni + ii;
}

inline std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(uint ii=0; ii < num - 1; ++ii){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

#endif