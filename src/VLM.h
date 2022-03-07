/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VLM_H
#define VLM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"
#include "Third_Party/Eigen/LU"
#include "Third_Party/Eigen/Geometry"
#include "Var.h"

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif


/*A structure for the coordinates to define a panel.*/
typedef struct Panel
{ /*A and B are 1/4 chord bounds for the vortex.
	C is the control Point location*/
	Panel(StateVecD A, StateVecD B, StateVecD C,
			StateVecD p1, StateVecD p2, StateVecD p3, StateVecD p4)
	: A(A), B(B), C(C), p1(p1), p2(p2), p3(p3), p4(p4) 
	{norm =  (C-A).cross(C-B).normalized();}

	StateVecD A, B, C;
	StateVecD p1, p2, p3, p4;
	StateVecD norm;
}Panel;

typedef class VLM
{
	public:
		VLM()
		{
			coords = Eigen::Matrix<real,2,1>(-1,-1);
			panels = Eigen::Matrix<int,2,1>(-1,-1);
			AoA = 0;
			sweep = 0;
			taper = 0;
			flap = Eigen::Vector3i(0,0,0);
			beta = 0;

			/*Want freestream to be aligned with jet axis*/
			Freestream = Eigen::Matrix<real,3,1>(0,1,0);
		}

		void Init(std::string input)
		{
			// input.append("VLM.dat");
			std::ifstream fin(input, std::ios::in);

			if(fin.is_open())
			{
				string line;
				while (getline(fin,line))
				{
					Get_Vector(line,"VLM wing dimensions", coords);
					Get_Vector(line,"VLM panel counts",panels);
					Get_Number(line,"VLM angle alpha (degree)",AoA);
					Get_Number(line,"VLM sweep angle (degree)",sweep);
					Get_Number(line,"VLM taper ratio",taper);
					Get_Number(line,"VLM flap start (panels)",flap[0]);
					Get_Number(line,"VLM flap width (panels)",flap[1]);
					Get_Number(line,"VLM flap depth (panels)",flap[2]);
					Get_Number(line,"VLM flap angle (degree)",beta);
				}
			}
			else 
			{
				std::cerr << "Couldn't open " << input << " to read VLM settings. Stopping." << std::endl;
				exit(-1);
			}
			fin.close();

			if(coords[0] == -1 || coords[1] == -1)
			{
				cout << "VLM coordinates not defined." << endl;
				exit(-1);
			}

			if(panels[0] == -1 || panels[1] == -1)
			{
				cout << "VLM panels numbers not defined." << endl;
				exit(-1);
			}

			npanels = 2*panels[0]*panels[1];

			AoA *= M_PI/180.0;
			sweep *= M_PI/180.0;
			beta *= M_PI/180.0;

			// coords[0] = 5;
			// coords[1] = 2;
			// panels(0) = 10;
			// panels(1) = 5;
			// AoA = 12.0 * M_PI/180.0;
			// sweep = 30 * M_PI/180.0;
			
			/*Initialise matrices*/
			panelData.reserve(npanels);
			aInf = Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic>(npanels,npanels);
			gamma = Eigen::Matrix<real,Eigen::Dynamic,1>(npanels);
			RHS = Eigen::Matrix<real,Eigen::Dynamic,1>(npanels);

			MakeMatrix();
		}

		void GetGamma(StateVecD inf)
		{
			Freestream = inf;
			/*Find the influence matrix, and invert it to find gamma*/

			/*Find the influence matrix aInf*/
 			for(int i=0; i < npanels; ++i)
 			{
 				for(int j=0; j < npanels; ++j)
 				{
 					/*aInf[i,j] = the influence of vortex j on control point i*/
 					StateVecD A, B, C;
 					A = panelData[i].A;
 					B = panelData[i].B;
 					C = panelData[j].C;

					aInf(j,i) = FindInfluence(A,B,C).dot(panelData[i].norm);
				}

			RHS(i) = -Freestream.dot(panelData[i].norm);
			}
			/*Find Gamma */
			Eigen::FullPivLU<Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic>> lu(aInf);
			if (lu.isInvertible())
 				gamma = lu.solve(RHS);
 			else
 				std::cerr << "Influence matrix is singular" << std::endl;

 			// for (int i= 0; i < gamma.rows(); ++i)
 			// {
 			// 	std::cout << gamma[i] << std::endl;
 			// }
 		}

		const StateVecD getVelocity(const StateVecD pos)
		{	/*Find velocity for a particle at its position*/
			StateVecD vel = StateVecD::Zero();

			for(int i = 0; i<npanels; ++i)
				vel += gamma(i)*FindInfluence(panelData[i].A,panelData[i].B,pos);
			
			return vel+Freestream;
		}

		void write_VLM_Panels(string &prefix)
		{	
			string file1 = prefix;
			file1.append("_Panels.dat");
			ofstream fp(file1, std::ios::out);
			if(fp.is_open())
			{
				fp << "TITLE=\"VLM Panels\"" << std::endl;
		  		fp << "VARIABLES = \"X\", \"Y\", \"Z\""<< std::endl;
		  		for (auto p:panelData)
				{
					fp << "ZONE" << std::endl;
					// fp << "VARLOCATION=([1-3]=NODAL)" << endl;
			        fp << p.p1(0) << " " << p.p1(1) << " " << p.p1(2) << std::endl;
			        fp << p.p2(0) << " " << p.p2(1) << " " << p.p2(2) << std::endl;
			        fp << p.p3(0) << " " << p.p3(1) << " " << p.p3(2) << std::endl;
			        fp << p.p4(0) << " " << p.p4(1) << " " << p.p4(2) << std::endl;
			        fp << p.p1(0) << " " << p.p1(1) << " " << p.p1(2) << std::endl;  
			  	}
			 	fp.close();
			}
			else
			{
				std::cerr << "Failed to open VLM_Panels.plt. Attempted path:" << std::endl;
				std::cerr << file1 << std::endl;
				exit(-1);
			}
	  	
		  	string file2 = prefix;
		  	file2.append("_Vortices.dat");
		  	std::ofstream fq(file2, std::ios::out);
		  	if(fq.is_open())
			{
		  		fq << "TITLE=\"VLM Vortices and Control Points\"" << std::endl;
		  		fq << "VARIABLES = \"X\", \"Y\", \"Z\""<< std::endl;
		  		for (auto p:panelData)
				{
					fq << "ZONE" << std::endl;
					fq << p.A(0) << " " << p.A(1) << " " << p.A(2) << std::endl;
					fq << p.B(0) << " " << p.B(1) << " " << p.B(2) << std::endl;
					fq << "ZONE" << std::endl;
					fq << p.C(0) << " " << p.C(1) << " " << p.C(2) << std::endl;
				}
				fq.close();
			}
			else
			{
				std::cerr << "Failed to open VLM_Vortices.plt. Attempted path:" << std::endl;
				std::cerr << file2 << std::endl;
				exit(-1);
			}
			
		}

	protected:
		
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
			size_t end = line.find_first_of("#",pos+1); /* Check if a comment exists on the line */

			if (end != string::npos)
			{
				string value = line.substr(pos + 1, (end-pos+2) );
				return ltrim(rtrim(value));
			}

			string value = line.substr(pos + 1);
			return ltrim(rtrim(value));
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

		void MakeMatrix(void)
		{
			/*Define the steps for each dimension*/
			real dr0 = coords[1]/real(panels[1]-1);
			real dx0 = (coords[0]/real(panels[0]-1))*cos(sweep);
			real dy0 = cos(AoA)*dr0;
			real dz0 = sin(AoA)*dr0;
			
			real x0 = 0.0;
			real y0 = cos(AoA)*coords[1];
			real z0 = sin(AoA)*coords[1];
			
			real xend = coords[0]*cos(sweep);
			
			/*Find 1/4 panel points (assuming symmetry) */
			/*Find 3/4 Control Points*/
			for(int i = -panels[0]; i < panels[0]; ++i)
			{	/*i = count in z, (spanwise)*/
				/*X coordinates*/
				real x1 = real(i)*dx0;
				real x2 = real(i+1)*dx0;
				real xhalf = (real(i)+0.5)*dx0;
				
				/*x and y values at j = 0 to apply sweep*/
				real yi = tan(sweep)*(fabs(x1)-x0);
				real zi = sin(AoA)*sin(sweep)*(fabs(x1)-x0);
				real yip1 = tan(sweep)*(fabs(x2)-x0);
				real zip1 = sin(AoA)*sin(sweep)*(fabs(x2)-x0);

				/*Vertex 1 deltas (i)*/
				real frac1 = (fabs(x1)-x0)/(xend-x0);
				real dy1 =  dy0*(1-frac1) + dy0*taper*frac1;
				real dz1 =  dz0*(1-frac1) + dz0*taper*frac1;
				
				/*Vertex 2 deltas (i+1)*/
				real frac2 = (fabs(x2)-x0)/(xend-x0);
				real dy2 =  dy0*(1-frac2) + dy0*taper*frac2;
				real dz2 =  dz0*(1-frac2) + dz0*taper*frac2;

				/*Halfway values (i+1/2) for control point*/
				real frac3 = (fabs(xhalf)-x0)/(xend-x0);
				real yhalf = tan(sweep)*(fabs(xhalf)-x0);
				real zhalf = sin(AoA)*sin(sweep)*(fabs(xhalf)-x0);
				real dzhalf = dz0*(1-frac3) + dz0*taper*frac3;
				real dyhalf =  dy0*(1-frac3) + dy0*taper*frac3;

				int jend;
				int doflap = 0; /*0 = no flap, 1 = flap*/

				
				if( -i > flap(0) && -i <= flap(1)+flap(0))
				{	/*Left wing*/
					jend = panels(1)-flap(2);
					doflap = 1;
				}
				else if ( i >= flap(0) && i < flap(1)+flap(0))
				{	/*Right wing*/
					jend = panels(1)-flap(2);
					doflap = 1;
				}
				else
				{	/*No flap*/
					jend = panels(1);
				}

				for (int j=0; j < jend; ++j)
				{	/*j = count in x and y, (chordwise)*/
					/*Panel Verticies (for visual)*/
					StateVecD p1(x1, dy1*(real(j))-y0+yi    , z0-zi-dz1*(real(j)));
					StateVecD p2(x1, dy1*(real(j+1))-y0+yi  , z0-zi-dz1*(real(j+1)));
					StateVecD p3(x2, dy2*(real(j+1))-y0+yip1, z0-zip1-dz2*(real(j+1)));
					StateVecD p4(x2, dy2*(real(j))-y0+yip1  , z0-zip1-dz2*(real(j)));

					/*1/4 chord points*/
					StateVecD A(x1, dy1*(0.25+real(j))-y0+yi  , z0-zi-dz1*(0.25+real(j)));
					StateVecD B(x2, dy2*(0.25+real(j))-y0+yip1, z0-zip1-dz2*(0.25+real(j)));
					/*3/4 control point*/
					StateVecD C(xhalf, dyhalf*(0.75+real(j))-y0+yhalf,z0-zhalf-dzhalf*(0.75+real(j)));
					panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
				}

				if(doflap == 1)
				{
					/*Flap delta parameters*/
					real dzflap1 = (sin(beta+AoA))*dr0*(1-frac1) + (sin(beta+AoA))*dr0*taper*frac1;
					real dyflap1 = (cos(beta+AoA))*dr0*(1-frac1) + (cos(beta+AoA))*dr0*taper*frac1;
					real dzflap2 = (sin(beta+AoA))*dr0*(1-frac2) + (sin(beta+AoA))*dr0*taper*frac2;
					real dyflap2 = (cos(beta+AoA))*dr0*(1-frac2) + (cos(beta+AoA))*dr0*taper*frac2;
					real dzflaphalf = (sin(beta+AoA))*dr0*(1-frac3) + (sin(beta+AoA))*dr0*taper*frac3;
					real dyflaphalf = (cos(beta+AoA))*dr0*(1-frac3) + (cos(beta+AoA))*dr0*taper*frac3;

					/*Flap start position*/
					real zflap1 = z0-zi-real(panels(1)-flap(2))*dz1;
					real yflap1 = real(panels(1)-flap(2))*dy1 + yi - y0;
					real zflap2 = z0 - zip1-real(panels(1)-flap(2))*dz2;
					real yflap2 = real(panels(1)-flap(2))*dy2 + yip1 - y0;
					real zflaphalf = z0-zhalf-real(panels(1)-flap(2))*dzhalf;
					real yflaphalf = real(panels(1)-flap(2))*dyhalf + yhalf - y0;

					for (int j = 0; j < flap(2); ++j)
					{
						/*j = count in x and y, (chordwise)*/
						/*Panel Verticies (for visual)*/
						StateVecD p1(x1, yflap1+dyflap1*(real(j))  , zflap1-dzflap1*(real(j)));
						StateVecD p2(x1, yflap1+dyflap1*(real(j+1)), zflap1-dzflap1*(real(j+1)));
						StateVecD p3(x2, yflap2+dyflap2*(real(j+1)), zflap2-dzflap2*(real(j+1)));
						StateVecD p4(x2, yflap2+dyflap2*(real(j))  , zflap2-dzflap2*(real(j)));

						/*1/4 chord points*/
						StateVecD A(x1, yflap1+dyflap1*(0.25+real(j)), zflap1-dzflap1*(0.25+real(j)));
						StateVecD B(x2, yflap2+dyflap2*(0.25+real(j)), zflap2-dzflap2*(0.25+real(j)));
						/*3/4 control point*/
						StateVecD C(xhalf, yflaphalf+dyflaphalf*(0.75+real(j)),zflaphalf-dzflaphalf*(0.75+real(j)));
						panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
					}
				}
			}
		
			/*End initialisation*/
		}

		const StateVecD FindInfluence(StateVecD const& A, StateVecD const& B, StateVecD const& C)
 		{	
			StateVecD r0, r1, r2, inf;

			/*Bounded Vortex*/
			StateVecD coefAB;

			r0 = B-A;
			r1 = C-A;
			r2 = C-B;

			coefAB = (1/(4*M_PI))*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(r0.dot(r1.normalized()-r2.normalized()));


			/*Horseshoe vortex from point A*/
			StateVecD coefA;
		
			inf = A + Freestream;
			r2 = C - A;
			r1 = C - inf;
			r0 = A - inf;

			coefA = (1/(4*M_PI))*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(1-r0.dot(r2)/(r0.norm()*r2.norm()));

			/*Horseshoe vortex from point B*/
			StateVecD coefB;

			/*Vector B to infinity*/
			inf = B + Freestream;
			r1 = C - B;
			r2 = C - inf;
			r0 = inf - B;

			coefB = (1/(4*M_PI))*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(r0.dot(r1)/(r1.norm()*r0.norm())+1.0);

			return (coefAB+coefB+coefA);			
 		}

 		std::vector<Panel> panelData;

		int npanels, nverts;
		real AoA, sweep, taper, beta;
		
		Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic> aInf;

		Eigen::Matrix<real,Eigen::Dynamic,1> gamma;
		Eigen::Matrix<real,Eigen::Dynamic,1> RHS;

		Eigen::Matrix<real,2,1> coords;
		Eigen::Vector2i panels;
		
		Eigen::Vector3i flap;

		StateVecD Freestream;
}VLM;

#endif