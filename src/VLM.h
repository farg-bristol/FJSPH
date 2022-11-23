/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VLM_H
#define VLM_H

#include "Var.h"
#include "Third_Party/Eigen/Core"
#include "Third_Party/Eigen/StdVector"
#include "Third_Party/Eigen/LU"
#include "Third_Party/Eigen/Geometry"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>

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
			Freestream = Eigen::Matrix<real,3,1>(1,0,0);
		}

		inline void Init(std::string input)
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
			
			/*Initialise matrices*/
			panelData.reserve(npanels);
			aInf = Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic>(npanels,npanels);
			gamma = Eigen::Matrix<real,Eigen::Dynamic,1>(npanels);
			RHS = Eigen::Matrix<real,Eigen::Dynamic,1>(npanels);

			MakeMatrix();
		}

		inline void GetGamma(StateVecD inf)
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

		inline StateVecD getVelocity(StateVecD const& pos)
		{	/*Find velocity for a particle at its position*/
			// Use a Kahan sum to avoid truncation errors in the accumulation
			StateVecD vel = StateVecD::Zero();
			StateVecD c = StateVecD::Zero();// A running compensation for lost low-order bits.

			for (int ii = 0; ii < npanels; ++ii)  
			{
				StateVecD y = gamma(ii)*FindInfluenceConditioned(panelData[ii].A,panelData[ii].B,pos) - c;         // c is zero the first time around.
				StateVecD t = vel + y;              // Alas, sum is big, y small, so low-order digits of y are lost.
				c = (t - vel) - y;            // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
				vel = t;                      // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
				// Next time around, the lost low part will be added to y in a fresh attempt.
			}

			// for(int i = 0; i<npanels; ++i)
			// 	vel += gamma(i)*FindInfluenceConditioned(panelData[i].A,panelData[i].B,pos);
			
			return vel+Freestream;
		}

		inline void write_VLM_Panels(string &prefix)
		{	
			string file1 = prefix;
			file1.append("_Panels.dat");
			ofstream fp(file1, std::ios::out);
			if(fp.is_open())
			{
				fp << "TITLE=\"VLM Panels\"" << std::endl;
		  		fp << "VARIABLES = \"X\", \"Y\", \"Z\""<< std::endl;
		  		for (auto const& p:panelData)
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
		  		for (auto const& p:panelData)
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

		
		template<typename T>
		inline void Get_Number(string const& line, string const& param, T &value)
		{
			size_t pos = line.find(":");
			string substr;
			if(pos != string::npos)
				substr = line.substr(0,pos);

			substr = ltrim(substr);
			if(substr == param)
			{
				string temp = Get_Parameter_Value(line);
				std::istringstream iss(temp);
				iss >> value;
			}
		}


		inline void Get_Vector(string const& line, string const& param, 
					Eigen::Matrix<real,3,1>/* vec<real,3> */ &value)
		{
			size_t pos = line.find(":");
			string substr;
			if(pos != string::npos)
				substr = line.substr(0,pos);

			substr = ltrim(substr);
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
				
				value = /* vec<real,3> */ Eigen::Matrix<real,3,1>(a,b,c);
			}
		}

		inline void Get_Vector(string const& line, string const& param, 
					Eigen::Matrix<real,2,1>/* vec<real,2> */ &value)
		{
			size_t pos = line.find(":");
			string substr;
			if(pos != string::npos)
				substr = line.substr(0,pos);

			substr = ltrim(substr);
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
				
				value = /* vec<real,2> */ Eigen::Matrix<real,2,1>(a,b);
			}
		}

		inline void Get_Vector(string const& line, string const& param, 
					Eigen::Matrix<int,3,1>/* vec<int,3> */ &value)
		{
			size_t pos = line.find(":");
			string substr;
			if(pos != string::npos)
				substr = line.substr(0,pos);

			substr = ltrim(substr);
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
				
				value = /* vec<int,3> */ Eigen::Matrix<int,3,1>(a,b,c);
			}
		}

		inline void Get_Vector(string const& line, string const& param, 
						Eigen::Matrix<int,2,1>& value)
		{
			size_t pos = line.find(":");
			string substr;
			if(pos != string::npos)
				substr = line.substr(0,pos);

			substr = ltrim(substr);
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
				
				value = Eigen::Matrix<int,2,1>(a,b) /* vec<int,2>(a,b) */;
			}
		}

		inline void MakeMatrix(void)
		{
			/*Define the steps for each dimension*/
			real dr0 = coords[1]/real(panels[1]);
			real dy0 = (coords[0]/real(panels[0])); // /cos(sweep);
			real dx0 = dr0; // cos(AoA)*dr0;
			real dz0 = 0.0; // sin(AoA)*dr0;
			
			real x0 = 0.0;
			real y0 = 0.0; // cos(AoA)*coords[1];
			real z0 = 0.0; // sin(AoA)*coords[1];
			
			real yend = coords[0]; // *cos(sweep); /* End is now actually the span, and hypotenuse adjusted for the sweep */
			
			/*Find 1/4 panel points (assuming symmetry) */
			/*Find 3/4 Control Points*/
			for(int i = -panels[0]; i < panels[0]; ++i)
			{	/*i = count in z, (spanwise)*/
				/*X coordinates*/
				real y1 = real(i)*dy0;
				real y2 = real(i+1)*dy0;
				real yhalf = (real(i)+0.5)*dy0;
				
				/*x and y values at j = 0 to apply sweep*/
				real xi = tan(sweep)*(fabs(y1)-y0);
				real zi = 0.0; // sin(sweep)*(fabs(y1)-y0);
				real xip1 = tan(sweep)*(fabs(y2)-y0);
				real zip1 = 0.0; // sin(sweep)*(fabs(y2)-y0);

				/*Vertex 1 deltas (i)*/
				real frac1 = (fabs(y1)-y0)/(yend-y0);
				real fac1 = ((1-frac1) + taper*frac1);
				real dx1 = dx0*fac1;
				real dz1 = dz0*fac1;
				
				/*Vertex 2 deltas (i+1)*/
				real frac2 = (fabs(y2)-y0)/(yend-y0);
				real fac2 = ((1-frac2) + taper*frac2);
				real dx2 =  dx0*fac2;
				real dz2 =  dz0*fac2;

				/*Halfway values (i+1/2) for control point*/
				real frac3 = (fabs(yhalf)-y0)/(yend-y0);
				real fac3 = ((1-frac3) + taper*frac3);
				real xhalf = tan(sweep)*(fabs(yhalf)-y0);
				real zhalf = 0.0; // sin(sweep)*(fabs(yhalf)-y0);
				real dzhalf = dz0*fac3;
				real dxhalf = dx0*fac3;

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
					StateVecD p1(x0 + xi   + dx1*(real(j))  , y1, z0 + zi   + dz1*(real(j))  );
					StateVecD p2(x0 + xi   + dx1*(real(j+1)), y1, z0 + zi   + dz1*(real(j+1)));
					StateVecD p3(x0 + xip1 + dx2*(real(j+1)), y2, z0 + zip1 + dz2*(real(j+1)));
					StateVecD p4(x0 + xip1 + dx2*(real(j))  , y2, z0 + zip1 + dz2*(real(j))  );

					/*1/4 chord points*/
					StateVecD A(x0 + xi   + dx1*(0.25+real(j)), y1, z0 + zi   + dz1*(0.25+real(j)));
					StateVecD B(x0 + xip1 + dx2*(0.25+real(j)), y2, z0 + zip1 + dz2*(0.25+real(j)));
					/*3/4 control point*/
					StateVecD C(x0 + xhalf + dxhalf*(0.75+real(j)), yhalf, z0 + zhalf + dzhalf*(0.75+real(j)));
					panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
				}

				if(doflap == 1)
				{
					/*Flap delta parameters*/
					real dzflap1 = -(sin(beta))*dr0*fac1;
					real dxflap1 = (cos(beta))*dr0*fac1;
					real dzflap2 = -(sin(beta))*dr0*fac2;
					real dxflap2 = (cos(beta))*dr0*fac2;
					real dzflaphalf = -(sin(beta))*dr0*fac3;
					real dxflaphalf = (cos(beta))*dr0*fac3;

					/*Flap start position*/
					real fac = real(panels(1)-flap(2));
					real zflap1 = z0 + zi   + fac*dz1;
					real xflap1 = x0 + xi   + fac*dx1;
					real zflap2 = z0 + zip1 + fac*dz2;
					real xflap2 = x0 + xip1 + fac*dx2;
					real zflaphalf = z0 + zhalf + fac*dzhalf;
					real xflaphalf = x0 + xhalf + fac*dxhalf;

					for (int j = 0; j < flap(2); ++j)
					{
						/*j = count in x and y, (chordwise)*/
						/*Panel Verticies (for visual)*/
						StateVecD p1(xflap1+dxflap1*(real(j))  , y1, zflap1 + dzflap1*(real(j))  );
						StateVecD p2(xflap1+dxflap1*(real(j+1)), y1, zflap1 + dzflap1*(real(j+1)));
						StateVecD p3(xflap2+dxflap2*(real(j+1)), y2, zflap2 + dzflap2*(real(j+1)));
						StateVecD p4(xflap2+dxflap2*(real(j))  , y2, zflap2 + dzflap2*(real(j))  );

						/*1/4 chord points*/
						StateVecD A(xflap1+dxflap1*(0.25+real(j)), y1, zflap1 + dzflap1*(0.25+real(j)));
						StateVecD B(xflap2+dxflap2*(0.25+real(j)), y2, zflap2 + dzflap2*(0.25+real(j)));
						/*3/4 control point*/
						StateVecD C(xflaphalf+dxflaphalf*(0.75+real(j)), yhalf, zflaphalf + dzflaphalf*(0.75+real(j)));
						panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
					}
				}
			}
		
			/*End initialisation*/
		}

		/* May need some conditioning for points that are nearly colinear */
		inline const StateVecD FindInfluence(StateVecD const& A, StateVecD const& B, StateVecD const& C)
 		{	
			StateVecD r0, r1, r2, inf;

			/*Bounded Vortex*/
			StateVecD coefAB;

			r0 = B-A;
			r1 = C-A;
			r2 = C-B;

			coefAB = (1/(4*M_PI))*((r1.cross(r2))/((r1.cross(r2)).squaredNorm()))*
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


		/* May need some conditioning for points that are nearly colinear */
		/* Current test - If any vortex is ill conditioned, ignore the entire horseshoe */
		/* Feels more appropriate to do this, but could be excessive. */
		inline const StateVecD FindInfluenceConditioned(StateVecD const& A, StateVecD const& B, StateVecD const& C)
 		{	
			StateVecD r0, r1, r2, inf;

			StateVecD coefAB, coefA, coefB;

			real const cutoff = 1e-10;

			r0 = B-A;
			r1 = C-A;
			r2 = C-B;

			// Bounded Vortex
			StateVecD numer = r1.cross(r2);
			double denom = numer.squaredNorm();
			// Cut off if denomenator is too small
			if(denom > cutoff)
				coefAB = (1/(4*M_PI))*(numer/denom)*
					(r0.dot(r1.normalized()-r2.normalized()));


			// Horseshoe vortex from point A to infinity		
			inf = A + Freestream;
			r2 = C - A;
			r1 = C - inf;
			r0 = A - inf;

			numer = r1.cross(r2);
			denom = numer.squaredNorm();
			// Cut off if denomenator is too small
			if(denom > cutoff)
				coefA = (1/(4*M_PI))*r0.norm()*(numer/denom)*
						(1-r0.dot(r2)/(r0.norm()*r2.norm()));

			// Horseshoe vortex from point B to infinity
			inf = B + Freestream;
			r1 = C - B;
			r2 = C - inf;
			r0 = inf - B;

			numer = r1.cross(r2);
			denom = numer.squaredNorm();
			// Cut off if denomenator is too small
			if(denom > cutoff)
				coefB = (1/(4*M_PI))*r0.norm()*(numer/denom)*
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