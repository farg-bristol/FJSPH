/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VLM_H
#define VLM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include "../Eigen/Core"
#include "../Eigen/StdVector"
#include "../Eigen/LU"
#include "../Eigen/Geometry"
// #include "IO.h"

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif


/*A structure for the coordinates to define a panel.*/
typedef struct Panel
{ /*A and B are 1/4 chord bounds for the vortex.
	C is the control Point location*/
	Panel(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C,
			Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d p4)
	: A(A), B(B), C(C), p1(p1), p2(p2), p3(p3), p4(p4) 
	{norm =  (C-A).cross(C-B).normalized();}

	Eigen::Vector3d A, B, C;
	Eigen::Vector3d p1, p2, p3, p4;
	Eigen::Vector3d norm;
}Panel;

int getI(std::ifstream& In)
{
	std::string line;
	getline(In,line);
	int i = stoi(line);
	return i;
}

double getD(std::ifstream& In)
{
	std::string line;
	getline(In,line);
	double d = stod(line);
	return d; 
}

std::string getS(std::ifstream& In)
{
	std::string line;
	getline(In,line);
	return line; 
}

Eigen::Vector2i getIVec(std::ifstream& In)
{
	std::string line;
	getline(In,line);
	std::istringstream sline(line);
	// cout << sline.str() << endl;
	Eigen::Vector2i x;
	sline >> x[0]; sline >> x[1];

	return x;
}

/*Function for a 2D Vector*/
Eigen::Vector2d getvec(std::ifstream& In)
{
	std::string line;
	getline(In,line);
	std::istringstream sline(line);
	
	Eigen::Vector2d x;
	sline >> x[0]; sline >> x[1]; 
		
	return x;
}

Eigen::Vector3i get3dVector(std::ifstream& In)
{
	std::string line;
	getline(In,line);
	std::istringstream sline(line);
	
	Eigen::Vector3i x;
	sline >> x[0]; sline >> x[1]; sline >> x[2]; 
		
	return x;
}


typedef class VLM
{
	public:
		VLM()
		{
			std::ifstream filein("VLM.dat", std::ios::in);

			/*Define x and y end coordinates of the wing*/
			coords = getvec(filein);

			/*Split it up into this many panels (Will be doubled on the other side)*/
			panels = getIVec(filein);
			npanels = 2*panels[0]*panels[1];

			/*Define Angle of attack*/
			AoA = getD(filein) * M_PI/180.0;
			
			/*Sweep*/
			sweep = getD(filein) * M_PI/180.0;
			taper = getD(filein);

			/*Flap Properties*/
			flap = get3dVector(filein);
			beta = getD(filein) * M_PI/180.0;

			// coords[0] = 5;
			// coords[1] = 2;
			// panels(0) = 10;
			// panels(1) = 5;
			// AoA = 12.0 * M_PI/180.0;
			// sweep = 30 * M_PI/180.0;
			
			/*Initialise matrices*/
			panelData.reserve(npanels);
			aInf = Eigen::MatrixXd(npanels,npanels);
			gamma = Eigen::VectorXd(npanels);
			RHS = Eigen::VectorXd(npanels);

			
			// panelxyz.reserve(nverts);

			/*Want freestream to be aligned with jet axis*/
			Freestream[0]= 0;
			Freestream[1]= 1;
			Freestream[2]= 0;

			Freestream = Freestream.normalized();

			MakeMatrix();

		}

		void GetGamma(Eigen::Vector3d inf)
		{
			Freestream = inf;
			/*Find the inflence matrix, and invert it to find gamma*/

			/*Find the influence matrix aInf*/
 			for(int i=0; i < npanels; ++i)
 			{
 				for(int j=0; j < npanels; ++j)
 				{
 					/*aInf[i,j] = the influence of vortex j on control point i*/
 					Eigen::Vector3d A, B, C;
 					A = panelData[i].A;
 					B = panelData[i].B;
 					C = panelData[j].C;

					aInf(j,i) = FindInfluence(A,B,C).dot(panelData[i].norm);
				}

			RHS(i) = -Freestream.dot(panelData[i].norm);
			}
			/*Find Gamma */
 			gamma = aInf.fullPivLu().solve(RHS);

 			// for (int i= 0; i < gamma.rows(); ++i)
 			// {
 			// 	std::cout << gamma[i] << std::endl;
 			// }
 		}

		Eigen::Vector3d getVelocity(Eigen::Vector3d pos)
		{	/*Find velocity for a particle at its position*/
			Eigen::Vector3d vel = Eigen::Vector3d::Zero();

			for(int i = 0; i<npanels; ++i)
				vel += gamma(i)*FindInfluence(panelData[i].A,panelData[i].B,pos);
			
			return vel+Freestream;
		}

		void write_VLM_Panels(std::string folder)
		{
			std::string file1 = folder;
			file1.append("/VLM_Panels.plt");
			std::ofstream fp(file1, std::ios::out);
	  		fp << "TITLE=\"VLM Panels\"" << std::endl;
	  		fp << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\""<< std::endl;
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
		  	std::string file2 = folder;
		  	file2.append("/VLM_Vortices.plt");
		  	std::ofstream fq(file2, std::ios::out);
	  		fq << "TITLE=\"VLM Vortices and Control Points\"" << std::endl;
	  		fq << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\""<< std::endl;
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

	protected:
		void MakeMatrix(void)
		{
			/*Define the steps for each dimension*/
			double dr0 = coords[1]/double(panels[1]);
			double dz0 = coords[0]/double(panels[0])*cos(sweep);
			double dy0 = cos(AoA)*dr0;
			double dx0 = sin(AoA)*dr0;

			double y0 = cos(AoA)*coords[1];
			double x0 = sin(AoA)*coords[1];
			double z0 = 0.0;
			double zend = coords[0]*cos(sweep);
			
			/*Find 1/4 panel points (assuming symmetry) */
			/*Find 3/4 Control Points*/
			for(int i = -panels[0]; i < panels[0]; ++i)
			{	/*i = count in z, (spanwise)*/
				/*Z coordinates*/
				double z1 = double(i)*dz0;
				double z2 = double(i+1)*dz0;
				double zhalf = (double(i)+0.5)*dz0;
				
				/*x and y values at j = 0 to apply sweep*/
				double yi = tan(sweep)*(fabs(z1)-z0);
				double xi = sin(AoA)*sin(sweep)*(fabs(z1)-z0);
				double yip1 = tan(sweep)*(fabs(z2)-z0);
				double xip1 = sin(AoA)*sin(sweep)*(fabs(z2)-z0);

				/*Vertex 1 deltas (i)*/
				double frac1 = (fabs(z1)-z0)/(zend-z0);
				double dy1 =  dy0*(1-frac1) + dy0*taper*frac1;
				double dx1 =  dx0*(1-frac1) + dx0*taper*frac1;
				
				/*Vertex 2 deltas (i+1)*/
				double frac2 = (fabs(z2)-z0)/(zend-z0);
				double dy2 =  dy0*(1-frac2) + dy0*taper*frac2;
				double dx2 =  dx0*(1-frac2) + dx0*taper*frac2;

				/*Halfway values (i+1/2) for control point*/
				double frac3 = (fabs(zhalf)-z0)/(zend-z0);
				double yhalf = tan(sweep)*(fabs(zhalf)-z0);
				double xhalf = sin(AoA)*sin(sweep)*(fabs(zhalf)-z0);
				double dxhalf = dx0*(1-frac3) + dx0*taper*frac3;
				double dyhalf =  dy0*(1-frac3) + dy0*taper*frac3;

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
					Eigen::Vector3d p1(x0-xi-dx1*(double(j)),dy1*(double(j))-y0+yi,z1);
					Eigen::Vector3d p2(x0-xi-dx1*(double(j+1)),dy1*(double(j+1))-y0+yi,z1);
					Eigen::Vector3d p3(x0-xip1-dx2*(double(j+1)),dy2*(double(j+1))-y0+yip1,z2);
					Eigen::Vector3d p4(x0-xip1-dx2*(double(j)),dy2*(double(j))-y0+yip1,z2);

					/*1/4 chord points*/
					Eigen::Vector3d A(x0-xi-dx1*(0.25+double(j)),dy1*(0.25+double(j))-y0+yi,z1);
					Eigen::Vector3d B(x0-xip1-dx2*(0.25+double(j)),dy2*(0.25+double(j))-y0+yip1,z2);
					/*3/4 control point*/
					Eigen::Vector3d C(x0-xhalf-dxhalf*(0.75+double(j)),dyhalf*(0.75+double(j))-y0+yhalf,zhalf);
					panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
				}

				if(doflap == 1)
				{
					/*Flap delta parameters*/
					double dxflap1 = (sin(beta+AoA))*dr0*(1-frac1) + (sin(beta+AoA))*dr0*taper*frac1;
					double dyflap1 = (cos(beta+AoA))*dr0*(1-frac1) + (cos(beta+AoA))*dr0*taper*frac1;
					double dxflap2 = (sin(beta+AoA))*dr0*(1-frac2) + (sin(beta+AoA))*dr0*taper*frac2;
					double dyflap2 = (cos(beta+AoA))*dr0*(1-frac2) + (cos(beta+AoA))*dr0*taper*frac2;
					double dxflaphalf = (sin(beta+AoA))*dr0*(1-frac3) + (sin(beta+AoA))*dr0*taper*frac3;
					double dyflaphalf = (cos(beta+AoA))*dr0*(1-frac3) + (cos(beta+AoA))*dr0*taper*frac3;

					/*Flap start position*/
					double xflap1 = x0-xi-double(panels(1)-flap(2))*dx1;
					double yflap1 = double(panels(1)-flap(2))*dy1 + yi - y0;
					double xflap2 = x0 - xip1-double(panels(1)-flap(2))*dx2;
					double yflap2 = double(panels(1)-flap(2))*dy2 + yip1 - y0;
					double xflaphalf = x0-xhalf-double(panels(1)-flap(2))*dxhalf;
					double yflaphalf = double(panels(1)-flap(2))*dyhalf + yhalf - y0;

					for (int j = 0; j < flap(2); ++j)
					{
						/*j = count in x and y, (chordwise)*/
						/*Panel Verticies (for visual)*/
						Eigen::Vector3d p1(xflap1-dxflap1*(double(j)),yflap1+dyflap1*(double(j)),z1);
						Eigen::Vector3d p2(xflap1-dxflap1*(double(j+1)),yflap1+dyflap1*(double(j+1)),z1);
						Eigen::Vector3d p3(xflap2-dxflap2*(double(j+1)),yflap2+dyflap2*(double(j+1)),z2);
						Eigen::Vector3d p4(xflap2-dxflap2*(double(j)),yflap2+dyflap2*(double(j)),z2);

						/*1/4 chord points*/
						Eigen::Vector3d A(xflap1-dxflap1*(0.25+double(j)),yflap1+dyflap1*(0.25+double(j)),z1);
						Eigen::Vector3d B(xflap2-dxflap2*(0.25+double(j)),yflap2+dyflap2*(0.25+double(j)),z2);
						/*3/4 control point*/
						Eigen::Vector3d C(xflaphalf-dxflaphalf*(0.75+double(j)),yflaphalf+dyflaphalf*(0.75+double(j)),zhalf);
						panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
					}
				}
			}
		
			/*End initialisation*/
		}

		Eigen::Vector3d FindInfluence(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
 		{	
			Eigen::Vector3d r0, r1, r2, inf;

			/*Bounded Vortex*/
			Eigen::Vector3d coefAB;

			r0 = B-A;
			r1 = C-A;
			r2 = C-B;

			coefAB = (1/(4*M_PI))*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(r0.dot(r1.normalized())-r0.dot(r2.normalized()));


			/*Horseshoe vortex from point A*/
			Eigen::Vector3d coefA;
		
			inf = A + Freestream;
			r2 = C - A;
			r1 = C - inf;
			r0 = A - inf;

			coefA = (1/(4*M_PI))*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
					(1-r0.dot(r2)/(r0.norm()*r2.norm()));

			/*Horseshoe vortex from point B*/
			Eigen::Vector3d coefB;

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
		double AoA, sweep, taper, beta;
		
		Eigen::MatrixXd aInf;

		Eigen::VectorXd gamma;
		Eigen::VectorXd RHS;

		Eigen::Vector2d coords;
		Eigen::Vector2i panels;
		
		Eigen::Vector3i flap;

		Eigen::Vector3d Freestream;
}VLM;

#endif