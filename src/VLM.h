/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef VLM_H
#define VLM_H

#include "Var.h"
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>

/*A structure for the coordinates to define a panel.*/
typedef struct Panel
{ /*A and B are 1/4 chord bounds for the vortex.
	C is the control Point location*/
#if SIMDIM == 3
	Panel(StateVecD A, StateVecD B, StateVecD C,
			StateVecD p1, StateVecD p2, StateVecD p3, StateVecD p4)
	: A(A), B(B), C(C), p1(p1), p2(p2), p3(p3), p4(p4) 
	{norm =  (C-A).cross(C-B).normalized();}
#endif

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

			write_traj = 0;
			maxX = 1000;
			streamDx = 0.1;
			maxIters = 10000;
			minDx = 1e-5;
		}

		void Init(std::string input);

		void GetGamma(StateVecD inf);

		const StateVecD getVelocity(StateVecD const& pos);

		void write_VLM_Panels(string &prefix);

		void Plot_Streamlines(std::string& prefix);

		int write_traj;
	protected:
		
		void MakeMatrix(void);

		/* May need some conditioning for points that are nearly colinear */
		const StateVecD FindInfluence(StateVecD const& A, StateVecD const& B, StateVecD const& C);

		/* May need some conditioning for points that are nearly colinear */
		/* Current test - If any vortex is ill conditioned, ignore the entire horseshoe */
		/* Feels more appropriate to do this, but could be excessive. */
		const StateVecD FindInfluenceConditioned(StateVecD const& A, StateVecD const& B, StateVecD const& C);

 		std::vector<Panel> panelData;

		Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic> aInf;

		Eigen::Matrix<real,Eigen::Dynamic,1> gamma;
		Eigen::Matrix<real,Eigen::Dynamic,1> RHS;

		Eigen::Matrix<real,3,1> Freestream;
		Eigen::Matrix<real,2,1> coords;
		Eigen::Vector2i panels;
		Eigen::Vector3i flap;

		real AoA, sweep, taper, beta;
		real maxX, streamDx, minDx;
		size_t maxIters;
		int npanels, nverts;
}VLM;

#endif