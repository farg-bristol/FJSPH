/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*********                   FOR THREE DIMENSION CODE                           *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"
#include "IO.h"
#include "Cross.h"

using namespace std;

/*Make a guess on how big the array is going to be (doesn't need to be totally exact)*/
int ParticleCount(SIM &svar)
{
	int partCount = 0;
	ldouble stepx = svar.Pstep*svar.Bstep;
	ldouble stepy = svar.Pstep*svar.Bstep;
	int Ny = ceil(svar.Box(1)/stepy);
	int Nx = ceil(svar.Box(0)/stepx);
	int Nz = ceil(svar.Box(2)/stepx);


	switch(svar.Bcase)
	{
		case 0:
			partCount += svar.simPts; /*Simulation pn*/
			break;
		case 1:
			partCount = Nx*Nz + 2*Nx*Ny + 2*Ny*Nz; /*Boundary particles*/
			partCount += svar.simPts; /*Simulation pn*/
			break;
		case 3:
		{	/*Boundary Points*/
			ldouble holeS = svar.Start(0); /*Distance from Box start to hole*/
			ldouble holeD = svar.Start(1)+4*svar.Pstep; /*Diameter of hole (or width)*/
			
			ldouble stepb = (svar.Pstep*svar.Bstep);
			int Nx = ceil((holeS)/stepb);
			stepb = holeS/(1.0*Nx);
			int Nz = ceil((svar.Box(2)/stepb));

			/* Trivial values*/			
			int preHole = (Nx+1)*Nz;
			int postHole = (ceil((svar.Box(0)-holeS-holeD)/stepx)+1)*Nz;

			/*Find how many points would be there without the hole*/
            int minusHole = (ceil(holeD/stepb)+1)*Nz;

            /*Find how many points would be in the circle*/
            int holeCount = 0;
            for (ldouble x = -0.5*holeD; x <= 0.5*holeD; x+=stepb)
            {
            	for(ldouble z = -0.5*holeD; z <= 0.5*holeD; z+=stepb)
            	{
	            	if((x*x)/(0.25*holeD*holeD) + (z*z)/(0.25*holeD*holeD) < 1 )
	                        ++holeCount;
            	}
            }

            /*Find the points on the side of the pipe (Assume a circle)*/
            float dtheta = atan((stepb)/(0.5*holeD));
            int holeWall = ceil(2*M_PI/dtheta)*Ny;

            int holeBound = minusHole - holeCount + holeWall;

			/*Simulation Points*/
			int simCount = 0;
			holeD = svar.Start(1);
			for (ldouble x = -0.5*holeD; x <= 0.5*holeD; x+=stepb)
            {
            	for(ldouble z = -0.5*holeD; z <= 0.5*holeD; z+=stepb)
            	{
	            	if((x*x)/(0.25*holeD*holeD) + (z*z)/(0.25*holeD*holeD) < 1 )
	                        ++simCount;
            	}
            }

			/*Need to add the pn already present*/
			int simPts = simCount*svar.nmax + simCount*ceil(svar.Box[1]/svar.Pstep);

			partCount = preHole + holeBound + postHole + simPts;
			
			break;
		}

		case 4:
		{
			unsigned int count = 0;
			/*Create the simulation pn*/
			for(double z=0; z < 1.5; z+=svar.Pstep) 
			{
				for (double y = -sqrt(z); y <= sqrt(z); y+=svar.Pstep)
				{				
					for (double x = -sqrt(z); x <= sqrt(z); x+=svar.Pstep)
					{	
						if (pow(x,2) + pow(y,2) <= z)
						{
							++count;
						}				
					}	
				}
			}

			partCount = count;
			int decision = 0;
			cout << "Particle count: " << count << endl;
			cout << "Continue?" << endl;
			cin >>  decision;

			if (decision == 1)
				break;
			else if (decision == 0)
				exit(0);

		}

			
	}
	
	return partCount;
}

void InitSPH(SIM &svar, FLUID &fvar, CROSS &cvar, State &pn, State &pnp1)
{


	switch (svar.Bcase)
	{
		default:
			cout << "Initialising simulation with " << svar.simPts << " points" << endl;
			break;

		case 3:
			cout << "Initialising simulation..." << endl;
			break;
	}
	
	
	//Structure input initialiser
	//Particle(StateVecD x, StateVecD v, StateVecD vh, StateVecD f, float rho, float Rrho, bool bound) :
	StateVecD v = StateVecD::Zero();  
	StateVecD f = StateVecD::Zero(); 
	ldouble rho=fvar.rho0;  
	 
/************** Create the boundary pn  *****************/ 	 
	ldouble stepx = svar.Pstep*svar.Bstep;
	ldouble stepy = svar.Pstep*svar.Bstep;
	ldouble stepz = svar.Pstep*svar.Bstep;
	
	int Ny = ceil(svar.Box(1)/stepy);
	stepy = svar.Box(1)/Ny;	/*Find exact step to meet dimensions*/
	
	int Nx = ceil(svar.Box(0)/stepx);
	stepx = svar.Box(0)/Nx;

	int Nz = ceil(svar.Box(2)/stepz);
	stepz = svar.Box(2)/Nz;
	
	switch (svar.Bcase) 
	{
		case 0:
		{ /*No boundary*/
			break;
		}
		case 1: /*Rectangle*/
		{
			for(int j = 0; j <= Ny ; ++j) 
			{
				for (int k=0; k<=Nz; ++k) 
				{	/*Create Left and right boundary faces*/
					StateVecD xi(0.0,j*stepy,k*stepz);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
					xi(0) = svar.Box(0);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
				}
				
			}
			for(int i = 1; i<Nx ; ++i) 
			{
				for (int j=1; j<Ny; ++j)
				{	/*Create top and bottom boundary faces*/
					StateVecD xi(i*stepx,j*stepy,0);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
					// xi(2)= Box(2); //Top boundary (Typically omitted)
					// pn.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,0));
				}
				
			}
			for(int i= 1; i<Nx; ++i) 
			{
				for(int k = 0; k <= Nz; ++k) 
				{	/*Create far and near boundary*/
					StateVecD xi(i*stepx, 0, k*stepz);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
					xi(1) = svar.Box(1);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
				}
			}		
			break;
		}
		// case 2: /*Bowl*/
		// {	/*DONT USE WITH SIMDIM = 3*/
			
		// 	ldouble r= svar.Box(0);
		// 	ldouble dtheta = atan((svar.Pstep*svar.Bstep)/r);
		// 	for (ldouble theta = 0.0; theta < M_PI; theta+=dtheta)
		// 	{
		// 		StateVecD xi(-r*cos(theta)-svar.Start(0),r*(1-sin(theta))-svar.Start(1));
		// 		pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
		// 	}
		// 	break;
		// }
		case 3:
		{	/*Jet in Crossflow*/
			ldouble holeS = svar.Start(0); /*Distance from Box start to hole*/
			ldouble holeD = svar.Start(1)+4*fvar.H; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			int Nb = ceil((holeS)/stepb);
			stepb = holeS/(1.0*Nb);
			
			// cout << "Bound x-cond: " << svar.Box(0) << " Y-cond: " << svar.Box(1) << " z-cond: " << svar.Box(2) << endl;

			/*Pre Hole*/
			for (ldouble x=0.0; x<holeS; x+= stepb)
			{
				for (ldouble z = -0.5*svar.Box(2); z <= 0.5*svar.Box(2); z+=stepb)
				{
					StateVecD xi(x,0.0,z);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
				}
			}

			/*Plate around the hole*/
			for (ldouble x=holeS; x<holeS+holeD; x+= stepb)
			{
				for (ldouble z = -0.5*svar.Box(2); z <= 0.5*svar.Box(2); z+=stepb)
				{
					ldouble r = (x-holeS-0.5*holeD); /*Normalise the circle around 0,0*/
					if((r*r)/(0.25*holeD*holeD) + (z*z)/(0.25*holeD*holeD) >= 1.0 )
					{	/*If the point is outside the hole diameter, add it*/
						StateVecD xi(x,0.0,z);
						pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
					}
				}
			}

			/*Post Hole*/
			for (ldouble x=holeS+holeD; x<svar.Box(0); x+= stepb)
			{
				for (ldouble z = -0.5*svar.Box(2); z <= 0.5*svar.Box(2); z+=stepb)
				{
					StateVecD xi(x,0.0,z);
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
				}
			}


			/*Create a bit of the pipe downward.*/
			double r = 0.5*holeD;
        	double dtheta = atan((stepb)/(r));
			for (ldouble y = -stepb; y >= -svar.Box(1)-stepb; y-=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*(1-cos(theta))+holeS, y, r*sin(theta));
					pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
				}
				
			}
				
			break;
		}
		case 4:
		{
			// ifstream afile("23015.csv");
			break;
			
		}
		default: 
		{
			cerr << "Boundary case is not within the design. 0 <= Bcase <= 4." << endl;
			exit(-1);
		}
	}
	
	svar.bndPts = pn.size();
	
/***********  Create the simulation pn  **************/

	
	switch(svar.Bcase)
	{
		case 3: 
		{	/*Crossflow case*/
		// cout << "Got to making simulation points..." << endl;
			svar.simPts = 0;
			svar.totPts = pn.size();
			/*Update n+1 before adding sim pn*/
			for (auto p: pn)
				pnp1.emplace_back(p);

			for (ldouble y = 0.0; y > -svar.Box[1]; y-=svar.Pstep)
			{
				// cout << "In add points for-loop" << endl;
				AddPoints(y, svar, fvar, cvar, pn, pnp1);
			}

			break;
		}

		case 4:
		{
			/*Create the simulation pn*/
			for(double z=0; z < 1.5; z+=svar.Pstep) 
			{
				
				for (double y = -sqrt(z); y <= sqrt(z); y+=svar.Pstep)
				{				
					for (double x = -sqrt(z); x <= sqrt(z); x+=svar.Pstep)
					{	
						if (pow(x,2) + pow(y,2) <= z)
						{
							StateVecD xi(svar.Start(0)+x,svar.Start(1)+y,svar.Start(2)+z);		
							pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
						}				
					}	
				}
			}

			for (auto p: pn)
				pnp1.emplace_back(p);

			svar.simPts = pnp1.size();

			break;
		}

		default:
		{
			/*Create the simulation pn*/
			for( int i=0; i< svar.xyPART(0); ++i) 
			{
				for(int j=0; j< svar.xyPART(1); ++j)
				{				
					for (int k=0; k < svar.xyPART(2); ++k )
					{
						StateVecD xi(svar.Start(0)+i*svar.Pstep,svar.Start(1)+j*svar.Pstep,svar.Start(2)+k*svar.Pstep);		
						pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
					}
						
				}
			}

			for (auto p: pn)
				pnp1.emplace_back(p);

			break;
		}
	}
		
	// svar.simPts+=10*10;

	svar.totPts = pn.size();
	// cout << "Boundary pn: " << svar.bndPts << endl;
	// cout << "Sim pn:" << svar.simPts << endl;
 // 	cout << "Total Points: " << svar.totPts << endl;
	
	// if(svar.totPts!=svar.bndPts+svar.simPts)
	// {
	// 	cerr<< "Mismatch of particle count." << endl;
	// 	cerr<< "Particle array size doesn't match defined values." << endl;
	// 	write_settings(svar,fvar);
	// 	exit(-1);
	// }
	// cout << "Total pn: " << svar.totPts << endl;
	
	// cout << "Refresh round pn: " << svar.nrefresh << endl;
}


#endif