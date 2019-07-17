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
		{	
			ldouble holeD = svar.Jet(0)+4*svar.Pstep; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);

            /*Find the points on the side of the pipe (Assume a circle)*/
            float dtheta = atan((stepb)/(0.5*holeD));
            Ny = ceil(svar.Jet(1)/stepb);
            int holeWall = ceil(2*M_PI/dtheta)*Ny;

			/*Simulation Points*/
			int simCount = 0;
			ldouble jetR = 0.5*(svar.Jet(0));
			
			/*Do the centerline of points*/
			for (ldouble z = -jetR; z <= jetR; z+= svar.dx)
				simCount++;

			for (ldouble x = svar.dx; x < jetR ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				for (ldouble z = -jetR; z <= jetR; z+= svar.dx)
				{	/*If the point is inside the hole diameter, add it*/
					if(((x*x)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
		    			simCount += 2;
				}
			}

			/*Need to add the pn already present*/
			int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

			partCount = holeWall + simPts;
			
			break;
		}

		case 4:
		{	
			ldouble holeD = svar.Jet(0)+4*svar.Pstep; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);

            /*Find the points on the side of the pipe (Assume a circle)*/
            double dtheta = atan((stepb)/(0.5*holeD));
            Ny = ceil(svar.Jet(1)/stepb);
            int holeWall = ceil(2*M_PI/dtheta)*Ny;

			/*Simulation Points*/
			int simCount = 0;
			ldouble jetR = 0.5*(svar.Jet(0));
			
			/*Do the centerline of points*/
			for (ldouble z = -jetR; z <= jetR; z+= svar.dx)
				simCount++;

			for (ldouble x = svar.dx; x < jetR ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				for (ldouble z = -jetR; z <= jetR; z+= svar.dx)
				{	/*If the point is inside the hole diameter, add it*/
					if(((x*x)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
		    			simCount += 2;
				}
			}

			/*Need to add the pn already present*/
			int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

			partCount = holeWall + simPts;
			
			break;
		}

		case 5:
		{
			uint simCount = 0;
			ldouble radius = 0.5*svar.Start(1);
			for (ldouble y = -radius; y <= radius; y+=svar.dx)
			{	
				ldouble xradius = sqrt(radius*radius - y*y);
				for (ldouble z = -xradius; z <= xradius; z+= svar.dx)
				{ /*Do the centerline of points*/
					++simCount;
				}

				for (ldouble x = svar.dx; x <= xradius ; x+=svar.dx)
				{ /*Do the either side of the centerline*/
					for (ldouble z = -xradius; z <= xradius; z+= svar.dx)
					{
						if(((x*x) + (z*z) + (y*y)) <= (radius*radius))
			    		{   /*If the point is inside the hole diameter, add it*/
							++simCount;
							++simCount;
						}
					}	
				}
			}
			partCount = simCount;
			break;
		}		
	}
	
	return partCount;
}

void InitSPH(SIM &svar, FLUID &fvar, CROSS &cvar, State &pn, State &pnp1)
{
	switch (svar.Bcase)
	{
		case 3:
			cout << "Initialising simulation..." << endl;
			break;
		case 4:
			cout << "Initialising simulation..." << endl;
			break;
		default:
			cout << "Initialising simulation with " << svar.simPts << " points" << endl;
			break;

		
	}
	
	
	//Structure input initialiser
	//Particle(StateVecD x, StateVecD v, StateVecD vh, StateVecD f, float rho, float Rrho, bool bound) :
	StateVecD v = StateVecD::Zero();   
	ldouble rho=fvar.rho0;
	ldouble press = fvar.B*(pow(rho/fvar.rho0,fvar.gam)-1);  
	 
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
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
					xi(0) = svar.Box(0);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
				}
				
			}
			for(int i = 1; i<Nx ; ++i) 
			{
				for (int j=1; j<Ny; ++j)
				{	/*Create top and bottom boundary faces*/
					StateVecD xi(i*stepx,j*stepy,0);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
					// xi(2)= Box(2); //Top boundary (Typically omitted)
					// pn.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,0));
				}
				
			}
			for(int i= 1; i<Nx; ++i) 
			{
				for(int k = 0; k <= Nz; ++k) 
				{	/*Create far and near boundary*/
					StateVecD xi(i*stepx, 0, k*stepz);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
					xi(1) = svar.Box(1);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
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
			
			ldouble holeD = svar.Jet(1)+4*fvar.H; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			double r = 0.5*holeD;
        	double dtheta = atan((stepb)/(r));
			for (ldouble y = -stepb; y >= -svar.Jet(1)-stepb; y-=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(-r*(1-cos(theta)), y, r*sin(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
				}
				
			}
				
			break;
		}
		case 4:
		{	/*Jet in VLM*/
			ldouble holeD = svar.Jet(1)+4*fvar.H; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			double r = 0.5*holeD;
        	double dtheta = atan((stepb)/(r));
			for (ldouble y = -stepb; y >= -svar.Jet(1)-stepb; y-=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
				}
				
			}
				
			break;
		}
		case 5: 
		{	/*Droplet case - no boundary*/
			break;
		}
		default: 
		{
			cerr << "Boundary case is not within the design. 0 <= Bcase <= 5." << endl;
			exit(-1);
		}
	}
	
	svar.bndPts = pn.size();
	
/***********  Create the simulation pn  **************/

	
	switch(svar.Bcase)
	{
		case 3: 
		{	/*Crossflow case*/
			svar.simPts = 0;
			svar.totPts = pn.size();
			/*Update n+1 before adding sim pn*/
			for (auto p: pn)
				pnp1.emplace_back(p);

			for (ldouble y = 0.0; y > -svar.Jet[1]; y-=svar.dx)
			{
				// cout << "In add points for-loop" << endl;
				AddPoints(y, svar, fvar, cvar, pn, pnp1);
			}

			break;
		}

		case 4:
		{
			/*VLM case*/
			svar.simPts = 0;
			svar.totPts = pn.size();
			/*Update n+1 before adding sim pn*/
			for (auto p: pn)
				pnp1.emplace_back(p);

			for (ldouble y = 0.0; y > -svar.Jet[1]; y-=svar.dx)
			{
				// cout << "In add points for-loop" << endl;
				AddPoints(y, svar, fvar, cvar, pn, pnp1);
			}

			break;
		}
		case 5:
		{
			/*Droplet case*/
			svar.simPts = 0;
			svar.totPts = pn.size();
			/*Update n+1 before adding sim pn*/
			for (auto p: pn)
				pnp1.emplace_back(p);

			CreateDroplet(svar,fvar,pn,pnp1);

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
						pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
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
	
	if(svar.totPts!=svar.bndPts+svar.simPts)
	{
		cerr<< "Mismatch of particle count." << endl;
		cerr<< "Particle array size doesn't match defined values." << endl;
		write_settings(svar,fvar);
		exit(-1);
	}
	// cout << "Total pn: " << svar.totPts << endl;
	
	// cout << "Refresh round pn: " << svar.nrefresh << endl;
}


#endif