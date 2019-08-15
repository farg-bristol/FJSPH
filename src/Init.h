/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*********                   FOR THREE DIMENSION CODE                           *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"
#include "IO.h"
#include "Add.h"

using namespace std;

/*Make a guess on how big the array is going to be (doesn't need to be totally exact)*/
int ParticleCount(SIM &svar)
{
	int partCount = 0;
	ldouble step = svar.Pstep*svar.Bstep;
	int Ny = ceil(svar.Box(1)/step);
	int Nx = ceil(svar.Box(0)/step);

	#if(SIMDIM == 3)
		int Nz = ceil(svar.Box(2)/step);
	#endif

	switch(svar.Bcase)
	{
		case 0:
			partCount += svar.simPts; /*Simulation pn*/
			break;
		case 1:
			#if (SIMDIM == 3)
				partCount = Nx*Nz + 2*Nx*Ny + 2*Ny*Nz; /*Boundary particles*/
			#else
				partCount = 2*Ny + Nx; /*Boundary particles*/
			#endif
			partCount += svar.simPts; /*Simulation pn*/

			break;
		case 3:
		{	
			#if SIMDIM == 3
			{
				ldouble holeD = svar.Jet(0)+4*svar.Pstep; /*Diameter of hole (or width)*/

	            /*Find the points on the side of the pipe (Assume a circle)*/
	            float dtheta = atan((step)/(0.5*holeD));
	            Ny = ceil(svar.Jet(1)/step);
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
			}
			#else
			{
	            /*Find the points on the side of the pipe*/
	            Ny = ceil(svar.Jet(1)/step);
	            int holeWall = 2*Ny;

				/*Simulation Points*/
				int simCount = 0;
				ldouble jetR = 0.5*(svar.Jet(0));
				for (ldouble x = -jetR; x <= jetR; x+= svar.dx)
					simCount++;

				/*Need to add the pn already present*/
				int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

				partCount = holeWall + simPts;
			}
			#endif
			break;
		}

		case 4:
		{	
			#if SIMDIM == 3
				ldouble holeD = svar.Jet(0)+4*svar.Pstep; /*Diameter of hole (or width)*/

	            /*Find the points on the side of the pipe (Assume a circle)*/
	            double dtheta = atan((step)/(0.5*holeD));
	            Ny = ceil(svar.Jet(1)/step);
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
			#else
				cout << "VLM case is not available in 2D. Stopping." << endl;
				exit(-1);
			#endif
			break;
		}

		case 5:
		{	
			#if SIMDIM == 3
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
			#else
				uint simCount = 0;
				ldouble radius = 0.5*svar.Start(1);
				for (ldouble y = -radius; y <= radius; y+=svar.dx)
				{	
					++simCount;
					

					for (ldouble x = svar.dx; x <= radius ; x+=svar.dx)
					{ /*Do the either side of the centerline*/
							if(((x*x) + (y*y)) <= (radius*radius))
				    		{   /*If the point is inside the hole diameter, add it*/
								++simCount;
								++simCount;
							}	
					}
				}
				partCount = simCount;
			#endif

			break;
		}		
	}
	
	return partCount;
}

void InitSPH(SIM &svar, FLUID &fvar, CROSS &cvar, MESH &cells, State &pn, State &pnp1)
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
	ldouble step = svar.Pstep*svar.Bstep;
	int Ny = ceil(svar.Box(1)/step);	
	int Nx = ceil(svar.Box(0)/step);
	#if (SIMDIM == 3)
		int Nz = ceil(svar.Box(2)/step);
	#endif
	
 	if(svar.Bcase == 0 || svar.Bcase == 5)
	{ /*No boundary*/

	}
	else if(svar.Bcase == 1) /*Rectangle*/
	{	
		#if SIMDIM == 3
			for(int j = 0; j <= Ny ; ++j) 
			{
				for (int k=0; k<=Nz; ++k) 
				{	/*Create Left and right boundary faces*/
					StateVecD xi(0.0,j*step,k*step);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
					xi(0) = svar.Box(0);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
				}
				
			}
			for(int i = 1; i<Nx ; ++i) 
			{
				for (int j=1; j<Ny; ++j)
				{	/*Create top and bottom boundary faces*/
					StateVecD xi(i*step,j*step,0);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
					// xi(2)= Box(2); //Top boundary (Typically omitted)
					// pn.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,0));
				}
				
			}
			for(int i= 1; i<Nx; ++i) 
			{
				for(int k = 0; k <= Nz; ++k) 
				{	/*Create far and near boundary*/
					StateVecD xi(i*step, 0, k*step);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
					xi(1) = svar.Box(1);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
				}
			}
		#else 
			for(int i = 0; i <= Ny ; ++i) {
			StateVecD xi(-svar.Start(0),i*step-svar.Start(1));
			pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
			}
	/*			// Optional lid
			for(int i = 1; i <Nx ; ++i) {
				StateVecD xi(i*stepx,Box(1));
				particles.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,Bound));	
			}
			StateVecD x(stepx-svar.Start(0),(Ny+0.5)*stepy);
			pn.emplace_back(Particle(x,v,f,rho,fvar.Boundmass,true));
			x(0) = svar.Box(0) -stepx;
			pn.emplace_back(Particle(x,v,f,rho,fvar.Boundmass,true));
	*/
			for(int i= Ny; i>0; --i) {
				StateVecD xi(svar.Box(0)-svar.Start(0),i*step-svar.Start(1));
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));	
			}
			for(int i = Nx; i > 0; --i) {
				StateVecD xi(i*step-svar.Start(0),-svar.Start(1));
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
			}
		#endif		
	}
	else if(svar.Bcase == 3 || svar.Bcase == 4)
	{	/*Jet in Crossflow*/
		#if SIMDIM == 3
			ldouble holeD = svar.Jet(1)+4*fvar.H; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			double r = 0.5*holeD;
	    	double dtheta = atan((stepb)/(r));
			for (ldouble y = 0; y >= -svar.Jet(1)-stepb; y-=stepb)			
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
		#else
			ldouble jetR = 0.5*(svar.Jet(1)+4*fvar.H); /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			/*Create a bit of the pipe downward.*/
			for (ldouble y = -stepb; y >= -svar.Jet(1)-stepb; y-=stepb)			{
				StateVecD xi(-jetR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
			}

			for (ldouble y = -stepb; y >= -svar.Jet(1)-stepb; y-=stepb)
			{
				StateVecD xi(jetR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
			}
		#endif
	}
	else if(svar.Bcase > 5) 
	{
		cerr << "Boundary case is not within the design. 0 <= Bcase <= 5." << endl;
		exit(-1);
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
			/*VLM case*/ /*NOT AVAILABLE IN 2D*/
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
			#if SIMDIM == 3
				/*Create the simulation pn*/
				for( int i=0; i< svar.xyPART(0); ++i) 
				{
					for(int j=0; j< svar.xyPART(1); ++j)
					{				
						for (int k=0; k < svar.xyPART(2); ++k )
						{
							StateVecD xi(svar.Start(0)+i*svar.Pstep,
								svar.Start(1)+j*svar.Pstep,svar.Start(2)+k*svar.Pstep);		
							pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));
						}		
					}
				}
			#else
				/*Create the simulation pn*/
				for( int i=0; i< svar.xyPART(0); ++i) 
				{
					for(int j=0; j< svar.xyPART(1); ++j)
					{				
						StateVecD xi(svar.Start(0)+i*svar.Pstep,svar.Start(1)+j*svar.Pstep);		
						pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0));		
					}
				}
			#endif

			

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
		Write_settings(svar,fvar);
		exit(-1);
	}
	// cout << "Total pn: " << svar.totPts << endl;
	
	// cout << "Refresh round pn: " << svar.nrefresh << endl;
}


#endif