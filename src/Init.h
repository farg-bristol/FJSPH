/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*********                   FOR THREE DIMENSION CODE                           *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"
#include "IO.h"
#include "Add.h"

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

		if(svar.Bcase == 0)
			partCount += svar.simPts; /*Simulation pn*/
		else if (svar.Bcase == 1)
		{
			#if (SIMDIM == 3)
				partCount = Nx*Nz + 2*Nx*Ny + 2*Ny*Nz; /*Boundary particles*/
			#else
				partCount = 2*Ny + Nx; /*Boundary particles*/
			#endif
			partCount += svar.simPts; /*Simulation pn*/
		}
		else if (svar.Bcase == 2)
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
	            Ny = ceil((svar.Jet(1)*3)/step);
	            int holeWall = 2*Ny;

				/*Simulation Points*/
				int simCount = 0;
				ldouble jetR = 2*(0.5*(svar.Jet(0)));
				for (ldouble x = -jetR; x <= jetR; x+= svar.dx)
					simCount++;

				/*Need to add the pn already present*/
				int simPts = simCount*svar.nmax + simCount*ceil(svar.Jet[1]/svar.dx);

				partCount = holeWall + simPts;
			}
			#endif
		}
		else if(svar.Bcase == 3 || svar.Bcase == 4 || svar.Bcase == 6)
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
		}
		else if (svar.Bcase == 5)
		{	
			#if SIMDIM == 3
				uint simCount = 0;
				ldouble radius = 0.5*svar.Start(0);
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
				ldouble radius = 0.5*svar.Start(0);
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
		}
		else if (svar.Bcase == 7) 
		{	/*Piston driven flow*/

			#if SIMDIM == 2
				/*Create the reservoir tank*/
				ldouble tankW = svar.Start(0);
				ldouble tankD = svar.Start(1);
				ldouble stepb = (svar.Pstep*svar.Bstep);

				uint pisCnt = ceil((tankW + 8*svar.dx-4*svar.Pstep)/stepb);
				
				/*Create the reservoir tank*/
				uint wall = ceil((tankD+4*svar.dx+6*svar.Pstep)/stepb);
				
				/*Create the tapering section*/
				ldouble theta = atan(svar.Jet(1)/(0.5*tankW-0.5*svar.Jet(0)));
				ldouble stepy = stepb*sin(theta);
				uint taper = ceil((svar.Jet(1))/stepy);

				/*Create the exit bit.*/
				uint exit = ceil(svar.Jet(1)/stepb);

				/*Simulation Particles*/
				uint simCount = ceil(tankW/svar.dx)*ceil(tankD/svar.dx);	

				partCount = pisCnt + 2*(wall+taper+exit) + simCount;
			#endif
		}	

	return partCount;
}

void InitSPH(SIM &svar, FLUID &fvar, AERO &avar, State &pn, State &pnp1)
{
	if (svar.Bcase == 3 || svar.Bcase == 4 || svar.Bcase == 6)
		cout << "Initialising simulation..." << endl;
	else if (svar.Bcase == 5 || svar.Bcase == 7)
		cout << "Initialising simulation with " << svar.finPts << " points" << endl;
	else
		cout << "Initialising simulation with " << svar.simPts << " points" << endl;	
	
	
	//Structure input initialiser
	//Particle(StateVecD x, StateVecD v, StateVecD vh, StateVecD f, float rho, float Rrho, bool bound) :
	StateVecD v = StateVecD::Zero();   
	// ldouble rho=fvar.rho0;
	// ldouble press = fvar.B*(pow(rho/fvar.rho0,fvar.gam)-1);  
	
	ldouble press = fvar.pPress;
	ldouble rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

/************** Create the boundary pn  *****************/ 	 
	ldouble step = svar.Pstep*svar.Bstep;
	int Ny = ceil((svar.Box(1)+4*svar.Pstep)/step);	
	int Nx = ceil((svar.Box(0)+4*svar.Pstep)/step);
	#if (SIMDIM == 3)
		int Nz = ceil(svar.Box(2)/step);
	#endif
	uint pID = 0;

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
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0,pID));
					pID++;
					xi(0) = svar.Box(0);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0,pID));
					pID++;
				}
				
			}
			for(int i = 1; i<Nx ; ++i) 
			{
				for (int j=1; j<Ny; ++j)
				{	/*Create top and bottom boundary faces*/
					StateVecD xi(i*step,j*step,0);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,0,pID));
					pID++;
					// xi(2)= Box(2); //Top boundary (Typically omitted)
					// pn.emplace_back(Particle(xi,v,f,rho,Rrho,Boundmass,0));
				}
				
			}
			for(int i= 1; i<Nx; ++i) 
			{
				for(int k = 0; k <= Nz; ++k) 
				{	/*Create far and near boundary*/
					StateVecD xi(i*step, 0, k*step);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
					xi(1) = svar.Box(1);
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}
			}
		#else /*SIMDIM == 2*/
			for(int i = 0; i <= Ny ; ++i) 
			{	/* Left Wall*/
				StateVecD xi(-svar.Start(0)-2*svar.Pstep,i*step-svar.Start(1)-2*svar.Pstep);
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
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
			for(int i= Ny; i>0; --i) 
			{	/*Right Wall*/
				StateVecD xi(Nx*step-svar.Start(0)-2*svar.Pstep,i*step-svar.Start(1)-2*svar.Pstep);
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));	
				pID++;
			}
			for(int i = Nx; i > 0; --i) 
			{	/*Floor*/
				StateVecD xi(i*step-svar.Start(0)-2*svar.Pstep,-svar.Start(1)-2*svar.Pstep);
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}
		#endif		
	}
	else if (svar.Bcase == 2)
	{	/*Converging nozzle for jet spray*/
		#if SIMDIM == 3
			ldouble holeD = svar.Jet(0)+8*svar.dx; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			

			/*Create a bit of the pipe downward.*/
			double r = 0.5*holeD;
	    	double dtheta = atan((stepb)/(r));
			for (ldouble y = 0; y > -svar.Jet(1); y-=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}

			/*Interpolate between the big and small diameters*/
			for (ldouble y = -svar.Jet(1)*1; y > -svar.Jet(1)*2; y-=stepb)			
			{	
				ldouble x = holeD + (y-2*svar.Jet(1))*(r-holeD)/(-svar.Jet(1));
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(x*sin(theta), y, x*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}

			/*Create the wide bit of the nozzle*/
			r = holeD;
	    	dtheta = atan((stepb)/(r));
			for (ldouble y = -svar.Jet(1)*2; y >= -svar.Jet(1)*3; y-=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}

		#else

			ldouble stepb = (svar.Pstep*svar.Bstep);
			ldouble jetR = 0.5*(svar.Jet(0)+4*svar.dx); /*Diameter of hole (or width)*/
			ldouble resR = jetR*2.0;

			/*Create the wide bit of the nozzle*/
			for (ldouble y = -svar.Jet(1)*2; y >= -svar.Jet(1)*3; y-=stepb)
			{
				StateVecD xi(-resR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}	

			/*Create the tapering section*/
			for (ldouble y = -svar.Jet(1); y > -svar.Jet(1)*2; y-=stepb)
			{
				/*Interpolate between resR and jetR*/
				ldouble x = resR + (y-2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				StateVecD xi(-x,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (ldouble y = 0; y > -svar.Jet(1); y-=stepb)
			{
				StateVecD xi(-jetR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the wide bit of the nozzle*/
			for (ldouble y = -svar.Jet(1)*2; y >= -svar.Jet(1)*3; y-=stepb)
			{
				StateVecD xi(resR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the tapering section*/
			for (ldouble y = -svar.Jet(1); y > -svar.Jet(1)*2; y-=stepb)
			{
				/*Interpolate between resR and jetR*/
				ldouble x = resR + (y-2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				StateVecD xi(x,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (ldouble y = 0; y > -svar.Jet(1); y-=stepb)
			{
				StateVecD xi(jetR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}
		#endif
	}
	else if(svar.Bcase == 3 || svar.Bcase == 4 || svar.Bcase == 6)
	{	/*Jet in Crossflow*/
		#if SIMDIM == 3
			ldouble holeD = svar.Jet(0)+8*svar.dx; /*Diameter of hole (or width)*/
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
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}
		#else
			ldouble jetR = 0.5*(svar.Jet(0)+4*svar.dx); /*Radius of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			for (ldouble y = -stepb; y >= -svar.Jet(1)-stepb; y-=stepb)			{
				StateVecD xi(-jetR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			for (ldouble y = -stepb; y >= -svar.Jet(1)-stepb; y-=stepb)
			{
				StateVecD xi(jetR,y);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}
		#endif
	}
	else if (svar.Bcase == 7) 
	{	
		#if SIMDIM == 3
			/*Piston driven flow*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			ldouble tankW = svar.Start(0)+8*svar.dx;
			ldouble tankD = svar.Start(1);

			ldouble jetR = 0.5*(svar.Jet(0)+8*svar.dx); /*Diameter of hole (or width)*/
			ldouble resR = 0.5*tankW;
			v = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
			

			/*Create the piston*/
			for (ldouble z = -resR; z <= resR; z+= svar.dx)
			{ /*Do the centerline of points*/
				ldouble y = -(2*svar.Jet(1)+tankD+fvar.H);
				StateVecD xi(0.0,y,z);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,BOUND,pID));
				++pID;
			}

			for (ldouble x = svar.dx; x < resR ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				for (ldouble z = -resR; z <= resR; z+= svar.dx)
				{
					if(((x*x)/(resR*resR) + (z*z)/(resR*resR)) <= 1.0 )
		    		{   /*If the point is inside the hole diameter, add it*/
		    			ldouble y = -(2*svar.Jet(1)+tankD+fvar.H);
						StateVecD temp(x,y,z);
						StateVecD xi = svar.Rotate*temp;
						xi+= svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,BOUND,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;

						temp(0) = -x;
						xi = svar.Rotate*temp;
						xi+= svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,BOUND,pID));
						++pID;
					}
				}
			}


			/*Create the reservior*/
	    	ldouble dtheta = atan((stepb)/(resR));
			for (ldouble y = (-svar.Jet(1)*2-tankW); y <= -svar.Jet(1)*2; y+=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(resR*sin(theta), y, resR*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}

			/*Interpolate between the big and small diameters*/
			for (ldouble y = -svar.Jet(1)*2; y < -svar.Jet(1)*1; y+=stepb)			
			{	
				ldouble x = resR + (y-2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				dtheta = atan((stepb)/(x));
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(x*sin(theta), y, x*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}

			/*Create the exit section.*/
			dtheta = atan((stepb)/(jetR));
			for (ldouble y = -svar.Jet(1); y < 0; y+=stepb)			
			{	
				for(ldouble theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(jetR*sin(theta), y, jetR*cos(theta));
					/*Apply Rotation...*/
					xi = svar.Rotate*xi;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
					pID++;
				}	
			}
		#endif

		#if SIMDIM == 2
			/*Piston driven flow*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			ldouble tankW = svar.Start(0)+8*svar.dx;
			ldouble tankD = svar.Start(1)+4*svar.dx;

			ldouble jetR = 0.5*(svar.Jet(0)+4*svar.dx); /*Diameter of hole (or width)*/
			ldouble resR = 0.5*tankW;

			v = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
			/*Create the piston*/
			for(ldouble x = -(resR-fvar.H); x < (resR -fvar.H); x+=stepb)
			{
				StateVecD xi(x,-(2*svar.Jet(1)+tankD+3*fvar.H));
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
				svar.psnPts = pID;
			}

			v=StateVecD::Zero();
			/*Create the reservoir tank*/
			for (ldouble y = -(svar.Jet(1)*2+tankD+3*fvar.H); y <= -svar.Jet(1)*2; y+=stepb)
			{
				StateVecD xi(-resR,y);
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			ldouble theta = atan(svar.Jet(1)/(resR-jetR));
			ldouble stepy = stepb*sin(theta);
			/*Create the tapering section*/
			for (ldouble y = -svar.Jet(1)*2; y < -svar.Jet(1); y+=stepy)
			{
				/*Interpolate between resR and jetR*/
				ldouble x = resR - (y+2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				StateVecD xi(-x,y);
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (ldouble y = -svar.Jet(1); y < 0; y+=stepb)
			{
				StateVecD xi(-jetR,y);
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (ldouble y = 0; y > -svar.Jet(1); y-=stepb)
			{
				StateVecD xi(jetR,y);
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the tapering section*/
			for (ldouble y = -svar.Jet(1); y > -svar.Jet(1)*2; y-=stepy)
			{
				/*Interpolate between resR and jetR*/
				ldouble x = resR - (y+2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				StateVecD xi(x,y);
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}

			/*Create the other wall.*/
			for (ldouble y = -svar.Jet(1)*2; y >= -(svar.Jet(1)*2+tankD+3*fvar.H); y-=stepb)
			{
				StateVecD xi(resR,y);
				xi = svar.Rotate*xi;
				pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,BOUND,pID));
				pID++;
			}
		#endif


	}
	else if(svar.Bcase > 7) 
	{
		cerr << "Boundary case is not within the design. 0 <= Bcase <= 7." << endl;
		exit(-1);
	}
	
	svar.bndPts = pn.size();
	
/***********  Create the simulation pn  **************/
	if(svar.Bcase == 2)
	{
		/*Crossflow case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		/*Update n+1 before adding sim pn*/
		for (auto p: pn)
			pnp1.emplace_back(p);

		for (ldouble y = -svar.Jet[1]*2; y > -svar.Jet[1]*3; y-=svar.dx)
		{
			// cout << "In add points for-loop" << endl;
			AddPoints(y, svar, fvar, avar, pn, pnp1);
		}
		svar.clear = -svar.Jet[1] + 4*svar.dx;
	}
	else if (svar.Bcase == 3 || svar.Bcase == 4 || svar.Bcase == 6)
	{
		/*Crossflow case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		/*Update n+1 before adding sim pn*/
		for (auto p: pn)
			pnp1.emplace_back(p);

		for (ldouble y = 0.0; y > -svar.Jet[1]; y-=svar.dx)
		{
			// cout << "In add points for-loop" << endl;
			AddPoints(y, svar, fvar, avar, pn, pnp1);
		}

		for(size_t ii = svar.totPts-svar.nrefresh; ii < svar.totPts; ++ii)
		{	/*Fill the vector of the last particles*/
			svar.back.emplace_back(ii);
		}

		svar.clear = -svar.Jet[1] + 4*svar.dx;
	}
	else if(svar.Bcase == 5)
	{
		/*Droplet case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		/*Update n+1 before adding sim pn*/
		for (auto p: pn)
			pnp1.emplace_back(p);

		CreateDroplet(svar,fvar,pn,pnp1);
	}
	else if (svar.Bcase == 7)
	{
		svar.simPts = 0;
		#if SIMDIM == 3
			/*Create fluid in the reservoir*/
			ldouble tankW = svar.Start(0);
			ldouble tankD = svar.Start(1);

			press = fvar.pPress;
			rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

			/*Create the simulation pn*/
			for (ldouble y = -(tankD + 2*svar.Jet(1)); y < 2*svar.Jet(1); y+=svar.dx)
			{
				for (ldouble z = -tankW; z <= tankW; z+= svar.dx)
				{ /*Do the centerline of points*/
					StateVecD xi(0.0,y,z);
					xi = svar.Rotate*xi;
					pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,PIPE,pID));
					++pID;
					++svar.simPts;
				}

				for (ldouble x = svar.dx; x < tankW ; x+=svar.dx)
				{ /*Do the either side of the centerline*/
					for (ldouble z = -tankW; z <= tankW; z+= svar.dx)
					{
						if(((x*x)/(tankW*tankW) + (z*z)/(tankW*tankW)) <= 1.0 )
			    		{   /*If the point is inside the hole diameter, add it*/
							StateVecD temp(x,y,z);
							StateVecD xi = svar.Rotate*temp;
						
							pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,PIPE,pID));
							++pID;
							++svar.simPts;
							temp(0) = -x;
							xi = svar.Rotate*temp;
							pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,PIPE,pID));
							++pID;
							++svar.simPts;
						}
					}
				}
			}

		#endif
		#if SIMDIM == 2
			/*Create fluid in the reservoir*/
			ldouble tankW = svar.Start(0);
			ldouble tankD = svar.Start(1);

			press = fvar.pPress;
			rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

			/*Create the simulation pn*/
			for(ldouble y = -(tankD + 2*svar.Jet(1)+4*svar.dx); y < -(2*svar.Jet(1)+4*svar.dx); y+=svar.dx)
			{
				for(ldouble x = -(0.5*tankW); x < 0.5*tankW; x+=svar.dx) 
				{		
					StateVecD xi(x,y);		
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,PIPE,pID));
					pn.back().cellRho = fvar.rhog;
					pn.back().cellP = 20000-fvar.gasPress;
					pID++;
					++svar.simPts;
				}
			}
		#endif

		for (auto p: pn)
			pnp1.emplace_back(p);

		svar.clear = 0.0;
	}
	else
	{	
		press =fvar.pPress;
		rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);
		#if SIMDIM == 3
			/*Create the simulation pn*/
			for( int i=0; i< svar.xyPART(0); ++i) 
			{
				for(int j=0; j< svar.xyPART(1); ++j)
				{				
					for (int k=0; k < svar.xyPART(2); ++k )
					{
						StateVecD xi(svar.Start(0)+i*svar.dx,
							svar.Start(1)+j*svar.dx,svar.Start(2)+k*svar.dx);		
						pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,FREE,pID));
						pID++;
					}		
				}
			}
		#else
			/*Create the simulation pn*/
			for( int i=0; i< svar.xyPART(0); ++i) 
			{
				for(int j=0; j< svar.xyPART(1); ++j)
				{				
					StateVecD xi(svar.Start(0)+i*svar.dx,svar.Start(1)+j*svar.dx);		
					pn.emplace_back(Particle(xi,v,rho,fvar.Boundmass,press,FREE,pID));
					pn.back().cellRho = fvar.rhog;
					pn.back().cellP = 20000-fvar.gasPress;
					pID++;
				}
			}
		#endif

		for (auto p: pn)
			pnp1.emplace_back(p);
	}
	
		
	// svar.simPts+=10*10;

	svar.totPts = pn.size();
	// cout << "Boundary pn: " << svar.bndPts << endl;
	// cout << "Sim pn:" << svar.simPts << endl;
 // 	cout << "Total Points: " << svar.totPts << endl;
#ifdef DEBUG
	dbout << "Initialising simulation complete." << endl;
	dbout<< "Total Points: " << svar.totPts << endl;
	dbout << "Boundary Points: " << svar.bndPts << endl;
	dbout << "Simulation Points: " << svar.simPts << endl;
#endif
	
	if(svar.totPts!=svar.bndPts+svar.simPts)
	{
		cerr<< "Mismatch of particle count." << endl;
		cerr<< "Particle array size doesn't match defined values." << endl;
		cerr<< "Total Points: " << svar.totPts << "  Boundary Points: " << svar.bndPts
			<< "  Simulation Points: " << svar.simPts << endl;
		exit(-1);
	}
	// cout << "Total pn: " << svar.totPts << endl;
	
	// cout << "Refresh round pn: " << svar.nrefresh << endl;
}

#endif