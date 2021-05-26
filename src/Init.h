/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*********                   FOR THREE DIMENSION CODE                           *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"
#include "IO.h"
#include "Add.h"

void InitSPH(SIM& svar, FLUID const& fvar, AERO const& avar, State& pn, State& pnp1)
{
	if (svar.Bcase == 3 || svar.Bcase == 4 || svar.Bcase == 6 || svar.Bcase == 7)
		cout << "Initialising simulation..." << endl;
	else if (svar.Bcase == 5 || svar.Bcase == 8)
		cout << "Initialising simulation with " << svar.finPts << " points" << endl;
	
	
	
	//Structure input initialiser
	//Particle(StateVecD x, StateVecD v, StateVecD vh, StateVecD f, float rho, float Rrho, bool bound) :
	StateVecD v = StateVecD::Zero();   
	// real rho=fvar.rho0;
	// real press = fvar.B*(pow(rho/fvar.rho0,fvar.gam)-1);  
	
	real press = fvar.pPress;
	real rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

/************** Create the boundary pn  *****************/ 	 
	real step = svar.Pstep*svar.Bstep;
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
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
					xi(0) = svar.Box(0);
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}
				
			}
			for(int i = 1; i<Nx ; ++i) 
			{
				for (int j=1; j<Ny; ++j)
				{	/*Create top and bottom boundary faces*/
					StateVecD xi(i*step,j*step,0);
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
					// xi(2)= Box(2); //Top boundary (Typically omitted)
					// pn.emplace_back(Particle(xi,v,f,rho,Rrho,bndM,0));
				}
				
			}
			for(int i= 1; i<Nx; ++i) 
			{
				for(int k = 0; k <= Nz; ++k) 
				{	/*Create far and near boundary*/
					StateVecD xi(i*step, 0, k*step);
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
					xi(1) = svar.Box(1);
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}
			}
		#else /*SIMDIM == 2*/
			for(int i = 0; i <= Ny ; ++i) 
			{	/* Left Wall*/
				StateVecD xi(-svar.Start(0)-2*svar.Pstep,i*step-svar.Start(1)-2*svar.Pstep);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}
	/*			// Optional lid
			for(int i = 1; i <Nx ; ++i) {
				StateVecD xi(i*stepx,Box(1));
				particles.emplace_back(Particle(xi,v,f,rho,Rrho,bndM,Bound));	
			}
			StateVecD x(stepx-svar.Start(0),(Ny+0.5)*stepy);
			pn.emplace_back(Particle(x,v,f,rho,fvar.bndM,true));
			x(0) = svar.Box(0) -stepx;
			pn.emplace_back(Particle(x,v,f,rho,fvar.bndM,true));
	*/
			for(int i= Ny; i>0; --i) 
			{	/*Right Wall*/
				StateVecD xi(Nx*step-svar.Start(0)-2*svar.Pstep,i*step-svar.Start(1)-2*svar.Pstep);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));	
				pID++;
			}
			for(int i = Nx; i > 0; --i) 
			{	/*Floor*/
				StateVecD xi(i*step-svar.Start(0)-2*svar.Pstep,-svar.Start(1)-2*svar.Pstep);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}
		#endif		
	}
	else if (svar.Bcase == 2)
	{	/*Converging nozzle for jet spray*/
		#if SIMDIM == 3
			real holeD = svar.Jet(0)+4*svar.dx; /*Diameter of hole (or width)*/
			real stepb = (svar.Pstep*svar.Bstep);
			

			/*Create a bit of the pipe downward.*/
			real r = 0.5*holeD;
	    	real dtheta = atan((stepb)/(r));
			for (real y = 0; y > -svar.Jet(1); y-=stepb)			
			{	
				for(real theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}	
			}

			/*Interpolate between the big and small diameters*/
			for (real y = -svar.Jet(1)*1; y > -svar.Jet(1)*2; y-=stepb)			
			{	
				real x = holeD + (y-2*svar.Jet(1))*(r-holeD)/(-svar.Jet(1));
				for(real theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(x*sin(theta), y, x*cos(theta));
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}	
			}

			/*Create the wide bit of the nozzle*/
			r = holeD;
	    	dtheta = atan((stepb)/(r));
			for (real y = -svar.Jet(1)*2; y >= -svar.Jet(1)*3; y-=stepb)			
			{	
				for(real theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}	
			}

		#else
			real stepb = (svar.Pstep*svar.Bstep);
			real jetR = 0.5*(svar.Jet(0)+4*svar.dx); /*Diameter of hole (or width)*/
			real resR = jetR*2.0;

			/*Create the wide bit of the nozzle*/
			for (real y = -svar.Jet(1)*2; y >= -svar.Jet(1)*3; y-=stepb)
			{
				StateVecD xi(-resR,y);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}	

			/*Create the tapering section*/
			for (real y = -svar.Jet(1); y > -svar.Jet(1)*2; y-=stepb)
			{
				/*Interpolate between resR and jetR*/
				real x = resR + (y-2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				StateVecD xi(-x,y);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (real y = 0; y > -svar.Jet(1); y-=stepb)
			{
				StateVecD xi(-jetR,y);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the wide bit of the nozzle*/
			for (real y = -svar.Jet(1)*2; y >= -svar.Jet(1)*3; y-=stepb)
			{
				StateVecD xi(resR,y);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the tapering section*/
			for (real y = -svar.Jet(1); y > -svar.Jet(1)*2; y-=stepb)
			{
				/*Interpolate between resR and jetR*/
				real x = resR + (y-2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
				StateVecD xi(x,y);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (real y = 0; y > -svar.Jet(1); y-=stepb)
			{
				StateVecD xi(jetR,y);
				pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}
		#endif
	}
	else if(svar.Bcase == 3)
	{	/*Jet*/
		int const interval = 20000;
		#if SIMDIM == 3
			real holeD = svar.Jet(0)+7*svar.dx; /*Diameter of hole (or width)*/
			real stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			real r = 0.5*holeD;
	    	real dtheta = atan((stepb)/(r));

	    	int ncirc = floor(2.0*M_PI/dtheta);
	    	dtheta = 2.0*M_PI/real(ncirc);

			for(real ii = 0; ii < 4; ii+=1.0)
			{
				r = 0.5*holeD + ii*stepb;
				for (real y = 0; y >= -svar.Jet(1)-stepb; y-=stepb)			
				{	
					for(real theta = 0; theta < 2*M_PI; theta += dtheta)
					{
						StateVecD perturb(random(interval), random(interval), random(interval));
						StateVecD xi(r*sin(theta), y, r*cos(theta));
						xi += perturb;
						pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
						pID++;
					}	
				}
			}
		#else

			real jetR = 0.5*(svar.Jet(0)+5.5*svar.dx); /*Radius of hole (or width)*/
			real stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			for(real ii = 0; ii < 4; ii+=1.0)
			{
				real x = jetR + ii*stepb;
				for (real y = 0.0; y >= -svar.Jet(1) - 2.0 * stepb; y-=stepb)			
				{
					StateVecD perturb(random(interval), random(interval));
					StateVecD xi(-x,y);
					xi += perturb;
					xi += svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}

				for (real y = 0.0; y >= -svar.Jet(1) - 2.0 * stepb; y -= stepb)
				{
					StateVecD perturb(random(interval), random(interval));
					StateVecD xi(x,y);
					xi += perturb;
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}
			}
		#endif
	}
	else if (svar.Bcase == 5) 
	{	
		// #if SIMDIM == 3
		// 	/*Piston driven flow*/
		// 	real stepb = (svar.Pstep*svar.Bstep);
		// 	real tankW = svar.Start(0)+8*svar.dx;
		// 	real tankD = svar.Start(1);

		// 	real jetR = 0.5*(svar.Jet(0)+8*svar.dx); /*Diameter of hole (or width)*/
		// 	real resR = 0.5*tankW;
		// 	v = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
			

		// 	/*Create the piston*/
		// 	for (real z = -resR; z <= resR; z+= svar.dx)
		// 	{ /*Do the centerline of points*/
		// 		real y = -(2*svar.Jet(1)+tankD+fvar.H);
		// 		StateVecD xi(0.0,y,z);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.BOUND_,pID));
		// 		++pID;
		// 	}

		// 	for (real x = svar.dx; x < resR ; x+=svar.dx)
		// 	{ /*Do the either side of the centerline*/
		// 		for (real z = -resR; z <= resR; z+= svar.dx)
		// 		{
		// 			if(((x*x)/(resR*resR) + (z*z)/(resR*resR)) <= 1.0 )
		//     		{   /*If the point is inside the hole diameter, add it*/
		//     			real y = -(2*svar.Jet(1)+tankD+fvar.H);
		// 				StateVecD xi(x,y,z);
		// 				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.BOUND_,pID));
		// 				++pID;
		// 				++svar.simPts;
		// 				++svar.nrefresh;

		// 				xi(0) = -x;
		// 				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.BOUND_,pID));
		// 				++pID;
		// 			}
		// 		}
		// 	}


		// 	/*Create the reservior*/
	 //    	real dtheta = atan((stepb)/(resR));
		// 	for (real y = (-svar.Jet(1)*2-tankW); y <= -svar.Jet(1)*2; y+=stepb)			
		// 	{	
		// 		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		// 		{
		// 			StateVecD xi(resR*sin(theta), y, resR*cos(theta));
		// 			/*Apply Rotation...*/
		// 			pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 			pID++;
		// 		}	
		// 	}

		// 	/*Interpolate between the big and small diameters*/
		// 	for (real y = -svar.Jet(1)*2; y < -svar.Jet(1)*1; y+=stepb)			
		// 	{	
		// 		real x = resR + (y-2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
		// 		dtheta = atan((stepb)/(x));
		// 		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		// 		{
		// 			StateVecD xi(x*sin(theta), y, x*cos(theta));
		// 			/*Apply Rotation...*/
		// 			pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 			pID++;
		// 		}	
		// 	}

		// 	/*Create the exit section.*/
		// 	dtheta = atan((stepb)/(jetR));
		// 	for (real y = -svar.Jet(1); y < 0; y+=stepb)			
		// 	{	
		// 		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		// 		{
		// 			StateVecD xi(jetR*sin(theta), y, jetR*cos(theta));
		// 			/*Apply Rotation...*/
		// 			pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 			pID++;
		// 		}	
		// 	}
		// #endif

		// #if SIMDIM == 2
		// 	/*Piston driven flow*/
		// 	real stepb = (svar.Pstep*svar.Bstep);
		// 	real tankW = svar.Start(0)+8*svar.dx;
		// 	real tankD = svar.Start(1)+4*svar.dx;

		// 	real jetR = 0.5*(svar.Jet(0)+4*svar.dx); /*Diameter of hole (or width)*/
		// 	real resR = 0.5*tankW;

		// 	v = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
		// 	/*Create the piston*/
		// 	for(real x = -(resR-fvar.H); x < (resR -fvar.H); x+=stepb)
		// 	{
		// 		StateVecD xi(x,-(2*svar.Jet(1)+tankD+3*fvar.H));
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 		svar.psnPts = pID;
		// 	}

		// 	v=StateVecD::Zero();
		// 	/*Create the reservoir tank*/
		// 	for (real y = -(svar.Jet(1)*2+tankD+3*fvar.H); y <= -svar.Jet(1)*2; y+=stepb)
		// 	{
		// 		StateVecD xi(-resR,y);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	real theta = atan(svar.Jet(1)/(resR-jetR));
		// 	real stepy = stepb*sin(theta);
		// 	/*Create the tapering section*/
		// 	for (real y = -svar.Jet(1)*2; y < -svar.Jet(1); y+=stepy)
		// 	{
		// 		/*Interpolate between resR and jetR*/
		// 		real x = resR - (y+2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
		// 		StateVecD xi(-x,y);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the exit bit.*/
		// 	for (real y = -svar.Jet(1); y < 0; y+=stepb)
		// 	{
		// 		StateVecD xi(-jetR,y);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the exit bit.*/
		// 	for (real y = 0; y > -svar.Jet(1); y-=stepb)
		// 	{
		// 		StateVecD xi(jetR,y);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the tapering section*/
		// 	for (real y = -svar.Jet(1); y > -svar.Jet(1)*2; y-=stepy)
		// 	{
		// 		/*Interpolate between resR and jetR*/
		// 		real x = resR - (y+2*svar.Jet(1))*(jetR-resR)/(-svar.Jet(1));
		// 		StateVecD xi(x,y);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the other wall.*/
		// 	for (real y = -svar.Jet(1)*2; y >= -(svar.Jet(1)*2+tankD+3*fvar.H); y-=stepb)
		// 	{
		// 		StateVecD xi(resR,y);
		// 		pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}
		// #endif

		// Point in cell.



	}
	else if(svar.Bcase > 5) 
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
		for (auto& pi : pn)
			pnp1.emplace_back(pi);

		for (real y = -svar.Jet[1]*2; y > -svar.Jet[1]*3; y-=svar.dx)
		{
			// cout << "In add points for-loop" << endl;
			AddPoints(y, svar, fvar, avar, pn, pnp1, PartState.PIPE_);
		}
	}
	else if (svar.Bcase == 3)
	{
		/*Jet case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		/*Update n+1 before adding sim pn*/
		for (auto& pi : pn)
			pnp1.emplace_back(pi);

		real y = -svar.Jet[1]+4*svar.dx;
		for (y = -svar.dx; y >  -svar.Jet[1]; y-= svar.dx)
		{
			#if SIMDIM == 3
				Add_Radial_Points(y, svar, fvar, avar, pn, pnp1, PartState.PIPE_);
			#else
				AddPoints(y, svar, fvar, avar, pn, pnp1, PartState.PIPE_);
			#endif
		}

		for (size_t ii = svar.totPts - svar.nrefresh; ii < svar.totPts; ++ii)
		{ /*Fill the vector of the last particles*/
			pn[ii].b = PartState.PIPE_;
			pnp1[ii].b = PartState.PIPE_;
			svar.back.emplace_back(ii);
		}

		Add_Buffer(svar,fvar,pn,pnp1);

	}
	else if(svar.Bcase == 4)
	{
		/*Droplet case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		/*Update n+1 before adding sim pn*/
		for (auto& pi : pn)
			pnp1.emplace_back(pi);

		// CreateRDroplet(svar,fvar,pn,pnp1);
		CreateDroplet(svar,fvar,pn,pnp1);
	}
	else if (svar.Bcase == 5)
	{	/*Piston driven flow*/
		svar.simPts = 0;
		// #if SIMDIM == 3
		// 	/*Create fluid in the reservoir*/
		// 	real tankW = svar.Start(0);
		// 	real tankD = svar.Start(1);

		// 	press = fvar.pPress;
		// 	rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

		// 	/*Create the simulation pn*/
		// 	for (real y = -(tankD + 2*svar.Jet(1)); y < 2*svar.Jet(1); y+=svar.dx)
		// 	{
		// 		for (real z = -tankW; z <= tankW; z+= svar.dx)
		// 		{ /*Do the centerline of points*/
		// 			StateVecD xi(0.0,y,z);
		// 			xi = svar.Rotate*xi;
		// 			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PIPE,pID));
		// 			++pID;
		// 			++svar.simPts;
		// 		}

		// 		for (real x = svar.dx; x < tankW ; x+=svar.dx)
		// 		{ /*Do the either side of the centerline*/
		// 			for (real z = -tankW; z <= tankW; z+= svar.dx)
		// 			{
		// 				if(((x*x)/(tankW*tankW) + (z*z)/(tankW*tankW)) <= 1.0 )
		// 	    		{   /*If the point is inside the hole diameter, add it*/
		// 					StateVecD temp(x,y,z);
		// 					StateVecD xi = svar.Rotate*temp;
						
		// 					pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PIPE,pID));
		// 					++pID;
		// 					++svar.simPts;
		// 					temp(0) = -x;
		// 					xi = svar.Rotate*temp;
		// 					pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PIPE,pID));
		// 					++pID;
		// 					++svar.simPts;
		// 				}
		// 			}
		// 		}
		// 	}

		// #endif
		// #if SIMDIM == 2
		// 	/*Create fluid in the reservoir*/
		// 	real tankW = svar.Start(0);
		// 	real tankD = svar.Start(1);

		// 	press = fvar.pPress;
		// 	rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

		// 	/*Create the simulation pn*/
		// 	for(real y = -(tankD + 2*svar.Jet(1)+4*svar.dx); y < -(2*svar.Jet(1)+4*svar.dx); y+=svar.dx)
		// 	{
		// 		for(real x = -(0.5*tankW); x < 0.5*tankW; x+=svar.dx) 
		// 		{		
		// 			StateVecD xi(x,y);		
		// 			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PIPE,pID));
		// 			pn.back().cellRho = fvar.rhog;
		// 			pn.back().cellP = 20000-fvar.gasPress;
		// 			pID++;
		// 			++svar.simPts;
		// 		}
		// 	}
		// #endif


		// StateVecD xi = StateVecD::Zero();

		// xi(0)-=0.4;
		// xi(1)+=0.1;


		// pn.emplace_back(Particle(xi,StateVecD::Zero(),fvar.rho0,fvar.simM,0.0,PartState.FREE_,pID));

		// for (auto p: pn)
		// 	pnp1.emplace_back(p);
		// svar.simPts++;


		CreateDroplet(svar,fvar,pn,pnp1);
		// Perturb so that the points are at numerical zero. 
		// real pturb = 1e-1;
		for (auto& pi : pnp1)
		{
// #if SIMDIM == 3
// 			pi.xi += StateVecD(-pturb,-pturb,-pturb);
// #else
// 			pi.xi += StateVecD(-pturb,-pturb);
// #endif
			pi.cellID = 0;
			pi.cellV = avar.vInf;
			pi.cellP = avar.pRef;
			// pi.cellRho = fvar.rho0 * pow((fvar.gasPress/fvar.B + 1),1/fvar.gam);
		}

		for (auto& pi : pn)
		{
// #if SIMDIM == 3
// 			pi.xi += StateVecD(-pturb,-pturb,-pturb);
// #else
// 			pi.xi += StateVecD(-pturb,-pturb);
// #endif
			pi.cellID = 0;
			pi.cellV = avar.vInf;
			pi.cellP = avar.pRef;
			// pi.cellRho = fvar.rho0 * pow((fvar.gasPress/fvar.B + 1),1/fvar.gam);
		}

		svar.clear = 0.0;
		
	}
	else
	{	
		press = fvar.pPress;
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
						pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.FREE_,pID));
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
					pn.emplace_back(Particle(xi,v,rho,fvar.bndM,press,PartState.FREE_,pID));
					// pn.back().cellP = 20000-fvar.gasPress;
					pID++;
				}
			}
		#endif

		for (auto& pi : pn)
			pnp1.emplace_back(pi);
	}


		
	// svar.simPts+=10*10;

	svar.totPts = pn.size();

	for(size_t ii = 0; ii < svar.totPts; ++ii)
	{
		pn[ii].xi = (svar.Rotate * pn[ii].xi) + svar.Start;
		pnp1[ii].xi = (svar.Rotate * pnp1[ii].xi) + svar.Start;
	}
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