/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef CROSS_H
#define CROSS_H

#include "Var.h"
#include "IOFunctions.h"
#include <random>
#include <stdint.h>
#include <time.h>

using std::cout;
using std::endl;

StateVecD getVelocity(StateVecD const& vavg, real const& zsq, real const& r)
{
	// return 1.092*vavg*(1.25-(zsq)/(r*r));
	return vavg;
}

void AddPoints(const real y, SIM &svar, const FLUID &fvar, const AERO &avar, State &pn, State &pnp1, size_t const& pState)
{	
	// cout << "Adding points..." << endl;
	uint pID = svar.totPts;
	
	svar.nrefresh = 0;	
	real jetR = 0.5*(svar.Jet(0));
	real r = jetR + 2*svar.dx;
	real resR = 2*jetR;
	
	StateVecD vavg;
	if(svar.Bcase == 2)
	{
		vavg = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
		jetR *= 2;
	}
	else
	{
		vavg = avar.vJet;  /*Jet velocity*/
	}
	

	/*Squeeze particles together to emulate increased pressure*/
	real press =fvar.pPress;
	real rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

	#if SIMDIM == 3
		/*Create the simulation particles*/
		StateVecD xi(0.0,y,0.0);
		xi = svar.Rotate*xi;
		xi += svar.Start;
		StateVecD v = getVelocity(vavg,0.0,r);
		pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
		pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
		++pID;
		++svar.simPts;
		++svar.nrefresh;

		for (real z = svar.dx; z <= jetR; z+= svar.dx)
		{ /*Do the centerline of points*/
			StateVecD xi(0.0,y,z);
			xi = svar.Rotate*xi;
			xi += svar.Start;
			StateVecD v = getVelocity(vavg,z*z,r);
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
		}

		for (real z = -svar.dx; z >= -jetR; z-= svar.dx)
		{ /*Do the centerline of points*/
			StateVecD xi(0.0,y,z);
			xi = svar.Rotate*xi;
			xi += svar.Start;
			StateVecD v = getVelocity(vavg,z*z,r);
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
		}

		for (real x = svar.dx; x <= jetR; x+= svar.dx)
		{ /*Do the centerline of points*/
			StateVecD xi(x,y,0.0);
			xi = svar.Rotate*xi;
			xi += svar.Start;
			StateVecD v = getVelocity(vavg,x*x,r);
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
		}

		for (real x = -svar.dx; x >= -jetR; x-= svar.dx)
		{ /*Do the centerline of points*/
			StateVecD xi(x,y,0.0);
			xi = svar.Rotate*xi;
			xi += svar.Start;
			StateVecD v = getVelocity(vavg,x*x,r);
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
		}

		for (real x = svar.dx; x <= jetR ; x+=svar.dx)
		{ /*Do the either side of the centerline*/
			for (real z = svar.dx; z <= jetR; z+= svar.dx)
			{
				if( (x*x + z*z)/(jetR*jetR) <= 1.0 )
	    		{   /*If the point is inside the hole diameter, add it*/
					StateVecD temp(x,y,z);
					StateVecD xi = svar.Rotate*temp;
					xi+= svar.Start;
					StateVecD v = getVelocity(vavg,(z*z+x*x),r);
					pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					++pID;

					temp(0) = -x;
					xi = svar.Rotate*temp;
					xi+= svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					++pID;

					temp(2) = -z;
					xi = svar.Rotate*temp;
					xi+= svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					++pID;

					temp(0) = x;
					xi = svar.Rotate*temp;
					xi+= svar.Start;
					pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
					++pID;
					svar.simPts+=4;
					svar.nrefresh+=4;
				}
			}
		}

	#else 

		StateVecD xi1(0,y);
		xi1 = svar.Rotate*xi1;
		xi1 += svar.Start;
		StateVecD v1 = getVelocity(vavg,0.0,r);
		pn.emplace_back(Particle(xi1,v1,rho,fvar.simM,press,pState,pID));
		pnp1.emplace_back(Particle(xi1,v1,rho,fvar.simM,press,pState,pID));
		++pID;
		++svar.simPts;
		++svar.nrefresh;

		/*Create the simulation particles*/
		for (real x = svar.dx; x <= jetR; x+= svar.dx)
		{ /*Do the centerline of points*/
			StateVecD xi(x,y);
			xi = svar.Rotate*xi;
			xi += svar.Start;
			StateVecD v = getVelocity(vavg,x*x,r);
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;

			xi = StateVecD(-x,y);
			xi = svar.Rotate*xi;
			xi += svar.Start;
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,pState,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
		}
	#endif

	svar.totPts += svar.nrefresh;
	++svar.addcount;
	// cout << "New points: " << svar.nrefresh << "  totPts: " <<
	// svar.totPts << " simPts: "<< svar.simPts <<  endl;
}

void CreateDroplet(SIM &svar, const FLUID &fvar, State &pn, State &pnp1)
{
	uint pID = svar.totPts;
	StateVecD v = StateVecD::Zero();
	real rho = fvar.simM/pow(svar.dx,SIMDIM);
	// real rho = fvar.rho0;
	real press = fvar.pPress;
	// real press = 0.0;
	svar.nrefresh = 0;	
	real radius = 0.500001*svar.Jet(0);

#if SIMDIM == 3
		
		for (real y = 0; y <= radius; y+=svar.dx)
		{	
			real xradius = sqrt(radius*radius - y*y);

			StateVecD xi_y(0.0,y,0.0);
			xi_y = svar.Rotate*xi_y;
			xi_y += svar.Start;
			pn.emplace_back(Particle(xi_y,v,rho,fvar.simM,press,PartState.FREE_,pID));
			pnp1.emplace_back(Particle(xi_y,v,rho,fvar.simM,press,PartState.FREE_,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;

			for (real z = svar.dx; z <= xradius; z+= svar.dx)
			{ /*Do the centerline of points*/
				StateVecD xi(0.0,y,z);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;

				
				xi = StateVecD(0.0,y,-z);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;

			}

			for (real x = svar.dx; x <= xradius ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				
				StateVecD xi_z(x,y,0.0);
				xi_z = svar.Rotate*xi_z;
				xi_z += svar.Start;
				pn.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;

				
				xi_z = StateVecD(-x,y,0.0);
				xi_z = svar.Rotate*xi_z;
				xi_z += svar.Start;
				pn.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;


				for (real z = svar.dx; z <= xradius; z+= svar.dx)
				{
					if(((x*x) + (z*z) + (y*y)) <= (radius*radius) )
		    		{   /*If the point is inside the hole diameter, add it*/
						StateVecD xi(x,y,z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;
						
						xi = StateVecD(-x,y,z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;

						xi = StateVecD(-x,y,-z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;

						xi = StateVecD(x,y,-z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;
					}
				}	
			}
		}

		for (real y = -svar.dx; y >= -radius; y-=svar.dx)
		{	
			real xradius = sqrt(radius*radius - y*y);

			StateVecD xi_y(0.0,y,0.0);
			xi_y = svar.Rotate*xi_y;
			xi_y += svar.Start;
			pn.emplace_back(Particle(xi_y,v,rho,fvar.simM,press,PartState.FREE_,pID));
			pnp1.emplace_back(Particle(xi_y,v,rho,fvar.simM,press,PartState.FREE_,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;

			for (real z = svar.dx; z <= xradius; z+= svar.dx)
			{ /*Do the centerline of points*/
				StateVecD xi(0.0,y,z);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;

				xi = StateVecD(0.0,y,-z);
				xi = svar.Rotate*xi;
				xi += svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;

			}

			for (real x = svar.dx; x <= xradius ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				
				StateVecD xi_z(x,y,0.0);
				xi_z = svar.Rotate*xi_z;
				xi_z += svar.Start;
				pn.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;

				xi_z = StateVecD(-x,y,0.0);
				xi_z = svar.Rotate*xi_z;
				xi_z += svar.Start;
				pn.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi_z,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;


				for (real z = svar.dx; z <= xradius; z+= svar.dx)
				{
					if(((x*x) + (z*z) + (y*y)) <= (radius*radius) )
		    		{   /*If the point is inside the hole diameter, add it*/
						StateVecD xi(x,y,z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;

						xi = StateVecD(-x,y,z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;

						xi = StateVecD(-x,y,-z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;

						xi = StateVecD(x,y,-z);
						xi = svar.Rotate*xi;
						xi += svar.Start;
						pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
						++pID;
						++svar.simPts;
						++svar.nrefresh;
					}
				}	
			}
		}
#else
		/*Want to perturb points on the order of machine error */
		int const interval = 2000;
		for (real y = 0; y <= radius; y+=svar.dx)
		{	
			/*Do the centerline of points*/
			srand(time(NULL));
			StateVecD perturb((rand()%interval-interval/2) * MERROR, (rand()%interval-interval/2) * MERROR);
			
			StateVecD xi(0.0,y);
			
			xi = svar.Rotate*(xi+perturb);
			xi += svar.Start;
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
			
			for (real x = svar.dx; x <= radius ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				if(((x*x) + (y*y)) <= (radius*radius) )
	    		{   /*If the point is inside the hole diameter, add it*/
	    			srand(time(NULL));
					StateVecD perturb2((rand()%interval-interval/2) * MERROR, (rand()%interval-interval/2) * MERROR);

					StateVecD xi2(x,y);
					xi2 = svar.Rotate*(xi2+perturb2);
					xi2 += svar.Start;
					pn.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					pnp1.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					++pID;
					++svar.simPts;
					++svar.nrefresh;

					perturb2 = StateVecD((rand()%interval-interval/2) * MERROR, (rand()%interval-interval/2) * MERROR);
					xi2 = StateVecD(-x,y);
					xi2 = svar.Rotate*(xi2+perturb2);
					xi2 += svar.Start;
					pn.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					pnp1.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					++pID;
					++svar.simPts;
					++svar.nrefresh;
				}	
			}
		}

		for (real y = -svar.dx; y >= -radius; y-=svar.dx)
		{	
			/*Do the centerline of points*/
			srand(time(NULL));
			StateVecD perturb((rand()%interval-interval/2) * MERROR, (rand()%interval-interval/2) * MERROR);

			StateVecD xi(0.0,y);
			xi = svar.Rotate*(xi+perturb);
			xi += svar.Start;
			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
			
			for (real x = svar.dx; x <= radius ; x+=svar.dx)
			{ /*Do the either side of the centerline*/
				if(((x*x) + (y*y)) <= (radius*radius) )
	    		{   /*If the point is inside the hole diameter, add it*/
					srand(time(NULL));
					StateVecD perturb2((rand()%interval-interval/2) * MERROR, (rand()%interval-interval/2) * MERROR);

					StateVecD xi2(x,y);
					xi2 = svar.Rotate*(xi2+perturb2);
					xi2 += svar.Start;
					pn.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					pnp1.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					++pID;
					++svar.simPts;
					++svar.nrefresh;

					perturb2 = StateVecD((rand()%interval-interval/2) * MERROR, (rand()%interval-interval/2) * MERROR);
					xi2 = StateVecD(-x,y);
					xi2 = svar.Rotate*(xi2+perturb2);
					xi2 += svar.Start;
					pn.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					pnp1.emplace_back(Particle(xi2,v,rho,fvar.simM,press,PartState.FREE_,pID));
					++pID;
					++svar.simPts;
					++svar.nrefresh;
				}	
			}
		}
#endif

	svar.totPts += svar.nrefresh;
}

void CreateRDroplet(SIM &svar, const FLUID &fvar, State &pn, State &pnp1)
{

#if SIMDIM == 3
	// /*Create particles in a circle*/
	// for(real phi = 0.0; phi < 0.25*M_PI - 0.1*dthe_; phi += dthe_)
	// {
	// 	real z = radius*sin(phi);
	// 	real rad = sqrt(radius*radius - z*z);

	// 	/*First create outer radius*/
	// 	uint kk = 0;
	// 	for(real theta = 0.0; theta < 2*M_PI - 0.1*dthe_; theta += dthe_)
	// 	{
	// 		real x = radius*sin(theta);
	// 		real y = radius*cos(theta);
	// 		StateVecD xi(x,y);

	// 		pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 		pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 		++pID;
	// 		++svar.simPts;
	// 		++svar.nrefresh;
	// 		++kk;
	// 	}

	// 	int ii = 1;
	// 	if (nrad > 5)
	// 	{
	// 		for(real rad = radius - dx_r; rad > 0.85*radius; rad -= dx_r)
	// 		{
	// 			real dtheta = 2.0*M_PI/real(kk);

	// 			real tstart = 0.0;
	// 			if((ii+1) % 2 == 0)
	// 			{
	// 				tstart += dtheta/2;
	// 			}

	// 			kk = 0;
	// 			for(real theta = tstart; theta < 2*M_PI-0.2*dtheta; theta += dtheta)
	// 			{
	// 				real x = rad*sin(theta);
	// 				real y = rad*cos(theta);
	// 				StateVecD xi(x,y);

	// 				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 				pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 				++pID;
	// 				++svar.simPts;
	// 				++svar.nrefresh;
	// 				++kk;
	// 			}
	// 			++ii;
	// 		}
	// 	}

		
	// 	for(real rad = radius-dx_r*ii; rad > 0.99*dx_r; rad-= dx_r)
	// 	{
	// 		real dtheta = atan(svar.dx/rad);
	// 		real ncirc = floor(abs(2.0*M_PI/dtheta));
	// 		dtheta = 2.0*M_PI/(ncirc);
			

	// 		real tstart = 0.0;
	// 		if(ii % 2 == 0)
	// 		{
	// 			tstart += dtheta/2;
	// 		}

	// 		for(real theta = tstart; theta < 2*M_PI-0.2*dtheta; theta += dtheta)
	// 		{
	// 			real x = rad*sin(theta);
	// 			real y = rad*cos(theta);
	// 			StateVecD xi(x,y);

	// 			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 			++pID;
	// 			++svar.simPts;
	// 			++svar.nrefresh;
	// 		}
	// 		++ii;

	// 	}

	// 	StateVecD xi = StateVecD::Zero();
	// 	pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 	pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	// 	++pID;
	// 	++svar.simPts;
	// 	++svar.nrefresh;
	// }
#else 
	uint pID = svar.totPts;
	StateVecD v = StateVecD::Zero();
	real rho = fvar.simM/pow(svar.dx,SIMDIM);
	// real rho = fvar.rho0;
	real press = fvar.pPress;
	// real press = 0.0;
	svar.nrefresh = 0;	
	real radius = 0.5*svar.Jet(0);

	real dthe_ = atan(svar.dx/radius);
	uint ncirc = floor(abs(2.0*M_PI/dthe_));
	dthe_ = 2.0*M_PI/ncirc;

	real dx_r = radius*tan(dthe_);
	uint nrad = floor(radius/(sqrt(3)*dx_r/2));
	dx_r = radius/real(nrad);
	/*Create particles in a circle*/
	/*First create outer radius*/
	uint kk = 0;
	for(real theta = 0.0; theta < 2*M_PI - 0.1*dthe_; theta += dthe_)
	{
		real x = radius*sin(theta);
		real y = radius*cos(theta);
		StateVecD xi(x,y);

		pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
		pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
		++pID;
		++svar.simPts;
		++svar.nrefresh;
		++kk;
	}

	int ii = 1;
	if (nrad > 5)
	{
		for(real rad = radius - dx_r; rad > 0.85*radius; rad -= dx_r)
		{
			real dtheta = 2.0*M_PI/real(kk);

			real tstart = 0.0;
			if((ii+1) % 2 == 0)
			{
				tstart += dtheta/2;
			}

			kk = 0;
			for(real theta = tstart; theta < 2*M_PI-0.2*dtheta; theta += dtheta)
			{
				real x = rad*sin(theta);
				real y = rad*cos(theta);
				StateVecD xi(x,y);

				pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
				++pID;
				++svar.simPts;
				++svar.nrefresh;
				++kk;
			}
			++ii;
		}
	}

	
	for(real rad = radius-dx_r*ii; rad > 0.99*dx_r; rad-= dx_r)
	{
		real dtheta = atan(svar.dx/rad);
		real ncirc = floor(abs(2.0*M_PI/dtheta));
		dtheta = 2.0*M_PI/(ncirc);
		

		real tstart = 0.0;
		if(ii % 2 == 0)
		{
			tstart += dtheta/2;
		}

		for(real theta = tstart; theta < 2*M_PI-0.2*dtheta; theta += dtheta)
		{
			real x = rad*sin(theta);
			real y = rad*cos(theta);
			StateVecD xi(x,y);

			pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
			++pID;
			++svar.simPts;
			++svar.nrefresh;
		}
		++ii;

	}

	StateVecD xi = StateVecD::Zero();
	pn.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	pnp1.emplace_back(Particle(xi,v,rho,fvar.simM,press,PartState.FREE_,pID));
	++pID;
	++svar.simPts;
	++svar.nrefresh;

#endif

}


namespace PoissonSample
{
	class PRNG
	{
	public:
		PRNG(): gen_(std::random_device()()), dis_( 0.0, 1.0 )
		{
			// prepare PRNG
		}

		real randomReal()
		{
			return static_cast<real>( dis_( gen_ ) );
		}

		int randomInt( int maxValue )
		{
			std::uniform_int_distribution<> disInt( 0, maxValue );
			return disInt( gen_ );
		}

	private:
		std::mt19937 gen_;
		std::uniform_real_distribution<real> dis_;
	};

	bool isInCircle(const StateVecD& centre, const StateVecD& p, const real radius)
	{
		const StateVecD f = p - centre;
		return f.squaredNorm() <= radius;
	}

	StateVecI imageToGrid(const StateVecD& P, const real cellSize )
	{
		#if SIMDIM == 2
		return StateVecI( (int)(P(0)/cellSize), (int)(P(1)/cellSize));
		#else
		return StateVecI((int)(P(0)/cellSize), (int)(P(1)/cellSize), (int)(P(2)/cellSize));
		#endif
	}

	struct Grid
	{
		Grid( uint w, real minDist, real cellSize ):
			 minDist_(minDist), cellSize_( cellSize ), w_(w), h_(w)
		#if SIMDIM == 3
		, d_(w)
		#endif
		{
			#if SIMDIM == 2
				grid_ = std::vector<std::vector<StateVecD>>(w,std::vector<StateVecD>(w,StateVecD::Zero()));
			#endif
			#if SIMDIM == 3
				grid_ = std::vector<std::vector<std::vector<StateVecD>>>
				(w,std::vector<std::vector<StateVecD>>(w,std::vector<StateVecD>(w,StateVecD::Zero())));
			#endif
		}

		void insert(const StateVecD& p)
		{
			const StateVecI g = imageToGrid(p, cellSize_);
			// #pragma omp critical
			// {
			// cout << "Grid position: " << p(0) << "  " << p(1);
			// #if SIMDIM == 3
			// cout << "  " << p(2);
			// #endif 
			// cout << endl;
			// cout << "Grid index: " << g(0) << "  " << g(1);
			// #if SIMDIM == 3
			// cout << "  " << g(2);
			// #endif 
			
			// cout << endl << endl;
			// cout << grid_.size() << endl;
			if(grid_.size() == 0)
			{
				cout << "Ghost particle grid size is zero." << endl;
				exit(-1);
			}

			if(g(0) >  static_cast<int>(grid_.size()))
			{
				std::cout << "Tried to access grid_ out of bounds in i direction." << std::endl;
				std::cout << g(0) << "  " << g(1) << std::endl;
				exit(-1);
			}
			else if ( g(1) > static_cast<int>(grid_[g(0)].size()))
			{
				std::cout << "Tried to access grid_ out of bounds in j direction." << std::endl;
				exit(-1);
			}
			// }

			#if SIMDIM == 2
				grid_[g(0)][g(1)] = p;
			#endif
			#if SIMDIM == 3
				grid_[g(0)][g(1)][g(2)] = p;
			#endif
		}

		bool isInNeighbourhood(const StateVecD& point)
		{
			StateVecI g = imageToGrid(point, cellSize_);

			// number of adjucent cells to look for neighbour points
			const int D = 5;

			// scan the neighbourhood of the point in the grid
			for ( int ii = g(0) - D; ii < g(0) + D; ii++ )
			{
				for ( int jj = g.y() - D; jj < g.y() + D; jj++ )
				{	
					#if SIMDIM == 2
						if ( ii >= 0 && ii < int(w_) && jj >= 0 && jj < int(h_) )
						{
							const StateVecD P = grid_[ii][jj];

							if ( (P-point).norm() < minDist_ ) { return true; }
						}
					#endif

					#if SIMDIM == 3
						for (int kk = g.z() - D; kk < g.z() + D; kk++)
						{
							if ( ii >= 0 && ii < int(w_) && jj >= 0 && jj < int(h_)
							&& kk >= 0 && kk < int(d_) )
							{
								const StateVecD P = grid_[ii][jj][kk];

								if ( (P-point).norm() < minDist_ ) { return true; }
							}
						}
					#endif
				}
			}

			return false;
		}

	private:

		real minDist_, cellSize_;
		uint w_;
		uint h_;
		
		
		#if SIMDIM == 2
			std::vector<std::vector<StateVecD>> grid_;
		#endif

		#if SIMDIM == 3
			uint d_;
			std::vector<std::vector<std::vector<StateVecD>>> grid_;
		#endif
	};

	StateVecD popRandom(std::vector<StateVecD>& points, PRNG& generator)
	{
		const int idx = generator.randomInt( points.size()-1 );
		const StateVecD p = points[idx];
		points.erase( points.begin() + idx );
		return p;
	}

	StateVecD generateRandomPointAround( const StateVecD& p, real minDist, PRNG& generator )
	{
		// start with non-uniform distribution
		const real R1 = generator.randomReal();
		const real R2 = generator.randomReal();
		
		// radius should be between MinDist and 2 * MinDist
		const real radius = minDist * ( R1 + 1.0 );

		// random angle
		const real angle1 = 2 * M_PI * R2;

		#if SIMDIM == 3
			const real R3 = generator.randomReal();
			const real angle2 = 2 * M_PI * R3;
		#endif

		// the new point is generated around the point (x, y)
		#if SIMDIM == 2
		return StateVecD(p(0) + radius*cos(angle1), p(1) + radius*sin(angle1));
		#endif
		#if SIMDIM == 3
		return StateVecD(p(0) + radius*cos(angle1)*sin(angle2), 
						 p(1) + radius*sin(angle1)*sin(angle2),
						 p(2) + radius*cos(angle2));
		#endif
	}

	/**
		Return a vector of generated points
		sampleLimit - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
	**/
	std::vector<Part> generatePoissonPoints(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells,
	 uint const& host, State const& pnp1, outl const& outlist, StateVecD const& norm, StateVecD const& avgV)
	{
		/*Variables for the poisson disk sampling*/
		real radius = fvar.sr;
		uint sampleLimit = 30;
		PRNG generator;
		uint numPoints = svar.nfull;
		uint pID = svar.totPts;

		/*Properties for new particles*/
		StateVecD vel= avar.vInf;
		StateVecD Vdiff = StateVecD::Zero();
		real Pbase = 0.0;
		real press = 0;
		real rho = fvar.rho0;

		if(svar.Asource == 1)
		{
			vel = pnp1[host].cellV;
			Vdiff =  vel - avgV;
			Pbase = pnp1[host].cellP - avar.pRef;
		}
		else if (svar.Asource == 2)
		{
			vel = (pnp1[host].cellV+cells.cPertnp1[pnp1[host].cellID]);
			Vdiff =  vel - /*pi.v*/ avgV;
			Pbase = pnp1[host].cellP - avar.pRef;
		}
#if SIMDIM == 3
		else if(svar.Asource == 3)
		{	
			vel = svar.vortex.getVelocity(pnp1[host].xi);
			Vdiff = vel - avgV;
			Pbase = 0.5*avar.rhog*(pow(avar.vRef,2.0)-pow(vel.norm(),2.0));
		}
#endif
		else
		{
			Vdiff = vel - avgV;
			Pbase = 0.5*avar.rhog*Vdiff.squaredNorm();
		}

		real theta = acos(-norm.normalized().dot(Vdiff.normalized()));
		
		real Cp = 0.0;

		if(abs(theta) < 2.4877)
		{
			Cp = 1.0 - 2.5*pow(sin(abs(theta)),2.0);
		}
		else
		{
			Cp = 0.075;
		}


		press = Pbase + 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;
		rho = fvar.rho0 * pow((press/fvar.B + 1.0),1.0/fvar.gam);


		// const real rho = pnp1[host].cellRho;
		// const real mass = fvar.rhog* pow(svar.Pstep, SIMDIM);

		
		const real mass = pnp1[host].m;

		const real deltax = pow(mass/rho, 1.0/real(SIMDIM));
		
		// #pragma omp critical
		// cout << "CellID: " << pnp1[host].cellID << "  " << press 
		// << "  " << rho << "  " << mass << "  " << deltax << endl;

 		const real minDist = /*svar.Pstep*/ deltax;
		std::vector<Part> samplePoints;
		std::vector<StateVecD> processList;
		std::vector<Part> airP;

		// create the grid
		#if SIMDIM == 2
			const StateVecD origin(pnp1[host].xi(0)-2*fvar.H,pnp1[host].xi(1)-2*fvar.H);
		#endif
		#if SIMDIM == 3
			const StateVecD origin(pnp1[host].xi(0)-2*fvar.H,pnp1[host].xi(1)-2*fvar.H,
									pnp1[host].xi(2)-2*fvar.H);
		#endif

		const real cellSize = minDist / sqrt(2.0);
		const uint gridW = (uint)ceil(4*fvar.H / cellSize);
		// cout << "GridW: " << gridW << " minDist: " << minDist << " cellSize: " << cellSize << endl; 
		Grid grid(gridW, minDist, cellSize);

		/*Fill out the prexisting particles to add points around*/
		for(auto jj:outlist[host])
		{
			samplePoints.push_back(pnp1[jj.first]);
			grid.insert(pnp1[jj.first].xi-origin);
		}

		/*Try and add a point where it won't conflict*/
		#if SIMDIM == 2
			const StateVecD circCent(2*fvar.H, 2*fvar.H);
		#endif
		#if SIMDIM == 3
			const StateVecD circCent(2*fvar.H, 2*fvar.H, 2*fvar.H);
		#endif
		
		StateVecD firstPoint;
	 	
		do {/*Generate a random number between 0 and 1, then normalise it to the grid.*/
			#if SIMDIM == 2
			firstPoint = StateVecD(generator.randomReal()*4*fvar.H, generator.randomReal()*4*fvar.H);
			#endif
			#if SIMDIM == 3
			firstPoint = StateVecD(generator.randomReal()*4*fvar.H, generator.randomReal()*4*fvar.H,
						generator.randomReal()*4*fvar.H);
			#endif	

		} while (isInCircle(circCent,firstPoint,radius) != 1 || grid.isInNeighbourhood(firstPoint) != 0);
		
		// update containers
		processList.push_back(firstPoint);
		
		grid.insert(firstPoint);

		// generate new points for each point in the queue
		while (!processList.empty() && samplePoints.size() < floor(float(0.95*numPoints)))
		{
			const StateVecD point = popRandom(processList, generator);

			for (uint ii = 0; ii < sampleLimit; ii++)
			{
				const StateVecD newPoint = generateRandomPointAround(point, minDist, generator);
				const bool canFitPoint = isInCircle(circCent,newPoint,radius);

				if (canFitPoint == true && grid.isInNeighbourhood(newPoint) == false)
				{
					processList.push_back(newPoint);
					// StateVecD Point = newPoint+origin;
					samplePoints.push_back(Part(newPoint+origin,StateVecD::Zero(),0.0,0.0,0.0,PartState.GHOST_,0));
					airP.emplace_back(Part(newPoint+origin,vel,press,rho,mass,PartState.GHOST_,pID));
					pID++;
					grid.insert(newPoint);
					// cout << "New Point: " << point(0) << " " << point(1) << endl;
					continue;
				}
			}
		}
		return airP;
	}
}

#endif
