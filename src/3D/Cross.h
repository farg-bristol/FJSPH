/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef CROSS_H
#define CROSS_H

#include "Var.h"
#include "Neighbours.h"

using namespace std;

void AddPoints(ldouble y, SIM &svar, FLUID &fvar, CROSS &cvar, State &pn, State &pnp1)
{
	// cout << "Adding points..." << endl;
	StateVecD v = cvar.vJet;  /*Jet velocity*/
	
	// ldouble rho=fvar.rho0;
	// ldouble jetS = svar.Jet(0)+2*fvar.H;
	svar.nrefresh = 0;	
	ldouble jetR = 0.5*(svar.Jet(0));

	/*Squeeze particles together to emulate increased pressure*/
	
	
	ldouble press =fvar.pPress;
	ldouble rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

	/*Create the simulation particles*/
	for (ldouble z = -jetR; z <= jetR; z+= svar.dx)
	{ /*Do the centerline of points*/
		StateVecD xi(0.0,y,z);
		xi = svar.Rotate*xi;
		xi += svar.Start;
		pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
		pnp1.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
		++svar.simPts;
		++svar.nrefresh;
	}

	for (ldouble x = svar.dx; x < jetR ; x+=svar.dx)
	{ /*Do the either side of the centerline*/
		for (ldouble z = -jetR; z <= jetR; z+= svar.dx)
		{
			if(((x*x)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
    		{   /*If the point is inside the hole diameter, add it*/
				StateVecD temp(x,y,z);
				StateVecD xi = svar.Rotate*temp;
				xi+= svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
				++svar.simPts;
				++svar.nrefresh;

				temp(0) = -x;
				xi = svar.Rotate*temp;
				xi+= svar.Start;
				pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
				pnp1.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
				++svar.simPts;
				++svar.nrefresh;

			}
		}
		
	}

	svar.totPts += svar.nrefresh;
	++svar.addcount;
	// cout << "New points: " << svar.nrefresh << "  totPts: " <<
	// svar.totPts << " simPts: "<< svar.simPts <<  endl;
}

void CloseBoundary(SIM &svar, FLUID &fvar, State &pn, State &pnp1)
{
	std::cout << "Closing boundary..." << endl;
	StateVecD v = StateVecD::Zero();
	ldouble rho=fvar.rho0;
	ldouble press = fvar.B*(pow(rho/fvar.rho0,fvar.gam)-1);
	
	ldouble holeS = svar.Start(0); /*Distance from Box start to hole*/
	ldouble holeD = svar.Start(1); /*Diameter of hole (or width)*/
	ldouble stepb = (svar.Pstep*svar.Bstep);
	int Nb = ceil((holeD)/stepb);
	stepb = holeD/(1.0*Nb);
	State temp;

	
	ldouble jetR = 0.5*(holeD);
	for(ldouble x = holeS; x<holeS+holeD; x+=stepb)
	{
		for (ldouble z = -jetR; z <= jetR; z+= svar.Pstep)
		{
			double r = (x-holeS-jetR); /*Normalise the circle around 0,0*/
			if(((r*r)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
    		{   /*If the point is inside the hole diameter, add it*/
				StateVecD xi(x,0.0,z);
				temp.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,0));
			}
		}
	}
	
	
	/*Insert particles at the end of boundary particles*/
	pn.insert(pn.begin()+svar.bndPts,temp.begin(),temp.end());
	pnp1.insert(pnp1.begin()+svar.bndPts,temp.begin(),temp.end());
	svar.bndPts+=temp.size(); /*Adjust counter values*/
	svar.totPts +=temp.size();
}

void CreateDroplet(SIM &svar, FLUID &fvar, State &pn, State &pnp1)
{
	StateVecD v = StateVecD::Zero();
	ldouble rho = fvar.Simmass/pow(svar.dx,3.0);
	ldouble press = fvar.B*(pow(rho/fvar.rho0,fvar.gam)-1);
	svar.nrefresh = 0;	

	ldouble radius = 0.5*svar.Start(1);
	for (ldouble y = -radius; y <= radius; y+=svar.dx)
	{	
		ldouble xradius = sqrt(radius*radius - y*y);
		for (ldouble z = -xradius; z <= xradius; z+= svar.dx)
		{ /*Do the centerline of points*/
			StateVecD xi(0.0,y,z);
			pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
			pnp1.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
			++svar.simPts;
			++svar.nrefresh;
		}

		for (ldouble x = svar.dx; x <= xradius ; x+=svar.dx)
		{ /*Do the either side of the centerline*/
			for (ldouble z = -xradius; z <= xradius; z+= svar.dx)
			{
				if(((x*x) + (z*z) + (y*y)) <= (radius*radius) )
	    		{   /*If the point is inside the hole diameter, add it*/
					StateVecD xi(x,y,z);
					pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
					pnp1.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
					++svar.simPts;
					++svar.nrefresh;
					xi(0) = -x;
					pn.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
					pnp1.emplace_back(Particle(xi,v,rho,fvar.Simmass,press,1));
					++svar.simPts;
					++svar.nrefresh;
				}
			}	
		}
	}

}

#endif
