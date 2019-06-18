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
	StateVecD f = StateVecD::Zero();
	ldouble rho=fvar.rho0;
	ldouble jetS = svar.Start(0)+2*fvar.H;
	svar.nrefresh = 0;	
	ldouble jetR = 0.5*(svar.Start(1));

	/*Create the simulation particles*/
	for (ldouble z = -jetR; z <= jetR; z+= svar.Pstep)
	{ /*Do the centerline of points*/
		StateVecD xi(-(jetS+jetR),y,z);
		pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
		pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
		++svar.simPts;
		++svar.nrefresh;
	}

	for (ldouble x = svar.Pstep; x < jetR ; x+=svar.Pstep)
	{ /*Do the either side of the centerline*/
		for (ldouble z = -jetR; z <= jetR; z+= svar.Pstep)
		{
			double r = (x-jetS-jetR); /*Normalise the circle around 0,0*/
			// cout << (r*r)/(jetR*jetR) + (z*z)/(jetR*jetR) << endl;

			if(((x*x)/(jetR*jetR) + (z*z)/(jetR*jetR)) <= 1.0 )
    		{   /*If the point is inside the hole diameter, add it*/
				StateVecD xi(r,y,z);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,1));
				pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,1));
				++svar.simPts;
				++svar.nrefresh;
				xi(0) = r-2*x;
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
				pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
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
	StateVecD f = StateVecD::Zero();
	ldouble rho=fvar.rho0;

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
				temp.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
			}
		}
	}
	
	
	/*Insert particles at the end of boundary particles*/
	pn.insert(pn.begin()+svar.bndPts,temp.begin(),temp.end());
	pnp1.insert(pnp1.begin()+svar.bndPts,temp.begin(),temp.end());
	svar.bndPts+=temp.size(); /*Adjust counter values*/
	svar.totPts +=temp.size();
}

#endif
