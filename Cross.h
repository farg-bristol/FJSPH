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
	ldouble jetE = jetS+svar.Start(1);
	svar.nrefresh = 0;


	if (cvar.acase == 5)
	{
		/*Create the simulation particles*/
		for( ldouble x = jetS; x<jetE - (jetE-jetS)/2; x+=svar.Pstep)
		{ /*Do the left set of points*/
			StateVecD xi(x,y);
			int pos = svar.totPts+svar.nrefresh;
			pn.insert(pn.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,2));
			pnp1.insert(pnp1.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,2));
			++svar.simPts;
			++svar.nrefresh;
		}

		ldouble start2 = pn.back().xi[0]+svar.Pstep;
		for( ldouble x = start2; x<=jetE; x+=svar.Pstep)
		{	/*Do the right set of points*/
			StateVecD xi(x,y);
			int pos = svar.totPts+svar.nrefresh;
			pn.insert(pn.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,1));
			pnp1.insert(pnp1.begin()+pos,Particle(xi,v,f,rho,fvar.Simmass,1));
			++svar.simPts;
			++svar.nrefresh;
		}
	}
	else
	{
		/*Create the simulation particles*/
		for( ldouble x = jetS; x<jetE - (jetE-jetS)/2; x+=svar.Pstep)
		{ /*Do the left set of points*/
			StateVecD xi(x,y);
			pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
			pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,2));
			++svar.simPts;
			++svar.nrefresh;
		}

		ldouble start2 = pn.back().xi[0]+svar.Pstep;
		for( ldouble x = start2; x<=jetE; x+=svar.Pstep)
		{	/*Do the right set of points*/
			StateVecD xi(x,y);
			pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,1));
			pnp1.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,1));
			++svar.simPts;
			++svar.nrefresh;
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

	for(ldouble x = holeS; x<holeS+holeD; x+=stepb)
	{
		StateVecD xi(x,0.0);
		temp.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
	}
	/*Insert particles at the end of boundary particles*/
	pn.insert(pn.begin()+svar.bndPts,temp.begin(),temp.end());
	pnp1.insert(pnp1.begin()+svar.bndPts,temp.begin(),temp.end());
	svar.bndPts+=temp.size(); /*Adjust counter values*/
	svar.totPts +=temp.size();
}

#endif
