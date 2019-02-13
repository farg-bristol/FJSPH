/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"
#include "IO.h"
#include "Cross.h"

using namespace std;

/*Make a guess on how big the array is going to be (doesn't need to be totally exact)*/
int ParticleCount(SIM &svar, FLUID &fvar)
{
	int partCount = 0;
	ldouble stepx = svar.Pstep*svar.Bstep;
	ldouble stepy = svar.Pstep*svar.Bstep;
	int Ny = ceil(svar.Box(1)/stepy);
	int Nx = ceil(svar.Box(0)/stepx);
	switch(svar.Bcase)
	{
		case 0:
			break;
		case 1:
			partCount = 2*Ny + Nx; /*Boundary particles*/
			partCount += svar.simPts; /*Simulation particles*/
			break;
		case 3:
		{	/*Boundary Points*/
			ldouble holeS = svar.Start(0); /*Distance from Box start to hole*/
			ldouble holeD = svar.Start(1)+4*svar.Pstep; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			int Nb = ceil((holeS)/stepb);
			stepb = holeS/(1.0*Nb);
			/*Actual values*/			
			int preHole = Nb+1 + ceil(svar.Box[1]/stepb);
			int postHole = ceil(svar.Box[1]/stepb) + ceil((svar.Box[0]-holeS-holeD)/stepb)+1;

			/*Simulation Points*/
			ldouble jetS = svar.Start(0)+2*fvar.H;
			ldouble jetE = svar.Start(0)+svar.Start(1);
			int addRound = floor((jetE-jetS)/svar.Pstep);
			/*Need to add the particles already present*/
			int simPts = addRound*svar.nmax + addRound*ceil(svar.Box[1]/svar.Pstep);

			partCount = preHole + postHole + simPts;
			
			break;
		}

			
	}
	
	return partCount;
}

void InitSPH(SIM &svar, FLUID &fvar, CROSS &cvar, State &pn, State &pnp1)
{


	switch (svar.Bcase)
	{
		default:
			cout << "Initialising simulation with " << svar.simPts << " particles" << endl;
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
	 
/************** Create the boundary particles  *****************/ 	 
	ldouble stepx = svar.Pstep*svar.Bstep;
	ldouble stepy = svar.Pstep*svar.Bstep;
	
	int Ny = ceil(svar.Box(1)/stepy);
	stepy = svar.Box(1)/Ny;	/*Find exact step to meet dimensions*/
	
	int Nx = ceil(svar.Box(0)/stepx);
	stepx = svar.Box(0)/Nx;

	
	switch (svar.Bcase) 
	{
		case 0:
		{ /*No boundary*/
			break;
		}
		case 1: /*Rectangle*/
		{
			for(int i = 0; i <= Ny ; ++i) {
				StateVecD xi(-svar.Start(0),i*stepy-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}
/*			
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
				StateVecD xi(svar.Box(0)-svar.Start(0),i*stepy-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,true));	
			}
			for(int i = Nx; i > 0; --i) {
				StateVecD xi(i*stepx-svar.Start(0),-svar.Start(1));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}
			break;
		}
		case 2: /*Bowl*/
		{
			
			ldouble r= svar.Box(0);
			ldouble dtheta = atan((svar.Pstep*svar.Bstep)/r);
			for (ldouble theta = 0.0; theta < M_PI; theta+=dtheta)
			{
				StateVecD xi(-r*cos(theta),r*(1-sin(theta)));
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}
			break;
		}
		case 3:
		{	/*Jet in Crossflow*/
			ldouble holeS = svar.Start(0); /*Distance from Box start to hole*/
			ldouble holeD = svar.Start(1)+4*svar.Pstep; /*Diameter of hole (or width)*/
			ldouble stepb = (svar.Pstep*svar.Bstep);
			int Nb = ceil((holeS)/stepb);
			stepb = holeS/(1.0*Nb);
			
			for (ldouble x=0.0; x<holeS; x+= stepb)
			{
				StateVecD xi(x,0.0);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}

			/*Create a bit of the pipe downward.*/
			for (ldouble y = -stepb; y >= -svar.Box[1]-stepb; y-=stepb)			{
				StateVecD xi(pn.back().xi[0],y);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}

			for (ldouble y = pn.back().xi[1]; y < 0.0 ; y+=stepb)
			{
				StateVecD xi(holeS+holeD,y);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}


			for(ldouble x = holeS+holeD; x<=svar.Box(0); x+=stepb)
			{
				StateVecD xi(x,0.0);
				pn.emplace_back(Particle(xi,v,f,rho,fvar.Boundmass,0));
			}
			break;
		}
		default: 
		{
			cerr << "Boundary case is not within the design. 0 <= Bcase <= 4." << endl;
			exit(-1);
		}
	}
	
	svar.bndPts = pn.size();
	
/***********  Create the simulation particles  **************/

	
	switch(svar.Bcase)
	{
		case 3: 
		{	/*Crossflow case*/
			svar.simPts = 0;
			svar.totPts = pn.size();
			/*Update n+1 before adding sim particles*/
			for (auto p: pn)
				pnp1.emplace_back(p);

			for (double y = 0.0; y > -svar.Box[1]; y-=svar.Pstep)
				AddPoints(y, svar, fvar, cvar, pn, pnp1);

			break;
		}

		default:
		{
			int which = 2 ;
			if (which ==1) 
			{	/*Hexahedral start*/
				svar.xyPART(0) = 75;
				svar.xyPART(1) = 150;
				svar.simPts = 75*150;
				stepx = sqrt(2)*svar.Pstep;
				stepy = stepx/2;

				for(int j=0; j< svar.xyPART(1); ++j)
				{
					for( int i=0; i< svar.xyPART(0); ++i) 
					{
						double indent = 0.0;
						if (i != svar.xyPART(0) && (j+1)%2 == 0) indent = stepy;

						double x = /*svar.Start(0)+*/i*stepx+indent;
						double y = /*svar.Start(1)+*/j*stepy;
						StateVecD xi(x,y);		
						pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
					}
				}
				// fvar.height0 = (svar.xyPART(1)-1)*svar.Pstep/sqrt(2);
			}
			else if (which == 2)
			{	/*Square lattice start*/
				for( int i=0; i< svar.xyPART(0); ++i) 
				{
					for(int j=0; j< svar.xyPART(1); ++j)
					{				
							StateVecD xi(i*svar.Pstep,j*svar.Pstep);		
							pn.emplace_back(Particle(xi,v,f,rho,fvar.Simmass,0));
					}
				}
				// fvar.height0 = (svar.xyPART(1)-1)*svar.Pstep;
			}

			for (auto p: pn)
				pnp1.emplace_back(p);

			break;
		}
	}
		
	// svar.simPts+=10*10;

	svar.totPts = pn.size();
	
	if(svar.totPts!=svar.bndPts+svar.simPts)
	{
		cerr<< "Mismatch of particle count." << endl;
		cerr<< "Particle array size doesn't match defined values." << endl;
		write_settings(svar,fvar);
		exit(-1);
	}
	cout << "Total Particles: " << svar.totPts << endl;
	// cout << "Boundary Particles: " << svar.bndPts << endl;
	// cout << "Sim Particles:" << svar.simPts << endl;
	// cout << "Refresh round Particles: " << svar.nrefresh << endl;
}
#endif