/*********   WCSPH (Weakly Compressible Smoothed Particles Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*********                   FOR THREE DIMENSION CODE                           *************/

#ifndef INIT_H
#define INIT_H

#include "Var.h"
#include "IO.h"
#include "Add.h"

void Init_Surface(SIM const& svar, MESH const& cells, vector<SURF>& surf_marks)
{
	surf_marks = vector<SURF>(svar.markers.size());

	vector<vector<size_t>> faceIDs(svar.markers.size());
	vector<vector<int>> markers(svar.markers.size());
	/* for each surface, find how many faces are in it. */
	for(std::pair<size_t,int> const& marker:cells.smarkers)
	{
		auto index = find(svar.markers.begin(),svar.markers.end(),marker.second);
		if(index != svar.markers.end())
		{
			size_t mark = index - svar.markers.begin();
			faceIDs[mark].emplace_back(marker.first);
			markers[mark].emplace_back(marker.second);

			// surf_faces[mark].back().faceID  = marker.first;
			// if()
			// cout << mark << "  " << markers.first << endl;
			// cout << surf_faces[mark].back().faceID << "  " << markers.first << endl;
			// surf_faces[mark].back().marker  = marker.second;
		}
		else
		{
			cout << "Couldn't find the marker in the index" << endl;
		}
		
		
	}

	for(size_t ii = 0; ii < svar.markers.size(); ii++)
	{
		surf_marks[ii].name = svar.bnames[ii];
		surf_marks[ii].marker = svar.markers[ii];
		surf_marks[ii].output = svar.bwrite[ii];
		
		size_t nFaces = faceIDs[ii].size();
		surf_marks[ii].faceIDs = faceIDs[ii];
		surf_marks[ii].face_count = vector<uint>(nFaces,0);
		surf_marks[ii].face_beta = vector<real>(nFaces,0.0);
		surf_marks[ii].face_area = vector<real>(nFaces,0.0);

	}

	// /* allocate the impacte vector */
	// for(size_t ii = 0; ii < nSurf; ++ii)
	// {
	// 	surfs[ii].alloc(nFaces[ii]);
	// }
}

void InitSPH(SIM& svar, FLUID const& fvar, AERO const& avar, SPHState& pn, SPHState& pnp1)
{

	cout << "Initialising simulation..." << endl;

	//Structure input initialiser
	//SPHPart(StateVecD x, StateVecD v, StateVecD vh, StateVecD f, float rho, float Rrho, bool bound) :
	StateVecD v = StateVecD::Zero();   
	// real rho=fvar.rho0;
	// real press = fvar.B*(pow(rho/fvar.rho0,fvar.gam)-1);  
	
	real press = fvar.pPress;
	real rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);

/************** Create the boundary pn  *****************/ 	 
	real step = svar.Pstep*svar.Bstep;
	int Ny = ceil((svar.bound_box(1)+4*svar.Pstep)/step);	
	int Nx = ceil((svar.bound_box(0)+4*svar.Pstep)/step);
	#if (SIMDIM == 3)
		int Nz = ceil(svar.bound_box(2)/step);
	#endif
	uint pID = 0;

	/* Bcase == 0 -> no boundary */
	if(svar.Bcase == 1) /*Rectangle*/
	{	
		#if SIMDIM == 3
			for(int j = 0; j <= Ny ; ++j) 
			{
				for (int k=0; k<=Nz; ++k) 
				{	/*Create Left and right boundary faces*/
					StateVecD xi(0.0,j*step,k*step);
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
					xi(0) = svar.bound_box(0);
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}
				
			}
			for(int i = 1; i<Nx ; ++i) 
			{
				for (int j=1; j<Ny; ++j)
				{	/*Create top and bottom boundary faces*/
					StateVecD xi(i*step,j*step,0);
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
					// xi(2)= bound_box(2); //Top boundary (Typically omitted)
					// pn.emplace_back(SPHPart(xi,v,f,rho,Rrho,bndM,0));
				}
				
			}
			for(int i= 1; i<Nx; ++i) 
			{
				for(int k = 0; k <= Nz; ++k) 
				{	/*Create far and near boundary*/
					StateVecD xi(i*step, 0, k*step);
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
					xi(1) = svar.bound_box(1);
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}
			}
		#else /*SIMDIM == 2*/
			for(int i = 0; i <= Ny ; ++i) 
			{	/* Left Wall*/
				StateVecD xi(-svar.bound_start(0)-2*svar.Pstep,i*step-svar.bound_start(1)-2*svar.Pstep);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}
	/*			// Optional lid
			for(int i = 1; i <Nx ; ++i) {
				StateVecD xi(i*stepx,bound_box(1));
				particles.emplace_back(SPHPart(xi,v,f,rho,Rrho,bndM,Bound));	
			}
			StateVecD x(stepx-svar.bound_start(0),(Ny+0.5)*stepy);
			pn.emplace_back(SPHPart(x,v,f,rho,fvar.bndM,true));
			x(0) = svar.bound_box(0) -stepx;
			pn.emplace_back(SPHPart(x,v,f,rho,fvar.bndM,true));
	*/
			for(int i= Ny; i>0; --i) 
			{	/*Right Wall*/
				StateVecD xi(Nx*step-svar.bound_start(0)-2*svar.Pstep,i*step-svar.bound_start(1)-2*svar.Pstep);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));	
				pID++;
			}
			for(int i = Nx; i > 0; --i) 
			{	/*Floor*/
				StateVecD xi(i*step-svar.bound_start(0)-2*svar.Pstep,-svar.bound_start(1)-2*svar.Pstep);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}
		#endif		
	}
	else if(svar.Bcase == 2)
	{	/*Jet*/
		int const interval = 20000;
		#if SIMDIM == 3
			real holeD = svar.jet_diam+6.75*svar.dx; /*Diameter of hole (or width)*/
			real stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			real r = 0.5*holeD;
	    	real dtheta = atan((stepb)/(r));

	    	int ncirc = floor(2.0*M_PI/dtheta);
	    	dtheta = 2.0*M_PI/real(ncirc);

			for(real ii = 0; ii < 4; ii+=1.0)
			{
				r = 0.5*holeD + ii*stepb;
				for (real y = 0; y >= -svar.jet_depth-stepb; y-=stepb)			
				{	
					for(real theta = 0; theta < 2*M_PI; theta += dtheta)
					{
						StateVecD perturb(random(interval), random(interval), random(interval));
						StateVecD xi(r*sin(theta), y, r*cos(theta));
						xi += perturb;
						pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
						pID++;
					}	
				}
			}
		#else

			real jetR = 0.5*(svar.jet_diam+2*svar.dx); /*Radius of hole (or width)*/
			real stepb = (svar.Pstep*svar.Bstep);
			
			/*Create a bit of the pipe downward.*/
			for(real ii = 0; ii < 4; ii+=1.0)
			{
				real x = jetR + ii*stepb;
				for (real y = 0.0; y >= -svar.jet_depth - 2.0 * stepb; y-=stepb)			
				{
					StateVecD perturb(random(interval), random(interval));
					StateVecD xi(-x,y);
					xi += perturb;
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}

				for (real y = 0.0; y >= -svar.jet_depth - 2.0 * stepb; y -= stepb)
				{
					StateVecD perturb(random(interval), random(interval));
					StateVecD xi(x,y);
					xi += perturb;
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}
			}
		#endif
	}
	else if (svar.Bcase == 3)
	{	/*Converging nozzle for jet spray*/
		#if SIMDIM == 3
			real holeD = svar.jet_diam+4*svar.dx; /*Diameter of hole (or width)*/
			real stepb = (svar.Pstep*svar.Bstep);
			

			/*Create a bit of the pipe downward.*/
			real r = 0.5*holeD;
	    	real dtheta = atan((stepb)/(r));
			for (real y = 0; y > -svar.jet_depth; y-=stepb)			
			{	
				for(real theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}	
			}

			/*Interpolate between the big and small diameters*/
			for (real y = -svar.jet_depth*1; y > -svar.jet_depth*2; y-=stepb)			
			{	
				real x = holeD + (y-2*svar.jet_depth)*(r-holeD)/(-svar.jet_depth);
				for(real theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(x*sin(theta), y, x*cos(theta));
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}	
			}

			/*Create the wide bit of the nozzle*/
			r = holeD;
	    	dtheta = atan((stepb)/(r));
			for (real y = -svar.jet_depth*2; y >= -svar.jet_depth*3; y-=stepb)			
			{	
				for(real theta = 0; theta < 2*M_PI; theta += dtheta)
				{
					StateVecD xi(r*sin(theta), y, r*cos(theta));
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
					pID++;
				}	
			}

		#else
			real stepb = (svar.Pstep*svar.Bstep);
			real jetR = 0.5*(svar.jet_diam+4*svar.dx); /*Diameter of hole (or width)*/
			real resR = jetR*2.0;

			/*Create the wide bit of the nozzle*/
			for (real y = -svar.jet_depth*2; y >= -svar.jet_depth*3; y-=stepb)
			{
				StateVecD xi(-resR,y);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}	

			/*Create the tapering section*/
			for (real y = -svar.jet_depth; y > -svar.jet_depth*2; y-=stepb)
			{
				/*Interpolate between resR and jetR*/
				real x = resR + (y-2*svar.jet_depth)*(jetR-resR)/(-svar.jet_depth);
				StateVecD xi(-x,y);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (real y = 0; y > -svar.jet_depth; y-=stepb)
			{
				StateVecD xi(-jetR,y);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the wide bit of the nozzle*/
			for (real y = -svar.jet_depth*2; y >= -svar.jet_depth*3; y-=stepb)
			{
				StateVecD xi(resR,y);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the tapering section*/
			for (real y = -svar.jet_depth; y > -svar.jet_depth*2; y-=stepb)
			{
				/*Interpolate between resR and jetR*/
				real x = resR + (y-2*svar.jet_depth)*(jetR-resR)/(-svar.jet_depth);
				StateVecD xi(x,y);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}

			/*Create the exit bit.*/
			for (real y = 0; y > -svar.jet_depth; y-=stepb)
			{
				StateVecD xi(jetR,y);
				pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
				pID++;
			}
		#endif
	}
	else if (svar.Bcase == 5) 
	{	
		// #if SIMDIM == 3
		// 	/*Piston driven flow*/
		// 	real stepb = (svar.Pstep*svar.Bstep);
		// 	real tankW = svar.bound_start(0)+8*svar.dx;
		// 	real tankD = svar.bound_start(1);

		// 	real jetR = 0.5*(svar.jet_diam+8*svar.dx); /*Diameter of hole (or width)*/
		// 	real resR = 0.5*tankW;
		// 	v = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
			

		// 	/*Create the piston*/
		// 	for (real z = -resR; z <= resR; z+= svar.dx)
		// 	{ /*Do the centerline of points*/
		// 		real y = -(2*svar.jet_depth+tankD+fvar.H);
		// 		StateVecD xi(0.0,y,z);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.simM,press,PartState.BOUND_,pID));
		// 		++pID;
		// 	}

		// 	for (real x = svar.dx; x < resR ; x+=svar.dx)
		// 	{ /*Do the either side of the centerline*/
		// 		for (real z = -resR; z <= resR; z+= svar.dx)
		// 		{
		// 			if(((x*x)/(resR*resR) + (z*z)/(resR*resR)) <= 1.0 )
		//     		{   /*If the point is inside the hole diameter, add it*/
		//     			real y = -(2*svar.jet_depth+tankD+fvar.H);
		// 				StateVecD xi(x,y,z);
		// 				pn.emplace_back(SPHPart(xi,v,rho,fvar.simM,press,PartState.BOUND_,pID));
		// 				++pID;
		// 				++svar.simPts;
		// 				++svar.nrefresh;

		// 				xi(0) = -x;
		// 				pn.emplace_back(SPHPart(xi,v,rho,fvar.simM,press,PartState.BOUND_,pID));
		// 				++pID;
		// 			}
		// 		}
		// 	}


		// 	/*Create the reservior*/
	 //    	real dtheta = atan((stepb)/(resR));
		// 	for (real y = (-svar.jet_depth*2-tankW); y <= -svar.jet_depth*2; y+=stepb)			
		// 	{	
		// 		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		// 		{
		// 			StateVecD xi(resR*sin(theta), y, resR*cos(theta));
		// 			/*Apply Rotation...*/
		// 			pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 			pID++;
		// 		}	
		// 	}

		// 	/*Interpolate between the big and small diameters*/
		// 	for (real y = -svar.jet_depth*2; y < -svar.jet_depth*1; y+=stepb)			
		// 	{	
		// 		real x = resR + (y-2*svar.jet_depth)*(jetR-resR)/(-svar.jet_depth);
		// 		dtheta = atan((stepb)/(x));
		// 		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		// 		{
		// 			StateVecD xi(x*sin(theta), y, x*cos(theta));
		// 			/*Apply Rotation...*/
		// 			pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 			pID++;
		// 		}	
		// 	}

		// 	/*Create the exit section.*/
		// 	dtheta = atan((stepb)/(jetR));
		// 	for (real y = -svar.jet_depth; y < 0; y+=stepb)			
		// 	{	
		// 		for(real theta = 0; theta < 2*M_PI; theta += dtheta)
		// 		{
		// 			StateVecD xi(jetR*sin(theta), y, jetR*cos(theta));
		// 			/*Apply Rotation...*/
		// 			pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 			pID++;
		// 		}	
		// 	}
		// #endif

		// #if SIMDIM == 2
		// 	/*Piston driven flow*/
		// 	real stepb = (svar.Pstep*svar.Bstep);
		// 	real tankW = svar.Start(0)+8*svar.dx;
		// 	real tankD = svar.Start(1)+4*svar.dx;

		// 	real jetR = 0.5*(svar.jet_diam+4*svar.dx); /*Diameter of hole (or width)*/
		// 	real resR = 0.5*tankW;

		// 	v = (avar.vJet*pow(jetR,2))/(0.6*pow(resR,2));
		// 	/*Create the piston*/
		// 	for(real x = -(resR-fvar.H); x < (resR -fvar.H); x+=stepb)
		// 	{
		// 		StateVecD xi(x,-(2*svar.jet_depth+tankD+3*fvar.H));
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 		svar.psnPts = pID;
		// 	}

		// 	v=StateVecD::Zero();
		// 	/*Create the reservoir tank*/
		// 	for (real y = -(svar.jet_depth*2+tankD+3*fvar.H); y <= -svar.jet_depth*2; y+=stepb)
		// 	{
		// 		StateVecD xi(-resR,y);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	real theta = atan(svar.jet_depth/(resR-jetR));
		// 	real stepy = stepb*sin(theta);
		// 	/*Create the tapering section*/
		// 	for (real y = -svar.jet_depth*2; y < -svar.jet_depth; y+=stepy)
		// 	{
		// 		/*Interpolate between resR and jetR*/
		// 		real x = resR - (y+2*svar.jet_depth)*(jetR-resR)/(-svar.jet_depth);
		// 		StateVecD xi(-x,y);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the exit bit.*/
		// 	for (real y = -svar.jet_depth; y < 0; y+=stepb)
		// 	{
		// 		StateVecD xi(-jetR,y);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the exit bit.*/
		// 	for (real y = 0; y > -svar.jet_depth; y-=stepb)
		// 	{
		// 		StateVecD xi(jetR,y);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the tapering section*/
		// 	for (real y = -svar.jet_depth; y > -svar.jet_depth*2; y-=stepy)
		// 	{
		// 		/*Interpolate between resR and jetR*/
		// 		real x = resR - (y+2*svar.jet_depth)*(jetR-resR)/(-svar.jet_depth);
		// 		StateVecD xi(x,y);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
		// 		pID++;
		// 	}

		// 	/*Create the other wall.*/
		// 	for (real y = -svar.jet_depth*2; y >= -(svar.jet_depth*2+tankD+3*fvar.H); y-=stepb)
		// 	{
		// 		StateVecD xi(resR,y);
		// 		pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.BOUND_,pID));
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
	if(svar.Scase == 2 )
	{	
		/*Cylinder case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		/* Just a cylinder case */
		real y = 0.0;
		for (y = 0.0; y >  -svar.jet_depth; y-= svar.dx)
		{
			#if SIMDIM == 3
				Add_Radial_Points(y, svar, fvar, avar, pn, PartState.FREE_);
			#else
				AddPoints(y, svar, fvar, avar, pn, PartState.FREE_);
			#endif
		}
	}
	else if(svar.Scase == 3)
	{
		/*Droplet case*/
		svar.simPts = 0;
		svar.totPts = pn.size();

		#if SIMDIM == 2
		CreateRDroplet(svar,fvar,pn);
		#else
		CreateDroplet(svar,fvar,pn);
		#endif
	}
	else if (svar.Scase == 4)
	{
		/*Jet case*/
		svar.simPts = 0;
		svar.totPts = pn.size();
		
		if (svar.Bcase == 3)
		{	/* Converging jet */
			for (real y = -svar.jet_depth*2; y > -svar.jet_depth*3; y-=svar.dx)
			{
				// cout << "In add points for-loop" << endl;
				AddPoints(y, svar, fvar, avar, pn, PartState.PIPE_);
			}

			Add_Buffer(svar,fvar,pn);
		}
		else
		{	/* Normal Jet */
			real y = 0.0;
			for (y = 0.0; y >  -svar.jet_depth; y -= svar.dx)
			{
				#if SIMDIM == 3
					Add_Radial_Points(y, svar, fvar, avar, pn, PartState.PIPE_);
				#else
					AddPoints(y, svar, fvar, avar, pn, PartState.PIPE_);
				#endif
			}

			Add_Buffer(svar,fvar,pn);
		}
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
						StateVecD xi(i*svar.dx,	j*svar.dx, k*svar.dx);		
						pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.FREE_,pID));
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
					StateVecD xi(i*svar.dx,j*svar.dx);		
					pn.emplace_back(SPHPart(xi,v,rho,fvar.bndM,press,PartState.FREE_,pID));
					// pn.back().cellP = 20000-fvar.gasPress;
					pID++;
				}
			}
		#endif

	}


		
	// svar.simPts+=10*10;

	svar.totPts = pn.size();
	svar.partID = pn.size();

	

	for(size_t ii = 0; ii < svar.bndPts; ++ii)
		pn[ii].xi = (svar.Rotate * pn[ii].xi) + svar.bound_start;

	for(size_t ii = svar.bndPts; ii < svar.totPts; ++ii)
	{
		pn[ii].xi = (svar.Rotate * pn[ii].xi) + svar.sim_start;
		// if(svar.Asource == 0)
		// {
		// 	pn[ii].cellV = avar.vInf;
		// }
	}
	pnp1 = pn;

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
		cerr<< "Particles array size doesn't match defined values." << endl;
		cerr<< "Total Points: " << svar.totPts << "  Boundary Points: " << svar.bndPts
			<< "  Simulation Points: " << svar.simPts << endl;
		exit(-1);
	}
	// cout << "Total pn: " << svar.totPts << endl;
	
	// cout << "Refresh round pn: " << svar.nrefresh << endl;
}

#endif