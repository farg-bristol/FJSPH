/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <chrono>

#include "Var.h"
#include "IO.h"
#include "Neighbours.h"
#include "Resid.h"
#include "Crossing.h"
#include "Newmark_Beta.h"
// #include "Runge_Kutta.h"


///**************** Integration loop **************///
real Integrate(KDTREE& TREE, SIM& svar, const FLUID& fvar, const AERO& avar, 
	MESH& cells, DELTAP& dp, State& pn, State& pnp1, State& airP, outl& outlist)
{
	// cout << "Entered Newmark_Beta" << endl;
	uint   k = 0;	
	real logbase = 0.0;
	real error1 = 0.0;
	real error2 = 0.0;

	// Find maximum safe timestep
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});

	vector<Particle>::iterator maxUi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](Particle p1, Particle p2){return p1.v.norm()< p2.v.norm();});

	real maxf = maxfi->f.squaredNorm();
	real maxU = maxUi->v.norm();
	real dtv = fvar.HSQ * fvar.rho0/fvar.mu;
	real dtf = sqrt(fvar.H/maxf);
	real dtc = fvar.H/(maxU);

	// Only use if -fno-finite-math-only is on
	// if (std::isinf(maxf))
	// {
	// 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
	// 	exit(-1);
	// }

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	svar.dt = 0.125*std::min(dtf,std::min(dtc,dtv));
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

// #ifdef DEBUG
	cout << "time: " << svar.t << " dt: " << svar.dt << "  dtv: " << dtv <<  "  dtf: " << dtf << "  dtc: " << dtc << " Maxf: " << maxf << endl;
// #endif

	if (svar.dt > (svar.frame+1)*svar.framet-svar.t)
	{
		svar.dt = (svar.frame+1)*svar.framet-svar.t;
	}
	
	// Check if the particle has moved to a new cell
	if (svar.Asource == 1 || svar.Asource == 2)
	{
		// cout << "Finding cells" << endl;
		FindCell(svar,fvar.sr,TREE,cells,pnp1,pn);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.simPts -= svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
		}	
	}

	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;
	// const size_t piston = svar.psnPts;

	// cout << "Cells found" << endl;

	vector<size_t> cellsused; // Cells that contain a particle
	if(svar.Asource == 2)
	{
		vector<size_t> tempcell;
		#pragma omp parallel shared(pnp1)
		{
			std::vector<size_t> local;
			#pragma omp for schedule(static) nowait
			for (size_t ii = start; ii < end; ++ii)
			{	
				if(pnp1[ii].b == PartState.FREE_)
					local.emplace_back(pnp1[ii].cellID);
			}

			#pragma omp for schedule(static) ordered
			for(int i=0; i<omp_get_num_threads(); i++)
			{
				#pragma omp ordered
				tempcell.insert(tempcell.end(),local.begin(),local.end());
			}
		}

		// Sort the vector and remove unique values.
		std::unordered_set<size_t> s;
		for(size_t const& ii:tempcell)
			s.insert(ii);

		cellsused.assign(s.begin(),s.end());
		std::sort(cellsused.begin(),cellsused.end());
	}


	// Check if a particle is running low on neighbours, and add ficticious particles
	std::vector<std::vector<Part>> air;
	if(svar.ghost == 1)
	{
		air = vector<vector<Part>>(start,vector<Part>());
		#pragma omp parallel shared(pnp1, outlist)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (size_t ii = start; ii < end; ++ii)
			{
				std::vector<Part> temp;
				if(pnp1[ii].surf == 1 && pnp1[ii].b == PartState.FREE_ && outlist[ii].size() > 0.4*avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,cells,ii,pnp1,outlist,dp.norm[ii],dp.avgV[ii]);
				}
				local.emplace_back(temp);
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i<omp_get_num_threads(); i++)
	    	{
	    		#pragma omp ordered
	    		air.insert(air.end(),local.begin(),local.end());
	    	}
		}
	}
	airP.clear();

	// #pragma omp parallel for schedule(static) ordered
	for(size_t ii = 0; ii < air.size(); ++ii)
		for(size_t jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	#pragma omp parallel for shared(outlist) schedule(static)
	for(size_t ii = start; ii < end; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}

	vector<StateVecD> xih(end);
	
	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	const real B = fvar.B;
	const real gam = fvar.gam;

	StateVecD dropVel = StateVecD::Zero();
	StateVecD Force = StateVecD::Zero();

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
	
	// airP.clear();
	std::vector<std::vector<Part>> neighb;
	neighb.reserve(end);
	

	#pragma omp parallel shared(pnp1, outlist, air)
	{
		std::vector<std::vector<Part>> local;
		if(svar.ghost == 1 )
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				if(air[ii].size()>0)
				{
					temp.insert(temp.end(), air[ii].begin(), air[ii].end());
				}
				local.push_back(temp);
			}
		}
		else
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				local.push_back(temp);
			}
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),local.begin(),local.end());
    	}
	}

	
	/*Get preliminary new state to find neighbours, then freeze*/
	Get_Resid(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,
		cells,cellsused,neighb,outlist,dp,pn,pnp1,airP,Force,dropVel);

	uint nUnstab = 0;

	void(Check_Error(TREE,svar,fvar,start,end,error1,error2,logbase,
							cellsused,outlist,xih,pn,pnp1,k,nUnstab));
	k++;

	// State st_2 = pn;

	// void(Get_First_RK(TREE,svar,fvar,avar,start,end,cells,cellsused,
	// 							neighb,outlist,dp,logbase,pn,st_2,error1));

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(svar,fvar,start,end,pnp1,outlist,dp);

	Detect_Surface(svar,fvar,avar,start,end,dp,outlist,cells,pnp1);

	// Apply_XSPH(fvar,start,end,outlist,dp,pnp1);
	Particle_Shift(svar,fvar,start,end,outlist,dp,pnp1);
	
	air.clear();
	if(svar.ghost == 1)
	{
		air = vector<vector<Part>>(start,vector<Part>());
		#pragma omp parallel shared(pnp1, outlist)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (size_t ii = start; ii < end; ++ii)
			{
				std::vector<Part> temp;
				if(pnp1[ii].surf == 1 && pnp1[ii].b == PartState.FREE_ && outlist[ii].size() > 0.4*avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,cells,ii,pnp1,outlist,dp.norm[ii],dp.avgV[ii]);
				}
				local.emplace_back(temp);
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i < omp_get_num_threads(); i++)
	    	{
	    		#pragma omp ordered
	    		air.insert(air.end(),local.begin(),local.end());
	    	}
		}
	}
	airP.clear();

	// #pragma omp parallel for schedule(static) ordered
	for(size_t ii = 0; ii < air.size(); ++ii)
		for(size_t jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	neighb.clear();
	neighb.reserve(end);
	
	dropVel = StateVecD::Zero();

	#pragma omp parallel shared(pnp1, outlist, air)
	{
		std::vector<std::vector<Part>> local;
		if(svar.ghost == 1 )
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				if(air[ii].size()>0)
				{
					temp.insert(temp.end(), air[ii].begin(), air[ii].end());
				}
				local.push_back(temp);
			}
		}
		else
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				local.push_back(temp);
			}
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),local.begin(),local.end());
    	}
	}

	if(svar.Asource == 2)
	{
		cellsused.clear();
		vector<size_t> tempcell;
		#pragma omp parallel shared(pnp1)
		{
			std::vector<size_t> local;
			#pragma omp for schedule(static) nowait
			for (size_t ii = start; ii < end; ++ii)
			{	
				if(pnp1[ii].b == PartState.FREE_)
					local.emplace_back(pnp1[ii].cellID);
			}

			#pragma omp for schedule(static) ordered
			for(int i=0; i<omp_get_num_threads(); i++)
			{
				#pragma omp ordered
				tempcell.insert(tempcell.end(),local.begin(),local.end());
			}
		}

		// Sort the vector and remove unique values.
		std::unordered_set<size_t> s;
		for(size_t const& ii:tempcell)
			s.insert(ii);

		cellsused.assign(s.begin(),s.end());
		std::sort(cellsused.begin(),cellsused.end());

		#pragma omp parallel for schedule(static)
		for(size_t const& ii : cellsused)
		{
			// Work out the mass and volume fractions
			if(cells.fNum[ii] != 0)
			{
				real fVol = real(cells.fNum[ii]) * avar.pVol;

				real aFrac = (cells.cVol[ii]-fVol)/cells.cVol[ii];

				
				if (aFrac < 0.2)
					continue;
				
				real aMass = cells.cMass[ii]*aFrac;

				// Do the momentum exchange
				StateVecD newPert = (cells.fMass[ii]/aMass)*(cells.vFnp1[ii]-cells.vFn[ii])/real(cells.fNum[ii]);
				// StateVecD diffusion = 0.2*cells.cPertn[ii];
				cells.cPertnp1[ii] = cells.cPertn[ii]*0.75 - newPert;

// 				#pragma omp critical
// 				{
// 				cout << "Cell " << ii << ":" << endl;

// 				cout << "pertnp1: " << cells.cPertnp1[ii](0) << "  "
// 						<< cells.cPertnp1[ii](1) 
// #if SIMDIM == 3
// 				 		<< "  " << cells.cPertnp1[ii](2)
// #endif
// 						<< endl;

// 				cout << "pertn:   " << cells.cPertn[ii](0) << "  "
// 				 		<< cells.cPertn[ii](1) 
// #if SIMDIM == 3
// 				 		<< "  " << cells.cPertn[ii](2)
// #endif
// 						<< endl;
// 				cout << "Update: " << newPert(0) << "  " << newPert(1) 
// #if SIMDIM == 3
// 				<< "  " << newPert(2)
// #endif
// 				 << endl; 

// 				cout << "Fuel count: " << cells.fNum[ii] << endl;
// 				cout << "Mass fraction: " << cells.fMass[ii]/aMass << endl;
// 				// cout << "Fuel Volume: " << fVol << " Fuel Mass: " << cells.fMass[ii] << endl;
// 				// // cout << "Cell Volume: " << cells.cVol[ii] << " Air fraction: " << aFrac << "  Air Mass: " << aMass << endl;
// 				cout << "Fuel Vel difference: " << cells.vFnp1[ii](0)-cells.vFn[ii](0) << "  " << cells.vFnp1[ii](1)-cells.vFn[ii](1)
// #if SIMDIM == 3
// 					<< "  " << cells.vFnp1[ii](2)-cells.vFn[ii](2)
// #endif
// 				 	<< endl << endl;
// 				}
			}
			// else
			// {
			// 	cout << "Cell with no fuel in it considered" << endl;
			// }
		}
	}


	/*Do time integration*/
	Newmark_Beta(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,cells,cellsused,neighb,
		outlist,dp,logbase,k,error1,error2,pn,pnp1,airP,Force,dropVel);

	// error1 = Runge_Kutta4(TREE,svar,fvar,avar,start,end,cells,cellsused,neighb,outlist,
	// 		dp,logbase,pn,st_2,pnp1,airP,Force,dropVel);

	/*Add time to global*/
	svar.t+=svar.dt;


	// cout << "Timestep Params: " << maxf << " " << fvar.Cs + svar.maxmu << " " << dtf << " " << dtcv << endl;
	// cout << "New Time: " << svar.t << endl;

	/*Check if more particles need to be created*/
	if(svar.Bcase == 2 || svar.Bcase == 3)
	{
		if(svar.Bclosed == 0)
		{
			uint nAdd = 0;
			for (auto& back:svar.back)
			{	
				if(svar.totPts < svar.nmax)
				{
					/*Check that the starting area is clear first...*/
					StateVecD vec = svar.Transp*(pnp1[back].xi-svar.Start);
					real clear;
					if(svar.Bcase == 2)
						clear =  svar.dx-svar.Jet(1)*3;
					else
						clear = svar.dx-svar.Jet(1);
					
					if(vec[1] > clear)
					{	/*Create a new particle behind the last one*/
						StateVecD xi = vec;
						xi(1) -= svar.dx;
						xi = svar.Rotate*xi;
						xi += svar.Start;

						pn.emplace_back(Particle(xi,pnp1[back],svar.totPts,PartState.BACK_));
						pnp1.emplace_back(Particle(xi,pnp1[back],svar.totPts,PartState.BACK_));

						pnp1[back].b = PartState.START_;
						/*Update the back vector*/
						back = svar.totPts;
						svar.simPts++;
						svar.totPts++;
						nAdd++;
					}
				}
				else
				{
					svar.Bclosed = 1;
				}	
			}
			
			if(nAdd != 0)
			{
				TREE.NP1.index->buildIndex();
				FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
			}
		}
	}

	/****** UPDATE TIME N ***********/
	if(svar.totPts != pnp1.size())
	{
		cout << "Size mismatch. Total points not equal to array size. Stopping" << endl;
		exit(-1);
	}

	// cout << "Updating. Time: " << svar.t << "  dt: " << svar.dt << endl;
	pn = pnp1;

	if(svar.Asource == 2)
	{
		cells.cPertn = cells.cPertnp1;
		
		// Work out the mass and volume fractions
		real tMom = 0.0;
		real aMomT = 0.0;
		// vector<real> aMom(cells.size(),0.0);
		real fMom = 0.0;
		// cout << cells.size() << endl;

		#pragma omp parallel for schedule(static)
		for(size_t ii = 0; ii < cells.size(); ii++)
		{
			real fVol = real(cells.fNum[ii]) * avar.pVol;

			real aFrac = (cells.cVol[ii]-fVol)/cells.cVol[ii];

			if(aFrac > 1)
				aFrac = 1;
			else if (aFrac < 0)
				aFrac = 0;

			real aMass = cells.cMass[ii]*aFrac;

			// cout << aMass << " " << cells.cVel[ii](0) << " " << cells.cPertnp1[ii](0) << endl;
			// aMom[ii] = ( aMass*(cells.cVel[ii]+cells.cPertnp1[ii]).norm());
			#pragma omp critical
			{
				aMomT += aMass*(cells.cVel[ii]+cells.cPertnp1[ii]).norm();
				tMom += aMass*(cells.cVel[ii]+cells.cPertnp1[ii]).norm();
			}
			
		}

		for(size_t jj = 0; jj < pnp1.size(); ++jj)
		{
			fMom += pnp1[jj].m * pnp1[jj].v.norm();
		}

		tMom += fMom;
		// cout << "fVol: " << fVol << " aFrac: " << aFrac  << " cMass: " << cells.cMass[0] << " aMass: " << aMass <<  endl;
		// cout << "cVel: "  << cells.cVel[0](0) << " " << cells.cVel[0](1) << " " << cells.cVel[0](2) << endl;

		// cout << pnp1[0].cellV(0) << "  " << pnp1[0].cellV(1) << "  " << pnp1[0].cellV(2) << endl;
		// cout << pn[0].cellV(0) << "  " << pn[0].cellV(1) << "  " << pn[0].cellV(2) << endl;
		// pertLog << svar.t << " " << tMom-(svar.tMom) << " " << aMomT-svar.aMom << " " << fMom << endl;
	}

// 	if (svar.Bcase == 4)
// 	{
// 		// Calculate the force expected for a droplet of the same size.
// 		svar.Force = Force;

// #if SIMDIM == 3
// 		real radius = svar.diam* std::cbrt(3.0/(4.0*M_PI));
// #else
// 		real radius = svar.diam / std::sqrt(M_PI);
// #endif
// 		// Average velocity of the particles
// 		dropVel /= real(svar.totPts);
// 		StateVecD Vdiff = avar.vInf - dropVel;

// 		real Re = 2.0*avar.rhog*Vdiff.norm()*radius/avar.mug;
// 		real Cds = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));

// #if SIMDIM == 3
// 		real Adrop = M_PI*radius*radius;
// #else
// 		real Adrop = radius;
// #endif

// 		// cout << Re << "  " << Cds << " " << Adrop << " " << Vdiff.norm() << "  " << svar.mass << endl;

// 		// Undeformed drag
// 		StateVecD dropForce = (0.5*avar.rhog*Vdiff.norm()*Vdiff*Cds*Adrop);

// 		// Calculate deformed drag
// 		AERO bigdrop;
// 		bigdrop.rhog = avar.rhog;
// 		bigdrop.mug = avar.mug;
// 		bigdrop.aPlate = radius*2;
// 		bigdrop.nfull = avar.nfull;

// 		// real temp = radius / std::cbrt(3.0/(4.0*M_PI));
// 		GetYcoef(bigdrop,fvar,svar.diam);

// 		real ymax = Vdiff.squaredNorm()*bigdrop.ycoef;

// 		Re = 2.0*avar.rhog*Vdiff.norm()*bigdrop.L/avar.mug;
// 		Cds = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.000001));

// 		// cout << bigdrop.ycoef << endl;
// 		// cout << radius << "  "  << bigdrop.L << "  " << ymax << endl;
// 		// cout << bigdrop.ycoef << endl;
		
// 		// if (ymax > 1.0)
// 		// 	ymax = 1.0;


// 		real const Cdl = Cds*(1+2.632*ymax);

// 		#if SIMDIM == 3 
// 			Adrop = M_PI*pow((bigdrop.L + bigdrop.Cb*bigdrop.L*ymax),2);

// 		#endif
// 		#if SIMDIM == 2
// 			Adrop = 2*(bigdrop.L + bigdrop.Cb*bigdrop.L*ymax);
// 		#endif

// 		StateVecD dropDefForce =  0.5*avar.rhog*Vdiff.norm()*Vdiff*Cdl*Adrop;	

// 		StateVecD test = GisslerForce(bigdrop, Vdiff, 1.0, 1.0, 0);

// 		cout << "Time:  " << svar.t << "  " << dropForce.norm() << "  " << dropDefForce.norm() << "  " << test.norm() << "  " 
// 			<< svar.Force.norm() << "  " << svar.Force(1) << endl;
// 	}

	return error1;
}

void First_Step(KDTREE& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	MESH& cells, outl& outlist, DELTAP& dp, State& pnp1, State& pn, State& airP)
{
	const size_t start = svar.bndPts;
	const size_t end = svar.totPts;
	// cout << "Calculating first step" << endl;

	#if DEBUG 
		dbout << "Starting first step. ";
		dbout << "  Start index: " << start << "  End index: " << end << endl;
	#endif

	if (svar.Asource == 1 || svar.Asource == 2)
	{
		FindCell(svar,fvar.sr,TREE,cells,pnp1,pn);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
		}
	}

	// cout << "Retrieved the cells" << endl;

	std::vector<std::vector<Part>> neighb;
	neighb.reserve(end);
	for(size_t ii = 0; ii < start; ++ii)
		neighb.emplace_back();

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(svar,fvar,start,end,pnp1,outlist,dp);
	cout << "Got past the first prestep" << endl;
	
	Detect_Surface(svar,fvar,avar,start,end,dp,outlist,cells,pnp1);
	
	vector<vector<Part>> air;
	if(svar.ghost == 1)
	{
		air = vector<vector<Part>>(start,vector<Part>());
		#pragma omp parallel shared(pnp1, outlist)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (size_t ii = start; ii < end; ++ii)
			{
				std::vector<Part> temp;
				if(pnp1[ii].surf == 1 && pnp1[ii].b == PartState.FREE_ && outlist[ii].size() > 0.4*avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,cells,ii,pnp1,outlist,dp.norm[ii],dp.avgV[ii]);
				}
				local.emplace_back(temp);
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i<omp_get_num_threads(); i++)
	    	{
	    		#pragma omp ordered
	    		air.insert(air.end(),local.begin(),local.end());
	    	}
		}
	}
	airP.clear();

	// #pragma omp parallel for schedule(static) ordered
	for(size_t ii = 0; ii < air.size(); ++ii)
		for(size_t jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	neighb.clear();
	neighb.reserve(end);
	

	#pragma omp parallel shared(pnp1, outlist, air)
	{
		std::vector<std::vector<Part>> local;
		if(svar.ghost == 1 )
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				if(air[ii].size()>0)
				{
					temp.insert(temp.end(), air[ii].begin(), air[ii].end());
				}
				local.push_back(temp);
			}
		}
		else
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				local.push_back(temp);
			}
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),local.begin(),local.end());
    	}
	}

	uint k = 0;	
	real logbase = 0.0;
	real error1 = 0.0;
	real error2 = 0.0;

	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	const real B = fvar.B;
	const real gam = fvar.gam;

	StateVecD dropVel = StateVecD::Zero();
	StateVecD Force = StateVecD::Zero();
	
	/*find force at time n*/
	vector<StateVecD> res;
	vector<StateVecD> Af;
	vector<real> Rrho;
	vector<real> curve;

	Forces(svar,fvar,avar,cells,pnp1,neighb,outlist,dp,res,Rrho,Af,Force,curve); 

	// #pragma omp parallel for shared(res,Rrho) schedule(static)
	// for(size_t ii = 0; ii < end; ++ii)
	// {
	// 	pn[ii].f = res[ii];
	// 	pn[ii].Rrho = Rrho[ii];
	// }

	// State st_2 = pn;
	vector<size_t> cellsused;
	
	real dt = svar.dt;
	real dt2 = dt*dt;
	// void(Get_First_RK(TREE,svar,fvar,avar,start,end,cells,cellsused,
	// 							neighb,outlist,dp,logbase,pn,st_2,error1));

	// /*Previous State for error calc*/
	vector<StateVecD> xih(svar.totPts);
	vector<StateVecD> xi(svar.totPts); /*Keep original positions to reset after finding forces*/
	#pragma omp parallel for shared(pnp1)
	for (size_t  ii=0; ii < end; ++ii)
		xih[ii] = pnp1[ii].xi;

	

	#pragma omp parallel for shared(res, Rrho, Af/*, wDiff, norm, curve*/)
	for(size_t ii = 0; ii < end; ++ii)
	{
		
		pnp1[ii].Af = Af[ii]; 
		pnp1[ii].rho = pn[ii].rho+dt*(b*Rrho[ii]);
		pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
		// pnp1[ii].p = fvar.Cs*fvar.Cs * (pnp1[ii].rho - fvar.rho0);
		
		// cout << res[ii](0) << "  " << res[ii](1) << endl;			

		if (pnp1[ii].b == PartState.START_ || pnp1[ii].b == PartState.BACK_)
		{
			xih[ii] = pn[ii].xi + dt*pn[ii].v;
		}
		else if (pnp1[ii].b > PartState.START_)
		{
			// pnp1[ii].v = pn[ii].v + dt*b*res[ii]; 
			pnp1[ii].f = res[ii];
			pnp1[ii].Rrho = Rrho[ii];
			xih[ii] = pn[ii].xi + dt*pn[ii].v + dt2*(d*res[ii]);
		}

		// pnp1[ii].theta = kern[ii]/ 5.51645e+009;
		pnp1[ii].theta = dp.kernsum[ii];
		pnp1[ii].s = outlist[ii].size();

		// pnp1[ii].theta = wDiff[ii];
		pnp1[ii].nNeigb = real(outlist[ii].size());
		pnp1[ii].bNorm = dp.norm[ii];

		// Force += res[ii];
		// dropVel += pnp1[ii].v;

		xi[ii] = pnp1[ii].xi;

	}


	/*Find maximum safe timestep*/
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});

	vector<Particle>::iterator maxUi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.v.norm()< p2.v.norm();});

	real maxf = maxfi->f.squaredNorm();
	real maxU = maxUi->v.norm();
	real dtv = fvar.HSQ * fvar.rho0/fvar.mu;
	real dtf = sqrt(fvar.H/maxf);
	real dtc = fvar.H/(maxU);

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	dt = 0.125*std::min(dtf,std::min(dtc,dtv));
	svar.dt = dt;
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(svar,fvar,start,end,pnp1,outlist,dp);


	Detect_Surface(svar,fvar,avar,start,end,dp,outlist,cells,pnp1);
	
	air.clear();
	if(svar.ghost == 1)
	{
		air = vector<vector<Part>>(start,vector<Part>());
		#pragma omp parallel shared(pnp1, outlist)
		{
			std::vector<std::vector<Part>> local;
			#pragma omp for schedule(static) nowait
			for (size_t ii = start; ii < end; ++ii)
			{
				std::vector<Part> temp;
				if(pnp1[ii].surf == 1 && pnp1[ii].b == PartState.FREE_ && outlist[ii].size() > 0.4*avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,cells,ii,pnp1,outlist,dp.norm[ii],dp.avgV[ii]);
				}
				local.emplace_back(temp);
			}

			#pragma omp for schedule(static) ordered
	    	for(int i=0; i<omp_get_num_threads(); i++)
	    	{
	    		#pragma omp ordered
	    		air.insert(air.end(),local.begin(),local.end());
	    	}
		}
	}
	airP.clear();

	// #pragma omp parallel for schedule(static) ordered
	for(size_t ii = 0; ii < air.size(); ++ii)
		for(size_t jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	neighb.clear();
	neighb.reserve(end);

	#pragma omp parallel shared(pnp1, outlist, air)
	{
		std::vector<std::vector<Part>> local;
		if(svar.ghost == 1 )
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				if(air[ii].size()>0)
				{
					temp.insert(temp.end(), air[ii].begin(), air[ii].end());
				}
				local.push_back(temp);
			}
		}
		else
		{
			#pragma omp for schedule(static) nowait
			for (size_t ii = 0; ii < end; ++ii)
			{
				std::vector<Part> temp;
				temp.reserve(outlist[ii].size());
				for(auto jj:outlist[ii])
					temp.push_back(Part(pnp1[jj])); 

				local.push_back(temp);
			}
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),local.begin(),local.end());
    	}
	}


	/*************************** Runge Kutta ***************************************/
	// void(Check_Error(TREE,svar,fvar,start,end,error1,error2,logbase,
	// 						cellsused,outlist,xih,pn,pnp1,k));

	// State st_4 = pn;
	

	// error1 = Runge_Kutta4(TREE,svar,fvar,avar,start,end,cells,cellsused,neighb,outlist,
	// 		dp,logbase,pn,st_2,st_4,airP,Force,dropVel);

	// #pragma omp parallel for shared(res, Rrho, Af/*, wDiff, norm, curve*/)
	// for(size_t ii = 0; ii < end; ++ii)
	// {
	// 	pnp1[ii].f = st_4[ii].f;
	// 	pnp1[ii].Rrho = st_4[ii].Rrho;
	// 	pnp1[ii].rho = st_4[ii].rho;
	// 	pnp1[ii].v = st_4[ii].v;
	// 	pnp1[ii].p = st_4[ii].p;
	// }

	/***************************** Newmark Beta *************************************/
	Newmark_Beta(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,cells,cellsused,neighb,
		outlist,dp,logbase,k,error1,error2,pn,pnp1,airP,Force,dropVel);


	/*Reset positions*/
	for(size_t ii = 0; ii < end; ii++)
	{
		pnp1[ii].xi = xi[ii];
	}

	// for(size_t ii = start; ii < end; ++ii)
	// {
	// 	cout << "L Mag: " << dp.L[ii].determinant() << "  gradRho: " << dp.gradRho[ii].norm() << "  norm: " << dp.norm[ii].norm() 
	// 	<< "  avgV: " << dp.avgV[ii].norm() << "  lam: " << dp.lam[ii] << "  kernsum: " << dp.kernsum[ii] << endl;
	// }

#if DEBUG 	
		dbout << "Exiting first step. Error: " << error1 << endl;
#endif

	if (svar.Bcase == 4)
	{
		// Calculate the force expected for a droplet of the same size.
		svar.Force = Force;

		real radius = 0.5*svar.diam;

		// Average velocity of the particles
		
		StateVecD Vdiff = avar.vInf;

		real Re = avar.rhog*Vdiff.norm()*svar.diam/avar.mug;
		real Cds_undef,Cds_def;

		Cds_undef = GetCd(Re);

		
		// if( Re <= 1000.0)
		// 	Cds = (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
		// else 
		// 	Cds = 0.424;

#if SIMDIM == 3
		real Adrop_und = M_PI*radius*radius;
#else
		real Adrop_und = 2*radius;
#endif
		cout << endl << "Droplet calc parameters: " << endl;
		cout << Re << "  " << Cds_undef << " " << Adrop_und << " " << Vdiff.norm() << "  " << svar.mass << endl;

		// Undeformed drag
		// StateVecD dropForce = (0.5*avar.rhog*Vdiff.norm()*Vdiff*Cds*Adrop)*sqrt(radius);

		// Calculate deformed drag
		AERO bigdrop;
		bigdrop.rhog = avar.rhog;
		bigdrop.mug = avar.mug;
		bigdrop.nfull = avar.nfull;

		// real temp = radius / std::cbrt(3.0/(4.0*M_PI));
		GetYcoef(bigdrop,fvar,svar.diam);
		
		real ymax = Vdiff.squaredNorm()*bigdrop.ycoef;
		if (ymax > 1)
			ymax = 1;

		Re = 2.0*avar.rhog*Vdiff.norm()*bigdrop.L/avar.mug;
		Cds_def = (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re));

		// if( Re <= 1000.0)
		// 	Cds = (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
		// else 
		// 	Cds = 0.424;

		// cout << bigdrop.ycoef << endl;
		// cout << radius << "  "  << bigdrop.L << "  " << ymax << endl;
		// cout << bigdrop.ycoef << endl;
		
		// if (ymax > 1.0)
		// 	ymax = 1.0;


		real const Cdl = Cds_def*(1+2.632*ymax);

		#if SIMDIM == 3 
			real Adrop = M_PI*pow((bigdrop.L + bigdrop.Cb*bigdrop.L*ymax),2);

		#endif
		#if SIMDIM == 2
			real Adrop = 2*(bigdrop.L + bigdrop.Cb*bigdrop.L*ymax);
		#endif

		StateVecD dropDefCd =  0.5*avar.rhog*Vdiff.norm()*Vdiff*Cdl*Adrop;	

		// StateVecD sumForce = StateVecD::Zero();

		// for(size_t ii = 0; ii < end; ++ii)
		// {
		// 	sumForce += pnp1[ii].m*pnp1[ii].f;	
		// }

		Force = StateVecD::Zero();

		for(size_t ii = 0; ii < end; ii++)
		{
			Force += pnp1[ii].Af;
		}

		// real totMass = real(svar.totPts) * pnp1[0].m;

		// cout << totMass << "  " << pnp1[0].m << endl;

		// Force = pnp1[0].m*Force;

		// cout << sumForce[0] << "  " << sumForce[1] << "  " << sumForce[2] << endl;
		// real f_undef = 0.5*avar.rhog*Vdiff.squaredNorm()*Adrop_und*Cds_undef;
		real calcCd = Force[1]/(0.5*avar.rhog*Vdiff.squaredNorm()*Adrop_und);
		real defCd = dropDefCd[1]/(0.5*avar.rhog*Vdiff.squaredNorm()*Adrop_und);
		// StateVecD test = GisslerForce(bigdrop, Vdiff, 1.0, 1.0, 0);
		cout << std::fixed;
		cout << endl << "Drop Cd: " << Cds_undef << "  Drop deformed Cd: "
		 << defCd << "  Calculated Cd: " << calcCd/*Force(1)*/  /*<< " Magnitude: " << Force.norm()*/   << endl;
		cout << "Difference factor: " << (calcCd/Cds_undef-1.0)*100 << "%" << endl << endl;

		cout << std::scientific;
	}
}

#endif
