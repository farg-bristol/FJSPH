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
#include "Shifting.h"
#ifdef RK
#include "Runge_Kutta.h"
#else
#include "Newmark_Beta.h"
#endif



///**************** Integration loop **************///
real Integrate(KDTREE& TREE, SIM& svar, const FLUID& fvar, const AERO& avar, 
	MESH& cells, DELTAP& dp, State& pn, State& pnp1, State& airP, outl& outlist)
{
	// cout << "Entered Newmark_Beta" << endl;
	#ifndef RK
	uint   k = 0;	//iteration number
	real error2 = 1.0;
	#endif
	real logbase = 0.0;
	real error1 = 0.0;
	

	// Find maximum safe timestep
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](Particle const& p1, Particle const& p2){return p1.f.norm()< p2.f.norm();});

	vector<Particle>::iterator maxUi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](Particle const& p1, Particle const& p2){return p1.v.norm()< p2.v.norm();});

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
	#ifdef RK
		svar.dt = 0.125*std::min(dtf,std::min(dtc,dtv));
	#else
		svar.dt = 0.175*std::min(dtf,std::min(dtc,dtv));
	#endif
	/***********************************************************************************/
	/***********************************************************************************/
	/***********************************************************************************/

	if(svar.dt < svar.dt_min)
		svar.dt = svar.dt_min;

	if (svar.dt > (svar.frame+1)*svar.framet-svar.t)
	{
		svar.dt = (svar.frame+1)*svar.framet-svar.t + 1e-10;
	}


	// #ifdef DEBUG
	cout << "time: " << svar.t << " dt: " << svar.dt << "  dtv: " << dtv <<  "  dtf: " << dtf << "  dtc: " << dtc << " Maxf: " << maxf << endl;
	// #endif

	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);

	// Check if the particle has moved to a new cell
	if (svar.Asource == 1 || svar.Asource == 2)
	{
		// cout << "Finding cells" << endl;
		FindCell(svar,fvar.sr,TREE,cells,dp,pn,pnp1);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.simPts -= svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
			dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
		}	
	}
	
	Detect_Surface(svar,fvar,avar,svar.bndPts,svar.totPts,dp,outlist,cells,pnp1);

	if(svar.ghost == 1)
	{	/* Poisson points */
		PoissonGhost(svar,fvar,avar,cells,TREE.NP1,outlist,pn,pnp1);
		dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
	}
	else if (svar.ghost == 2)
	{	/* Lattice points */
		LatticeGhost(svar,fvar,cells,TREE,outlist,pn,pnp1);
		dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
	}

	size_t const start = svar.bndPts;
	size_t end = svar.totPts;
	size_t end_ng = svar.bndPts + svar.simPts;
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

	vector<StateVecD> xih(end-start);
	
	#ifndef RK
	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	#endif

	const real B = fvar.B;
	const real gam = fvar.gam;

	StateVecD dropVel = StateVecD::Zero();
	StateVecD Force = StateVecD::Zero();

	#pragma omp parallel for
	for(size_t ii = start; ii < end; ++ii )
	{
		xih[ii-start] = pnp1[ii].xi;
	}
	
	/*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
	#ifdef RK
	int errstate = (Get_First_RK(TREE, svar, fvar, avar, start, end, B, gam, cells, cellsused,
								 neighb, outlist, dp, logbase, pn, pnp1, error1));

	if (errstate)
		cout << "First step indicates instability. Caution..." << endl;
	#else
		Do_NB_Iter(TREE,svar,fvar,avar,start,end_ng,a,b,c,d,B,gam,
			cells,cellsused,outlist,dp,pn,pnp1,Force,dropVel);
	
		uint nUnstab = 0;

		void(Check_Error(TREE,svar,fvar,start,end_ng,error1,error2,logbase,
								cellsused,outlist,xih,pn,pnp1,k,nUnstab));
		k++; //Update iteration count

	#endif

	// cout << "Error: " << error1 << endl;

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(svar,fvar,end,pnp1,outlist,dp);

	if (svar.Asource == 1 || svar.Asource == 2)
	{
		// cout << "Finding cells" << endl;
		FindCell(svar,fvar.sr,TREE,cells,dp,pn,pnp1);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.simPts -= svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			end = svar.totPts;
			end_ng = svar.bndPts + svar.simPts;
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
			dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
		}	
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);

	// Apply_XSPH(fvar,start,end,outlist,dp,pnp1);
	#ifndef NOALE
		if(svar.ghost > 0)
			Particle_Shift_Ghost(svar,fvar,start,end_ng,outlist,dp,pnp1);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end_ng,outlist,dp,pnp1);
	#endif
	
	dropVel = StateVecD::Zero();

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
	#ifdef RK
		State st_2 = pnp1;

		error1 = Runge_Kutta4(TREE,svar,fvar,avar,start,end,B,gam,cells,cellsused,neighb,outlist,
					dp,logbase,pn,st_2,pnp1,Force,dropVel);
	#else 
		Newmark_Beta(TREE,svar,fvar,avar,start,end_ng,a,b,c,d,B,gam,cells,cellsused,
			outlist,dp,logbase,k,error1,error2,xih,pn,pnp1,Force,dropVel);
	#endif	

	/*Add time to global*/
	svar.t+=svar.dt;

	// cout << "Timestep Params: " << maxf << " " << fvar.Cs + svar.maxmu << " " << dtf << " " << dtcv << endl;
	// cout << "New Time: " << svar.t << endl;

	airP.clear();
	for(size_t ii = end_ng; ii < end; ++ii)
	{
		airP.emplace_back(pnp1[ii]);
	}

	/* Remove the ghost particles */
	// if(svar.ghost != 0)
	// {
	// 	pnp1.erase(pnp1.begin()+svar.bndPts+svar.simPts, pnp1.begin()+svar.bndPts+svar.simPts+svar.gstPts);
	// 	pn.erase(pn.begin()+svar.bndPts+svar.simPts, pn.begin()+svar.bndPts+svar.simPts+svar.gstPts);
	// 	dp.erase(svar.bndPts+svar.simPts, svar.bndPts+svar.simPts+svar.gstPts);
	// 	outlist.erase(outlist.begin()+svar.bndPts+svar.simPts, outlist.begin()+svar.bndPts+svar.simPts+svar.gstPts);
	// 	svar.totPts = svar.bndPts + svar.simPts;
	// 	svar.gstPts = 0;
	// }

	/*Check if more particles need to be created*/
	if(svar.Bcase == 2 || svar.Bcase == 3)
	{
		uint nAdd = 0;
		uint partID = end_ng;
		/* Check if any buffer particles have become pipe particles */
		for (size_t ii = 0; ii < svar.back.size(); ++ii)
		{
			size_t const& pID = svar.back[ii];
			/*Check that the starting area is clear first...*/
			StateVecD const vec = svar.Transp*(pnp1[pID].xi-svar.Start);
			real clear;
			if(svar.Bcase == 2)
				clear =  -(3.0 * svar.Jet(1)-svar.dx) ;
			else
				clear = -(svar.Jet(1)-svar.dx);
			

			if(vec[1] > clear)
			{
				/* Particle has cleared the back zone */
				pnp1[pID].b = PartState.PIPE_;

				/* Update the back vector */
				pnp1[svar.buffer[ii][0]].b = PartState.BACK_;
				svar.back[ii] = svar.buffer[ii][0];

				/* Update the buffer vector */
				svar.buffer[ii][0] = svar.buffer[ii][1];
				svar.buffer[ii][1] = svar.buffer[ii][2];
				svar.buffer[ii][2] = svar.buffer[ii][3];

				/* Create a new particle */
				if(svar.totPts < svar.finPts)
				{
					StateVecD xi = svar.Transp*(pnp1[svar.buffer[ii][3]].xi-svar.Start);

					xi[1] -= svar.dx;
					xi = svar.Rotate*xi + svar.Start;

					pnp1.insert(pnp1.begin() + partID,
					Particle(xi,pnp1[svar.buffer[ii][3]],PartState.BUFFER_,partID));
					svar.buffer[ii][3] = partID;
					if(svar.Asource == 0)
					{
						pnp1[partID].cellRho = pnp1[pID].cellRho;
					}

					svar.simPts++;
					svar.totPts++;
					partID++;
					nAdd++;
				}
			}

		}

		if(nAdd != 0)
		{
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
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
	size_t const start = svar.bndPts;
	size_t end = svar.totPts;
	size_t end_ng = svar.bndPts + svar.simPts;

	#if DEBUG 
		dbout << "Starting first step. ";
		dbout << "  Start index: " << start << "  End index: " << end << endl;
	#endif

	dSPH_PreStep(svar,fvar,end,pnp1,outlist,dp);
	if (svar.Asource == 1 || svar.Asource == 2)
	{
		FindCell(svar,fvar.sr,TREE,cells,dp,pnp1,pn);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			svar.simPts = svar.totPts-svar.bndPts;
			end = svar.totPts;
			end_ng = end;
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
			dSPH_PreStep(svar,fvar,end,pnp1,outlist,dp);
			
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);

	if(svar.ghost == 1)
	{	/* Poisson points */
		PoissonGhost(svar,fvar,avar,cells,TREE.NP1,outlist,pn,pnp1);
		dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
		// Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);
	}
	else if (svar.ghost == 2)
	{	/* Lattice points */
		LatticeGhost(svar,fvar,cells,TREE,outlist,pn,pnp1);
		dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
		// Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);
	}

	end = svar.totPts;

	#ifndef RK
	uint k = 0;	
	real error2 = 0.0;
	#endif
	real logbase = 0.0;
	real error1 = 0.0;
	

	#ifndef RK
	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	#endif

	const real B = fvar.B;
	const real gam = fvar.gam;

	StateVecD dropVel = StateVecD::Zero();
	StateVecD Force = StateVecD::Zero();
	
	/*find force at time n*/
	vector<StateVecD> res;
	vector<StateVecD> Af;
	vector<real> Rrho;
	vector<real> curve;

	
	// #pragma omp parallel for shared(res,Rrho) schedule(static)
	// for(size_t ii = 0; ii < end; ++ii)
	// {
	// 	pn[ii].f = res[ii];
	// 	pn[ii].Rrho = Rrho[ii];
	// }

	// State st_2 = pn;
	vector<size_t> cellsused;

	// Find maximum safe timestep
	// vector<Particle>::iterator maxfi = std::max_element(pnp1.begin() + svar.bndPts, pnp1.end(),
	// 								[](Particle p1, Particle p2) { return p1.f.norm() < p2.f.norm(); });

	// vector<Particle>::iterator maxUi = std::max_element(pnp1.begin() + svar.bndPts, pnp1.end(),
	// 								[](Particle p1, Particle p2) { return p1.v.norm() < p2.v.norm(); });

	// real maxf = maxfi->f.squaredNorm();
	// real maxU = maxUi->v.norm();
	// real dtv = fvar.HSQ * fvar.rho0 / fvar.mu;
	// real dtf = sqrt(fvar.H / maxf);
	// real dtc = fvar.H / (maxU);

	// /***********************************************************************************/
	// /***********************************************************************************/
	// /***********************************************************************************/
	// svar.dt = 0.125 * std::min(dtf, std::min(dtc, dtv));
	// /***********************************************************************************/
	// /***********************************************************************************/
	// /***********************************************************************************/

	// if (svar.dt > (svar.frame + 1) * svar.framet - svar.t)
	// {
	// 	svar.dt = (svar.frame + 1) * svar.framet - svar.t;
	// }

	// cout << "time: " << svar.t << " dt: " << svar.dt << "  dtv: " << dtv << "  dtf: " << dtf << "  dtc: " << dtc << " Maxf: " << maxf << endl;

	// /*Previous State for error calc*/
	vector<StateVecD> xih(end-start);
	// vector<StateVecD> xi(end); /*Keep original positions to reset after finding forces*/
	#pragma omp parallel for shared(pnp1)
	for (size_t  ii=start; ii < end; ++ii)
	{
		xih[ii-start] = pnp1[ii].xi;
		// xi[ii] = pnp1[ii].xi;
	}

	/*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
	#ifdef RK
		int errstate = (Get_First_RK(TREE, svar, fvar, avar, start, end, B, gam, cells, cellsused,
									neighb, outlist, dp, logbase, pn, pnp1, error1));

		if(errstate)
			cout << "First step indicates instability. Caution..." << endl;
	#else

		Do_NB_Iter(TREE, svar, fvar, avar, start, end_ng, a, b, c, d, B, gam,
				cells, cellsused, outlist, dp, pn, pnp1, Force, dropVel);

		uint nUnstab = 0;

		void(Check_Error(TREE, svar, fvar, start, end_ng, error1, error2, logbase,
						cellsused, outlist, xih, pn, pnp1, k, nUnstab));
		k++; //Update iteration count

	#endif

	/*Find maximum safe timestep*/
	// maxUi = std::max_element(pnp1.begin(),pnp1.end(),
	// 	[](Particle p1, Particle p2){return p1.v.norm()< p2.v.norm();});

	// maxU = maxUi->v.norm();
	// dtc = fvar.H/(maxU);

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	// svar.dt = 0.125*std::min(dtf,std::min(dtc,dtv));
	
	// if (svar.dt > (svar.frame + 1) * svar.framet - svar.t)
	// {
	// 	svar.dt = (svar.frame + 1) * svar.framet - svar.t;
	// }

	// dt = dt;
	// dt2 = dt * dt;

	// cout << "time: " << svar.t << " dt: " << svar.dt << "  dtv: " << dtv << "  dtf: " << dtf << "  dtc: " << dtc << " Maxf: " << maxf << endl;
	/***********************************************************************************/
	/***********************************************************************************/
	/***********************************************************************************/

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(svar,fvar,end,pnp1,outlist,dp);

	if (svar.Asource == 1 || svar.Asource == 2)
	{
		FindCell(svar,fvar.sr,TREE,cells,dp,pnp1,pn);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			svar.simPts = svar.totPts-svar.bndPts;
			end = svar.totPts;
			end_ng = end;
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
			dSPH_PreStep(svar,fvar,end,pnp1,outlist,dp);
			
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);
	
	/*Do time integration*/
	#ifdef RK
		State st_2 = pnp1;

		error1 = Runge_Kutta4(TREE, svar, fvar, avar, start, end, B, gam, cells, cellsused, neighb, outlist,
							  dp, logbase, pn, st_2, pnp1, Force, dropVel);
	#else

		Newmark_Beta(TREE, svar, fvar, avar, start, end_ng, a, b, c, d, B, gam, cells, cellsused,
				outlist, dp, logbase, k, error1, error2, xih, pn, pnp1, Force, dropVel);
	#endif

	/*Reset positions*/
	// for(size_t ii = 0; ii < end; ii++)
	// {
	// 	pnp1[ii].xi = xi[ii];
	// }

	// for(size_t ii = start; ii < end; ++ii)
	// {
	// 	cout << "L Mag: " << dp.L[ii].determinant() << "  gradRho: " << dp.gradRho[ii].norm() << "  norm: " << dp.norm[ii].norm() 
	// 	<< "  avgV: " << dp.avgV[ii].norm() << "  lam: " << dp.lam[ii] << "  kernsum: " << dp.kernsum[ii] << endl;
	// }

	airP.clear();
	for(size_t ii = end_ng; ii < end; ++ii)
	{
		airP.emplace_back(pnp1[ii]);
	}

	/* Remove the ghost particles */
	if(svar.ghost != 0)
	{
		pnp1.erase(pnp1.begin()+svar.bndPts+svar.simPts, pnp1.begin()+svar.bndPts+svar.simPts+svar.gstPts);
		pn.erase(pn.begin()+svar.bndPts+svar.simPts, pn.begin()+svar.bndPts+svar.simPts+svar.gstPts);
		dp.erase(svar.bndPts+svar.simPts, svar.bndPts+svar.simPts+svar.gstPts);
		outlist.erase(outlist.begin()+svar.bndPts+svar.simPts, outlist.begin()+svar.bndPts+svar.simPts+svar.gstPts);
		svar.totPts = svar.bndPts + svar.simPts;
		svar.gstPts = 0;

	}

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
