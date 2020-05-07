/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef NEWMARK_BETA_H
#define NEWMARK_BETA_H


#include <chrono>
#include <unordered_set>
#include "Var.h"
#include "IO.h"
#include "Neighbours.h"
#include "Resid.h"
#include "Crossing.h"

///**************** Integration loop **************///
real Newmark_Beta(KDTREE& TREE, SIM& svar, const FLUID& fvar, const AERO& avar, 
	MESH& cells, State& pn, State& pnp1, State& airP, outl& outlist)
{
	// cout << "Entered Newmark_Beta" << endl;
	uint   k = 0;	
	real errsum = 1.0;
	real logbase = 0.0;
	real error1 = 0.0;
	real error2 = 0.0;

	// Find maximum safe timestep
	vector<Particle>::iterator maxfi = std::max_element(pnp1.begin(),pnp1.end(),
		[](Particle p1, Particle p2){return p1.f.norm()< p2.f.norm();});
	real maxf = maxfi->f.norm();
	real dtf = sqrt(fvar.H/maxf);
	real dtcv = fvar.H/(fvar.Cs+svar.maxmu);

	// Only use if -fno-finite-math-only is on
	// if (std::isinf(maxf))
	// {
	// 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
	// 	exit(-1);
	// }

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
	svar.dt = 0.6*std::min(dtf,dtcv);
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
#ifdef DEBUG
	cout << "time: " << svar.t << " dt: " << svar.dt << "  dtf: " << dtf << "  dtcv: " << dtcv << " maxmu: " << svar.maxmu 
	<< " Maxf: " << maxf << endl;
#endif

	// cout << svar.framet * (svar.frame+1) << "  " << svar.t  << endl;

	if (svar.dt > svar.framet)
	{
		svar.dt = svar.framet;
	}
	else if (svar.dt > (svar.frame+1)*svar.framet-svar.t)
	{
		svar.dt = (svar.frame+1)*svar.framet-svar.t;
		// cout << "New dt: " << svar.dt << endl;
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
				if(pnp1[ii].b == FREE)
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
				if(pnp1[ii].b == FREE && outlist[ii].size() > 0.4*avar.nfull 
					&& outlist[ii].size() < avar.nfull)
				{
					temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,ii,pnp1,outlist);
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


	while (log10(sqrt(errsum/(real(svar.totPts)))) - logbase > -7.0)
	{
		/****** UPDATE TREE ***********/
		TREE.NP1.index->buildIndex();
		FindNeighbours(TREE.NP1, fvar, pnp1, outlist);
		
		// airP.clear();
		std::vector<std::vector<Part>> neighb;
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
				
		// cout << "K: " << k << endl;
		// cout << "Calculating forces" << endl;
		vector<StateVecD> res(end,StateVecD::Zero());
		vector<StateVecD> Af(end,StateVecD::Zero());
		vector<real> Rrho(end,0.0);
		// vector<real> wDiff(end,0.0);
		// vector<StateVecD> norm(end, StateVecD::Zero());
		// vector<real> curve(end,0.0);

		// cout << "Neighbour list size: " << neighb.size() << "  Outlist size: " << outlist.size() << endl;
		// cout << "Start: " << start << "  End: " << end << endl;
		Force = StateVecD::Zero();
		svar.AForce = StateVecD::Zero();
 		Forces(svar,fvar,avar,cells,pnp1,neighb,outlist,res,Rrho,Af/*,wDiff,norm,curve*/); /*Guess force at time n+1*/

	
		/*Previous State for error calc*/
		#pragma omp parallel for shared(pnp1)
		for (size_t  ii=0; ii < end; ++ii)
			xih[ii] = pnp1[ii].xi;

		const real dt = svar.dt;
		const real dt2 = dt*dt;

		/*Update the state at time n+1*/
		#pragma omp parallel shared(pn,res,Rrho,svar,fvar) /*reduction(+:Force,dropVel)*/
		{
			if(svar.Asource == 2)
			{
				#pragma omp for schedule(static) nowait
				for(size_t const& ii : cellsused)
				{
						cells.fNum[ii] = 0;
						cells.fMass[ii] = 0.0;
						cells.vFn[ii] = StateVecD::Zero();
						cells.vFnp1[ii] = StateVecD::Zero();
				}
			}
			// #pragma omp for schedule(static) nowait
			// for (size_t ii = 0; ii < piston; ++ii)
			// {	/****** PISTON PARTICLES *************/
			// 	pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+Rrho[ii]);
			// 	pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
			// 	pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
			// 	pnp1[ii].f = res[ii];
			// 	pnp1[ii].Rrho = Rrho[ii];
			// }

			#pragma omp for schedule(static) nowait
			for (size_t ii=0/*piston*/; ii < start; ++ii)
			{	/****** BOUNDARY PARTICLES ***********/
				pnp1[ii].rho = pn[ii].rho+0.5*dt*(pn[ii].Rrho+Rrho[ii]);
				pnp1[ii].p = B*(pow(pnp1[ii].rho/fvar.rho0,gam)-1);
				// pnp1[ii].p = fvar.pPress;
				// pnp1[ii].rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);
				pnp1[ii].f = res[ii];
				pnp1[ii].Rrho = Rrho[ii];
			}


			#pragma omp for schedule(static) nowait 
			for (size_t ii=start; ii < end; ++ii)
			{	/****** FLUID PARTICLES **************/
				/*Check if the particle is clear of the starting area*/
				if(pnp1[ii].b == START || pnp1[ii].b == BACK)
				{   /* START = pipe particle receiving prescribed motion            */
					/* BACK = the latest particle in that column                    */
					/* PIPE = in the pipe, with free motion                         */
					/* FREE = free of the pipe and receives an aero force           */
					StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
					if(vec(1) > svar.clear)
					{	/*Tag it as clear if it's higher than the plane of the exit*/
						pnp1[ii].b=PIPE;
					}
					else
					{	/*For the particles marked 1, perform a prescribed motion*/
						pnp1[ii].xi = pn[ii].xi + dt*pn[ii].v;
						// pnp1[ii].f = res[ii];
						// pnp1[ii].Rrho = Rrho[ii];
						// pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*pnp1[ii].Rrho);
						// pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					}
				}

				if(pnp1[ii].b == PIPE)
				{	/*Do a check to see if it needs to be given an aero force*/
					StateVecD vec = svar.Transp*(pnp1[ii].xi-svar.Start);
					if(vec(1) > 0.0)
					{
						pnp1[ii].b = FREE;
						if(svar.Asource == 1 || svar.Asource == 2)
						{
							/*retrieve the cell it's in*/
							FirstCell(svar,end,ii,TREE.CELL, cells, pnp1, pn);
						}
					}	
				}

				if(pnp1[ii].b > START)
				{	/*For any other particles, intergrate as normal*/
					pnp1[ii].xi = pn[ii].xi+dt*pn[ii].v+dt2*(c*pn[ii].f+d*res[ii]);
					pnp1[ii].v =  pn[ii].v +dt*(a*pn[ii].f+b*res[ii]);
					pnp1[ii].f = res[ii];
					pnp1[ii].Af = Af[ii];
					pnp1[ii].Rrho = Rrho[ii];
					
					
					pnp1[ii].rho = pn[ii].rho+dt*(a*pn[ii].Rrho+b*Rrho[ii]);
					pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					
					
					Force += res[ii]*pnp1[ii].m;
					dropVel += pnp1[ii].v;

					
					// pnp1[ii].theta = wDiff[ii];
					// pnp1[ii].nNeigb = real(outlist[ii].size());

// #if SIMDIM == 2
// 					if(real(outlist[ii].size()) < 5.26 * wDiff[ii] +30)
// 					{
// 						pnp1[ii].normal = norm[ii];
// 						pnp1[ii].surf = curve[ii];
// 						pnp1[ii].theta = 1.0;
// 					}
// 					else
// 					{
// 						pnp1[ii].normal = StateVecD::Zero();
// 						pnp1[ii].surf = 0.0;
// 						pnp1[ii].theta = 0.0;
// 					}
// #else
// 					if(real(outlist[ii].size()) < 297.25 * wDiff[ii] + 34)
// 					{
// 						pnp1[ii].theta = wDiff[ii];
// 						pnp1[ii].normal = norm[ii];
// 						pnp1[ii].surf = curve[ii];
// 					}
// 					else
// 					{
// 						pnp1[ii].normal = StateVecD::Zero();
// 						pnp1[ii].surf = 0.0;
// 						pnp1[ii].theta = 0.0;
// 					}
// #endif

					
					if(svar.Asource == 2 && pnp1[ii].b == FREE)
					{
						#pragma omp atomic
						cells.fNum[pnp1[ii].cellID]++;
						#pragma omp atomic
						cells.fMass[pnp1[ii].cellID] += pnp1[ii].m;

						#pragma omp critical
						{
							cells.vFn[pnp1[ii].cellID] += pn[ii].v;
							cells.vFnp1[pnp1[ii].cellID] += pnp1[ii].v;
						}	
					}
					

				}
				
			} /*End fluid particles*/


			if(svar.Asource == 2)
			{
				#pragma omp for schedule(static) nowait
				for(size_t const& ii : cellsused)
				{
					// Work out the mass and volume fractions
					if(cells.fNum[ii] != 0)
					{
						real fVol = real(cells.fNum[ii]) * avar.pVol;

						real aFrac = (cells.cVol[ii]-fVol)/cells.cVol[ii];

						
						if (aFrac < 0.1)
							continue;
						
						real aMass = cells.cMass[ii]*aFrac;

						// Do the momentum exchange
						StateVecD newPert = (cells.fMass[ii]/aMass)*(cells.vFnp1[ii]-cells.vFn[ii])/real(cells.fNum[ii]);
						// StateVecD diffusion = 0.2*cells.cPertn[ii];
						cells.cPertnp1[ii] = cells.cPertn[ii]*0.9 - newPert;

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

		}/*End pragma omp parallel*/

		/****** FIND ERROR ***********/
		errsum = 0.0;
		#pragma omp parallel for reduction(+:errsum) schedule(static) shared(pnp1,xih)
		for (size_t ii = start; ii < end; ++ii)
		{
			StateVecD r = pnp1[ii].xi-xih[ii];
			errsum += r.squaredNorm();
		}
	
		if(k == 0)
			logbase=log10(sqrt(errsum/(real(svar.totPts))));


		error1 = log10(sqrt(errsum/(real(svar.totPts)))) - logbase;
		// cout << RestartCount << "  " << k << "  " << error1  << "  " << svar.dt << endl;
		// cout << k << "  " << error1 << "  " << svar.dt << endl;

		if (error1-error2 > 0.0 /*|| std::isnan(error1)*/)
		{	/*If simulation starts diverging, then reduce the timestep and try again.*/
			
			pnp1 = pn;
			
			TREE.NP1.index->buildIndex();
			FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

			if(svar.Asource == 2)
			{
				cellsused.clear();
				#pragma omp parallel shared(pnp1)
				{
					std::vector<size_t> local;
					#pragma omp for schedule(static) nowait
					for (size_t ii = start; ii < end; ++ii)
					{	
						if(pnp1[ii].b == FREE)
							local.emplace_back(pnp1[ii].cellID);
					}

					#pragma omp for schedule(static) ordered
					for(int i=0; i<omp_get_num_threads(); i++)
					{
						#pragma omp ordered
						cellsused.insert(cellsused.end(),local.begin(),local.end());
					}
				}

				// Sort the vector and remove unique values.
				std::unordered_set<size_t> s;
				for(size_t const& ii:cellsused)
					s.insert(ii);

				cellsused.assign(s.begin(),s.end());
				std::sort(cellsused.begin(),cellsused.end());
			}


			

			svar.dt = 0.1*svar.dt;
			// cout << "Unstable timestep. New dt: " << svar.dt << endl;
			k = 0;
			error1 = 0.0;
			// RestartCount++;
		}	/*Check if we've exceeded the maximum iteration count*/
		else if (k > svar.subits)
			break;
		else
		{	/*Otherwise, roll forwards*/
			++k;
		}

		error2 = error1;
		
	} /*End of subits*/
	// cout << "Out of the parallel loop." << endl;
	// RestartCount = 0;
	/*Add time to global*/

	svar.t+=svar.dt;
	// cout << "Timestep Params: " << maxf << " " << fvar.Cs + svar.maxmu << " " << dtf << " " << dtcv << endl;
	// cout << "New Time: " << svar.t << endl;

	/*Check if more particles need to be created*/
	if(svar.Bcase == 2 || svar.Bcase == 3)
	{
		if(svar.Bclosed == 0)
		{
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

						pn.emplace_back(Particle(xi,pnp1[back],svar.totPts));
						pnp1.emplace_back(Particle(xi,pnp1[back],svar.totPts));

						pnp1[back].b = START;
						/*Update the back vector*/
						back = svar.totPts;
						svar.simPts++;
						svar.totPts++;
					}
				}
				else
				{
					svar.Bclosed = 1;
				}	
			}
			
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
			aMomT += aMass*(cells.cVel[ii]+cells.cPertnp1[ii]).norm();
			tMom += aMass*(cells.cVel[ii]+cells.cPertnp1[ii]).norm();
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
		pertLog << svar.t << " " << tMom-(svar.tMom) << " " << aMomT-svar.aMom << " " << fMom << endl;
	}
	else if (svar.Bcase == 4)
	{
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
	}

	return log10(sqrt(errsum/(real(svar.totPts))))-logbase;
}



void First_Step(KDTREE& TREE, SIM& svar, const FLUID& fvar, const AERO& avar, 
	MESH& cells, outl& outlist, State& pnp1, State& pn, State& airP)
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


	/*Check if a particle is running low on neighbours, and add ficticious particles*/
	vector<vector<Part>> air;
	#pragma omp parallel shared(svar, pnp1, outlist)
	{
		std::vector<std::vector<Part>> localN;
		vector<vector<Part>> localA;
		#pragma omp for schedule(static) nowait 
		for (size_t ii = start; ii < end; ++ii)
		{
			std::vector<Part> temp;
			if(svar.ghost == 1 && pnp1[ii].b == FREE && outlist[ii].size() < avar.nfull &&
				outlist[ii].size() > 0.4*avar.nfull)
				temp = PoissonSample::generatePoissonPoints(svar,fvar,avar,ii,pnp1,outlist);

			localA.emplace_back(temp);

			for(auto j:outlist[ii])
				temp.emplace_back(Part(pnp1[j]));

			localN.emplace_back(temp);
		}

		#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++)
    	{
    		#pragma omp ordered
    		neighb.insert(neighb.end(),localN.begin(),localN.end());
    		air.insert(air.end(),localA.begin(),localA.end());
    	}
	}
	airP.clear();

	for(size_t ii = 0; ii < air.size(); ++ii)
		for(size_t jj = 0; jj < air[ii].size(); ++jj)
			airP.emplace_back(PartToParticle(air[ii][jj]));

	#pragma omp parallel for shared(outlist)
	for(size_t ii = start; ii < end; ++ii)
	{
		pnp1[ii].theta = outlist[ii].size(); 
	}

	/*Previous State for error calc*/
	vector<StateVecD> xih(svar.totPts);
	#pragma omp parallel for shared(pnp1)
	for (size_t  ii=0; ii < end; ++ii)
		xih[ii] = pnp1[ii].xi;
	
	// vector<StateVecD> vPert(end,StateVecD::Zero());
	vector<StateVecD> res(end,StateVecD::Zero());
	vector<StateVecD> Af(end,StateVecD::Zero());
	vector<real> Rrho(svar.totPts,0.0);
	// vector<real> wDiff(end,0.0);
	// vector<StateVecD> norm(end,StateVecD::Zero());
	// vector<real> curve(end,0.0);

	Forces(svar,fvar,avar,cells,pnp1,neighb,outlist,res,Rrho,Af/*,wDiff,norm,curve*/); /*Guess force at time n+1*/

	/*Find maximum safe timestep*/
	vector<StateVecD>::iterator maxfi = std::max_element(res.begin(),res.end(),
		[](StateVecD p1, StateVecD p2){return p1.norm() < p2.norm();});
	real maxf = maxfi->norm();
	real dtf = sqrt(fvar.H/maxf);
	real dtcv = fvar.H/(fvar.Cs+svar.maxmu);
	const real dt = 0.3*std::min(dtf,dtcv);

	#pragma omp parallel for shared(res, Rrho, Af/*, wDiff, norm, curve*/)
	for(size_t ii = 0; ii < end; ++ii)
	{
		pnp1[ii].f = res[ii];
		pnp1[ii].Rrho = Rrho[ii];
		pnp1[ii].Af = Af[ii]; 
		pnp1[ii].rho = pn[ii].rho+dt*(svar.gamma*Rrho[ii]);
		pnp1[ii].p = fvar.B*(pow(pnp1[ii].rho/fvar.rho0,fvar.gam)-1);
					
		xih[ii] = pnp1[ii].xi + dt*pnp1[ii].v;

		// pnp1[ii].theta = wDiff[ii];
		pnp1[ii].nNeigb = real(outlist[ii].size());
		


// #if SIMDIM == 2
// 		if(real(outlist[ii].size()) < 5.26 * wDiff[ii] +30)
// 		{
// 			pnp1[ii].normal = norm[ii];
// 			pnp1[ii].surf = curve[ii];
// 			pnp1[ii].theta = 1.0;
// 		}
// 		else
// 		{
// 			pnp1[ii].normal = StateVecD::Zero();
// 			pnp1[ii].surf = 0.0;
// 			pnp1[ii].theta = 0.0;
// 		}
// #else
// 		if(real(outlist[ii].size()) < 297.25 * wDiff[ii] + 34)
// 		{
// 			pnp1[ii].theta = wDiff[ii];
// 			pnp1[ii].normal = norm[ii];
// 			pnp1[ii].surf = curve[ii];
// 		}
// 		else
// 		{
// 			pnp1[ii].normal = StateVecD::Zero();
// 			pnp1[ii].surf = 0.0;
// 			pnp1[ii].theta = 0.0;
// 		}
// #endif
	}


#if DEBUG 
	real errsum = 0.0;

	for (size_t ii = start; ii < end; ++ii)
	{
		StateVecD r = xih[ii]-pnp1[ii].xi;
		errsum += r.squaredNorm();
	}

	real error=log10(sqrt(errsum/(real(svar.totPts))));

	
		dbout << "Exiting first step. Error: " << error << endl;
#endif
}

#endif