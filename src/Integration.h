/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <chrono>

#include "Var.h"
#include "IO.h"
#include "Neighbours.h"
#include "Resid.h"
#include "Containment.h"
#include "Shifting.h"
#include "IPT.h"
#ifdef RK
#include "Runge_Kutta.h"
#else
#include "Newmark_Beta.h"
#endif



///**************** Integration loop **************///
real Integrate(KDTREE& TREE, SIM& svar, const FLUID& fvar, const AERO& avar, 
	MESH& cells, SURFS& surf_marks, DELTAP& dp, SPHState& pn, SPHState& pnp1, vector<IPTState>& iptdata, OUTL& outlist)
{
	size_t const start = svar.bndPts;
	size_t end_ng = svar.bndPts + svar.simPts;
	
	// cout << "Entered Newmark_Beta" << endl;
	#ifndef RK
	uint   k = 0;	//iteration number
	real error2 = 1.0;
	#endif
	real logbase = 0.0;
	real error1 = 0.0;
	

	// Find maximum safe timestep
	vector<SPHPart>::iterator maxfi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](SPHPart const& p1, SPHPart const& p2){return p1.acc.norm()< p2.acc.norm();});

	vector<SPHPart>::iterator maxUi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](SPHPart const& p1, SPHPart const& p2){return p1.v.norm()< p2.v.norm();});

	real maxf = maxfi->acc.norm();
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

	if(svar.dt < svar.dt_min)
		svar.dt = svar.dt_min;
	else if(svar.dt > svar.dt_max)
		svar.dt = svar.dt_max;

	if (svar.dt > (svar.frame+1)*svar.framet-svar.t)
		svar.dt = (svar.frame+1)*svar.framet-svar.t + svar.dt_min;

	#ifdef DEBUG
	vector<SPHPart>::iterator maxAfi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
	[](SPHPart const& p1, SPHPart const& p2){return p1.Af.norm()< p2.Af.norm();});
	real maxAf = maxAfi->Af.norm();
	cout << "time: " << svar.t << " dt: " << svar.dt <<  "  dtf: " << dtf << "  dtc: " << dtc << " Maxf: " << maxf << " MaxAf: " << maxAf << endl;
	#endif

	


	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);

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
			dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);
		}	
	}
	#if SIMDIM == 3
	else if (svar.Asource == 3)
	{
		/* Find the VLM influence for the particles */
		for(size_t ii = start; ii < end_ng ; ii++)
		{
			if(dp.lam_ng[ii] < 0.75  && pnp1[ii].b == FREE)
			{
				pnp1[ii].cellV = svar.vortex.getVelocity(pnp1[ii].xi);
				pnp1[ii].cellID = 1;
			}
			else
				pnp1[ii].cellID = 0;
		}
	}
	#endif
	else
	{
		for(size_t ii = start; ii < end_ng ; ii++)
		{
			if(dp.lam_ng[ii] < 0.75  && pnp1[ii].b == FREE)
			{
				pnp1[ii].cellV = avar.vInf;
				pnp1[ii].cellID = 1;
			}
			else
				pnp1[ii].cellID = 0;
		}
	}

	Detect_Surface(svar,fvar,avar,svar.bndPts,svar.totPts,dp,outlist,cells,pnp1);

	if(svar.ghost == 1)
	{	/* Poisson points */
		PoissonGhost(svar,fvar,avar,cells,TREE.NP1,outlist,pn,pnp1);
		dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);
	}
	else if (svar.ghost == 2)
	{	/* Lattice points */
		LatticeGhost(svar,fvar,avar,cells,TREE,outlist,pn,pnp1);
		dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);
	}

	size_t end = svar.totPts;

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,dp,pnp1);
	#endif

	#ifdef ALE
		if(svar.ghost > 0)
			Particle_Shift_Ghost(svar,fvar,start,end_ng,outlist,dp,pnp1);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end_ng,outlist,dp,pnp1);
	#endif

	vector<size_t> cellsused; // Cells that contain a particle
	if(svar.Asource == 2)
	{
		vector<size_t> tempcell;
		#pragma omp parallel default(shared) // shared(pnp1)
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

	#pragma omp parallel for default(shared)
	for(size_t ii = start; ii < end; ++ii )
	{
		xih[ii-start] = pnp1[ii].xi;
	}
	
	/*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
	#ifdef RK
		error1 = (Get_First_RK(TREE, svar, fvar, avar, start, end, B, gam, cells, cellsused,
								 outlist, dp, logbase, pn, pnp1, error1));


	#else
		Do_NB_Iter(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,
			cells,cellsused,outlist,dp,pn,pnp1,Force,dropVel);
	
		uint nUnstab = 0;

		void(Check_Error(TREE,svar,fvar,start,end,error1,error2,logbase,
								cellsused,outlist,xih,pn,pnp1,k,nUnstab));
		k++; //Update iteration count

	#endif

	// cout << "Error: " << error1 << endl;

	/****** UPDATE TREE ***********/
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	dSPH_PreStep(fvar,end,pnp1,outlist,dp);

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
			dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);
			
		}	
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,dp,pnp1);
	#endif
	
	// Apply_XSPH(fvar,start,end,outlist,dp,pnp1);
	
	#ifdef ALE
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
		#pragma omp parallel default(shared) // shared(pnp1)
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

		#pragma omp parallel for schedule(static) default(shared)
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
		SPHState st_2 = pnp1;

		error1 = Runge_Kutta4(TREE,svar,fvar,avar,start,end,B,gam,cells,cellsused,outlist,
					dp,logbase,pn,st_2,pnp1,Force,dropVel);
	#else 
		Newmark_Beta(TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,cells,cellsused,
			outlist,dp,logbase,k,error1,error2,xih,pn,pnp1,Force,dropVel);
	#endif	

	

	// cout << "Timestep Params: " << maxf << " " << fvar.Cs + svar.maxmu << " " << dtf << " " << dtcv << endl;
	// cout << "New Time: " << svar.t << endl;

	// airP.clear();
	// for(size_t ii = end_ng; ii < end; ++ii)
	// {
	// 	airP.emplace_back(pnp1[ii]);
	// }

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

	if(svar.ghost != 0)
		Check_If_Ghost_Needs_Removing(svar,fvar,TREE.NP1,pn,pnp1);

	/*Check if more particles need to be created*/
	if(svar.Scase == 4)
	{
		uint nAdd = 0;
		size_t& partID = svar.partID;
		/* Check if any buffer particles have become pipe particles */
		for (size_t ii = 0; ii < svar.back.size(); ++ii)
		{
			size_t const& pID = svar.back[ii];
			/*Check that the starting area is clear first...*/
			StateVecD const vec = svar.Transp*(pnp1[pID].xi-svar.sim_start);
			real clear;
			if(svar.Bcase == 3)
				clear =  -(3.0 * svar.jet_depth-svar.dx) ;
			else
				clear = -(svar.jet_depth-svar.dx);
			

			if(vec[1] > clear)
			{
				/* particle has cleared the back zone */
				pnp1[pID].b = PIPE;

				/* Update the back vector */
				pnp1[svar.buffer[ii][0]].b = BACK;
				svar.back[ii] = svar.buffer[ii][0];

				/* Update the buffer vector */
				svar.buffer[ii][0] = svar.buffer[ii][1];
				svar.buffer[ii][1] = svar.buffer[ii][2];
				svar.buffer[ii][2] = svar.buffer[ii][3];

				/* Create a new particle */
				if(svar.totPts < svar.finPts)
				{
					StateVecD xi = svar.Transp*(pnp1[svar.buffer[ii][3]].xi-svar.sim_start);

					xi[1] -= svar.dx;
					xi = svar.Rotate*xi + svar.sim_start;
					pnp1.insert(pnp1.begin() + end_ng,
					SPHPart(xi,pnp1[svar.buffer[ii][3]],BUFFER,partID));
					svar.buffer[ii][3] = end_ng;
					if(svar.Asource == 0 || svar.Asource == 3)
					{
						pnp1[partID].cellRho = pnp1[pID].cellRho;
						pnp1[partID].cellP = pnp1[pID].cellP;
					}

					svar.simPts++;
					svar.totPts++;
					end_ng++;
					end++;
					partID++;
					nAdd++;
				}
			}

		}

		/* Check if any particles need to be deleted, and turned to particle tracking */
		
		uint nDel = 0;
		if(svar.using_ipt && svar.Asource != 0)
		{
			vector<size_t> to_del;
			vector<IPTPart> IPT_nm1, IPT_n, IPT_np1;
			#pragma omp parallel for default(shared)
			for(size_t ii = 0; ii < end_ng; ii++)
			{
				if(pnp1[ii].xi[0] > svar.max_x_sph)
				{
					/* SPHPart is downstream enough to convert to IPT */
					#pragma omp critical
					{
						to_del.emplace_back(ii);
						IPT_nm1.emplace_back(IPTPart(pnp1[ii],svar.t,svar.IPT_diam,svar.IPT_area));
					}
				}
			}

			/* Delete the old particles */
			std::sort(to_del.begin(), to_del.end());
			for(vector<size_t>::reverse_iterator itr = to_del.rbegin(); itr!=to_del.rend(); ++itr)
			{
				pnp1.erase(pnp1.begin() + *itr);
				svar.totPts--;
				svar.simPts--;
				end_ng--;
				end--;
				nDel++;
			}

			/* Need to shift back vector and buffer vector for correct inlet */
			for(size_t& back:svar.back)
			{
				back -= nDel;
			}

			for(vector<size_t>& buffer:svar.buffer)
				for(size_t& part:buffer)
				{
					part -= nDel;
				}


			/* Do particle tracking on the particles converted */
		
			IPT_n = IPT_nm1;
			IPT_np1 = IPT_nm1;
			#pragma omp parallel for default(shared)
			for(size_t ii = 0; ii < IPT_np1.size(); ++ii)
			{
				IPT::Integrate(svar, fvar,avar, cells,ii,IPT_nm1[ii],IPT_n[ii],IPT_np1[ii]
							,surf_marks,iptdata);
			}
			
		}

		if(nAdd != 0 || nDel != 0)
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

		#pragma omp parallel for schedule(static) default(shared)
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

		for(size_t jj = 0; jj < end_ng; ++jj)
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

	/*Add time to global*/
	svar.t+=svar.dt;

	/* Find max x and y coordinates of the fluid */
		// Find maximum safe timestep
	vector<SPHPart>::iterator maxXi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](SPHPart const& p1, SPHPart const& p2){return p1.xi[0] < p2.xi[0];});

	vector<SPHPart>::iterator maxYi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](SPHPart const& p1, SPHPart const& p2){return p1.xi[1] < p2.xi[1];});

	real maxX = maxXi->xi[0];
	real maxY = maxYi->xi[1];

	dambreak << svar.t << "  " << maxX << "  " << maxY << endl;

	return error1;
}

void First_Step(KDTREE& TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	MESH& cells, OUTL& outlist, DELTAP& dp, SPHState& pnp1, SPHState& pn, vector<IPTState>& iptdata)
{
	size_t const start = svar.bndPts;
	size_t end = svar.totPts;
	size_t end_ng = svar.bndPts + svar.simPts;

	#if DEBUG 
		dbout << "Starting first step. ";
		dbout << "  Start index: " << start << "  End index: " << end << endl;
	#endif

	dSPH_PreStep(fvar,end,pnp1,outlist,dp);
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
			dSPH_PreStep(fvar,end,pnp1,outlist,dp);
			
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);

	if(svar.ghost == 1)
	{	/* Poisson points */
		PoissonGhost(svar,fvar,avar,cells,TREE.NP1,outlist,pn,pnp1);
		dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);
		// Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);
	}
	else if (svar.ghost == 2)
	{	/* Lattice points */
		LatticeGhost(svar,fvar,avar,cells,TREE,outlist,pn,pnp1);
		dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,dp);
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
	// 	pn[ii].acc = res[ii];
	// 	pn[ii].Rrho = Rrho[ii];
	// }

	// SPHState st_2 = pn;
	vector<size_t> cellsused;

	// Find maximum safe timestep
	// vector<SPHPart>::iterator maxfi = std::max_element(pnp1.begin() + svar.bndPts, pnp1.end(),
	// 								[](SPHPart p1, SPHPart p2) { return p1.acc.norm() < p2.acc.norm(); });

	// vector<SPHPart>::iterator maxUi = std::max_element(pnp1.begin() + svar.bndPts, pnp1.end(),
	// 								[](SPHPart p1, SPHPart p2) { return p1.v.norm() < p2.v.norm(); });

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

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,dp,pnp1);
	#endif

	// /*Previous SPHState for error calc*/
	vector<StateVecD> xih(end-start);
	// vector<StateVecD> xi(end); /*Keep original positions to reset after finding forces*/
	#pragma omp parallel for default(shared) // shared(pnp1)
	for (size_t  ii=start; ii < end; ++ii)
	{
		xih[ii-start] = pnp1[ii].xi;
		// xi[ii] = pnp1[ii].xi;
	}

	/*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
	#ifdef RK
		error1 = (Get_First_RK(TREE, svar, fvar, avar, start, end, B, gam, cells, cellsused,
									 outlist, dp, logbase, pn, pnp1, error1));
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
	// 	[](SPHPart p1, SPHPart p2){return p1.v.norm()< p2.v.norm();});

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

	dSPH_PreStep(fvar,end,pnp1,outlist,dp);

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
			dSPH_PreStep(fvar,end,pnp1,outlist,dp);
			
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,dp,outlist,cells,pnp1);

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,dp,pnp1);
	#endif
	
	/*Do time integration*/
	#ifdef RK
		SPHState st_2 = pnp1;

		error1 = Runge_Kutta4(TREE, svar, fvar, avar, start, end, B, gam, cells, cellsused, outlist,
							  dp, logbase, pn, st_2, pnp1, Force, dropVel);
	#else
		Newmark_Beta(TREE, svar, fvar, avar, start, end, a, b, c, d, B, gam, cells, cellsused,
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

	// airP.clear();
	// for(size_t ii = end_ng; ii < end; ++ii)
	// {
	// 	airP.emplace_back(pnp1[ii]);
	// }

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

	Check_If_Ghost_Needs_Removing(svar,fvar,TREE.NP1,pn,pnp1);

	#if DEBUG 	
		dbout << "Exiting first step. Error: " << error1 << endl;
	#endif

	if (svar.Scase == 3)
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

		// real temp = radius / std::cbrt(3.0/(4.0*M_PI));
		bigdrop.GetYcoef(fvar,svar.diam);
		
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
		#else
			real Adrop = 2*(bigdrop.L + bigdrop.Cb*bigdrop.L*ymax);
		#endif

		StateVecD dropDefCd =  0.5*avar.rhog*Vdiff.norm()*Vdiff*Cdl*Adrop;	

		// StateVecD sumForce = StateVecD::Zero();

		// for(size_t ii = 0; ii < end; ++ii)
		// {
		// 	sumForce += pnp1[ii].m*pnp1[ii].acc;	
		// }

		Force = StateVecD::Zero();

		for(size_t ii = 0; ii < end; ii++)
		{
			Force += pnp1[ii].Af*pnp1[ii].m;
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
