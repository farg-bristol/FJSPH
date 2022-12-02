/*********   WCSPH (Weakly Compressible Smoothed SPHPart Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "Integration.h"

#include "Add.h"
#include "Containment.h"
#include "Geometry.h"
#include "Helper_Functions.h"
#include "IPT.h"
#include "Neighbours.h"
#include "Resid.h"
#include "Shifting.h"
#include "shapes/inlet.h"
#include "Runge_Kutta.h"
#include "Newmark_Beta.h"
#include <chrono>
using namespace std::chrono;
///**************** Integration loop **************///
real Integrate(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	VLM& vortex, MESH& cells, SURFS& surf_marks, LIMITS& limits, OUTL& outlist, 
	SPHState& pn, SPHState& pnp1, vector<IPTState>& iptdata)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	size_t const start = svar.bndPts;
	size_t end_ng = svar.bndPts + svar.simPts;
	
	// cout << "Entered Newmark_Beta" << endl;
	#ifndef RK
	uint   k = 0;	//iteration number
	real error2 = 1.0;
	#endif
	real logbase = 0.0;
	real error1 = 0.0;
	real npd = 1.0;	

	// Find maximum safe timestep
	real maxf = 0.0;
	real maxU = 0.0;
	real minST = 9999999.0;
	#pragma omp parallel for reduction(max: maxf, maxU) reduction(min: minST)
	for(size_t ii = start; ii < end_ng; ++ii)
	{
		maxf = std::max(maxf, pnp1[ii].acc.norm());
		maxU = std::max(maxU, pnp1[ii].v.norm());
		minST = std::min(minST, sqrt(pnp1[ii].rho*svar.dx*svar.dx/(2.0*M_PI*fvar.sig*fabs(pnp1[ii].curve))));
	}

	real dtv = 0.125 * fvar.HSQ * fvar.rho0/fvar.mu; /*Viscosity timestep constraint*/
	real dtf = 0.25 * sqrt(fvar.H/maxf); 			 /*Force timestep constraint*/
	real dtc = 2 * fvar.H/(maxU);					 /*Velocity constraint*/
	real dtc2 = 1.5 * fvar.H/fvar.Cs; /* Acoustic constraint */ /* 2* can't be used without delta-SPH it seems. Divergent in tensile instability */

	real dtst = 0.067 * minST;

	// Only use if -fno-finite-math-only is on
	// if (std::isinf(maxf))
	// {
	// 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
	// 	exit(-1);
	// }

	/***********************************************************************************/
	/***********************************************************************************/
	/***********************************************************************************/
		real safe_dt = 0.75*std::min(std::min(dtf,dtc2),std::min(dtc,std::min(dtv,dtst)));
		svar.dt = svar.cfl*safe_dt;
	/***********************************************************************************/
	/***********************************************************************************/
	/***********************************************************************************/

	if(svar.dt < svar.dt_min)
		svar.dt = svar.dt_min;
	else if(svar.dt > svar.dt_max)
		svar.dt = svar.dt_max;

	if (svar.dt > svar.tframem1 + svar.framet - svar.t)
		svar.dt = svar.tframem1 + svar.framet - svar.t + svar.dt_min;

	
	if(svar.speedTest)
	{
		// Bound SPH timestep by the mesh size if it's larger than that to traverse a cell, to prevent skipping.
		if(svar.dt * pnp1[0].v.norm() > cells.minlength)
		{
			svar.dt = cells.minlength/pnp1[0].v.norm();
		}		
	}

	SPH_TREE.index->buildIndex();
	FindNeighbours(SPH_TREE, fvar, pnp1, outlist);

	dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,npd);

	// Check if the particle has moved to a new cell
	if (svar.Asource == meshInfl)
	{
		// cout << "Finding cells" << endl;
		vector<size_t> toDelete = FindCell(svar,avar,CELL_TREE,cells,pn,pnp1);
		    
		if(!toDelete.empty())
		{	
			size_t nDel = toDelete.size();
			std::sort(toDelete.begin(),toDelete.end());
			vector<size_t> nshift(limits.size(),0);
			for(vector<size_t>::reverse_iterator index = toDelete.rbegin(); 
				index!=toDelete.rend(); ++index)
			{
				pnp1.erase(pnp1.begin()+ *index);
				pn.erase(pn.begin() + *index);
				for(size_t block = 0; block < limits.size(); ++block)
				{
					if(*index < limits[block].index.second)
					{
						nshift[block]++;
					}
				}
			}

			for(size_t block = 0; block < limits.size(); ++block)
			{
				for(auto& back:limits[block].back)
					back -= nDel;

				for(auto& buffer:limits[block].buffer)
					for(auto& p:buffer)
						p -= nDel;
			}
			//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += nDel;
			svar.simPts -= nDel;
			svar.totPts -= nDel;
			end_ng -= nDel;
			SPH_TREE.index->buildIndex();
			FindNeighbours(SPH_TREE, fvar, pnp1, outlist);
			dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,npd);
		}	
	}
	#if SIMDIM == 3
	else if (svar.Asource == VLMInfl)
	{
		/* Find the VLM influence for the particles */
		#pragma omp parallel for
		for(size_t ii = start; ii < end_ng ; ii++)
		{
			if(pnp1[ii].lam_ng < avar.cutoff  && pnp1[ii].b == FREE)
			{
				pnp1[ii].cellV = vortex.getVelocity(pnp1[ii].xi);
				pnp1[ii].cellID = 1;
			}
			else
			{
				pnp1[ii].cellV = StateVecD::Zero();
				pnp1[ii].cellID = -1;
			}
		}
	}
	#endif
	else
	{
		#pragma omp parallel for
		for(size_t ii = start; ii < end_ng ; ii++)
		{
			if(pnp1[ii].lam_ng < avar.cutoff  && pnp1[ii].b == FREE)
			{
				pnp1[ii].cellV = avar.vInf;
				pnp1[ii].cellID = 1;
			}
			else
			{
				pnp1[ii].cellV = StateVecD::Zero();
				pnp1[ii].cellID = -1;
			}
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,vortex,pnp1);

	if(svar.ghost == 1)
	{	/* Poisson points */
		PoissonGhost(svar,fvar,avar,cells,SPH_TREE,outlist,pn,pnp1);
		dSPH_PreStep(fvar,end_ng,pnp1,outlist,npd);
	}
	else if (svar.ghost == 2)
	{	/* Lattice points */
		LatticeGhost(svar,fvar,avar,cells,SPH_TREE,outlist,pn,pnp1);
		dSPH_PreStep(fvar,end_ng,pnp1,outlist,npd);
	}

	size_t end = svar.totPts;

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,pnp1);
	#endif

	#ifdef ALE
		if(svar.ghost > 0)
			Particle_Shift_Ghost(svar,fvar,start,end_ng,outlist,pnp1);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end_ng,outlist,pnp1);
	#endif

    Check_Pipe_Outlet(CELL_TREE,svar,avar,cells,limits,pn,pnp1,end,end_ng);

	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	
	const real B = fvar.B;
	const real gam = fvar.gam;

	StateVecD dropVel = StateVecD::Zero();
	StateVecD Force = StateVecD::Zero();

	vector<StateVecD> xih(end-start);
	#pragma omp parallel for default(shared)
	for(size_t ii = start; ii < end; ++ii )
	{
		xih[ii-start] = pnp1[ii].xi;
	}
	
	/*Get preliminary new state to find neighbours and d-SPH values, then freeze*/
	switch(svar.solver_type)
	{ 
		case runge_kutta:
		{
			error1 = Get_First_RK(svar, fvar, avar, start, end, B, gam, npd, cells, 
									limits,outlist, logbase, pn, pnp1, error1);
			break;
		}
		case newmark_beta:
		{
			Do_NB_Iter(CELL_TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,npd,
				cells,limits,outlist,pn,pnp1,Force,dropVel);
		
			uint nUnstab = 0;

			void(Check_Error(SPH_TREE,svar,fvar,start,end,error1,error2,logbase,
									outlist,xih,pn,pnp1,k,nUnstab));
			k++; //Update iteration count
			break;
		}
	}

	// cout << "Error: " << error1 << endl;

	/****** UPDATE TREE ***********/
	SPH_TREE.index->buildIndex();
	FindNeighbours(SPH_TREE, fvar, pnp1, outlist);

	dSPH_PreStep(fvar,end,pnp1,outlist,npd);

	if (svar.Asource == meshInfl)
	{
		// cout << "Finding cells" << endl;
		FindCell(svar,avar,CELL_TREE,cells,pn,pnp1);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.simPts -= svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			end = svar.totPts;
			end_ng = svar.bndPts + svar.simPts;
			SPH_TREE.index->buildIndex();
			FindNeighbours(SPH_TREE, fvar, pnp1, outlist);
			dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,npd);
			
		}	
	}
	#if SIMDIM == 3
	else if (svar.Asource == VLMInfl)
	{
		/* Find the VLM influence for the particles */
		#pragma omp parallel for
		for(size_t ii = start; ii < end_ng ; ii++)
		{
			if(pnp1[ii].lam_ng < avar.cutoff  && pnp1[ii].b == FREE)
			{
				pnp1[ii].cellV = vortex.getVelocity(pnp1[ii].xi);
				pnp1[ii].cellID = 1;
			}
			else
			{
				pnp1[ii].cellV = StateVecD::Zero();
				pnp1[ii].cellID = -1;
			}
		}
	}
	#endif

	Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,vortex,pnp1);

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,pnp1);
	#endif
	
	// Apply_XSPH(fvar,start,end,outlist,pnp1);
	
	#ifdef ALE
		if(svar.ghost > 0)
			Particle_Shift_Ghost(svar,fvar,start,end_ng,outlist,pnp1);
		else
			Particle_Shift_No_Ghost(svar,fvar,start,end_ng,outlist,pnp1);
	#endif
	
    Check_Pipe_Outlet(CELL_TREE,svar,avar,cells,limits,pn,pnp1,end,end_ng);
    
	dropVel = StateVecD::Zero();

	/*Do time integration*/
	switch (svar.solver_type)
	{
		case runge_kutta:
		{
			SPHState st_2 = pnp1;

			error1 = Runge_Kutta4(CELL_TREE,svar,fvar,avar,start,end,B,gam,npd,cells,
						limits,outlist,logbase,pn,st_2,pnp1,Force,dropVel);
			break;
		}
		case newmark_beta:
		{
			Newmark_Beta(SPH_TREE,CELL_TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,npd,cells,
				limits,outlist,logbase,k,error1,error2,xih,pn,pnp1,Force,dropVel);
			break;
		}
	}

	if(svar.ghost != 0)
		Check_If_Ghost_Needs_Removing(svar, fvar, SPH_TREE, pn, pnp1);

	uint nAdd = update_buffer_region(svar, limits, pnp1, end, end_ng);

	/* Check if any particles need to be deleted, and turned to particle tracking */
	vector<size_t> to_del;
	vector<size_t> ndel(limits.size(),0);
	for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
	{
		// Need to check if the delete plane is active.
		if(limits[block].delconst == default_val)
			continue;

		#pragma omp parallel for default(shared)
		for (size_t ii=limits[block].index.first; ii < limits[block].index.second; ++ii)
		{
			// if(pnp1[ii].xi[0] > svar.max_x_sph - fvar.H*fvar.Hfac)
			// {
				/* Create a buffer outlet zone, to avoid truncation of the evolved physics */
				// pnp1[ii].b = OUTLET;
				if(pnp1[ii].xi.dot(limits[block].delete_norm) > limits[block].delconst)
				{
					/* SPHPart is downstream enough to convert to IPT */
					#pragma omp critical
					{
						to_del.emplace_back(ii);
					}
					#pragma omp atomic
					ndel[block]++;
				}
			// }

			// if(pnp1[ii].b == LOST)
			// {
			// 	to_del.emplace_back(ii);
			// }
		}
	}

	if(svar.using_ipt && svar.Asource != constVel)
	{
		vector<IPTPart> IPT_nm1, IPT_n, IPT_np1;
		for(size_t const& ii:to_del)
		{
			IPT_nm1.emplace_back(IPTPart(pnp1[ii],svar.t,svar.IPT_diam,svar.IPT_area));
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

	/* Delete the old particles */
	uint nDel = 0;
	if(!to_del.empty())
	{
		std::sort(to_del.begin(), to_del.end());
		for(vector<size_t>::reverse_iterator itr = to_del.rbegin(); itr!=to_del.rend(); ++itr)
		{
			pnp1.erase(pnp1.begin() + *itr);
			svar.totPts--;
			svar.simPts--;
			end_ng--;
			end--;
			nDel++;
			svar.delNum++;
		}

		size_t delshift = 0;
		for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
		{
			// Shift the blocks, which will happen cumulatively
			limits[block].index.first -= delshift;
			delshift += ndel[block];
			limits[block].index.second -= delshift;

			/* Need to shift back vector and buffer vector for correct inlet */
			for(size_t& back:limits[block].back)
			{
				back -= delshift;
			}

			for(vector<size_t>& buffer:limits[block].buffer)
				for(size_t& part:buffer)
				{
					part -= delshift;
				}
		}
	}

	#ifndef DEBUG
	if(svar.speedTest == 0)
	{
	#endif
		real maxAf = 0.0;
		real maxRhoi = 0.0;
		#ifdef ALE
		real maxShift = 0.0;
		#endif
		#pragma omp parallel for reduction(max: maxAf, maxRhoi)
		for(size_t ii = start; ii < end_ng; ++ii)
		{
			maxAf = std::max(maxAf, pnp1[ii].Af.norm());
			maxRhoi = std::max(maxRhoi, abs(pnp1[ii].rho-fvar.rho0));
			#ifdef ALE
			maxShift = std::max(maxShift,pnp1[ii].vPert.norm());
			#endif
		}
		real maxRho = 100*maxRhoi/fvar.rho0;
		real cfl = svar.dt/safe_dt;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(t2-t1).count();
		#ifdef ALE
		printf("%9.3e | %7.2e | %4.2f | %8.4f | %7.3f | %9.3e | %9.3e | %9.3e | %14ld|\n",
				svar.t, svar.dt, cfl, error1, maxRho, maxf, maxAf, maxShift, duration);

		#else
		printf("%9.3e | %7.2e | %4.2f | %9.4f | %3u | %8.3f | %9.3e | %9.3e | %14ld|\n",
				svar.t, svar.dt, cfl, error1, k, maxRho, maxf, maxAf, duration);
		#endif
		// cout << "time: " << svar.t << " dt: " << svar.dt <<  " CFL: " << cfl 
		// 	<< " dRho: " << maxRho << " Maxf: " << maxf << " MaxAf: " << maxAf 
		// 	#ifdef ALE
		// 	<< " MaxdU: " << maxShift
		// 	#endif
		// 	 << endl;
	#ifndef DEBUG
	}
	#endif

	if(pnp1.size() == 0 || svar.simPts == 0)
	{ 
		// cout << "No more particles. Stopping." << endl;
		return 0;
	}

	if(nAdd != 0 || nDel != 0)
	{
		SPH_TREE.index->buildIndex();
		FindNeighbours(SPH_TREE, fvar, pnp1, outlist);
	}


	/****** UPDATE TIME N ***********/
	if(svar.totPts != pnp1.size())
	{
		cout << "Size mismatch. Total points not equal to array size. Stopping" << endl;
		exit(-1);
	}

	// cout << "Updating. Time: " << svar.t << "  dt: " << svar.dt << endl;
	if(pn.size() != pnp1.size())
		pn.resize(pnp1.size());

	copy_omp(pnp1,pn);

	/*Add time to global*/
	svar.t+=svar.dt;

	#ifdef DAMBREAK
	/* Find max x and y coordinates of the fluid */
		// Find maximum safe timestep
	vector<SPHPart>::iterator maxXi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](SPHPart const& p1, SPHPart const& p2){return p1.xi[0] < p2.xi[0];});

	vector<SPHPart>::iterator maxYi = std::max_element(pnp1.begin()+svar.bndPts,pnp1.end(),
		[](SPHPart const& p1, SPHPart const& p2){return p1.xi[1] < p2.xi[1];});

	real maxX = maxXi->xi[0];
	real maxY = maxYi->xi[1];

	dambreak << svar.t << "  " << maxX+0.5*svar.Pstep << "  " << maxY+0.5*svar.Pstep << endl;
	#endif

	return error1;
}

void First_Step(Sim_Tree& SPH_TREE, Vec_Tree const& CELL_TREE, SIM& svar, FLUID const& fvar, AERO const& avar, 
	VLM& vortex, MESH& cells, LIMITS const& limits, OUTL& outlist, SPHState& pnp1, SPHState& pn, vector<IPTState>& iptdata)
{
	size_t const start = svar.bndPts;
	size_t end = svar.totPts;
	size_t end_ng = svar.bndPts + svar.simPts;
	real npd = 1.0;
	#if DEBUG 
		fprintf(dbout, "Starting first step. Start index: %zu End index: %zu\n",start, end);
	#endif

	dSPH_PreStep(fvar,end,pnp1,outlist,npd);
	if (svar.Asource == meshInfl)
	{
		FindCell(svar,avar,CELL_TREE,cells,pnp1,pn);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			svar.simPts = svar.totPts-svar.bndPts;
			end = svar.totPts;
			end_ng = end;
			SPH_TREE.index->buildIndex();
			FindNeighbours(SPH_TREE, fvar, pnp1, outlist);
			dSPH_PreStep(fvar,end,pnp1,outlist,npd);
			
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,vortex,pnp1);

	if(svar.ghost == 1)
	{	/* Poisson points */
		PoissonGhost(svar,fvar,avar,cells,SPH_TREE,outlist,pn,pnp1);
		dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,npd);
		// Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,pnp1);
	}
	else if (svar.ghost == 2)
	{	/* Lattice points */
		LatticeGhost(svar,fvar,avar,cells,SPH_TREE,outlist,pn,pnp1);
		dSPH_PreStep(fvar,svar.totPts,pnp1,outlist,npd);
		// Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,pnp1);
	}

	end = svar.totPts;

	uint k = 0;	
	real error2 = 0.0;
	real logbase = 0.0;
	real error1 = 0.0;

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

	// Find maximum safe timestep
	real maxf = 0.0;
	real maxU = 0.0;
	real minST = 9999999.0;
	#pragma omp parallel for reduction(max: maxf, maxU) reduction(min: minST)
	for(size_t ii = start; ii < end_ng; ++ii)
	{
		maxf = std::max(maxf, pnp1[ii].acc.norm());
		maxU = std::max(maxU, pnp1[ii].v.norm());
		minST = std::min(minST, sqrt(pnp1[ii].rho*svar.dx*svar.dx/(2.0*M_PI*fvar.sig*fabs(pnp1[ii].curve))));
	}

	real dtv = 0.125 * fvar.HSQ * fvar.rho0/fvar.mu; /*Viscosity timestep constraint*/
	real dtf = 0.25 * sqrt(fvar.H/maxf); 			 /*Force timestep constraint*/
	real dtc = 2 * fvar.H/(maxU);					 /*Velocity constraint*/
	real dtc2 = 1.5 * fvar.H/fvar.Cs; /* Acoustic constraint */ /* 2* can't be used without delta-SPH it seems. Divergent in tensile instability */

	real dtst = 0.067 * minST;

	// Only use if -fno-finite-math-only is on
	// if (std::isinf(maxf))
	// {
	// 	std::cerr << "Forces are quasi-infinite. Stopping..." << std::endl;
	// 	exit(-1);
	// }

	/***********************************************************************************/
	/***********************************************************************************/
	/***********************************************************************************/
		real safe_dt = 0.75*std::min(std::min(dtf,dtc2),std::min(dtc,std::min(dtv,dtst)));
		svar.dt = svar.cfl*safe_dt;
	/***********************************************************************************/
	/***********************************************************************************/
	/***********************************************************************************/

	if(svar.dt < svar.dt_min)
		svar.dt = svar.dt_min;
	else if(svar.dt > svar.dt_max)
		svar.dt = svar.dt_max;

	if (svar.dt > svar.tframem1 + svar.framet - svar.t)
		svar.dt = svar.tframem1 + svar.framet - svar.t + svar.dt_min;


	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,pnp1);
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
	switch(svar.solver_type)
	{ 
		case runge_kutta:
		{
			error1 = Get_First_RK(svar, fvar, avar, start, end, B, gam, npd, cells, 
									limits,outlist, logbase, pn, pnp1, error1);
			break;
		}
		case newmark_beta:
		{
			Do_NB_Iter(CELL_TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,npd,
				cells,limits,outlist,pn,pnp1,Force,dropVel);
		
			uint nUnstab = 0;

			void(Check_Error(SPH_TREE,svar,fvar,start,end,error1,error2,logbase,
									outlist,xih,pn,pnp1,k,nUnstab));
			k++; //Update iteration count
			break;
		}
	}

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
	SPH_TREE.index->buildIndex();
	FindNeighbours(SPH_TREE, fvar, pnp1, outlist);

	dSPH_PreStep(fvar,end,pnp1,outlist,npd);

	if (svar.Asource == meshInfl)
	{
		FindCell(svar,avar,CELL_TREE,cells,pnp1,pn);
		if (svar.totPts != pnp1.size())
		{	//Rebuild the neighbour list
			// cout << "Updating neighbour list" << endl;
			// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			svar.delNum += svar.totPts-pnp1.size();
			svar.totPts = pnp1.size();
			svar.simPts = svar.totPts-svar.bndPts;
			end = svar.totPts;
			end_ng = end;
			SPH_TREE.index->buildIndex();
			FindNeighbours(SPH_TREE, fvar, pnp1, outlist);
			dSPH_PreStep(fvar,end,pnp1,outlist,npd);
			
		}
	}

	Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,vortex,pnp1);

	#ifndef NOFROZEN
	dissipation_terms(fvar,start,end,outlist,pnp1);
	#endif
	
	/*Do time integration*/
	switch (svar.solver_type)
	{
		case runge_kutta:
		{
			SPHState st_2 = pnp1;

			error1 = Runge_Kutta4(CELL_TREE,svar,fvar,avar,start,end,B,gam,npd,cells,
						limits,outlist,logbase,pn,st_2,pnp1,Force,dropVel);
			break;
		}
		case newmark_beta:
		{
			Newmark_Beta(SPH_TREE,CELL_TREE,svar,fvar,avar,start,end,a,b,c,d,B,gam,npd,cells,
				limits,outlist,logbase,k,error1,error2,xih,pn,pnp1,Force,dropVel);
			break;
		}
	}

	Check_If_Ghost_Needs_Removing(svar,fvar,SPH_TREE,pn,pnp1);

	#if DEBUG 	
		fprintf(dbout,"Exiting first step. Error: %f\n", error1);
	#endif
}
