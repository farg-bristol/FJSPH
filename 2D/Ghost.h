#ifndef GHOST_H
#define GHOOST_H

#include "Var.h"
#include "Neighbours.h"

void Ghost_Particles(Sim_Tree &NP1_INDEX, ldouble numpartdens,
	SIM &svar, FLUID &fvar, CROSS &cvar, State &pn, State &pnp1, outl &outlist)
{
/***** Find particles outside of the influence of the Liquid *******/

	/*Default to an air particle that is outside the influence of the simulation*/
	for (size_t i=svar.totPts; i< svar.totPts+svar.aircount; ++i)
		pnp1[i].b = 4;

	for (size_t i=svar.bndPts; i< svar.totPts; ++i)
	{
		for (size_t j=0; j < outlist[i].size(); ++j)
		{	/*If it's inside the list of neighbours of the liquid, keep it.*/
			if (pnp1[outlist[i][j]].b == 4) pnp1[outlist[i][j]].b = 3;
		}
	}

	for (size_t i=svar.totPts; i< svar.totPts+svar.aircount; ++i)
	{	/*If it's still outside the influence of the liquid, then delete it.*/
		if(pnp1[i].b == 4)
		{
			pnp1[i] = pnp1.back();
			pnp1.pop_back();
		}
	}

/****** Find particles with reduced density to create particles around ********/


/********** Create more air particles *********/
	std::vector<StateVecD> temp; /*Temporary storage for air particles*/
	for (size_t i=svar.bndPts; i< svar.totPts; ++i)
	{	/*Find the surface of fluid particles.*/
		Particle pi = pnp1[i];
		StateVecD SurfC= StateVecD::Zero();
		pi.b = 1;

		/*Surface tension factor*/
		const static ldouble lam = (6.0/81.0*pow((2.0*fvar.H),4.0)/pow(M_PI,4.0)*
							(9.0/4.0*pow(M_PI,3.0)-6.0*M_PI-4.0));
		for (size_t j=0; j < outlist[i].size(); ++j)
		{
			Particle pj = pnp1[outlist[i][j]];
			/*Check if the position is the same, and skip the particle if yes*/
			if(pi.xi == pj.xi)
				continue;

			StateVecD Rij = pj.xi-pi.xi;
			ldouble r = Rij.norm();

			/*Surface Tension as described by Nair & Poeschel (2017)*/
			ldouble fac=1.0;
			if(pj.b==0)
	            fac=(1+0.5*cos(M_PI*(fvar.contangb/180)));
	        //cout << lam(svar.H) << endl;
			ldouble sij = 0.5*pow(numpartdens,-2.0)*(fvar.sig/lam)*fac;
			SurfC -= (Rij/r)*sij*cos((3.0*M_PI*r)/(4.0*fvar.H));
		}

		if (SurfC.norm()>0.05)
		{
			/*Create particles... */
			namespace pds = thinks::poisson_disk_sampling;
			ldouble radius = svar.Pstep;
			std::array<ldouble,2> xmin = {pn[i].xi[0]-2.0*fvar.H, pn[i].xi[1]-2.0*fvar.H};
			std::array<ldouble,2> xmax = {pn[i].xi[0]+2.0*fvar.H,pn[i].xi[1]+2.0*fvar.H};

			std::vector<StateVecD> samples =
			pds::PoissonDiskSampling<ldouble,2,StateVecD,EVecTraits>(radius,xmin,xmax);

			// cout << samples.size() << endl;
			for (auto j:samples)
			{
				StateVecD xi(j[0],j[1]);
				temp.emplace_back(xi);
			}
		}

	}

	if(temp.size()!=0)
	{
		Temp_Tree temp_index(2,temp,10);
		temp_index.index->buildIndex();
		ldouble search_radius = svar.Pstep*svar.Pstep;
		//std::vector<size_t> delete_list;
		//cout << temp.size() << endl;
    nanoflann::SearchParams params;

		for (auto i=temp.begin(); i!=temp.end(); )
		{	/*Check for duplicate particles and delete them when too close.*/
			StateVecD xi = *i;
			std::vector<std::pair<size_t,ldouble>> matches;
			temp_index.index->radiusSearch(&xi[0],search_radius,matches,params);
			//cout << matches.size() << endl;
			if (matches.size()!=1)
			{
				for (auto j:matches)
				{
					if (j.second == 0.0)
						continue;

					temp[j.first] = temp.back();
					temp.pop_back();
				}
				temp_index.index->buildIndex();
			}
			else
				++i;
		}

		for (size_t i=svar.bndPts; i < svar.totPts; ++i)
		{	/*Check for particles inside the fluid*/
			std::vector<std::pair<size_t,ldouble>> matches;
			temp_index.index->radiusSearch(&pn[i].xi[0],search_radius,matches,params);
			//cout << matches.size()<< endl;
			if (matches.size()!=0)
			{
				for (auto j:matches)
				{
					temp[j.first] = temp.back();
					temp.pop_back();
				}
				temp_index.index->buildIndex();
			}
		}

		// search_radius = 2*svar.Pstep*svar.Pstep;
		// for (size_t i=0; i< svar.bndPts; ++i)
		// {	/*Check for particles next to the boudnary*/
		// 	std::vector<std::pair<size_t,double>> matches;
		// 	temp_index.index->radiusSearch(&pn[i].xi[0],search_radius,matches,params);
		// 	//cout << matches.size()<< endl;
		// 	if (matches.size()!=0)
		// 	{
		// 		for (auto j:matches)
		// 		{
		// 			temp[j.first] = temp.back();
		// 			temp.pop_back();
		// 		}
		// 		temp_index.index->buildIndex();
		// 	}
		// }

		/*Place particles in the simulation vector*/
		svar.aircount = temp.size();

		ldouble rho = 1.225;
		StateVecD f = StateVecD::Zero();
		ldouble airmass = fvar.Simmass*rho/fvar.rho0;
		for (auto i:temp)
		{
			pnp1.emplace_back(Particle(i,cvar.vInf,f,rho,airmass,1));
			pn.emplace_back(Particle(i,cvar.vInf,f,rho,airmass,1));
		}

		NP1_INDEX.index->buildIndex();
		FindNeighbours(NP1_INDEX,fvar,pnp1,outlist);
	}
}

#endif