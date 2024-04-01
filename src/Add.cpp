/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "Add.h"
#include "Containment.h"
#include "Neighbours.h"
#include <random>
#include <stdint.h>
#include <time.h>

namespace PoissonSample
{
	class PRNG
	{
	public:
		PRNG(): gen_(std::random_device()()), dis_( 0.0, 1.0 )
		{
			// prepare PRNG
		}

		real randomReal()
		{
			return static_cast<real>( dis_( gen_ ) );
		}

		int randomInt( int maxValue )
		{
			std::uniform_int_distribution<> disInt( 0, maxValue );
			return disInt( gen_ );
		}

	private:
		std::mt19937 gen_;
		std::uniform_real_distribution<real> dis_;
	};

	bool isInCircle(StateVecD const& centre, StateVecD const& p, real const radius)
	{
		const StateVecD f = p - centre;
		return f.squaredNorm() <= radius;
	}

	StateVecI imageToGrid(StateVecD const& P, real const cellSize )
	{
		#if SIMDIM == 2
		return StateVecI( (int)(P(0)/cellSize), (int)(P(1)/cellSize));
		#else
		return StateVecI((int)(P(0)/cellSize), (int)(P(1)/cellSize), (int)(P(2)/cellSize));
		#endif
	}

	struct Grid
	{
		Grid( uint w, real minDist, real cellSize ):
			 minDist_(minDist), cellSize_( cellSize ), w_(w), h_(w)
		#if SIMDIM == 3
		, d_(w)
		#endif
		{
			#if SIMDIM == 2
				grid_ = std::vector<std::vector<StateVecD>>(w,std::vector<StateVecD>(w,StateVecD::Zero()));
			#endif
			#if SIMDIM == 3
				grid_ = std::vector<std::vector<std::vector<StateVecD>>>
				(w,std::vector<std::vector<StateVecD>>(w,std::vector<StateVecD>(w,StateVecD::Zero())));
			#endif
		}

		void insert( StateVecD const& p)
		{
			const StateVecI g = imageToGrid(p, cellSize_);
			if(grid_.size() == 0)
			{
				cout << "Ghost particle grid size is zero." << endl;
				exit(-1);
			}

			if(g(0) >  static_cast<int>(grid_.size()))
			{
				std::cout << "Tried to access grid_ out of bounds in i direction." << std::endl;
				std::cout << g(0) << "  " << g(1) << std::endl;
				exit(-1);
			}
			else if ( g(1) > static_cast<int>(grid_[g(0)].size()))
			{
				std::cout << "Tried to access grid_ out of bounds in j direction." << std::endl;
				exit(-1);
			}
			// }

			#if SIMDIM == 2
				grid_[g(0)][g(1)] = p;
			#endif
			#if SIMDIM == 3
				grid_[g(0)][g(1)][g(2)] = p;
			#endif
		}

		bool isInNeighbourhood( StateVecD const& point)
		{
			StateVecI g = imageToGrid(point, cellSize_);

			// number of adjucent cells to look for neighbour points
			const int D = 5;

			// scan the neighbourhood of the point in the grid
			for ( int ii = g(0) - D; ii < g(0) + D; ii++ )
			{
				for ( int jj = g.y() - D; jj < g.y() + D; jj++ )
				{	
					#if SIMDIM == 2
						if ( ii >= 0 && ii < int(w_) && jj >= 0 && jj < int(h_) )
						{
							const StateVecD P = grid_[ii][jj];

							if ( (P-point).norm() < minDist_ ) { return true; }
						}
					#endif

					#if SIMDIM == 3
						for (int kk = g.z() - D; kk < g.z() + D; kk++)
						{
							if ( ii >= 0 && ii < int(w_) && jj >= 0 && jj < int(h_)
							&& kk >= 0 && kk < int(d_) )
							{
								const StateVecD P = grid_[ii][jj][kk];

								if ( (P-point).norm() < minDist_ ) { return true; }
							}
						}
					#endif
				}
			}

			return false;
		}

	private:

		real minDist_, cellSize_;
		uint w_;
		uint h_;
		
		
		#if SIMDIM == 2
			std::vector<std::vector<StateVecD>> grid_;
		#endif

		#if SIMDIM == 3
			uint d_;
			std::vector<std::vector<std::vector<StateVecD>>> grid_;
		#endif
	};

	StateVecD popRandom(std::vector<StateVecD>& points, PRNG& generator)
	{
		const int idx = generator.randomInt( points.size()-1 );
		const StateVecD p = points[idx];
		points.erase( points.begin() + idx );
		return p;
	}

	StateVecD generateRandomPointAround( StateVecD const& p, real minDist, PRNG& generator )
	{
		// start with non-uniform distribution
		const real R1 = generator.randomReal();
		const real R2 = generator.randomReal();
		
		// radius should be between MinDist and 2 * MinDist
		const real radius = minDist * ( R1 + 1.0 );

		// random angle
		const real angle1 = 2 * M_PI * R2;

		#if SIMDIM == 3
			const real R3 = generator.randomReal();
			const real angle2 = 2 * M_PI * R3;
		#endif

		// the new point is generated around the point (x, y)
		#if SIMDIM == 2
		return StateVecD(p(0) + radius*cos(angle1), p(1) + radius*sin(angle1));
		#endif
		#if SIMDIM == 3
		return StateVecD(p(0) + radius*cos(angle1)*sin(angle2), 
						 p(1) + radius*sin(angle1)*sin(angle2),
						 p(2) + radius*cos(angle2));
		#endif
	}

	/**
		Return a vector of generated points
		sampleLimit - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
	**/
	SPHState generatePoissonPoints(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells,
	 uint const& host, SPHState const& pnp1, OUTL const& outlist/* , StateVecD const& norm, StateVecD const& avgV */)
	{
		/*Variables for the poisson disk sampling*/
		real radius = 1.2*fvar.sr;
		uint sampleLimit = 30;
		PRNG generator;
		#if SIMDIM == 3
			uint numPoints = 257;
		#else
			uint numPoints = 48;
		#endif
		uint pID = svar.totPts;

		/*Properties for new particles*/

		/* Based off the aerodynamic mesh */
		// StateVecD vel= avar.vInf;
		// StateVecD Vdiff = StateVecD::Zero();
		// real Pbase = 0.0;
		// real press = 0;
		// real rho = fvar.rho0;

// 		if(svar.Asource == 1)
// 		{
// 			vel = pnp1[host].cellV;
// 			Vdiff =  vel - avgV;
// 			Pbase = pnp1[host].cellP - avar.pRef;
// 		}
// 		else if (svar.Asource == 2)
// 		{
// 			vel = (pnp1[host].cellV+cells.cPertnp1[pnp1[host].cellID]);
// 			Vdiff =  vel - /*pi.v*/ avgV;
// 			Pbase = pnp1[host].cellP - avar.pRef;
// 		}
// #if SIMDIM == 3
// 		else if(svar.Asource == 3)
// 		{	
// 			vel = svar.vortex.getVelocity(pnp1[host].xi);
// 			Vdiff = vel - avgV;
// 			Pbase = 0.5*avar.rhog*(pow(avar.vRef,2.0)-pow(vel.norm(),2.0));
// 		}
// #endif
// 		else
// 		{
// 			Vdiff = vel - avgV;
// 			Pbase = 0.5*avar.rhog*Vdiff.squaredNorm();
// 		}

// 		real theta = acos(-norm.normalized().dot(Vdiff.normalized()));
		
// 		real Cp = 0.0;

// 		if(abs(theta) < 2.4877)
// 		{
// 			Cp = 1.0 - 2.5*pow(sin(abs(theta)),2.0);
// 		}
// 		else
// 		{
// 			Cp = 0.075;
// 		}


		// press = Pbase + 0.5*avar.rhog*Vdiff.squaredNorm()*Cp;
		// rho = fvar.rho0 * pow((press/fvar.B + 1.0),1.0/fvar.gam);


		// const real rho = pnp1[host].cellRho;
		// const real mass = fvar.rhog* pow(svar.Pstep, SIMDIM);

				// StateVecD vel= avar.vInf;
		// StateVecD Vdiff = StateVecD::Zero();
		// real Pbase = 0.0;
		// real press = 0;
		// real rho = fvar.rho0;

		/* Host particle properties */
		StateVecD vel= pnp1[host].v;
		real press = pnp1[host].p;
		real rho = pnp1[host].rho;

		real const& mass = pnp1[host].m;

		real const deltax = pow(mass/rho, 1.0/real(SIMDIM));
		
		// #pragma omp critical
		// cout << "CellID: " << pnp1[host].cellID << "  " << press 
		// << "  " << rho << "  " << mass << "  " << deltax << endl;

 		const real minDist = /*svar.Pstep*/ deltax;
		SPHState samplePoints;
		std::vector<StateVecD> processList;
		SPHState airP;

		// create the grid
		#if SIMDIM == 2
			const StateVecD origin(pnp1[host].xi(0)-2*fvar.H,pnp1[host].xi(1)-2*fvar.H);
		#endif
		#if SIMDIM == 3
			const StateVecD origin(pnp1[host].xi(0)-2*fvar.H,pnp1[host].xi(1)-2*fvar.H,
									pnp1[host].xi(2)-2*fvar.H);
		#endif

		const real cellSize = minDist / sqrt(2.0);
		const uint gridW = (uint)ceil(4*fvar.H / cellSize);
		// cout << "GridW: " << gridW << " minDist: " << minDist << " cellSize: " << cellSize << endl; 
		Grid grid(gridW, minDist, cellSize);

		/*Fill out the prexisting particles to add points around*/
		for(auto jj:outlist[host])
		{
			samplePoints.push_back(pnp1[jj.first]);
			grid.insert(pnp1[jj.first].xi-origin);
		}

		/*Try and add a point where it won't conflict*/
		#if SIMDIM == 2
			const StateVecD circCent(2*fvar.H, 2*fvar.H);
		#endif
		#if SIMDIM == 3
			const StateVecD circCent(2*fvar.H, 2*fvar.H, 2*fvar.H);
		#endif
		
		StateVecD firstPoint;
	 	
		do {/*Generate a random number between 0 and 1, then normalise it to the grid.*/
			#if SIMDIM == 2
			firstPoint = StateVecD(generator.randomReal()*4*fvar.H, generator.randomReal()*4*fvar.H);
			#endif
			#if SIMDIM == 3
			firstPoint = StateVecD(generator.randomReal()*4*fvar.H, generator.randomReal()*4*fvar.H,
						generator.randomReal()*4*fvar.H);
			#endif	

		} while (isInCircle(circCent,firstPoint,radius) != 1 || grid.isInNeighbourhood(firstPoint) != 0);
		
		// update containers
		processList.push_back(firstPoint);
		
		grid.insert(firstPoint);

		// generate new points for each point in the queue
		while (!processList.empty() && samplePoints.size() < floor(float(0.95*numPoints)))
		{
			const StateVecD point = popRandom(processList, generator);

			for (uint ii = 0; ii < sampleLimit; ii++)
			{
				const StateVecD newPoint = generateRandomPointAround(point, minDist, generator);
				const bool canFitPoint = isInCircle(circCent,newPoint,radius);

				if (canFitPoint == true && grid.isInNeighbourhood(newPoint) == false)
				{
					processList.push_back(newPoint);
					// StateVecD Point = newPoint+origin;
					samplePoints.push_back(SPHPart(newPoint+origin,StateVecD::Zero(),0.0,0.0,0.0,GHOST,0));
					airP.emplace_back(SPHPart(newPoint+origin,vel,press,rho,mass,GHOST,pID));
					pID++;
					grid.insert(newPoint);
					// cout << "New Point: " << point(0) << " " << point(1) << endl;
					continue;
				}
			}
		}
		return airP;
	}
}

void PoissonGhost(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, Sim_Tree& NP1_INDEX, OUTL& outlist, SPHState& pn, SPHState& pnp1)
{
	size_t const& start = svar.bndPts;
	size_t const& end = svar.totPts;

	size_t nGhost = 0;

	for(size_t ii = start; ii < end; ii++)
	{
		if(pnp1[ii].surf == 1)
		{
			/* Create a lattice grid (perturbed, so its not a perfect grid) around the point. */
			SPHState ghost_particles = PoissonSample::generatePoissonPoints(svar,fvar,avar,cells,ii,pnp1,outlist);		

			if(!ghost_particles.empty())
			{
				/* Add ghost particles to the vector */
				pnp1.insert(pnp1.end(),ghost_particles.begin(),ghost_particles.end());
				pn.insert(pn.end(),ghost_particles.begin(),ghost_particles.end());
				/* Rebuild the tree, including the ghost particles just made, so no overlap */
				NP1_INDEX.index->buildIndex();
				outlist = find_neighbours(NP1_INDEX, fvar, pnp1);
				nGhost += ghost_particles.size();
			}
		}
	}

	svar.gstPts = nGhost;
	svar.totPts = svar.bndPts + svar.simPts + svar.gstPts;  

	/*Write boundary particles*/
	// Write_Binary_Timestep(svar,pnp1,svar.bndPts+svar.simPts,svar.bndPts+svar.simPts+svar.gstPts,"Ghost",3,svar.ghostFile);	
}

void check_if_too_close(Sim_Tree const& NP1_INDEX, real const& sr, SPHState const& pnp1, StateVecD const& xi, 
		StateVecD const& vel, real const& dens, real const& mass, real const& press, size_t const& pID, SPHState& ghost_particles, size_t& nGhost)
{
	/* Check if point is too close to an existing point... */
	vector<size_t> ret_indexes(1);
	vector<real> out_dists_sqr(1);

	nanoflann::KNNResultSet<real> resultSet(1);
	resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
	
	NP1_INDEX.index->findNeighbors(resultSet, &xi[0], flann_params);

	if(out_dists_sqr[0] > sr)
	{
		/* it's far enough away, so point can be added to the array */
		ghost_particles.emplace_back(SPHPart(xi, vel, dens, mass, press, GHOST, pID));
		nGhost++;
	}
}

void LatticeGhost(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, 
			Sim_Tree& SPH_TREE, OUTL& outlist, SPHState& pn, SPHState& pnp1, LIMITS const& limits)
{
	// size_t const& start = svar.bndPts;
	// size_t const& end = svar.bndPts+svar.simPts;

	size_t pID = svar.totPts;
	size_t nGhost = 0;

	real const sr = 0.99*svar.dx*svar.dx;
	int const interval = 1000;

	for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; ++block)
	{
		for(size_t ii = limits[block].index.first; ii < limits[block].index.second; ii++)
		{
			if(pnp1[ii].surf == 1)
			{
				/* Check if it's neighbourhood is fully supported?  */
				/* Find how many are within 2H+dx */
				real const search_radius = (2*fvar.H+svar.dx)*(2*fvar.H+svar.dx);

				std::vector<neighbour_index> matches = radius_search(SPH_TREE, pnp1[ii].xi, search_radius);

				// cout << matches.size() << endl;

				#if SIMDIM == 3
				if(matches.size() >  250/* ? not sure yet */)
					continue;
				#else
				if(matches.size() >  90/* ? not sure yet */)
					continue;
				#endif
				/* Create a lattice grid (perturbed, so its not a perfect grid) around the point. */
				SPHState ghost_particles;

				/* Velocity of the air particle */
				StateVecD vel; 
				if(svar.Asource == constVel)
					vel = (avar.gasM/ pnp1[ii].m) * avar.vInf;
				else
					vel = (avar.gasM/ pnp1[ii].m) * pnp1[ii].cellV;

				real press = pnp1[ii].cellP;
				real rho = fvar.rho0*pow((press/fvar.B) + 1.0, 1.0/fvar.gam);
				real mass = fvar.simM;
				/* How far away does it need to be? */
				/* The fluid surface cannot interact with a 'surface' particle */
				/* Origin */
				StateVecD const& origin = pnp1[ii].xi;

				for(real x = 0.0; x < (2*fvar.H+svar.dx); x += svar.dx)
				{
					for(real y = 0.0; y < (2*fvar.H+svar.dx); y += svar.dx)
					{
						#if SIMDIM == 3
						for(real z = 0.0; z < (2*fvar.H+svar.dx); z += svar.dx)
						{
							StateVecD perturb(random(interval), random(interval), random(interval));
							StateVecD xi = StateVecD(x,y,z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);	

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(-x,y,z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(x,-y,z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(x,y,-z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(-x,-y,z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(-x,y,-z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(x,-y,-z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval), random(interval));
							xi = StateVecD(-x,-y,-z)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);	
						}
						#else
							StateVecD perturb(random(interval), random(interval));
							StateVecD xi = StateVecD(x,y)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);	

							perturb = StateVecD(random(interval),random(interval));
							xi = StateVecD(-x,y)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval));
							xi = StateVecD(x,-y)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);

							perturb = StateVecD(random(interval),random(interval));
							xi = StateVecD(-x,-y)+origin+perturb;
							check_if_too_close(SPH_TREE,sr,pnp1,xi,vel,rho,mass,press,pID,ghost_particles,nGhost);
						#endif
						
					}
				}

				if(!ghost_particles.empty())
				{
					/* Add ghost particles to the vector */
					pnp1.insert(pnp1.end(),ghost_particles.begin(),ghost_particles.end());
					pn.insert(pn.end(),ghost_particles.begin(),ghost_particles.end());

					/* Rebuild the tree, including the ghost particles just made, so no overlap */
					SPH_TREE.index->buildIndex();
				}
			}
		}
	}

	svar.gstPts += nGhost;
	svar.totPts = svar.bndPts + svar.simPts + svar.gstPts; 
	
	outlist = find_neighbours(SPH_TREE, fvar, pnp1);

	/* Find the aerodynamic cell the ghost particle is in, and apply its settings*/
	// for(size_t ii = end; ii < end+svar.gstPts; ++ii)
	// {
	// 	uint to_del = 0;
	// 	FirstCell(svar,SPH_TREE.CELL, cells, pnp1[ii], to_del);

	// 	if(pnp1[ii].cellID != 0)
	// 	{
	// 		/* Make the particle properties those of the cell. */
	// 		pnp1[ii].v = pnp1[ii].cellV;
	// 		pnp1[ii].p = pnp1[ii].cellP;
	// 		pnp1[ii].rho = fvar.rho0 * pow((pnp1[ii].p / fvar.B + 1), 1 / fvar.gam);
	// 	}
	// 	else
	// 	{
	// 		/* Otherwise set to a resting settings */
	// 		pnp1[ii].p = 0.0;
	// 		pnp1[ii].rho = fvar.rho0;
	// 	}

	// }

	/*Write boundary particles*/
	// Write_Binary_Timestep(svar,pnp1,svar.bndPts+svar.simPts,svar.bndPts+svar.simPts+svar.gstPts,"Ghost",3,svar.ghostFile);
}

/* Function to check if a ghost particle needs to be removed from the simulation, because */
/* it's left the support of the fluid */
void Check_If_Ghost_Needs_Removing(SIM& svar, FLUID const& fvar, Sim_Tree& NP1_INDEX, SPHState& pn, SPHState& pnp1)
{
	size_t const& start = svar.bndPts + svar.simPts;
	size_t const& end = svar.totPts;

	vector<size_t> to_del;
	real const search_radius = (2*fvar.H+svar.dx)*(2*fvar.H+svar.dx);

	#pragma omp parallel for default(shared)
	for(size_t ii = start; ii < end; ++ii)
	{
		uint has_interaction = 0;
		std::vector<neighbour_index> matches = radius_search(NP1_INDEX, pnp1[ii].xi, search_radius);

		for(neighbour_index const& jj: matches)
		{
			if(pnp1[jj.first].b == FREE || pnp1[jj.first].b == PIPE)
			{
				has_interaction = 1;
				break;
			}
		}	

		if(has_interaction == 0)
		{
			#pragma omp critical
			to_del.emplace_back(ii);
		}
	}

	/* remove the ghost particles that no longer interact */

	if(!to_del.empty())
	{
		/* Make sure it's sorted first, otherwise it will mess things up */
		std::sort(to_del.begin(), to_del.end());
		for(vector<size_t>::reverse_iterator itr = to_del.rbegin(); itr!=to_del.rend(); ++itr)
		{
			pnp1.erase(pnp1.begin() + *itr);
			svar.totPts--;
			svar.gstPts--;
		}
	}
}
