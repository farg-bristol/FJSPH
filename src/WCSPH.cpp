/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/***			        + XSPH Correction (Also described in Monaghan)                    ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)   (Currently inactive)           ***/
/*** Smoothing Kernel: Wendland's C2 ***/
/*** Integrator: Newmark-Beta ****/
/*** Variable Timestep Criteria: CFL + Monaghan, J.J. (1989) conditions ***/

#include <chrono>

#include "Var.h"
#include "IO.h"
#include "Neighbours.h"
#include "Kernel.h"
#include "Init.h"
#include "Add.h"
#include "Resid.h"
#include "Crossing.h"
#include "Integration.h"

using namespace std::chrono;
using namespace nanoflann;

int main(int argc, char *argv[])
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
	Eigen::initParallel();
	// omp_set_num_threads(1);

    real duration;
    real error = 0;
    cout.width(13);
    cout << std::scientific << std::left << std::setprecision(4);
	
	write_header();
    
    /******* Define the global simulation parameters ******/
	SIM svar;
	FLUID fvar;
	AERO avar;
	outl outlist;
	MESH cells;

	GetInput(argc,argv,svar,fvar,avar);
	
	if(MakeOutputDir(argc,argv,svar))
	{
		cout << "Couldn't make output directory. Please check permissions." << endl;
		exit(-1);
	}

	if(svar.Asource == 1 || svar.Asource == 2)
	{
		#if SIMDIM == 3
		Read_TAUMESH_FACE(svar,cells,fvar,avar);
		// Read_TAUMESH(svar,cells,fvar);
		#else
		Read_TAUMESH_EDGE(svar,cells,fvar,avar);
		#endif
	}
	else if (svar.Asource == 4)
	{
		// Create single cell
		Make_Cell(fvar,avar,cells);
	}	

	cout << std::setprecision(5);

	
	cout << "Adjusted Start Coordinates: " << endl;
	cout << svar.Start(0) << "  " << svar.Start(1);
#if SIMDIM == 3
	cout << "  " << svar.Start(2); 
#endif
	cout << endl << endl;
	/*Make a guess of how many there will be...*/
	// int partCount = ParticleCount(svar);
    ///****** Initialise the particles memory *********/
	State pn;	    /*Particles at n   */
	State pnp1; 	/*Particles at n+1 */
	State airP;
	DELTAP dp;

	if(svar.Bcase == 2 || svar.Bcase == 3)
	{
		cout << "Final particle count:  " << svar.nmax << endl;
	}

	svar.finPts = svar.nmax;
	pn.reserve(svar.nmax);
  	pnp1.reserve(svar.nmax);

	// if (svar.Bcase == 6){
	// 	Write_Mesh_Data(svar,cells);
	// }	

  	if(svar.restart == 1)
  	{
  		Restart(svar,fvar,avar,pn,pnp1,cells);
  	}
  	else
  	{
	  	svar.t = 0.0;				/*Total simulation time*/
  		InitSPH(svar,fvar,avar,pn,pnp1);

  		// Redefine the mass and spacing to make sure the required mass is conserved
  		if(svar.Bcase == 4 || svar.Bcase == 3)
	  		Set_Mass(svar,fvar,avar,pn,pnp1);
  	}

	// Define svar.clear to state a particle is clear of the starting area
	if (svar.Bcase == 2 || svar.Bcase == 3)
		svar.clear = -svar.Jet[1] + 4 * svar.dx;
	else
		svar.clear = 0.0;

	// Check if cells have been initialsed before making a tree off it
	if(cells.cCentre.size() == 0)
		cells.cCentre.emplace_back(StateVecD::Zero());
	if (cells.bVerts.size() == 0)
		cells.bVerts.emplace_back(StateVecD::Zero());
	if(cells.fNum.size() == 0)
	{
		cells.fNum.emplace_back();
		cells.fMass.emplace_back();
		cells.vFn.emplace_back();
		cells.vFnp1.emplace_back();
		cells.cVol.emplace_back();
		cells.cMass.emplace_back();
		cells.cRho.emplace_back();
		cells.cPertn.emplace_back();
		cells.cPertnp1.emplace_back();
	}
	
	cout << "Starting counts: " << endl;
	cout << "Boundary: " << svar.bndPts << "  Sim: " << svar.simPts << endl << endl;
	
	///********* Tree algorithm stuff ************/
	KDTREE TREE(pnp1,cells);
	// Sim_Tree NP1_INDEX(SIMDIM,pnp1,20);
	
	TREE.CELL.index->buildIndex();
	TREE.NP1.index->buildIndex();
	FindNeighbours(TREE.NP1, fvar, pnp1, outlist);

	
	// real aMom = 0.0;
	// for(size_t ii = 0; ii < cells.size(); ii++)
	// {
	// 	// Find the total momentum of the cells
	// 	aMom += cells.cMass[ii]*(cells.cVel[ii]).norm();
	// }

	// svar.aMom = aMom;
	// svar.tMom = aMom;

	///*** Perform an iteration to populate the vectors *****/
	if(svar.restart == 0)
	{
		First_Step(TREE,svar,fvar,avar,cells,outlist,dp,pnp1,pn,airP);
	}
	else
	{
		size_t const start = svar.bndPts;
		size_t const end = svar.totPts;
		dSPH_PreStep(svar,fvar,start,end,pnp1,outlist,dp);

		Detect_Surface(svar,fvar,avar,start,end,dp,outlist,cells,pnp1);

		// Apply_XSPH(fvar,start,end,outlist,dp,pnp1);
		#ifdef ALE
			Particle_Shift(svar,fvar,start,end,outlist,dp,pnp1);
		#endif
	}

	///*************** Open simulation files ***************/
	std::fstream f1,f2,f3,fb,fg;
	string framef = svar.outfolder;
	framef.append("frame.info");
	if(svar.restart == 1)
		f2.open(framef, std::ios::out | std::ios::app);
	else
		f2.open(framef, std::ios::out);

	f1 << std::scientific << std::setprecision(6);
	f3 << std::scientific << std::setw(10);

	// pertLog << std::scientific << std::setprecision(8);
	uint ghost_strand = 0;
	if(svar.boutform == 0)
		ghost_strand = 2;
	else
		ghost_strand = 3;

	Write_First_Step(f1,fb,fg,svar,pnp1,airP);
	
	if(svar.restart == 0)
	{
		if(svar.Bcase == 3 && (svar.Asource == 1 || svar.Asource == 2))
		{// Check if the pipe is inside the mesh
			cout << "Checking Pipe..." << endl;
			real holeD = svar.Jet(0)+8*svar.dx; /*Diameter of hole (or width)*/
			real r = 0.5*holeD;
#if SIMDIM == 3
			real stepb = (svar.Pstep*svar.Bstep);
	    	real dtheta = atan((stepb)/(r));
			for(real theta = 0; theta < 2*M_PI; theta += dtheta)
			{
				StateVecD xi(r*sin(theta), 0.0, r*cos(theta));
				/*Apply Rotation...*/
				xi = svar.Rotate*xi;
				xi += svar.Start;
			    if(!Check_Pipe(svar,TREE.CELL, cells, xi))
			    {

			    	cout << "Some of the pipe is outside of the simulation mesh." << endl;
			    	cout << "Fuel will be excessively close to begin." << endl;
			    	exit(-1);
			    }
	    	}
#else
	    	for(real x = -r; x <= r; x+=svar.Pstep)
	    	{
	    		StateVecD xi(x,0.0);
	    		xi = svar.Rotate*xi;
				xi += svar.Start;

				// cout << "Checking point: " << xi(0) << "  " << xi(1) << endl;
			    if(!Check_Pipe(svar,TREE.CELL, cells, xi))
			    {
			    	cout << "Some of the pipe is outside of the simulation mesh." << endl;
			    	cout << "Fuel will be excessively close to begin." << endl;
			    	exit(-1);
			    }
	    	}
#endif
		}

#if SIMDIM == 3
			if(svar.Asource == 3)
				svar.vortex.write_VLM_Panels(svar.outfolder);		
#endif
	}
	/*Timing calculation + error sum output*/
	t2 = high_resolution_clock::now();
	duration = duration_cast<microseconds>(t2-t1).count()/1e6;
	
	if(svar.restart != 1)
		svar.frame = 0;
	else
	{
		cout << "Restarting simulation..." << endl;
	}

	cout << "Frame: " << svar.frame << "  Sim Time: " << svar.t << "  Compute Time: "
	<< duration <<"  Error: " << error << endl;

	// pertLog << "Time,  Total Momentum,   Air Momentum,   Fuel Momentum,  Perturbation velocity magnitude" << endl;

	if(svar.restart == 0)
	{
		f2 << "Frame: " << svar.frame << endl;
		f2 << "Total Points: " << svar.totPts << " Boundary Points: " << svar.bndPts 
			<< " Fluid Points: " << svar.simPts << endl;
		f2 << "Sim Time:  " << std::scientific << std::setprecision(6) << svar.t 
		<< " Comp Time: " << std::fixed << duration << " Error: " << 0 << " Sub-iterations: " 
	    << 0 << endl;
		f2 << "Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;

		// Write a settings file in the solution folder.
		Write_Input_TECIO(svar,fvar,avar);
	}

	///************************* MAIN LOOP ********************/
	
#ifdef DEBUG
	const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;
	const real B = fvar.B;
	const real gam = fvar.gam;
	dbout << "Newmark Beta integration parameters" << endl;
	dbout << "a: " << a << "  b: " << b << endl;
	dbout << "c: " << c << "  d: " << d << endl;
	dbout << "B: " << B << "  gam: " << gam << endl << endl; 
#endif

	for (uint frame = 0; frame < svar.Nframe; ++frame)
	{
		int stepits=0;
		real stept=0.0;
		
		while (stept + MERROR < svar.framet)
		{
		    error = Integrate(TREE,svar,fvar,avar,cells,dp,pn,pnp1,airP,outlist);
		    stept+=svar.dt;
		    ++stepits;
		    svar.iter++;
		}
		++svar.frame;

		t2= high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;

		/*Write each frame info to file*/
		f2 << endl;
		f2 << "Frame: " << svar.frame << endl;

		f2 << "Total Points: " << svar.totPts << " Boundary Points: " << svar.bndPts 
			<< " Fluid Points: " << svar.simPts << endl;
		f2 << "Sim Time:  " << std::scientific << std::setprecision(6) << svar.t 
		<< " Comp Time: " << std::fixed << duration << " Error: " << error << " Sub-iterations: " 
	    << stepits << endl;
	    f2 << "Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;

		if(svar.outframe !=0)
		{
			if (frame % svar.outframe == 0 )
			{	/*Output to console every 20 or so steps*/
			  	cout << "Frame: " << svar.frame << "  Sim Time: " << svar.t << "  Compute Time: "
			  	<< duration <<"  Error: " << error << endl;
			  	cout << "Boundary particles:  " << svar.bndPts << " Sim particles: " << svar.totPts-svar.bndPts
			  	<< " Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;
			}
		}

		if (svar.totPts-svar.bndPts == 0) 
		{
			cout << "No more points in the simulation space. Ending...." << endl;
			break;
		}

		Write_Timestep(f1,fb,fg,ghost_strand,svar,pnp1,airP);
	}

	/*Wrap up simulation files and close them*/
	if (f1.is_open())
		f1.close();
	if (f2.is_open())
		f2.close();
	if (f3.is_open())
		f3.close();
	if (fb.is_open())
		fb.close();
	if (fg.is_open())
		fg.close();
	
	// if(svar.outtype == 0)
	// {
	// 	if(tecFileWriterClose(&svar.fuelFile))
	// 		exit(-1);

	// 	if(svar.boutform == 1)
	// 		if(tecFileWriterClose(&svar.boundFile))
	// 			exit(-1);

	// 	if(svar.gout == 1)
	// 		if(tecFileWriterClose(&svar.ghostFile))
	// 			exit(-1);
	// }

	if(svar.outtype == 0)
	{
		/*Combine the szplt files*/

		string outfile = svar.outfolder;
		outfile.append("Fuel.szplt");
		Combine_SZPLT(outfile);

		if(svar.boutform == 1)
		{
			outfile = svar.outfolder;
			outfile.append("Boundary.szplt");
			Combine_SZPLT(outfile);
		}

		if(svar.gout == 1)
		{
			outfile = svar.outfolder;
			outfile.append("Ghost.szplt");
			Combine_SZPLT(outfile);
		}
	}
	
	cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.t << " seconds" << endl;

	return 0;
}
