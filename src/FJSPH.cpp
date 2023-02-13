/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
/*** Force Calculation: On Simulating Free Surface Flows using SPH. Monaghan, J.J. (1994) ***/
/*** Continuity:        Delta-SPH. Marrone et al. (2011)                                  ***/
/*** Viscosity:         Laminar Viscosity as described by Morris (1997)                   ***/
/*** Surface Tension:   Tartakovsky and Panchenko (2016)                                  ***/
/*** Smoothing Kernel:  Wendland's C2                                                     ***/
/*** Integrator:        Newmark-Beta or 4th order Runge-Kutta                             ***/
/*** Variable Timestep Criteria: CFL + Monaghan, J.J. (1989) conditions                   ***/

#include <chrono>

#include "Var.h"
#include "Add.h"
#include "BinaryIO.h"
#include "CDFIO.h"
#include "Containment.h"
#include "FOAMIO.h"
#include "Geometry.h"
#include "H5IO.h"
#include "IPT.h"
#include "Init.h"
#include "Integration.h"
#include "IO.h"
#include "IOFunctions.h"
#include "Neighbours.h"
#include "Shifting.h"
#include "Droplet.h"
#include "Speedtest.h"
#include "VLM.h"

#ifdef DEBUG
	/*Open debug file to write to*/
	FILE* dbout;
#endif

using namespace std::chrono;
using namespace nanoflann;

int main(int argc, char *argv[])
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;
	Eigen::initParallel();
	// omp_set_num_threads(1);
	srand(unsigned(time(NULL)));
	
	#ifdef DEBUG
	dbout = fopen("WCSPH.log","w");
	#endif

	real duration;
    real error = 0;
    cout.width(13);
    cout << std::scientific << std::left << std::setprecision(4);
	
	write_header();

    /******* Define the global simulation parameters ******/
	SIM svar;
	FLUID fvar;
	AERO avar;
	LIMITS limits;
	OUTL outlist;
	MESH cells;
	VLM vortex;
	SURFS surf_marks;

	///****** Initialise the particles memory *********/
	SPHState pn;	    /*Particles at n   */
	SPHState pnp1; 	/*Particles at n+1 */
	vector<IPTState> iptdata; /* Particle tracking data */
	// SPHState airP;
	// DELTAP dp;

	GetInput(argc,argv,svar,fvar,avar,vortex);
	
	// if(MakeOutputDir(argc,argv,svar))
	// {
	// 	cout << "Couldn't make output directory. Please check permissions." << endl;
	// 	exit(-1);
	// }

	if(svar.dropDragSweep)
	{
		Droplet_Drag_Sweep(svar, fvar, avar);
		return 0; /* End after this */
	}

	if(svar.speedTest)
	{
		Speed_Test_Sweep(svar,fvar,avar);
		return 0;
	}

	if(svar.Asource == meshInfl)
	{
		#if SIMDIM == 3
			if(svar.CDForFOAM == 0)
			{
				vector<uint> empty;
				TAU::Read_BMAP(svar);
				TAU::Read_TAUMESH_FACE(svar,cells,fvar,avar);
				TAU::Read_SOLUTION(svar, fvar, avar, svar.offset_axis, cells, empty);
			}
			else
			{
				FOAM::Read_FOAM(svar,cells);
			}
		#else
			if (svar.CDForFOAM == 0)
			{
				vector<uint> used_verts;
				TAU::Read_BMAP(svar);
				TAU::Read_TAUMESH_EDGE(svar,cells,fvar,avar,used_verts);
				TAU::Read_SOLUTION(svar, fvar, avar, svar.offset_axis, cells, used_verts);
			}
			else
			{
				cout << "OpenFOAM mesh input is currently unsupported in two-dimensions" << endl;
				exit(-1);
			}
		#endif

		Init_Surface(svar,cells,surf_marks);
	}

	cout << std::setprecision(5);

	pn.reserve(svar.finPts);
  	pnp1.reserve(svar.finPts);
	
	if(!svar.restart)
	{
	  	svar.t = 0.0;				/*Total simulation time*/
		Init_Particles(svar,fvar,avar,pn,pnp1,limits);

  		// Redefine the mass and spacing to make sure the required mass is conserved
  		// if(svar.Scase == SPHERE || svar.Scase == JET)
	  	// 	Set_Mass(svar,fvar,avar,pn,pnp1);
		if(svar.Asource != meshInfl)
		{
			for(size_t ii = 0; ii < pnp1.size(); ++ii)
			{
				pn[ii].cellRho = avar.rhog;
				pn[ii].cellP = avar.pRef;
				pnp1[ii].cellRho = avar.rhog;
				pnp1[ii].cellP = avar.pRef;
				pn[ii].cellV = avar.vInf;
				pnp1[ii].cellV = avar.vInf;
			}
		}

		string file = svar.output_prefix + "_boundary.szplt.sz*";
		string cmd = "exec rm -f " + file;
		if(system(cmd.c_str()))
	    {
	    	cout << "No prexisting boundary files deleted." << endl;
	    }

		file = svar.output_prefix + "_fluid.szplt.sz*";
		cmd = "exec rm -f " + file;
		if(system(cmd.c_str()))
	    {
	    	cout << "No prexisting fluid files deleted." << endl;
	    }
  	}
	else
	{
		Init_Particles_Restart(svar,fvar,limits);
		/* Read the files */
		Read_HDF5(svar,fvar,avar,vortex,pn,pnp1,limits);
	}

	// Check if cells have been initialsed before making a tree off it
	if(cells.cCentre.size() == 0)
		cells.cCentre.emplace_back(StateVecD::Zero());

	if(cells.fNum.size() == 0)
	{
		cells.fNum.emplace_back();
		cells.fMass.emplace_back();
		cells.vFn.emplace_back();
		cells.vFnp1.emplace_back();
		cells.cRho.emplace_back();
		cells.cPertn.emplace_back();
		cells.cPertnp1.emplace_back();
	}
	
	cout << "Starting counts: " << endl;
	cout << "Boundary: " << svar.bndPts << "  Sim: " << svar.simPts << endl << endl;
	
	///********* Tree algorithm stuff ************/
	Sim_Tree SPH_TREE(SIMDIM,pnp1,20);
	Vec_Tree CELL_TREE(SIMDIM,cells.cCentre,10);
	SPH_TREE.index->buildIndex();
	CELL_TREE.index->buildIndex();
	FindNeighbours(SPH_TREE, fvar, pnp1, outlist);

	
	///*************** Open simulation files ***************/
	FILE* ff = NULL;
	FILE* f2 = NULL;
	FILE* f3 = NULL;
	FILE* fb = NULL;
	FILE* fg = NULL;
	string framef = svar.output_prefix;
	framef.append("_frame.info");
	if(svar.restart == 1)
		f2 = fopen(framef.c_str(), "a");
	else
		f2 = fopen(framef.c_str(), "w");

	if(svar.single_file)
		Write_Headers(ff,fb,fg,svar,fvar,avar);

	if(!svar.restart)
	{
		Write_Timestep(ff,fb,fg,svar,fvar,avar,limits,pnp1);
	}

	///*** Perform an iteration to populate the vectors *****/
	if(!svar.restart)
	{
		First_Step(SPH_TREE,CELL_TREE,svar,fvar,avar,vortex,cells,limits,outlist,pnp1,pn,iptdata);
	}
	else
	{
		if(svar.ghost == 1)
		{	/* Poisson points */
			PoissonGhost(svar,fvar,avar,cells,SPH_TREE,outlist,pn,pnp1);
		}
		else if (svar.ghost == 2)
		{	/* Lattice points */
			LatticeGhost(svar,fvar,avar,cells,SPH_TREE,outlist,pn,pnp1,limits);
		}

		size_t const start = svar.bndPts;
		size_t end = svar.totPts;
		size_t end_ng = svar.bndPts + svar.simPts;
		real npd = 1.0;
		dSPH_PreStep(fvar,end,pnp1,outlist,npd);

		if (svar.Asource == meshInfl)
		{
			// cout << "Finding cells" << endl;
			FindCell(svar,avar,CELL_TREE,cells,pn,pnp1);
			// if (svar.totPts != pnp1.size())
			// {	//Rebuild the neighbour list
			// 	// cout << "Updating neighbour list" << endl;
			// 	// cout << "Old: " << svar.totPts << "  New: " << pnp1.size() << endl;
			// 	svar.delNum += svar.totPts-pnp1.size();
			// 	svar.simPts -= svar.totPts-pnp1.size();
			// 	svar.totPts = pnp1.size();
			// 	end = svar.totPts;
			// 	end_ng = svar.bndPts + svar.simPts;
			// 	SPH_TREE.index->buildIndex();
			// 	FindNeighbours(SPH_TREE, fvar, pnp1, outlist);
			// 	dSPH_PreStep(svar,fvar,svar.totPts,pnp1,outlist,dp);
			// }	
		}

		Detect_Surface(svar,fvar,avar,start,end_ng,outlist,cells,vortex,pn);

		// Apply_XSPH(fvar,start,end,outlist,dp,pnp1);
		#ifdef ALE
			if(svar.ghost > 0)
				Particle_Shift_Ghost(svar,fvar,start,end_ng,outlist,pn);
			else
				Particle_Shift_No_Ghost(svar,fvar,start,end_ng,outlist,pn);
		#endif

		/* Update shifting velocity and surface data */
		// SPHState::const_iterator first = pnp1.begin();
		// SPHState::const_iterator last = pnp1.begin() + end_ng;
		// pn = SPHState(first,last);
	}

	// Append_Restart_Prefix(svar);

	if(svar.using_ipt)
	{
		IPT::Init_IPT_Files(svar);
		IPT::Write_Data(svar,cells,iptdata);
	}

	// Make sure that inlets are injecting inside the mesh domain.
	#if SIMDIM == 3
	if(!svar.restart && svar.Asource == 3)
		vortex.write_VLM_Panels(svar.output_prefix);		
	#endif

	/*Timing calculation + error sum output*/
	t2 = high_resolution_clock::now();
	duration = duration_cast<seconds>(t2-t1).count();
	
	if(!svar.restart)
	{
		svar.frame = 0;
		fprintf(f2,"Frame: %u\n",svar.frame);
		fprintf(f2, "Total Points: %zu Boundary Points: %zu Fluid Points: %zu\n", 
				svar.totPts, svar.bndPts, svar.simPts);
		fprintf(f2, "Sim Time:  %.7g Comp Time: %.6e Error: %.6f Sub-iterations: %d\n", svar.t, duration, 0.0, 0);
		fprintf(f2, "Deleted particles: %zu Internal collisions: %zu\n", svar.delNum, svar.intNum);
	}
	else
	{
		cout << "Restarting simulation..." << endl;	
		cout << "Frame: " << svar.frame << "  Sim Time: " << svar.t << "  Compute Time: "
				<< duration <<"  Error: " << error << endl;
	}

	///************************* MAIN LOOP ********************/
	
	#ifdef DEBUG
		const real a = 1 - svar.gamma;
		const real b = svar.gamma;
		const real c = 0.5*(1-2*svar.beta);
		const real d = svar.beta;
		const real B = fvar.B;
		const real gam = fvar.gam;
		fprintf(dbout,"Newmark Beta integration parameters\n");
		fprintf(dbout, "a: %f b: %f\n", a, b);
		fprintf(dbout, "c: %f d: %f\n", c, d);
		fprintf(dbout, "B: %f gam: %f\n\n", B, gam); 
	#endif

	
	for (uint frame = 0; frame < svar.Nframe; ++frame)
	{
		int stepits=0;
		real stept=0.0;
		while (stept + 0.1*svar.dt_min < svar.framet)
		{
			if(stepits % 50 == 0)
				#ifdef ALE
				printf("\nTime      | Timestep | CFL  | RMS error | its | dRho (%%) | Max-F     | Max-Af    | Max Shift | Step time (ms)| \n");
				#else
				printf("\nTime      | Timestep | CFL  | RMS error | its | dRho (%%) | Max-F     | Max-Af    | Step time (ms)| \n");
				#endif

			error = Integrate(SPH_TREE,CELL_TREE,svar,fvar,avar,vortex,cells,surf_marks,limits,outlist,pn,pnp1,iptdata);
			stept+=svar.dt;
			++stepits;
		}
		++svar.frame;

		t2= high_resolution_clock::now();
		duration = duration_cast<microseconds>(t2-t1).count()/1e6;

		/*Write each frame info to file*/
		fprintf(f2,"\nFrame: %u\n",svar.frame);
		fprintf(f2, "Total Points: %zu Boundary Points: %zu Fluid Points: %zu\n", 
				svar.totPts, svar.bndPts, svar.simPts);
		fprintf(f2, "Sim Time:  %.7g Comp Time: %.6e Error: %.6f Sub-iterations: %d\n", svar.t, duration, error, stepits);
		fprintf(f2, "Deleted particles: %zu Internal collisions: %zu\n", svar.delNum, svar.intNum);

		cout << "Frame: " << svar.frame << "  Sim Time: " << svar.t << "  Compute Time: "
			 << duration <<"  Error: " << error << endl;
		cout << "Boundary particles:  " << svar.bndPts << " Sim particles: " << svar.totPts-svar.bndPts
			 << " Deleted particles: " << svar.delNum << " Internal collisions: " << svar.intNum <<  endl;
			
		if (svar.totPts-svar.bndPts == 0) 
		{
			cout << "No more points in the simulation space. Ending...." << endl;
			break;
		}

		Write_Timestep(ff,fb,fg,svar,fvar,avar,limits,pnp1);
		if(svar.using_ipt)
			IPT::Write_Data(svar,cells,iptdata);

		svar.tframem1 += svar.framet; /* March frame time forward */
	}

	/*Wrap up simulation files and close them*/
	if (f2 != NULL)
		fclose(f2);
	if (f3 != NULL)
		fclose(f3);
	if (ff != NULL)
		fclose(ff);
	if (fb != NULL)
		fclose(fb);
	if (fg != NULL)
		fclose(fg);
	
	if(svar.out_encoding == 1)
	{
		/*Combine the szplt files*/

		string outfile = svar.output_prefix + "_fluid.szplt";
		Combine_SZPLT(outfile);

		outfile = svar.output_prefix + "_boundary.szplt";
		Combine_SZPLT(outfile);
		
		if(svar.gout == 1)
		{
			outfile = svar.output_prefix + "_ghost.szplt";
			Combine_SZPLT(outfile);
		}
	}
	
	cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.t << " seconds" << endl;

	return 0;
}
