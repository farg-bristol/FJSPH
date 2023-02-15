/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "AsciiIO.h"
#include "BinaryIO.h"
#include "Geometry.h"
#include "H5IO.h"
#include "IO.h"
#include "IOFunctions.h"
#include "Kernel.h"
#include "VLM.h"
#include <ctime>
#include <Eigen/LU>
#include <filesystem>
#include <regex>

using std::filesystem::directory_iterator;

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void Set_Values(SIM& svar, FLUID& fvar, AERO& avar, VLM& vortex)
{
	/*Universal parameters based on input values*/
	// StateVecD angle = svar.Angle;
	// angle *= M_PI/180.0;
	// svar.Rotate = GetRotationMat(angle);
	// svar.Transp = svar.Rotate.transpose();
	
	svar.offset_vec *= svar.scale;
	svar.max_x *= svar.scale;
	svar.max_x_sph *= svar.scale;

	// if(avar.vJetMag != -1)
	// {
	// 	avar.vStart = StateVecD::Zero();
	// 	avar.vStart[1] = avar.vJetMag;
	// 	avar.vStart = svar.Rotate*avar.vStart;
	// }

	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Cole's Eq*/

	/*Pipe Pressure calc*/
	fvar.rhoJ = density_equation(fvar.pPress,fvar.B,fvar.gam,fvar.Cs,fvar.rho0,fvar.backP);
	
	/* Upper and lower limits for density */
	if(fvar.rhoMax == 1500 && fvar.rhoMin == 500)
	{	/* If limits are undefined, use a variation around the base density */
		fvar.rhoMax = fvar.rho0*(1.0 + fvar.rhoVar*0.01);
		fvar.rhoMin = fvar.rho0*(1.0 - fvar.rhoVar*0.01);
	}
	
	// svar.nrad = 1;
	// if(svar.Scase == BOX)
	// {
	// 	svar.dx = svar.Pstep;
	// }
	// else if(svar.Scase == CYLINDER || svar.Scase == JET || svar.Scase == CONVJET)
	// {	/* Cylinder case */
	// 	// Defining dx to fit the pipe, then find rest spacing
	// 	if(svar.jet_diam/svar.Pstep < 1.5)
	// 	{	/*spacing is too close to full size to use standard adjustment*/
	// 		printf("Warning: particle spacing if of the same order of magintude of the jet diameter.\n");
	// 		printf("Consider a different size for accuracy.\n");
	// 	}
	// 	else
	// 	{
	// 		svar.nrad = ceil(abs(0.5*svar.jet_diam/svar.Pstep));
	// 		svar.dx = 0.5*(svar.jet_diam)/real(svar.nrad);
	// 	}	

	// 	// avar.dPipe = svar.Jet(0);
	// }
	// else if (svar.Scase == SPHERE)
	// {	/* Droplet case */
	// 	if(svar.diam/svar.Pstep < 1.5)
	// 	{	/*spacing is too close to full size to use standard adjustment*/
	// 		svar.nrad = 1;
	// 		svar.dx = svar.diam;
	// 	}
	// 	else
	// 	{
	// 		real radius = 0.5*svar.diam;
		
	// 		svar.nrad = ceil(abs(radius/svar.Pstep));
	// 		svar.dx = radius/real(svar.nrad);
	// 	}
	// }
	
	svar.dx = svar.Pstep * pow(fvar.rhoJ/fvar.rho0,1.0/SIMDIM);
 	
	// Correct the droplet to have the same volume as the original
	// if(svar.Scase == SPHERE)
	// {
	// 	printf("Droplet Diameter: ", svar.diam, endl;		
	// }

	#if SIMDIM == 3
		avar.pVol = 4.0/3.0 * M_PI * pow(svar.Pstep*0.5,SIMDIM);
	#else
		avar.pVol = M_PI* svar.Pstep*svar.Pstep/4.0;
	#endif

 	svar.beta = 0.25; svar.gamma = 0.5; /*Newmark Beta parameters*/

 	/*Mass from spacing and density*/
	fvar.simM = fvar.rho0*pow(svar.Pstep,SIMDIM); 
	fvar.bndM = fvar.simM;
	avar.gasM = avar.rhog*pow(svar.Pstep,SIMDIM);

	avar.sos = sqrt(avar.T*avar.Rgas*avar.gamma);
	avar.isossqr = 1.0/(avar.sos*avar.sos);
	
	if(avar.MRef != -1)
	{
		/* Mach has been defined, not velocity */
		avar.vRef = avar.sos * avar.MRef;
	}
	else
	{
		avar.MRef = avar.vRef/avar.sos;
	}
		
	svar.addcount = 0;
	
	if(svar.dt_min > 0)
		svar.dt = svar.dt_min;
	else
		svar.dt = 2E-010; 			/*Initial timestep*/

  	fvar.H = fvar.Hfac*svar.Pstep;
  	fvar.HSQ = fvar.H*fvar.H; 
	fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/

	fvar.dCont = 2.0 * fvar.delta * fvar.H * fvar.Cs;
  	// fvar.dMom = fvar.alpha * fvar.H * fvar.Cs * fvar.rho0;
  	fvar.dMom = 2.0*(SIMDIM + 2.0);
    fvar.nu = fvar.mu/fvar.rho0;

	#if SIMDIM == 2
		#ifdef CUBIC
			fvar.correc = 10.0 / (7.0 * M_PI * fvar.H * fvar.H);
		#else
			fvar.correc = 7.0 / (4.0 * M_PI * fvar.H * fvar.H);
		#endif
	#endif
	#if SIMDIM == 3
		#ifdef CUBIC
			fvar.correc = (1.0/(M_PI*fvar.H*fvar.H*fvar.H));
		#else
			fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
		#endif
	#endif
		
	fvar.Wdx = Kernel(svar.Pstep,fvar.H,fvar.correc);

	#if SIMDIM == 3
		if (svar.Asource == VLMInfl)
		{
			vortex.Init(svar.vlm_file);
			vortex.GetGamma(avar.vInf);
			vortex.write_VLM_Panels(svar.output_prefix);
			if(vortex.write_traj)
				vortex.Plot_Streamlines(svar.output_prefix);
		}
	#endif

	avar.GetYcoef(fvar, /*fvar.H*/ svar.Pstep);
	real nfull = get_n_full(svar.Pstep, fvar.H);
	avar.nfull = nfull;
	avar.infull = 1.0/avar.nfull;
	avar.interp_fac = 1.0/avar.iinterp_fac;

	#if SIMDIM == 3
		avar.aPlate = svar.Pstep * svar.Pstep;
		// avar.aPlate = fvar.H*fvar.H;
	#else
		avar.aPlate = svar.Pstep /**svar.Pstep*/ /** pow(avar.L,0.5)*/;
		// avar.aPlate = fvar.H;
	#endif

	/* Particle tracking values */
	svar.IPT_diam = pow((6.0*fvar.simM)/(M_PI*fvar.rho0),1.0/3.0);
	svar.IPT_area = M_PI * svar.IPT_diam*svar.IPT_diam/4.0; 
}

void print_vector(string const& pretext, StateVecD const& vec)
{
	#if SIMDIM == 3
	printf("%s: %f %f %f\n",pretext.c_str(), vec[0],vec[1],vec[2]);
	#else
	printf("%s: %f %f\n",pretext.c_str(), vec[0],vec[1]);
	#endif
}

void Print_Settings(FILE* out, SIM const& svar, FLUID const& fvar, AERO const& avar)
{
    fprintf(out, "Input values, after modification and interpretation, are: \n\n");

    #pragma omp parallel
    {
        #pragma omp single
        fprintf(out, "                         Number of threads: %d\n", omp_get_num_threads());
    }
    
    /* File Inputs */
    fprintf(out, " Input files --------------------------------: -\n");
    fprintf(out, "              Input fluid definition filename: %s\n", svar.fluidfile.c_str());
    fprintf(out, "           Input boundary definition filename: %s\n", svar.boundfile.c_str());
    fprintf(out, "                           SPH restart prefix: %s\n", svar.restart_prefix.c_str());
    fprintf(out, "                      VLM definition filename: %s\n\n", svar.vlm_file.c_str());
    
    /* TAU files */
    fprintf(out, " TAU files ----------------------------------: -\n");
    fprintf(out, "                   Primary grid face filename: %s\n", svar.taumesh.c_str());
    fprintf(out, "                    Boundary mapping filename: %s\n", svar.taubmap.c_str());
    fprintf(out, "                          Restart-data prefix: %s\n", svar.tausol.c_str());
    fprintf(out, "                                   Grid scale: %f\n", svar.scale);
    fprintf(out, "                         Angle alpha (degree): %f\n", svar.angle_alpha);
    fprintf(out, "           2D offset vector (0 / x=1,y=2,z=3): %d\n\n", svar.offset_axis);

    /* OpenFOAM files */
    fprintf(out, " OpenFOAM files -----------------------------: -\n");
    fprintf(out, "                     OpenFOAM input directory: %s\n", svar.foamdir.c_str());
    fprintf(out, "                  OpenFOAM solution directory: %s\n", svar.foamsol.c_str());
    fprintf(out, "                        OpenFOAM binary (0/1): %d\n", svar.isBinary);
    fprintf(out, "                           Label size (32/64): %d\n", svar.labelSize);
    fprintf(out, "                          Scalar size (32/64): %d\n", svar.scalarSize);
    fprintf(out, "                       OpenFOAM buoyant (0/1): %d\n\n", svar.buoyantSim);

    /* File outputs */
    fprintf(out, " Output parameters --------------------------: -\n");
    fprintf(out, "                 Single file for output (0/1): %u\n", svar.single_file);
	fprintf(out, "                   Write Tecplot output (0/1): %u\n", svar.write_tecio);
	fprintf(out, "                    Write H5Part output (0/1): %u\n", svar.write_h5part);
    fprintf(out, "                          Output files prefix: %s\n", svar.output_prefix.c_str());
    fprintf(out, "                      SPH frame time interval: %f\n", svar.framet);
    fprintf(out, "                              SPH frame count: %u\n", svar.Nframe);
    fprintf(out, "       SPH output encoding (0=ascii/1=binary): %u\n", svar.out_encoding);
    fprintf(out, "                       SPH ghost output (0/1): %u\n", svar.gout);
    fprintf(out, "                                Variable list: %s\n\n", svar.output_names.c_str());
    
    /* Fluid data */
    fprintf(out, " Fluid parameters ---------------------------: -\n");
    fprintf(out, "                            Reference density: %g\n", avar.rhog);
    fprintf(out, "                  Reference dispersed density: %g\n", fvar.rho0);
    fprintf(out, "                Reference dispersed viscosity: %g\n", fvar.mu);
    fprintf(out, "                    Reference surface tension: %g\n", fvar.sig);
    fprintf(out, "            SPH surface tension contact angle: %g\n", fvar.contangb);
    fprintf(out, "              SPH artificial viscosity factor: %g\n", fvar.alpha);
    fprintf(out, "                        SPH delta coefficient: %g\n", fvar.delta);
    fprintf(out, "                SPH maximum shifting velocity: %g\n", svar.maxshift);
    fprintf(out, "                           SPH speed of sound: %g\n", fvar.Cs);
    fprintf(out, "                      SPH background pressure: %g\n", fvar.backP);
    fprintf(out, "                        SPH starting pressure: %g\n", fvar.pPress);
    fprintf(out, "                    SPH density variation (%%): %g\n", fvar.rhoVar);
    fprintf(out, "                          SPH maximum density: %g\n", fvar.rhoMax);
    fprintf(out, "                          SPH minimum density: %g\n", fvar.rhoMin);
    fprintf(out, "              Init hydrostatic pressure (0/1): %d\n", svar.init_hydro_pressure);
    fprintf(out, "                           Hydrostatic height: %g\n\n", svar.hydro_height);

    /* Simulation settings */
    fprintf(out, " Time integration settings ------------------: -\n");
    fprintf(out, "                       SPH integration solver: %s\n", svar.solver_name.c_str());
    fprintf(out, "     SPH boundary solver (0=pressure/1=ghost): %u\n", svar.bound_solver);
    fprintf(out, "                  SPH solver minimum residual: %g\n", svar.minRes);
    fprintf(out, "                         SPH maximum timestep: %g\n", svar.dt_max);
    fprintf(out, "                         SPH minimum timestep: %g\n", svar.dt_min);
    fprintf(out, "                              SPH maximum CFL: %g\n", svar.cfl_max);
    fprintf(out, "                              SPH minimum CFL: %g\n", svar.cfl_min);
    fprintf(out, "                            SPH CFL condition: %g\n", svar.cfl);
    fprintf(out, "                        SPH unstable CFL step: %g\n", svar.cfl_step);
    fprintf(out, "             SPH Newmark-Beta iteration limit: %u\n", svar.subits);
    fprintf(out, "                 SPH unstable CFL count limit: %u\n", svar.nUnstable_Limit);
    fprintf(out, "                   SPH stable CFL count limit: %u\n", svar.nStable_Limit);
    fprintf(out, "        SPH stable CFL count iteration factor: %g\n", svar.subits_factor);
    fprintf(out, "                   SPH maximum particle count: %zu\n\n", svar.finPts);
    
    fprintf(out, " Geometric parameters -----------------------: -\n");
    #if SIMDIM == 3
    fprintf(out, "                           SPH gravity vector: %g, %g, %g\n", svar.grav[0],svar.grav[1],svar.grav[2]);
    #else
    fprintf(out, "                           SPH gravity vector: %g, %g\n", svar.grav[0],svar.grav[1]);
    #endif
    fprintf(out, "                          SPH initial spacing: %g\n", svar.Pstep);
    fprintf(out, "                  SPH boundary spacing factor: %g\n", svar.Bstep);
    fprintf(out, "                    SPH restart fit tolerance: %f\n", svar.restart_tol);
    fprintf(out, "                  SPH smoothing length factor: %g\n", fvar.Hfac);
    #if SIMDIM == 3
    fprintf(out, "                 SPH global offset coordinate: %g, %g, %g\n\n", svar.offset_vec[0],svar.offset_vec[1],svar.offset_vec[2]);
    #else
    fprintf(out, "                 SPH global offset coordinate: %g, %g\n\n", svar.offset_vec[0],svar.offset_vec[1]);
    #endif

    fprintf(out, " Aerodynamic coupling parameters ------------: -\n");
    fprintf(out, "                         SPH aerodynamic case: %s\n", avar.aero_case.c_str());
    fprintf(out, "                 SPH aerodynamic cutoff value: %f\n", avar.cutoff);
	fprintf(out, "          SPH aerodynamic iterpolation factor: %f\n", avar.iinterp_fac);
	fprintf(out, "           SPH aerodynamic full neighbourhood: %f\n", avar.nfull);
    fprintf(out, " SPH interpolation factor (0=ncount/1=lambda): %d\n", avar.use_lam);
    fprintf(out, "        SPH SP diameter definition (0=dx/1=h): %d\n", avar.use_dx);
    fprintf(out, "                SPH use TAB deformation (0/1): %d\n", avar.useDef);
    fprintf(out, "                SPH use ghost particles (0/1): %d\n", svar.ghost);
    #if SIMDIM == 3
    fprintf(out, "                      SPH freestream velocity: %g, %g, %g\n", avar.vInf[0],avar.vInf[1],avar.vInf[2]);
    #else
    fprintf(out, "                      SPH freestream velocity: %g, %g\n\n", avar.vInf[0],avar.vInf[1]);
    #endif
	
    /* Aerodynamic data */
    fprintf(out, " Carrier phase parameters -------------------: -\n");
    fprintf(out, "                           Reference velocity: %g\n", avar.vRef);
    fprintf(out, "                           Reference pressure: %g\n", avar.pRef);
    fprintf(out, "                        Reference Mach number: %g\n", avar.MRef);
    fprintf(out, "                            Reference density: %g\n", avar.rhog);
    fprintf(out, "                        Reference temperature: %g\n", avar.T);
    fprintf(out, "                           Gas constant gamma: %g\n", avar.gamma);
    fprintf(out, "               Sutherland reference viscosity: %g\n\n", avar.mug);

    /* Particle tracking settings */
    fprintf(out, " Particle tracking parameters ---------------: -\n");
    fprintf(out, "                      Transition to IPT (0/1): %d\n", svar.using_ipt);
    fprintf(out, "                Velocity equation order (1/2): %d\n", svar.eqOrder);
    fprintf(out, "         SPH tracking conversion x coordinate: %g\n", svar.max_x_sph);
    fprintf(out, "              Maximum x trajectory coordinate: %g\n", svar.max_x);
    fprintf(out, "              Particle scatter output (0/1/2): %u\n", svar.partout);
    fprintf(out, "               Particle streak output (0/1/2): %u\n", svar.streakout);
    fprintf(out, "    Particle cell intersection output (0/1/2): %u\n\n", svar.cellsout);

    /* Droplet drag sweep settings */
    fprintf(out, " Droplet drag parameters --------------------: -\n");
    fprintf(out, "                  Do droplet drag sweep (0/1): %d\n", svar.dropDragSweep);
    fprintf(out, "                          Do speed test (0/1): %d\n", svar.speedTest);
    fprintf(out, "                         Speed test run count: %d\n", svar.nRuns);

    string str;
    for(auto const& x : svar.nacross)
    str.append(std::to_string(x) + ",");
    str.pop_back(); //Remove the last comma
    fprintf(out, "                          Droplet resolutions: %s\n", str.c_str());

    str.clear();
    for(auto const& x : svar.diameters)
    str.append(std::to_string(x) + ",");
    str.pop_back(); //Remove the last comma
    fprintf(out, "                            Droplet diameters: %s\n", str.c_str());

        str.clear();
    for(auto const& x : svar.velocities)
    str.append(std::to_string(x) + ",");
    str.pop_back(); //Remove the last comma
    fprintf(out, "                           Droplet velocities: %s\n", str.c_str());
    
    str.clear();
    for(auto const& x : svar.Reynolds)
    str.append(std::to_string(x) + ",");
    str.pop_back(); //Remove the last comma
    fprintf(out, "                     Droplet Reynolds numbers: %s\n\n", str.c_str());
	fflush(out);
}


void GetInput(int argc, char **argv, SIM& svar, FLUID& fvar, AERO& avar, VLM& vortex)
{
	if (argc > 2) 
	{	/*Check number of input arguments*/
		printf("WARNING: only a maximum of one input arguments accepted,\n");
		printf("\t1: Input file directory\n");
		printf("Other inputs will be ignored.\n");
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	printf("ERROR: No inputs provided. Stopping... \n");
    	exit(-1);    	
    }

	/*Get parameters if it has been provided*/
	// printf(argv[1], endl;
	string file = argv[1];
	    	
	svar.infile = file;

	FILE* fin = fopen(argv[1], "r");
	if(fin == NULL)
	{
		printf("ERROR: Failed to open parameter file: %s\n", argv[1]);
		exit(-1);
	}
	else
	{
		printf("Para file open, reading contents...\n");
	}

    char* cline = NULL;
	size_t cline_buf_size = 0;
	ssize_t cline_size = getline(&cline, &cline_buf_size, fin);
    while (cline_size >= 0)
    {
		string line(cline); // Recast as a string to manipulate.
		line = ltrim(line);
        size_t end = line.find_first_of('#');
		if(end != std::string::npos)
			line = line.substr(0,end+1);

        /* File Inputs */
		Get_String(line, "Input fluid definition filename", svar.fluidfile);
		Get_String(line, "Input boundary definition filename", svar.boundfile);
        Get_String(line, "Primary grid face filename", svar.taumesh);
        Get_String(line, "Boundary mapping filename", svar.taubmap);
        Get_String(line, "Restart-data prefix", svar.tausol);
        Get_String(line, "SPH restart prefix", svar.restart_prefix);
        Get_Number(line, "Grid scale", svar.scale);
		Get_Number(line, "Angle alpha (degree)", svar.angle_alpha);
        Get_Number(line, "2D offset vector (0 / x=1,y=2,z=3)",svar.offset_axis);
        Get_String(line, "OpenFOAM input directory", svar.foamdir);
        Get_String(line, "OpenFOAM solution directory", svar.foamsol);
        Get_Number(line, "OpenFOAM binary (0/1)", svar.isBinary);
        Get_Number(line, "Label size (32/64)", svar.labelSize);
        Get_Number(line, "Scalar size (32/64)", svar.scalarSize);
        Get_Number(line, "OpenFOAM buoyant (0/1)", svar.buoyantSim);
		#if SIMDIM == 3
		Get_String(line, "VLM definition filename", svar.vlm_file); 
		#endif 

        /* File outputs */
		Get_Number(line, "Write Tecplot output (0/1)", svar.write_tecio);
		Get_Number(line, "Write H5Part output (0/1)", svar.write_h5part);
		Get_Number(line, "Single file for output (0/1)",svar.single_file);
		Get_String(line, "Output files prefix", svar.output_prefix);
        Get_Number(line, "SPH frame time interval", svar.framet);
        Get_Number(line, "SPH frame count", svar.Nframe);
        Get_Number(line, "SPH output encoding (0=ascii/1=binary)", svar.out_encoding);
        Get_Number(line, "SPH ghost output (0/1)", svar.gout);
        Get_String(line, "Variable list", svar.output_names);
        
        // Get_String(line, "Particle surface impact filename", svar.surfacefile);

        /* Fluid data */
        Get_Number(line, "Reference dispersed density", fvar.rho0);
        Get_Number(line, "Sutherland reference viscosity", avar.mug);
        Get_Number(line, "Reference dispersed viscosity", fvar.mu);
        Get_Number(line, "Reference surface tension", fvar.sig);
		Get_Number(line, "SPH surface tension contact angle", fvar.contangb);
        Get_Number(line, "Init hydrostatic pressure (0/1)", svar.init_hydro_pressure);
		Get_Number(line, "Hydrostatic height", svar.hydro_height);

		/* Aerodynamic data */

		/* Simulation settings */
		Get_String(line, "SPH integration solver", svar.solver_name);
		Get_Number(line, "SPH boundary solver (0=pressure/1=ghost)", svar.bound_solver);
		Get_Number(line, "SPH solver minimum residual", svar.minRes);
		Get_Number(line, "SPH maximum timestep", svar.dt_max);
		Get_Number(line, "SPH minimum timestep", svar.dt_min);
        Get_Number(line, "SPH maximum CFL", svar.cfl_max);
        Get_Number(line, "SPH minimum CFL", svar.cfl_min);
        Get_Number(line, "SPH CFL condition", svar.cfl);
        Get_Number(line, "SPH unstable CFL step", svar.cfl_step);
        Get_Number(line, "SPH unstable CFL count limit", svar.nUnstable_Limit);
        Get_Number(line, "SPH stable CFL count limit", svar.nStable_Limit);
        Get_Number(line, "SPH stable CFL count iteration factor", svar.subits_factor);
        Get_Number(line, "SPH maximum shifting velocity", svar.maxshift);

		Get_Number(line, "SPH background pressure", fvar.backP);
		Get_Number(line, "SPH starting pressure", fvar.pPress);
		Get_Number(line, "SPH density variation (%)", fvar.rhoVar);
		Get_Number(line, "SPH maximum density", fvar.rhoMax);
		Get_Number(line, "SPH minimum density", fvar.rhoMin);
		Get_Number(line, "SPH delta coefficient", fvar.delta);

		Get_Number(line, "SPH artificial viscosity factor", fvar.alpha);
		Get_Number(line, "SPH speed of sound", fvar.Cs);
        Get_Number(line, "SPH Newmark-Beta iteration limit", svar.subits);
		Get_Vector(line, "SPH gravity vector", svar.grav);

        Get_Number(line, "SPH initial spacing", svar.Pstep);
        Get_Number(line, "SPH boundary spacing factor", svar.Bstep);
        Get_Number(line, "SPH smoothing length factor", fvar.Hfac);
        Get_Number(line, "SPH use ghost particles (0/1)", svar.ghost);
        Get_Vector(line, "SPH global offset coordinate", svar.offset_vec);
		Get_Number(line, "SPH maximum particle count",svar.finPts);
        Get_Number(line, "SPH restart fit tolerance", svar.restart_tol);

		/* Aerodynamic coupling settings */
        Get_String(line, "SPH aerodynamic case", avar.aero_case);
        Get_Number(line, "SPH SP diameter definition (0=dx/1=h)", avar.use_dx);
        Get_Number(line, "SPH use TAB deformation (0/1)", avar.useDef);
        Get_Number(line, "SPH interpolation factor (0=ncount/1=lambda)", avar.use_lam);
        Get_Number(line, "SPH aerodynamic cutoff value", avar.cutoff);
        Get_Number(line, "SPH aerodynamic interpolation factor", avar.iinterp_fac);
        Get_Vector(line, "SPH freestream velocity", avar.vInf);
        Get_Number(line, "Reference velocity", avar.vRef);
        Get_Number(line, "Reference pressure", avar.pRef);
        Get_Number(line, "Reference Mach number", avar.MRef);
        Get_Number(line, "Reference density", avar.rhog);
        Get_Number(line, "Reference temperature", avar.T);
		Get_Number(line, "Gas constant gamma", avar.gamma);

        /* Particle tracking settings */
		Get_Number(line, "Transition to IPT (0/1)", svar.using_ipt);
        Get_Number(line, "Velocity equation order (1/2)", svar.eqOrder);
        Get_Number(line, "SPH tracking conversion x coordinate", svar.max_x_sph);
        Get_Number(line, "Maximum x trajectory coordinate", svar.max_x);
		Get_Number(line, "Particle scatter output (0/1/2)", svar.partout);
        Get_Number(line, "Particle streak output (0/1/2)", svar.streakout);
        Get_Number(line, "Particle cell intersection output (0/1/2)", svar.cellsout);

		/* Droplet drag sweep settings */
		Get_Number(line, "Do droplet drag sweep (0/1)", svar.dropDragSweep);
        Get_Array(line,  "Droplet resolutions", svar.nacross);
		Get_Array(line,  "Droplet diameters", svar.diameters);
		Get_Array(line,  "Droplet velocities", svar.velocities);
		Get_Array(line,  "Droplet Reynolds numbers", svar.Reynolds);

		/* Speed test settings */
		Get_Number(line, "Do speed test (0/1)", svar.speedTest);
		Get_Number(line, "Speed test run count", svar.nRuns);

		cline_size = getline(&cline, &cline_buf_size, fin);
    }

	free(cline);
	cline = NULL;
    fclose(fin);

	int fault = 0;

	/* Need to check if inputs are correct */
	if(svar.taumesh.empty())
    {
		if(svar.foamdir.empty())
		{
			#if SIMDIM == 3
			if(svar.vlm_file.empty())
				svar.Asource = constVel;
			else
			{
				svar.Asource = VLMInfl;

				if(svar.vlm_file == "(thisfile)")
					svar.vlm_file = svar.infile;
			}

			#else
			svar.Asource = constVel;
			#endif
		}
		else
		{
			if(svar.foamsol.empty())
			{
				printf("OpenFOAM solution directory not defined.\n");
				fault = 1;
			}

			svar.CDForFOAM = 1;
			svar.Asource = meshInfl;			
		}
    }
	else
	{
		svar.CDForFOAM = 0;
		svar.Asource = meshInfl;
		if(svar.taubmap.empty())
		{
			printf("Input TAU bmap file not defined.\n");
			fault = 1;
		}
		else if(svar.taubmap == "(thisfile)")
		{
			/* The para file is the boundary definition file, so set it so */
			svar.taubmap = svar.infile;
		}
		else
		{
			// Check files exist
			if(!std::filesystem::exists(svar.taubmap))
			{   
				printf("Input TAU boundary file \"%s\" not found.\n",svar.taubmap.c_str());
				fault = 1;
			}

			if(!std::filesystem::exists(svar.taumesh))
			{   
				printf("Input TAU mesh file \"%s\" not found.\n",svar.taumesh.c_str());
				fault = 1;
			}
		}
		
		if(svar.tausol.empty())
		{
			printf("Input TAU solution file not defined.\n");
			fault = 1;
		}
		else
		{			
			if(!std::filesystem::exists(svar.tausol))
			{   
				printf("Input TAU solution file \"%s\" not found.\n",svar.tausol.c_str());
				fault = 1;
			}
		}

	}
	
	/* Check for restart data now boundary case is known */	
	if(!svar.restart_prefix.empty())
	{
		svar.restart = 1;
		
		if(!std::filesystem::exists(svar.restart_prefix + "_particles.h5"))
		{   
			printf("SPH restart file \"%s\" not found.\n",(svar.restart_prefix + "_particles.h5").c_str());
			fault = 1;
		}
		// Check_If_Restart_Possible(svar);
		// svar.output_prefix = svar.restart_prefix;
	}

	if(svar.fluidfile.empty())
	{
		printf("Input SPH fluid file not defined.\n");
		fault = 1;
	}
	else
	{
		if(!std::filesystem::exists(svar.fluidfile))
		{   
			printf("SPH fluid file \"%s\" not found.\n",svar.fluidfile.c_str());
			fault = 1;
		}
	}

	if(svar.boundfile.empty())
	{
		printf("Input SPH boundary file not defined.\n");
		fault = 1;
	}
	else
	{
		if(!std::filesystem::exists(svar.boundfile))
		{   
			printf("SPH boundary file \"%s\" not found.\n",svar.boundfile.c_str());
			fault = 1;
		}
	}

	if(!svar.solver_name.empty())
	{
		if(svar.solver_name == "Newmark-Beta")
		{
			svar.solver_type = newmark_beta;
		}
		else if(svar.solver_name == "Runge-Kutta")
		{
			svar.solver_type = runge_kutta;
		}
		else
		{
			printf("ERROR: Unrecognised solver name. Choose from the following options:\n");
			printf("\t1. Newmark-Beta\n\t2. Runge-Kutta\n");
		fault = 1;
		}
	}

	if(svar.Pstep < 0)
	{
		printf("ERROR: SPH initial spacing has not been defined.\n");
		fault = 1;
	}

	if(svar.Nframe < 0)
	{
		printf("Number of frames to output not defined.\n");
		fault = 1;
	}

	if (svar.framet < 0)
	{
		printf("Frame time interval has not been defined.\n");
		fault = 1;
	}

	if(svar.init_hydro_pressure)
	{
		if(svar.hydro_height < 0)
		{
			printf("Hydrostatic height has not been defined.\n");
			fault = 1;
		}
	}

	/* Aerodynamic settings */
	if(avar.aero_case == "(none)")
	{
		avar.acase = none;
	}
	else if(avar.aero_case == "Gissler")
	{
		avar.acase = Gissler;
	}
	else if (avar.aero_case == "Induced_pressure")
	{
		avar.acase = Induced_Pressure;
	}
	else if (avar.aero_case == "Skin_friction")
	{
		avar.acase = SkinFric;
	}
	else
	{
		printf("Aerodynamic coupling model is not defined or correct.\n");
		fault = 1;
	}

	if(svar.dt_min > 0.0)
	{
		svar.dt = svar.dt_min;
	}

	if(avar.iinterp_fac < 0 || avar.iinterp_fac > 1.0)
	{
		printf("ERROR: SPH aerodynamic interpolation factor must be between 0 < x < 1.\n");
		fault = 1;
	}

	if(svar.dropDragSweep)
	{
		if(svar.nacross.empty())
		{
			if(svar.diameters.empty())
			{
				printf("\n");
				fault = 1;
			}
			else
			{
				for(real const& dx : svar.diameters)
				{
					/* Droplet case */
					if(svar.diam/dx < 1.5)
					{	/*spacing is too close to full size to use standard adjustment*/
						svar.nacross.emplace_back(1);
					}
					else
					{
						svar.nacross.emplace_back(static_cast<size_t>(ceil(abs(svar.diam/dx))));
					}
				}
			}
		}

		fvar.pPress = 0;
	}

	/* Particle Tracking Settings */
	if(svar.using_ipt)
	{
		if(svar.eqOrder > 2 || svar.eqOrder < 1)
		{
			printf("Equation order not 1 or 2. Please choose between these.\n");
			exit(-1);
		}

		if(svar.max_x < svar.max_x_sph)
		{
			printf("WARNING: Maximum x coordinate for particle tracking is less than that for SPH.\n");
			printf("         No particle tracking will be performed.\n");
			svar.using_ipt = 0;
		}

		// svar.streakdir = svar.outfile;
		// size_t pos = svar.streakfile.find_last_of(".");
		// svar.streakfile.insert(pos,"_streak");
	}
	
    if(svar.offset_axis != 0)
    {
        if(SIMDIM != 2)
		{
			printf("WARNING: trying to use 3D code with a 2D settings file.\n");
			svar.offset_axis = 0;
		}
    }
    else if (svar.offset_axis > 3)
    {
        printf("ERROR: 2D offset axis option out of bounds\n");
        fault = 1;
    }

	#if SIMDIM == 2
	if(svar.offset_axis == 0)
	{
		printf("ERROR: Offset axis has not been defined.\n");
		fault = 1;
	}
	#endif
	
	if(fault)
	{
		printf("Check of the input settings finished with errors. Stopping\n");
		exit(-1);
	}

	Check_Output_Variables(svar);

  	Set_Values(svar,fvar,avar,vortex);

	Print_Settings(stdout,svar,fvar,avar);
	#ifdef DEBUG
	Print_Settings(dbout,svar,fvar,avar);
	#endif

} /*End of GetInput()*/


void Write_Tec_Headers(FILE* ff, FILE* fb, FILE* fg, SIM& svar, FLUID const& fvar, AERO const& avar, std::string const& prefix)
{
	if(svar.out_encoding == 1)
	{
		Init_Binary_PLT(svar,fvar,avar,prefix,"_fluid.szplt","Simulation Particles",svar.fluidFile);
	
		Init_Binary_PLT(svar,fvar,avar,prefix,"_boundary.szplt","Boundary Particles",svar.boundFile);
	}
	else
	{
		/* Write first timestep */
		string mainfile = prefix;
		mainfile.append("_fluid.dat");
		ff = fopen(mainfile.c_str(), "a");
		if(ff != NULL)
			Write_ASCII_header(ff, svar, "Simulation Particles");
		else
		{
			printf("Failed to open %s. Stopping.\n",mainfile.c_str());
			exit(-1);
		}

		/*If the boundary exists, write it.*/
		string bfile = prefix;
		bfile.append("_boundary.dat");
		fb = fopen(bfile.c_str(), "a");
		if(fb != NULL)
			Write_ASCII_header(fb, svar, "Boundary Particles");
		else
		{
			printf("Error opening %s file. Stopping\n", bfile.c_str());
			exit(-1);
		}	
	}
}

void Write_h5part_Headers(SIM& svar, FLUID const& fvar, AERO const& avar, std::string const& prefix)
{
	open_h5part_files(svar,fvar,avar,prefix,svar.ffile,svar.bfile);
}

void Write_Timestep(FILE* ff, FILE* fb, FILE* fg, SIM& svar, FLUID const& fvar, AERO const& avar, 
			LIMITS const& limits, SPHState const& pnp1)
{
	if(svar.write_tecio)
	{
		std::ostringstream oss;
		oss << svar.t;
		if(!svar.single_file)
		{
			string prefix = svar.output_prefix + "_time_" + oss.str();
			Write_Tec_Headers(ff,fb,fg,svar,fvar,avar,prefix);
		}
			
		if (svar.out_encoding == 1)
		{
			for(size_t bound = 0; bound < svar.nbound; bound++)
			{
				string title = "Boundary_" + std::to_string(bound) + "_" + limits[bound].name + 
					"_time_" + oss.str() + "s";
				Write_Binary_Timestep(svar,fvar.rho0,pnp1,limits[bound],title.c_str(),
					static_cast<int32_t>(bound+1),svar.boundFile);
			}

			for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
			{
				string title = "Fluid_" + std::to_string(block) + "_" + limits[block].name + 
					"_time_" + oss.str() + "s";
				Write_Binary_Timestep(svar,fvar.rho0,pnp1,limits[block],title.c_str(),
					static_cast<int32_t>(block+1),svar.fluidFile);
			}
		} 
		else
		{
			for(size_t bound = 0; bound < svar.nbound; bound++)
			{
				string title = "Boundary_" + std::to_string(bound) + "_" + limits[bound].name + 
					"_time_" + oss.str() + "s";
				Write_ASCII_Timestep(svar,fvar.rho0,pnp1,limits[bound].index.first,
						limits[bound].index.second,title.c_str(),bound+1,fb);
			}

			for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
			{
				string title = "Fluid_" + std::to_string(block) + "_" + limits[block].name + 
					"_time_" + oss.str() + "s";
			
				Write_ASCII_Timestep(svar,fvar.rho0,pnp1,limits[block].index.first,
						limits[block].index.second,title.c_str(),block+1,ff);
			}
		}

		if(!svar.single_file)
		{	// Close files after use
			if (ff != NULL)
				fclose(ff);
			if (fb != NULL)
				fclose(fb);
			if (fg != NULL)
				fclose(fg);

			if(svar.out_encoding == 1)
			{
				/*Combine the szplt files*/
				if(tecFileWriterClose(&svar.fluidFile))
					printf("Failed to close fluid file.\n");
					
				if(tecFileWriterClose(&svar.boundFile))
					printf("Failed to close boundary file.\n");
			}
			
		}
	}

	if(svar.write_h5part)
	{
		string prefix = svar.output_prefix;
		if(!svar.single_file)
		{
			std::ostringstream oss;
			oss << svar.t;
			prefix += "_" + std::to_string(svar.frame);
		}

		Write_h5part_Headers(svar,fvar,avar,prefix);
		write_h5part_data(svar.ffile, svar.bfile, svar, fvar, pnp1);
		// Files are closed in the write function
	}
	
}

void Remove_Old_Files(SIM const& svar)
{
	string dir;
	string prefix = svar.output_prefix;
	if(svar.output_prefix.find_last_of("/\\") != string::npos)
	{
		dir = svar.output_prefix.substr(0,svar.output_prefix.find_last_of("/\\"));
		// prefix = svar.output_prefix.substr(svar.output_prefix.find_last_of("/\\")+1);
	}
	else
	{
		dir = ".";
		prefix = dir + "/" + svar.output_prefix;
	}

	for(int ii = prefix.size()-1; ii >= 0; ii--)
	{
		if(prefix[ii] == '.')
			prefix.insert(ii,"\\");

		if(prefix[ii] == '/')
			prefix.insert(ii,"\\");
	}

	int fault = 0;
	if(svar.write_tecio)
	{	
		std::regex file_expr("(" + prefix + ")(.*)(\\.szplt\\.sz)(.*)$");
		// std::regex file_expr("(" + prefix + ")(.*)(\\.szplt)$");
		
		for(auto const& file : directory_iterator(dir))
		{
			// printf("%s: ",file.path().string().c_str());
			// cout << file.path().string() << endl;
			// Check it has the right prefix
			if(std::regex_match(file.path().string(),file_expr))
			{
				// printf("File matches the expression \n");
				// Remove the file
				try{std::filesystem::remove(file.path());}
				catch(const std::filesystem::__cxx11::filesystem_error& e)
				{
					printf("Could not remove tecplot file. Trying to remove other files.\n");
					fault = 1;
				}
			}
			// else
			// 	printf("File doesn't match the expression\n");
		}		
	}

	if(fault)
	{
		printf("Failed to remove some old simulation files, so cannot continue.\n");
		exit(-1);
	}

	if(svar.write_h5part)
	{
		std::regex fluid_expr("(" + prefix + ")(.*)(\\.h5part)");
		
		for(auto const& file : directory_iterator(dir))
		{
			// printf("%s: ",file.path().filename().string().c_str());
			// cout << file.path().string() << endl;
			// Check it has the right prefix
			
			if(std::regex_match(file.path().string(),fluid_expr))
			{
				// printf("File matches the expression \n");
				// Remove the file
				try{std::filesystem::remove(file.path());}
				catch(const std::filesystem::__cxx11::filesystem_error& e)
				{
					printf("Could not remove fluid h5part file. Trying to remove other files.\n");
					fault = 1;
				}
			}
			// else
			// 	printf("File doesn't match the expression\n");
		
		}	
	}

	if(fault)
	{
		printf("Failed to remove some files. Try closing paraview, as it doesn't correctly free the files.\n");
		exit(-1);
	}
}

void Check_Output_Variables(SIM& svar)
{
	// Set the initial output options
	svar.outvar = std::vector<uint>(29,0);
	for(size_t ii = 0; ii < 8; ++ii)
		svar.outvar[ii] = 1;

	size_t nVars = SIMDIM*3 + 5;
	svar.var_types = std::vector<int32_t>(nVars,2);
	size_t nOptVars = 0;

	if(!svar.output_names.empty())
	{
		// Split the names using underscore as a delimiter
		std::vector<std::string> vals;
		auto start = 0U;
		auto end = svar.output_names.find_first_of("_");
		while(end != std::string::npos)
		{
			vals.emplace_back(svar.output_names.substr(start, end - start));
			start = end + 1;
			end = svar.output_names.find_first_of("_",start);
		}
		vals.emplace_back(svar.output_names.substr(start)); // catch the final value

		// Now check the variables for which ones want to be output
		for(std::string const& var : vals )
		{
			if(var == "dens")
			{
				svar.outvar[8] = 1;
				nOptVars+=1;
			}
			else if(var == "densVar")
			{
				svar.outvar[9] = 1;
				nOptVars += 1;
			}
			else if(var == "vmag")
			{
				svar.outvar[10] = 1;
				nOptVars+=1;
			}
			else if(var == "surf")
			{
				svar.outvar[11] = 1;
				nOptVars+=1;
			}
			else if(var == "surfZ")
			{
				svar.outvar[12] = 1;
				nOptVars+=1;
			}
			else if(var == "aero-mag")
			{
				svar.outvar[13] = 1;
				nOptVars+=1;
			}
			else if(var == "aero-vec")
			{
				svar.outvar[14] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "curv")
			{
				svar.outvar[15] = 1;
				nOptVars+=1;
			}
			else if(var == "occl")
			{
				svar.outvar[16] = 1;
				nOptVars+=1;
			}
			else if(var == "cellP")
			{
				svar.outvar[17] = 1;
				nOptVars+=1;
			}
			else if(var == "cellRho")
			{
				svar.outvar[18] = 1;
				nOptVars+=1;
			}
			else if(var == "cellV-mag")
			{
				svar.outvar[19] = 1;
				nOptVars+=1;
			}
			else if(var == "cellV-vec")
			{
				svar.outvar[20] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "dsphG")
			{
				svar.outvar[21] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "lam")
			{
				svar.outvar[22] = 1;
				nOptVars+=1;
			}
			else if(var == "lam-nb")
			{
				svar.outvar[23] = 1;
				nOptVars+=1;
			}
			else if(var == "colour")
			{
				svar.outvar[24] = 1;
				nOptVars+=1;
			}
			else if(var == "colour-G")
			{
				svar.outvar[25] = 1;
				nOptVars+=1;
			}
			else if(var == "norm")
			{
				svar.outvar[26] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "shiftV-mag")
			{
				svar.outvar[27] = 1;
				nOptVars+=1;
			}
			else if(var == "shiftV-vec")
			{
				svar.outvar[28] = 1;
				nOptVars+=SIMDIM;
			}
			else
			{
				std::printf("Unrecognised output variable \"%s\" defined\n", var.c_str());
			}
		}
	}
	
	svar.var_types.resize(nVars+nOptVars,0);

	// Define the variable names
	string axis1 = "X"; string axis1s = "x";
	string axis2 = "Y"; string axis2s = "y";
	#if SIMDIM == 2
	if(svar.offset_axis != 0)
	{
		if(svar.offset_axis == 1)
		{
			axis1 = "Y"; axis1s = "y";
			axis2 = "Z"; axis2s = "z";
		}
		else if(svar.offset_axis == 2)
		{
			axis1 = "X"; axis1s = "x";
			axis2 = "Z"; axis2s = "z";
		}
		else if(svar.offset_axis == 3)
		{
			axis1 = "X"; axis1s = "x";
			axis2 = "Y"; axis2s = "y";
		}
	}
	#endif

	size_t varCount = 0;
	svar.var_names.clear();
	if(svar.outvar[0])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += axis1 + "," + axis2;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Z";
		#endif
	}

	if(svar.outvar[1])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",A" + axis1s + ",A" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Az";
		#endif
	}

	if(svar.outvar[2])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",V" + axis1s + ",V" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Vz";
		#endif
	}

	if(svar.outvar[3])
	{
		svar.var_names += ",Pressure";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[4])
	{
		svar.var_names += ",dRho";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[5])
	{
		svar.var_names += ",partID";
		svar.var_types[varCount] = 3; varCount++;
	}

	if(svar.outvar[6])
	{
		svar.var_names += ",cellID";
		svar.var_types[varCount] = 3; varCount++;
	}

	if(svar.outvar[7])
	{
		svar.var_names += ",bound";
		svar.var_types[varCount] = 5; varCount++;
	}
	
	// Add optionals
	if(svar.outvar[8])
	{
		svar.var_names += ",Dens";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[9])
	{
		svar.var_names += ",DensVar";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[10])
	{
		svar.var_names += ",Vmag";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[11])
	{
		svar.var_names += ",Surf";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[12])
	{
		svar.var_names += ",Surfzone";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[13])
	{
		svar.var_names += ",Aero-mag";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[14])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",Af" + axis1s + ",Af" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Afz";
		#endif
	}

	if(svar.outvar[15])
	{
		svar.var_names += ",curvature";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[16])
	{
		svar.var_names += ",occl";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[17])
	{
		svar.var_names += ",cellP";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[18])
	{
		svar.var_names += ",cellRho";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[19])
	{
		svar.var_names += ",cellV-mag";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[20])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",cellV" + axis1s + ",cellV" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",cellVz";
		#endif
	}

	if(svar.outvar[21])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",dsphG" + axis1s + ",dsphG" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",dsphGz";
		#endif
	}
	
	if(svar.outvar[22])
	{
		svar.var_names += ",lam";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[23])
	{
		svar.var_names += ",lam-nb";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[24])
	{
		svar.var_names += ",colour";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[25])
	{
		svar.var_names += ",colourG";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[26])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",surf-norm" + axis1s + ",surf-norm" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",surf-normz";
		#endif
	}
	
	if(svar.outvar[27])
	{
		svar.var_names += ",shiftV-mag";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[28])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",shiftV" + axis1s + ",shiftV" + axis2s;
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",shiftVz";
		#endif
	}
}

