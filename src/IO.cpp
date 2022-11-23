/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "IO.h"
#include "AsciiIO.h"
#include "BinaryIO.h"
#include "Geometry.h"
#include "IOFunctions.h"
#include "Kernel.h"
#include <ctime>
#include <Eigen/LU>
// #include "Restart.h"
// #include "TauIO.h"

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void Set_Values(SIM& svar, FLUID& fvar, AERO& avar)
{
	/*Universal parameters based on input values*/
	svar.Angle *= M_PI/180.0;
	svar.Rotate = GetRotationMat(svar.Angle);
	svar.Transp = svar.Rotate.transpose();
	
	svar.sim_start *= svar.scale;
	svar.bound_start *= svar.scale;
	svar.max_x *= svar.scale;
	svar.max_x_sph *= svar.scale;

	if(avar.vJetMag != -1)
	{
		avar.vStart = StateVecD::Zero();
		avar.vStart[1] = avar.vJetMag;
		avar.vStart = svar.Rotate*avar.vStart;
	}

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
	// 		cout << "Warning: particle spacing if of the same order of magintude of the jet diameter." << endl;
	// 		cout << "Consider a different size for accuracy." << endl;
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
	// 	cout << "Droplet Diameter: " << svar.diam << endl;		
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
		
	avar.qInf = 0.5*avar.MRef*avar.MRef*avar.gamma*avar.pRef;

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

	// if(svar.Scase == BOX)
	// {
	// 	#if SIMDIM == 2
	// 	svar.simPts = svar.xyPART[0]*svar.xyPART[1];
	// 	#else
	// 	svar.simPts = svar.xyPART[0]*svar.xyPART[1]*svar.xyPART[2]; /*total sim particles*/
	// 	#endif
	// }

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
		if (svar.Asource == 3)
		{
			svar.vortex.Init(svar.vlm_file);
			svar.vortex.GetGamma(avar.vInf);
			svar.vortex.write_VLM_Panels(svar.output_prefix);
		}
	#endif

	avar.GetYcoef(fvar, /*fvar.H*/ svar.Pstep);

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

void Print_Settings(char** argv, SIM const& svar, FLUID const& fvar, AERO const& avar)
{
	#ifdef DEBUG
		dbout << "Tait Gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
		dbout << "Pipe rho: " << fvar.rhoJ << endl;
		dbout << "Number of fluid particles along diameter: " << 2*svar.nrad+1 << endl;
		dbout << "Pipe step (dx): " << svar.dx << endl;
		dbout << "Particle mass: " << fvar.simM << endl;
		dbout << "Freestream initial spacing: " << svar.Pstep << endl;
		dbout << "Support Radius: " << fvar.H << endl;
		dbout << "Gas Dynamic pressure: " << avar.qInf << endl << endl;
	#endif

	cout << "****** SIMULATION SETTINGS *******" << endl;
	#pragma omp parallel
	{
		#pragma omp single
		cout << "Number of threads: " << omp_get_num_threads() << endl;
	}
	
	cout << "Frame time interval: " << svar.framet << endl;
	cout << "Number of frames: " << svar.Nframe << endl;
	cout << "Output variables: " << svar.output_names << endl;
	// cout << "Tait gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
	fprintf(stdout, "Tait gamma: %.1f  Tait B: %g\n",fvar.gam,fvar.B);
	// cout << "Newmark-Beta parameters: " << svar.beta << ", " << svar.gamma << endl;
	fprintf(stdout, "Newmark-Beta parameters: %.2f, %.2f\n",svar.beta,svar.gamma);
	cout << "Boundary case: " << svar.Bcase << endl;
	cout << "Aerodynamic source: " << svar.Asource << endl;
	cout << "Aerodynamic case: " << avar.acase << endl;
	cout << endl;

	cout << "****** PIPE SETTINGS *******" << endl;
	cout << "Pipe pressure: " << fvar.pPress << endl;
	cout << "Pipe density: " << fvar.rhoJ << endl;
	cout << "Pipe step (dx): " << svar.dx << endl;
	cout << "Pipe diameter: " << svar.jet_diam << endl;
	cout << "Number of fluid particles along diameter: " << 2*svar.nrad+1 << endl;
	cout << "Pipe start position: " << svar.sim_start(0) << "  " << svar.sim_start(1);
	#if SIMDIM == 3
		cout << "  " << svar.sim_start(2);
	#endif
	cout << endl;
	cout << "Pipe start rotation: " << svar.Angle(0) << "  " << svar.Angle(1);
	#if SIMDIM == 3
		cout << "  " << svar.Angle(2);
	#endif
	cout << endl;
	cout << "Jet velocity: " << avar.vJetMag << endl;

	cout << endl;		
	cout << "****** RESTING FUEL SETTINGS *******" << endl;
	cout << "Resting density: " << fvar.rho0 << endl;
	cout << "Maxmimum density: " << fvar.rhoMax << endl;
	cout << "Minimum density: " << fvar.rhoMin << endl;
	cout << "Particle spacing: " << svar.Pstep << endl;
	cout << "Support radius: " << fvar.H << endl;
	cout << "Particle mass: " << fvar.simM << endl;
	cout << "Liquid viscosity: " << fvar.mu << endl;
	cout << "Speed of sound: " << fvar.Cs << endl;
	cout << endl;

	cout << "****** FREESTREAM SETTINGS ******" << endl;
	if(svar.Asource != 1 && svar.Asource != 2)
	{
		cout << "Gas Velocity: " << avar.vInf(0) << "  " << avar.vInf(1);
		#if SIMDIM == 3
			cout << "  " << avar.vInf(2);
		#endif
		cout << endl;
	}

	cout << "Gas density: " << avar.rhog << endl;
	cout << "Gas viscosity: " << avar.mug << endl;
	cout << "Aerodynamic length: " << avar.L << endl;
	cout << endl;

	if(svar.Asource == 1 || svar.Asource == 2)
	{
		cout << "Reference velocity: " << avar.vRef << endl;
		cout << "Reference pressure: " << avar.pRef << endl;
		cout << "Reference temperature: " << avar.T << endl;
		cout << "Gas dynamic pressure: " << avar.qInf << endl << endl;
	}
	
	cout << "******** FILE SETTINGS ********" << endl;	
	// cout << "Working directory: " << cCurrentPath << endl;
	cout << "Input file: " << argv[1] << endl;
	if(svar.Asource == 1 || svar.Asource == 2)
	{
		if(svar.CDForFOAM == 0)
		{
			cout << "Mesh filename: " << svar.taumesh << endl;
			cout << "Solution filename: " << svar.tausol << endl;
		}
		else
		{
			cout << "OpenFOAM root directory: " << svar.foamdir << endl;
			cout << "OpenFOAM solution time: " << svar.foamsol << endl;
		}
	}

	cout << "Output prefix: " << svar.output_prefix << endl << endl;

}


void GetInput(int argc, char **argv, SIM& svar, FLUID& fvar, AERO& avar)
{
	if (argc > 2) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only a maximum of one input arguments accepted,\n";
		cout << "1: Input file directory\n";
		cout << "Other inputs will be ignored." << endl << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tERROR: No inputs provided. Stopping... \n";
    	exit(-1);    	
    }

	/*Get parameters if it has been provided*/
	// cout << argv[1] << endl;
	string file = argv[1];
	    	
	svar.infile = file;

	ifstream fin(argv[1]);
	if(!fin.is_open())
	{
		cout << "ERROR: Failed to open parameter file: " << argv[1] << endl;
		exit(-1);
	}
	else
	{
		cout << "para file open, reading contents..." << endl;
	}

    string line;
    while (getline(fin,line))
    {
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
        Get_Number(line, "2D offset vector (0 / x=1,y=2,z=3)",svar.offset_axis);
        Get_String(line, "OpenFOAM folder", svar.foamdir);
        Get_String(line, "OpenFOAM solution folder", svar.foamsol);
        Get_Number(line, "OpenFOAM binary (0/1)", svar.isBinary);
        Get_Number(line, "Label size (32/64)", svar.labelSize);
        Get_Number(line, "Scalar size (32/64)", svar.scalarSize);
        Get_Number(line, "OpenFOAM buoyant (0/1)", svar.buoyantSim);
		#if SIMDIM == 3
		Get_String(line, "VLM definition filename", svar.vlm_file); 
		#endif 

        /* File outputs */
		Get_Number(line, "Single file for output (0/1)",svar.single_file);
		Get_String(line, "Output files prefix", svar.output_prefix);
        Get_Number(line, "SPH frame time interval", svar.framet);
        Get_Number(line, "SPH frame count", svar.Nframe);
        Get_Number(line, "SPH output encoding (0=ascii/1=binary)", svar.out_encoding);
        Get_Number(line, "SPH ghost output (0/1)", svar.gout);
        Get_String(line, "Variable list", svar.output_names);
        
        // Get_String(line, "Particle surface impact filename", svar.surfacefile);

        /* Fluid data */
        Get_Number(line, "Reference density", avar.rhog);
        Get_Number(line, "Reference dispersed density", fvar.rho0);
        Get_Number(line, "Sutherland reference viscosity", avar.mug);
        Get_Number(line, "Reference dispersed viscosity", fvar.mu);
        Get_Number(line, "Reference surface tension", fvar.sig);
		Get_Number(line, "SPH surface tension contact angle", fvar.contangb);
        Get_Number(line, "Init hydrostatic pressure (0/1)", svar.init_hydro_pressure);
		Get_Number(line, "Hydrostatic height", svar.hydro_height);

		/* Aerodynamic data */
        Get_Number(line, "Reference velocity", avar.vRef);
        Get_Number(line, "Reference pressure", avar.pRef);
        Get_Number(line, "Reference Mach number", avar.MRef);
        Get_Number(line, "Reference density", avar.qInf);
        Get_Number(line, "Reference temperature", avar.T);
		Get_Number(line, "Gas constant gamma", avar.gamma);

		/* Simulation settings */
		Get_Number(line, "SPH boundary solver (0=pressure/1=ghost)", svar.bound_solver);
		Get_Number(line, "SPH maximum timestep", svar.dt_max);
		Get_Number(line, "SPH minimum timestep", svar.dt_min);
        Get_Number(line, "SPH CFL condition", svar.cfl);
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
        Get_String(line, "SPH aerodynamic case", avar.aero_case);
        Get_Number(line, "SPH SP diameter definition (0=dx/1=h)", avar.use_dx);
        Get_Number(line, "SPH use TAB deformation (0/1)", avar.useDef);
        Get_Number(line, "SPH use ghost particles (0/1)", svar.ghost);
        Get_Vector(line, "SPH start coordinates", svar.sim_start);
		Get_Vector(line, "SPH boundary start coordinates", svar.bound_start);
		Get_Number(line, "SPH maximum particle count",svar.finPts);
        Get_Number(line, "SPH aerodynamic cutoff value", avar.cutoff);

        /* Starting area conditions */
        Get_String(line, "SPH start geometry type", svar.start_type);
		Get_String(line, "SPH boundary geometry type", svar.bound_type);
        Get_Vector(line, "SPH box resolution", svar.xyPART);
        Get_Vector(line, "SPH box lengths", svar.sim_box);
		Get_Vector(line, "SPH boundary box lengths", svar.bound_box);
        Get_Vector(line, "SPH rotation angles (deg)", svar.Angle);
        Get_Number(line, "SPH jet diameter", svar.jet_diam);
        Get_Number(line, "SPH jet depth", svar.jet_depth);
        Get_Number(line, "SPH jet velocity", avar.vJetMag);
        Get_Number(line, "SPH sphere diameter", svar.diam);
		Get_Vector(line, "SPH starting velocity", avar.vStart);
        Get_Vector(line, "SPH freestream velocity", avar.vInf);
        Get_Number(line, "SPH restart fit tolerance", svar.restart_tol);

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

		Get_Number(line, "Do speed test (0/1)", svar.speedTest);
		Get_Number(line, "Speed test run count", svar.nRuns);
    }

    fin.close();

	/* Need to check if inputs are correct */
	
	if(svar.taumesh.empty())
    {
		
		if(svar.foamdir.empty())
		{
			#if SIMDIM == 3
			if(svar.vlm_file.empty())
				svar.Asource = 0;
			else
			{
				svar.Asource = 3;

				if(svar.vlm_file == "(thisfile)")
					svar.vlm_file = svar.infile;
			}

			#else
			svar.Asource = 0;
			#endif
		}
		else
		{
			if(svar.foamsol.empty())
			{
				cout << "OpenFOAM solution directory not defined." << endl;
				exit(-1);
			}

			svar.CDForFOAM = 1;
			svar.Asource = 1;			
		}
    }
	else
	{
		svar.CDForFOAM = 0;
		svar.Asource = 1;
		if(svar.taubmap.empty())
		{
			cout << "Input TAU bmap file not defined." << endl;
			exit(-1);
		}
		else if(svar.taubmap == "(thisfile)")
		{
			/* The para file is the boundary definition file, so set it so */
			svar.taubmap = svar.infile;
		}

		if(svar.tausol.empty())
		{
			cout << "Input TAU solution file not defined." << endl;
			exit(-1);
		}

	}
	
	/* Check starting geometry conditions */
	if (svar.sim_start[0] == 999999 || svar.sim_start[1] == 999999
		#if SIMDIM == 3
			|| svar.sim_start[2] == 999999
		#endif
	)
	{
		cout << "Some or all of the SPH start coordinates have not been defined." << endl;
		exit(-1);
	}

	// if(svar.start_type == "box")
	// {
	// 	svar.Scase = BOX;

	// 	if(svar.xyPART[0] < 0 || svar.xyPART[1] < 0
	// 		#if SIMDIM == 3
	// 			|| svar.xyPART[2] < 0
	// 		#endif
	// 		)
	// 	{
	// 		if(svar.Pstep < 0)
	// 		{
	// 			cout << "SPH box resolution or initial spacing has not been defined." << endl;
	// 			exit(-1);
	// 		}	
	// 	}
	// 	else
	// 	{
	// 		real dx = svar.sim_box[0]/svar.xyPART[0];
	// 		real dy = svar.sim_box[1]/svar.xyPART[1];
	// 		svar.Pstep = std::min(dx,dy);
	// 	}
		
	// }
	// else if(svar.start_type == "cylinder")
	// {
	// 	svar.Scase = CYLINDER;

	// 	if(svar.jet_diam < 0)
	// 	{
	// 		cout << "SPH cylinder diameter has not been defined." << endl;
	// 		exit(-1);
	// 	}

	// 	if(svar.jet_depth < 0)
	// 	{
	// 		cout << "SPH cylinder depth has not been defined." << endl;
	// 		exit(-1);
	// 	}

	// }	
	// else if(svar.start_type == "sphere")
	// {
	// 	svar.Scase = SPHERE;
	// 	if(svar.diam < 0)
	// 	{
	// 		cout << "SPH sphere diameter has not been defined." << endl;
	// 		exit(-1);
	// 	}
	// }
	// else if (svar.start_type == "jet")
	// {
	// 	svar.Scase = JET;

	// 	if(svar.jet_diam < 0)
	// 	{
	// 		cout << "SPH cylinder diameter has not been defined." << endl;
	// 		exit(-1);
	// 	}

	// 	if(svar.jet_depth < 0)
	// 	{
	// 		cout << "SPH cylinder depth has not been defined." << endl;
	// 		exit(-1);
	// 	}
	// }
	// else if (svar.start_type == "taper_jet")
	// {
	// 	svar.Scase = CONVJET;

	// 	if(svar.jet_diam < 0)
	// 	{
	// 		cout << "SPH cylinder diameter has not been defined." << endl;
	// 		exit(-1);
	// 	}

	// 	if(svar.jet_depth < 0)
	// 	{
	// 		cout << "SPH cylinder depth has not been defined." << endl;
	// 		exit(-1);
	// 	}
	// }
	// else
	// {
	// 	cout << "SPH starting geometry is not defined correctly if at all." << endl;
	// 	exit(-1);
	// }


	// if(svar.bound_type == "(none)")
	// {
	// 	svar.Bcase = NONE;
	// }
	// else if (svar.bound_type == "box")
	// {
	// 	svar.Bcase = BOX;

	// 	if (svar.bound_start[0] == 999999 || svar.bound_start[1] == 999999
	// 		#if SIMDIM == 3
	// 		|| svar.bound_start[2] == 999999
	// 		#endif
	// 	)
	// 	{
	// 		cout << "Some or all of the SPH boundary start coordinates have not been defined." << endl;
	// 		exit(-1);
	// 	}

	// 	if(svar.bound_box[0] < 0 || svar.bound_box[1] < 0
	// 		#if SIMDIM == 3
	// 			|| svar.bound_box[2] < 0
	// 		#endif
	// 		)
	// 	{
	// 		cout << "Some or all of the SPH boundary box lengths have not been defined." << endl;
	// 		exit(-1);
	// 	}
	// }
	// else if (svar.bound_type == "jet")
	// {
	// 	svar.Bcase = JET;
	// 	if(svar.Scase != JET)
	// 	{
	// 		cout << "Boundary geometry is not compatible with starting geometry." << endl;
	// 		cout << "Starting geometry must be a cylinder to use a jet boundary." << endl;
	// 		exit(-1);
	// 	}

	// 	if(avar.vJetMag < 0)
	// 	{
	// 		cout << "Jet velocity has not been defined" << endl;
	// 		exit(-1);
	// 	}

	// 	svar.bound_start = svar.sim_start;
	// }
	// else if (svar.bound_type == "taper_jet")
	// {
	// 	svar.Bcase = CONVJET;
	// 	if(svar.Scase != CONVJET)
	// 	{
	// 		cout << "Boundary geometry is not compatible with starting geometry." << endl;
	// 		cout << "Starting geometry must be a cylinder to use a jet boundary." << endl;
	// 		exit(-1);
	// 	}

	// 	if(avar.vJetMag < 0)
	// 	{
	// 		cout << "Jet velocity has not been defined" << endl;
	// 		exit(-1);
	// 	}

	// 	svar.bound_start = svar.sim_start;
	// }
	// else
	// {
	// 	cout << "SPH boundary geometry is not defined correctly if at all." << endl;
	// 	exit(-1);
	// }

	// if(svar.Bcase != NONE)
	// {
	// 	if(svar.bound_solver > 1)
	// 	{
	// 		cout << "Boundary solver not defined correctly." << endl;
	// 		exit(-1);
	// 	}
	// }

	/* Check for restart data now boundary case is known */	
	if(!svar.restart_prefix.empty())
	{
		svar.restart = 1;
		// Check_If_Restart_Possible(svar);
		// svar.output_prefix = svar.restart_prefix;
	}


	if(svar.Pstep < 0)
	{
		cout << "SPH initial spacing has not been defined." << endl;
		exit(-1);
	}

	if(svar.Nframe < 0)
	{
		cout << "Number of frames to output not defined." << endl;
		exit(-1);
	}

	if (svar.framet < 0)
	{
		cout << "Frame time interval has not been defined." << endl;
		exit(-1);
	}

	if(svar.init_hydro_pressure)
	{
		if(svar.hydro_height < 0)
		{
			cout << "Hydrostatic height has not been defined." << endl;
			exit(-1);
		}
	}

	/* Aerodynamic settings */
	if(avar.aero_case == "Gissler")
	{
		avar.acase = 1;
	}
	else if (avar.aero_case == "Induced_pressure")
	{
		avar.acase = 2;
	}
	else if (avar.aero_case == "Skin_friction")
	{
		avar.acase = 3;
	}
	else if(avar.aero_case == "(none)")
	{
		avar.acase = 0;
	}
	else
	{
		cout << "Aerodynamic coupling model is not defined or correct." << endl;
		exit(-1);
	}

	if(svar.dt_min > 0.0)
	{
		svar.dt = svar.dt_min;
	}

	if(svar.dropDragSweep)
	{
		if(svar.nacross.empty())
		{
			if(svar.diameters.empty())
			{
				cout << "" << endl;
				exit(-1);
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
		// if(svar.partout == 1)
		// {
		//     // if(svar.outdir.empty())
		//     // {
		//     //     cout << "Output particle directory not defined." << endl;
		//     //     exit(-1);
		//     // }
		// 	string partf = svar.output_prefix + "_parts.dat";

		// 	svar.partfile.open(partf,std::ios::out);
		// 	if(!svar.cellfile)
		// 	{
		// 		cout << "Failed to open the particle tracking scatter file." << endl;
		// 		exit(-1);
		// 	}
		// }

		// if(svar.streakout == 1)
		// {
		//     // if(svar.streakdir.empty())
		//     // {
		//     //     cout << "Output particle streaks directory not defined." << endl;
		//     //     exit(-1);
		//     // }
		// 	string streakf = svar.output_prefix + "_streaks.dat";
		// 	svar.streakfile.open(streakf,std::ios::out);
		// 	if(!svar.cellfile)
		// 	{
		// 		cout << "Failed to open the particle tracking streaks file." << endl;
		// 		exit(-1);
		// 	}
		// }

		// if(svar.cellsout == 1)
		// {
		//     // if(svar.celldir.empty())
		//     // {
		//     //     cout << "Output cell intersections directory not defined." << endl;
		//     //     exit(-1);
		//     // }
		// 	string cellf = svar.output_prefix + "_cells.dat";
		// 	svar.cellfile.open(cellf,std::ios::out);
		// 	if(!svar.cellfile)
		// 	{
		// 		cout << "Failed to open the particle tracking cell file." << endl;
		// 		exit(-1);
		// 	}
		// }

		if(svar.eqOrder > 2 || svar.eqOrder < 1)
		{
			cout << "Equation order not 1 or 2. Please choose between these." << endl;
			exit(-1);
		}

		if(svar.max_x < svar.sim_start[0])
		{
			cout << "ERROR: Maximum x coordinate for particle tracking or SPH is less than the starting position." << endl;
			exit(-1);
		}

		if(svar.max_x < svar.max_x_sph)
		{
			cout << "WARNING: Maximum x coordinate for particle tracking is less than that for SPH." << endl;
			cout << "         No particle tracking will be performed." << endl;
			svar.using_ipt = 0;
		}

		if(svar.max_x_sph < svar.sim_start[0])
		{
			cout << "WARNING: Maximum x coordinate for SPH particles is less than the starting position." << endl;
			cout << "         Particles will immediately transition to tracking." << endl;
		}

		// svar.streakdir = svar.outfile;
		// size_t pos = svar.streakfile.find_last_of(".");
		// svar.streakfile.insert(pos,"_streak");
	}
	
    if(svar.offset_axis != 0)
    {
        if(SIMDIM != 2)
		{
			cout << "WARNING: trying to use 3D code with a 2D settings file." << endl;
			svar.offset_axis = 0;
		}
    }
    else if (svar.offset_axis > 3)
    {
        cout << "2D offset axis option out of bounds" << endl;
        exit(-1);
    }

	#if SIMDIM == 2
	if(svar.offset_axis == 0)
	{
		cout << "ERROR: Offset axis has not been defined." << endl;
		exit(-1);
	}
	#endif
	
	Check_Output_Variables(svar);

  	Set_Values(svar,fvar,avar);

	Print_Settings(argv,svar,fvar,avar);

} /*End of GetInput()*/


void Write_Headers(FILE* f1, FILE* fb, FILE* fg, SIM& svar)
{
	if(svar.out_encoding == 1)
	{
		Init_Binary_PLT(svar,"_fluid.szplt","Simulation Particles",svar.fluidFile);
	
		Init_Binary_PLT(svar,"_boundary.szplt","Boundary Particles",svar.boundFile);
	}
	else
	{
		/* Write first timestep */
		string mainfile = svar.output_prefix;
		mainfile.append("_fluid.dat");
		f1 = fopen(mainfile.c_str(), "a");
		if(f1 != NULL)
			Write_ASCII_header(f1, svar, "Simulation Particles");
		else
		{
			cerr << "Failed to open " << mainfile << ". Stopping." << endl;
			exit(-1);
		}

		/*If the boundary exists, write it.*/
		string bfile = svar.output_prefix;
		bfile.append("_boundary.dat");
		fb = fopen(bfile.c_str(), "a");
		if(fb != NULL)
			Write_ASCII_header(fb, svar, "Boundary Particles");
		else
		{
			cerr << "Error opening " << bfile << " file. Stopping" << endl;
			exit(-1);
		}
		
	}
}

void Write_Timestep(FILE* f1, FILE* fb, FILE* fg, SIM& svar, LIMITS const& limits,
				 SPHState const& pnp1/* , DELTAP const& dp *//* , SPHState const& airP */)
{
	if (svar.out_encoding == 1)
	{
		for(size_t bound = 0; bound < svar.nbound; bound++)
		{
			string title = "Boundary_" + std::to_string(bound) + "_" + limits[bound].name;
			Write_Binary_Timestep(svar,pnp1,limits[bound],title.c_str(),
                static_cast<int32_t>(bound+1),svar.boundFile);
		}

		for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
		{
			string title = "Fluid_" + std::to_string(block) + "_" + limits[block].name ;
			Write_Binary_Timestep(svar,pnp1,limits[block],title.c_str(),
                static_cast<int32_t>(block+1),svar.fluidFile);
		}
	} 
	else
	{
		for(size_t bound = 0; bound < svar.nbound; bound++)
		{
			string title = "Boundary_" + std::to_string(bound) + "_" + limits[bound].name;
			Write_ASCII_Timestep(svar,pnp1,limits[bound].index.first,
					limits[bound].index.second,title.c_str(),bound+1,fb);
		}

		for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
		{
			string title = "Fluid_" + std::to_string(block) + "_" + limits[block].name ;
		
			Write_ASCII_Timestep(svar,pnp1,limits[block].index.first,
					limits[block].index.second,title.c_str(),block+1,f1);
		}
	}
}

void Append_Restart_Prefix(SIM const& svar)
{
	if(svar.restart)
	{
		if(svar.output_prefix == svar.restart_prefix)
			return;
	}

	/* Append a restart prefix to the end of the para file */
	ofstream para(svar.infile, std::ios_base::app | std::ios_base::out);
	if(!para.is_open())
	{
		cout << "Failed to reopen para file" << endl;
		exit(-1);
	}
	std::time_t now = std::time(0);
	char* time = ctime(&now);
	para << "\n";
	para << "    solver at " << time;
	para << "                                               " << 
			"SPH restart prefix: " << svar.output_prefix << endl;

	para.close();
	
}

void Check_Output_Variables(SIM& svar)
{
	// Set the initial output options
	svar.outvar = std::vector<uint>(28,0);
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
			else if(var == "vmag")
			{
				svar.outvar[9] = 1;
				nOptVars+=1;
			}
			else if(var == "surf")
			{
				svar.outvar[10] = 1;
				nOptVars+=1;
			}
			else if(var == "surfZ")
			{
				svar.outvar[11] = 1;
				nOptVars+=1;
			}
			else if(var == "aero-mag")
			{
				svar.outvar[12] = 1;
				nOptVars+=1;
			}
			else if(var == "aero-vec")
			{
				svar.outvar[13] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "curv")
			{
				svar.outvar[14] = 1;
				nOptVars+=1;
			}
			else if(var == "occl")
			{
				svar.outvar[15] = 1;
				nOptVars+=1;
			}
			else if(var == "cellP")
			{
				svar.outvar[16] = 1;
				nOptVars+=1;
			}
			else if(var == "cellRho")
			{
				svar.outvar[17] = 1;
				nOptVars+=1;
			}
			else if(var == "cellV-mag")
			{
				svar.outvar[18] = 1;
				nOptVars+=1;
			}
			else if(var == "cellV-vec")
			{
				svar.outvar[19] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "dsphG")
			{
				svar.outvar[20] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "lam")
			{
				svar.outvar[21] = 1;
				nOptVars+=1;
			}
			else if(var == "lam-nb")
			{
				svar.outvar[22] = 1;
				nOptVars+=1;
			}
			else if(var == "colour")
			{
				svar.outvar[23] = 1;
				nOptVars+=1;
			}
			else if(var == "colour-G")
			{
				svar.outvar[24] = 1;
				nOptVars+=1;
			}
			else if(var == "norm")
			{
				svar.outvar[25] = 1;
				nOptVars+=SIMDIM;
			}
			else if(var == "shiftV-mag")
			{
				svar.outvar[26] = 1;
				nOptVars+=1;
			}
			else if(var == "shiftV-vec")
			{
				svar.outvar[27] = 1;
				nOptVars+=SIMDIM;
			}
			else
			{
				std::cout << "Unrecognised output variable \"" << var << "\" defined" << std::endl;
			}
		}
	}
	
	svar.var_types.resize(nVars+nOptVars,0);

	// Define the variable names
	size_t varCount = 0;
	svar.var_names.clear();
	if(svar.outvar[0])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += "X,Y";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Z";
		#endif
	}

	if(svar.outvar[1])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",Ax,Ay";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Az";
		#endif
	}

	if(svar.outvar[2])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",Vx,Vy";
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
		svar.var_names += ",Vmag";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[10])
	{
		svar.var_names += ",Surf";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[11])
	{
		svar.var_names += ",Surfzone";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[12])
	{
		svar.var_names += ",Aero-mag";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[13])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",Afx,Afy";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",Afz";
		#endif
	}

	if(svar.outvar[14])
	{
		svar.var_names += ",curvature";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[15])
	{
		svar.var_names += ",occl";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[16])
	{
		svar.var_names += ",cellP";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[17])
	{
		svar.var_names += ",cellRho";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[18])
	{
		svar.var_names += ",cellV-mag";
		svar.var_types[varCount] = realType; varCount++;
	}

	if(svar.outvar[19])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",cellVx,cellVy";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",cellVz";
		#endif
	}

	if(svar.outvar[20])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",dsphGx,dsphGy";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",dsphGz";
		#endif
	}
	
	if(svar.outvar[21])
	{
		svar.var_names += ",lam";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[22])
	{
		svar.var_names += ",lam-nb";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[23])
	{
		svar.var_names += ",colour";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[24])
	{
		svar.var_names += ",colourG";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[25])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",surf-normx,surf-normy";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",surf-normz";
		#endif
	}
	
	if(svar.outvar[26])
	{
		svar.var_names += ",shiftV-mag";
		svar.var_types[varCount] = realType; varCount++;
	}
	
	if(svar.outvar[27])
	{
		svar.var_types[varCount] = realType; varCount++;
		svar.var_types[varCount] = realType; varCount++;
		svar.var_names += ",shiftVx,shiftVy";
		#if SIMDIM == 3
			svar.var_types[varCount] = realType; varCount++;
			svar.var_names += ",shiftVz";
		#endif
	}

}

void Restart_Simulation(SIM& svar, FLUID const& fvar, AERO const& avar, 
		MESH const& cells, SPHState& pn, SPHState& pnp1, LIMITS& limits)
{
	// Now get the data from the files. Start with the boundary
	if(svar.out_encoding == 1)
	{
		Restart_Binary(svar,fvar,pn,limits);
	}
	else
	{
		// TODO: ASCII Restart.
		// Particle numbers can be found from frame file.
		// Find EOF, then walk back from there how many particles.
		ASCII_Restart(svar,fvar,pn);
		cout << "ASCII file restart not yet implemented." << endl;
		exit(-1);
	}

	// Go through the particles giving them the properties of the cell
	size_t pID = 0;
	vector<vector<size_t>> buffer(limits.size());
	for(size_t block = 0; block < limits.size(); ++block)
	{
		#pragma omp parallel for reduction(max:pID)
		for(size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
		{
			pn[ii].xi *= svar.scale;
			
			/* Set density based on pressure. More information this way since density has much less variation*/
			pn[ii].rho = std::max(fvar.rhoMin,std::min(fvar.rhoMax, 
				density_equation(pn[ii].p, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP)));

			if(pn[ii].b == BACK)
			{
				#pragma omp critical
				limits[block].back.emplace_back(ii);
			}
			else if(pn[ii].b == BUFFER)
			{
				#pragma omp critical
				buffer[block].emplace_back(ii);
			}

			if(svar.Asource == 0)
			{
				pn[ii].cellRho = avar.rhog;
				pn[ii].cellP = avar.pRef;
			}
			
			// Initialise the rest of the values to 0
			pn[ii].surf = 0;
			pn[ii].surfzone = 0;
			pn[ii].nFailed = 0;

			pn[ii].s = 0.0;
			pn[ii].woccl = 0.0;
			pn[ii].pDist = 0.0;
			pn[ii].deltaD = 0.0;
			pn[ii].internal = 0;

			pn[ii].bNorm = StateVecD::Zero();
			pn[ii].y = 0.0;
			pn[ii].vPert = StateVecD::Zero();

			pID = std::max(pID, pn[ii].partID);
		}	

		size_t nBuff = limits[block].hcpl == 1 ? 5 : 4;
		if(buffer[block].size() != limits[block].back.size()*nBuff)
		{
			cout << "Mismatch of back vector size and buffer vector size" << endl;
			cout << "Buffer vector should be 4/5x back vector size." << endl;
			cout << "Sizes are:  Buffer: " << buffer[block].size() << 
					"  Back: " << limits[block].back.size() << endl;
			exit(-1);
		}

		std::sort(limits[block].back.begin(),limits[block].back.end());
		std::sort(buffer[block].begin(),buffer[block].end());
		limits[block].buffer = vector<vector<size_t>>(limits[block].back.size(),
					vector<size_t>(nBuff, 0));
	}
	svar.partID = pID+1;

	/* Put the particles into the buffer vector */
	real eps = svar.restart_tol*svar.dx; /* Tolerance value */
	eps = eps*eps;
	for(size_t block = 0; block < limits.size(); ++block)
	{
		size_t nBack = limits[block].back.size();
		size_t nBuffer = buffer[block].size();
		size_t nDeep = limits[block].hcpl == 1 ? 5 : 4;
		size_t nFound = 0;
		vector<vector<size_t>> found(nBack,vector<size_t>(nDeep,0));
		for(size_t ii = 0; ii < nBack; ++ii)
		{	
			// Use a plane test this time. 
			StateVecD insnorm = limits[block].insert_norm.normalized();
			StateVecD xi = pn[limits[block].back[ii]].xi;
			
			for(size_t index = 0; index < nBuffer; ++index)
			{
				/* Check which particle in the back vector it corresponds to by checking which  */
				/* particle it lies behind */
				StateVecD const& test = pn[buffer[block][index]].xi;

				// If point lies within a given radius tolerance of the expected position, assign.
				if(((xi - insnorm * svar.dx)-test).squaredNorm() < eps)
				{
					limits[block].buffer[ii][0] = buffer[block][index];
					found[ii][0]++;
					nFound++;
				}

				if(((xi - 2.0 * insnorm * svar.dx)-test).squaredNorm() < eps)
				{
					limits[block].buffer[ii][1] = buffer[block][index];
					found[ii][1]++;
					nFound++;
				}
				
				if(((xi - 3.0 * insnorm * svar.dx)-test).squaredNorm() < eps)
				{
					limits[block].buffer[ii][2] = buffer[block][index];
					found[ii][2]++;
					nFound++;
				}
				
				if(((xi - 4.0 * insnorm * svar.dx)-test).squaredNorm() < eps)
				{
					limits[block].buffer[ii][3] = buffer[block][index];
					found[ii][3]++;
					nFound++;
				}
				
				if(((xi - 5.0 * insnorm * svar.dx)-test).squaredNorm() < eps)
				{
					limits[block].buffer[ii][4] = buffer[block][index];
					found[ii][4]++;
					nFound++;
				}
			}
		}	

		for(size_t ii = 0; ii < nBack; ++ii)
		{
			for(size_t jj = 0; jj < nDeep; ++jj)
			{
				if(found[ii][jj] == 0)
				{
					cout << "Buffer particle " << ii << ", " << jj << " has not been assigned. " << endl;
				}

				if(found[ii][jj] > 1)
				{
					cout << "Buffer particle " << ii << ", " << jj << " has been assigned more than once. " << endl;
				}
			}
		}

		if(nFound != nBuffer)
		{
			cout << "Some buffer particles have not been attached to a back particle. Cannot continue" << endl;
			cout << "Particles located: " << nFound << "  Particles failed: " << nBuffer - nFound << endl;
			exit(-1);
		}
	}
	pnp1 = pn;
}
