/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H

#include "Third_Party/Eigen/LU"
#include "Var.h"
#include "IOFunctions.h"
#include "BinaryIO.h"
#include "Geometry.h"
#include <ctime>
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

	if(avar.vJetMag != -1)
	{
		avar.vStart = StateVecD::Zero();
		avar.vStart[1] = avar.vJetMag;
		avar.vStart = svar.Rotate*avar.vStart;
	}

	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Cole's Eq*/

	/*Pipe Pressure calc*/
	fvar.rhoJ = density_equation(fvar.pPress,fvar.B,fvar.gam,fvar.Cs,fvar.rho0);
	
	/* Upper and lower limits for density */
	if(fvar.rhoMax == 1500 && fvar.rhoMin == 500)
	{	/* If limits are undefined, use a variation around the base density */
		fvar.rhoMax = fvar.rho0*(1.0 + fvar.rhoVar*0.01);
		fvar.rhoMin = fvar.rho0*(1.0 - fvar.rhoVar*0.01);
	}
	
	svar.nrad = 1;
	if(svar.Scase == 1)
	{
		svar.dx = svar.Pstep;
	}
	else if(svar.Scase == 2 || svar.Scase == 4)
	{	/* Cylinder case */
		// Defining dx to fit the pipe, then find rest spacing
		if(svar.jet_diam/svar.Pstep < 1.5)
		{	/*spacing is too close to full size to use standard adjustment*/
			cout << "Warning: particle spacing if of the same order of magintude of the jet diameter." << endl;
			cout << "Consider a different size for accuracy." << endl;
		}
		else
		{
			svar.nrad = ceil(abs(0.5*svar.jet_diam/svar.Pstep));
			svar.dx = 0.5*(svar.jet_diam)/real(svar.nrad);
		}	

		// avar.dPipe = svar.Jet(0);
	}
	else if (svar.Scase == 3)
	{	/* Droplet case */
		if(svar.diam/svar.Pstep < 1.5)
		{	/*spacing is too close to full size to use standard adjustment*/
			svar.nrad = 1;
			svar.dx = svar.diam;
		}
		else
		{
			real radius = 0.5*svar.diam;
		
			svar.nrad = ceil(abs(radius/svar.Pstep));
			svar.dx = radius/real(svar.nrad);
		}
	}
	// else
	// {
	// 	if(svar.jet_diam/svar.Pstep < 1.5)
	// 	{	/*spacing is too close to full size to use standard adjustment*/
	// 		svar.nrad = 1;
	// 		svar.dx = svar.jet_diam;
	// 	}
	// 	else
	// 	{
	//  		svar.nrad = ceil(abs(0.5*svar.jet_diam/svar.Pstep));
	//  		svar.dx = 0.5*(svar.jet_diam)/real(svar.nrad);
	// 	}
	// }

	// svar.dx = svar.Pstep;	
	svar.Pstep = svar.dx * pow(fvar.rhoJ/fvar.rho0,1.0/SIMDIM);
 	
	// Correct the droplet to have the same volume as the original
	if(svar.Scase == 3)
	{
		cout << "Droplet Diameter: " << svar.diam << endl;		
	}


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

	if(svar.Scase == 1)
	{
		#if SIMDIM == 2
		svar.simPts = svar.xyPART[0]*svar.xyPART[1];
		#else
		svar.simPts = svar.xyPART[0]*svar.xyPART[1]*svar.xyPART[2]; /*total sim particles*/
		#endif
	}

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
	cout << "Output type: " << svar.outform << endl;
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
	if(svar.Asource != 1 || svar.Asource != 2)
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

    string line;
    while (getline(fin,line))
    {
		line = ltrim(line);
        size_t end = line.find_first_of('#');
		if(end != std::string::npos)
			line = line.substr(0,end+1);

        /* File Inputs */
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
		#if SIMDIM == 3
		Get_String(line, "VLM definition filename", svar.vlm_file); 
		#endif 

        /* File outputs */
		Get_String(line, "Output files prefix", svar.output_prefix);
        Get_Number(line, "SPH frame time interval", svar.framet);
        Get_Number(line, "SPH frame count", svar.Nframe);
        Get_Number(line, "SPH output information (0/1/2/3/4/5)", svar.outform);
        Get_Number(line, "SPH output encoding (0=ascii/1=binary)", svar.out_encoding);
        Get_Number(line, "SPH ghost output (0/1)", svar.gout);
        
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
		Get_Number(line, "SPH starting pressure", fvar.pPress);
		Get_Number(line, "SPH density variation (%)", fvar.rhoVar);
		Get_Number(line, "SPH maximum density", fvar.rhoMax);
		Get_Number(line, "SPH minimum density", fvar.rhoMin);
		Get_Number(line, "SPH artificial viscosity factor", fvar.alpha);
		Get_Number(line, "SPH speed of sound", fvar.Cs);
        Get_Number(line, "SPH Newmark-Beta iteration limit", svar.subits);
		Get_Vector(line, "SPH gravity vector", svar.grav);

        Get_Number(line, "SPH initial spacing", svar.Pstep);
        Get_Number(line, "SPH boundary spacing factor", svar.Bstep);
        Get_Number(line, "SPH smoothing length factor", fvar.Hfac);
        Get_String(line, "SPH aerodynamic case", avar.aero_case);
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
        Get_Array(line, "Droplet resolutions", svar.nacross);
		Get_Array(line, "Droplet diameters", svar.diameters);
		Get_Array(line, "Droplet velocities", svar.velocities);
		Get_Array(line, "Droplet Reynolds numbers", svar.Reynolds);
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

	if(svar.start_type == "box")
	{
		svar.Scase = 1;

		if(svar.xyPART[0] < 0 || svar.xyPART[1] < 0
			#if SIMDIM == 3
				|| svar.xyPART[2] < 0
			#endif
			)
		{
			if(svar.Pstep < 0)
			{
				cout << "SPH box resolution or initial spacing has not been defined." << endl;
				exit(-1);
			}	
		}
		else
		{
			real dx = svar.sim_box[0]/svar.xyPART[0];
			real dy = svar.sim_box[1]/svar.xyPART[1];
			svar.Pstep = std::min(dx,dy);
		}
		
	}
	else if(svar.start_type == "cylinder")
	{
		svar.Scase = 2;

		if(svar.jet_diam < 0)
		{
			cout << "SPH cylinder diameter has not been defined." << endl;
			exit(-1);
		}

		if(svar.jet_depth < 0)
		{
			cout << "SPH cylinder depth has not been defined." << endl;
			exit(-1);
		}

	}	
	else if(svar.start_type == "sphere")
	{
		svar.Scase = 3;
		if(svar.diam < 0)
		{
			cout << "SPH sphere diameter has not been defined." << endl;
			exit(-1);
		}
	}
	else if (svar.start_type == "jet")
	{
		svar.Scase = 4;

		if(svar.jet_diam < 0)
		{
			cout << "SPH cylinder diameter has not been defined." << endl;
			exit(-1);
		}

		if(svar.jet_depth < 0)
		{
			cout << "SPH cylinder depth has not been defined." << endl;
			exit(-1);
		}
	}
	else
	{
		cout << "SPH starting geometry is not defined correctly if at all." << endl;
		exit(-1);
	}


	if(svar.bound_type == "(none)")
	{
		svar.Bcase = 0;
	}
	else if (svar.bound_type == "box")
	{
		svar.Bcase = 1;

		if (svar.bound_start[0] == 999999 || svar.bound_start[1] == 999999
			#if SIMDIM == 3
			|| svar.bound_start[2] == 999999
			#endif
		)
		{
			cout << "Some or all of the SPH boundary start coordinates have not been defined." << endl;
			exit(-1);
		}

		if(svar.bound_box[0] < 0 || svar.bound_box[1] < 0
			#if SIMDIM == 3
				|| svar.bound_box[2] < 0
			#endif
			)
		{
			cout << "Some or all of the SPH boundary box lengths have not been defined." << endl;
			exit(-1);
		}
	}
	else if (svar.bound_type == "jet")
	{
		svar.Bcase = 2;
		if(svar.Scase != 4)
		{
			cout << "Boundary geometry is not compatible with starting geometry." << endl;
			cout << "Starting geometry must be a cylinder to use a jet boundary." << endl;
			exit(-1);
		}

		if(avar.vJetMag < 0)
		{
			cout << "Jet velocity has not been defined" << endl;
			exit(-1);
		}

		svar.bound_start = svar.sim_start;
	}
	else if (svar.bound_type == "taper_jet")
	{
		svar.Bcase = 3;
		if(svar.Scase != 4)
		{
			cout << "Boundary geometry is not compatible with starting geometry." << endl;
			cout << "Starting geometry must be a cylinder to use a jet boundary." << endl;
			exit(-1);
		}

		if(avar.vJetMag < 0)
		{
			cout << "Jet velocity has not been defined" << endl;
			exit(-1);
		}

		svar.bound_start = svar.sim_start;
	}
	else
	{
		cout << "SPH boundary geometry is not defined correctly if at all." << endl;
		exit(-1);
	}

	if(svar.Bcase != 0)
	{
		if(svar.bound_solver > 1)
		{
			cout << "Boundary solver not defined correctly." << endl;
			exit(-1);
		}
	}

	/* Check for restart data now boundary case is known */	
	if(!svar.restart_prefix.empty())
	{
		svar.restart = 1;
		Check_If_Restart_Possible(svar);
		svar.output_prefix = svar.restart_prefix;
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
						svar.nacross.emplace_back(ceil(abs(svar.diam/dx)));
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
	
  	Set_Values(svar,fvar,avar);

	Print_Settings(argv,svar,fvar,avar);

} /*End of GetInput()*/

/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/
void Write_ASCII_header(std::fstream& fp, SIM &svar)
{
	string variables;

#if SIMDIM == 2
	variables = "\"X\", \"Z\"";  
	if (svar.outform == 1)
	{
		variables = "\"X\", \"Z\", \"rho\", \"P\", \"m\", \"v\", \"a\"";
	}
	else if (svar.outform == 2)
	{
		variables = "\"X\", \"Z\", \"rho\", \"P\", \"m\", \"v_x\", \"v_z\", \"a_x\", \"a_z\"";
	}
	else if (svar.outform == 3)
	{
		variables = 
"\"X\", \"Z\", \"rho\", \"P\", \"m\", \"v_x\", \"v_z\", \"a_x\", \"a_z\", \"Cell_Vx\", \"Cell_Vz\", \"Cell_P\", \"Cell_ID\"";
	}
	else if (svar.outform == 4)
	{
		variables = "\"X\", \"Z\", \"rho\", \"P\", \"m\", \"v\", \"a\", \"b\", \"Neighbours\", \"Aero\"";
	}

#endif

#if SIMDIM == 3
	variables = "\"X\", \"Y\", \"Z\"";  
	if (svar.outform == 1)
	{
		variables = "\"X\", \"Y\", \"Z\", \"rho\", \"P\", \"m\", \"v\", \"a\"";
	}
	else if (svar.outform == 2)
	{
		variables = "\"X\", \"Y\", \"Z\", \"rho\", \"P\", \"m\", \"v_x\", \"v_y\", \"v_z\", \"a_x\", \"a_y\", \"a_z\"";
	}
	else if (svar.outform == 3)
	{
		variables = 
"\"X\", \"Y\", \"Z\", \"rho\", \"P\", \"m\", \"v_x\", \"v_y\", \"v_z\", \"a_x\", \"a_y\", \"a_z\", \"Cell_Vx\", \"Cell_Vy\", \"Cell_Vz\", \"Cell_P\", \"Cell_ID\"";
	}
	else if (svar.outform == 4)
	{
		variables = "\"X\", \"Y\", \"Z\", \"rho\", \"P\", \"m\", \"v\", \"a\", \"b\", \"Neighbours\", \"Aero\"";
	}

#endif

	fp << "VARIABLES = " << variables << "\n";
}


void Write_ASCII_Timestep(std::fstream& fp, SIM& svar, SPHState const& pnp1, 
	const uint bwrite, const uint start, const uint end, string const& name)
{
	// if(bwrite == 1)
	//  	fp <<  "ZONE T=\"" << name << "\"";
	// else
		fp <<  "ZONE T=\"" << name << "\"";

    fp <<", I=" << end - start << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
    fp << std::left << std::scientific << std::setprecision(6);
    const static uint width = 15;
    
    if(svar.outform == 0)
    {
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
	        	fp << setw(width) << pnp1[ii].xi(dim)/svar.scale; 
			fp << "\n";  
	  	}

	  	// if (airP.size() > 0 )
	  	// {
		  // 	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
		  //   ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
		  // 	for(auto p:airP)
		  // 	{
		  // 		for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.xi(i) << " "; 
				// fp << "\n"; 
		  // 	}
	  	// }
	  	fp << std::flush;
	}	
    if(svar.outform == 1)
    {
		for (size_t ii = start; ii < end; ++ii)
		{	
			for(uint dim = 0; dim < SIMDIM; ++dim)
	        	fp << setw(width) << pnp1[ii].xi(dim)/svar.scale;

			fp << setw(width) << pnp1[ii].rho << setw(width) << pnp1[ii].p;
			fp << setw(width) << pnp1[ii].m;
	        fp << setw(width) << pnp1[ii].v.norm();
	        fp << setw(width) << pnp1[ii].acc.norm() << "\n";
	  	}

	  // 	if (airP.size() > 0 )
	  // 	{
		 //  	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
		 //    ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";

		 //  	for(auto p:airP)
		 //  	{
		 //  		for(uint i = 0; i < SIMDIM; ++i)
		 //        	fp << p.xi(i) << " "; 

		 //        fp << p.v.norm() << " ";
		 //        fp << p.acc.norm() << " ";
		 //        fp << p.rho << " "  << p.p  << "\n";
		 //  	}
	 	// }
	  	fp << std::flush; 	
    }
    else if (svar.outform == 2)
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
	        	fp << setw(width) << pnp1[ii].xi(dim)/svar.scale;
	        
	        fp << setw(width) << pnp1[ii].rho << setw(width) << pnp1[ii].p;
			fp << setw(width) << pnp1[ii].m;

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << pnp1[ii].v(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
				fp << setw(width) << pnp1[ii].acc(i); 

			fp << "\n";
	  	}  

	  	// if (airP.size() > 0 )
	  	// {
		  // 	fp <<  "ZONE T=\"Air Data\"" <<", I=" << airP.size() << ", F=POINT" <<
		  //   ", STRANDID=2, SOLUTIONTIME=" << svar.t  << "\n";
		    
		  // 	for(auto p:airP)
		  // 	{
		  // 		for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.xi(i) << " "; 

		  //     	fp << p.acc.norm() << " " << p.Af.norm() << " " << p.Sf.norm() << " ";
		  //       for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.cellV(i) << " "; 

		  //       fp << p.b << " " << p.theta  << "\n"; 
		  // 	}
	  	// }	
	  	fp << std::flush;
    }
    else if (svar.outform == 3)
    {
    	for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
	        	fp << setw(width) << pnp1[ii].xi(dim)/svar.scale;
	        
	        fp << setw(width) << pnp1[ii].rho << setw(width) << pnp1[ii].p;
			fp << setw(width) << pnp1[ii].m;

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << pnp1[ii].v(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << pnp1[ii].acc(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
				fp << setw(width) << pnp1[ii].cellV(i);

			fp << setw(width) << pnp1[ii].cellP;
	        fp << setw(width) << pnp1[ii].cellID << "\n"; 
	  	}
	  	fp << std::flush;
    }
    else if (svar.outform == 4)
    {
    	for (size_t ii = start; ii < end; ++ii)
		{	
			for(uint dim = 0; dim < SIMDIM; ++dim)
	        	fp << setw(width) << pnp1[ii].xi(dim)/svar.scale;

			fp << setw(width) << pnp1[ii].rho << setw(width) << pnp1[ii].p;
			fp << setw(width) << pnp1[ii].m;
	        fp << setw(width) << pnp1[ii].v.norm();
	        fp << setw(width) << pnp1[ii].acc.norm();
	        fp << setw(width) << pnp1[ii].b << setw(width) << pnp1[ii].s;
	        fp << setw(width) << pnp1[ii].Af.norm() << "\n";
	  	}
	  	fp << std::flush;

    }
}

void Write_First_Step(std::fstream& f1, std::fstream& fb, std::fstream& fg, 
						SIM& svar, SPHState const& pnp1/* , SPHState const& airP */)
{
	if(svar.out_encoding == 1 )
	{
		/*Write sim particles*/
		Init_Binary_PLT(svar,"_fuel.szplt","Simulation Particles",svar.fuelFile);
		if(svar.restart == 0)
		{
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1,svar.fuelFile);
		}

		if (svar.Bcase != 0)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"_boundary.szplt","Boundary Particles",svar.boundFile);

			if(svar.restart == 0)
				Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2,svar.boundFile); 
		}	

		// if (svar.ghost == 1 && svar.gout == 1)
		// {
		// 	/*Write boundary particles*/
		// 	Init_Binary_PLT(svar,"_ghost.szplt","Ghost Particles",svar.ghostFile);		
		// }	
	}
	else
	{

		if (svar.Bcase != 0 && svar.restart == 0)
		{	/*If the boundary exists, write it.*/
			string bfile = svar.output_prefix;
			bfile.append("_boundary.plt");
			fb.open(bfile, std::ios::out);
			if(fb.is_open())
			{
				Write_ASCII_header(fb,svar);
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
			}
			else
			{
				cerr << "Error opening boundary file." << endl;
				exit(-1);
			}
		}

		/* Write first timestep */
		string mainfile = svar.output_prefix;
		mainfile.append("_fuel.plt");
		f1.open(mainfile, std::ios::out);
		if(f1.is_open())
		{
			Write_ASCII_header(f1,svar);
			Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,"Fuel");
		}
		else
		{
			cerr << "Failed to open fuel.plt. Stopping." << endl;
			exit(-1);
		}

		// if(svar.ghost == 1 && svar.gout == 1)
		// {
		// 	string ghostfile = svar.output_prefix;
		// 	ghostfile.append("_ghost.plt");
		// 	fg.open(ghostfile,std::ios::out);
		// 	if(fg.is_open())
		// 	{
		// 		Write_ASCII_header(fg,svar);
		// 		Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
		// 	}
		// }
	}

	if(svar.restart == 0)
	{
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
}

void Write_Timestep(std::fstream& f1, std::fstream& fb, std::fstream& fg, 
				SIM& svar, SPHState const& pnp1/* , SPHState const& airP */)
{
	if (svar.out_encoding == 1)
	{
		if(svar.Bcase != 0)
		{	/*Write boundary particles*/
			Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2,svar.boundFile); 
		}
		Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1,svar.fuelFile); /*Write sim particles*/
		// if(svar.ghost != 0 && svar.gout == 1 && airP.size() != 0)
		// 	Write_Binary_Timestep(svar,airP,0,airP.size(),"Ghost",3,svar.ghostFile);
	} 
	else
	{
		Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,"Fuel");
		if(svar.Bcase != 0)
		{
			SPHState empty;
			Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
		}

		// if(svar.ghost != 0 && svar.gout == 1 && airP.size() != 0)
		// 	Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
	}
}

#endif