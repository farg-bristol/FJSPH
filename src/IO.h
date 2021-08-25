/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H

#include "Third_Party/Eigen/LU"
#include "Var.h"
#include "IOFunctions.h"
#include "CDFIO.h"
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

	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/

	/*Pipe Pressure calc*/
	fvar.rhoJ = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);

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
	avar.qInf = 0.5*avar.rhog*avar.vRef*avar.vRef;

	svar.addcount = 0;
  	svar.dt = 2E-010; 			/*Initial timestep*/
  	fvar.H = 2.0*svar.Pstep;
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
			svar.vortex.Init(svar.infile);
			svar.vortex.GetGamma(avar.vInf);
		}
	#endif

		GetYcoef(avar, fvar, /*fvar.H*/ svar.Pstep);
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
        if(line[0] == '#') /* Skip commented line */
            continue;

        /* File Inputs */
        Get_String(line, "Primary grid face filename", svar.taumesh);
        Get_String(line, "Boundary mapping filename", svar.taubmap);
        Get_String(line, "Restart-data prefix", svar.tausol);
        Get_String(line, "SPH restart-data prefix", svar.restart_prefix);
        Get_Number(line, "Grid scale", svar.scale);
        Get_Number(line, "2D offset vector (0 / x=1,y=2,z=3)",svar.offset_axis);
        Get_String(line, "OpenFOAM folder", svar.foamdir);
        Get_String(line, "OpenFOAM solution folder", svar.foamsol);
        Get_Number(line, "OpenFOAM binary (0/1)", svar.isBinary);
        Get_Number(line, "Label size (32/64)", svar.labelSize);
        Get_Number(line, "Scalar size (32/64)", svar.scalarSize);

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
        Get_Number(line, "Reference velocity", avar.vRef);
        Get_Number(line, "Reference pressure", avar.pRef);
        Get_Number(line, "Reference Mach number", avar.MRef);
        Get_Number(line, "Reference density", avar.qInf);
        Get_Number(line, "Reference temperature", avar.T);

		// Get_Number(line, "SPH surface tension contact angle", fvar.contangb);

        /* Simulation settings */
		Get_Number(line, "SPH maximum timestep", svar.dt_max);
		Get_Number(line, "SPH minimum timestep", svar.dt_min);
		Get_Number(line, "SPH starting pressure", fvar.pPress);
		Get_Number(line, "SPH speed of sound", fvar.Cs);
		Get_Number(line, "SPH artificial viscosity factor", fvar.alpha);
        Get_Number(line, "SPH Newmark-Beta iteration limit", svar.subits);
        Get_Number(line, "SPH initial spacing", svar.Pstep);
        Get_Number(line, "SPH boundary spacing factor", svar.Bstep);
        Get_String(line, "SPH aerodynamic case", avar.aero_case);
        Get_Number(line, "SPH use ghost particles (0/1)", svar.ghost);
        Get_Vector(line, "SPH start coordinates", svar.sim_start);
		Get_Vector(line, "SPH boundary start coordinates", svar.bound_start);
		Get_Number(line, "SPH maximum particle count",svar.finPts);

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
        Get_Number(line, "Velocity equation order (1/2)", svar.eqOrder);
        Get_Number(line, "SPH tracking conversion x coordinate", svar.max_x_sph);
        Get_Number(line, "Maximum x trajectory coordinate", svar.max_x);
		Get_Number(line, "Particle scatter output (0/1)", svar.partout);
        Get_Number(line, "Particle streak output (0/1)", svar.streakout);
        Get_Number(line, "Particle cell intersection output (0/1)", svar.cellsout);
    }

    fin.close();

	/* Need to check if inputs are correct */
	
	if(svar.taumesh.empty())
    {
		
		if(svar.foamdir.empty())
		{
			svar.Asource = 0;	
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

		if(svar.tausol.empty())
		{
			cout << "Input TAU solution file not defined." << endl;
			exit(-1);
		}

	}

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
			cout << "Some or all of the SPH box resolutions have not been defined." << endl;
			exit(-1);
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
	/* Particle Tracking Settings */
    // if(svar.partout == 1)
    // {
    //     if(svar.outdir.empty())
    //     {
    //         cout << "Output particle directory not defined." << endl;
    //         exit(-1);
    //     }
    // }

    // if(svar.streakout == 1)
    // {
    //     if(svar.streakdir.empty())
    //     {
    //         cout << "Output particle streaks directory not defined." << endl;
    //         exit(-1);
    //     }
    // }

    // if(svar.cellsout == 1)
    // {
    //     if(svar.celldir.empty())
    //     {
    //         cout << "Output cell intersections directory not defined." << endl;
    //         exit(-1);
    //     }
    // }

    if(svar.eqOrder > 2 || svar.eqOrder < 1)
    {
        cout << "Equation order not 1 or 2. Please choose between these." << endl;
        exit(-1);
    }

    // svar.streakdir = svar.outfile;
    // size_t pos = svar.streakfile.find_last_of(".");
    // svar.streakfile.insert(pos,"_streak");
    
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
						SIM& svar, SPHState const& pnp1, SPHState const& airP)
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

		if (svar.ghost == 1 && svar.gout == 1)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"_ghost.szplt","Ghost Particles",svar.ghostFile);		
		}	
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

		if(svar.ghost == 1 && svar.gout == 1)
		{
			string ghostfile = svar.output_prefix;
			ghostfile.append("_ghost.plt");
			fg.open(ghostfile,std::ios::out);
			if(fg.is_open())
			{
				Write_ASCII_header(fg,svar);
				Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
			}
		}
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
		para << "                                          " << 
				"SPH restart-data prefix: " << svar.output_prefix << endl;

		para.close();
	}
}

void Write_Timestep(std::fstream& f1, std::fstream& fb, std::fstream& fg, 
				SIM& svar, SPHState const& pnp1, SPHState const& airP)
{
	if (svar.out_encoding == 1)
	{
		if(svar.Bcase != 0)
		{	/*Write boundary particles*/
			Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2,svar.boundFile); 
		}
		Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1,svar.fuelFile); /*Write sim particles*/
		if(svar.ghost != 0 && svar.gout == 1 && airP.size() != 0)
			Write_Binary_Timestep(svar,airP,0,airP.size(),"Ghost",3,svar.ghostFile);
	} 
	else
	{
		Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,"Fuel");
		if(svar.Bcase != 0)
		{
			SPHState empty;
			Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
		}

		if(svar.ghost != 0 && svar.gout == 1 && airP.size() != 0)
			Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
	}
}


/* Particle tracking functions */
namespace IPT
{
	void Write_ASCII_variables(ofstream& fp)
	{
		string variables = 
	"\"X\", \"Y\", \"Z\", \"t\", \"dt\", \"v\", \"a\", \"ptID\", \"Cell_V\", \"Cell_Rho\", \"Cell_ID\"";

		fp << "VARIABLES = " << variables << "\n";
	}

	void Write_ASCII_Scatter_Header(ofstream& fp)
	{
		fp << "TITLE =\"IPT particle scatter data\" \n";
		fp << "F=POINT\n";
		fp << std::left << std::scientific << std::setprecision(6);
	}

	void Write_ASCII_Streaks_Header(ofstream& fout)
	{
		fout << "TITLE = \"IPT Streaks\"\n";
		Write_ASCII_variables(fout);
	}

	void Write_ASCII_Cells_Header(ofstream& fout)
	{
		fout << "TITLE = \"IPT intersecting cells\"\n";
		fout << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
	}
	
	void Init_IPT_Files(SIM& svar)
	{
		string partf, cellf, streakf, surfacef;

		partf = svar.output_prefix;
		partf.append("_IPT_scatter.dat");

		streakf = svar.output_prefix;
		streakf.append("_IPT_streaks.dat");

		cellf = svar.output_prefix;
		cellf.append("_IPT_cells.dat");

		surfacef = svar.output_prefix;
		surfacef.append("_surface_impacts.dat");

		if(svar.partout == 1)
		{
			svar.partfile.open(partf,std::ios::in);
			
			if(svar.partfile.is_open())
			{
				Write_ASCII_Scatter_Header(svar.partfile);
			}
			else
			{
				cout << "Couldn't open the IPT particle scatter output file." << endl;
				exit(-1);
			}
		}

		if(svar.streakout == 1)
		{
			svar.streakfile.open(streakf,std::ios::in);

			if(svar.streakfile.is_open())
			{
				Write_ASCII_Streaks_Header(svar.streakfile);
			}
			else
			{
				cout << "Couldn't open the IPT streaks output file." << endl;
				exit(-1);
			}
		}

		if(svar.cellsout == 1)
		{
			svar.cellfile.open(cellf,std::ios::in);

			if(svar.cellfile.is_open())
			{
				Write_ASCII_Cells_Header(svar.cellfile);
			}
			else
			{
				cout << "Couldn't open the IPT cell intersection output file." << endl;
				exit(-1);
			}
		}

		// svar.surfacefile.open(surfacef,std::ios::in);
	}

	void Write_ASCII_Point(ofstream& fp, real const& scale, IPTPart const& pnp1)
	{
		const static uint width = 15;

		for(uint dim = 0; dim < 3; ++dim)
			fp << setw(width) << pnp1.xi[dim]/scale;

		fp << setw(width) << pnp1.t; 

		fp << setw(width) << pnp1.dt; 
		
		fp << setw(width) << pnp1.v.norm(); 

		fp << setw(width) << pnp1.acc; 

		fp << setw(width) << pnp1.partID;

		fp << setw(width) << pnp1.cellV.norm();

		fp << setw(width) << pnp1.cellRho;
		fp << setw(width) << pnp1.cellID << "\n"; 
	}

	void Write_ASCII_Timestep(ofstream& fp, SIM const& svar, IPTState const& pnp1)
	{
		fp <<  "ZONE T=\"" << "IPT scatter data" << "\"";
		fp <<", I=" << pnp1.size() << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
		fp << std::left << std::scientific << std::setprecision(6);
		
		for (size_t ii = 0; ii < pnp1.size(); ++ii)
		{
			Write_ASCII_Point(fp, svar.scale, pnp1[ii]);
		}
		fp << std::flush;
	}

	void Write_ASCII_Streaks( SIM& svar, IPTState const& t_pnp1, IPTPart const& pnp1)
	{
		ofstream& fp = svar.streakfile;
		if(!fp.is_open())
		{
			cout << "The streak output file is not open. Cannot write." << endl;
			exit(-1);
		}
		
		size_t nTimes = t_pnp1.size();
		size_t time = 0;
		
		fp <<  "ZONE T=\"" << "Ash particle " << pnp1.partID << "\"";
		
		if(pnp1.failed == 0)
			fp <<", I=" << nTimes+1 << ", J=1, K=1, DATAPACKING=POINT"<< "\n";
		else
			fp <<", I=" << nTimes << ", J=1, K=1, DATAPACKING=POINT"<< "\n";

		fp << std::left << std::scientific << std::setprecision(6);

		for (time = 0; time < nTimes; ++time)
		{ /* Inner loop to write the times of the particle */
			Write_ASCII_Point(fp, svar.scale, t_pnp1[time]);
		}

		if(pnp1.failed == 0)
			Write_ASCII_Point(fp, svar.scale, pnp1);

	}


	void Write_ASCII_Cells(SIM& svar, MESH const& cells, IPTState const& t_pnp1)
	{
		size_t nTimes = t_pnp1.size();
		size_t time = 0;

		/* Write file for each particle. Gonna get way too confusing otherwise */
		ofstream& fout = svar.cellfile;

		if(!fout.is_open())
		{
			cout << "Couldn't open the cell intersection output file." << endl;
			exit(-1);
		}

		int TotalNumFaceNodes = 0;

		vector<size_t> vertIndexes;
		vector<size_t> faceIndexes;
		vector<int> cellIndexes;

		vector<StateVecD> usedVerts;
		vector<vector<size_t>> faces;
		vector<int> left;
		vector<int> right;


		/* Get cell properties for the intersected cells (delete duplicates later)*/
		for (time = 0; time < nTimes; ++time)
		{
			int cellID = t_pnp1[time].cellID;
			
			/* Find how many vertices the cell has */
			vertIndexes.insert(vertIndexes.end(), cells.elems[cellID].begin(), cells.elems[cellID].end());
			
			faceIndexes.insert(faceIndexes.end(), cells.cFaces[cellID].begin(), cells.cFaces[cellID].end());

			cellIndexes.emplace_back(cellID);
		}

		/* Delete repeats of vertex mentions*/
		std::sort(vertIndexes.begin(),vertIndexes.end());
		vertIndexes.erase( std::unique( vertIndexes.begin(), vertIndexes.end()), vertIndexes.end() );

		std::sort(faceIndexes.begin(),faceIndexes.end());
		faceIndexes.erase( std::unique( faceIndexes.begin(), faceIndexes.end()), faceIndexes.end() );

		for(size_t const& vert : vertIndexes)
		{
			usedVerts.emplace_back(cells.verts[vert]);
		}

		/* Get faces properties used in the cell */
		for(size_t const& face: faceIndexes )
		{
			faces.emplace_back(cells.faces[face]);
			left.emplace_back(cells.leftright[face].first);
			right.emplace_back(cells.leftright[face].second);
		}

		/* Go through the faces indexes, finding the appropriate index to change it to */

		for(vector<size_t>& face : faces)
		{
			for(size_t& vert : face)
			{
				vector<size_t>::iterator index = std::find(vertIndexes.begin(), vertIndexes.end(), vert);
				if(index != vertIndexes.end())
				{
					vert = index - vertIndexes.begin() + 1 ;
				}
				else
				{
					cout << "Couldn't find the vertex used in the prepared list." << endl;
					exit(-1);
				}                
			}
			TotalNumFaceNodes += face.size();
		}

		/* Find the cell used and change it's index */
		for(int& leftCell:left)
		{
			vector<int>::iterator index = std::find(cellIndexes.begin(), cellIndexes.end(), leftCell);
			if(index != cellIndexes.end())
			{
				leftCell = (index - cellIndexes.begin())+1;
			}
			else
			{   /* Cell isn't intersected, so it won't be drawn. Boundary face. */
				leftCell = 0;
			} 
		}

		for(int& rightCell:right)
		{
			vector<int>::iterator index = std::find(cellIndexes.begin(), cellIndexes.end(), rightCell);
			if(index != cellIndexes.end())
			{
				rightCell = index - cellIndexes.begin() + 1;
			}
			else
			{   /* Cell isn't intersected, so it won't be drawn. Boundary face. */
				rightCell = 0;
			}
		}

		fout << "ZONE T=\"particle " << t_pnp1[0].partID << " intersecting cells\"" << endl;
		fout << "ZONETYPE=FEPOLYHEDRON" << endl;
		fout << "NODES=" << usedVerts.size() << " ELEMENTS=" << cellIndexes.size() << 
				" FACES=" << faces.size() << endl;
		fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
		fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

		size_t w = 15;
		size_t preci = 6;
		fout << std::left << std::scientific << std::setprecision(preci);

		size_t newl = 0;
		fout << std::setw(1);
		for(size_t DIM = 0; DIM < 3; ++DIM)
		{
			for(size_t ii = 0; ii < usedVerts.size(); ++ii)
			{
				fout << std::setw(w) << usedVerts[ii][DIM];
				newl++;

				if(newl>4)
				{
					fout << endl;
					fout << " ";
					newl=0;
				}
			}
		}
		fout << endl;

		fout << std::left << std::fixed;
		w = 9;
		/*Inform of how many vertices in each face*/
		fout << "#node count per face" << endl;
		newl = 0;
		for (size_t ii = 0; ii < faces.size(); ++ii)
		{
			fout << std::setw(w) << faces[ii].size();
			newl++;

			if(newl>4)
			{
				fout << endl;
				newl=0;
			}
		}
		fout << endl;
		/*Write the face data*/
		fout << "#face nodes" << endl;
		for (size_t ii = 0; ii < faces.size(); ++ii)
		{
			for(auto const& vertex:faces[ii])
			{	/*Write face vertex indexes*/
				fout << std::setw(w) << vertex;
				// if (vertex > fdata.nPnts)
				// {
				// 	cout << "Trying to write a vertex outside of the number of points." << endl;
				// }
			}
			fout << endl;
		}

		/*Write face left and right*/
		newl = 0;
		fout << "#left elements" << endl;
		for (size_t ii = 0; ii < left.size(); ++ii)
		{
			fout << std::setw(w) << left[ii] ;
			newl++;

			if(newl>4)
			{
				fout << endl;
				newl=0;
			}
		}
		fout << endl;

		newl = 0;
		fout << "#right elements" << endl;
		for (size_t ii = 0; ii < right.size(); ++ii)
		{
			fout << std::setw(w) << right[ii] ;
			newl++;

			if(newl>4)
			{
				fout << endl;
				newl=0;
			}
		}

		fout.close();
	}


	void Write_ASCII_Impacts(SIM const& svar, FLUID const& fvar, MESH const& cells, 
				vector<SURF> const& surfs, vector<vector<real>> const& beta_data)
	{
		// uint nSurf = svar.markers.size();
		uint nSurf = 0;
		vector<int> marks;
		vector<string> names;
		vector<SURF> surfaces_to_write;
		for(size_t ii = 0; ii < surfs.size(); ii++)
		{
			if(surfs[ii].output == 1)
			{
				marks.emplace_back(svar.markers[ii]);
				names.emplace_back(svar.bnames[ii]);
				surfaces_to_write.emplace_back(surfs[ii]);
			}
			nSurf += surfs[ii].output;
		}

		nSurf = 0;
		for(size_t ii = 0; ii < surfs.size(); ii++)
		{
			if(surfs[ii].output == 1)
			{
				surfaces_to_write[nSurf].face_count = surfs[ii].face_count;
				surfaces_to_write[nSurf].face_beta = surfs[ii].face_beta;
				surfaces_to_write[nSurf].face_area = surfs[ii].face_area;
			}
			nSurf += surfs[ii].output;
		}
		

		vector<vector<vector<size_t>>> faces(nSurf);
		vector<std::pair<size_t,int>> smarkers = cells.smarkers;

		vector<vector<size_t>> vertIndexes(nSurf);
		vector<vector<StateVecD>> usedVerts(nSurf);
		
		/* Do I need to sort the markers? */
		std::sort(smarkers.begin(), smarkers.end(),
		[](std::pair<size_t,int> const& p1, std::pair<size_t,int> const& p2){return p1.second > p2.second;});

		
		for(size_t ii = 0; ii < surfaces_to_write.size(); ++ii)
		{

			for(size_t jj = 0; jj < surfaces_to_write[ii].faceIDs.size(); ++jj )
			{
				size_t faceID = surfaces_to_write[ii].faceIDs[jj];
				faces[ii].emplace_back(cells.faces[faceID]);
			}
		
			/* Get the vertex indexes, delete duplicates, reorder, and recast the face indexes */
			for (vector<size_t> const& face : faces[ii])
				vertIndexes[ii].insert(vertIndexes[ii].end(), face.begin(), face.end());
			

			std::sort(vertIndexes[ii].begin(),vertIndexes[ii].end());
			vertIndexes[ii].erase(std::unique( vertIndexes[ii].begin(), vertIndexes[ii].end()), vertIndexes[ii].end() );

			for(size_t const& vert : vertIndexes[ii])
				usedVerts[ii].emplace_back(cells.verts[vert]);
		

			for(vector<size_t>& face : faces[ii])
			{
				for(size_t& vert : face)
				{
					vector<size_t>::iterator index = std::find(vertIndexes[ii].begin(), vertIndexes[ii].end(), vert);
					if(index != vertIndexes[ii].end())
					{
						vert = index - vertIndexes[ii].begin() + 1 ;
					}
					else
					{
						cout << "Couldn't find the vertex used in the prepared list." << endl;
						exit(-1);
					}                
				}
			}
		}

		
		string file = svar.output_prefix;
		file.append("_surface_impacts.dat");

		ofstream fout(file);

		if(!fout.is_open())
		{
			cout << "Couldn't open the surface impact output file." << endl;
			exit(-1);
		}

		fout << "TITLE=\"Surface collision metrics\"\n";
		fout << "VARIABLES= \"X\" \"Y\" \"Z\" \"Number of Impacts\" \"beta\" \"average area\"\n";
		
		for(size_t ii = 0; ii < faces.size(); ++ii)
		{
			fout << "ZONE T=\"" << names[ii] << "\"\n";
			fout << "N=" << usedVerts[ii].size() << ", E=" << faces[ii].size();
			fout << ", F=FEBLOCK ET=QUADRILATERAL, VARLOCATION=([1-3]=NODAL,[4-7]=CELLCENTERED)\n\n";
			
			
			size_t w = 15;
			size_t preci = 6;
			fout << std::left << std::scientific << std::setprecision(preci);

			size_t newl = 0;
			fout << std::setw(1);
			for(size_t DIM = 0; DIM < SIMDIM; ++DIM)
			{
				for(size_t jj = 0; jj < usedVerts[ii].size(); ++jj)
				{
					fout << std::setw(w) << usedVerts[ii][jj][DIM];
					newl++;

					if(newl>4)
					{
						fout << endl;
						fout << " ";
						newl=0;
					}
				}
			}
			fout << endl << endl;

			/* Variable data goes here */
			fout << "#face impact count" << endl;
			for(size_t jj = 0; jj < surfaces_to_write[ii].face_count.size(); ++jj)
			{
				fout << std::setw(w) << surfaces_to_write[ii].face_count[jj];
				newl++;

				if(newl>4)
				{
					fout << endl;
					fout << " ";
					newl=0;
				}
			}
			fout << endl << endl;

			// for(size_t jj = 0; jj < surfs[surf].mass.size(); ++jj)
			// {
			//     fout << std::setw(w) << surfs[surf].mass[jj];
			//     newl++;

			//     if(newl>4)
			//     {
			//         fout << endl;
			//         fout << " ";
			//         newl=0;
			//     }
			// }
			// fout << endl;
			fout << "#face beta value" << endl; 
			for(size_t jj = 0; jj < surfaces_to_write[ii].face_beta.size(); ++jj)
			{
				fout << std::setw(w) << surfaces_to_write[ii].face_beta[jj];
				newl++;

				if(newl>4)
				{
					fout << endl;
					fout << " ";
					newl=0;
				}
			}
			fout << endl << endl;

			fout << "#face area value" << endl;
			for(size_t jj = 0; jj < surfaces_to_write[ii].face_area.size(); ++jj)
			{
				fout << std::setw(w) << surfaces_to_write[ii].face_area[jj];
				newl++;

				if(newl>4)
				{
					fout << endl;
					fout << " ";
					newl=0;
				}
			}
			fout << endl << endl;


			/*Write the face data*/
			fout << "#face nodes" << endl;
			for (size_t jj = 0; jj < faces[ii].size(); ++jj)
			{
				for(auto const& vertex:faces[ii][jj])
				{	/*Write face vertex indexes*/
					fout << std::setw(w) << vertex;
					// if (vertex > fdata.nPnts)
					// {
					// 	cout << "Trying to write a vertex outside of the number of points." << endl;
					// }
				}

				if(faces[ii][jj].size() == 3)
					fout << std::setw(w) << faces[ii][jj].back();

				fout << endl;
			}
		}
		fout.close();
	}   
}

#endif