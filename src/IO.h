/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IO_H
#define IO_H


#include "Eigen/Core"
#include "Eigen/StdVector"
#include "Eigen/LU"
#include "Var.h"
#include "IOFunctions.h"
#include "CDFIO.h"
// #include "TauIO.h"


/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void Read_SIM_Var(string& infolder, SIM& svar, FLUID& fvar, AERO& avar)
{
	svar.infolder = infolder;
	string file = infolder;
	file.append("Settings");
#ifdef DEBUG
	dbout << "Reading settings file. Path:" << endl << file << endl;
#endif
	std::ifstream in(file);
  	if(!in.is_open()) 
  	{	
  		cerr << "Error opening the settings file." << endl;
	    exit(-1);
  	}

	/*Simulation parameters*/
	cout << "Input file opened. Reading settings..." << endl;
  	uint lineno = 0;
  	svar.scale = 1.0;
	svar.framet = getDouble(in, lineno, "Frame time");
	svar.Nframe = getInt(in, lineno, "Number of frames");
	svar.outframe = getInt(in, lineno, "Output frame info");
	svar.outtype = getInt(in, lineno, "Output data type");
	svar.outform = getInt(in, lineno, "Output contents");
	svar.boutform = getInt(in, lineno, "Boundary time output");
	svar.gout = getInt(in, lineno, "Output ghost particles to file");
	svar.subits = getInt(in, lineno, "Max sub iterations");
	svar.nmax = getInt(in, lineno, "Max particle add rounds");
	/*Get post processing options*/
	svar.cellSize = getDouble(in, lineno, "Post processing mesh size");
	svar.postRadius = getDouble(in, lineno, "Post processing support radius");
	/*Particle settings*/	
	svar.Pstep = getDouble(in, lineno, "Particle initial spacing");
	svar.Bstep = getDouble(in, lineno, "Boundary spacing factor");
	svar.Bcase = getInt(in, lineno, "Simulation initial case");
	svar.Asource = getInt(in, lineno, "Simulation aerodynamic solution source");
	avar.acase = getInt(in, lineno, "Simulation aerodynamic case");
	svar.ghost = getInt(in, lineno, "Ghost particles?");
	svar.Start = getDVector(in, lineno, "Starting position");
	if(svar.Bcase < 2)
	{
		svar.xyPART = getIVector(in, lineno, "Particles in each coordinate");
		svar.Box= getDVector(in, lineno, "Box dimensions");
		fvar.pPress = getDouble(in, lineno, "Pipe pressure");
		if(svar.Box(0) < 0 || svar.Box(1) < 0 
		#if SIMDIM == 3
			|| svar.Box(2) < 0
		#endif
		)
		{
			cout << "Box dimensions are negative. Please check the input and try again." << endl;
			cout << "Line " << lineno << endl;
			exit(-1);
		}
	}
	else if(svar.Bcase > 1 && svar.Bcase < 8)
	{	
		StateVecD angles = getDVector(in, lineno, "Starting angle");
		svar.Angle = angles;
		angles = angles*M_PI/180;
		svar.Rotate = GetRotationMat(angles);
		svar.Transp = svar.Rotate.transpose();
		svar.Jet = getvector(in, lineno, "Jet dimensions");
		fvar.pPress = getDouble(in, lineno, "Pipe pressure");
		avar.vJet = StateVecD::Zero(); avar.vInf = StateVecD::Zero();
		avar.vJet(1) = getDouble(in, lineno, "Jet velocity");  
		avar.vJet = svar.Rotate*avar.vJet;
		avar.vInf = getDVector(in, lineno, "Freestream velocity");
		if(avar.acase == 2 || avar.acase == 3)
		{
			avar.a = getDouble(in, lineno, "a");
			avar.h1 = getDouble(in, lineno, "h1");
			avar.b = getDouble(in, lineno, "b");
			avar.h2 = getDouble(in, lineno, "h2");
		}
		if(avar.acase > 6)
		{
			cout << "Aerodynamic case is not in design. Stopping..." << endl;
			exit(-1);
		}
	}
	else
	{
		cout << "Boundary case not within design. Stopping." << endl;
		exit(-1);
	}
	

#ifdef DEBUG
	dbout << "Closing settings file"  << endl;
#endif
	in.close();	
}

void Read_FLUID_Var(string& infolder, SIM& svar, FLUID& fvar, AERO& avar)
{
	string file = infolder;
	file.append("Fluid");
#ifdef DEBUG
	dbout << "Reading fluid file. Path:" << endl << file << endl;
#endif
	std::ifstream fluid(file);
	if (!fluid.is_open())
	{	
		cerr << "Error opening the fluid file." << endl;
	    exit(-1);
	}

	/*Fluid parameters read*/
	uint lineno = 0;
	Eigen::Vector2d nb = getvector(fluid, lineno, "Newmark-Beta terms");
	svar.beta = nb[0];	svar.gamma = nb[1];
	fvar.Hfac = getDouble(fluid, lineno, "Smoothing length factor"); /*End of state read*/

  	fvar.alpha = getDouble(fluid, lineno, "Artificial visc");
	fvar.contangb = getDouble(fluid, lineno, "contact angle");
	fvar.rho0 = getDouble(fluid, lineno, "Fluid density rho0");
	avar.rhog = getDouble(fluid, lineno, "Air density rhog");
	fvar.Cs = getDouble(fluid, lineno, "Speed of sound");
	fvar.mu = getDouble(fluid, lineno, "Fluid viscosity");
	avar.mug = getDouble(fluid, lineno, "Air viscosity");
	fvar.sig = getDouble(fluid, lineno, "Surface Tension");
	svar.outfolder = getString(fluid,lineno, "Output Folder name");
	svar.outdir = svar.outfolder;
	if(svar.Asource != 0)
	{
		svar.meshfile = getString(fluid, lineno, "Mesh input file");
  		svar.solfile = getString(fluid,lineno, "Mesh solution file");
  		svar.scale = getDouble(fluid,lineno, "Mesh scale");
  		avar.vRef = getDouble(fluid, lineno, "Gas ref Vel");
  		avar.pRef = getDouble(fluid, lineno, "Get ref Press");
  		avar.T = getDouble(fluid, lineno, "Gas ref Temp");
	}
		

#ifdef DEBUG
	dbout << "Closing fluid file"  << endl;
#endif
	fluid.close();

	if(svar.Asource == 3)
	{
		avar.vRef=avar.vInf.norm();
	}
}


void GetInput(int argc, char **argv, SIM& svar, FLUID& fvar, AERO& avar)
{
	svar.restart = 0;
	StateVecD angle;
	if (argc > 3) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only a maximum of two input arguments accepted,\n";
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
	    	
	if(file == "-r")
	{
		cout << "Restart option selected." << endl;
		file = argv[2];
    	svar.restart = 1; 
	}

	if(file.back() != '/')
    	file.append("/");

  	Read_SIM_Var(file,svar,fvar,avar);
	
	/*Get fluid properties from fluid.dat*/
	Read_FLUID_Var(file,svar,fvar,avar);
	

  	/*Universal parameters based on input values*/
	fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
	fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/

	/*Pipe Pressure calc*/
	real rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);

	svar.nrad = 1;
	if(svar.Bcase == 0 || svar.Bcase == 1)
	{
		svar.dx = svar.Pstep;
	}
	else if(svar.Bcase == 2 || svar.Bcase == 3)
	{
		// Defining dx to fit the pipe, then find rest spacing
		if(svar.Jet(0)/svar.Pstep < 1.5)
		{	/*spacing is too close to full size to use standard adjustment*/
			cout << "Warning: particle spacing if of the same order of magintude of the jet diameter." << endl;
			cout << "Consider a different size for accuracy." << endl;
		}
		else
		{
			svar.nrad = ceil(abs(0.5*svar.Jet(0)/svar.Pstep));
			svar.dx = 0.5*(svar.Jet(0))/real(svar.nrad);
		}	
 	}
 	else if (svar.Bcase == 4)
 	{
 		if(svar.Jet(0)/svar.Pstep < 1.5)
		{	/*spacing is too close to full size to use standard adjustment*/
 			svar.nrad = 1;
 			svar.dx = svar.Jet(0);
		}
		else
		{
	 		svar.nrad = ceil(abs(svar.Jet(0)*M_PI/svar.Pstep));
	 		svar.dx = svar.Jet(0)*M_PI/real(svar.nrad);
 		}
 	}
 	else
 	{
 		if(svar.Jet(0)/svar.Pstep < 1.5)
		{	/*spacing is too close to full size to use standard adjustment*/
 			svar.nrad = 1;
 			svar.dx = svar.Jet(0);
		}
		else
		{
	 		svar.nrad = ceil(abs(0.5*svar.Jet(0)/svar.Pstep));
	 		svar.dx = 0.5*(svar.Jet(0))/real(svar.nrad);
 		}
 	}


 	// svar.dx = svar.Pstep;	
	svar.Pstep = svar.dx * pow(rho/fvar.rho0,1.0/SIMDIM);

	// Correct the droplet to have the same volume as the original
	if(svar.Bcase == 4)
	{
		cout << "Droplet Diameter: " << svar.Jet(0) << endl;
		svar.diam = svar.Jet(0);
	}


#if SIMDIM == 3
 	avar.pVol = 4.0/3.0 * M_PI * pow(svar.Pstep,SIMDIM);
#else
 	avar.pVol = M_PI* svar.Pstep*svar.Pstep;
#endif

 	/*Mass from spacing and density*/
	fvar.simM = fvar.rho0*pow(svar.Pstep,SIMDIM); 
	fvar.bndM = fvar.simM;
	avar.gasM = avar.rhog*pow(svar.Pstep,SIMDIM);
	avar.qInf = 0.5*avar.rhog*avar.vRef*avar.vRef;

	svar.addcount = 0;
  	svar.dt = 2E-010; 			/*Initial timestep*/
  	svar.t = 0.0;				/*Total simulation time*/
  	fvar.H = fvar.Hfac*svar.Pstep;
  	fvar.HSQ = fvar.H*fvar.H; 

	fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/
	svar.Bclosed = 0; 		/*Boundary begins open*/
	svar.psnPts = 0; 		/*Start with no pitson points*/
  	svar.delNum = 0;
  	svar.intNum = 0;

  	fvar.dCont = 2.0 * fvar.delta * fvar.H * fvar.Cs;
  	fvar.dMom = fvar.alpha * fvar.H * fvar.Cs * fvar.rho0;

#if SIMDIM == 2
	fvar.correc = (7/(4*M_PI*fvar.H*fvar.H));
	svar.simPts = svar.xyPART[0]*svar.xyPART[1];
#endif
#if SIMDIM == 3
	fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
	svar.simPts = svar.xyPART[0]*svar.xyPART[1]*svar.xyPART[2]; /*total sim particles*/
#endif

#ifdef DEBUG
	dbout << "Tait Gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
	dbout << "Pipe rho: " << rho << endl;
	dbout << "Number of fluid particles along diameter: " << 2*svar.nrad+1 << endl;
	dbout << "Pipe step (dx): " << svar.dx << endl;
	dbout << "Particle mass: " << fvar.simM << endl;
	dbout << "Freestream initial spacing: " << svar.Pstep << endl;
	dbout << "Support Radius: " << fvar.H << endl;
	dbout << "Gas Dynamic pressure: " << avar.qInf << endl << endl;
#endif

	cout << "****** SIMULATION SETTINGS *******" << endl;
	cout << "Tait gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
	cout << "Newmark-Beta parameters: " << svar.beta << ", " << svar.gamma << endl;
	cout << "Boundary case: " << svar.Bcase << endl;
	cout << "Aerodynamic source: " << svar.Asource << endl;
	cout << "Aerodynamic case: " << avar.acase << endl;
	cout << "Speed of sound: " << fvar.Cs << endl;
	cout << endl;
	cout << "****** PIPE SETTINGS *******" << endl;
	cout << "Pipe pressure: " << fvar.pPress << endl;
	cout << "Pipe density: " << rho << endl;
	cout << "Pipe step (dx): " << svar.dx << endl;
	cout << "Pipe diameter: " << svar.Jet(0) << endl;
	cout << "Number of fluid particles along diameter: " << 2*svar.nrad+1 << endl;
	cout << "Pipe start position: " << svar.Start(0) << "  " << svar.Start(1);
	#if SIMDIM == 3
	cout << "  " << svar.Start(2);
	#endif
	cout << endl;
	cout << "Pipe start rotation: " << svar.Angle(0) << "  " << svar.Angle(1);
	#if SIMDIM == 3
	cout << "  " << svar.Angle(2);
	#endif
	cout << endl;
	cout << "Jet velocity: " << avar.vJet.norm() << endl;

	cout << endl;		
	cout << "****** RESTING FUEL SETTINGS *******" << endl;
	cout << "Particle spacing: " << svar.Pstep << endl;
	cout << "Support radius: " << fvar.H << endl;
	cout << "Particle mass: " << fvar.simM << endl;
	cout << "Liquid viscosity: " << fvar.mu << endl;
	cout << "Resting density: " << fvar.rho0 << endl;
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

	if(svar.Asource == 1 || svar.Asource == 2)
	{
		cout << "Mesh filename: " << svar.meshfile << endl;
		cout << "Solution filename: " << svar.solfile << endl;
		cout << "Reference velocity: " << avar.vRef << endl;
		cout << "Reference pressure: " << avar.pRef << endl;
		cout << "Reference temperature: " << avar.T << endl;
		cout << "Gas dynamic pressure: " << avar.qInf << endl << endl;
	}
	cout << "Gas density: " << avar.rhog << endl;
	cout << "Gas viscosity: " << avar.mug << endl << endl;
	

	cout << "******** FILE SETTINGS ********" << endl;
	cout << "Input folder: " << svar.infolder << endl;
	if(svar.Asource == 1 || svar.Asource == 2)
	{
		cout << "Input mesh face file: " << svar.meshfile  << endl;
		cout << "Input solution file: " << svar.solfile << endl;
	}

	cout << "Output folder: " << svar.outfolder << endl << endl;
	



	#if SIMDIM == 3
		if(svar.Asource == 3)
		{
			svar.vortex.Init(svar.infolder);	
			svar.vortex.GetGamma(avar.vInf);	
		}
	#endif

	
	GetYcoef(avar, fvar, /*fvar.H*/ svar.Pstep);
#if SIMDIM == 3
	avar.aPlate = svar.Pstep*svar.Pstep;
		// avar.aPlate = fvar.H*fvar.H;
#else
	avar.aPlate = svar.Pstep/**svar.Pstep*/ /** pow(avar.L,0.5)*/;
		// avar.aPlate = fvar.H;
#endif
	// if(svar.restart!=1)
	// 	Write_Input(svar,fvar,avar,angle);
}


void Restart(SIM& svar, FLUID& fvar, AERO& avar, State& pn, State& pnp1, MESH& cells)
{	
	// Read the values from the solution folder, then check. 

	// Check that a restart can be performed. 
	if(svar.outform == 0 || svar.outform == 1 || svar.outform == 4)
	{
		cout << "Output type cannot be restarted from. Make sure that you have the correct output selected." << endl;
		exit(-1);
	}

	if(svar.outform == 2)
	{
		if(svar.Bcase == 6)
		{
			cout << "No cell information. Cannot restart" << endl;
			exit(-1);
		}
	}
	// if(svar.boutform == 0)
	// {
	// 	cout << "No boundary time output, so can't be certain of forces from the boundary." << endl;
	// 	exit(-1);
	// }

	if(svar.boutform == 0 && (svar.Bcase != 4 && svar.Bcase !=0))
	{
		cout << "Time data for the boundary has not been output. Cannot restart." << endl;
		exit(-1);
	}

#ifdef DEBUG
	dbout << "Reading frame info file for restart information." << endl;
#endif

	// Re-read the settings and fluid file
	// cout << svar.outfolder << endl;

	string outdir = svar.outfolder;
	if(outdir.back() != '/')
    	outdir.append("/");


  	Read_SIM_Var(outdir,svar,fvar,avar);

  	Read_FLUID_Var(outdir,svar,fvar,avar);

  	svar.outfolder = outdir;

	// Now get the data from the files. Start with the boundary
	if(svar.outtype == 0)
	{
		State boundary, fuel;

		string file;
		INTEGER4 I;
		INTEGER4 frameNo;
		

		// Read the fuel
		void* fuelHandle = NULL;
		string fuelf = outdir;
		fuelf.append("Fuel.szplt");

		// cout << fuelf << endl;
		// fuelf = "Droplet2D/solution/Fuel.szplt";
		I = tecFileReaderOpen(fuelf.c_str(),&fuelHandle);
		if(I == -1)
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << file << endl;
			exit(-1);
		}

		cout << "Checking Fuel file..." << endl;
		CheckContents(fuelHandle,svar);

	    // Read the actual data.
	    frameNo = svar.frame+1;

		cout << endl << "Attempting to read the fuel..." << endl;
		Read_Binary_Timestep(fuelHandle,svar,frameNo,fuel);

		if (svar.Bcase != 4 && svar.Bcase !=0)
		{
			void* boundHandle = NULL;
			string boundf = outdir;
			boundf.append("Boundary.szplt");

			I = tecFileReaderOpen(boundf.c_str(),&boundHandle);
			if(I == -1)
			{
				cout << "Error opening szplt file. Path:" << endl;
				cout << file << endl;
				exit(-1);
			}

			cout << "Checking Boundary file..." << endl;
			CheckContents(boundHandle,svar);

			cout << "Attempting to read the boundary..." << endl;
			Read_Binary_Timestep(boundHandle,svar,frameNo,boundary);
		
		}

		pn = boundary;
		svar.bndPts = boundary.size();
		svar.simPts = fuel.size();
		pn.insert(pn.end(),fuel.begin(),fuel.end());	
		svar.totPts = pn.size();
		
		if(svar.simPts + svar.bndPts != svar.totPts)
		{
			cout << "Mismatch of array sizes. Total array is not the sum of the others" << endl;
			exit(-1);
		}

		// if(svar.totPts != totPts || svar.bndPts != bndPts || svar.simPts != simPts)
		// {
		// 	cout << "Mismatch of the particle numbers in the frame file and data file." << endl;
		// }
	}
	else
	{
		// TODO: ASCII Restart.
		// Particle numbers can be found from frame file.
		// Find EOF, then walk back from there how many particles.
	}

	// Go through the particles giving them the properties of the cell
	#pragma omp parallel for 
	for(size_t ii = 0; ii < svar.totPts; ++ii)
	{
		pn[ii].partID = ii;
		pn[ii].p =  fvar.B*(pow(pn[ii].rho/fvar.rho0,fvar.gam)-1);

		if((svar.Asource == 1 || svar.Asource == 2) && svar.outform == 5)
		{
			if(pn[ii].b == PartState.FREE_)
			{
				pn[ii].cellV = cells.cVel[pn[ii].cellID];
				pn[ii].cellP = cells.cP[pn[ii].cellID];
			}
		}

		if(pn[ii].b == PartState.BACK_)
		{
			svar.back.emplace_back(ii);
		}
	}
	
	pnp1 = pn;
	// Define svar.clear to state a particle is clear of the starting area
	if (svar.Bcase == 2 || svar.Bcase == 3 || svar.Bcase == 5)
		svar.clear = -svar.Jet[1] + 4*svar.dx;
	else
		svar.clear = 0.0;

	// int width = 15;

	// for(size_t ii = 0; ii < svar.totPts; ++ii)
	// {
	// 	cout << setw(5) << pnp1[ii].partID << setw(width) << pnp1[ii].xi(0) << setw(width) << pnp1[ii].xi(1) << setw(width) <<
	// 	pnp1[ii].rho << setw(width) << pnp1[ii].Rrho << setw(width) << pnp1[ii].m << setw(width) << pnp1[ii].v(0) << 
	// 	setw(width) << pnp1[ii].v(1) << setw(width) << pnp1[ii].f(0) << setw(width) << pnp1[ii].f(1) << setw(4) << 
	// 	pnp1[ii].b << setw(width) << pnp1[ii].cellV(0) << setw(width) << pnp1[ii].cellV(1) << setw(width) << pnp1[ii].cellP
	// 	<< setw(width) << pnp1[ii].cellID << endl;
	// }
}


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


void Write_ASCII_Timestep(std::fstream& fp, SIM &svar, const State &pnp1, 
	const uint bwrite, const uint start, const uint end, const string& name)
{
	if(bwrite == 1)
	 	fp <<  "ZONE T=\"" << name << "\"";
	else
		fp <<  "ZONE T=\"" << name << "\"";

    fp <<", I=" << end - start << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
    fp << std::left << std::scientific << std::setprecision(6);
    const static uint width = 15;
    
    if(svar.outform == 0)
    {
		for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
		{
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale; 
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
	}	
    if(svar.outform == 1)
    {
		for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
		{	
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;

			fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;
	        fp << setw(width) << p->v.norm();
	        fp << setw(width) << p->f.norm() << "\n";
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
		 //        fp << p.f.norm() << " ";
		 //        fp << p.rho << " "  << p.p  << "\n";
		 //  	}
	 	// }
	  	fp << std::flush; 	
    }
    else if (svar.outform == 2)
	{
		for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
		{
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;
	        
	        fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->v(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
				fp << setw(width) << p->f(i); 

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

		  //     	fp << p.f.norm() << " " << p.Af.norm() << " " << p.Sf.norm() << " ";
		  //       for(uint i = 0; i < SIMDIM; ++i)
		  //       	fp << p.cellV(i) << " "; 

		  //       fp << p.b << " " << p.theta  << "\n"; 
		  // 	}
	  	// }	
	  	fp << std::flush;
    }
    else if (svar.outform == 3)
    {
    	for (auto p=std::next(pnp1.begin(),svar.bndPts); p!=std::next(pnp1.begin(),svar.bndPts+svar.simPts); ++p)
		{
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;
	        
	        fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->v(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->f(i); 

	        for(uint i = 0; i < SIMDIM; ++i)
				fp << setw(width) << p->cellV(i);

			fp << setw(width) << p->cellP;
	        fp << setw(width) << p->cellID << "\n"; 
	  	}
    }
    else if (svar.outform == 4)
    {
    	for (auto p=std::next(pnp1.begin(),start); p!=std::next(pnp1.begin(),end); ++p)
		{	
			for(uint i = 0; i < SIMDIM; ++i)
	        	fp << setw(width) << p->xi(i)/svar.scale;

			fp << setw(width) << p->rho << setw(width) << p->p;
			fp << setw(width) << p->m;
	        fp << setw(width) << p->v.norm();
	        fp << setw(width) << p->f.norm();
	        fp << setw(width) << p->b << setw(width) << p->theta;
	        fp << setw(width) << p->Af.norm() << "\n";
	  	}

    }
}

void Write_First_Step(std::fstream& f1, std::fstream& fb, std::fstream& fg, 
						SIM& svar, State const& pnp1, State const& airP)
{
	if(svar.outtype == 0 )
	{
		/*Write sim particles*/
		
		Init_Binary_PLT(svar,"Fuel.szplt","Simulation Particles",svar.fuelFile);
		if(svar.restart == 0)
		{
			Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1,svar.fuelFile);
		}

		if (svar.Bcase != 0 && svar.Bcase !=4)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"Boundary.szplt","Boundary Particles",svar.boundFile);

			if(svar.boutform == 0 && svar.restart == 0)
			{   //Don't write a strand, so zone is static. 
				Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",0,svar.boundFile); 
				int32_t i = tecFileWriterClose(&svar.boundFile);			
				if(i == -1)
					exit(-1);
			}
			else if(svar.restart == 0)
				Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2,svar.boundFile); 
		}	

		if (svar.ghost == 1 && svar.gout == 1)
		{
			/*Write boundary particles*/
			Init_Binary_PLT(svar,"Ghost.szplt","Ghost Particles",svar.ghostFile);		
		}	
	}
	else if (svar.outtype == 1)
	{

		if (svar.Bcase != 0 && svar.Bcase != 4 &&  svar.restart == 0)
		{	/*If the boundary exists, write it.*/
			string bfile = svar.outfolder;
			bfile.append("Boundary.plt");
			fb.open(bfile, std::ios::out);
			if(fb.is_open())
			{
				Write_ASCII_header(fb,svar);
				Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
				if(svar.boutform == 0)
					fb.close();
			}
			else
			{
				cerr << "Error opening boundary file." << endl;
				exit(-1);
			}
		}

		/* Write first timestep */
		string mainfile = svar.outfolder;
		mainfile.append("Fuel.plt");
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
			string ghostfile = svar.outfolder;
			ghostfile.append("Ghost.plt");
			fg.open(ghostfile,std::ios::out);
			if(fg.is_open())
			{
				Write_ASCII_header(fg,svar);
				Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
			}
		}
	}
	else
	{
		cerr << "Output type ambiguous. Please select 0 or 1 for output data type." << endl;
		exit(-1);
	}
}

void Write_Timestep(std::fstream& f1, std::fstream& fb, std::fstream& fg, uint ghost_strand, 
				SIM& svar, State const& pnp1, State const& airP)
{
	if (svar.outtype == 0)
	{
		if(svar.Bcase != 0 && svar.Bcase != 4 && svar.boutform == 1)
		{	/*Write boundary particles*/
			Write_Binary_Timestep(svar,pnp1,0,svar.bndPts,"Boundary",2,svar.boundFile); 
		}
		Write_Binary_Timestep(svar,pnp1,svar.bndPts,svar.totPts,"Fuel",1,svar.fuelFile); /*Write sim particles*/
		if(svar.ghost == 1 && svar.gout == 1 && airP.size() != 0)
			Write_Binary_Timestep(svar,airP,0,airP.size(),"Ghost",ghost_strand,svar.ghostFile);
	} 
	else if (svar.outtype == 1)
	{
		Write_ASCII_Timestep(f1,svar,pnp1,0,svar.bndPts,svar.totPts,"Fuel");
		if(svar.Bcase != 0 && svar.Bcase != 4 && svar.boutform == 1)
		{
			State empty;
			Write_ASCII_Timestep(fb,svar,pnp1,1,0,svar.bndPts,"Boundary");
		}

		if(svar.ghost == 1 && svar.gout == 1)
			Write_ASCII_Timestep(fg,svar,airP,0,0,airP.size(),"Ghost");
	}
}


#endif