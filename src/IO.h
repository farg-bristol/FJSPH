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
void Restart(SIM& svar, State& pn, MESH& cells)
{	
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

	if(svar.boutform == 0)
	{
		cout << "Time data for the boundary has not been output. Cannot restart." << endl;
		exit(-1);
	}

#ifdef DEBUG
	dbout << "Reading frame info file for restart information." << endl;
#endif

	// find the last timestep by reading the frame file
	// string file = svar.outfolder;
// 	file.append("frame.info");

// 	// Find end of file
// 	std::ifstream f1(file,std::ios::in);
// 	if(!f1.is_open())
// 	{
// 		cout << "Failed to open frame.info file" << endl;
// 		exit(-1);
// 	}

// 	uint lsize = 0;
// 	while(f1.ignore(std::numeric_limits<std::streamsize>::max(),'\n'))	
// 		lsize++;

// 	// cout << "Found file size: " << lsize 
// 	// 	 << ". Now getting last frame." << endl;
// 	f1.clear();
// 	GotoLine(f1, lsize-4);

// 	string line;
// 	getline(f1,line);
// 	std::stringstream sstr(line);
// 	string temp;
// 	uint frame;
// 	while(!sstr.eof())
// 	{
// 		sstr >> temp;
// 		if(std::stringstream(temp) >> frame)
// 		{
// 			break;
// 		}
// 		temp = "";
// 	}

// #ifdef DEBUG
// 	dbout << "Final frame: " << frame << endl;
// #endif
// 	cout << "Found the final frame: " << frame << endl;
// 	svar.frame = frame;

// #ifdef DEBUG
// 	dbout << "Getting particle numbers" << endl;
// #endif

// 	// Get the particle counts to checksum.
// 	uint bndPts, simPts, totPts;
	
// 	getline(f1,line);
// 	std::stringstream sline(line);
// 	vector<uint> pts;
// 	temp = "";
// 	int var;
// 	while(!sline.eof())
// 	{
// 		sline >> temp;
// 		if(std::stringstream(temp) >> var)
// 		{
// 			pts.emplace_back(var);
// 		}
// 		temp = "";
// 	}

// 	totPts = pts[0];
// 	bndPts = pts[1];
// 	simPts = pts[2];
// 	cout << line << endl;
// #ifdef DEBUG
// 	dbout << line << endl;
// #endif

	// Now get the data from the files. Start with the boundary
	if(svar.outtype == 0)
	{
		string file;
		INTEGER4 I;
		INTEGER4 frameNo;
		void* boundHandle = NULL;
		string boundf = svar.outfolder;
		boundf.append("Boundary.szplt");

		// Read the fuel
		void* fuelHandle = NULL;
		string fuelf = svar.outfolder;
		fuelf.append("Fuel.szplt");

		
		I = tecFileReaderOpen(boundf.c_str(),&boundHandle);
		if(I == -1)
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << file << endl;
			exit(-1);
		}

		I = tecFileReaderOpen(fuelf.c_str(),&fuelHandle);
		if(I == -1)
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << file << endl;
			exit(-1);
		}

		cout << "Checking Boundary file..." << endl;
		CheckContents(boundHandle,svar);
		cout << "Checking Fuel file..." << endl;
		CheckContents(fuelHandle,svar);

	    // Read the actual data.
	    
	    State boundary, fuel;
	    frameNo = svar.frame+1;

	    cout << "Attempting to read the boundary..." << endl;
		Read_Binary_Timestep(boundHandle,svar,frameNo,boundary);
		
		cout << endl << "Attempting to read the fuel..." << endl;
		Read_Binary_Timestep(fuelHandle,svar,frameNo,fuel);

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
		if(svar.Bcase == 6 && svar.outform == 5)
		{
			if(pn[ii].b == FREE)
			{
				pn[ii].cellV = cells.cVel[pn[ii].cellID];
				pn[ii].cellP = cells.cP[pn[ii].cellID];
			}
		}

		if(pn[ii].b == BACK)
		{
			svar.back.emplace_back(ii);
		}
	}
}

void GetInput(int argc, char **argv, SIM& svar, FLUID& fvar, AERO& avar)
{
	svar.restart = 0;
	StateVecD angle;
	if (argc > 3) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only two input arguments accepted,\n";
		cout << "1: Input file directoy   2: Output file directory.\n";
		cout << "Other inputs will be ignored." << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tERROR: No inputs provided. Stopping... \n";
    	exit(-1);    	
    }
    else if (argc > 1)
    {	/*Get parameters if it has been provided*/
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

	    svar.infolder = file;
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
  		svar.Bcase = getInt(in, lineno, "Simulation boundary case");
  		avar.acase = getInt(in, lineno, "Simulation Aerodynamic case");
  		svar.ghost = getInt(in, lineno, "Ghost particles?");
  		svar.Start = getDVector(in, lineno, "Starting position");
  		if(svar.Bcase < 2)
  		{
  			svar.xyPART = getIVector(in, lineno, "Particles in each coordinate");
  			svar.Box= getDVector(in, lineno, "Box dimensions");
  			(void)getDouble(in, lineno, "Skipped");
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
  			angle = angles;
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
	  		if(avar.acase > 5)
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

	/*Get fluid properties from fluid.dat*/
	svar.scale = 1.0;
	string file = svar.infolder;
	file.append("Fluid");
#ifdef DEBUG
	dbout << "Reading fluid file. Path:" << endl << file << endl;
#endif
	std::ifstream fluid(file);
	if (fluid.is_open())
	{	/*Fluid parameters read*/
		uint lineno = 0;
		Eigen::Vector2d nb = getvector(fluid, lineno, "Newmark-Beta terms");
		svar.beta = nb[0];	svar.gamma = nb[1];
		fvar.Hfac = getDouble(fluid, lineno, "Smoothing length factor"); /*End of state read*/

	  	fvar.alpha = getDouble(fluid, lineno, "Artificial visc");
  		fvar.contangb = getDouble(fluid, lineno, "contact angle");
  		fvar.rho0 = getDouble(fluid, lineno, "Fluid density rho0");
  		fvar.rhog = getDouble(fluid, lineno, "Air density rhog");
  		fvar.Cs = getDouble(fluid, lineno, "Speed of sound");
  		fvar.mu = getDouble(fluid, lineno, "Fluid viscosity");
  		fvar.mug = getDouble(fluid, lineno, "Air viscosity");
  		fvar.sig = getDouble(fluid, lineno, "Surface Tension");
  		if(svar.Bcase == 6 || svar.Bcase == 4 || svar.Bcase == 7)
  		{
	  		fvar.gasVel = getDouble(fluid, lineno, "Gas ref Vel");
	  		fvar.gasPress = getDouble(fluid, lineno, "Get ref Press");
	  		fvar.T = getDouble(fluid, lineno, "Gas ref Temp");

	  		if(svar.Bcase == 6)
	  		{
		  		svar.meshfile = getString(fluid, lineno, "Mesh input file");
		  		svar.solfile = getString(fluid,lineno, "Mesh solution file");
		  		svar.scale = getDouble(fluid,lineno, "Mesh scale");
	  		}
  		}
  		svar.outfolder = getString(fluid,lineno, "Output Folder name");

#ifdef DEBUG
		dbout << "Closing fluid file"  << endl;
#endif
  		fluid.close();

  		if(svar.Bcase == 4)
  		{
  			fvar.gasVel=avar.vInf(1);
  		}
	}
	else 
	{
		cerr << "Fluid file not found. Assuming standard set of parameters (water)." << endl;
		svar.beta = 0.25;	svar.gamma = 0.5;
		fvar.H= 2.0*svar.Pstep;
	  	fvar.alpha = 0.1;
  		fvar.contangb = 150.0;
  		fvar.rho0 = 1000.0;
  		fvar.Cs = 100.0;
  		fvar.mu = 1.0;
  		fvar.sig = 0.0728;
  		fvar.rhog = 1.225;
  		fvar.mug = 18.5E-06;
	}

	  	/*Universal parameters based on input values*/
		fvar.gam = 7.0;  							 /*Factor for Tait's Eq*/
		fvar.B = fvar.rho0*pow(fvar.Cs,2)/fvar.gam;  /*Factor for Tait's Eq*/

		/*Pipe Pressure calc*/
		real rho = fvar.rho0*pow((fvar.pPress/fvar.B) + 1.0, 1.0/fvar.gam);

		svar.nrad = 0;

		if(svar.Bcase == 3 || svar.Bcase == 4 || svar.Bcase == 6)
		{
			// Defining dx to fit the pipe, then find rest spacing
			svar.nrad = ceil(abs(0.5*svar.Jet(0)/svar.Pstep));
			svar.dx = 0.5*(svar.Jet(0))/real(svar.nrad);
	 	}
	 	else
	 	{
	 		svar.nrad = ceil(abs(0.5*svar.Start(0)/svar.Pstep));
	 		svar.dx = 0.5*(svar.Start(0))/real(svar.nrad);
	 	}

	 // 	svar.dx = svar.Pstep;	
		svar.Pstep = svar.dx * pow(rho/fvar.rho0,1.0/SIMDIM);

		// Correct the droplet to have the same volume as the original
		if(svar.Bcase == 5)
		{
			int partCount = ParticleCount(svar);
			cout << "Predicted particle count: " << partCount << endl;

			
#if SIMDIM == 3
			real dvol = (4.0*M_PI/3.0)*pow(0.5*svar.Start(0),3);
#else
			real dvol = M_PI*pow(0.5*svar.Start(0),2);
#endif

			// Simulation mass
			// real volume = pow(svar.dx,SIMDIM)*partCount;
			// svar.mass = volume * rho;

			// cout << "SPH Volume: " << volume << "  Droplet expected volume: " << dvol << endl;
			// cout << "SPH Mass: " << svar.mass << "  Droplet expected mass: " << dvol*rho << endl;	

			// Want to match the volume, so adjust step size to match. 

			real newstep = 0.98*pow(dvol/real(partCount),1.0/SIMDIM);
			// cout << "Old step: " << svar.dx << " New step: " << newstep << endl;
			svar.dx = newstep;
			// real radius = svar.dx * real(svar.nrad);

			// cout << "Old diameter: " << svar.Start(0) << "  New diameter: " << radius*2.0 << endl;
			svar.diam = svar.Start(0);
			svar.Start(0) = 2.0*newstep*real(svar.nrad);

			// svar.nrad = ceil(abs(0.5*svar.Start(0)/svar.Pstep));
	 		// svar.dx = 0.5*(svar.Start(0))/real(svar.nrad);

			svar.Pstep = svar.dx * pow(rho/fvar.rho0,1.0/SIMDIM);

			// int newCount = ParticleCount(svar);

			// cout << "New particle count: " << newCount << endl;
		}


#if SIMDIM == 3
	 	avar.pVol = 4.0/3.0 * M_PI * pow(svar.Pstep,SIMDIM);
#else
	 	avar.pVol = M_PI* svar.Pstep*svar.Pstep;
#endif

	 	/*Mass from spacing and density*/
		fvar.simM = fvar.rho0*pow(svar.Pstep,SIMDIM); 
		fvar.bndM = fvar.simM;
		fvar.gasM = fvar.rhog*pow(svar.Pstep,SIMDIM);
		fvar.gasDynamic = 0.5*fvar.rhog*fvar.gasVel;

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
		dbout << "Number of fluid particles along diameter: " << 2*nrad+1 << endl;
		dbout << "Pipe step (dx): " << svar.dx << endl;
		dbout << "Particle mass: " << fvar.simM << endl;
		dbout << "Freestream initial spacing: " << svar.Pstep << endl;
		dbout << "Support Radius: " << fvar.H << endl;
		dbout << "Gas Dynamic pressure: " << fvar.gasDynamic << endl << endl;
#endif

		cout << "****** SIMULATION SETTINGS *******" << endl;
		cout << "Tait gamma: " << fvar.gam << "  Tait B: " << fvar.B << endl;
		cout << "Newmark-Beta parameters: " << svar.beta << ", " << svar.gamma << endl;
		cout << "Boundary case: " << svar.Bcase << endl;
		cout << "Aerodynamic case: " << avar.acase << endl;
		cout << "Speed of sound: " << fvar.Cs << endl;
		cout << endl;
		cout << "****** PIPE SETTINGS *******" << endl;
		cout << "Pipe density: " << rho << endl;
		cout << "Pipe step (dx): " << svar.dx << endl;
		cout << "Pipe diameter: " << svar.Jet(0) << endl;
		cout << "Number of fluid particles along diameter: " << 2*svar.nrad+1 << endl;
		cout << "Pipe start position: " << svar.Start(0) << "  " << svar.Start(1);
		#if SIMDIM == 3
		cout << "  " << svar.Start(2);
		#endif
		cout << endl;
		cout << "Pipe start rotation: " << angle(0) << "  " << angle(1);
		#if SIMDIM == 3
		cout << "  " << angle(2);
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
		if(svar.Bcase != 6)
		{
			cout << "Gas Velocity: " << avar.vInf(0) << "  " << avar.vInf(1);
#if SIMDIM == 3
			cout << "  " << avar.vInf(2);
#endif
			cout << endl;
		}

		if(svar.Bcase == 6)
		{
			cout << "Mesh filename: " << svar.meshfile << endl;
			cout << "Solution filename: " << svar.solfile << endl;
			cout << "Reference velocity: " << fvar.gasVel << endl;
			cout << "Reference pressure: " << fvar.gasPress << endl;
			cout << "Reference temperature: " << fvar.T << endl;
			cout << "Gas dynamic pressure: " << fvar.gasDynamic << endl << endl;
		}
		cout << "Gas density: " << fvar.rhog << endl;
		cout << "Gas viscosity: " << fvar.mug << endl;
		

		#if SIMDIM == 3
			if(svar.Bcase == 4)
			{
				svar.vortex.Init(svar.infolder);	
				svar.vortex.GetGamma(avar.vInf);	
			}
		#endif


		GetAero(avar, fvar, /*fvar.H*/svar.Pstep);
#if SIMDIM == 3
		avar.aPlate = svar.Pstep*svar.Pstep;
#else
		avar.aPlate = svar.Pstep;
#endif
		// if(svar.restart!=1)
		// 	Write_Input(svar,fvar,avar,angle);
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




#endif