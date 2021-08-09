#ifndef RESTART_H
#define RESTART_H
	

#include "Var.h"
#include "BinaryIO.h"
#include "IOFunctions.h"
#include "Geometry.h"

/*************************************************************************/
/************************** RESTART OUTPUTS ******************************/
/*************************************************************************/
void Write_Input_ASCII(SIM const& svar, FLUID const& fvar, AERO const& avar)
{
	// Open the file.
	string file = svar.output_prefix;
	file.append("_restart");
	uint width = 60;
	std::ofstream restf(file);
	restf << std::scientific << std::setprecision(15);
	restf << svar.framet << setw(width) << "#Frame time interval" << endl;
	restf << svar.Nframe << setw(width) << "#Number of frames" << endl;
	restf << svar.outframe << setw(width) << "#Output frame info" << endl;
	restf << svar.out_encoding << setw(width) << "#Output data type" << endl;
	restf << svar.outform << setw(width) << "#Output contents" << endl;
	restf << svar.gout << setw(width) << "#Ghost particle output" << endl;
	restf << svar.subits << setw(width) << "#Maximum sub iterations" << endl;
	restf << svar.finPts << setw(width) << "#Maximum number of particles" << endl;
	restf << svar.cellSize << setw(width) << "#Post processing mesh size" << endl;
	restf << svar.postRadius << setw(width) << "#Post processing support radius" << endl;
	restf << svar.dx << setw(width) << "#Particle actual spacing" << endl;
	restf << svar.Pstep << setw(width) << "#Particle resting spacing" << endl;
	restf << svar.Bstep << setw(width) << "#Boundary spacing factor" << endl;
	restf << svar.Bcase << setw(width) << "#Boundary case" << endl;
	restf << svar.Asource << setw(width) << "#Aerodynamic source" << endl; 
	restf << avar.acase << setw(width) << "#Aerodynamic case" << endl;
	restf << svar.ghost << setw(width) << "#Ghost particles?" << endl;
	restf << svar.Start(0) << " " << svar.sim_start(1);
	#if SIMDIM == 3
		restf << " " << svar.Start(2) << setw(28);
	#else
		restf << setw(44);
	#endif
	restf << "#Fluid start position" << endl;

	if(svar.Bcase < 2)
	{
		restf << svar.xyPART(0) << " " << svar.xyPART(1);
		#if SIMDIM == 3
			restf << " " << svar.xyPART(2) << setw(28);
		#else
			restf << setw(44);
		#endif
		restf << "#Particle counts" << endl;

		restf << svar.Box(0) << " " << svar.Box(1);
		#if SIMDIM == 3
			restf << " " << svar.Box(2);
		#endif
		restf << setw(width) << "#Box dimensions" << endl << endl;
		restf << fvar.pPress << setw(width) << "#Pipe Pressure" << endl;
	}
	else
	{
		restf << svar.Angle(0) << " " << svar.Angle(1);
		#if SIMDIM == 3
			restf << " " << svar.Angle(2) << setw(28);
		#else
			restf << setw(44);
		#endif
		restf << "#Fluid start rotation" << endl;

		restf << svar.jet_diam << " " << svar.jet_depth << setw(width) << "#Jet dimensions" << endl;
		restf << fvar.pPress << setw(width) << "#Pipe Pressure" << endl;

		restf << vJetMag << setw(width) << "#Jet velocity" << endl;
		restf << avar.vInf(0) << " " << avar.vInf(1);
		#if SIMDIM == 3
			restf << " " << avar.vInf(2)<< setw(28);
		#else
			restf << setw(44);
		#endif
		restf << "#Freestream velocity" << endl;

		// if(avar.acase == 2 || avar.acase == 3)
  		// {
  		// 	restf << avar.a << setw(width) << "#a" << endl;
  		// 	restf << avar.h1 << setw(width) << "h1" << endl;
  		// 	restf << avar.b << setw(width) << "#b" << endl;
  		// 	restf << avar.h2 << setw(width) << "#h2" << endl;
  		// }
	}

	restf << endl;

	restf << fvar.alpha << setw(width) << "#Artificial viscosity" << endl;
	restf << fvar.contangb << setw(width) << "#Contact angle" << endl;
	restf << fvar.rho0 << setw(width) << "#Fluid density" << endl;
	restf << avar.rhog << setw(width) << "#Gas density" << endl;
	restf << fvar.Cs << setw(width) << "#Speed of sound" << endl;
	restf << fvar.mu << setw(width) << "#Fluid viscosity" << endl;
	restf << avar.mug << setw(width) << "#Gas viscosity" << endl;
	restf << fvar.sig << setw(width) << "#Surface Tension" << endl;
	restf << svar.outdir << setw(width) << "#Output Folder" << endl;
	if(svar.Asource > 0)
	{
		restf << svar.meshfile << "   #Mesh edge/face file" << endl;
		restf << svar.solfile << "   #Mesh Solution file" << endl;
		restf << svar.scale << setw(width) << "#Mesh scale" << endl;
		restf << avar.vRef << setw(width) << "#Gas reference velocity" << endl;
		restf << avar.pRef << setw(width) << "#Gas reference pressure" << endl;
		restf << avar.T << setw(width) << "#Gas reference temperature" << endl;
	}
	

	/*End of restfings write*/
	restf.close();
}


void Write_Input_TECIO(SIM const& svar, FLUID const& fvar, AERO const& avar)
{
	int32_t FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
    int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT

    string file = svar.output_prefix;
	file.append("_restart.szplt");

	vector<int32_t> varTypes;
	#if FOD == 1
	int32_t realType = 2;
	#else
		int32_t realType = 1;
	#endif

	string var = "";
	var.append("framet,"); varTypes.emplace_back(realType);
	var.append("Nframe,"); varTypes.emplace_back(3);
	var.append("outframe,"); varTypes.emplace_back(5);
	var.append("outtype,"); varTypes.emplace_back(5);
	var.append("outform,"); varTypes.emplace_back(5);
	var.append("gout,"); varTypes.emplace_back(5);
	var.append("subits,"); varTypes.emplace_back(5);
	var.append("nmax,"); varTypes.emplace_back(3);
	var.append("cellSize,"); varTypes.emplace_back(realType);
	var.append("PostRadius,"); varTypes.emplace_back(realType);
	var.append("dx,"); varTypes.emplace_back(realType);
	var.append("Pstep,"); varTypes.emplace_back(realType);
	var.append("Bstep,"); varTypes.emplace_back(realType);
	var.append("Bcase,"); varTypes.emplace_back(3);
	var.append("Asource,"); varTypes.emplace_back(3);
	var.append("acase,"); varTypes.emplace_back(3);
	var.append("ghost,"); varTypes.emplace_back(3);
	var.append("start,"); varTypes.emplace_back(realType);
	if(svar.Bcase < 2)
	{
		var.append("xyPART,"); varTypes.emplace_back(3);
		var.append("Box,"); varTypes.emplace_back(realType);
		var.append("pPress,"); varTypes.emplace_back(realType);
	}
	else
	{
		var.append("nrad,"); varTypes.emplace_back(3);
		var.append("Angle,"); varTypes.emplace_back(realType);
		var.append("Jet,"); varTypes.emplace_back(realType);
		var.append("pPress,"); varTypes.emplace_back(realType);
		var.append("vJet,"); varTypes.emplace_back(realType);
		var.append("vInf,"); varTypes.emplace_back(realType);

		if(avar.acase == 2 || avar.acase == 3)
  		{
			var.append("a,"); varTypes.emplace_back(realType);
			var.append("h1,"); varTypes.emplace_back(realType);
			var.append("b,"); varTypes.emplace_back(realType);
			var.append("h2,"); varTypes.emplace_back(realType);
		}
	}

	var.append("alpha,"); varTypes.emplace_back(realType);
	var.append("contangb,"); varTypes.emplace_back(realType);
	var.append("rho0,"); varTypes.emplace_back(realType);
	var.append("rhog,"); varTypes.emplace_back(realType);
	var.append("Cs,"); varTypes.emplace_back(realType);
	var.append("mu,"); varTypes.emplace_back(realType);
	var.append("mug,"); varTypes.emplace_back(realType);
	var.append("sig"); varTypes.emplace_back(realType);
	if(svar.Asource > 0)
	{
		var.append(",scale,"); varTypes.emplace_back(realType);
		var.append("vRef,"); varTypes.emplace_back(realType);
		var.append("pRef,"); varTypes.emplace_back(realType);
		var.append("T"); varTypes.emplace_back(realType);
	}


	string zone = "Restart File";
	void* fileHandle = NULL;
	if(tecFileWriterOpen(file.c_str(),zone.c_str(),var.c_str(),fileFormat,FileType,1,NULL,&fileHandle))
    {
    	cout << "Failed to open " << file << endl;
    	exit(-1);
    }

#ifdef DEBUG
	if(tecFileSetDiagnosticsLevel(fileHandle, 1))
	{
		cerr << "Failed to set debug option for output file: " << file << endl;
		exit(-1);
	}
#endif

    int32_t frame;
    zone = "Restart data";
    int64_t size = SIMDIM;
	// Get zone data types
	vector<int32_t> shareVarFromZone(varTypes.size(),0);
    vector<int32_t> valueLocation(varTypes.size(),1);
    vector<int32_t> passiveVarList(varTypes.size(),0);

	if(tecZoneCreateIJK(fileHandle,zone.c_str(),size,1,1,&varTypes[0],
		&shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&frame))
	{
		cerr << "Failed to create IJK zone." << endl;
		exit(-1);
	}
	
	int32_t varnum = 1;
	vector<uint8_t> uintVar(size,0);
	vector<int32_t> intVar(size,0);

	Write_Real_Value(fileHandle, frame, varnum, size, svar.framet, "frame time interval");
	Write_Int_Value(fileHandle, frame, varnum, size, svar.Nframe, "number of frames");
	Write_UInt_Value(fileHandle, frame, varnum, size, svar.outframe, "output frame info");
	Write_UInt_Value(fileHandle, frame, varnum, size, svar.outtype, "output data type");
	Write_UInt_Value(fileHandle, frame, varnum, size, svar.outform, "output content");
	Write_UInt_Value(fileHandle, frame, varnum, size, svar.gout, "ghost output");
	Write_UInt_Value(fileHandle, frame, varnum, size, svar.subits, "sub iterations");
	Write_Int_Value(fileHandle, frame, varnum, size, svar.finPts, "max particles");
	Write_Real_Value(fileHandle, frame, varnum, size, svar.cellSize, "cell size");
	Write_Real_Value(fileHandle, frame, varnum, size, svar.postRadius, "post process support radius");
	Write_Real_Value(fileHandle, frame, varnum, size, svar.dx, "dx");
	Write_Real_Value(fileHandle, frame, varnum, size, svar.Pstep, "Pstep");
	Write_Real_Value(fileHandle, frame, varnum, size, svar.Bstep, "Bstep");
	Write_Int_Value(fileHandle, frame, varnum, size, svar.Bcase, "boundary case");
	Write_Int_Value(fileHandle, frame, varnum, size, svar.Asource, "aerodynamic source");
	Write_Int_Value(fileHandle, frame, varnum, size, avar.acase, "aerodynamic case");
	Write_Int_Value(fileHandle, frame, varnum, size, svar.ghost, "ghost particle option");

	vector<real> vecreal(SIMDIM);
	for(size_t dim = 0; dim < SIMDIM; dim++)
	{
		vecreal[dim] = svar.Start(dim);
	}

	Write_Real_Vector(fileHandle, frame, varnum, size, vecreal, "start coordinates");


	if(svar.Bcase < 2)
	{
		vector<int32_t> vecint(SIMDIM);
		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			vecint[dim] = static_cast<int32_t>(svar.xyPART(dim));
		}
		
		Write_Int_Vector(fileHandle, frame, varnum, size, vecint, "particle dimension counts");

		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			vecreal[dim] = svar.Box(dim);
		}

		Write_Real_Vector(fileHandle, frame, varnum, size, vecreal, "box dimensions");

		Write_Real_Value(fileHandle, frame, varnum, size, fvar.pPress, "starting pressure");
	}
	else
	{
		Write_Int_Value(fileHandle, frame, varnum, size, svar.nrad, "diameter particle number");

		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			vecreal[dim] = svar.Angle(dim);
		}

		Write_Real_Vector(fileHandle, frame, varnum, size, vecreal, "starting angles");


		vector<real> vec2(3,0.0);
		for(size_t dim = 0; dim < 2; dim++)
		{
			vec2[dim] = svar.Jet(dim);
		}

		Write_Real_Vector(fileHandle, frame, varnum, size, vec2, "jet dimensions");

		Write_Real_Value(fileHandle, frame, varnum, size, fvar.pPress, "starting pressure");


		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			vecreal[dim] = avar.vJet(dim);
		}

		Write_Real_Vector(fileHandle, frame, varnum, size, vecreal, "jet velocity");

		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			vecreal[dim] = avar.vInf(dim);
		}
		
		Write_Real_Vector(fileHandle, frame, varnum, size, vecreal, "freestream velocity");

		if(avar.acase == 2 || avar.acase == 3)
  		{
			Write_Real_Value(fileHandle, frame, varnum, size, avar.a, "parameter a");
			Write_Real_Value(fileHandle, frame, varnum, size, avar.h1, "parameter h1");
			Write_Real_Value(fileHandle, frame, varnum, size, avar.b, "parameter b");
			Write_Real_Value(fileHandle, frame, varnum, size, avar.h2, "parameter h2");
  		}
	}

	Write_Real_Value(fileHandle, frame, varnum, size, fvar.alpha, "artificial dissipation factor");
	Write_Real_Value(fileHandle, frame, varnum, size, fvar.contangb, "contact angle");
	Write_Real_Value(fileHandle, frame, varnum, size, fvar.rho0, "liquid density");
	Write_Real_Value(fileHandle, frame, varnum, size, avar.rhog, "gas density");
	Write_Real_Value(fileHandle, frame, varnum, size, fvar.Cs, "speed of sound");
	Write_Real_Value(fileHandle, frame, varnum, size, fvar.mu, "liquid viscosity");
	Write_Real_Value(fileHandle, frame, varnum, size, avar.mug, "gas viscosity");
	Write_Real_Value(fileHandle, frame, varnum, size, fvar.sig, "surface tension");
	
	if(svar.Asource > 0)
	{
		Write_Real_Value(fileHandle, frame, varnum, size, svar.scale, "scale");
		Write_Real_Value(fileHandle, frame, varnum, size, avar.vRef, "reference velocity");
		Write_Real_Value(fileHandle, frame, varnum, size, avar.pRef, "reference pressure");
		Write_Real_Value(fileHandle, frame, varnum, size, avar.T, "reference temperature");
	}

	if(tecFileWriterClose(&fileHandle))
	{
		cerr << "Failed to close file." << endl;
	}
}

/*************************************************************************/
/************************** RESTART INPUTS *******************************/
/*************************************************************************/
void Read_Restart(string& infolder, SIM& svar, FLUID& fvar, AERO& avar)
{
	string file = infolder;
	file.append("Restart");
#ifdef DEBUG
	dbout << "Reading restart file. Path:" << endl << file << endl;
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
	svar.finPts = getInt(in, lineno, "Max number of particles");
	/*Get post processing options*/
	svar.cellSize = getDouble(in, lineno, "Post processing mesh size");
	svar.postRadius = getDouble(in, lineno, "Post processing support radius");
	/*Particle settings*/
	svar.dx = getDouble(in, lineno, "Particle actual spacing");	
	svar.Pstep = getDouble(in, lineno, "Particle resting spacing");
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
		avar.vJetMag = avar.vJet(1);
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
	
	string temp;
	getline(in,temp);
	lineno++;
	
	fvar.alpha = getDouble(in, lineno, "Artificial viscosity factor");
	fvar.contangb = getDouble(in, lineno, "Surface tension contact angle");
	fvar.rho0 = getDouble(in, lineno, "Fluid density rho0");
	avar.rhog = getDouble(in, lineno, "Air density rhog");
	fvar.Cs = getDouble(in, lineno, "Speed of sound");
	fvar.mu = getDouble(in, lineno, "Fluid viscosity");
	avar.mug = getDouble(in, lineno, "Air viscosity");
	fvar.sig = getDouble(in, lineno, "Surface Tension");
	svar.outfolder = getString(in,lineno, "Output folder name");
	svar.outdir = svar.outfolder;
	if(svar.Asource != 0)
	{
		svar.meshfile = getString(in, lineno, "Mesh input file");
  		svar.solfile = getString(in,lineno, "Mesh solution file");
  		svar.scale = getDouble(in,lineno, "Mesh scale");
  		avar.vRef = getDouble(in, lineno, "Gas ref Vel");
  		avar.pRef = getDouble(in, lineno, "Get ref Press");
  		avar.T = getDouble(in, lineno, "Gas ref Temp");
  	}

#ifdef DEBUG
	dbout << "Closing settings file"  << endl;
#endif
	in.close();	
}

void Read_Input_TECIO(string& infolder, SIM& svar, FLUID& fvar, AERO& avar)
{
	void* inputHandle = NULL;

	string inFile = svar.outfolder;
	inFile.append("Restart.szplt");

	if(tecFileReaderOpen(inFile.c_str(),&inputHandle))
	{
		cout << "Error opening szplt file. Path:" << endl;
		cout << inFile << endl;
		exit(-1);
	}

	int32_t frame = 1;

	int64_t iMax, jMax, kMax;
	if(tecZoneGetIJK(inputHandle, frame, &iMax, &jMax, &kMax))
	{
		cout << "Failed to read restart IJK info. Stopping." << endl;
		exit(-1);
	}

	if(iMax != SIMDIM)
	{
		cerr << "Dimension mismatch in restart file to compiled dimension. Stopping" << endl;
		exit(-1);
	}

	int32_t varnum = 1;
	int var;
	Read_Real_Value(inputHandle, frame, varnum, iMax, svar.framet, "frame time interval");
	Read_Int_Value(inputHandle, frame, varnum, iMax, var, "number of frames");
	svar.Nframe = var;
	Read_UInt_Value(inputHandle, frame, varnum, iMax, svar.outframe, "output frame info");
	Read_UInt_Value(inputHandle, frame, varnum, iMax, svar.outtype, "output data type");
	Read_UInt_Value(inputHandle, frame, varnum, iMax, svar.outform, "output content");
	Read_UInt_Value(inputHandle, frame, varnum, iMax, svar.gout, "ghost output");
	Read_UInt_Value(inputHandle, frame, varnum, iMax, svar.subits, "sub iterations");
	Read_Int_Value(inputHandle, frame, varnum, iMax, var, "max particles");
	svar.finPts = var;
	Read_Real_Value(inputHandle, frame, varnum, iMax, svar.cellSize, "cell size");
	Read_Real_Value(inputHandle, frame, varnum, iMax, svar.postRadius, "post process support radius");
	Read_Real_Value(inputHandle, frame, varnum, iMax, svar.dx, "dx");
	Read_Real_Value(inputHandle, frame, varnum, iMax, svar.Pstep, "Pstep");
	Read_Real_Value(inputHandle, frame, varnum, iMax, svar.Bstep, "Bstep");
	Read_Int_Value(inputHandle, frame, varnum, iMax, svar.Bcase, "boundary case");
	Read_Int_Value(inputHandle, frame, varnum, iMax, svar.Asource, "aerodynamic source");
	Read_Int_Value(inputHandle, frame, varnum, iMax, avar.acase, "aerodynamic case");
	Read_Int_Value(inputHandle, frame, varnum, iMax, svar.ghost, "ghost particle option");

	vector<real> vecreal;
	vector<int> vecint;

	Read_Real_Vector(inputHandle, frame, varnum, iMax, vecreal, "start coordinates");
	
	for(size_t dim = 0; dim < SIMDIM; dim++)
	{
		svar.Start(dim) = vecreal[dim];
	}

	if(svar.Bcase < 2)
	{
		Read_Int_Vector(inputHandle, frame, varnum, iMax, vecint, "particle dimension counts");

		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			svar.xyPART(dim) = vecint[dim];
		}

		Read_Real_Vector(inputHandle, frame, varnum, iMax, vecreal, "box dimensions");
		
		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			svar.Box(dim) = vecreal[dim];
		}

		Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.pPress, "starting pressure");
	}
	else
	{
		Read_Int_Value(inputHandle, frame, varnum, iMax, var, "diameter particle number");
		svar.nrad = var;

		Read_Real_Vector(inputHandle, frame, varnum, iMax, vecreal, "starting angles");
		
		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			svar.Angle(dim) = vecreal[dim];
		}

		Read_Real_Vector(inputHandle, frame, varnum, iMax, vecreal, "jet dimensions");
	
		for(size_t dim = 0; dim < 2; dim++)
		{
			svar.Jet(dim) = vecreal[dim];
		}

		
		Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.pPress, "starting pressure");

		Read_Real_Value(inputHandle, frame, varnum, iMax, avar.vJetMag, "jet velocity");
		

		
		Read_Real_Vector(inputHandle, frame, varnum, iMax, vecreal, "freestream velocity");
		
		for(size_t dim = 0; dim < SIMDIM; dim++)
		{
			avar.vInf(dim) = vecreal[dim];
		}

		// if(avar.acase == 2 || avar.acase == 3)
  		// {
		// 	Read_Real_Value(inputHandle, frame, varnum, iMax, avar.a, "parameter a");
		// 	Read_Real_Value(inputHandle, frame, varnum, iMax, avar.h1, "parameter h1");
		// 	Read_Real_Value(inputHandle, frame, varnum, iMax, avar.b, "parameter b");
		// 	Read_Real_Value(inputHandle, frame, varnum, iMax, avar.h2, "parameter h2");
  		// }
	}

	Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.alpha, "artificial dissipation factor");
	Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.contangb, "contact angle");
	Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.rho0, "liquid density");
	Read_Real_Value(inputHandle, frame, varnum, iMax, avar.rhog, "gas density");
	Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.Cs, "speed of sound");
	Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.mu, "liquid viscosity");
	Read_Real_Value(inputHandle, frame, varnum, iMax, avar.mug, "gas viscosity");
	Read_Real_Value(inputHandle, frame, varnum, iMax, fvar.sig, "surface tension");
	
	if(svar.Asource > 0)
	{
		Read_Real_Value(inputHandle, frame, varnum, iMax, svar.scale, "scale");
		Read_Real_Value(inputHandle, frame, varnum, iMax, avar.vRef, "reference velocity");
		Read_Real_Value(inputHandle, frame, varnum, iMax, avar.pRef, "reference pressure");
		Read_Real_Value(inputHandle, frame, varnum, iMax, avar.T, "reference temperature");
	}

	if(tecFileReaderClose(&inputHandle))
	{
		cerr << "Failed to close file" << endl;
	}
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

	if(svar.Asource != 0)
	{
		if(svar.outform != 3 || svar.outform != 5)
		{
			cout << "No cell information. Cannot restart" << endl;
			exit(-1);
		}
	}

#ifdef DEBUG
	dbout << "Reading frame info file for restart information." << endl;
#endif

	// Re-read the settings and fluid file
	// cout << svar.outfolder << endl;

	string outdir = svar.output_prefix;

	// Now get the data from the files. Start with the boundary
	if(svar.out_encoding == 0)
	{
		State boundary, fuel, rest;	

		// Read the fuel
		void* fuelHandle = NULL;
		void* boundHandle = NULL;

		int32_t fuelFrames, boundFrames;
		double fuelTime, boundTime;

		string fuelf = outdir;
		fuelf.append("_fuel.szplt");

		if(tecFileReaderOpen(fuelf.c_str(),&fuelHandle))
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << fuelf << endl;
			exit(-1);
		}

		// Check how many frames are in the fuel file.
		cout << "Checking Fuel file..." << endl;
		CheckContents(fuelHandle,svar,fuelFrames,fuelTime);

		if (svar.Bcase!= 3 && svar.Bcase != 4 && svar.Bcase !=0)
		{
			string boundf = outdir;
			boundf.append("Boundary.szplt");

			if(tecFileReaderOpen(boundf.c_str(),&boundHandle))
			{
				cout << "Error opening szplt file. Path:" << endl;
				cout << boundf << endl;
				exit(-1);
			}

			cout << "Checking Boundary file..." << endl;
			CheckContents(boundHandle,svar,boundFrames,boundTime);

			if(fuelFrames!= boundFrames)
			{
				cout << "Caution! Number of frames is not consistent between fuel and boundary files." << endl;
			}

			if(fuelTime != boundTime)
			{
				cout << "Caution! Frame times are not consistent between fuel and boundary files." << endl;

				if(fuelTime > boundTime)
				{
					double time = 0.0;
					for(int32_t frame = fuelFrames-1; frame > 1; frame--)
					{
						if(tecZoneGetSolutionTime(fuelHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from fuel file." << endl;
							continue;
						}
						
						if(time == boundTime)
						{
							cout << "Found the correct frame" << endl;
							fuelFrames = frame;
							fuelTime = time;
							break;
						}
					}
				}
				else
				{
					double time = 0.0;
					for(int32_t frame = boundFrames-1; frame >= 1; frame--)
					{
						if(tecZoneGetSolutionTime(boundHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from boundary file." << endl;
							continue;
						}
						

						if(time == fuelTime)
						{
							cout << "Found the correct frame" << endl;
							boundFrames = frame;
							boundTime = time;
							break;
						}
					}
				}

				if(fuelTime != boundTime)
				{
					cout << "Could not find a consistent time in each file. Stopping." << endl;
					exit(-1);
				}
			}

			cout << "Attempting to read the boundary..." << endl;
			Read_Binary_Timestep(boundHandle,svar,boundFrames,boundary);
		}

	    svar.frame = fuelFrames-1;

	    // Read the actual data.
		cout  << "Attempting to read the fuel..." << endl;
		Read_Binary_Timestep(fuelHandle,svar,fuelFrames,fuel);


		pn = boundary;
		svar.bndPts = boundary.size();
		svar.simPts = fuel.size();
		pn.insert(pn.end(),fuel.begin(),fuel.end());	
		svar.totPts = pn.size();
		
		if(svar.simPts + svar.bndPts != svar.totPts)
		{
			cout << "Mismatch of array sizes. Total array is not the sum of sim and boundary arrays" << endl;
			exit(-1);
		}

		// Check the frame timings to ensure dt does not go negative.
		double frametime = svar.framet * real(fuelFrames-1);

		if(svar.t > frametime)
		{
			// There's a problem with the frame times.
			cout << "Time is further ahead than what is supposed by the frame timings." << endl;
			cout << "Need to adjust frame time to match the current time" << endl;

			real framet = svar.t/real(fuelFrames-1.0);
			cout << "Old frame time: " << svar.framet << "  New frame time: " << framet << endl << endl;
			svar.framet = framet;			
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

		cout << "ASCII file restart not yet implemented." << endl;
		exit(-1);
	}

	// Go through the particles giving them the properties of the cell
	vector<size_t> buffer;
	#pragma omp parallel for 
	for(size_t ii = 0; ii < svar.totPts; ++ii)
	{
		pn[ii].partID = ii;
		pn[ii].p =  fvar.B*(pow(pn[ii].rho/fvar.rho0,fvar.gam)-1);

		if(pn[ii].b == PartState.BACK_)
		{
			#pragma omp critical
			svar.back.emplace_back(ii);
		}
		else if(pn[ii].b == PartState.BUFFER_)
		{
			#pragma omp critical
			buffer.emplace_back(ii);
		}
		
		// Initialise the rest of the values to 0
		pn[ii].s = 0.0;
		pn[ii].woccl = 0.0;
		pn[ii].pDist = 0.0;
		pn[ii].internal = 0;

		pn[ii].vPert = StateVecD::Zero();
	}	

	/* Put the particles into the buffer vector */
	svar.buffer = vector<vector<size_t>>(svar.back.size(),vector<size_t>(4));
	real eps = 0.01*svar.dx; /* Tolerance value */
	for(size_t ii = 0; ii < svar.back.size(); ++ii)
	{
		StateVecD const test = svar.Transp*(pn[svar.back[ii]].xi-svar.sim_start);
		
		for(size_t index = 0; index < buffer.size(); ++index)
		{
			/* Check which particle in the back vector it corresponds to by checking which  */
			/* particle it lies behind */
			StateVecD const xi = svar.Transp*(pn[buffer[index]].xi-svar.sim_start);
			// cout << test[0] << "  " << test[1] << "  " << xi[0] << "  " << xi[1] << endl;

			if(xi[0] < test[0] + eps && xi[0] > test[0] - eps)
			{	/* X coordinate is within bounds, so should lie behind this point */

				/* Start wit the furthest away, and go closer */
				if (xi[1] < test[1] - 4.0*svar.dx + eps)
				{
					svar.buffer[ii][3] = buffer[index];
				}
				else if (xi[1] < test[1] - 3.0*svar.dx + eps)
				{
					svar.buffer[ii][2] = buffer[index];
				}
				else if (xi[1] < test[1] - 2.0*svar.dx + eps)
				{
					svar.buffer[ii][1] = buffer[index];
				}
				else if (xi[1] < test[1] - svar.dx + eps)
				{
					svar.buffer[ii][0] = buffer[index];
				}
				else
				{
					cout << "Couldn't identify where to place the buffer particle" << endl;
					exit(-1);
				}
			}
		}
	}	
	
	pnp1 = pn;


/*	int width = 15;

	for(size_t ii = 0; ii < svar.totPts; ++ii)
	{
		cout << setw(5) << pnp1[ii].partID << setw(width) << pnp1[ii].xi(0) << setw(width) << pnp1[ii].xi(1) << setw(width) <<
		pnp1[ii].rho << setw(width) << pnp1[ii].Rrho << setw(width) << pnp1[ii].m << setw(width) << pnp1[ii].v(0) << 
		setw(width) << pnp1[ii].v(1) << setw(width) << pnp1[ii].f(0) << setw(width) << pnp1[ii].f(1) << setw(4) << 
		pnp1[ii].b << setw(width) << pnp1[ii].cellV(0) << setw(width) << pnp1[ii].cellV(1) << setw(width) << pnp1[ii].cellP
		<< setw(width) << pnp1[ii].cellID << endl;
	}*/
}

#endif