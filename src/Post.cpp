
#include "PostVar.h"
#include "PostFunctions.h"

#include <tecio/TECIO.h>

void* Read_Binary_Info(SIM& svar, const string file)
{
	INTEGER4 I;
	void* inputHandle = NULL;
	I = tecFileReaderOpen(file.c_str(),&inputHandle);
	if(I == -1)
	{
		cout << "Error opening szplt file. Path:" << endl;
		cout << file << endl;
		exit(-1);
	}

	INTEGER4 numVars;
	I = tecDataSetGetNumVars(inputHandle, &numVars);
	cout << "Number of variables in output file: " << numVars << endl;
	if(numVars == 2 || numVars == 3)
    {
    	cout << "Only basic data has been output. No post processing can be done." << endl;
    	exit(-1);
    }

    /*Check if the data has components or not*/
    uint dataType = 0;
    if (DIM == 2)
    {
	    if (numVars == 7)
	    	dataType = 1;
	    else if (numVars == 10)
	    	dataType = 2;
		else if (numVars == 15)
			dataType = 3;
		else if (numVars == 9)
			dataType = 4;
		else if (numVars == 11)
			dataType = 5;
	}
	else if (DIM == 3)
	{
	    if(numVars == 8)
	    	dataType = 1;
	    else if(numVars == 13)	
	    	dataType = 2;
	    else if(numVars == 19)
	    	dataType = 3;
	    else if (numVars == 10)
	    	dataType = 4;
	    else if(numVars == 14)
	    	dataType = 5;
	}

    if(dataType!= svar.outform)
    {
    	cout << "Data is different from that in the settings file." << endl;
    	cout << "Changing to the file data type: " << dataType << endl;
    	svar.outform = dataType;
    }

    INTEGER4 numZones;
    I = tecDataSetGetNumZones(inputHandle, &numZones);
    cout << "Number of zones (timesteps) in output file: " << numZones << endl;

    if(file.find("Boundary.szplt") != std::string::npos) 
    {
    	if(svar.boutform == 1)
    	{
    		
	    	if(static_cast<uint>(numZones) != svar.Nframe+1)
		    {
	cout << "Warning: number of written timesteps is different from the frame count." << endl;
		    }
		    
		    
		    svar.Nframe = numZones;
    	}
    }
    else
    {
	    
    	if(static_cast<uint>(numZones) != svar.Nframe+1)
	    {
	cout << "Warning: number of written timesteps is different from the frame count." << endl;
	    }
	    
	    
	    svar.Nframe = numZones;
	}
	
    return inputHandle;
}

void Get_Sim_Info(int argc, char **argv, SIM& svar)
{

#ifdef DEBUG
	dbout << "Entering Get_Sim_Info." << endl;
#endif

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tERROR: No inputs provided. Stopping... \n";
    	exit(-1);    	
    }
    else if (argc > 2) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only one input arguments accepted,\n";
		cout << "1: Input file directoy\n";
		cout << "Other inputs will be ignored." << endl;
	}

	string file = argv[1];
	if(file.back() != '/')
    	file.append("/");

    svar.infolder = file;

	/*Read the settings file in the output folder*/
	file.append("Settings");

	std::ifstream sett(file);
	if(!sett.is_open())
	{
		cout << "Failed to open settings file." << endl;
		exit(-1);
	}

	cout << "Reading settings file." << endl;

	/*Skip lines*/
	string line;
	uint lineno=0;
	for(uint ii = 0; ii < 3; ii++)
	{
		getline(sett,line);
		lineno++;
	}

	svar.outtype = getInt(sett,lineno,"Output data type");
	svar.outform = getInt(sett,lineno,"Output data content");
	svar.boutform = getInt(sett,lineno,"Boundary Time Output");

	for(uint ii = 0; ii < 3; ii++)
	{
		getline(sett,line);
		lineno++;
	}

	// read in the settings
	svar.maxCells = getInt(sett,lineno,"Mesh cells along longest axis");
	svar.H = getDouble(sett,lineno,"Mesh support radius");
	svar.HSQ = svar.H*svar.H; 
	svar.sr = 4*svar.HSQ; 	/*KDtree search radius*/

	for(uint ii = 0; ii < 3; ii++)
	{
		getline(sett,line);
		lineno++;
	}

	svar.Bcase = getInt(sett,lineno,"Boundary case");


	getline(sett,line);
	lineno++;
	getline(sett,line);
	lineno++;

	/*Verify what dimension the system is*/
	getline(sett,line);
	lineno++;
	std::stringstream sline(line);
	
	vector<real> vec;
	real value;
	string temp;
	for(uint ii = 0; ii < 3; ii++)
	{
		sline >> temp;
		if(std::stringstream(temp)>> value)
		{
			vec.emplace_back(value);
		}
		temp = "";
	}

	if(vec.size() == 0)
	{
		cout << "Error: Dimension is zero. Check the input order." << endl;
		cout << "Line " << lineno << ": " << line << endl; 
		exit(-1);
	}	

	DIM = vec.size();

	sett.close();

	// Get information of the output folder. 
	file = svar.infolder;
	file.append("Fluid");

	lineno = 0;
	std::ifstream fluid(file);
	if(!fluid.is_open())
	{
		cout << "Failed to open fluid file." << endl;
		exit(-1);
	}

	uint nskip = 0;

	if(svar.Bcase == 6)
	{
		nskip = 16;
	}
	else
	{
		nskip = 10;
	}

	for(uint ii = 0; ii < nskip; ii++)
	{
		getline(fluid,line);
		lineno++;
	}

	svar.outfolder = getString(fluid,lineno,"Output Folder");
	fluid.close();

	file = svar.infolder;
	file.append(svar.outfolder);
	file.append("/");
	file.append("frame.info");

	std::ifstream f1(file);
	
	cout << "Simulation Parameters found..." << endl;
	cout << "Dimension: " << DIM << endl;
	cout << "Output data format: " ;
	if(svar.outtype == 0)
		cout << "Binary" << endl;
	else
		cout << "ACSII" << endl;

	cout << "Output data case: " << svar.outform << endl;
	cout << "Boundary case: " << svar.Bcase << endl;
 	cout << "Boundary time output: ";
	if(svar.boutform == 1)
		cout << "Yes" << endl;
	else
		cout << "No" << endl;

	cout << "Mesh cell max dimension size: " << svar.maxCells << endl;
	cout << "Mesh support radius: " << svar.H << endl;
	cout << "Output folder: " << svar.outfolder << endl;

#ifdef DEBUG
	dbout << "Exiting Get_Sim_Info." << endl;
#endif
			
}

vector<StateVecD> Read_Binary_Vector(void* inputHandle, const INTEGER4& frame, 
		INTEGER4& varCount, const INTEGER8& iMax)
{
	INTEGER4 I;
	vector<real> xVar(iMax);
	vector<real> yVar(iMax);
	vector<real> zVar;
    // cout << "Trying to get vector x-component. Var: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &xVar[0]);
	++varCount;
    // cout << "Trying to get vector y-component. Var: " << varCount << endl;
	I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &yVar[0]);
	++varCount;

	if(DIM == 3)
	{
		zVar = vector<real>(iMax);
		// cout << "Trying to get vector z-component. Var: " << varCount << endl;
		I = tecZoneVarGetDoubleValues(inputHandle, frame, varCount, 1, iMax, &zVar[0]);
		++varCount;
	}


	vector<StateVecD> vec(iMax,StateVecD(DIM));
	#pragma omp parallel for
	for(uint ii = 0; ii < xVar.size(); ++ii)
	{
		vec[ii](0) = xVar[ii];
		vec[ii](1) = yVar[ii];
		
		if(DIM == 3)
			vec[ii](2) = zVar[ii];
	}

	if(I == -1)
	{
		cout << "Errors occured obtaining vector. Variable count: " << varCount << endl;
	}

	return vec;
}


Point Read_Binary_Timestep(void* inputHandle, const INTEGER4 zoneNum, const uint dataType)
{
// variables = "x y z rho P m v_x v_y v_z a_x a_y a_z Cell_ID Cell_Vx Cell_Vy Cell_Vz Cell_P Cell_Rho";
    // cout << "Reading zone number: " << zoneNum << endl;

    INTEGER4 I;
    INTEGER8 iMax, jMax, kMax;
    I = tecZoneGetIJK(inputHandle, zoneNum, &iMax, &jMax, &kMax);
    if(I == -1)
    {
    	cout << "Reading zone data caused a tecplot error. Stopping." << endl;
    	exit(-1);
    }

    // cout << "Number of particles in zone: " << iMax << endl;

    Point tsData(iMax,dataType);
    I = tecZoneGetSolutionTime(inputHandle, zoneNum, &tsData.time);


    INTEGER4 varCount = 1;
    tsData.verts = Read_Binary_Vector(inputHandle, zoneNum, varCount, iMax); 
	

	/*Get density, pressure and mass*/
	I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.rho[0]);
	++varCount; //Skip Rrho
	++varCount;    
    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.m[0]);
	++varCount;  

	/*Check if the data has components or not*/
	if(dataType == 2 || dataType == 3 || dataType == 5)
	{	/*Data has components*/

	    /*Get velocity components*/
		tsData.vel = Read_Binary_Vector(inputHandle, zoneNum, varCount, iMax); 

		/*Get acceleration components*/
		tsData.acc = Read_Binary_Vector(inputHandle, zoneNum, varCount, iMax); 

		if(dataType == 3)
		{	/*Get the cell information for the points*/
			++varCount;
			tsData.cellV = Read_Binary_Vector(inputHandle, zoneNum, varCount, iMax); 

			/*Get density and pressure*/
			I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.cellRho[0]);
		    ++varCount;
		    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.cellP[0]);
			++varCount;
		}
	}
	else
	{
	    I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.vnorm[0]);
		++varCount;
		I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.anorm[0]);
		++varCount;

		if(dataType == 4)
		{
			++varCount;
			I = tecZoneVarGetDoubleValues(inputHandle, zoneNum, varCount, 1, iMax, &tsData.Af[0]);
			++varCount;
		}
	}	
    
	/*End of that timestep*/
	return tsData;
}


typedef class TECMESH {
	public:
		TECMESH(size_t DIM) : dim(DIM) {}

		void FindGridSize(SIM& svar)
		{
			cout << "Finding particle limits" << endl;
			// Get min/max from searching through data file
			// Start with boundary file. 
			INTEGER4 I = 0;
			INTEGER8 iMax, jMax, kMax;
			string file = svar.infolder;
			file.append(svar.outfolder);
			file.append("/");

			std::pair<StateVecD,StateVecD> minmax;
			/*Read in the timestep*/
			if(svar.Bcase != 0 && svar.Bcase !=5)
			{
				string boundFile = file;
				boundFile.append("Boundary.szplt");
				bFile = Read_Binary_Info(svar, boundFile);

				
			    I = tecZoneGetIJK(bFile, 1, &iMax, &jMax, &kMax);
			    if(I == -1)
			    {
			    	cout << "Reading zone data caused a tecplot error. Stopping." << endl;
			    	exit(-1);
			    }

				// Read the position vector
				INTEGER4 varCount = 1;
				vector<StateVecD> pos = Read_Binary_Vector(bFile, 1,varCount,iMax); 

				minmax = Find_MinMax(svar,pos);
				svar.minC = minmax.first;
				svar.maxC = minmax.second;
			}

			file.append("Fuel.szplt");
			fFile = Read_Binary_Info(svar, file);
		    I = tecZoneGetIJK(fFile, 1, &iMax, &jMax, &kMax);
		    if(I == -1)
		    {
		    	cout << "Reading zone data caused a tecplot error. Stopping." << endl;
		    	exit(-1);
		    }

		    /*Populate the Min/Max vectors so they aren't checked against ~0*/
		    // Read the position vector
			INTEGER4 varCount = 1;
			vector<StateVecD> pos = Read_Binary_Vector(fFile, 1,varCount,iMax); 
			minmax = Find_MinMax(svar,pos);

			if(svar.Bcase != 0 && svar.Bcase !=5)
			{
				// Compare against the boundary
				/*Check if they are more than the previous step*/
				for (size_t ii = 0; ii < dim; ++ii)
				{
					if(minmax.first(ii) < svar.minC(ii))
						svar.minC(ii) = minmax.first(ii);

					if(minmax.second(ii) > svar.maxC(ii))
						svar.maxC(ii) = minmax.second(ii);
				}
			}
			else
			{
				svar.minC = minmax.first;
				svar.maxC = minmax.second;
			}

			for(INTEGER4 frame = 2; frame <= static_cast<INTEGER4>(svar.Nframe); frame++)
			{
			    I = tecZoneGetIJK(fFile, frame, &iMax, &jMax, &kMax);
			    if(I == -1)
			    {
			    	cout << "Reading zone data caused a tecplot error. Stopping." << endl;
			    	exit(-1);
			    }

				// Read the position vector
				INTEGER4 varCount = 1;
				vector<StateVecD> pos = Read_Binary_Vector(fFile, frame,varCount,iMax); 

				minmax = Find_MinMax(svar,pos);

				/*Check if they are more than the previous step*/
				for (size_t ii = 0; ii < dim; ++ii)
				{
					if(minmax.first(ii) < svar.minC(ii))
						svar.minC(ii) = minmax.first(ii);

					if(minmax.second(ii) > svar.maxC(ii))
						svar.maxC(ii) = minmax.second(ii);
				}

			}

			/*data needs attaining from frame info*/
			minC = svar.minC;
			maxC = svar.maxC;
			step = StateVecD(dim);

			cout << "Maximum dimension: " << svar.maxC(0) << "  " << svar.maxC(1);
			if(dim == 3)
				cout << "  " << svar.maxC(2);
			cout << endl;
			cout << "Minimum dimension: " << svar.minC(0) << "  " << svar.minC(1);
			if(dim == 3)
				cout << "  " << svar.minC(2);
			cout << endl;

			// Find the cell step size
			// Find the max dimension
			real maxDiff;
			maxDiff = (fabs(maxC(0)-minC(0))+4*svar.H);
			if((fabs(maxC(1)-minC(1))+4*svar.H) > maxDiff)
				maxDiff = (fabs(maxC(1)-minC(1))+4*svar.H);

			if(dim == 3)
			{
				if((fabs(maxC(2)-minC(2))+4*svar.H) > maxDiff)
					maxDiff = (fabs(maxC(2)-minC(2))+4*svar.H);
			}

			real cellSize = maxDiff/(1.0*svar.maxCells);

			cout << "Maximum length: " << maxDiff << " Cell Size: " << cellSize << endl; 
			nX = ceil((fabs(maxC(0)-minC(0))+4*svar.H)/cellSize);
			step(0) = (fabs(maxC(0)-minC(0))+4*svar.H)/real(nX-1);
			nY = ceil((fabs(maxC(1)-minC(1))+4*svar.H)/cellSize);
			step(1) = (fabs(maxC(1)-minC(1))+4*svar.H)/real(nY-1);


			long unsigned int verts, cells;
			if (dim == 3)
			{
				nZ = ceil(fabs((maxC(2)-minC(2))+4*svar.H)/cellSize);
				step(2) = (fabs(maxC(2)-minC(2))+4*svar.H)/real(nZ-1);
				verts = nX*nY*nZ;
				cells = (nX-1)*(nY-1)*(nZ-1);
				nConns = 8*nCells;
				nFaces = 6;
			}
			else if (dim == 2)
			{
				verts = nX*nY;
				cells = (nX-1)*(nY-1);
				nConns = 4*nCells;
				nFaces = 4;
			}

			

			cout << "Dimensions numbers: " << nX << "  " << nY;
			if(dim == 3)
				cout << "  " << nZ;
			cout << endl;
			cout << "Dimension steps: " << step(0) << "  " << step(1);
			if(dim == 3)
				cout  << "  " << step(2);
			cout << endl;

			cout << "Total size: " << verts << "  " << cells << endl;
			nVerts = verts;
			nCells = cells;

			if(verts > std::numeric_limits<INTEGER4>::max())
			{
				cout << "Number of vertices exceeds the limit of the variable container." << endl;
				cout << "Please check your cell size and try again." << endl;
				exit(-1);
			}

			if(verts >  std::numeric_limits<size_t>::max())
			{
				cout << "Number of vertices exceeds the limit of the variable container." << endl;
				cout << "Please check your cell size and try again." << endl;
				exit(-1);
			}

			cout << "Grid Size Found. Vertices: " << nVerts << " Cells: " << nCells << endl;
			
			dataType = svar.outform;
			init(nVerts); 			
		}

		void Create_Grid(const SIM& svar)
		{

			cout << "Creating grid..." << endl;
			verts = vector<StateVecD>(nVerts,StateVecD(dim));
			if (dim == 3)
			{	
				xC = new real[verts.size()];
				yC = new real[verts.size()];
				zC = new real[verts.size()];

				for(INTEGER4 jj = 0; jj < nY; jj++)
					for(INTEGER4 ii = 0; ii < nX; ii++)
						for(INTEGER4 kk = 0; kk < nZ; kk++)
						{
							xC[index3(ii,jj,kk)] = minC(0)-2*svar.H + step(0)*ii;
							yC[index3(ii,jj,kk)] = minC(1)-2*svar.H + step(1)*jj;
							zC[index3(ii,jj,kk)] = minC(2)-2*svar.H + step(2)*kk;
						}

				// cout << maxC(0)+2*svar.H << "  " << maxC(1)+2*svar.H << endl;
				// cout << verts[nVerts-1](0) << "  " << verts[nVerts-1](1) << endl;
			}
			else if (dim == 2)
			{
				xC = new real[verts.size()];
				yC = new real[verts.size()];
				
				for(INTEGER4 jj = 0; jj < nY; jj++)
					for(INTEGER4 ii = 0; ii < nX; ii++)
					{
						xC[index2(ii,jj)] = minC(0)-2*svar.H + step(0)*ii;
						yC[index2(ii,jj)] = minC(1)-2*svar.H + step(1)*jj;
					}
			}
			
			if(verts.size() != static_cast<uint>(nVerts))
			{
				cout << "Mismatch of post process mesh dimensions." << endl;
				cout << "Mesh actual size: " << verts.size() << " Predicted: " << nVerts << endl;
				exit(-1);
			}
			// FindCellNeighbours()

			cout << "Grid made. Vertices: " << nVerts << " Cells: " << nCells << endl;
		}

		void Interp_Data(const Point& fluid, const SIM& svar)
		{
			INTEGER4 I;
			INTEGER4 DIsDouble                = 0;
		    if(sizeof(real) == 8)
				DIsDouble            = 1;

		    const INTEGER4 strandID = 1;     
		    const INTEGER4 parentZn                 = 0;
		    const INTEGER4 isBlock                  = 1;      /* Block */
		    const INTEGER4 nFConns                  = 0;
		    const INTEGER4 fNMode                   = 0;
		    const INTEGER4 shrConn                  = 0;
			const INTEGER4 zoneType = 0;

			const INTEGER4 iMax = nX;
			const INTEGER4 jMax = nY;
			INTEGER4 kMax;
			if (DIM == 3)
				kMax = nZ;
			else
				kMax = 1;
			

			// cout << fluid.size() << endl;
			/*Fluid from the timestep is read. build a neighbour tree.*/
			const Vec_Tree INDEX(DIM,fluid.verts,50);

			/*Now that there is the neighbour tree, interpolate the properties to the grid*/
			const nanoflann::SearchParams params;
			const real search_radius = svar.sr;

			#pragma omp parallel for
			for(uint ii=0; ii < verts.size(); ++ii)
			{
				StateVecD testp = verts[ii];
				std::vector<std::pair<size_t, real>> matches; /* Nearest Neighbour Search*/
				INDEX.index->radiusSearch(&testp[0], search_radius, matches, params);

				real ktemp = 0.0;
				real mass = 0.0;
				real dens = 0.0;

				if (dataType == 1 || dataType == 4)
				{
					real vels = 0.0;
					real accs = 0.0;
					real Af = 0.0;
					for(auto temp:matches)
					{
						uint jj = temp.first;
						/*Find the kernel*/
						real Rij = (fluid.verts[jj] - testp).norm();
						real k = W2Kernel(Rij,svar.H,1.0)/kernsum;

						ktemp += k;
						mass += k*fluid.m[jj];
						dens += k*fluid.rho[jj];
						vels += k*fluid.vnorm[jj];
						accs += k*fluid.anorm[jj];

						if(dataType == 4)
						{
							Af += k*fluid.Af[jj];
						}
					}
					if(dataType == 4)
						input(ii,ktemp,vels,accs,dens,mass,Af);
					else
						input(ii,ktemp,vels,accs,dens,mass);
				}

				else if (dataType == 2 || dataType == 5)
				{
					StateVecD vels = StateVecD::Zero(dim);
					StateVecD accs = StateVecD::Zero(dim);
					for(auto temp:matches)
					{
						uint jj = temp.first;
						/*Find the kernel*/
						real Rij = (fluid.verts[jj] - testp).norm();
						real k = W2Kernel(Rij,svar.H,1.0)/kernsum;

						ktemp += k;
						mass += k*fluid.m[jj];
						dens += k*fluid.rho[jj];
						vels += k*fluid.vel[jj];
						accs += k*fluid.acc[jj];
					}
					input(ii,ktemp,vels,accs,dens,mass);
				}
				else if (dataType == 3)
				{
					StateVecD vels = StateVecD::Zero(dim);
					StateVecD accs = StateVecD::Zero(dim);
					StateVecD cV = StateVecD::Zero(dim);
					real cR = 0.0;
					real cP = 0.0;
					for(auto temp:matches)
					{
						uint jj = temp.first;
						/*Find the kernel*/
						real Rij = (fluid.verts[jj] - testp).norm();
						real k = W2Kernel(Rij,svar.H,1.0)/kernsum;

						ktemp += k;
						vels += k*fluid.vel[jj];
						accs += k*fluid.acc[jj];
						dens += k*fluid.rho[jj];
						mass += k*fluid.m[jj];
						cV += k*fluid.cellV[jj];
						cR +=k*fluid.cellRho[jj];
						cP +=k*fluid.cellP[jj];
					}

					input(ii,ktemp,vels,accs,dens,mass,cV,cP,cR);
				}

			}

			
			/*Data has been interpolated to the grid, so now write to the file*/
			time = fluid.time;
			string zoneName = std::to_string(time);
			I = TECZNE142(zoneName.c_str(),
			                  &zoneType,
			                  &iMax,
			                  &jMax,
			                  &kMax,
			                  0,
			                  0,
			                  0,
			                  &time,
			                  &strandID,
			                  &parentZn,
			                  &isBlock,
			                  &nFConns,
			                  &fNMode,
			                  0,              /* TotalNumFaceNodes */
			                  0,              /* NumConnectedBoundaryFaces */
			                  0,              /* TotalNumBoundaryConnections */
			                  NULL,           /* PassiveVarList */
			                  NULL, 		  /* ValueLocation = Nodal */
			                  NULL,           /* SharVarFromZone */
			                  &shrConn);

			if(I == -1)
			{
				cout << "A tecplot error occured when making a new zone" << endl;
				exit(-1);
			}

			I = TECDAT142(&nVerts, xC, &DIsDouble);
			I = TECDAT142(&nVerts, yC, &DIsDouble);

			if (DIM == 3)
				I = TECDAT142(&nVerts, zC, &DIsDouble);
			

			// I = TECDAT142(&nVerts, kvalue, &DIsDouble);
			I = TECDAT142(&nVerts, rho, &DIsDouble);
			I = TECDAT142(&nVerts, m, &DIsDouble);

			
			if(dataType == 1)
			{
				I = TECDAT142(&nVerts, vnorm, &DIsDouble);
				I = TECDAT142(&nVerts, anorm, &DIsDouble);
			}
			else if(dataType == 2 || dataType == 3 || dataType == 5)
			{
				I = TECDAT142(&nVerts, vX, &DIsDouble);
				I = TECDAT142(&nVerts, vY, &DIsDouble);
				if (DIM == 3)
					I = TECDAT142(&nVerts, vZ, &DIsDouble);

				I = TECDAT142(&nVerts, aX, &DIsDouble);
				I = TECDAT142(&nVerts, aY, &DIsDouble);
				if (DIM == 3)
					I = TECDAT142(&nVerts, aZ, &DIsDouble);


				if (dataType == 3)
				{
					I = TECDAT142(&nVerts, cVX, &DIsDouble);
					I = TECDAT142(&nVerts, cVY, &DIsDouble);
					if (DIM == 3)
						I = TECDAT142(&nVerts, cVZ, &DIsDouble);

					I = TECDAT142(&nVerts, cellRho, &DIsDouble);
					I = TECDAT142(&nVerts, cellP, &DIsDouble);
				}
			}
			else if (dataType == 4)
			{
				I = TECDAT142(&nVerts, Aero, &DIsDouble);
			}

			INTEGER4 numZonesToRetain = 0;
		    I = TECFLUSH142(&numZonesToRetain,NULL);
		}

		void Write_Data(SIM& svar)
		{
			


			#if DEBUG
				const INTEGER4 Debug = 1;
			#else 
				const INTEGER4 Debug = 0;
			#endif
			INTEGER4 VIsDouble       = 0; 
			if(sizeof(real) == 8)
				VIsDouble            = 1;
		    
		    const INTEGER4 FileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
		    const INTEGER4 fileFormat = 1; // 0 == PLT, 1 == SZPLT
		    INTEGER4 I          = 0; /* Used to track return codes */
		    int* valueLocation = NULL;


		    string outfile = svar.infolder;
		    outfile.append(svar.outfolder);
		    outfile.append("/");
		    outfile.append("grid.szplt");
		    string variables;

		    if (DIM == 2)
		    {
		    	variables = "X Z";	
				if (dataType == 1)
				{
					variables = "X Z rho m v a";
					valueLocation = (int*)realloc(valueLocation, 6 * sizeof(int));
					nVar = 6;
				}
				else if (dataType == 2 || dataType == 5)
				{
					variables = "X Z rho m v_x v_z a_x a_z";
					valueLocation = (int*)realloc(valueLocation, 8 * sizeof(int));
					nVar = 8;
				}
				else if (dataType == 3)
				{
		variables = "X Z rho m v_x v_z a_x a_z Cell_Vx Cell_Vz Cell_Rho Cell_P Cell_ID";
					valueLocation = (int*)realloc(valueLocation, 13 * sizeof(int));
					nVar = 13;
				}
				else if (dataType == 4)
				{
					variables = "X Z rho m v a Aero";
					valueLocation = (int*)realloc(valueLocation, 7 * sizeof(int));
					nVar = 7;
				}
			}

			else if (DIM == 3)
			{
				variables = "X Y Z";  
				if (dataType == 1)
				{
					variables = "X Y Z rho m v a";
					valueLocation = (int*)realloc(valueLocation, 7 * sizeof(int));
					nVar = 7;
				}
				else if (dataType == 2 || dataType == 5)
				{
					variables = "X Y Z rho m v_x v_y v_z a_x a_y a_z";
					valueLocation = (int*)realloc(valueLocation, 11 * sizeof(int));
					nVar = 11;
				}
				else if (dataType == 3)
				{
					variables = 
	"X Y Z rho m v_x v_y v_z a_x a_y a_z Cell_Vx Cell_Vy Cell_Vz Cell_Rho Cell_P Cell_ID";
					valueLocation = (int*)realloc(valueLocation, 17 * sizeof(int));
					nVar = 17;
				}
				else if (dataType == 4)
				{
					variables = "X Y Z rho m v a Aero";
					valueLocation = (int*)realloc(valueLocation, 8 * sizeof(int));
					nVar = 8;
				}
			}

			for(uint ii = 0; ii < nVar; ++ii)
				valueLocation[ii] = 1;

			cout << "Opening the grid file." << endl;
			I = TECINI142((char*)"SPH Fluid Grid",
						  variables.c_str(),
						  outfile.c_str(),
						  (char*)".",
						  &fileFormat,
		                  &FileType,
		                  &Debug,
		                  &VIsDouble);
			if(I == -1)
			{
				cout << "Error opening output grid file for post processing." << endl;
				exit(-1);
			}

			GetKernelSum(svar,fFile);
			
			/*Read the data from the files into something.*/
			uint frameout = 50;
			if(svar.Bcase !=0 && svar.Bcase != 5 && svar.boutform == 1)
			{	/*The boundary is being written each timestep*/
				for( INTEGER4 zone = 1; zone <= static_cast<INTEGER4>(svar.Nframe); ++zone)
				{
					Point boundary = Read_Binary_Timestep(bFile, zone, dataType);
					Point fuel = Read_Binary_Timestep(fFile, zone, dataType);

					Point fluid = boundary;
					fluid.append(fuel);
					
					Interp_Data(fluid,svar);

					if(zone%frameout == 0)
					{
						cout << "Frames processed: " << zone << endl;
						cout << "Particle count on last timestep: " << fluid.size() << endl;
					}
				}
			}
			else
			{
				Point boundary;

				if(svar.Bcase != 0 && svar.Bcase !=5)
				{
					boundary = Read_Binary_Timestep(bFile, 1, dataType);
					/*Close the boundary file*/
					I = tecFileReaderClose(&bFile);
				}

				for( INTEGER4 zone = 1; zone <= static_cast<INTEGER4>(svar.Nframe); ++zone)
				{
					Point fuel = Read_Binary_Timestep(fFile, zone, dataType);

					Point fluid = boundary;
					fluid.append(fuel);
					
					Interp_Data(fluid, svar);
					if(zone%frameout == 0)
					{
						cout << "Frames processed: " << zone << endl;
						cout << "Particle count on last timestep: " << fluid.size() << endl;
					}
				}

			}

			TECEND142();
			free(valueLocation);
			
		}
		

		void DoPostProcessing(SIM& svar)
		{
			/*Check the folder for the grid szplt files*/

			string file = svar.infolder;
			file.append(svar.outfolder);
	  		file.append("/grid.szplt.szdat");
	  		struct stat info;
	  		if(stat( file.c_str(), &info ) == 0)
	  		{
		  		string cmd = "exec rm -r \"";
		  		cmd.append(svar.infolder);
		  		cmd.append(svar.outfolder);
		  		cmd.append("/");
		  		cmd.append("\"grid.szplt.sz*");
		  		if(system(cmd.c_str()))
		  		{
			    	cout << "System command failed to execute." << endl;
			    	exit(-1);
			    }
			}

			FindGridSize(svar);
			Create_Grid(svar);
			Write_Data(svar);
		}

		
	private:
		std::pair<StateVecD,StateVecD> Find_MinMax(SIM& svar, const vector<StateVecD>& xi)
		{
			/*Find the max and min positions*/
			auto xC = std::minmax_element(xi.begin(),xi.end(),
						[](StateVecD p1, StateVecD p2){return p1(0)< p2(0);});

			auto yC = std::minmax_element(xi.begin(),xi.end(),
						[](StateVecD p1, StateVecD p2){return p1(1)< p2(1);});

			StateVecD minC(dim);
			StateVecD maxC(dim);

			if(DIM == 3)
			{
				auto zC = std::minmax_element(xi.begin(),xi.end(),
						[](StateVecD p1, StateVecD p2){return p1(2)< p2(2);});

				// cout << "Min Pos: " << xC.first - xi.begin() << " " <<  yC.first - xi.begin()
				// << " " << zC.first - xi.begin() << endl;
				// cout << "Min Value: " << (*xC.first)(0) << " " << (*yC.first)(1) << " " << (*zC.first)(2) << endl; 
				minC(0) = (*xC.first)(0);
				minC(1) = (*yC.first)(1);
				minC(2) = (*zC.first)(2);


				// cout << "Max Pos: " << xC.second - xi.begin() << " " <<  yC.second - xi.begin()
				// << " " << zC.second - xi.begin() << endl;
				// cout << "Max Value: " << (*xC.second)(0) << " " << (*yC.second)(1) << " " << (*zC.second)(2) << endl; 

				maxC(0) = (*xC.second)(0);
				maxC(1) = (*yC.second)(1);
				maxC(2) = (*zC.second)(2);
			}
			else
			{
				minC(0) = xC.first->x();
				minC(1) = yC.first->y();

				maxC(0) = xC.second->x();
				maxC(1) = yC.second->y();
			}


				return std::pair<StateVecD,StateVecD>(minC,maxC);
		}

		void GetKernelSum(const SIM& svar, void* fluidFile)
		{

			Point fluid = Read_Binary_Timestep(fluidFile, 1, dataType);
			/*Fluid from the timestep is read. build a neighbour tree.*/
			const Vec_Tree INDEX(DIM,fluid.verts,50);

			/*Now that there is the neighbour tree, interpolate the properties to the grid*/
			const nanoflann::SearchParams params;
			const real search_radius = svar.sr;

			/*Find the maximum kernel sum to normalise by*/
			real ksum = 0.0;
			#pragma omp parallel for
			for(uint ii=0; ii < verts.size(); ++ii)
			{
				StateVecD testp = verts[ii];
				std::vector<std::pair<size_t, real>> matches; /* Nearest Neighbour Search*/
				INDEX.index->radiusSearch(&testp[0], search_radius, matches, params);

				real ktemp = 0.0;
				for(auto temp:matches)
				{
					uint jj = temp.first;
					/*Find the kernel*/
					real Rij = (fluid.verts[jj] - testp).norm();
					real k = W2Kernel(Rij,svar.H,1.0);

					ktemp += k;
				}

				if (ktemp > ksum)
					ksum = ktemp;
			}

			kernsum = ksum;
			cout << ksum << endl;
		}

		const uint index2(const uint ii, const uint jj)
		{
			return(ii + jj*nX);
		}
		const uint index3(const uint ii, const uint jj, const uint kk)
		{
			return ((kk*nY + jj)*nX + ii);
		}

		void init(const uint size)
		{
			rho = new real[size];
			m = new real[size];
			kvalue = new real[size];

			if (dataType == 1)
			{
				vnorm = new real[size]; anorm = new real[size]; 
			}
			else if(dataType == 2 || dataType == 3 || dataType == 5)
			{
				vX = new real[size]; vY = new real[size]; 
				aX = new real[size]; aY = new real[size];
				if (dim == 3)
				{
					vZ = new real[size]; aZ = new real[size]; 
				} 

				if(dataType == 3)
				{	
					cVX = new real[size]; cVY = new real[size];  
					cellP = new real[size]; cellRho = new real[size]; 
					if (DIM == 3)
						cVZ = new real[size];
				}
			}
		}
		// Input for data case 1
		void input(const uint ii, const real k, const real& vels, 
			const real& acc, const real dens, const real mass)
		{
			vnorm[ii] = vels; anorm[ii] = acc; kvalue[ii] = k;
			rho[ii] = dens; m[ii] = mass;
		}
		// Function for data case 2
		void input(const uint ii, const real k, const StateVecD& vels, 
		const StateVecD& acc, const real dens, const real mass)
		{
			vX[ii] = vels(0); vY[ii] = vels(1); aX[ii] = acc(0); aY[ii] = acc(1); 
			rho[ii] = dens; m[ii] = mass; kvalue[ii] = k;
			if (DIM == 3)
			{
				vZ[ii] = vels(2);  aZ[ii] = acc(2); 
			}
		}
		// Function for data case 3
		void input(const uint ii, const real k, const StateVecD& vels, const StateVecD& acc, 
		const real dens, const real mass, const StateVecD& cellV, 
		const real cP, const real cR)
		{
			kvalue[ii] = k;
			vX[ii] = vels(0); vY[ii] = vels(1); aX[ii] = acc(0); aY[ii] = acc(1); 
			rho[ii] = dens; m[ii] = mass; cellP[ii] = cP; cellRho[ii] = cR;
			cVX[ii] = cellV(0); cVY[ii] = cellV(1); 
			if (DIM == 3)
			{
				vZ[ii] = vels(2);  aZ[ii] = acc(2); cVZ[ii] = cellV(2);
			}	
		}

		// Function for data case 4
		void input(const uint ii, const real k, const real& vels, const real& acc,
				const real dens, const real mass, const real Af)
		{
			vnorm[ii] = vels; anorm[ii] = acc; kvalue[ii] = k;
			rho[ii] = dens; m[ii] = mass; Aero[ii] = Af;
		}

		std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
		{
		    file.seekg(file.beg);
		    for(uint ii=0; ii < num - 1; ++ii){
		        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
		    }
		    return file;
		}
		
		const size_t dim;

		void* bFile;
		void* fFile;

		StateVecD minC, maxC, step;
		INTEGER4 nVerts, nCells, nFaces, nConns, nX, nY, nZ;
		double time;
		INTEGER4 zone;
		uint dataType, nVar;
		vector<StateVecD> verts;

		real kernsum;
		real* xC;
		real* yC;
		real* zC;
		real* kvalue;
		real* vnorm;
		real* anorm;
		real* vX;
		real* vY;
		real* vZ;
		real* aX;
		real* aY;
		real* aZ;
		real* cVX;
		real* cVY;
		real* cVZ;
		real* rho;
		real* m;
		real* Aero;
		real* cellP;
		real* cellRho;
}TECMESH;



int main(int argc, char* argv[])
{
	write_header();

	SIM svar;

	Get_Sim_Info(argc,argv,svar);

  	// cout << "Found path: " << pathname << endl;
	TECMESH postgrid(DIM);

	postgrid.DoPostProcessing(svar);

	cout << "Post Processing complete!" << endl;
	exit(0);


	return 1;
}