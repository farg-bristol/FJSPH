#ifndef TAUIO_H
#define TAUIO_H

#include "Var.h"
#include "Crossing.h"
#include "CDFIO.h"

using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string;


void Check_if_convex(MESH& cells)
{
	/*Perform a check using the cell centre to test if a cell is defined as convex*/
	#pragma omp parallel for
	for(uint ii = 0; ii < cells.elems.size(); ++ii)
	{
		StateVecD testp = cells.cCentre[ii];

		#if SIMDIM == 3
		if(!Crossings3D(cells.verts,cells.cFaces[ii],testp))
		{
			cout << "Cell is not defined as convex. Please check the cell definition" << endl;
			cout << "Cell ID: " << ii << " No Faces: " << cells.cFaces[ii].size() << 
			 " No vertices: " << cells.elems[ii].size() << endl;
			exit(-1); 
		}

		#endif
		#if SIMDIM == 2
		if(!Crossings2D(cells.verts, cells.elems[ii],testp))
		{
			cout << "Cell is not defined as convex. Please check the cell definition" << endl;
			cout << "Cell ID: " << ii << " No vertices: " << cells.cVerts[ii].size() << endl;
			exit(-1); 
		}
		#endif

	}
}

void Write_Mesh_Data(SIM &svar, MESH &cells)
{
	string mesh = svar.outfolder;
	mesh.append("/Mesh.plt");
	std::ofstream fm(mesh, std::ios::out);

	if(!fm.is_open())
	{ 
		cout << "Failed to open mesh output file" << endl;
		exit(-1);
	}

	#if SIMDIM == 2
		fm << "TITLE = \"2D TAU Solution\"\n";
		fm << "VARIABLES = \"x (m)\" \"z (m)\" \"x_velocity\" \"z_velocity\"\n";
		fm << "ZONE T=\"2D Solution Plane\"\n";
		fm << "VARLOCATION=([1-2]=NODAL,[3-4]=CELLCENTERED)\n";
		fm << "N=" << cells.numPoint << ", E=" << cells.numElem << ", F=FEBLOCK ET=Quadrilateral\n\n";
	#endif

	#if SIMDIM == 3
		fm << "TITLE = \"3D TAU Solution\"\n";
		fm << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\" \"x_velocity\" \"y_velocity\" \"z_velocity\"\n";
		fm << "ZONE T=\"3D Solution Plane\"\n";
		fm << "VARLOCATION=([1-3]=NODAL,[4-6]=CELLCENTERED)\n";
		fm << "N=" << cells.numPoint << ", E=" << cells.numElem << ", F=FEBLOCK ET=Brick\n\n";
	#endif

	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	uint kk = 0;
		for(uint jj = 0; jj < cells.numPoint; ++jj)
		{
			fm << cells.verts[jj][ii] << " ";
			kk++;

			if(kk == 5)
			{
				fm << "\n";
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fm << "\n";
	}

	for(uint ii = 0; ii < SIMDIM; ++ii)
	{	uint kk = 0;
		for(uint jj = 0; jj < cells.numElem; ++jj)
		{
			fm << cells.cVel[jj][ii] << " ";
			kk++;

			if(kk == 5)
			{
				fm << "\n";
				kk = 0;
			}
		}

		if(kk % 5 != 0)
			fm << "\n";
	}

	for(uint ii = 0; ii <= cells.numElem; ++ii)
	{	
		for(auto elem:cells.elems[ii])
		{
			fm << elem+1 << " ";
		}
		fm << "\n";
	}
}

uint Identify_Vol(string line, ZONE& zn)
{
	uint endof3D = 0;
	if(line.find("hexa")!=string::npos)
	{
		zn.name = "hexa";
		zn.ETtype = "Brick";
		zn.ctype = 1;
		/*N verts = 8, N faces = 6, N edges = 12*/
		#if SIMDIM == 2
			zn.nCverts = 4;
			zn.nF = 0;
			zn.nFverts = 0;
		#else
			zn.nF = 6;
			zn.nCverts = 8;
			zn.nFverts = 4;
		#endif
	}
	else if(line.find("tetra")!=string::npos)
	{
		zn.name = "tetra";
		zn.ETtype = "Tetrahedron";
		zn.ctype = 2;
		/*N verts = 4, N faces = 4, N edges = 6*/
		#if SIMDIM == 2
			zn.nCverts = 3;
			zn.nF = 0;
			zn.nFverts = 0;
		#else
			zn.nF = 4;
			zn.nCverts = 4;
			zn.nFverts = 3;
		#endif
	}
	else if(line.find("prism")!=string::npos)
	{
		zn.name = "prism";
		
		zn.ETtype = "Brick";
		zn.ctype = 3;
		#if SIMDIM == 2
			cout << "Tau Mesh is 3D. Simulation is 2D. Stopping" << endl;
			exit(-1);
		#else
			zn.nF = 5;
			zn.nCverts = 6;
			zn.nFverts = 3; /*The end faces*/
			zn.nFverts2 = 4;	/*The centre faces*/
		#endif
	}
	else if(line.find("pyra")!=string::npos)
	{
		zn.name = "pyra";
		zn.ETtype = "Brick";
		zn.ctype = 4;
		#if SIMDIM == 2
			cout << "Tau Mesh is 3D. Simulation is 2D. Stopping" << endl;
			exit(-1);
		#else
			zn.nF = 5;
			zn.nCverts = 5;
			zn.nFverts = 3; /*The end faces*/
			zn.nFverts2 = 4;	/*The centre faces*/
		#endif
	}
	else if(line.find("symmetry")!=string::npos)
	{
		#if SIMDIM == 2
			zn.name = "Symmetry Plane";
			zn.ctype = 1;
			zn.nCverts = 4;
			zn.nF = 0;
			zn.nFverts = 0;
		#else
			endof3D = 1;
		#endif
	}
	else
	{
		endof3D = 1;
	}
	return endof3D;
}

uint Get_Zone_Info(ifstream& fin, ZONE& zn)
{
	string line,zone;
	getline(fin, line); /*Get Zone line data. Tells which type of volume*/
	zn.lineNo++;

	uint endof3D = Identify_Vol(line,zn);
	
	getline(fin, line); /*Get numbers*/
	zn.lineNo++;
	zn.nP = 0;
	zn.nE = 0;
	std::size_t ptr2;
	std::size_t ptr = line.find("N=");
	if(ptr!=string::npos)
	{
		ptr2 = line.find_first_not_of("0123456789",ptr+2);
		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
		zn.nP = stoi(temp);
	}
	ptr = line.find("E=");
	if(ptr!=string::npos)
	{
		ptr2 = line.find_first_not_of("0123456789",ptr+2);
		string temp = line.substr(ptr+2,ptr2-(ptr+2));
		
		zn.nE = stoi(temp);
	}

	if(line.find("Quad")!=string::npos)
	{
		zn.ctype = 1;
	}
	else if (line.find("Triangle")!=string::npos)
	{
		zn.ctype = 2;
		zn.nCverts = 3;
		zn.nF = 0;
		zn.nFverts = 0;
	}


	return endof3D;
}

// std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
//     file.seekg(std::ios::beg);
//     for(uint ii=0; ii < num - 1; ++ii){
//         file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
//     }
//     return file;
// }

void Skip_Variable(ifstream& fin, const int np)
{
	for(int ii = 0; ii < ceil(float(np)/5.0); ++ii)
	{
		fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	}
}

vector<StateVecD> Get_Vector(ifstream& fin, ZONE& zn, const uint yskip)
{
	string line;
	vector<StateVecD> var(zn.nP,StateVecD::Zero());

	for(uint dim =0; dim < SIMDIM; dim++)
	{	
		uint kk = 0;
		/*If in 2D, skip the y-dimension, and read the z-dim.*/
		if(yskip == 1)
		{
			if (dim == 1)
			{
				for(uint ii = 0; ii < ceil(float(zn.nP)/5.0); ++ii)
				{
					fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

					// getline(fin,line);
				}
				zn.lineNo+= ceil(float(zn.nP)/5.0);
			}
		}
		for(int ii = 0; ii < ceil(float(zn.nP)/5.0); ++ii)
		{
			getline(fin, line);
			zn.lineNo++;
			std::istringstream sline(line);

			for (uint jj = 0; jj < 5; ++jj)
			{	
				double temp;
				sline >> temp;

				if (!sline)
					break;
				else
				{
					var[kk](dim) = temp;
					++kk;
				}
			}
		}
		if (kk!= var.size())
		{
			cout << "Mismatch of array size.\n" << 
			" Not all of the array has been populated." << endl;
			cout << "populated: " << kk << " Array size: " << var.size() << endl;
			cout << var[kk] << endl;
		} 
	}

	return var;
}

template <class T>
void Get_Scalar_Data(ifstream& fin, ZONE& zn, T& var)
{
	string line;

	uint kk = 0;
	for(uint ii = 0; ii < ceil(float(zn.nP)/5.0); ++ii)
	{
		getline(fin, line);
		zn.lineNo++;
		std::istringstream sline(line);

		for (uint jj = 0; jj < 5; ++jj)
		{	
			double temp;
			sline >> temp;

			if (!sline)
				break;
			else
			{
				var[kk] = temp;
				++kk;
			}
		}
	}

	if (kk!= var.size())
	{
		cout << "Mismatch of array size.\n" << 
		" Not all of the array has been populated." << endl;
		cout << "populated: " << kk << " Array size: " << var.size() << endl;
		cout << var[kk] <<  "  " << zn.lineNo << endl;
	}
}

void Get_3DCells(ifstream& fin, /*const vector<StateVecD> verts,*/ ZONE& zn,
 		std::vector<std::vector<std::vector<uint>>>& cFaces, 
 		std::vector<std::vector<uint>>& cell, vector<vector<StateVecD>>& cVerts)
{
	string line;
	if(zn.ctype == 1 || zn.ctype == 2)
	{
		for(uint ii = 0; ii < zn.nE; ++ii)
		{
			cVerts.emplace_back();
			getline(fin, line);
			zn.lineNo++;
			std::istringstream sline(line);

			for (uint jj = 0; jj < zn.nCverts; ++jj)
			{	
				uint temp;
				sline >> temp;
				cell[ii][jj] = (temp-1);
				// cVerts[ii].emplace_back(verts[temp-1]);								
			}
		}

		std::vector<std::vector<uint>> facenum;
		if(zn.ctype == 1)
		{	/*Hexahedron*/
			facenum = {{1,5,6,2},{2,6,7,3},{0,3,7,4},{0,4,5,1},{0,1,2,3},{7,6,5,4}};

			// cFaces = std::vector<std::vector<std::vector<StateVecD>>>(nE,
			// 		std::vector<std::vector<StateVecD>>(6,
			// 			std::vector<StateVecD>(4)));									
		}
		else if (zn.ctype == 2)
		{	/*Tetrahedron*/
			facenum = {{0,1,2},{0,2,3},{0,3,1},{1,3,2}};
			// facenum = {{2,1,0},{3,2,0},{1,3,0},{3,2,1}};
			// cFaces = std::vector<std::vector<std::vector<StateVecD>>>(nE, /*Number of elems*/
			// 		std::vector<std::vector<StateVecD>>(4,     /*Number of faces per elem*/
			// 			std::vector<StateVecD>(4)));	    /*Number of vertices per face*/			
		}

		Get_Cell_Faces(cell,facenum,cFaces);
	}
	else if (zn.ctype == 3)
	{	/*Prism*/
		for(uint ii = 0; ii < zn.nE; ++ii)
		{
			cVerts.emplace_back();
			getline(fin, line);
			zn.lineNo++;
			std::istringstream sline(line);
			uint count = 0;
			for (uint jj = 0; jj < 6; ++jj)
			{	
				uint temp;
				sline >> temp;
				cell[ii][jj] = (temp-1);
				// cVerts[ii].emplace_back(verts[temp-1]);

				if(count == 2 || count ==  5)
				{
					sline >> temp;
				}
				count++;
			}
		}

		std::vector<std::vector<uint>> facenum = {{1,4,5,2},{0,3,5,2},{0,3,4,1},{0,1,2},{5,4,3}};
		Get_Cell_Faces(cell,facenum,cFaces);
	}	
	else if (zn.ctype == 4)
	{	/*Pyramid*/
		for(uint ii = 0; ii < zn.nE; ++ii)
		{
			cVerts.emplace_back();
			getline(fin, line);
			zn.lineNo++;
			std::istringstream sline(line);
			for (uint jj = 0; jj < 5; ++jj)
			{	
				uint temp;
				sline >> temp;
				cell[ii][jj] = (temp-1);
				// cVerts[ii].emplace_back(verts[temp-1]);
			}
		}

		std::vector<std::vector<uint>> facenum = {{0,3,2,1},{0,1,4},{0,4,3},{1,4,2},{2,3,4}};
		Get_Cell_Faces(cell,facenum,cFaces);
	}	
}

void Get_2DCells(ifstream& fin, const vector<StateVecD> verts, ZONE& zn,
 	 std::vector<std::vector<uint>>& elems, std::vector<std::vector<StateVecD>>& cVerts)
{
	string line;
	if(zn.ctype == 1 || zn.ctype == 2)
	{	/*Hexa or Tetrahedron*/
		for(uint ii = 0; ii < zn.nE; ++ii)
		{
			getline(fin, line);
			zn.lineNo++;
			std::istringstream sline(line);

			for (uint jj = 0; jj < zn.nCverts; ++jj)
			{	
				uint temp;
				sline >> temp;
				elems[ii][jj] = (temp-1);								
			}
		}
	}
	else
	{
		cout << "Not sure what to do with these volumes in 2D. Stopping" << endl;
		exit(-1);
	}

	for(uint ii = 0; ii < zn.nE; ++ii)
	{
		for (uint jj = 0; jj < zn.nCverts; ++jj)
		{	
			if (ii > cVerts.size())
			{
				cout << "Loop attempted to access out of bounds." << endl;
				exit(-1);
			}
			if(elems[ii][jj]>verts.size())
			{
				cout << "Value in element list exceeds vertex list size." << endl;
				cout << ii << "  " << jj << "  " << elems[ii][jj] << endl;
				exit(-1);
			}
			cVerts[ii][jj] = verts[elems[ii][jj]];								
		}
	}

}

uint Get_Zone(ifstream& fin, string input, MESH& cells, ZONE& zone)
{
	string line;

	uint endof3D = Get_Zone_Info(fin, zone);

	if(endof3D == 1)
		return endof3D;

	cout << "Reading zone: " << zone.name << endl;
	cout << "\tNumber of vertices: " << zone.nP << " Number of elements: " << zone.nE << endl;
	getline(fin,line);
	zone.lineNo++;
	// cout << "Line: " << zone.lineNo << endl;

	/*************** START OF VERTICES DATA *******************/
	uint varcount = 0;
	/*Get the position vectors*/
	#if SIMDIM == 2	/*Skip y component*/
		vector<StateVecD> verts = Get_Vector(fin, zone, 1);
	#else /*get y for 3D sim*/
		vector<StateVecD> verts = Get_Vector(fin, zone, 0);
	#endif
	
	cells.verts.insert(cells.verts.end(),verts.begin(),verts.end());

	varcount +=3;
	/*Skip variables aside from the velocity vectors*/
	for (uint ii = 3; ii < zone.velstart ; ++ii)
	{
		if(varcount == zone.cpstart)
		{	/*If Cp data is encountered, then read it in.*/
			vector<ldouble> scalar(zone.nP,0.0);
			Get_Scalar_Data(fin, zone, scalar);
			if(zone.pressOrcp == 1)
				cells.pointCp.insert(cells.pointCp.end(),scalar.begin(),scalar.end());
			else
				cells.pointP.insert(cells.pointP.end(),scalar.begin(),scalar.end());
		}
		else if(varcount == zone.densstart)
		{	/*If density data is encountered, then read it in.*/
			vector<ldouble> scalar(zone.nP,0.0);
			Get_Scalar_Data(fin, zone, scalar);
			cells.pointRho.insert(cells.pointRho.end(),scalar.begin(),scalar.end());
		}
		else
		{
			Skip_Variable(fin,zone.nP);
		}
		varcount++;
	}

	/*Read velocity data*/
	#if SIMDIM == 2
		if (zone.veltype == 1) /*Don't skip since there isnt a y vel*/
		{
			vector<StateVecD> vels = Get_Vector(fin, zone, 0);
			cells.pVel.insert(cells.pVel.end(),vels.begin(),vels.end());
			varcount += 2;
		} 
		else /*Skip y velocity component*/
		{
			vector<StateVecD> vels = Get_Vector(fin, zone, 1);
			cells.pVel.insert(cells.pVel.end(),vels.begin(),vels.end());
			varcount +=3;
		}
	#else /*Get the 3D velocity*/
		if(zone.veltype == 0)
		{	
			vector<StateVecD> vels = Get_Vector(fin, zone, 0);
			cells.pVel.insert(cells.pVel.end(),vels.begin(),vels.end());
			varcount +=3;
		}
		else 
		{
			cout << "y-velocity component not present. Stopping." << endl;
			exit(-1);
		}
	#endif

	uint velend;
	#if SIMDIM == 2
		if (zone.veltype == 1)
			velend = SIMDIM;
		else
			velend = 3;
	#else
		velend = SIMDIM;
	#endif 

	for (uint ii = 0; ii < zone.nvar - (zone.velstart+velend); ++ii)
	{
		if(varcount == zone.cpstart)
		{	/*If Cp data is encountered, then read it in.*/
			vector<ldouble> scalar(zone.nP,0.0);
			Get_Scalar_Data(fin, zone, scalar);
			if(zone.pressOrcp == 1)
				cells.pointCp.insert(cells.pointCp.end(),scalar.begin(),scalar.end());
			else
				cells.pointP.insert(cells.pointP.end(),scalar.begin(),scalar.end());
		}
		else if(varcount == zone.densstart)
		{	/*If density data is encountered, then read it in.*/
			vector<ldouble> scalar(zone.nP,0.0);
			Get_Scalar_Data(fin, zone, scalar);
			cells.pointRho.insert(cells.pointRho.end(),scalar.begin(),scalar.end());
		}
		else
		{
			Skip_Variable(fin,zone.nP);
		}
		varcount++;		
	}

	if(varcount != zone.nvar)
	{
		cout << "Some point data has been missed. \nCell data won't be read correctly. Stopping." << endl;
		exit(-1);
	}	

	cout << "\tVertex data complete. Reading cells..." << endl;
	// cout << zone.ctype << endl;
	// cout << "Line: " << zone.lineNo << endl;

	/************ BEGINNING OF CELL CONNECTIVITY ******************/
	#if SIMDIM == 3
		vector<vector<uint>> celldata(zone.nE,vector<uint>(zone.nCverts));
		vector<vector<vector<uint>>> cFaces;
		vector<vector<StateVecD>> cVerts;


		Get_3DCells(fin,/*verts,*/zone,cFaces,celldata,cVerts);

		cells.elems.insert(cells.elems.end(),celldata.begin(),celldata.end());
		cells.cFaces.insert(cells.cFaces.end(), cFaces.begin(), cFaces.end());
		cells.cVerts.insert(cells.cVerts.end(), cVerts.begin(), cVerts.end());
	#endif

	#if SIMDIM == 2
		vector<vector<uint>> celldata(zone.nE,vector<uint>(zone.nCverts));
		std::vector<std::vector<StateVecD>> cVerts(zone.nE,
						std::vector<StateVecD>(zone.nCverts));
		Get_2DCells(fin,verts,zone,celldata,cVerts);

		cells.elems.insert(cells.elems.end(),celldata.begin(),celldata.end());
		cells.cVerts.insert(cells.cVerts.end(), cVerts.begin(), cVerts.end());

	#endif

	cout << "\tCell data complete." << endl;
	// cout << "Line: " << zone.lineNo << endl;
	cout << "\tZone capture complete. Continuing..." << endl;
	// Write_Zone(input, zone, cells);
	return endof3D;
} 


// std::vector<ldouble> CpToPressure(const std::vector<ldouble>& Cp, const FLUID& fvar)
// {
// 	std::vector<ldouble> press(Cp.size());
// 	#pragma omp parallel for shared(Cp)
// 	for (uint ii = 0; ii < Cp.size(); ++ii)
// 	{
// 		press[ii] = Cp[ii]*fvar.gasDynamic /*+ fvar.gasPress*/;
// 	}
// 	return press;
// }

void NormalisePressure(MESH &cells, const FLUID& fvar)
{
	/*Check for density and pressure information*/
	/*Create the data based on the other.*/
	if(cells.pointP.size() == 0)
	{
		cells.pointP = CpToPressure(cells.pointCp,fvar);
	}
	else
	{
		/*Normalise pressure to be in terms of the Tait equation*/
	/*I have no idea if this is conservative...*/
		for(uint ii = 0; ii < cells.pointP.size(); ++ii)
		{
			cells.pointP[ii] -= fvar.gasPress;
			// cells.pointRho[ii] = fvar.rhog*pow((cells.pointP[ii]/fvar.B +1),1/fvar.gam);
		}
	}

	/*Check if the density has been initialised*/
	if(cells.pointRho.size()!=cells.pointP.size())
	{	/*Make them the same size if different*/
		cells.pointRho.resize(cells.pointP.size());
	}

	for(uint ii = 0; ii < cells.pointP.size(); ++ii)
	{
		cells.pointRho[ii] = fvar.rho0 * pow((cells.pointP[ii]/fvar.B + 1),1/fvar.gam);
	}	
}

void Read_TAUPLT(string input, MESH& cells, FLUID& fvar)
{
	std::ifstream fin(input, std::ios::in);

	if(!fin.is_open())
	{
		cout << "Couldn't open mesh file. Stopping." << endl;
		cout << "Path attempted: " << input << endl;
		exit(-1);
	}
	else 
	{
		cout << "Mesh file open, reading data..." << endl;
	}
	ZONE zone;
	std::string line;
	getline(fin,line);
	getline(fin,line);
	zone.lineNo += 2;

	zone.veltype = 0;
	if(line.find("\"x_velocity\"")!=string::npos)
	{
		if(line.find("\"y_velocity\"")!=string::npos)
		{
			if(line.find("\"z_velocity\"")!=string::npos)
			{
				cout << "All velocity components found!" << endl;
			}
			else
			{
				cout << "velocity components \"u\" and \"v\" found, but no \"w\"..." << endl;
				cout << "I can't work with this. Stopping..." << endl;
				exit(-1);
			}
		}
		else if (line.find("\"z_velocity\"")!=string::npos)
		{
			cout << "velocity components \"u\" and \"w\" found, but no \"v\"..." << endl;
			#if SIMDIM == 2
				cout << "I can work with this. Continuing..." << endl;
				zone.veltype = 1;
			#else 
				cerr << "SIMDIM is 3, and \"v\" velocity component missing. Stopping." << endl;
				exit(-1);
			#endif
		}
		else
		{
			cout << "Only velocity component \"u\" found, but no \"v\" and \"w\"..." << endl;
			cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
		}
	}
	else if (line.find("\"y_velocity\"")!=string::npos)
	{
		if(line.find("\"z_velocity\"")!=string::npos)
		{
			cout << "velocity components \"v\" and \"w\" found, but no \"u\"..." << endl;
			cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
		}
		else
		{
			cout << "Only velocity component \"v\" found, but no \"u\" and \"w\"..." << endl;
			cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
		}
	}
	else if (line.find("\"z_velocity\"")!=string::npos)
	{
		cout << "Only velocity component \"w\" found, but no \"u\" and \"v\"..." << endl;
		cout << "I can't work with this. Stopping..." << endl;
			exit(-1);
	}
	else
	{
		cout << "Warning: No velocity components provided.\n" ;
		cout << "I can't work with this. Stopping..." << endl;
		exit(-1);
	}

	zone.nvar = std::count(line.begin(),line.end(),'\"');
	zone.nvar /= 2;
	std::size_t ptr = line.find("\"x_velocity\"");
	zone.velstart = std::count(line.begin(),line.begin()+ptr, '\"');
	zone.velstart/=2;

	/*Check to see if there is pressure data available*/
	ptr = line.find("\"pressure\"");
	zone.pressOrcp = 1;
	zone.cpstart = 0;
	if (ptr != string::npos)
	{
		cout << "Pressure data directly available!" << endl;
		zone.pressOrcp = 0;
		zone.cpstart = std::count(line.begin(),line.begin()+ptr, '\"');
		zone.cpstart/=2;
	}
	else
	{
		ptr = line.find("\"cp\"");	
		if (ptr != string::npos)
		{
			zone.cpstart = std::count(line.begin(),line.begin()+ptr, '\"');
			zone.cpstart/=2; 
		}
		else 
		{
			cout << "Couldn't find any pressure data" << endl;
		}
	}

	/*Check to see if there is density data available*/
	ptr = line.find("\"density\"");
	zone.densstart = 0;
	if (ptr != string::npos)
	{
		cout << "Density data directly available!" << endl;
		
		zone.densstart = std::count(line.begin(),line.begin()+ptr, '\"');
		zone.densstart/=2;
	}

	/*Next bit depends on dimensions. If 3D, read the hexa data.*/
	/*If 2D, skip this and read the symmetry plane data.*/
	#if SIMDIM == 3
		uint endof3D = 0;
		while (endof3D == 0)
			endof3D = Get_Zone(fin, input, cells, zone);

		if(cells.elems.size()==0)
		{
			cout << "No mesh data has been taken. \n" << 
			"Couldn't determine which type of volumes used." << endl;
			exit(-1);
		}
	#endif
	
	
	#if SIMDIM == 2
		getline(fin, line); /*Get Zone line data. Tells which type of volume*/
		zone.lineNo++;
		string zonename;
		zonename = Identify_Vol(line,zone);
		uint nP = 0;
		uint nE = 0;
		getline(fin, line); /*Get numbers*/
		zone.lineNo++;
		
		std::size_t ptr2;
		ptr = line.find("N=");
		if(ptr!=string::npos)
		{
			ptr2 = line.find_first_not_of("0123456789",ptr+2);
			string temp = line.substr(ptr+2,ptr2-(ptr+2));
			
			nP = stoi(temp);
		}
		ptr = line.find("E=");
		if(ptr!=string::npos)
		{
			ptr2 = line.find_first_not_of("0123456789",ptr+2);
			string temp = line.substr(ptr+2,ptr2-(ptr+2));
			
			nE = stoi(temp);
		}

		
		/*Skip the 3D data and read a symmetry plane*/
		uint lineNo = ceil(float(nP)/5.0)*zone.nvar + nE+6;
		GotoLine(fin,lineNo);
		zone.lineNo = lineNo;

		Get_Zone(fin, input, cells, zone);
	#endif
	
	fin.close(); 
	
	/*Build Cell connectivity*/
	/*Find the cell neighbours so that they can be checked 
		when a particle isn't in the cell any more*/

	/* check if a cell has 2 vertices the same in it as the checked cell
		If yes, then it's a neighbour. */
	cout <<  "Total vertices: " << cells.verts.size() << 
			" Total elements: " << cells.elems.size() << endl; 
	cells.numPoint = zone.nP;
	cells.numElem = zone.nE;
	
	cout << "Averaging point data to the cell..." << endl;
	NormalisePressure(cells,fvar);

	
	/*Average data from the points to find the cell based data*/
	cells.SetCells();
	StateVecD zero = StateVecD::Zero();
	Average_Point_to_Cell(cells.pVel,cells.cVel, cells.elems, zero);
	Average_Point_to_Cell(cells.pointRho,cells.cellRho,cells.elems,0.0);
	if(zone.pressOrcp == 1)
	{
		Average_Point_to_Cell(cells.pointCp,cells.cellCp,cells.elems,0.0);
		cells.cellP = CpToPressure(cells.cellCp,fvar);
	}
	else
		Average_Point_to_Cell(cells.pointP,cells.cellP,cells.elems,0.0);
	
	/*Find cell centres*/
	Average_Point_to_Cell(cells.verts,cells.cCentre,cells.elems,zero);

	/*Do a convex check*/
	// Check_if_convex(cells);

	cout << "End of mesh data intake. Starting SPH initialisation..." << endl;
	// Average_Point_to_Cell(cells.pointMach,cells.cellMach, cells.elems, 0.0);

}

#endif