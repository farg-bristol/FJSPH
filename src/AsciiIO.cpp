/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "AsciiIO.h"
#include "Var.h"
#include <regex>

#ifdef __STDC_ALLOC_LIB__
#define __STDC_WANT_LIB_EXT2__ 1
#else
#define _POSIX_C_SOURCE 200809L
#endif

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/
void ASCII_Restart(SIM& svar, FLUID const& fvar, SPHState& pn)
{
    #ifdef DEBUG
		fprintf(dbout,"Reading ASCII files for restart information.\n");
	#endif

	string fuelf = svar.restart_prefix;
	fuelf.append("_fuel.dat");

	FILE* fin = fopen(fuelf.c_str(),"r");

	if(fin == NULL)
	{
		cout << "Failed to open ASCII restart file: " << fuelf << endl;
		cout << "Stopping." << endl;
		exit(-1);
	}

	size_t const len = 200;
	char line[len];
	// string line;
	
	if(fgets(line,len,fin))
	{
		// Read the contents of the line
		
		cout << line << endl;
	}
	


	#ifdef DEBUG
		fprintf(dbout,"Exiting ASCII restart function\n");
	#endif
}


/*************************************************************************/
/**************************** ASCII OUTPUTS ******************************/
/*************************************************************************/
void Write_ASCII_header(FILE* fp, SIM const& svar, char const* title)
{
	string var_names = svar.var_names;
	std::regex r(",");
	string var = std::regex_replace(var_names,r,"\",\"");
	fprintf(fp, "TITLE=\"%s\"\n",title);
	fprintf(fp, "VARIABLES=\"%s\"",var.c_str());
}

inline void Write_ASCII_Vector(FILE* fp, vector<real> const& vec, size_t const& start, size_t const& end)
{
	size_t newl = 0;
	for(size_t ii = start; ii < end; ++ii)
	{
		fprintf(fp,"%3.7e",vec[ii]);
		if(newl > 4)
		{
			fprintf(fp,"\n");
			newl=0;
		}
		else
			newl++;
	}
}

inline void Write_ASCII_Vector(FILE* fp, vector<StateVecD> const& vec, size_t const& start, size_t const& end)
{
	for(size_t ii = start; ii < end; ++ii)
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",vec[ii][dim]);
			fprintf(fp,"\n");  
	  	}
	}
}

void Write_ASCII_Timestep(SIM& svar, SPHState const& pnp1, /* DELTAP const& dp, */
	size_t const& start, size_t const& end, char const* name, uint const& strandID, FILE* fp)
{
	// if(bwrite == 1)
	//  	fp <<  "ZONE T=\"" << name << "\"";
	// else
	fprintf(fp,"ZONE T=\"%s\", I=%zu, DATAPACKING=BLOCK, STRANDID=%d, SOLUTIONTIME=%.7g\n",
			name,end-start,strandID,svar.t);

    // fp <<", I=" << end - start << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
    // fp << std::left << std::scientific << std::setprecision(6);
    // const static uint width = 15;
    
    if(svar.outvar[0])	// Position coordinates, scaled
    {
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].xi(dim)/svar.scale);
			fprintf(fp,"\n");  
	  	}
	}	

	if(svar.outvar[1])	// Acceleration vector
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].acc(dim));
			fprintf(fp,"\n");  
	  	}
	}

	if(svar.outvar[2])	// Velocity vector
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].v(dim));
			fprintf(fp,"\n");  
	  	}
	}

	if(svar.outvar[3])	// Pressure scalar
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].p);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}

	if(svar.outvar[4])	// Density gradient
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].Rrho);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}

	if(svar.outvar[5])	// Particle ID, integer
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%zu",pnp1[ii].partID);
			if(newl > 6)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[6])	// cellID
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%zu",pnp1[ii].cellID);
			if(newl > 6)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}

	if (svar.outvar[7])	// Boundary type
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp, "%u", pnp1[ii].b);
			if (newl > 6)
			{
				fprintf(fp, "\n");
				newl = 0;
			}
			else
				newl++;
		}
		fprintf(fp, "\n");
	}
	
	if(svar.outvar[8])	// Density 
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].rho);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[9])	// Velocity magnitude
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].v.norm());
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}

	if(svar.outvar[10])	// Surface flag
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3d",pnp1[ii].surf);
			if(newl > 8)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[11])	// Surface zone flag
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3d",pnp1[ii].surfzone);
			if(newl > 8)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[12])	// Aerodynamic force magnitude
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].Af.norm());
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[13])	// Aerodynamic force vector
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].Af(dim));
			fprintf(fp,"\n");  
	  	}
	}
	
	if(svar.outvar[14])	// Curvature scalar
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].curve);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[15])	// Occlusion factor
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].woccl);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[16])	// Cell pressure
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].cellP);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[17])	// Cell density
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].cellRho);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[18])	// Cell velocity magnitude
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].cellV.norm());
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
	
	if(svar.outvar[19])	// Cell velocity vector
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].cellV(dim));
			fprintf(fp,"\n");  
	  	} 
	}
	
	if(svar.outvar[20]) // dSPH density gradient
	{
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].gradRho(dim));
			fprintf(fp,"\n");  
	  	} 
	}
		
	if(svar.outvar[21]) // Lambda eigenvalue of L matrix
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].lam);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
			
	if(svar.outvar[22]) // Lambda eigenvalue without boundary
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].lam_nb);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
			
	if(svar.outvar[23]) // Surface colour function
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].colour);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}

	if(svar.outvar[24]) // Surface colour gradient
	{
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].colourG);
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
		
	if(svar.outvar[25]) // Surface normal
	{	// Cell velocity vector
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
				fprintf(fp,"%3.7e",pnp1[ii].norm(dim));
			fprintf(fp,"\n");  
	  	} 
	}
		
	if(svar.outvar[26])
	{	// Shifting velocity magnitude
		size_t newl = 0;
		for (size_t ii = start; ii < end; ++ii)
		{
			fprintf(fp,"%3.7e",pnp1[ii].vPert.norm());
			if(newl > 4)
			{
				fprintf(fp,"\n");
				newl=0;
			}
			else
				newl++;  
	  	}
		fprintf(fp,"\n");  
	}
		
	if(svar.outvar[27])
	{	// Shifting velocity vector
		for (size_t ii = start; ii < end; ++ii)
		{
			for(uint dim = 0; dim < SIMDIM; ++dim)
			fprintf(fp,"%3.7e",pnp1[ii].vPert(dim));
			fprintf(fp,"\n");  
	  	} 
	}
}

void Write_Cell_Data(MESH const& cdata)
{
	#ifdef DEBUG
		fprintf(dbout,"Entering Write_Cell_Data...\n");
	#endif

	cout << "Writing cell based data." << endl;

	std::ofstream fout("Cell.dat", std::ios::out);
	if (!fout.is_open())
	{
		cout << "Failed to open data file for writing mesh." << endl;
		exit(-1);
	}

	fout << "TITLE = \"3D Mesh Solution\"\n";
	fout << "VARIABLES = \"x (m)\" \"y (m)\" \"z (m)\"\n";
	fout << "ZONE T=\"Cell Data\"" << endl;
	fout << "N=" << cdata.nPnts << ", E=" << cdata.nElem << ", F=FEBLOCK, ET=BRICK" << endl
		 << endl;

	/*Write vertices*/
	fout << std::left << std::scientific << std::setprecision(6);
	fout << std::setw(1);
	for (uint ii = 0; ii < SIMDIM; ++ii)
	{
		uint kk = 0;
		for (uint jj = 0; jj < cdata.verts.size(); ++jj)
		{
			fout << std::setw(15) << cdata.verts[jj][ii];
			kk++;

			if (kk == 5)
			{
				fout << endl;
				fout << std::setw(1);
				kk = 0;
			}
		}

		if (kk % 5 != 0)
			fout << "\n";
	}
	
	/*Write element indexing*/
	fout << std::fixed;
	for (uint index = 0; index < cdata.cFaces.size(); ++index)
	{
		/* Build the element unique index list */
		vector<size_t> elem;
		for (auto const& faceID: cdata.cFaces[index])
		{
			vector<size_t> const face = cdata.faces[faceID];
			for (auto const& vert : face)
			{
				if (std::find(elem.begin(), elem.end(), vert) == elem.end())
				{ /*Vertex doesn't exist in the elems vector yet.*/
					#pragma omp critical
					{
						elem.emplace_back(vert);
					}
				}
			}
		}

		for (auto const& vert : elem)
		{
			fout << std::setw(6) << vert + 1;
		}
		fout << "\n";
	}

	fout.close();

	#ifdef DEBUG
		fprintf(dbout,"Exiting Write_Cell_Data...\n");
	#endif
}

void Write_Face_Data(MESH const& cells)
{
	
	std::ofstream f1("foam_mesh.dat", std::ios::out);
	f1 << "VARIABLES= \"X\", \"Y\", \"Z\"" << endl;
	uint w = 17;
	f1 << std::left << std::scientific << std::setprecision(8);

	/*Write zone header information*/
	f1 << "ZONE T=\"OpenFOAM MESH\"" << endl;
	f1 << "ZONETYPE=FEPOLYHEDRON" << endl;
	f1 << "NODES=" << cells.nPnts << " ELEMENTS=" << cells.nElem << " FACES=" << cells.nFace << endl;
	size_t TotalNumFaceNodes = cells.nFace * 3;
	f1 << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
	f1 << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

	
	/*Write vertices in block format*/
	size_t n = 0;

	for (size_t dim = 0; dim < SIMDIM; ++dim)
	{
		for (auto const &pnt : cells.verts)
		{
			f1 << std::setw(w) << pnt[dim];
			n++;

			if(n == 6)
			{
				f1 << endl;
				n = 0;
			}
		}
		f1 << endl;
	}

	/*Write how many vertices per face*/
	n = 0;
	for (size_t ii = 0; ii < cells.nFace; ++ii)
	{
		f1 << std::setw(5) << 3;
		n++;

		if(n == 6)
		{
			f1 << endl;
			n = 0;
		}
	}
	f1 << endl;

	/*Write the face vertex list*/
	n = 0;
	w = 10;
	for (auto const &face : cells.faces)
	{
		for (auto const &vertex : face)
		{ 
			f1 << std::setw(w) << vertex + 1;
			n++;

			if (n == 6)
			{
				f1 << endl;
				n = 0;
			}
		}
	}
	f1 << endl;

	/*Write left elements*/
	n = 0;
	for (auto const &lr : cells.leftright)
	{
		f1 << std::setw(w) << lr.first+1;
		n++;

		if (n == 6)
		{
			f1 << endl;
			n = 0;
		}
	}
	f1 << endl;

	/*Write right elements*/
	n = 0;
	for (auto const &lr : cells.leftright)
	{
		if(lr.second < 0)
		{
			f1 << std::setw(w) << 0;
		}
		else
		{
			f1 << std::setw(w) << lr.second+1;
		}
		n++;

		if(n == 6)
		{
			f1 << endl;
			n = 0;
		}
	}
	f1 << endl	<< endl;
	f1.close();
	// exit(0);
}
