/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IPT_H
#define IPT_H

#include "Var.h"
#include "Containment.h"
#include "IO.h"
#include "Aero.h"

/* Particle tracking functions */
namespace IPT
{
    

    namespace ASCII
    {
        string Get_Variables(int const& offset)
        {
            string variables;
                
            if(offset == 1)
                variables = "\"Y\", \"Z\"";
            else if(offset == 2)
                variables = "\"X\", \"Z\"";
            else if (offset == 3)
                variables = "\"X\", \"Y\"";
            else
                variables = "\"X\", \"Y\", \"Z\"";

            variables += ", \"t\", \"dt\", \"v\", \"a\", \"ptID\", \"Cell_V\", \"Cell_Rho\", \"Cell_ID\"";
            return variables;
        }

        string Cell_Variables(int const& offset)
        {
            string variables;
                
            if(offset == 1)
                variables = "\"Y\", \"Z\"";
            else if(offset == 2)
                variables = "\"X\", \"Z\"";
            else if (offset == 3)
                variables = "\"X\", \"Y\"";
            else
                variables = "\"X\", \"Y\", \"Z\"";

            return variables;
        }


        void Write_variables(ofstream& fp, int const& offset)
        {
            fp << "VARIABLES = " << Get_Variables(offset) << endl;
        }

        void Write_Scatter_Header(ofstream& fp, int const& offset)
        {
            fp << "TITLE =\"IPT particle scatter data\" \n";
            Write_variables(fp, offset);
            fp << "F=POINT" << endl;
            fp << std::left << std::scientific << std::setprecision(6);
        }

        void Write_Streaks_Header(ofstream& fout, int const& offset)
        {
            fout << "TITLE = \"IPT Streaks\"\n";
            Write_variables(fout, offset);
        }

        void Write_Cells_Header(ofstream& fout, int const& offset)
        {
            fout << "TITLE = \"IPT intersecting cells\"\n";

            string variables;
            if(offset == 1)
                variables = "\"Y\", \"Z\"";
            else if(offset == 2)
                variables = "\"X\", \"Z\"";
            else if (offset == 3)
                variables = "\"X\", \"Y\"";
            else
                variables = "\"X\", \"Y\", \"Z\"";

            fout << "VARIABLES= " << variables << endl;
        }
        
        void Init_IPT_Files(SIM& svar)
        {
            string partf, cellf, streakf, surfacef;

            partf = svar.output_prefix + "_IPT_scatter.dat";

            streakf = svar.output_prefix + "_IPT_streaks.dat";

            cellf = svar.output_prefix + "_IPT_cells.dat";

            surfacef = svar.output_prefix + "_surface_impacts.dat";

            if(svar.partout == 1)
            {
                if(svar.restart == 1)
                {
                    svar.partfile.open(partf,std::ios::app);
                    if(!svar.partfile)
                    {
                        cout << "Couldn't open the IPT particle scatter output file." << endl;
                        exit(-1);
                    }
                }
                else
                {
                    svar.partfile.open(partf,std::ios::out);

                    if(svar.partfile.is_open())
                    {
                        Write_Scatter_Header(svar.partfile, svar.offset_axis);
                    }
                    else
                    {
                        cout << "Couldn't open the IPT particle scatter output file." << endl;
                        exit(-1);
                    }
                }
            }

            if(svar.streakout == 1)
            {
                if(svar.restart == 1)
                {
                    svar.streakfile.open(streakf,std::ios::app);
                    if(!svar.streakfile)
                    {
                        cout << "Couldn't open the IPT streaks output file." << endl;
                        exit(-1);
                    }
                }
                else
                {
                    svar.streakfile.open(streakf,std::ios::out);

                    if(svar.streakfile.is_open())
                    {
                        Write_Streaks_Header(svar.streakfile, svar.offset_axis);
                    }
                    else
                    {
                        cout << "Couldn't open the IPT streaks output file." << endl;
                        exit(-1);
                    }
                }
            }

            if(svar.cellsout == 1)
            {
                if(svar.restart == 1)
                {
                    svar.cellfile.open(cellf,std::ios::app);
                    if(!svar.cellfile)
                    {
                        cout << "Couldn't open the IPT cell intersection output file." << endl;
                        exit(-1);
                    }
                }
                else
                {
                    svar.cellfile.open(cellf,std::ios::out);

                    if(svar.cellfile.is_open())
                    {
                        Write_Cells_Header(svar.cellfile, svar.offset_axis);
                    }
                    else
                    {
                        cout << "Couldn't open the IPT cell intersection output file." << endl;
                        exit(-1);
                    }
                }
            }
        }

        void Write_Point(ofstream& fp, real const& scale, IPTPart const& pnp1)
        {
            const static uint width = 15;

            for(uint dim = 0; dim < SIMDIM; ++dim)
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

        void Write_Timestep(ofstream& fp, SIM const& svar, IPTState const& pnp1)
        {
            fp <<  "ZONE T=\"" << "IPT scatter data" << "\"";
            fp <<", I=" << pnp1.size() << ", F=POINT" <<", STRANDID=5, SOLUTIONTIME=" << svar.t  << "\n";
            fp << std::left << std::scientific << std::setprecision(6);
            
            for (size_t ii = 0; ii < pnp1.size(); ++ii)
            {
                Write_Point(fp, svar.scale, pnp1[ii]);
            }
            fp << std::flush;
        }

        void Write_Streaks( SIM& svar, IPTState const& t_pnp1)
        {
            ofstream& fp = svar.streakfile;
            if(!fp.is_open())
            {
                cout << "The streak output file is not open. Cannot write." << endl;
                exit(-1);
            }
            
            size_t nTimes = t_pnp1.size();
            size_t time = 0;
            
            fp <<  "ZONE T=\"" << "Particle " << t_pnp1[0].partID << "\"";
            
            fp <<", I=" << nTimes << ", J=1, K=1, DATAPACKING=POINT"<< "\n";

            fp << std::left << std::scientific << std::setprecision(6);

            for (time = 0; time < nTimes; ++time)
            { /* Inner loop to write the times of the particle */
                Write_Point(fp, svar.scale, t_pnp1[time]);
            }
        }


        void Write_Cells(SIM& svar, MESH const& cells, IPTState const& t_pnp1, int32_t const& totalNumFaceNodes, 
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left, vector<int32_t> const& right)
        {
            ofstream& fout = svar.cellfile;

            if(!fout.is_open())
            {
                cout << "Couldn't open the cell intersection output file." << endl;
                exit(-1);
            }

            fout << "ZONE T=\"particle " << t_pnp1[0].partID << " intersecting cells\"" << endl;
            #if SIMDIM == 3
            fout << "ZONETYPE=FEPOLYHEDRON" << endl;
            #else
            fout << "ZONETYPE=FEPOLYGON" << endl;
            #endif
            fout << "NODES=" << usedVerts.size() << " ELEMENTS=" << cellIndexes.size() << 
                    " FACES=" << faces.size() << endl;
            fout << "TotalNumFaceNodes=" << totalNumFaceNodes << endl;
            fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

            size_t w = 15;
            size_t preci = 6;
            fout << std::left << std::scientific << std::setprecision(preci);

            size_t newl = 0;
            fout << std::setw(1);
            for(size_t DIM = 0; DIM < SIMDIM; ++DIM)
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

        void Write_Impacts(SIM const& svar, FLUID const& fvar, MESH const& cells, 
                    vector<SURF> const& surfaces_to_write, vector<string> const& names, 
                    vector<vector<real>> const& beta_data, vector<vector<StateVecD>> const& usedVerts, 
                    vector<vector<vector<size_t>>> const& faces)
        {
            string file = svar.output_prefix;
            file.append("_surface_impacts.dat");

            ofstream fout(file);

            if(!fout.is_open())
            {
                cout << "Couldn't open the surface impact output file." << endl;
                exit(-1);
            }

            fout << "TITLE=\"Surface collision metrics\"\n";
            #if SIMDIM == 3
            fout << "VARIABLES= \"X\" \"Y\" \"Z\" \"Number of Impacts\" \"beta\" \"average area\"\n";
            #else
            fout << "VARIABLES= \"X\" \"Z\" \"Number of Impacts\" \"beta\" \"average area\"\n";
            #endif
            for(size_t ii = 0; ii < faces.size(); ++ii)
            {
                fout << "ZONE T=\"" << names[ii] << "\"\n";
                fout << "N=" << usedVerts[ii].size() << ", E=" << faces[ii].size();
                #if SIMDIM == 3
                fout << ", F=FEBLOCK ET=QUADRILATERAL, VARLOCATION=([1-3]=NODAL,[4-7]=CELLCENTERED)\n\n";
                #else
                fout << ", F=FELINESEG, VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED)\n\n";
                #endif
                
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

    namespace BINARY
    {
        string Get_Variables(int const& offset)
        {
            string variables;
                
            if(offset == 1)
                variables = "Y,Z";
            else if(offset == 2)
                variables = "X,Z";
            else if (offset == 3)
                variables = "X,Y";
            else
                variables = "X,Y,Z";

            variables += 
            ",t,dt,v,a,ptID,Cell_V,Cell_Rho,Cell_ID";
            return variables;
        }

        string Cell_Variables(int const& offset)
        {
            string variables;
                
            if(offset == 1)
                variables = "Y,Z";
            else if(offset == 2)
                variables = "X,Z";
            else if (offset == 3)
                variables = "X,Y";
            else
                variables = "X,Y,Z";

            return variables;
        }

        void Open_File(string const& file, string const& title, string const& var, void* &fileHandle)
        {
            if(tecFileWriterOpen(file.c_str(),title.c_str(),var.c_str(),1,0,1,NULL,&fileHandle))
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
        }

        void Init_IPT_Files(SIM& svar)
        {
            string partf, cellf, streakf, surfacef;

            partf = svar.output_prefix + "_IPT_scatter.szplt";

            streakf = svar.output_prefix + "_streaks.szplt";

            cellf = svar.output_prefix + "_cells.szplt";

            surfacef = svar.output_prefix + "_surface_impacts.szplt";

            // int32_t fileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
            // int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT
            string part_variables = Get_Variables(svar.offset_axis);
            string cell_variables = Cell_Variables(svar.offset_axis);

            #if FOD == 1
                int32_t realType = 2;
            #else
                int32_t realType = 1;
            #endif

            #if SIMDIM == 3
            vector<int32_t> varTypes = {realType,realType,realType,realType,3,realType,realType,3};
            #else
            vector<int32_t> varTypes = {realType,realType,realType,realType,3,realType,realType,3};
            #endif

            if(svar.partout == 1)
            {
                string title = "IPT particle scatter data";
                Open_File(partf,title,part_variables,svar.partHandle);
            }

            if(svar.streakout == 1)
            {
                string title = "IPT particle streak data";
                Open_File(streakf,title,part_variables,svar.streakHandle);
            }

            if(svar.cellsout == 1)
            {
                string title = "IPT particle cell intersection data";
                Open_File(cellf,title,cell_variables,svar.cellHandle);
            }
        }

        void Write_Point(SIM const& svar, IPTPart const& pnp1)
        {
            int64_t const size = 1;

            double solTime = pnp1.t;     
            int32_t outputZone;

            #if FOD == 1
                int32_t realType = 2;
            #else
                int32_t realType = 1;
            #endif

            #if SIMDIM == 3
            vector<int32_t> varTypes = {realType,realType,realType,realType,3,realType,realType,3};
            #else
            vector<int32_t> varTypes = {realType,realType,realType,3,realType,realType,3};
            #endif
            vector<int32_t> shareVarFromZone(varTypes.size(),0);
            vector<int32_t> valueLocation(varTypes.size(),1);
            vector<int32_t> passiveVarList(varTypes.size(),0);
 
            string group = "IPT Particle " + std::to_string(pnp1.partID) + " scatter data";
            if(tecZoneCreateIJK(svar.partHandle,group.c_str(),size,1,1,&varTypes[0],
                &shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&outputZone))
            {
                cerr << "Failed to create IJK zone." << endl;
                exit(-1);
            }

            if(tecZoneSetUnsteadyOptions(svar.partHandle, outputZone, solTime, 4))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

            vector<real> x(size);
            int32_t var = 1;

            for(uint dim = 0; dim < SIMDIM; ++dim)
            {
                x[0] = pnp1.xi(dim)/svar.scale;

                string name = "position coordinate ";
                name.append(std::to_string(dim));
                Write_Real_Vector(svar.partHandle, outputZone, var, size, x, name);
            }

            vector<real> t(size);
            vector<real> dt(size);
            vector<real> vel(size);
            vector<real> acc(size);
            vector<real> cVel(size);
            vector<real> cRho(size);
            vector<int> pID(size);
            vector<int> cID(size);

            t[0] = pnp1.t;
            dt[0] = pnp1.dt;
            vel[0] = pnp1.v.norm();
            acc[0] = pnp1.acc;
            cVel[0] = pnp1.cellV.norm();
            cRho[0] = pnp1.cellRho;
            pID[0] = pnp1.partID;
            cID[0] = pnp1.cellID;

            Write_Real_Vector(svar.partHandle, outputZone, var, size, t, "particle time");
            Write_Real_Vector(svar.partHandle, outputZone, var, size, dt, "timestep");
            Write_Real_Vector(svar.partHandle, outputZone, var, size, vel, "velocity magnitude");
            Write_Real_Vector(svar.partHandle, outputZone, var, size, acc, "acceleration magnitude");
            Write_Int_Vector(svar.partHandle, outputZone, var, size, pID, "particle ID");
			Write_Real_Vector(svar.partHandle, outputZone, var, size, cVel, "cell velocity");
			Write_Real_Vector(svar.partHandle, outputZone, var, size, cRho, "cell density");
			Write_Int_Vector(svar.partHandle, outputZone, var, size, cID, "cell ID");

            if(tecFileWriterFlush(svar.partHandle,0,NULL))
            {
                cout << "Failed to flush data. Retrying..." << endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                // retry the flush
                if(tecFileWriterFlush(svar.partHandle,0,NULL))
                {
                    cerr << "Failed to flush data to particle file" << endl;
                    exit(-1);
                }
            }
        }

        void Write_State(SIM const& svar, IPTState const& pnp1, string const& zoneName, void* &fileHandle)
        {
            int64_t const size = pnp1.size();

            double solTime = svar.t;     
            int32_t outputZone;

            #if FOD == 1
                int32_t realType = 2;
            #else
                int32_t realType = 1;
            #endif

            #if SIMDIM == 3    
            // variables =  "X,Y,Z,t,dt,v,a,ptID,Cell_V,Cell_Rho,Cell_ID"
            vector<int32_t> varTypes = 
            {realType,realType,realType,realType,realType,realType,realType,3,realType,realType,3};
            #else
            // variables =  "X,Z,t,dt,v,a,ptID,Cell_V,Cell_Rho,Cell_ID"
            vector<int32_t> varTypes = 
                {realType,realType,realType,realType,realType,realType,3,realType,realType,3};
            #endif
            vector<int32_t> shareVarFromZone(varTypes.size(),0);
            vector<int32_t> valueLocation(varTypes.size(),1);
            vector<int32_t> passiveVarList(varTypes.size(),0);
 
            if(tecZoneCreateIJK(fileHandle,zoneName.c_str(),size,1,1,&varTypes[0],
                &shareVarFromZone[0],&valueLocation[0],&passiveVarList[0],0,0,0,&outputZone))
            {
                cerr << "Failed to create IJK zone." << endl;
                exit(-1);
            }

            if(tecZoneSetUnsteadyOptions(fileHandle, outputZone, solTime, 6))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

            vector<real> x(size);
            int32_t var = 1;

            for(uint dim = 0; dim < SIMDIM; ++dim)
            {
                #pragma omp parallel for
                for(uint ii = 0; ii < size; ++ii)
                    x[ii] = pnp1[ii].xi(dim)/svar.scale;

                string name = "position coordinate " + std::to_string(dim);
                Write_Real_Vector(fileHandle, outputZone, var, size, x, name);
            }

            vector<real> t(size);
            vector<real> dt(size);
            vector<real> vel(size);
            vector<real> acc(size);
            vector<real> cVel(size);
            vector<real> cRho(size);
            vector<int> pID(size);
            vector<int> cID(size);

            for (int ii = 0; ii < size; ++ii)
            {
                t[ii] = pnp1[ii].t;
                dt[ii] = pnp1[ii].dt;
                vel[ii] = pnp1[ii].v.norm();
                acc[ii] = pnp1[ii].acc;
                cVel[ii] = pnp1[ii].cellV.norm();
                cRho[ii] = pnp1[ii].cellRho;
                pID[ii] = pnp1[ii].partID;
                cID[ii] = pnp1[ii].cellID;
            }

            Write_Real_Vector(fileHandle, outputZone, var, size, t, "particle time");
            Write_Real_Vector(fileHandle, outputZone, var, size, dt, "timestep");
            Write_Real_Vector(fileHandle, outputZone, var, size, vel, "velocity magnitude");
            Write_Real_Vector(fileHandle, outputZone, var, size, acc, "acceleration magnitude");
            Write_Int_Vector(fileHandle, outputZone, var, size, pID, "particle ID");
			Write_Real_Vector(fileHandle, outputZone, var, size, cVel, "cell velocity");
			Write_Real_Vector(fileHandle, outputZone, var, size, cRho, "cell density");
			Write_Int_Vector(fileHandle, outputZone, var, size, cID, "cell ID");

        }

        void Write_Cells(SIM& svar, MESH const& cells, IPTState const& pnp1, int32_t const& totalNumFaceNodes, 
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left, vector<int32_t> const& right)
        {
            // int64_t const size = pnp1.size();

            double solTime = svar.t;     
            int32_t outputZone;

            #if FOD == 1
                int32_t realType = 2;
            #else
                int32_t realType = 1;
            #endif

            #if SIMDIM == 3
            vector<int32_t> varTypes = {realType,realType,realType,realType,3,realType,realType,3};
            #else
            vector<int32_t> varTypes = {realType,realType,realType,3,realType,realType,3};
            #endif
            vector<int32_t> shareVarFromZone(varTypes.size(),0);
            vector<int32_t> valueLocation(varTypes.size(),1);
            vector<int32_t> passiveVarList(varTypes.size(),0);

            #if SIMDIM == 3
            int32_t zoneType = 7; /* FE Polyhedron */
            #else
            int32_t zoneType = 6; /* FE Polygon */
            #endif
            int32_t nNodes = usedVerts.size();
            int32_t nFaces = faces.size();
            int32_t nCells = cellIndexes.size();

            string group = "IPT cell intersection data";
            if(tecZoneCreatePoly(svar.cellHandle,group.c_str(),zoneType,nNodes,nFaces,nCells,
                totalNumFaceNodes,&varTypes[0],&shareVarFromZone[0],&valueLocation[0],
                &passiveVarList[0],0,0,0,&outputZone))
            {
                cerr << "Failed to create polyhedral/polygonal zone." << endl;
                exit(-1);
            }

            if(tecZoneSetUnsteadyOptions(svar.cellHandle, outputZone, solTime, 7))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

            vector<real> x(nNodes);
            int32_t var = 1;

            for(uint dim = 0; dim < SIMDIM; ++dim)
            {
                #pragma omp parallel for
                for(int32_t ii = 0; ii < nNodes; ++ii)
                    x[ii] = usedVerts[ii][dim];

                string name = "position coordinate " + std::to_string(dim);
                Write_Real_Vector(svar.cellHandle, outputZone, var, nNodes, x, name);
            }

            vector<int32_t> faceCounts(faces.size());
            vector<int32_t> faceNodes(totalNumFaceNodes);
            size_t jj = 0;
            /*Inform of how many vertices in each face*/
            for (size_t ii = 0; ii < faces.size(); ++ii)
            {
                faceCounts[ii] =  faces[ii].size();
                for(auto const& vertex:faces[ii])
                {	/*Write face vertex indexes*/
                    faceNodes[jj] = vertex;
                    jj++;
                }
            }

            tecZoneWritePolyFaces32(svar.cellHandle,outputZone,0,nFaces,&faceCounts[0],
                            &faceNodes[0],&left[0],&right[0],1);
        }
    }

    void Get_Intersection_Data(SIM const& svar, MESH const& cells, IPTState const& t_pnp1,
        vector<StateVecD>& usedVerts, vector<int>& cellIndexes, vector<vector<int32_t>>& faces, 
        vector<int32_t>& left, vector<int32_t>& right,  int32_t& totalNumFaceNodes)
    {
        size_t nTimes = t_pnp1.size();
        size_t time = 0;


        vector<int32_t> vertIndexes;
        vector<int32_t> faceIndexes;

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
            faces.emplace_back(cells.faces[face].begin(),cells.faces[face].end());
            left.emplace_back(cells.leftright[face].first);
            right.emplace_back(cells.leftright[face].second);
        }

        /* Go through the faces indexes, finding the appropriate index to change it to */
        for(vector<int32_t>& face : faces)
        {
            for(int32_t& vert : face)
            {
                vector<int32_t>::iterator index = std::find(vertIndexes.begin(), vertIndexes.end(), vert);
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
            totalNumFaceNodes += face.size();
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
    }

    void Get_Impact_Data(SIM const& svar, FLUID const& fvar, MESH const& cells, 
                vector<SURF> const& surfs, vector<vector<real>> const& beta_data,
                vector<SURF>& surfaces_to_write, vector<string>& names, 
                vector<vector<StateVecD>>& usedVerts, vector<vector<vector<size_t>>>& faces)
    {
        // uint nSurf = svar.markers.size();
        uint nSurf = 0;
        vector<int> marks;
        // vector<string> names;
        // vector<SURF> surfaces_to_write;
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

        faces =  vector<vector<vector<size_t>>>(nSurf);
        vector<std::pair<size_t,int>> smarkers = cells.smarkers;

        vector<vector<size_t>> vertIndexes(nSurf);
        usedVerts = vector<vector<StateVecD>>(nSurf);
        
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
    }  
    
    void Init_IPT_Files(SIM& svar)
    {
        if(svar.out_encoding == 0)
            ASCII::Init_IPT_Files(svar);
        else
            BINARY::Init_IPT_Files(svar);
    }

    void Write_Point(SIM& svar, IPTPart& pnp1)
    {
        if(svar.out_encoding == 0)
            ASCII::Write_Point(svar.partfile, svar.scale, pnp1);
        else
            BINARY::Write_Point(svar,pnp1);

    }

    void Write_Data(SIM& svar, MESH& cells, vector<IPTState>& time_record)
    {
        if(svar.partout == 1)
        {
            // for(size_t ii = 0; ii < time_record.size(); ii++)
            // {
            //     if(svar.out_encoding == 0)
            //         ASCII::Write_Timestep(svar.partfile,svar,time_record[ii]);
            //     else
            //         BINARY::
            // }
        }
        if(svar.streakout == 1)
        {
            for(size_t ii = 0; ii < time_record.size(); ii++)
            {
                if(!time_record[ii].empty())
                {
                    if(svar.out_encoding == 0)
                        ASCII::Write_Streaks(svar,time_record[ii]);
                    else
                    {
                        string zone = "IPT Particle" + std::to_string(time_record[ii][0].partID) + " Streak";
                        BINARY::Write_State(svar,time_record[ii],zone,svar.streakHandle);
                    }
                }
            }

            if(!time_record.empty())
            {
                if(svar.out_encoding == 0)
                    svar.streakfile << std::flush;
                else
                {
                
                    if(tecFileWriterFlush(svar.streakHandle,0,NULL))
                    {
                        cout << "Failed to flush data. Retrying..." << endl;
                        std::this_thread::sleep_for(std::chrono::milliseconds(10));
                        // retry the flush
                        if(tecFileWriterFlush(svar.streakHandle,0,NULL))
                        {
                            cerr << "Failed to flush data to streak file" << endl;
                            exit(-1);
                        }
                    }
                }
            }
        }
        if(svar.cellsout == 1)
        {
            for(size_t ii = 0; ii < time_record.size(); ii++)
            {
                if(!time_record[ii].empty())
                {
                    vector<StateVecD> usedVerts;
                    vector<int32_t> cellIndexes, left, right;
                    vector<vector<int32_t>> faces;
                    int32_t totalNumFaceNodes = 0;
                    Get_Intersection_Data(svar, cells, time_record[ii], usedVerts, cellIndexes, faces, 
                                        left, right, totalNumFaceNodes);
                    if(svar.out_encoding == 0)
                    {
                        ASCII::Write_Cells(svar,cells,time_record[ii],totalNumFaceNodes,usedVerts,
                                        cellIndexes,faces,left, right);
                    }
                    else
                    {
                        BINARY::Write_Cells(svar,cells,time_record[ii],totalNumFaceNodes,usedVerts,
                                        cellIndexes,faces,left, right);
                    }
                }
            }

            if(!time_record.empty())
            {
                if(svar.out_encoding == 0)
                    svar.cellfile << std::flush;
                else
                {
                    if(tecFileWriterFlush(svar.cellHandle,0,NULL))
                    {
                        cout << "Failed to flush data. Retrying..." << endl;
                        std::this_thread::sleep_for(std::chrono::milliseconds(10));
                        // retry the flush
                        if(tecFileWriterFlush(svar.cellHandle,0,NULL))
                        {
                            cerr << "Failed to flush data to cell intersection file" << endl;
                            exit(-1);
                        }
                    }
                }
            }
        }

        if(!time_record.empty())
            time_record.clear();
    } 

    void Terminate_Particle(SIM& svar, MESH const& cells, vector<IPTPart>& time_record, size_t const& ii,
        IPTPart& pnp1, vector<SURF>& surface_data, vector<IPTState>& iptdata)
    {
        pnp1.going = 0;

        if(pnp1.failed != 1)
            time_record.emplace_back(pnp1);

        // cout << "Particle " << pnp1.partID << " stopped iterating"; 
        if(pnp1.cellID < 0)
        {
            /* Find which surface it has impacted  */
            size_t const& face = pnp1.faceID;

            auto const index = std::find_if(cells.smarkers.begin(), cells.smarkers.end(), 
            [face](std::pair<size_t,int> const& p1){return p1.first == face;});

            if(index != cells.smarkers.end())
            {
                /* Find which surface to add to, and which face of that*/
                auto const surface = std::find(svar.markers.begin(), svar.markers.end(), index->second);

                if(surface != svar.markers.end())
                {
                    size_t const surf = surface - svar.markers.begin();

                    // cout << " after hitting surface: " << svar.bnames[surf] << endl;
                    #pragma omp critical
                    {
                        surface_data[surf].marker_count++;
                        surface_data[surf].marker_pIDs.emplace_back(ii);
                        surface_data[surf].end_pos.emplace_back(pnp1.xi);
                        surface_data[surf].start_pos.emplace_back(time_record[0].xi);
                    }
                    auto const index2 =
                        std::find(surface_data[surf].faceIDs.begin(), surface_data[surf].faceIDs.end(), face);

                    if(index2 != surface_data[surf].faceIDs.end())
                    {
                        size_t fIndex = index2 - surface_data[surf].faceIDs.begin();
                        #pragma omp critical
                        {
                            surface_data[surf].face_count[fIndex]++;
                            surface_data[surf].impacted_face.emplace_back(fIndex);
                        }
                    }
                    else
                    {
                        cout << "Couldn't find which face of the surface the particle impacted" << endl;
                    }
                }
                else
                {
                    cout << "Couldn't find which surface index for the marker." << endl;
                }
                
            }
            else
            {
                cout << "Couldn't find which surface the particle impacted" << endl;
            }

            
        }
        else
        {
            // cout << " prematurely" << endl;
        }

        iptdata.emplace_back(time_record);
        // if(svar.streakout == 1)
        // {
        //     #pragma omp critical
        //     {
        //         Write_ASCII_Streaks(svar,time_record,pnp1);
        //     }
        // }

        // if(svar.cellsout == 1)
        // {
        //     #pragma omp critical
        //     {
        //         Write_ASCII_Cells(svar,cells,time_record);
        //     }
        // }
    }


    inline real AeroForce(StateVecD const& Vdiff, FLUID const& fvar, AERO const& avar, IPTPart const& pi)
    {
        real const Re = 2.0*pi.faceRho*Vdiff.norm()*pi.d/avar.mug;
        real const Cd = GetCd(Re);


        return Vdiff.norm() * (3.0*Cd*pi.faceRho)/(4.0 * pi.d * fvar.rho0);
    }

    inline void BFD1(FLUID const& fvar, AERO const& avar, StateVecD const& g, IPTPart const& pn, IPTPart& pnp1)
    {
        StateVecD const Vdiff = pnp1.faceV - pnp1.v ;

        real res = AeroForce(Vdiff, fvar, avar, pnp1);
        
        pnp1.acc = res;
        pnp1.v = (pn.v + pnp1.dt*res*pnp1.faceV + pnp1.dt*g)/(1.0 + pnp1.dt*res);
        pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); 
    }

    inline void BFD2(FLUID const& fvar, AERO const& avar, StateVecD const& g, real const& relax, real const& dt, real const& dtm1, 
                            IPTPart const& pnm1, IPTPart const& pn, IPTPart& pnp1)
    {
        StateVecD const Vdiff = pnp1.faceV - pnp1.v ;

        real res = AeroForce(Vdiff, fvar, avar, pnp1);

        pnp1.acc = res;
        /* Second order velocity calculation */
        pnp1.v = ((dt+dtm1)*(dt*dtm1*(relax*res*pnp1.faceV+g) + (dt+dtm1)* pn.v) - dt*dt*pnm1.v)
                    /(dtm1*((2*dt+dtm1)+dt*(dt+dtm1)*res));

                    
        // pnp1.xi = pn.xi-pnp1.dt*(pnp1.v + pn.v)+0.5*pn.dt*(pn.v + pnm1.v); 
        pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); /* Just a first order space integration. */
    }

    void Integrate(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, size_t const& ii, 
                            IPTPart& pnm1, IPTPart& pn, IPTPart& pnp1, vector<SURF>& marker_data, 
                            vector<IPTState>& iptdata)
    {
        /* Do the first step */
        IPTState time_record;
        time_record.emplace_back(pnp1);
        ofstream partfile;

        /* Find intersecting face for np1 */
        // pnp1.faceV = pnp1.cellV; /* Start by guessing the face as the cell*/
        // pnp1.faceRho = pnp1.cellRho;
        FindFace(svar, cells, pn, pnp1);
        /* Now c+1 is known, use the central difference face properties*/
        pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
        pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

        // int iter = 0;
        // int min_iter = 3;
        // real error = 1.0;
        real relaxation = 1.0; /* Relaxation factor for 2nd order */
        for(size_t iter = 0; iter < 10; iter++)
        // while(error > -7.0 && iter < min_iter)
        {
            relaxation = (real(iter)+1.0)/6.0;
            if (relaxation > 1.0)
                relaxation = 1.0;
            // cout << "Iteration: " << iter << endl; 
            // vec<real,3> temp = pnp1.xi;

            if (svar.eqOrder == 2)
            {
                real const dt = pnp1.dt;
                real const dtm1 = pnp1.dt;
                BFD2(fvar, avar, svar.grav, relaxation, dt, dtm1, pnm1, pn, pnp1);
            }
            else
                BFD1(fvar, avar, svar.grav, pn, pnp1);

            FindFace(svar, cells, pn, pnp1);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            // error = log10((pnp1.xi-temp).norm());

            // if(iter < svar.max_iter)
            //     iter++;
            // else
            //     break;
        }

        if(svar.partout == 1)
            Write_Point(svar, pnp1);

        if (pnp1.going == 0 || pnp1.v.norm() > 1e3
            || (pnp1.xi - pn.xi).norm() > cells.maxlength )
        {
            pnp1.failed = 1;
            Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data, iptdata);
            #pragma omp atomic
                svar.nFailed++;
        }
        else if (pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x)
        {      
            Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data, iptdata);
            #pragma omp atomic
                svar.nSuccess++;
        }

        /* Move forward in time */
        // cout << "Time: " << pnp1.t << " dt: " << pnp1.dt  << " x: " <<
        // pnp1.xi[0] << " y: " << pnp1.xi[1] << " z: " << pnp1.xi[2] << 
        // " old cell: " << pn.cellID << " new cell: " << pnp1.cellID << endl;

        /* March forward in time */
        if(svar.eqOrder == 2)
            pnm1 = pn;

        pn = pnp1;
        pnp1.t += pnp1.dt;
        

        /* Only need the information if streaks or cells are output */
        if(svar.cellsout == 1 || svar.streakout == 1)
            time_record.emplace_back(pnp1);

        while(pnp1.going != 0)
        {
        
            /* Find intersecting face for np1 */
            // pnp1.faceV = pnp1.cellV; /* Start by guessing the face as the cell*/
            // pnp1.faceRho = pnp1.cellRho;
            
            if(pnp1.d < 1.0e-5)
            {   /* If particle is small, use cell velocity to start */
                pnp1.v = pnp1.cellV;
            }
            else if (pnp1.d < 100.0e-6)
            {   /* If particle is heavier, use sliding scale up to 100 microns */
                pnp1.v = (1.0 - (pnp1.d-1e-05)/(9e-5)) * pnp1.cellV + 
                (pnp1.d-1e-05)/(9e-5) * pnp1.v;
            }


            FindFace(svar, cells, pn, pnp1);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            // int iter = 0;
            // int min_iter = 3;
            // real error = 1.0;
            real relaxation = 1.0;
            for(size_t iter = 0; iter < 10; iter++)
            // while(error > -7.0 && iter < min_iter)
            {
                relaxation = (real(iter)+1.0)/6.0;
                if (relaxation > 1.0)
                    relaxation = 1.0;

                // cout << "Iteration: " << iter << endl; 
                // vec<real,3> temp = pnp1.xi;
                
                if(svar.eqOrder == 2)
                {
                    real const dt = pnp1.dt;
                    real dtm1 = pnp1.dt;
                    if(iter > 0)
                        dtm1 = pn.dt;
                    BFD2(fvar, avar, svar.grav, relaxation, dt, dtm1, pnm1, pn, pnp1);
                }
                else
                    BFD1(fvar, avar, svar.grav, pn, pnp1);


                FindFace(svar, cells, pn, pnp1);
                /* Now c+1 is known, use the central difference face properties*/
                pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
                pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

                // error = log10((pnp1.xi-temp).norm());

                // if(iter < svar.max_iter)
                //     iter++;
                // else
                //     break;
            }

            if(svar.partout == 1)
                Write_Point(svar, pnp1);

            if ( pnp1.going == 0 || pnp1.v.norm() > 1e3
                    || (pnp1.xi - pn.xi).norm() > cells.maxlength )
            {
                if((pnp1.xi - pn.xi).norm() > cells.maxlength)
                {
                    cout << "Particle terminated because of max length" << endl;
                    cout << (pnp1.xi - pn.xi).norm() << "  " << cells.maxlength << endl;
                }
                pnp1.failed = 1;
                Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data, iptdata);
                #pragma omp atomic
                    svar.nFailed++;
                break;
            }
            else if ( pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x)
            {    
                cout << "Particle terminated because of a success" << endl;
                Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data, iptdata);
                #pragma omp atomic
                    svar.nSuccess++;
                break;
            }

            // /* Move forward in time */
            // cout << "Time: " << pnp1.t << " dt: " << pnp1.dt  << " x: " <<
            // pnp1.xi[0] << " y: " << pnp1.xi[1] << " z: " << pnp1.xi[2] << 
            // " old cell: " << pn.cellID << " new cell: " << pnp1.cellID << endl;

            /* March forward in time */
            if(svar.eqOrder == 2)
                pnm1 = pn;
            
            pn = pnp1;
            pnp1.t += pnp1.dt;
            
            /* Only need the information if streaks or cells are output */
            if(svar.cellsout == 1 || svar.streakout == 1)
                time_record.emplace_back(pnp1);

        }
        
    }
}
#endif