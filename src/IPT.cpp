/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/
#include "IPT.h"
#include "Aero.h"
#include "BinaryIO.h"
#include "Containment.h"
#include "IOFunctions.h"
#include <unordered_set>
#include <thread>

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


        void Write_variables(FILE* fout, int const& offset)
        {
            fprintf(fout,"VARIABLES = %s", Get_Variables(offset).c_str());
        }

        void Write_Scatter_Header(FILE* fout, int const& offset)
        {
            fprintf(fout, "TITLE =\"IPT particle scatter data\" \n");
            Write_variables(fout, offset);
            fprintf(fout,"F=POINT\n");
        }

        void Write_Streaks_Header(FILE* fout, int const& offset)
        {
            fprintf(fout,"TITLE = \"IPT Streaks\"\n");
            Write_variables(fout, offset);
        }

        void Write_Cells_Header(FILE* fout, int const& offset)
        {
            fprintf(fout,"TITLE = \"IPT intersecting cells\"\n");

            string variables;
            if(offset == 1)
                variables = "\"Y\", \"Z\"";
            else if(offset == 2)
                variables = "\"X\", \"Z\"";
            else if (offset == 3)
                variables = "\"X\", \"Y\"";
            else
                variables = "\"X\", \"Y\", \"Z\"";

            fprintf(fout,"VARIABLES= %s\n", variables.c_str());
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
                    svar.partfile = fopen(partf.c_str(),"a");
                    if(svar.partfile == NULL)
                    {
                        cout << "Couldn't open the IPT particle scatter output file." << endl;
                        exit(-1);
                    }
                }
                else
                {
                    svar.partfile = fopen(partf.c_str(),"w");

                    if(svar.partfile == NULL)
                    {
                        cout << "Couldn't open the IPT particle scatter output file." << endl;
                        exit(-1);
                    }
                    else
                    {
                        Write_Scatter_Header(svar.partfile, svar.offset_axis);
                    }
                }
            }

            if(svar.streakout == 1)
            {
                if(svar.restart == 1)
                {
                    svar.streakfile = fopen(streakf.c_str(),"a");
                    if(svar.streakfile == NULL)
                    {
                        cout << "Couldn't open the IPT streaks output file." << endl;
                        exit(-1);
                    }
                }
                else
                {
                    svar.streakfile = fopen(streakf.c_str(),"w");

                    if(svar.streakfile == NULL)
                    {
                        cout << "Couldn't open the IPT streaks output file." << endl;
                        exit(-1);
                    }
                    else
                    {
                        Write_Streaks_Header(svar.streakfile, svar.offset_axis);
                    }
                }
            }

            if(svar.cellsout == 1)
            {
                if(svar.restart == 1)
                {
                    svar.cellfile = fopen(cellf.c_str(),"a");
                    if(svar.cellfile == NULL)
                    {
                        cout << "Couldn't open the IPT cell intersection output file." << endl;
                        exit(-1);
                    }
                }
                else
                {
                    svar.cellfile = fopen(cellf.c_str(),"w");

                    if(svar.cellfile == NULL)
                    {
                        cout << "Couldn't open the IPT cell intersection output file." << endl;
                        exit(-1);
                    }
                    else
                    {
                        Write_Cells_Header(svar.cellfile, svar.offset_axis);
                    }
                }
            }
        }

        void Write_Point(FILE* fout, real const& scale, IPTPart const& pnp1)
        {
            for(uint dim = 0; dim < SIMDIM; ++dim)
                fprintf(fout," %2.7e", pnp1.xi[dim]/scale);

            fprintf(fout," %2.7e %2.7e %2.7e %2.7e %zu %2.7e %2.7e %6ld\n",
                pnp1.t,pnp1.dt,pnp1.v.norm(),pnp1.acc,pnp1.partID,pnp1.cellV.norm(),pnp1.cellRho,pnp1.cellID);
        }

        void Write_Timestep(FILE* fout, SIM const& svar, IPTState const& pnp1)
        {
            fprintf(fout,"ZONE T=\"IPT scatter data\"");
            fprintf(fout,", I=%zu, F=POINT, STRANDID=5, SOLUTIONTIME= %.7g\n", pnp1.size(), svar.t);

            for (size_t ii = 0; ii < pnp1.size(); ++ii)
                Write_Point(fout, svar.scale, pnp1[ii]);

            fflush(fout);
        }

        void Write_Streaks(SIM& svar, IPTState const& t_pnp1)
        {
            FILE* fout = svar.streakfile;
            if(fout == NULL)
            {
                cout << "The streak output file is not open. Cannot write." << endl;
                exit(-1);
            }
            
            size_t nTimes = t_pnp1.size();
            size_t time = 0;
            
            fprintf(fout, "ZONE T=\"Particle %zu\"\n",t_pnp1[0].partID);
            fprintf(fout,"I= %zu, J=1, K=1, DATAPACKING=POINT\n",nTimes);

            for (time = 0; time < nTimes; ++time)
            { /* Inner loop to write the times of the particle */
                Write_Point(fout, svar.scale, t_pnp1[time]);
            }
        }


        void Write_Cells(SIM& svar, MESH const& cells, IPTState const& t_pnp1, int32_t const& totalNumFaceNodes, 
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left, vector<int32_t> const& right)
        {
            FILE* fout = svar.cellfile;

            if(fout == NULL)
            {
                cout << "Couldn't open the cell intersection output file." << endl;
                exit(-1);
            }

            fprintf(fout,"ZONE T=\"particle %zu intersecting cells\"\n", t_pnp1[0].partID);
            #if SIMDIM == 3
            fprintf(fout,"ZONETYPE=FEPOLYHEDRON,");
            #else
            fprintf(fout,"ZONETYPE=FEPOLYGON,");
            #endif
            fprintf(fout," NODES=%zu, ELEMENTS=%zu, FACES=%zu, TotalNumFaceNodes=%d\n", 
                usedVerts.size(), cellIndexes.size(), faces.size(),totalNumFaceNodes);
            fprintf(fout, "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0, DATAPACKING=BLOCK\n");

            size_t newl = 0;
            for(size_t DIM = 0; DIM < SIMDIM; ++DIM)
            {
                for(size_t ii = 0; ii < usedVerts.size(); ++ii)
                {
                    fprintf(fout," %2.7e", usedVerts[ii][DIM]);

                    if(newl > 4)
                    {
                        fprintf(fout,"\n");
                        newl=0;
                    }
                    else
                        newl++;
                }
            }
            fprintf(fout,"\n");


            /*Inform of how many vertices in each face*/
            fprintf(fout,"#node count per face\n");
            newl = 0;
            for (size_t ii = 0; ii < faces.size(); ++ii)
            {
                fprintf(fout," %zu", faces[ii].size());

                if(newl > 8)
                {
                    fprintf(fout,"\n");
                    newl=0;
                }
                else
                    newl++;
            }
            fprintf(fout,"\n");

            /*Write the face data*/
            fprintf(fout,"#face nodes\n");
            for (size_t ii = 0; ii < faces.size(); ++ii)
            {
                for(auto const& vertex:faces[ii])
                {	/*Write face vertex indexes*/
                    fprintf(fout," %6d",vertex);
                }

                fprintf(fout,"\n");
            }
            fprintf(fout,"\n");

            /*Write face left and right*/
            newl = 0;
            fprintf(fout, "#left elements\n");
            for (size_t ii = 0; ii < left.size(); ++ii)
            {
                fprintf(fout," %6d", left[ii]);

                if(newl > 6)
                {
                    fprintf(fout,"\n");
                    newl=0;
                }
                else
                    newl++;
            }
            fprintf(fout,"\n");

            newl = 0;
            fprintf(fout, "#right elements\n");
            for (size_t ii = 0; ii < right.size(); ++ii)
            {
                fprintf(fout," %6d", right[ii]);

                if(newl > 6)
                {
                    fprintf(fout,"\n");
                    newl=0;
                }
                else
                    newl++;
            }

            fclose(fout);
        }

        void Write_Impacts(SIM const& svar, FLUID const& fvar, MESH const& cells, 
                    vector<SURF> const& surfaces_to_write, vector<string> const& names, 
                    vector<vector<real>> const& beta_data, vector<vector<StateVecD>> const& usedVerts, 
                    vector<vector<vector<size_t>>> const& faces)
        {
            string file = svar.output_prefix;
            file.append("_surface_impacts.dat");

            FILE* fout = fopen(file.c_str(), "w");

            if(fout == NULL)
            {
                cout << "Couldn't open the surface impact output file." << endl;
                exit(-1);
            }

            fprintf(fout,"TITLE=\"Surface collision metrics\"\n");
            #if SIMDIM == 3
            fprintf(fout,"VARIABLES= \"X\" \"Y\" \"Z\" \"Number of Impacts\" \"beta\" \"average area\"\n");
            #else
            fprintf(fout,"VARIABLES= \"X\" \"Z\" \"Number of Impacts\" \"beta\" \"average area\"\n");
            #endif
            for(size_t ii = 0; ii < faces.size(); ++ii)
            {
                fprintf(fout,"ZONE T=\"%s\"\n", names[ii].c_str());
                fprintf(fout,"N=%zu, E=%zu, ", usedVerts[ii].size(), faces[ii].size());
                #if SIMDIM == 3
                fprintf(fout,", F=FEBLOCK ET=QUADRILATERAL, VARLOCATION=([1-3]=NODAL,[4-7]=CELLCENTERED)\n\n");
                #else
                fprintf(fout,", F=FELINESEG, VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED)\n\n");
                #endif
     

                size_t newl = 0;
                for(size_t DIM = 0; DIM < SIMDIM; ++DIM)
                {
                    for(size_t jj = 0; jj < usedVerts[ii].size(); ++jj)
                    {
                        fprintf(fout,"%3.7e",usedVerts[ii][jj][DIM]);
                        if(newl > 4)
                        {
                            fprintf(fout,"\n");
                            newl=0;
                        }
                        else
                            newl++;  
                    }
                }
                fprintf(fout,"\n\n"); 

                /* Variable data goes here */
                fprintf(fout, "#face impact count\n");
                for(size_t jj = 0; jj < surfaces_to_write[ii].face_count.size(); ++jj)
                {
                    fprintf(fout," %8u",surfaces_to_write[ii].face_count[jj]);
                    if(newl > 6)
                    {
                        fprintf(fout,"\n");
                        newl=0;
                    }
                    else
                        newl++;
                }
                fprintf(fout,"\n\n"); 

                fprintf(fout,"#face beta value\n");
                for(size_t jj = 0; jj < surfaces_to_write[ii].face_beta.size(); ++jj)
                {
                    fprintf(fout," %2.7e",surfaces_to_write[ii].face_beta[jj]);
                    if(newl > 4)
                    {
                        fprintf(fout,"\n");
                        newl=0;
                    }
                    else
                        newl++;
                }
                fprintf(fout,"\n\n"); 

                fprintf(fout,"#face area value\n");
                for(size_t jj = 0; jj < surfaces_to_write[ii].face_area.size(); ++jj)
                {
                    fprintf(fout," %2.7e",surfaces_to_write[ii].face_area[jj]);
                    if(newl > 4)
                    {
                        fprintf(fout,"\n");
                        newl=0;
                    }
                    else
                        newl++;
                }
                fprintf(fout,"\n\n"); 


                /*Write the face data*/
                fprintf(fout,"#face nodes\n");
                for (size_t jj = 0; jj < faces[ii].size(); ++jj)
                {
                    for(auto const& vertex:faces[ii][jj])
                    {	/*Write face vertex indexes*/
                        fprintf(fout," %zu",vertex);
                    }

                    if(faces[ii][jj].size() == 3)
                        fprintf(fout," %zu",faces[ii][jj].back());

                    fprintf(fout,"\n");
                }
            }
            fclose(fout);
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

        void Open_File(string const& file, string const& title, string const& var, int32_t const& fileType, void* &fileHandle)
        {
            if(tecFileWriterOpen(file.c_str(),title.c_str(),var.c_str(),fileType,0,1,NULL,&fileHandle))
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

            cellf = svar.output_prefix + "_cells.plt";

            surfacef = svar.output_prefix + "_surface_impacts.szplt";

            // int32_t fileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
            // int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT
            string part_variables = Get_Variables(svar.offset_axis);
            string cell_variables = Cell_Variables(svar.offset_axis);

            if(svar.partout == 1)
            {
                string title = "IPT particle scatter data";
                Open_File(partf,title,part_variables,1,svar.partHandle);
            }

            if(svar.streakout == 1)
            {
                string title = "IPT particle streak data";
                Open_File(streakf,title,part_variables,1,svar.streakHandle);
            }

            if(svar.cellsout == 1)
            {
                string title = "IPT particle cell intersection data";
                Open_File(cellf,title,cell_variables,0,svar.cellHandle);
            }
        }

        void Write_Point(SIM const& svar, IPTPart const& pnp1)
        {
            int64_t const size = 1;

            double solTime = pnp1.t;     
            int32_t outputZone;

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
                Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, x, name);
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

            Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, t, "particle time");
            Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, dt, "timestep");
            Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, vel, "velocity magnitude");
            Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, acc, "acceleration magnitude");
            Write_Int_Vector(svar.partHandle, outputZone, var, 0, size, pID, "particle ID");
			Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, cVel, "cell velocity");
			Write_Real_Vector(svar.partHandle, outputZone, var, 0, size, cRho, "cell density");
			Write_Int_Vector(svar.partHandle, outputZone, var, 0, size, cID, "cell ID");

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
                Write_Real_Vector(fileHandle, outputZone, var, 0, size, x, name);
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

            Write_Real_Vector(fileHandle, outputZone, var, 0, size, t, "particle time");
            Write_Real_Vector(fileHandle, outputZone, var, 0, size, dt, "timestep");
            Write_Real_Vector(fileHandle, outputZone, var, 0, size, vel, "velocity magnitude");
            Write_Real_Vector(fileHandle, outputZone, var, 0, size, acc, "acceleration magnitude");
            Write_Int_Vector(fileHandle, outputZone, var, 0, size, pID, "particle ID");
			Write_Real_Vector(fileHandle, outputZone, var, 0, size, cVel, "cell velocity");
			Write_Real_Vector(fileHandle, outputZone, var, 0, size, cRho, "cell density");
			Write_Int_Vector(fileHandle, outputZone, var, 0, size, cID, "cell ID");

        }

        void Write_Cells(SIM& svar, MESH const& cells, IPTState const& pnp1, int32_t const& totalNumFaceNodes, 
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left, vector<int32_t> const& right)
        {
            // int64_t const size = pnp1.size();

            double solTime = svar.t;     
            int32_t outputZone;

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
                Write_Real_Vector(svar.cellHandle, outputZone, var, 0, nNodes, x, name);
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


        std::unordered_set<int32_t> vertIndexes;
        std::unordered_set<int32_t> faceIndexes;

        /* Get cell properties for the intersected cells (delete duplicates later)*/
        for (time = 0; time < nTimes; ++time)
        {
            int cellID = t_pnp1[time].cellID;
            
            /* Find how many vertices the cell has */
            for (size_t const& face: cells.cFaces[cellID])
                for (size_t const& vert : cells.faces[face])
                    vertIndexes.insert(static_cast<int32_t>(vert));

            for (size_t const& face : cells.cFaces[cellID])
                faceIndexes.insert(static_cast<int32_t>(face));

            cellIndexes.emplace_back(cellID);
        }


        for(size_t const vert : vertIndexes)
        {
            usedVerts.emplace_back(cells.verts[vert]);
        }

        /* Get faces properties used in the cell */
        for(size_t const face: faceIndexes )
        {
            faces.emplace_back();
            for (size_t const& vert : cells.faces[face])
                faces.back().emplace_back(static_cast<int32_t>(vert));
            left.emplace_back(cells.leftright[face].first);
            right.emplace_back(cells.leftright[face].second);
        }

        /* Go through the faces indexes, finding the appropriate index to change it to */
        for(vector<int32_t>& face : faces)
        {
            for(int32_t& vert : face)
            {
                std::unordered_set<int32_t>::iterator index = std::find(vertIndexes.begin(), vertIndexes.end(), vert);
                if(index != vertIndexes.end())
                {
                    vert = std::distance(vertIndexes.begin(),index) + 1 ;
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
                    fflush(svar.streakfile);
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
                    fflush(svar.cellfile);
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
                    #ifdef DEBUG
                    else
                    {
                        cout << "Couldn't find which face of the surface the particle impacted" << endl;
                    }
                    #endif
                }
                else
                {
                    cout << "Couldn't find which surface index for the marker." << endl;
                }
                
            }
            #ifdef DEBUG
            else
            {
                cout << "Couldn't find which surface the particle impacted" << endl;
            }
            #endif

            
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

    inline void BFD1(FLUID const& fvar, AERO const& avar, StateVecD const& g, real const& dt, 
        IPTPart const& pn, IPTPart& pnp1)
    {
        StateVecD const Vdiff = pnp1.faceV - pnp1.v ;

        real res = AeroForce(Vdiff, fvar, avar, pnp1);
        
        pnp1.acc = res;
        pnp1.v = (pn.v + dt*res*pnp1.faceV + dt*g)/(1.0 + dt*res);
        pnp1.xi = pn.xi+0.5*dt*(pnp1.v + pn.v); 
    }

    inline void BFD2(FLUID const& fvar, AERO const& avar, StateVecD const& g, real const& dt, real const& dtm1, 
                            IPTPart const& pnm1, IPTPart const& pn, IPTPart& pnp1)
    {
        StateVecD const Vdiff = pnp1.faceV - pnp1.v ;

        real res = AeroForce(Vdiff, fvar, avar, pnp1);

        pnp1.acc = res;
        /* Second order velocity calculation */
        pnp1.v = ((dt+dtm1)*(dt*dtm1*(res*pnp1.faceV+g) + (dt+dtm1)* pn.v) - dt*dt*pnm1.v)
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
        pnp1.dt = FindFace(svar, cells, pn, pnp1);
        /* Now c+1 is known, use the central difference face properties*/
        pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
        pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);


        #ifdef DEBUG
        size_t step = 0;
        ofstream timelog("ipt_timesteps.dat",std::ios::out);
        timelog << "TITLE=\"IPT Timestep logs\"" << endl;
        timelog << "VARIABLES = \"Iteration\", \"timestep\", \"Acceleration\", \"Error\"" << endl;
        timelog << "Zone T=\"Step " << step << "\"" << endl;
        timelog << std::left << std::scientific << std::setprecision(7);
        size_t const width = 16;
        #endif

        size_t iter = 0;
        size_t min_iter = 3;
        size_t max_iter = svar.subits;
        real error = 1.0;
        
        // real relaxation = 1.0; /* Relaxation factor for 2nd order */
        pnp1.relax = 0.4; // Relaxing the timestep over the iteractions works well. 
        // for(size_t iter = 0; iter < 100; iter++)
        while((error > -7.0 || iter < min_iter) && iter < max_iter)
        {
            // Try lower relaxation first then increase to a higher amount to avoid oscillation
            // Doesn't improve convergence so leave it
            // real fac = iter > 5 ? 1 : real(iter)/5.0;
            // pnp1.relax = 0.75 *(1.0-fac) + 0.4*fac;
            // relaxation = (real(iter)+1.0)/20.0;
            // relaxation = 0.4;
            // if (relaxation > 1.0)
            //     relaxation = 1.0;
            // cout << "Iteration: " << iter << endl; 

            real const dt = pnp1.dt;
            if (svar.eqOrder == 2)
            {
                real const dtm1 = pnp1.dt;
                BFD2(fvar, avar, svar.grav, dt, dtm1, pnm1, pn, pnp1);
            }
            else
                BFD1(fvar, avar, svar.grav, dt, pn, pnp1);

            real dt_temp = FindFace(svar, cells, pn, pnp1);
            pnp1.dt = (1.0 - svar.relax) * dt_temp + svar.relax * pnp1.dt;
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            error = log10(fabs(pnp1.dt-dt));
            #ifdef DEBUG
            cout << iter << "  "  << error << endl;
            timelog << setw(width) << iter << setw(width) << pnp1.dt << setw(width) 
                << pnp1.acc << setw(width) << error << endl;
            #endif

            iter++;
        }
        #ifdef DEBUG
        step++;
        #endif

        if(svar.partout == 1)
            Write_Point(svar, pnp1);

        if (pnp1.going == 0 || pnp1.v.norm() > 1e3
            || (pnp1.xi - pn.xi).norm() > cells.maxlength )
        {
            #ifdef DEBUG
            if((pnp1.xi - pn.xi).norm() > cells.maxlength)
            {
                cout << "Particle terminated because of max length" << endl;
                cout << (pnp1.xi - pn.xi).norm() << "  " << cells.maxlength << endl;
            }
            #endif
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
            #ifdef DEBUG
            timelog << "Zone T=\"Step " << step << "\"" << endl;
            #endif
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

            pnp1.dt = FindFace(svar, cells, pn, pnp1);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            // int iter = 0;
            // int min_iter = 3;
            // real error = 1.0;
            iter = 0;
            error = 1.0;
            // for(size_t iter = 0; iter < 10; iter++)
            while((error > -7.0 || iter < min_iter) && iter < max_iter)
            {
                // relaxation = (real(iter)+1.0)/6.0;
                // if (relaxation > 1.0)
                //     relaxation = 1.0;

                // cout << "Iteration: " << iter << endl; 
                // vec<real,3> temp = pnp1.xi;
                
                real const dt = pnp1.dt;
                if(svar.eqOrder == 2)
                {
                    real dtm1 = iter > 0 ? pn.dt : pnp1.dt;
                    BFD2(fvar, avar, svar.grav, dt, dtm1, pnm1, pn, pnp1);
                }
                else
                    BFD1(fvar, avar, svar.grav, dt, pn, pnp1);

                real dt_temp = FindFace(svar, cells, pn, pnp1);
   
                // Implement relaxation, since it appears to struggle with convergence
                if(iter > svar.nrelax)
                    pnp1.dt = (1.0 - svar.relax) * dt_temp + svar.relax * pnp1.dt;
                else
                    pnp1.dt = dt_temp;

                /* Now c+1 is known, use the central difference face properties*/
                pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
                pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

                error = log10(fabs(pnp1.dt-dt));

                #ifdef DEBUG
                // cout << iter << "  "  << error << endl;
                timelog << setw(width) << iter << setw(width) << pnp1.dt << setw(width) 
                    << pnp1.acc << setw(width) << error << endl;
                #endif
                iter++;
            }

            if(svar.partout == 1)
                Write_Point(svar, pnp1);

            if ( pnp1.going == 0 || pnp1.v.norm() > 1e3
                    || (pnp1.xi - pn.xi).norm() > cells.maxlength )
            {
                #ifdef DEBUG
                if((pnp1.xi - pn.xi).norm() > cells.maxlength)
                {
                    cout << "Particle terminated because of max length" << endl;
                    cout << (pnp1.xi - pn.xi).norm() << "  " << cells.maxlength << endl;
                }
                timelog.close();
                #endif
                pnp1.failed = 1;
                Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data, iptdata);
                #pragma omp atomic
                    svar.nFailed++;
                break;
            }
            else if ( pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x)
            {    
                #ifdef DEBUG
                    cout << "Particle terminated because of a success" << endl;
                    timelog.close();
                #endif
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
            
            #ifdef DEBUG
            step++;
            #endif
        }        
    }
}