/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "BinaryIO.h"
#include <TECIO.h>
#include <chrono>
#include <filesystem> //Requires C++17 standard, needed for restart functionality
#include <thread>

using std::string;
using std::filesystem::directory_iterator;

inline std::string space2underscore(std::string text)
{
    std::replace(text.begin(), text.end(), ' ', '_');
    return text;
}

/*************************************************************************/
/********************** BINARY OUTPUT FUNCTIONS **************************/
/*************************************************************************/
namespace SZLIO
{
    bool
    open_file(string const& file_name, string const& zone_name, string var_names, void*& file_handle)
    {
        int32_t FileType = 0;   /*0 = Full, 1 = Grid, 2 = Solution*/
        int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT

        int32_t result = tecFileWriterOpen(
            file_name.c_str(), zone_name.c_str(), var_names.c_str(), fileFormat, FileType, 1, NULL,
            &file_handle
        );

#ifdef DEBUG
        if (tecFileSetDiagnosticsLevel(fileHandle, 1))
        {
            cerr << "Failed to set debug option for output file: " << file << endl;
            exit(-1);
        }
#endif

        return result == 0;
    }

    int32_t create_zone(
        void* const& file_handle, int32_t const& strandID, std::string const& zone_header,
        vector<int32_t> const& var_types, int64_t const& size, double const& solution_time
    )
    {
        int32_t output_zone = 0;
        // Get zone data types
        vector<int32_t> shareVarFromZone(var_types.size(), 0);
        vector<int32_t> valueLocation(var_types.size(), 1);
        vector<int32_t> passiveVarList(var_types.size(), 0);

        if (tecZoneCreateIJK(
                file_handle, zone_header.c_str(), size, 1, 1, &var_types[0], &shareVarFromZone[0],
                &valueLocation[0], &passiveVarList[0], 0, 0, 0, &output_zone
            ))
        {
            cerr << "Failed to create IJK zone." << endl;
            exit(-1);
        }

        if (strandID != 0)
            if (tecZoneSetUnsteadyOptions(file_handle, output_zone, solution_time, strandID))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

        return output_zone;
    }

    void close_file(void* handle)
    {
        if (tecFileWriterClose(&handle))
            printf("Failed to close szplt file.\n");
    }

    void flush_file(void* handle)
    {
        if (tecFileWriterFlush(handle, 0, NULL))
        {
            cout << "Failed to flush data. Retrying..." << endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            // retry the flush
            if (tecFileWriterFlush(handle, 0, NULL))
            {
                cerr << "Failed to flush data to szplt file, stopping." << endl;
                exit(-1);
            }
        }
    }

    void add_file_aux_data(void* const& file, std::string const& name, std::string const& var)
    {
        if (!var.empty())
        {
            if (tecDataSetAddAuxData(
                    file, space2underscore(name).c_str(), space2underscore(var).c_str()
                ))
            {
                std::cout << "Failed to write file auxiliary data: " << name << std::endl;
                std::cout << "Stopping." << std::endl;
                exit(-1);
            }
        }
    }

    void add_file_aux_data(void* const& file, std::string const& name, size_t const& var)
    {
        std::ostringstream stream;
        stream << var;
        if (tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void add_file_aux_data(void* const& file, std::string const& name, int const& var)
    {
        std::ostringstream stream;
        stream << var;
        if (tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void add_file_aux_data(void* const& file, std::string const& name, uint const& var)
    {
        std::ostringstream stream;
        stream << var;
        if (tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void add_file_aux_data(void* const& file, std::string const& name, double const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var;
        if (tecDataSetAddAuxData(file, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void add_file_aux_data(void* const& file, std::string const& name, StateVecD const& var)
    {
        string name_ = name + " x";
        std::ostringstream stream;
        stream << std::setprecision(9) << var[0];
        if (tecDataSetAddAuxData(file, space2underscore(name_).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }

        name_ = name + " y";
        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[1];
        if (tecDataSetAddAuxData(file, space2underscore(name_).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }

#if SIMDIM == 3
        name_ = name + " z";
        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[2];
        if (tecDataSetAddAuxData(file, space2underscore(name_).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
#endif
    }

    void add_zone_aux_data(
        void* const& file, int32_t const& zone, std::string const& name, std::string const& var
    )
    {
        if (!var.empty())
        {
            if (tecZoneAddAuxData(
                    file, zone, space2underscore(name).c_str(), space2underscore(var).c_str()
                ))
            {
                std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
                std::cout << "Stopping." << std::endl;
                exit(-1);
            }
        }
    }

    void
    add_zone_aux_data(void* const& file, int32_t const& zone, std::string const& name, size_t const& var)
    {
        std::ostringstream stream;
        stream << var;
        if (tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void
    add_zone_aux_data(void* const& file, int32_t const& zone, std::string const& name, int const& var)
    {
        std::ostringstream stream;
        stream << var;
        if (tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void
    add_zone_aux_data(void* const& file, int32_t const& zone, std::string const& name, double const& var)
    {
        std::ostringstream stream;
        stream << std::setprecision(9) << var;
        if (tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }

    void add_zone_aux_data(
        void* const& file, int32_t const& zone, std::string const& name, StateVecD const& var
    )
    {
        string name_ = name + " x";
        std::ostringstream stream;
        stream << std::setprecision(9) << var[0];
        if (tecZoneAddAuxData(file, zone, space2underscore(name_).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }

        name_ = name + " y";
        stream.str("");
        stream.clear();
        stream << var[1];
        if (tecZoneAddAuxData(file, zone, space2underscore(name_).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }

#if SIMDIM == 3
        name_ = name + " z";
        stream.str("");
        stream.clear();
        stream << std::setprecision(9) << var[2];
        if (tecZoneAddAuxData(file, zone, space2underscore(name_).c_str(), stream.str().c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
#endif
    }

    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<real> const& varVec, string const& varName
    )
    {
        int retval;
        int32_t size = varVec.size();
#if FOD == 1
        retval = tecZoneVarWriteDoubleValues(fileHandle, outputZone, varCount, 0, size, &varVec[0]);
#else
        retval = tecZoneVarWriteFloatValues(fileHandle, outputZone, varCount, 0, size, &varVec[0]);
#endif

        if (retval)
        {
            printf(
                "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
            );
            exit(-1);
        }
        varCount++;
    }

    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<int32_t> const& varVec, string const& varName
    )
    {
        int32_t size = varVec.size();
        if (tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &varVec[0]))
        {
            printf(
                "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
            );
            exit(-1);
        }
        varCount++;
    }

    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<uint8_t> const& varVec, string const& varName
    )
    {
        int32_t size = varVec.size();
        if (tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &varVec[0]))
        {
            printf(
                "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
            );
            exit(-1);
        }
        varCount++;
    }

    // X/Y/Z vector
    void Write_Vector(
        void* const& fileHandle, int32_t const& outputZone, int32_t& varCount,
        vector<StateVecD> const& varVec, string const& varName
    )
    {
        int32_t size = varVec.size();
        vector<real> vec(size);
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (int ii = 0; ii < size; ++ii)
                vec[ii] = varVec[ii][dim];

            string name = varName + " " + std::to_string(dim);
            Write_Vector(fileHandle, outputZone, varCount, vec, name);
        }
    }
} // namespace SZLIO

/*************************************************************************/
/************************* IPT BINARY OUTPUTS ****************************/
/*************************************************************************/
namespace IPT
{
    namespace BINARY
    {
        string Get_Variables(int const& offset)
        {
            string variables;

            if (offset == 1)
                variables = "Y,Z";
            else if (offset == 2)
                variables = "X,Z";
            else if (offset == 3)
                variables = "X,Y";
            else
                variables = "X,Y,Z";

            variables += ",t,dt,v,a,ptID,Cell_V,Cell_Rho,Cell_ID";
            return variables;
        }

        string Cell_Variables(int const& offset)
        {
            string variables;

            if (offset == 1)
                variables = "Y,Z";
            else if (offset == 2)
                variables = "X,Z";
            else if (offset == 3)
                variables = "X,Y";
            else
                variables = "X,Y,Z";

            return variables;
        }

        void Open_File(
            string const& file, string const& title, string const& var, int32_t const& fileType,
            void*& fileHandle
        )
        {
            if (tecFileWriterOpen(
                    file.c_str(), title.c_str(), var.c_str(), fileType, 0, 1, NULL, &fileHandle
                ))
            {
                cout << "Failed to open " << file << endl;
                exit(-1);
            }
#ifdef DEBUG
            if (tecFileSetDiagnosticsLevel(fileHandle, 1))
            {
                cerr << "Failed to set debug option for output file: " << file << endl;
                exit(-1);
            }
#endif
        }

        void Init_IPT_Files(SIM& svar)
        {
            string partf, cellf, streakf, surfacef;

            partf = svar.io.output_prefix + "_IPT_scatter.szplt";

            streakf = svar.io.output_prefix + "_streaks.szplt";

            cellf = svar.io.output_prefix + "_cells.plt";

            surfacef = svar.io.output_prefix + "_surface_impacts.szplt";

            // int32_t fileType   = 0; /*0 = Full, 1 = Grid, 2 = Solution*/
            // int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT
            string part_variables = Get_Variables(svar.io.offset_axis);
            string cell_variables = Cell_Variables(svar.io.offset_axis);

            if (svar.ipt.part_out == 1)
            {
                string title = "IPT particle scatter data";
                Open_File(partf, title, part_variables, 1, svar.ipt.part_handle);
            }

            if (svar.ipt.streak_out == 1)
            {
                string title = "IPT particle streak data";
                Open_File(streakf, title, part_variables, 1, svar.ipt.streak_handle);
            }

            if (svar.ipt.cells_out == 1)
            {
                string title = "IPT particle cell intersection data";
                Open_File(cellf, title, cell_variables, 0, svar.ipt.cell_handle);
            }
        }

        void Write_Point(SIM const& svar, IPTPart const& pnp1)
        {
            int64_t const size = 1;

            double solTime = pnp1.t;
            int32_t outputZone;

#if SIMDIM == 3
            vector<int32_t> varTypes = {realType, realType, realType, realType,
                                        3,        realType, realType, 3};
#else
            vector<int32_t> varTypes = {realType, realType, realType, 3, realType, realType, 3};
#endif
            vector<int32_t> shareVarFromZone(varTypes.size(), 0);
            vector<int32_t> valueLocation(varTypes.size(), 1);
            vector<int32_t> passiveVarList(varTypes.size(), 0);

            string group = "IPT Particle " + std::to_string(pnp1.part_id) + " scatter data";
            if (tecZoneCreateIJK(
                    svar.ipt.part_handle, group.c_str(), size, 1, 1, &varTypes[0], &shareVarFromZone[0],
                    &valueLocation[0], &passiveVarList[0], 0, 0, 0, &outputZone
                ))
            {
                cerr << "Failed to create IJK zone." << endl;
                exit(-1);
            }

            if (tecZoneSetUnsteadyOptions(svar.ipt.part_handle, outputZone, solTime, 4))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

            int32_t var = 1;

            vector<StateVecD> xi(size);
            vector<real> t(size);
            vector<real> dt(size);
            vector<real> vel(size);
            vector<real> acc(size);
            vector<real> cVel(size);
            vector<real> cRho(size);
            vector<int> pID(size);
            vector<int> cID(size);

            xi[0] = pnp1.xi / svar.scale;
            t[0] = pnp1.t;
            dt[0] = pnp1.dt;
            vel[0] = pnp1.v.norm();
            acc[0] = pnp1.acc;
            cVel[0] = pnp1.cellV.norm();
            cRho[0] = pnp1.cellRho;
            pID[0] = pnp1.part_id;
            cID[0] = pnp1.cellID;

            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, xi, "position coordinate");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, t, "particle time");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, dt, "timestep");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, vel, "velocity magnitude");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, acc, "acceleration magnitude");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, pID, "particle ID");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, cVel, "cell velocity");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, cRho, "cell density");
            SZLIO::Write_Vector(svar.ipt.part_handle, outputZone, var, cID, "cell ID");

            SZLIO::flush_file(svar.ipt.part_handle);
        }

        void
        Write_State(SIM const& svar, IPTState const& pnp1, string const& zoneName, void*& fileHandle)
        {
            int64_t const size = pnp1.size();

            double solTime = svar.integrator.current_time;
            int32_t outputZone;

#if SIMDIM == 3
            // variables =  "X,Y,Z,t,dt,v,a,ptID,Cell_V,Cell_Rho,Cell_ID"
            vector<int32_t> varTypes = {realType, realType, realType, realType, realType, realType,
                                        realType, 3,        realType, realType, 3};
#else
            // variables =  "X,Z,t,dt,v,a,ptID,Cell_V,Cell_Rho,Cell_ID"
            vector<int32_t> varTypes = {realType, realType, realType, realType, realType,
                                        realType, 3,        realType, realType, 3};
#endif
            vector<int32_t> shareVarFromZone(varTypes.size(), 0);
            vector<int32_t> valueLocation(varTypes.size(), 1);
            vector<int32_t> passiveVarList(varTypes.size(), 0);

            if (tecZoneCreateIJK(
                    fileHandle, zoneName.c_str(), size, 1, 1, &varTypes[0], &shareVarFromZone[0],
                    &valueLocation[0], &passiveVarList[0], 0, 0, 0, &outputZone
                ))
            {
                cerr << "Failed to create IJK zone." << endl;
                exit(-1);
            }

            if (tecZoneSetUnsteadyOptions(fileHandle, outputZone, solTime, 6))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

            int32_t var = 1;

            vector<StateVecD> xi(size);
            vector<real> t(size);
            vector<real> dt(size);
            vector<real> vel(size);
            vector<real> acc(size);
            vector<real> cVel(size);
            vector<real> cRho(size);
            vector<int> pID(size);
            vector<int> cID(size);

#pragma omp parallel for
            for (int ii = 0; ii < size; ++ii)
            {
                xi[ii] = pnp1[ii].xi / svar.scale;
                t[ii] = pnp1[ii].t;
                dt[ii] = pnp1[ii].dt;
                vel[ii] = pnp1[ii].v.norm();
                acc[ii] = pnp1[ii].acc;
                cVel[ii] = pnp1[ii].cellV.norm();
                cRho[ii] = pnp1[ii].cellRho;
                pID[ii] = pnp1[ii].part_id;
                cID[ii] = pnp1[ii].cellID;
            }

            SZLIO::Write_Vector(fileHandle, outputZone, var, xi, "position coordinate");
            SZLIO::Write_Vector(fileHandle, outputZone, var, t, "particle time");
            SZLIO::Write_Vector(fileHandle, outputZone, var, dt, "timestep");
            SZLIO::Write_Vector(fileHandle, outputZone, var, vel, "velocity magnitude");
            SZLIO::Write_Vector(fileHandle, outputZone, var, acc, "acceleration magnitude");
            SZLIO::Write_Vector(fileHandle, outputZone, var, pID, "particle ID");
            SZLIO::Write_Vector(fileHandle, outputZone, var, cVel, "cell velocity");
            SZLIO::Write_Vector(fileHandle, outputZone, var, cRho, "cell density");
            SZLIO::Write_Vector(fileHandle, outputZone, var, cID, "cell ID");

            SZLIO::flush_file(fileHandle);
        }

        void Write_Cells(
            SIM& svar, MESH const& cells, IPTState const& pnp1, int32_t const& totalNumFaceNodes,
            vector<StateVecD> const& usedVerts, vector<int32_t> const& cellIndexes,
            vector<vector<int32_t>> const& faces, vector<int32_t> const& left,
            vector<int32_t> const& right
        )
        {
            // int64_t const size = pnp1.size();

            double solTime = svar.integrator.current_time;
            int32_t outputZone;

#if SIMDIM == 3
            vector<int32_t> varTypes = {realType, realType, realType, realType,
                                        3,        realType, realType, 3};
#else
            vector<int32_t> varTypes = {realType, realType, realType, 3, realType, realType, 3};
#endif
            vector<int32_t> shareVarFromZone(varTypes.size(), 0);
            vector<int32_t> valueLocation(varTypes.size(), 1);
            vector<int32_t> passiveVarList(varTypes.size(), 0);

#if SIMDIM == 3
            int32_t zoneType = 7; /* FE Polyhedron */
#else
            int32_t zoneType = 6; /* FE Polygon */
#endif
            int32_t nNodes = usedVerts.size();
            int32_t nFaces = faces.size();
            int32_t nCells = cellIndexes.size();

            string group = "IPT cell intersection data";
            if (tecZoneCreatePoly(
                    svar.ipt.cell_handle, group.c_str(), zoneType, nNodes, nFaces, nCells,
                    totalNumFaceNodes, &varTypes[0], &shareVarFromZone[0], &valueLocation[0],
                    &passiveVarList[0], 0, 0, 0, &outputZone
                ))
            {
                cerr << "Failed to create polyhedral/polygonal zone." << endl;
                exit(-1);
            }

            if (tecZoneSetUnsteadyOptions(svar.ipt.cell_handle, outputZone, solTime, 7))
            {
                cerr << "Failed to add unsteady options." << endl;
                exit(-1);
            }

            vector<real> x(nNodes);
            int32_t var = 1;

            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (int32_t ii = 0; ii < nNodes; ++ii)
                    x[ii] = usedVerts[ii][dim];

                string name = "position coordinate " + std::to_string(dim);
                SZLIO::Write_Vector(svar.ipt.cell_handle, outputZone, var, x, name);
            }

            vector<int32_t> faceCounts(faces.size());
            vector<int32_t> faceNodes(totalNumFaceNodes);
            size_t jj = 0;
            /*Inform of how many vertices in each face*/
            for (size_t ii = 0; ii < faces.size(); ++ii)
            {
                faceCounts[ii] = faces[ii].size();
                for (auto const& vertex : faces[ii])
                { /*Write face vertex indexes*/
                    faceNodes[jj] = vertex;
                    jj++;
                }
            }

            tecZoneWritePolyFaces32(
                svar.ipt.cell_handle, outputZone, 0, nFaces, &faceCounts[0], &faceNodes[0], &left[0],
                &right[0], 1
            );
        }
    } // namespace BINARY
} // namespace IPT