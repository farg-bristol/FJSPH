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
inline void add_file_aux_data(void* const& file, std::string const& name, std::string const& var)
{
    if (!var.empty())
    {
        if (tecDataSetAddAuxData(file, space2underscore(name).c_str(), space2underscore(var).c_str()))
        {
            std::cout << "Failed to write file auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }
}

inline void add_file_aux_data(void* const& file, std::string const& name, size_t const& var)
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

inline void add_file_aux_data(void* const& file, std::string const& name, int const& var)
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

inline void add_file_aux_data(void* const& file, std::string const& name, uint const& var)
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

inline void add_file_aux_data(void* const& file, std::string const& name, double const& var)
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

inline void add_file_aux_data(void* const& file, std::string const& name, StateVecD const& var)
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

inline void add_zone_aux_data(
    void* const& file, int32_t const& zone, std::string const& name, std::string const& var
)
{
    if (!var.empty())
    {
        if (tecZoneAddAuxData(file, zone, space2underscore(name).c_str(), space2underscore(var).c_str()))
        {
            std::cout << "Failed to write zone " << zone << " auxiliary data: " << name << std::endl;
            std::cout << "Stopping." << std::endl;
            exit(-1);
        }
    }
}

inline void
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

inline void
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

inline void
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

inline void
add_zone_aux_data(void* const& file, int32_t const& zone, std::string const& name, StateVecD const& var)
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

inline void Write_Real_Vector(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& start,
    int32_t const& size, vector<real> const& varVec, string const& varName
)
{
    int retval;
#if FOD == 1
    retval = tecZoneVarWriteDoubleValues(fileHandle, outputZone, varCount, 0, size, &varVec[start]);
#else
    retval = tecZoneVarWriteFloatValues(fileHandle, outputZone, varCount, 0, size, &varVec[start]);
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

inline void Write_State_Vector(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& start,
    int32_t const& size, vector<StateVecD> const& varVec, string const& varName
)
{
    vector<real> vec(size);
    for (uint dim = 0; dim < SIMDIM; ++dim)
    {
#pragma omp parallel for
        for (int ii = start; ii < start + size; ++ii)
            vec[ii - start] = varVec[ii][dim];

        string name = varName + " " + std::to_string(dim);
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, size, vec, name);
    }
}

inline void Write_Int_Vector(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& start,
    int32_t const& size, vector<int32_t> const& varVec, string const& varName
)
{
    if (tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &varVec[start]))
    {
        printf(
            "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
        );
        exit(-1);
    }
    varCount++;
}

inline void Write_UInt_Vector(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& start,
    int32_t const& size, vector<uint8_t> const& varVec, string const& varName
)
{
    if (tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &varVec[start]))
    {
        printf(
            "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
        );
        exit(-1);
    }
    varCount++;
}

inline void Write_Real_Value(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
    real const& value, string const& varName
)
{
    int retval;
    vector<real> rvec(size, 0.0);
    rvec[0] = value;
#if FOD == 1
    retval = tecZoneVarWriteDoubleValues(fileHandle, outputZone, varCount, 0, size, &rvec[0]);
#else
    retval = tecZoneVarWriteFloatValues(fileHandle, outputZone, varCount, 0, size, &rvec[0]);
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

inline void Write_Int_Value(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
    int32_t const& value, string const& varName
)
{
    int retval;
    vector<int32_t> rvec(size, 0.0);
    rvec[0] = value;

    retval = tecZoneVarWriteInt32Values(fileHandle, outputZone, varCount, 0, size, &rvec[0]);

    if (retval)
    {
        printf(
            "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
        );
        exit(-1);
    }
    varCount++;
}

inline void Write_UInt_Value(
    void* const& fileHandle, int32_t const& outputZone, int32_t& varCount, int32_t const& size,
    uint8_t const& value, string const& varName
)
{
    int retval;
    vector<uint8_t> rvec(size, 0.0);
    rvec[0] = value;

    retval = tecZoneVarWriteUInt8Values(fileHandle, outputZone, varCount, 0, size, &rvec[0]);

    if (retval)
    {
        printf(
            "Failed to write \"%s\". zone: %d. varCount: %d\n", varName.c_str(), outputZone, varCount
        );
        exit(-1);
    }
    varCount++;
}

/*************************************************************************/
/*************************** BINARY OUTPUTS ******************************/
/*************************************************************************/
inline void Write_Zone(
    SPHState const& pnp1, real const& rho_rest, real const& scale, OutputMap const& output_variables,
    size_t const& start, size_t const& end, void* const& fileHandle, int32_t const& outputZone
)
{
    int32_t varCount = 1;
    int32_t imax = end - start;

    vector<real> vec(imax);
    // Restart essential variables. These should not be zero, but could possibly in future
    if (output_variables.at("pos-vec").write)
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].xi(dim) / scale;

            string name = "Position coordinate " + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("vel-vec").write)
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].v(dim);

            string name = "Velocity " + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("acc-vec").write)
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].acc(dim);

            string name = "Acceleration " + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("press").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].p;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Pressure");
    }

    if (output_variables.at("dRho").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].Rrho;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Density gradient");
    }

    if (output_variables.at("part_id").write)
    {
        vector<int> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].part_id;
        Write_Int_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Particle ID");
    }

    if (output_variables.at("cellID").write)
    {
        vector<int> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].cellID;
        Write_Int_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Cell ID");
    }

    if (output_variables.at("bound").write)
    {
        vector<uint8_t> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].b;
        Write_UInt_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Particle status condition");
    }

    // Non essential variables, but to provide further information
    if (output_variables.at("dens").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].rho;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Density");
    }

    if (output_variables.at("densVar").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = 100.0 * (pnp1[ii].rho / rho_rest - 1.0);
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Density Variation");
    }

    if (output_variables.at("vmag").write)
    {
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].v.norm();
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Velocity magnitude");
    }

    if (output_variables.at("surf").write)
    {
        vector<uint8_t> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].surf;
        Write_UInt_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Surface flag");
    }

    if (output_variables.at("surfZ").write)
    {
        vector<uint8_t> uvec(imax);
        for (size_t ii = start; ii < end; ++ii)
            uvec[ii - start] = pnp1[ii].surfzone;
        Write_UInt_Vector(fileHandle, outputZone, varCount, 0, imax, uvec, "Surface zone flag");
    }

    if (output_variables.at("aero-mag").write)
    {
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].Af.norm();
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Aerodynamic force magnitude");
    }

    if (output_variables.at("aero-vec").write)
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].Af(dim);

            string name = "Aerodynamic force " + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("curv").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].curve;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Curvature");
    }

    if (output_variables.at("occl").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].woccl;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Occlusion factor");
    }

    if (output_variables.at("cellP").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellP;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Cell pressure");
    }

    if (output_variables.at("cellRho").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellP;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Cell density");
    }

    if (output_variables.at("cellV-mag").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].cellV.norm();
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Cell velocity magnitude");
    }

    if (output_variables.at("cellV-vec").write)
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].cellV(dim);

            string name = "Cell velocity " + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("dsphG-vec").write) // Delta SPH density gradient
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].gradRho(dim);

            string name = "dSPH density gradient" + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("lam").write) // Lambda eigenvalue of L matrix
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].lam;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Lambda eigenvalue");
    }

    if (output_variables.at("lam-nb").write) // Lambda eigenvalue of L matrix without boundary
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].lam_nb;
        Write_Real_Vector(
            fileHandle, outputZone, varCount, 0, imax, vec, "Lambda eigenvalue without boundary"
        );
    }

    if (output_variables.at("colour").write) // Surface colour function
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].colour;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Colour function");
    }

    if (output_variables.at("colour-G").write) // Gradient of the colour function for CSF
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].colourG;
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Colour function gradient");
    }

    if (output_variables.at("norm-vec").write) // Surface normal
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].norm(dim);

            string name = "Surface normal" + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }

    if (output_variables.at("shiftV-mag").write)
    {
#pragma omp parallel for
        for (size_t ii = start; ii < end; ++ii)
            vec[ii - start] = pnp1[ii].vPert.norm();
        Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, "Shifting velocity magnitude");
    }

    if (output_variables.at("shiftV-vec").write)
    {
        for (uint dim = 0; dim < SIMDIM; ++dim)
        {
#pragma omp parallel for
            for (size_t ii = start; ii < end; ++ii)
                vec[ii - start] = pnp1[ii].vPert(dim);

            string name = "Shifting velocity " + std::to_string(dim);
            Write_Real_Vector(fileHandle, outputZone, varCount, 0, imax, vec, name);
        }
    }
}

void Write_Binary_Timestep(
    SIM const& svar, real const& rho_rest, SPHState const& pnp1, bound_block const& limits,
    char const* group, int32_t const& strandID, void* const& fileHandle
)
{
    size_t start = limits.index.first;
    size_t end = limits.index.second;
    int32_t const size = end - start;

    double solTime = svar.integrator.current_time;
    int32_t outputZone;

    // Get zone data types
    vector<int32_t> varTypes = svar.io.var_types;
    vector<int32_t> shareVarFromZone(varTypes.size(), 0);
    vector<int32_t> valueLocation(varTypes.size(), 1);
    vector<int32_t> passiveVarList(varTypes.size(), 0);

    if (tecZoneCreateIJK(
            fileHandle, group, size, 1, 1, &varTypes[0], &shareVarFromZone[0], &valueLocation[0],
            &passiveVarList[0], 0, 0, 0, &outputZone
        ))
    {
        cerr << "Failed to create IJK zone." << endl;
        exit(-1);
    }

    if (strandID != 0)
        if (tecZoneSetUnsteadyOptions(fileHandle, outputZone, solTime, strandID))
        {
            cerr << "Failed to add unsteady options." << endl;
            exit(-1);
        }

    Write_Zone(pnp1, rho_rest, svar.scale, svar.io.output_variables, start, end, fileHandle, outputZone);

    // Write some auxiliary data, such as mass
    add_zone_aux_data(fileHandle, outputZone, "Particle mass", pnp1[start].m);
    // Add information about the insertion, aero, and deletion planes.
    add_zone_aux_data(fileHandle, outputZone, "Insertion normal x", limits.insert_norm[0]);
    add_zone_aux_data(fileHandle, outputZone, "Insertion normal y", limits.insert_norm[1]);
    add_zone_aux_data(fileHandle, outputZone, "Deletion normal x", limits.delete_norm[0]);
    add_zone_aux_data(fileHandle, outputZone, "Deletion normal y", limits.delete_norm[1]);
    add_zone_aux_data(fileHandle, outputZone, "Aerodynamic normal x", limits.aero_norm[0]);
    add_zone_aux_data(fileHandle, outputZone, "Aerodynamic normal y", limits.aero_norm[1]);
#if SIMDIM == 3
    add_zone_aux_data(fileHandle, outputZone, "Insertion normal z", limits.insert_norm[2]);
    add_zone_aux_data(fileHandle, outputZone, "Deletion normal z", limits.delete_norm[2]);
    add_zone_aux_data(fileHandle, outputZone, "Aerodynamic normal z", limits.aero_norm[2]);
#endif

    add_zone_aux_data(fileHandle, outputZone, "Insertion plane constant", limits.insconst);
    add_zone_aux_data(fileHandle, outputZone, "Deletion plane constant", limits.delconst);
    add_zone_aux_data(fileHandle, outputZone, "Aerodynamic plane constant", limits.aeroconst);

    // Add boundary time information
    if (limits.nTimes != 0)
    {
        add_zone_aux_data(fileHandle, outputZone, "Times count", limits.nTimes);
        for (size_t time = 0; time < limits.nTimes; time++)
        {
            add_zone_aux_data(
                fileHandle, outputZone, "Time " + std::to_string(time), limits.times[time]
            );
            add_zone_aux_data(
                fileHandle, outputZone, "Velocity " + std::to_string(time) + " x", limits.vels[time][0]
            );
            add_zone_aux_data(
                fileHandle, outputZone, "Velocity " + std::to_string(time) + " y", limits.vels[time][1]
            );
#if SIMDIM == 3
            add_zone_aux_data(
                fileHandle, outputZone, "Velocity " + std::to_string(time) + " z", limits.vels[time][2]
            );
#endif
        }
    }
    else
    {
        add_zone_aux_data(fileHandle, outputZone, "Times count", 0);
        add_zone_aux_data(fileHandle, outputZone, "Velocity x", limits.vels[0][0]);
        add_zone_aux_data(fileHandle, outputZone, "Velocity y", limits.vels[0][1]);
#if SIMDIM == 3
        add_zone_aux_data(fileHandle, outputZone, "Velocity z", limits.vels[0][2]);
#endif
    }

    add_zone_aux_data(
        fileHandle, outputZone, "Fixed velocity or dynamic inlet", limits.fixed_vel_or_dynamic
    );
    add_zone_aux_data(fileHandle, outputZone, "Lattice or HCPL packing", limits.particle_order);
    add_zone_aux_data(fileHandle, outputZone, "Boundary solver type", limits.bound_solver);
    add_zone_aux_data(fileHandle, outputZone, "Boundary is no slip", limits.no_slip);
    add_zone_aux_data(fileHandle, outputZone, "Block type", limits.block_type);

    // cout << "Flushing results." << endl;
    // INTEGER4 numZonesToRetain = 0;
    if (tecFileWriterFlush(fileHandle, 0, NULL))
    {
        cout << "Failed to flush data. Retrying..." << endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        // retry the flush
        if (tecFileWriterFlush(fileHandle, 0, NULL))
        {
            cerr << "Failed to flush data to file: " << fileHandle << endl;
            exit(-1);
        }
    }
}

void Init_Binary_PLT(
    SIM& svar, string const& prefix, string const& filename, string const& zoneName, void*& fileHandle
)
{
    int32_t FileType = 0;   /*0 = Full, 1 = Grid, 2 = Solution*/
    int32_t fileFormat = 1; // 0 == PLT, 1 == SZPLT

    string file = prefix + filename;

    if (tecFileWriterOpen(
            file.c_str(), zoneName.c_str(), svar.io.var_names.c_str(), fileFormat, FileType, 1, NULL,
            &fileHandle
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

    // Write crucial simulation auxiliary data. Potential to expand
    add_file_aux_data(fileHandle, "Number of boundary blocks", svar.n_bound_blocks);
    add_file_aux_data(fileHandle, "Number of fluid blocks", svar.n_fluid_blocks);
    add_file_aux_data(fileHandle, "Number of boundary points", svar.bound_points);
    add_file_aux_data(fileHandle, "Number of fluid points", svar.fluid_points);
    add_file_aux_data(fileHandle, "Current frame", svar.integrator.current_frame);

    add_file_aux_data(fileHandle, "Input para filename", svar.io.input_file);
    add_file_aux_data(fileHandle, "Input fluid definition filename", svar.io.input_fluid_file);
    add_file_aux_data(fileHandle, "Input boundary definition filename", svar.io.input_bound_file);

    // add_file_aux_data(fileHandle, "Primary grid face filename", svar.tau_mesh);
    // add_file_aux_data(fileHandle, "Boundary mapping filename", svar.tau_bmap);
    // add_file_aux_data(fileHandle, "Restart-data prefix", svar.tau_sol);
    add_file_aux_data(fileHandle, "Dimension offset vector", svar.io.offset_axis);
    add_file_aux_data(fileHandle, "Solution angle of attack", svar.io.angle_alpha);
    add_file_aux_data(fileHandle, "Grid scale", svar.scale);

    add_file_aux_data(fileHandle, "OpenFOAM input directory", svar.io.foam_dir);
    add_file_aux_data(fileHandle, "OpenFOAM solution directory", svar.io.foam_sol);
    add_file_aux_data(fileHandle, "OpenFOAM binary", svar.io.foam_is_binary);
    add_file_aux_data(fileHandle, "Label size", svar.io.foam_label_size);
    add_file_aux_data(fileHandle, "Scalar size", svar.io.foam_scalar_size);
    add_file_aux_data(fileHandle, "OpenFOAM buoyant", svar.io.foam_buoyant_sim);

    add_file_aux_data(fileHandle, "VLM definition filename", svar.io.vlm_file);

    add_file_aux_data(fileHandle, "Single file for output", svar.io.single_file);
    add_file_aux_data(fileHandle, "Output files prefix", svar.io.output_prefix);
    add_file_aux_data(fileHandle, "SPH restart prefix", svar.io.restart_prefix);
    add_file_aux_data(fileHandle, "SPH frame time interval", svar.integrator.frame_time_interval);
    add_file_aux_data(fileHandle, "SPH frame count", svar.integrator.max_frames);
    add_file_aux_data(fileHandle, "SPH output encoding", svar.io.out_encoding);
    add_file_aux_data(fileHandle, "Variable list", svar.io.output_names);

    /* Fluid data */
    add_file_aux_data(fileHandle, "Reference density", svar.air.rho_g);
    add_file_aux_data(fileHandle, "Reference dispersed density", svar.fluid.rho_rest);
    add_file_aux_data(fileHandle, "Sutherland reference viscosity", svar.air.mu_g);
    add_file_aux_data(fileHandle, "Reference dispersed viscosity", svar.fluid.mu);
    add_file_aux_data(fileHandle, "Reference surface tension", svar.fluid.sig);
    add_file_aux_data(fileHandle, "SPH surface tension contact angle", svar.fluid.contangb);
    add_file_aux_data(fileHandle, "Init hydrostatic pressure", svar.init_hydro_pressure);
    add_file_aux_data(fileHandle, "Hydrostatic height", svar.hydro_height);

    /* Aerodynamic data */
    add_file_aux_data(fileHandle, "Reference velocity", svar.air.v_ref);
    add_file_aux_data(fileHandle, "Reference pressure", svar.air.p_ref);
    add_file_aux_data(fileHandle, "Reference Mach number", svar.air.M_ref);
    add_file_aux_data(fileHandle, "Reference temperature", svar.air.temp_g);
    add_file_aux_data(fileHandle, "Gas constant gamma", svar.air.gamma);

    /* Simulation settings */
    add_file_aux_data(fileHandle, "SPH integration solver", svar.integrator.solver_name);
    add_file_aux_data(fileHandle, "SPH boundary solver", svar.integrator.bound_solver);
    add_file_aux_data(fileHandle, "SPH solver minimum residual", svar.integrator.min_residual);
    add_file_aux_data(fileHandle, "SPH maximum timestep", svar.integrator.delta_t_max);
    add_file_aux_data(fileHandle, "SPH minimum timestep", svar.integrator.delta_t_min);
    add_file_aux_data(fileHandle, "SPH maximum CFL", svar.integrator.cfl_max);
    add_file_aux_data(fileHandle, "SPH minimum CFL", svar.integrator.cfl_min);
    add_file_aux_data(fileHandle, "SPH CFL condition", svar.integrator.cfl);
    add_file_aux_data(fileHandle, "SPH unstable CFL step", svar.integrator.cfl_step);
    add_file_aux_data(fileHandle, "SPH unstable CFL count limit", svar.integrator.n_unstable_limit);
    add_file_aux_data(fileHandle, "SPH stable CFL count limit", svar.integrator.n_stable_limit);
    add_file_aux_data(
        fileHandle, "SPH stable CFL count iteration factor", svar.integrator.subits_factor
    );
    add_file_aux_data(fileHandle, "SPH maximum shifting velocity", svar.integrator.max_shift_vel);

    add_file_aux_data(fileHandle, "SPH background pressure", svar.fluid.press_back);
    add_file_aux_data(fileHandle, "SPH starting pressure", svar.fluid.press_pipe);
    add_file_aux_data(fileHandle, "SPH density variation", svar.fluid.rho_var);
    add_file_aux_data(fileHandle, "SPH maximum density", svar.fluid.rho_max);
    add_file_aux_data(fileHandle, "SPH minimum density", svar.fluid.rho_min);
    add_file_aux_data(fileHandle, "SPH delta coefficient", svar.fluid.dsph_delta);

    add_file_aux_data(fileHandle, "SPH artificial viscosity factor", svar.fluid.visc_alpha);
    add_file_aux_data(fileHandle, "SPH speed of sound", svar.fluid.speed_sound);
    add_file_aux_data(fileHandle, "SPH Newmark Beta iteration limit", svar.integrator.max_subits);
    add_file_aux_data(fileHandle, "SPH gravity vector", svar.grav);

    add_file_aux_data(fileHandle, "SPH initial spacing", svar.particle_step);
    add_file_aux_data(fileHandle, "SPH boundary spacing factor", svar.bound_step_factor);
    add_file_aux_data(fileHandle, "SPH smoothing length factor", svar.fluid.H_fac);
    add_file_aux_data(fileHandle, "SPH aerodynamic case", svar.air.aero_case);
    add_file_aux_data(fileHandle, "SPH SP diameter definition", svar.air.use_dx);
    add_file_aux_data(fileHandle, "SPH use TAB deformation", svar.air.use_TAB_def);
    add_file_aux_data(fileHandle, "SPH global offset coordinate", svar.offset_vec);
    add_file_aux_data(fileHandle, "SPH maximum particle count", svar.max_points);
    add_file_aux_data(fileHandle, "SPH aerodynamic cutoff value", svar.air.lam_cutoff);
    add_file_aux_data(fileHandle, "SPH freestream velocity", svar.air.v_inf);
    add_file_aux_data(fileHandle, "SPH restart fit tolerance", svar.io.restart_tol);

    /* Particle tracking settings */
    add_file_aux_data(fileHandle, "Transition to IPT", svar.ipt.using_ipt);
    add_file_aux_data(fileHandle, "Velocity equation order", svar.ipt.ipt_eq_order);
    add_file_aux_data(fileHandle, "SPH tracking conversion x coordinate", svar.ipt.max_x_sph);
    add_file_aux_data(fileHandle, "Maximum x trajectory coordinate", svar.ipt.max_x);
    add_file_aux_data(fileHandle, "Particle scatter output", svar.ipt.part_out);
    add_file_aux_data(fileHandle, "Particle streak output", svar.ipt.streak_out);
    add_file_aux_data(fileHandle, "Particle cell intersection output", svar.ipt.cells_out);
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

            vector<real> x(size);
            int32_t var = 1;

            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
                x[0] = pnp1.xi(dim) / svar.scale;

                string name = "position coordinate ";
                name.append(std::to_string(dim));
                Write_Real_Vector(svar.ipt.part_handle, outputZone, var, 0, size, x, name);
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
            pID[0] = pnp1.part_id;
            cID[0] = pnp1.cellID;

            Write_Real_Vector(svar.ipt.part_handle, outputZone, var, 0, size, t, "particle time");
            Write_Real_Vector(svar.ipt.part_handle, outputZone, var, 0, size, dt, "timestep");
            Write_Real_Vector(svar.ipt.part_handle, outputZone, var, 0, size, vel, "velocity magnitude");
            Write_Real_Vector(
                svar.ipt.part_handle, outputZone, var, 0, size, acc, "acceleration magnitude"
            );
            Write_Int_Vector(svar.ipt.part_handle, outputZone, var, 0, size, pID, "particle ID");
            Write_Real_Vector(svar.ipt.part_handle, outputZone, var, 0, size, cVel, "cell velocity");
            Write_Real_Vector(svar.ipt.part_handle, outputZone, var, 0, size, cRho, "cell density");
            Write_Int_Vector(svar.ipt.part_handle, outputZone, var, 0, size, cID, "cell ID");

            if (tecFileWriterFlush(svar.ipt.part_handle, 0, NULL))
            {
                cout << "Failed to flush data. Retrying..." << endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                // retry the flush
                if (tecFileWriterFlush(svar.ipt.part_handle, 0, NULL))
                {
                    cerr << "Failed to flush data to particle file" << endl;
                    exit(-1);
                }
            }
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

            vector<real> x(size);
            int32_t var = 1;

            for (uint dim = 0; dim < SIMDIM; ++dim)
            {
#pragma omp parallel for
                for (uint ii = 0; ii < size; ++ii)
                    x[ii] = pnp1[ii].xi(dim) / svar.scale;

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
                pID[ii] = pnp1[ii].part_id;
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
                Write_Real_Vector(svar.ipt.cell_handle, outputZone, var, 0, nNodes, x, name);
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