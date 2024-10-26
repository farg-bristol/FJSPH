/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include <chrono>

#include "BinaryIO.h"
#include "CDFIO.h"
#include "Containment.h"
#include "FOAMIO.h"
#include "Geometry.h"
#include "H5IO.h"
#include "IO.h"
#include "IOFunctions.h"
#include "IPT.h"
#include "Init.h"
#include "Integration.h"
#include "Neighbours.h"
#include "Shifting.h"
#include "Var.h"

#ifdef DEBUG
/*Open debug file to write to*/
FILE* dbout;
#endif

using namespace std::chrono;
using namespace nanoflann;

int main(int argc, char* argv[])
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    high_resolution_clock::time_point t2;
    Eigen::initParallel();
    // omp_set_num_threads(1);
    srand(unsigned(time(NULL)));

#ifdef DEBUG
    dbout = fopen("WCSPH.log", "w");
#endif

    real duration;
    real error = 0;
    cout.width(13);
    cout << std::scientific << std::left << std::setprecision(4);

    write_header();

    /******* Define the global simulation parameters ******/
    SIM svar;
    FLUID svar.fluid;
    AERO svar.air;
    LIMITS limits;
    OUTL outlist;
    MESH cells;
    VLM vortex;
    SURFS surf_marks;

    ///****** Initialise the particles memory *********/
    SPHState pn;              /*Particles at n   */
    SPHState pnp1;            /*Particles at n+1 */
    vector<IPTState> iptdata; /* Particle tracking data */
    // SPHState airP;
    // DELTAP dp;

    GetInput(argc, argv, svar, vortex);

    omp_set_num_threads(svar.numThreads);

    // if(MakeOutputDir(argc,argv,svar))
    // {
    // 	cout << "Couldn't make output directory. Please check permissions." << endl;
    // 	exit(-1);
    // }

    if (svar.Asource == meshInfl)
    {
#if SIMDIM == 3
        if (svar.mesh_source == TAU_CDF)
        {
            vector<uint> empty;
            TAU::Read_BMAP(svar);
            TAU::Read_tau_mesh_FACE(svar, cells, svar.fluid, svar.air);
            TAU::Read_SOLUTION(svar, svar.offset_axis, cells, empty);
        }
        else
        {
            FOAM::Read_FOAM(svar, cells);
        }
#else
        if (svar.mesh_source == OpenFOAM)
        {
            vector<uint> used_verts;
            TAU::Read_BMAP(svar);
            TAU::Read_tau_mesh_EDGE(svar, cells, svar.fluid, svar.air, used_verts);
            TAU::Read_SOLUTION(svar, svar.offset_axis, cells, used_verts);
        }
        else
        {
            cout << "OpenFOAM mesh input is currently unsupported in two-dimensions" << endl;
            exit(-1);
        }
#endif

        Init_Surface(svar, cells, surf_marks);
    }

    cout << std::setprecision(5);

    pn.reserve(svar.max_points);
    pnp1.reserve(svar.max_points);

    if (!svar.restart)
    {
        svar.current_time = 0.0; /*Total simulation time*/
        Init_Particles(svar, pn, pnp1, limits);

        // Redefine the mass and spacing to make sure the required mass is conserved
        // if(svar.Scase == SPHERE || svar.Scase == JET)
        // 	Set_Mass(svar,svar.fluid,svar.air,pn,pnp1);
        if (svar.Asource != meshInfl)
        {
            for (size_t ii = 0; ii < pnp1.size(); ++ii)
            {
                pn[ii].cellRho = svar.air.rho_g;
                pn[ii].cellP = svar.air.p_ref;
                pnp1[ii].cellRho = svar.air.rho_g;
                pnp1[ii].cellP = svar.air.p_ref;
                pn[ii].cellV = svar.air.v_inf;
                pnp1[ii].cellV = svar.air.v_inf;
            }
        }

        Remove_Old_Files(svar);
    }
    else
    {
        Init_Particles_Restart(svar, limits);
        /* Read the files */
        Read_HDF5(svar, vortex, pn, pnp1, limits);
    }

    // Check if cells have been initialsed before making a tree off it
    if (cells.cCentre.size() == 0)
        cells.cCentre.emplace_back(StateVecD::Zero());

    cout << "Starting counts: " << endl;
    cout << "Boundary: " << svar.bound_points << "  Sim: " << svar.fluid_points << endl << endl;

    ///********* Tree algorithm stuff ************/
    Sim_Tree SPH_TREE(SIMDIM, pnp1, 20);
    outlist = update_neighbours(SPH_TREE, svar.fluid, pnp1);

    Vec_Tree CELL_TREE(SIMDIM, cells.cCentre, 10);
    CELL_TREE.index->buildIndex();

    ///********* Declare the time integrator ************/
    Integrator integrator(svar.solver_type);

    ///*************** Open simulation files ***************/
    FILE* ff = NULL;
    FILE* f2 = NULL;
    FILE* f3 = NULL;
    FILE* fb = NULL;
    FILE* fg = NULL;
    string framef = svar.output_prefix;
    framef.append("_frame.info");
    if (svar.restart == 1)
        f2 = fopen(framef.c_str(), "a");
    else
        f2 = fopen(framef.c_str(), "w");

    if (svar.single_file)
    {
        if (svar.write_tecio)
            Write_Tec_Headers(ff, fb, fg, svar, svar.output_prefix);
        if (svar.write_h5part)
            Write_h5part_Headers(svar, svar.output_prefix);
    }

    if (!svar.restart)
    {
        Write_Timestep(ff, fb, fg, svar, limits, pnp1);
    }

    ///*** Perform an iteration to populate the vectors *****/
    if (!svar.restart)
    {
        integrator.integrate_no_update(
            SPH_TREE, CELL_TREE, svar, vortex, cells, limits, outlist, pn, pnp1
        );
    }
    else
    {
        size_t const start = svar.bound_points;
        size_t end = svar.total_points;
        size_t end_ng = svar.bound_points + svar.fluid_points;
        real npd = 1.0;
        dSPH_PreStep(svar.fluid, end, pnp1, outlist, npd);

        if (svar.Asource == meshInfl)
        {
            // cout << "Finding cells" << endl;
            FindCell(svar, CELL_TREE, cells, pn, pnp1);
        }

        Detect_Surface(svar, start, end_ng, outlist, cells, vortex, pn);

// Apply_XSPH(svar.fluid,start,end,outlist,dp,pnp1);
#ifdef ALE
        particle_shift(svar, start, end_ng, outlist, pn);
#endif

        /* Update shifting velocity and surface data */
        // SPHState::const_iterator first = pnp1.begin();
        // SPHState::const_iterator last = pnp1.begin() + end_ng;
        // pn = SPHState(first,last);
    }

    if (svar.using_ipt)
    {
        IPT::Init_IPT_Files(svar);
        IPT::Write_Data(svar, cells, iptdata);
    }

// Make sure that inlets are injecting inside the mesh domain.
#if SIMDIM == 3
    if (!svar.restart && svar.Asource == 3)
        vortex.write_VLM_Panels(svar.output_prefix);
#endif

    /*Timing calculation + error sum output*/
    t2 = high_resolution_clock::now();
    duration = duration_cast<seconds>(t2 - t1).count();

    if (!svar.restart)
    {
        svar.current_frame = 0;
        fprintf(f2, "Frame: %u\n", svar.current_frame);
        fprintf(
            f2, "Total Points: %zu Boundary Points: %zu Fluid Points: %zu\n", svar.total_points,
            svar.bound_points, svar.fluid_points
        );
        fprintf(
            f2, "Sim Time:  %.7g Comp Time: %.6e Error: %.6f Sub-iterations: %d\n", svar.current_time,
            duration, 0.0, 0
        );
        fprintf(
            f2, "Deleted particles: %zu Internal collisions: %zu\n", svar.delete_count,
            svar.internal_count
        );
    }
    else
    {
        cout << "Restarting simulation..." << endl;
        cout << "Frame: " << svar.current_frame << "  Sim Time: " << svar.current_time
             << "  Compute Time: " << duration << "  Error: " << error << endl;
    }

    ///************************* MAIN LOOP ********************/

#ifdef DEBUG
    const real a = 1 - svar.nb_gamma;
    const real b = svar.nb_gamma;
    const real c = 0.5 * (1 - 2 * svar.nb_beta);
    const real d = svar.nb_beta;
    const real B = svar.fluid.B;
    const real gam = svar.fluid.gam;
    fprintf(dbout, "Newmark Beta integration parameters\n");
    fprintf(dbout, "a: %f b: %f\n", a, b);
    fprintf(dbout, "c: %f d: %f\n", c, d);
    fprintf(dbout, "B: %f gam: %f\n\n", B, gam);
#endif

    for (svar.current_frame = 0; svar.current_frame < svar.max_frames; ++svar.current_frame)
    {
        int stepits = 0;
        real stept = 0.0;
        while (stept + 0.1 * svar.delta_t_min < svar.frame_time_interval)
        {
            if (stepits % 50 == 0)
#ifdef ALE
                printf("\nTime      | Timestep | CFL  | RMS error | its | dRho (%%) | Max-F     | "
                       "Max-Af    | Max Shift | Step time (ms)| \n");
#else
                printf("\nTime      | Timestep | CFL  | RMS error | its | dRho (%%) | Max-F     | "
                       "Max-Af    | Step time (ms)| \n");
#endif

            error = integrator.integrate(
                SPH_TREE, CELL_TREE, svar, vortex, cells, surf_marks, limits, outlist, pn, pnp1, iptdata
            );
            stept += svar.delta_t;
            ++stepits;
        }

        t2 = high_resolution_clock::now();
        duration = duration_cast<microseconds>(t2 - t1).count() / 1e6;

        /*Write each frame info to file*/
        fprintf(f2, "\nFrame: %u\n", svar.current_frame);
        fprintf(
            f2, "Total Points: %zu Boundary Points: %zu Fluid Points: %zu\n", svar.total_points,
            svar.bound_points, svar.fluid_points
        );
        fprintf(
            f2, "Sim Time:  %.7g Comp Time: %.6e Error: %.6f Sub-iterations: %d\n", svar.current_time,
            duration, error, stepits
        );
        fprintf(
            f2, "Deleted particles: %zu Internal collisions: %zu\n", svar.delete_count,
            svar.internal_count
        );

        cout << "Frame: " << svar.current_frame << "  Sim Time: " << svar.current_time
             << "  Compute Time: " << duration << "  Error: " << error << endl;
        cout << "Boundary particles:  " << svar.bound_points
             << " Sim particles: " << svar.total_points - svar.bound_points
             << " Deleted particles: " << svar.delete_count
             << " Internal collisions: " << svar.internal_count << endl;

        if (svar.total_points - svar.bound_points == 0)
        {
            cout << "No more points in the simulation space. Ending...." << endl;
            break;
        }

        Write_Timestep(ff, fb, fg, svar, limits, pnp1);
        Write_HDF5(svar, pnp1, limits);
        if (svar.using_ipt)
            IPT::Write_Data(svar, cells, iptdata);

        svar.last_frame_time += svar.frame_time_interval; /* March frame time forward */
    }

    /*Wrap up simulation files and close them*/
    if (f2 != NULL)
        fclose(f2);
    if (f3 != NULL)
        fclose(f3);
    if (ff != NULL)
        fclose(ff);
    if (fb != NULL)
        fclose(fb);
    if (fg != NULL)
        fclose(fg);

    if (svar.out_encoding == 1)
    {
        /*Combine the szplt files*/

        string outfile = svar.output_prefix + "_fluid.szplt";
        Combine_SZPLT(outfile);

        outfile = svar.output_prefix + "_boundary.szplt";
        Combine_SZPLT(outfile);
    }

    cout << "Simulation complete!" << endl;
    cout << "Time taken:\t" << duration << " seconds" << endl;
    cout << "Total simulation time:\t" << svar.current_time << " seconds" << endl;

    return 0;
}
