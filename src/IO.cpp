/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "IO.h"
#include "AsciiIO.h"
#include "BinaryIO.h"
#include "Geometry.h"
#include "H5IO.h"
#include "IOFunctions.h"
#include "Kernel.h"
#include "VLM.h"

#include "Third_Party/Eigen/LU"

#include <ctime>
#include <filesystem>
#include <regex>

using std::filesystem::directory_iterator;

/*************************************************************************/
/**************************** ASCII INPUTS *******************************/
/*************************************************************************/

/* Set simulation parameters based off the user file inputs */
void Set_Values(SIM& svar, FLUID& fvar, AERO& avar, VLM& vortex)
{
    svar.offset_vec *= svar.scale;
    svar.max_x *= svar.scale;
    svar.max_x_sph *= svar.scale;

    fvar.B = fvar.rho_rest * pow(fvar.speed_sound, 2) / fvar.gam; /*Factor for Cole's Eq*/

    /*Pipe Pressure calc*/
    fvar.rho_pipe = fvar.get_density(fvar.press_pipe);

    /* Upper and lower limits for density */
    if (fvar.rho_max == 1500 && fvar.rho_min == 500)
    { /* If limits are undefined, use a variation around the base density */
        fvar.rho_max = fvar.rho_rest * (1.0 + fvar.rho_var * 0.01);
        fvar.rho_min = fvar.rho_rest * (1.0 - fvar.rho_var * 0.01);
    }

    svar.dx = svar.particle_step * pow(fvar.rho_pipe / fvar.rho_rest, 1.0 / SIMDIM);

    svar.nb_beta = 0.25;
    svar.nb_gamma = 0.5; /*Newmark Beta parameters*/

    /*Mass from spacing and density*/
    fvar.sim_mass = fvar.rho_rest * pow(svar.particle_step, SIMDIM);
    fvar.bnd_mass = fvar.sim_mass;
    avar.mass_g = avar.rho_g * pow(svar.particle_step, SIMDIM);

    avar.sos = sqrt(avar.temp_g * avar.R_g * avar.gamma);
    avar.i_sos_sq = 1.0 / (avar.sos * avar.sos);

    if (avar.M_ref != -1)
    {
        /* Mach has been defined, not velocity */
        avar.v_ref = avar.sos * avar.M_ref;
    }
    else
    {
        avar.M_ref = avar.v_ref / avar.sos;
    }

    if (svar.delta_t_min > 0)
        svar.delta_t = svar.delta_t_min;
    else
        svar.delta_t = 2E-010; /*Initial timestep*/

    fvar.H = fvar.H_fac * svar.particle_step;
    fvar.H_sq = fvar.H * fvar.H;
    fvar.sr = 4 * fvar.H_sq; /*KDtree search radius*/

    fvar.dsph_cont = 2.0 * fvar.dsph_delta * fvar.H * fvar.speed_sound;
    // fvar.dsph_mom = fvar.visc_alpha * fvar.H * fvar.speed_sound * fvar.rho_rest;
    fvar.dsph_mom = 2.0 * (SIMDIM + 2.0);
    fvar.nu = fvar.mu / fvar.rho_rest;

#if SIMDIM == 2
#ifdef CUBIC
    fvar.W_correc = 10.0 / (7.0 * M_PI * fvar.H * fvar.H);
#else
    fvar.W_correc = 7.0 / (4.0 * M_PI * fvar.H * fvar.H);
#endif
#endif
#if SIMDIM == 3
#ifdef CUBIC
    fvar.W_correc = (1.0 / (M_PI * fvar.H * fvar.H * fvar.H));
#else
    fvar.W_correc = (21 / (16 * M_PI * fvar.H * fvar.H * fvar.H));
#endif
#endif

    fvar.W_dx = Kernel(svar.particle_step, fvar.H, fvar.W_correc);

#if SIMDIM == 3
    if (svar.Asource == VLMInfl)
    {
        vortex.Init(svar.vlm_file);
        vortex.GetGamma(avar.v_inf);
        vortex.write_VLM_Panels(svar.output_prefix);
        if (vortex.write_traj)
            vortex.Plot_Streamlines(svar.output_prefix);
    }
#endif

    avar.GetYcoef(fvar, /*fvar.H*/ svar.particle_step);
    real n_full = get_n_full(svar.particle_step, fvar.H);
    avar.n_full = n_full;
    avar.i_n_full = 1.0 / avar.n_full;
    avar.interp_fac = 1.0 / avar.i_interp_fac;

#if SIMDIM == 3
    avar.A_plate = svar.particle_step * svar.particle_step;
    // avar.A_plate = fvar.H*fvar.H;
#else
    avar.A_plate = svar.particle_step /**svar.particle_step*/ /** pow(avar.L,0.5)*/;
    // avar.A_plate = fvar.H;
#endif

    /* Particle tracking values */
    svar.ipt_diam = pow((6.0 * fvar.sim_mass) / (M_PI * fvar.rho_rest), 1.0 / 3.0);
    svar.ipt_area = M_PI * svar.ipt_diam * svar.ipt_diam / 4.0;
}

void print_vector(string const& pretext, StateVecD const& vec)
{
#if SIMDIM == 3
    printf("%s: %f %f %f\n", pretext.c_str(), vec[0], vec[1], vec[2]);
#else
    printf("%s: %f %f\n", pretext.c_str(), vec[0], vec[1]);
#endif
}

void Print_Settings(FILE* out, SIM const& svar, FLUID const& fvar, AERO const& avar)
{
    fprintf(out, "Input values, after modification and interpretation, are: \n\n");

#pragma omp parallel
    {
#pragma omp single
        fprintf(out, "                         Number of threads: %d\n", svar.numThreads);
    }

    /* File Inputs */
    fprintf(out, " Input files --------------------------------: -\n");
    fprintf(out, "              Input fluid definition filename: %s\n", svar.input_fluid_file.c_str());
    fprintf(out, "           Input boundary definition filename: %s\n", svar.input_bound_file.c_str());
    fprintf(out, "                           SPH restart prefix: %s\n", svar.restart_prefix.c_str());
    fprintf(out, "                      VLM definition filename: %s\n\n", svar.vlm_file.c_str());

    /* TAU files */
    fprintf(out, " TAU files ----------------------------------: -\n");
    fprintf(out, "                   Primary grid face filename: %s\n", svar.tau_mesh.c_str());
    fprintf(out, "                    Boundary mapping filename: %s\n", svar.tau_bmap.c_str());
    fprintf(out, "                          Restart-data prefix: %s\n", svar.tau_sol.c_str());
    fprintf(out, "                                   Grid scale: %f\n", svar.scale);
    fprintf(out, "                         Angle alpha (degree): %f\n", svar.angle_alpha);
    fprintf(out, "           2D offset vector (0 / x=1,y=2,z=3): %d\n\n", svar.offset_axis);

    /* OpenFOAM files */
    fprintf(out, " OpenFOAM files -----------------------------: -\n");
    fprintf(out, "                     OpenFOAM input directory: %s\n", svar.foam_dir.c_str());
    fprintf(out, "                  OpenFOAM solution directory: %s\n", svar.foam_sol.c_str());
    fprintf(out, "                        OpenFOAM binary (0/1): %d\n", svar.foam_is_binary);
    fprintf(out, "                           Label size (32/64): %d\n", svar.foam_label_size);
    fprintf(out, "                          Scalar size (32/64): %d\n", svar.foam_scalar_size);
    fprintf(out, "                       OpenFOAM buoyant (0/1): %d\n\n", svar.foam_buoyant_sim);

    /* File outputs */
    fprintf(out, " Output parameters --------------------------: -\n");
    fprintf(out, "                 Single file for output (0/1): %u\n", svar.single_file);
    fprintf(out, "                   Write Tecplot output (0/1): %u\n", svar.write_tecio);
    fprintf(out, "                    Write H5Part output (0/1): %u\n", svar.write_h5part);
    fprintf(out, "                          Output files prefix: %s\n", svar.output_prefix.c_str());
    fprintf(out, "                      SPH frame time interval: %f\n", svar.frame_time_interval);
    fprintf(out, "                              SPH frame count: %u\n", svar.max_frames);
    fprintf(out, "       SPH output encoding (0=ascii/1=binary): %u\n", svar.out_encoding);
    fprintf(out, "                                Variable list: %s\n\n", svar.output_names.c_str());

    /* Fluid data */
    fprintf(out, " Fluid parameters ---------------------------: -\n");
    fprintf(out, "                            Reference density: %g\n", avar.rho_g);
    fprintf(out, "                  Reference dispersed density: %g\n", fvar.rho_rest);
    fprintf(out, "                Reference dispersed viscosity: %g\n", fvar.mu);
    fprintf(out, "                    Reference surface tension: %g\n", fvar.sig);
    fprintf(out, "            SPH surface tension contact angle: %g\n", fvar.contangb);
    fprintf(out, "              SPH artificial viscosity factor: %g\n", fvar.visc_alpha);
    fprintf(out, "                        SPH delta coefficient: %g\n", fvar.dsph_delta);
    fprintf(out, "                SPH maximum shifting velocity: %g\n", svar.max_shift_vel);
    fprintf(out, "                           SPH speed of sound: %g\n", fvar.speed_sound);
    fprintf(out, "                      SPH background pressure: %g\n", fvar.press_back);
    fprintf(out, "                        SPH starting pressure: %g\n", fvar.press_pipe);
    fprintf(out, "   SPH maximum absolute density variation (%%): %g\n", fvar.rho_var);
    fprintf(out, "                          SPH maximum density: %g\n", fvar.rho_max);
    fprintf(out, "                          SPH minimum density: %g\n", fvar.rho_min);
    fprintf(out, " SPH density variation to reduce timestep (%%): %g\n", fvar.rho_max_iter);
    fprintf(out, "              Init hydrostatic pressure (0/1): %d\n", svar.init_hydro_pressure);
    fprintf(out, "                           Hydrostatic height: %g\n\n", svar.hydro_height);

    /* Simulation settings */
    fprintf(out, " Time integration settings ------------------: -\n");
    fprintf(out, "                       SPH integration solver: %s\n", svar.solver_name.c_str());
    fprintf(out, " SPH compressibility solver (0=WCSPH,1=ACSPH): %s\n", svar.compressibility_solver);
    fprintf(out, "  SPH equation of state (0=Cole/1=Isothermal): %u\n", fvar.pressure_rel);
    fprintf(out, "     SPH boundary solver (0=pressure/1=ghost): %u\n", svar.bound_solver);
    fprintf(out, "                  SPH solver minimum residual: %g\n", svar.min_residual);
    fprintf(out, "                         SPH maximum timestep: %g\n", svar.delta_t_max);
    fprintf(out, "                         SPH minimum timestep: %g\n", svar.delta_t_min);
    fprintf(out, "                              SPH maximum CFL: %g\n", svar.cfl_max);
    fprintf(out, "                              SPH minimum CFL: %g\n", svar.cfl_min);
    fprintf(out, "                            SPH CFL condition: %g\n", svar.cfl);
    fprintf(out, "                        SPH unstable CFL step: %g\n", svar.cfl_step);
    fprintf(out, "             SPH Newmark-Beta iteration limit: %u\n", svar.max_subits);
    fprintf(out, "                 SPH unstable CFL count limit: %u\n", svar.n_unstable_limit);
    fprintf(out, "                   SPH stable CFL count limit: %u\n", svar.n_stable_limit);
    fprintf(out, "        SPH stable CFL count iteration factor: %g\n", svar.subits_factor);
    fprintf(out, "                   SPH maximum particle count: %zu\n\n", svar.max_points);

    fprintf(out, " Geometric parameters -----------------------: -\n");
#if SIMDIM == 3
    fprintf(
        out, "                           SPH gravity vector: %g, %g, %g\n", svar.grav[0], svar.grav[1],
        svar.grav[2]
    );
#else
    fprintf(out, "                           SPH gravity vector: %g, %g\n", svar.grav[0], svar.grav[1]);
#endif
    fprintf(out, "                          SPH initial spacing: %g\n", svar.particle_step);
    fprintf(out, "                  SPH boundary spacing factor: %g\n", svar.bound_step_factor);
    fprintf(out, "                    SPH restart fit tolerance: %f\n", svar.restart_tol);
    fprintf(out, "                  SPH smoothing length factor: %g\n", fvar.H_fac);
#if SIMDIM == 3
    fprintf(
        out, "                 SPH global offset coordinate: %g, %g, %g\n\n", svar.offset_vec[0],
        svar.offset_vec[1], svar.offset_vec[2]
    );
#else
    fprintf(
        out, "                 SPH global offset coordinate: %g, %g\n\n", svar.offset_vec[0],
        svar.offset_vec[1]
    );
#endif

    fprintf(out, " Aerodynamic coupling parameters ------------: -\n");
    fprintf(out, "                         SPH aerodynamic case: %s\n", avar.aero_case.c_str());
    fprintf(out, "                 SPH aerodynamic cutoff value: %f\n", avar.lam_cutoff);
    fprintf(out, "          SPH aerodynamic iterpolation factor: %f\n", avar.i_interp_fac);
    fprintf(out, "           SPH aerodynamic full neighbourhood: %f\n", avar.n_full);
    fprintf(out, " SPH interpolation factor (0=ncount/1=lambda): %d\n", avar.use_lam);
    fprintf(out, "        SPH SP diameter definition (0=dx/1=h): %d\n", avar.use_dx);
    fprintf(out, "                SPH use TAB deformation (0/1): %d\n", avar.use_TAB_def);
#if SIMDIM == 3
    fprintf(
        out, "                      SPH freestream velocity: %g, %g, %g\n", avar.v_inf[0], avar.v_inf[1],
        avar.v_inf[2]
    );
#else
    fprintf(
        out, "                      SPH freestream velocity: %g, %g\n\n", avar.v_inf[0], avar.v_inf[1]
    );
#endif

    /* Aerodynamic data */
    fprintf(out, " Carrier phase parameters -------------------: -\n");
    fprintf(out, "                           Reference velocity: %g\n", avar.v_ref);
    fprintf(out, "                           Reference pressure: %g\n", avar.p_ref);
    fprintf(out, "                        Reference Mach number: %g\n", avar.M_ref);
    fprintf(out, "                            Reference density: %g\n", avar.rho_g);
    fprintf(out, "                        Reference temperature: %g\n", avar.temp_g);
    fprintf(out, "                           Gas constant gamma: %g\n", avar.gamma);
    fprintf(out, "               Sutherland reference viscosity: %g\n\n", avar.mu_g);

    /* Particle tracking settings */
    fprintf(out, " Particle tracking parameters ---------------: -\n");
    fprintf(out, "                      Transition to IPT (0/1): %d\n", svar.using_ipt);
    fprintf(out, "                Velocity equation order (1/2): %d\n", svar.ipt_eq_order);
    fprintf(out, "         SPH tracking conversion x coordinate: %g\n", svar.max_x_sph);
    fprintf(out, "              Maximum x trajectory coordinate: %g\n", svar.max_x);
    fprintf(out, "              Particle scatter output (0/1/2): %u\n", svar.part_out);
    fprintf(out, "               Particle streak output (0/1/2): %u\n", svar.streak_out);
    fprintf(out, "    Particle cell intersection output (0/1/2): %u\n\n", svar.cells_out);

    fflush(out);
}

void GetInput(int argc, char** argv, SIM& svar, FLUID& fvar, AERO& avar, VLM& vortex)
{
    if (argc > 2)
    { /*Check number of input arguments*/
        printf("WARNING: only a maximum of one input arguments accepted,\n");
        printf("\t1: Input file directory\n");
        printf("Other inputs will be ignored.\n");
    }

    if (argc == 1)
    { /*Check if input has been provided*/
        printf("ERROR: No inputs provided. Stopping... \n");
        exit(-1);
    }

    /*Get parameters if it has been provided*/
    // printf(argv[1], endl;
    string file = argv[1];

    svar.input_file = file;

    FILE* fin = fopen(argv[1], "r");
    if (fin == NULL)
    {
        printf("ERROR: Failed to open parameter file: %s\n", argv[1]);
        exit(-1);
    }
    else
    {
        printf("Para file open, reading contents...\n");
    }

    char* cline = NULL;
    size_t cline_buf_size = 0;
    ssize_t cline_size = getline(&cline, &cline_buf_size, fin);
    while (cline_size >= 0)
    {
        string line(cline); // Recast as a string to manipulate.
        line = ltrim(line);
        size_t end = line.find_first_of('#');
        if (end != std::string::npos)
            line = line.substr(0, end + 1);

        /* File Inputs */
        Get_Number(line, "Number of threads", svar.numThreads);
        Get_String(line, "Input fluid definition filename", svar.input_fluid_file);
        Get_String(line, "Input boundary definition filename", svar.input_bound_file);
        Get_String(line, "Primary grid face filename", svar.tau_mesh);
        Get_String(line, "Boundary mapping filename", svar.tau_bmap);
        Get_String(line, "Restart-data prefix", svar.tau_sol);
        Get_String(line, "SPH restart prefix", svar.restart_prefix);
        Get_Number(line, "Grid scale", svar.scale);
        Get_Number(line, "Angle alpha (degree)", svar.angle_alpha);
        Get_Number(line, "2D offset vector (0 / x=1,y=2,z=3)", svar.offset_axis);
        Get_String(line, "OpenFOAM input directory", svar.foam_dir);
        Get_String(line, "OpenFOAM solution directory", svar.foam_sol);
        Get_Bool(line, "OpenFOAM binary (0/1)", svar.foam_is_binary);
        Get_Number(line, "Label size (32/64)", svar.foam_label_size);
        Get_Number(line, "Scalar size (32/64)", svar.foam_scalar_size);
        Get_Bool(line, "OpenFOAM buoyant (0/1)", svar.foam_buoyant_sim);
#if SIMDIM == 3
        Get_String(line, "VLM definition filename", svar.vlm_file);
#endif

        /* File outputs */
        Get_Bool(line, "Write Tecplot output (0/1)", svar.write_tecio);
        Get_Bool(line, "Write H5Part output (0/1)", svar.write_h5part);
        Get_Number(line, "Single file for output (0/1)", svar.single_file);
        Get_String(line, "Output files prefix", svar.output_prefix);
        Get_Number(line, "SPH frame time interval", svar.frame_time_interval);
        Get_Number(line, "SPH frame count", svar.max_frames);
        Get_Number(line, "SPH output encoding (0=ascii/1=binary)", svar.out_encoding);
        Get_String(line, "Variable list", svar.output_names);

        // Get_String(line, "Particle surface impact filename", svar.ascii_surface_file);

        /* Fluid data */
        Get_Number(line, "Reference dispersed density", fvar.rho_rest);
        Get_Number(line, "Sutherland reference viscosity", avar.mu_g);
        Get_Number(line, "Reference dispersed viscosity", fvar.mu);
        Get_Number(line, "Reference surface tension", fvar.sig);
        Get_Number(line, "SPH surface tension contact angle", fvar.contangb);
        Get_Number(line, "Init hydrostatic pressure (0/1)", svar.init_hydro_pressure);
        Get_Number(line, "Hydrostatic height", svar.hydro_height);

        /* Aerodynamic data */

        /* Simulation settings */
        Get_String(line, "SPH integration solver", svar.solver_name);
        Get_Number(line, "SPH compressibility solver (0=WCSPH,1=ACSPH)", svar.compressibility_solver);
        Get_Number(line, "SPH equation of state (0=Cole/1=Isothermal)", fvar.pressure_rel);
        Get_Number(line, "SPH boundary solver (0=pressure/1=ghost)", svar.bound_solver);
        Get_Number(line, "SPH solver minimum residual", svar.min_residual);
        Get_Number(line, "SPH maximum timestep", svar.delta_t_max);
        Get_Number(line, "SPH minimum timestep", svar.delta_t_min);
        Get_Number(line, "SPH maximum CFL", svar.cfl_max);
        Get_Number(line, "SPH minimum CFL", svar.cfl_min);
        Get_Number(line, "SPH CFL condition", svar.cfl);
        Get_Number(line, "SPH unstable CFL step", svar.cfl_step);
        Get_Number(line, "SPH unstable CFL count limit", svar.n_unstable_limit);
        Get_Number(line, "SPH stable CFL count limit", svar.n_stable_limit);
        Get_Number(line, "SPH stable CFL count iteration factor", svar.subits_factor);
        Get_Number(line, "SPH maximum shifting velocity", svar.max_shift_vel);

        Get_Number(line, "SPH background pressure", fvar.press_back);
        Get_Number(line, "SPH starting pressure", fvar.press_pipe);
        Get_Number(line, "SPH maximum absolute density variation (%)", fvar.rho_var);
        Get_Number(line, "SPH density variation to reduce timestep (%)", fvar.rho_max_iter);
        Get_Number(line, "SPH maximum density", fvar.rho_max);
        Get_Number(line, "SPH minimum density", fvar.rho_min);
        Get_Number(line, "SPH delta coefficient", fvar.dsph_delta);

        Get_Number(line, "SPH artificial viscosity factor", fvar.visc_alpha);
        Get_Number(line, "SPH speed of sound", fvar.speed_sound);
        Get_Number(line, "SPH Newmark-Beta iteration limit", svar.max_subits);
        Get_Vector(line, "SPH gravity vector", svar.grav);

        Get_Number(line, "SPH initial spacing", svar.particle_step);
        Get_Number(line, "SPH boundary spacing factor", svar.bound_step_factor);
        Get_Number(line, "SPH smoothing length factor", fvar.H_fac);
        Get_Vector(line, "SPH global offset coordinate", svar.offset_vec);
        Get_Number(line, "SPH maximum particle count", svar.max_points);
        Get_Number(line, "SPH restart fit tolerance", svar.restart_tol);

        /* Aerodynamic coupling settings */
        Get_String(line, "SPH aerodynamic case", avar.aero_case);
        Get_Number(line, "SPH SP diameter definition (0=dx/1=h)", avar.use_dx);
        Get_Number(line, "SPH use TAB deformation (0/1)", avar.use_TAB_def);
        Get_Number(line, "SPH interpolation factor (0=ncount/1=lambda)", avar.use_lam);
        Get_Number(line, "SPH aerodynamic cutoff value", avar.lam_cutoff);
        Get_Number(line, "SPH aerodynamic interpolation factor", avar.i_interp_fac);
        Get_Vector(line, "SPH freestream velocity", avar.v_inf);
        Get_Number(line, "Reference velocity", avar.v_ref);
        Get_Number(line, "Reference pressure", avar.p_ref);
        Get_Number(line, "Reference Mach number", avar.M_ref);
        Get_Number(line, "Reference density", avar.rho_g);
        Get_Number(line, "Reference temperature", avar.temp_g);
        Get_Number(line, "Gas constant gamma", avar.gamma);

        /* Particle tracking settings */
        Get_Number(line, "Transition to IPT (0/1)", svar.using_ipt);
        Get_Number(line, "Velocity equation order (1/2)", svar.ipt_eq_order);
        Get_Number(line, "SPH tracking conversion x coordinate", svar.max_x_sph);
        Get_Number(line, "Maximum x trajectory coordinate", svar.max_x);
        Get_Number(line, "Particle scatter output (0/1/2)", svar.part_out);
        Get_Number(line, "Particle streak output (0/1/2)", svar.streak_out);
        Get_Number(line, "Particle cell intersection output (0/1/2)", svar.cells_out);

        cline_size = getline(&cline, &cline_buf_size, fin);
    }

    free(cline);
    cline = NULL;
    fclose(fin);

    int fault = 0;

    /* Need to check if inputs are correct */
    if (svar.tau_mesh.empty())
    {
        if (svar.foam_dir.empty())
        {
#if SIMDIM == 3
            if (svar.vlm_file.empty())
                svar.Asource = constVel;
            else
            {
                svar.Asource = VLMInfl;

                if (svar.vlm_file == "(thisfile)")
                    svar.vlm_file = svar.input_file;
            }

#else
            svar.Asource = constVel;
#endif
        }
        else
        {
            if (svar.foam_sol.empty())
            {
                printf("OpenFOAM solution directory not defined.\n");
                fault = 1;
            }

            svar.mesh_source = OpenFOAM;
            svar.Asource = meshInfl;
        }
    }
    else
    {
        svar.mesh_source = TAU_CDF;
        svar.Asource = meshInfl;
        if (svar.tau_bmap.empty())
        {
            printf("Input TAU bmap file not defined.\n");
            fault = 1;
        }
        else if (svar.tau_bmap == "(thisfile)")
        {
            /* The para file is the boundary definition file, so set it so */
            svar.tau_bmap = svar.input_file;
        }
        else
        {
            // Check files exist
            if (!std::filesystem::exists(svar.tau_bmap))
            {
                printf("Input TAU boundary file \"%s\" not found.\n", svar.tau_bmap.c_str());
                fault = 1;
            }

            if (!std::filesystem::exists(svar.tau_mesh))
            {
                printf("Input TAU mesh file \"%s\" not found.\n", svar.tau_mesh.c_str());
                fault = 1;
            }
        }

        if (svar.tau_sol.empty())
        {
            printf("Input TAU solution file not defined.\n");
            fault = 1;
        }
        else
        {
            if (!std::filesystem::exists(svar.tau_sol))
            {
                printf("Input TAU solution file \"%s\" not found.\n", svar.tau_sol.c_str());
                fault = 1;
            }
        }
    }

    /* Check for restart data now boundary case is known */
    if (!svar.restart_prefix.empty())
    {
        svar.restart = 1;

        if (!std::filesystem::exists(svar.restart_prefix + "_particles.h5"))
        {
            printf(
                "SPH restart file \"%s\" not found.\n", (svar.restart_prefix + "_particles.h5").c_str()
            );
            fault = 1;
        }
        // Check_If_Restart_Possible(svar);
        // svar.output_prefix = svar.restart_prefix;
    }

    if (svar.input_fluid_file.empty())
    {
        printf("Input SPH fluid file not defined.\n");
        fault = 1;
    }
    else
    {
        if (!std::filesystem::exists(svar.input_fluid_file))
        {
            printf("SPH fluid file \"%s\" not found.\n", svar.input_fluid_file.c_str());
            fault = 1;
        }
    }

    if (svar.input_bound_file.empty())
    {
        printf("Input SPH boundary file not defined.\n");
        fault = 1;
    }
    else
    {
        if (!std::filesystem::exists(svar.input_bound_file))
        {
            printf("SPH boundary file \"%s\" not found.\n", svar.input_bound_file.c_str());
            fault = 1;
        }
    }

    if (!svar.solver_name.empty())
    {
        if (svar.solver_name == "Newmark-Beta")
        {
            svar.solver_type = newmark_beta;
        }
        else if (svar.solver_name == "Runge-Kutta")
        {
            svar.solver_type = runge_kutta;
        }
        else
        {
            printf("ERROR: Unrecognised solver name. Choose from the following options:\n");
            printf("\t1. Newmark-Beta\n\t2. Runge-Kutta\n");
            fault = 1;
        }
    }

    if (svar.particle_step < 0)
    {
        printf("ERROR: SPH initial spacing has not been defined.\n");
        fault = 1;
    }

    if (svar.max_frames < 0)
    {
        printf("Number of frames to output not defined.\n");
        fault = 1;
    }

    if (svar.frame_time_interval < 0)
    {
        printf("Frame time interval has not been defined.\n");
        fault = 1;
    }

    if (svar.init_hydro_pressure)
    {
        if (svar.hydro_height < 0)
        {
            printf("Hydrostatic height has not been defined.\n");
            fault = 1;
        }
    }

    /* Aerodynamic settings */
    if (avar.aero_case == "(none)")
    {
        avar.acase = NoAero;
    }
    else if (avar.aero_case == "Gissler")
    {
        avar.acase = Gissler;
    }
    else if (avar.aero_case == "Induced_pressure")
    {
        avar.acase = InducedPressure;
    }
    else if (avar.aero_case == "Skin_friction")
    {
        avar.acase = SkinFric;
    }
    else
    {
        printf("Aerodynamic coupling model is not defined or correct.\n");
        fault = 1;
    }

    if (svar.delta_t_min > 0.0)
    {
        svar.delta_t = svar.delta_t_min;
    }

    if (avar.i_interp_fac < 0 || avar.i_interp_fac > 1.0)
    {
        printf("ERROR: SPH aerodynamic interpolation factor must be between 0 < x < 1.\n");
        fault = 1;
    }

    /* Particle Tracking Settings */
    if (svar.using_ipt)
    {
        if (svar.ipt_eq_order > 2 || svar.ipt_eq_order < 1)
        {
            printf("Equation order not 1 or 2. Please choose between these.\n");
            exit(-1);
        }

        if (svar.max_x < svar.max_x_sph)
        {
            printf("WARNING: Maximum x coordinate for particle tracking is less than that for SPH.\n");
            printf("         No particle tracking will be performed.\n");
            svar.using_ipt = 0;
        }

        // svar.streakdir = svar.outfile;
        // size_t pos = svar.streak_file.find_last_of(".");
        // svar.streak_file.insert(pos,"_streak");
    }

    if (svar.offset_axis != no_offset)
    {
        if (SIMDIM != 2)
        {
            printf("WARNING: trying to use 3D code with a 2D settings file.\n");
            svar.offset_axis = no_offset;
        }
    }
    else if (svar.offset_axis > z_axis)
    {
        printf("ERROR: 2D offset axis option out of bounds\n");
        fault = 1;
    }

#if SIMDIM == 2
    if (svar.offset_axis == no_offset)
    {
        printf("ERROR: Offset axis has not been defined.\n");
        fault = 1;
    }
#endif

    if (fault)
    {
        printf("Check of the input settings finished with errors. Stopping\n");
        exit(-1);
    }

    Check_Output_Variables(svar);

    Set_Values(svar, fvar, avar, vortex);

    Print_Settings(stdout, svar, fvar, avar);
#ifdef DEBUG
    Print_Settings(dbout, svar, fvar, avar);
#endif

} /*End of GetInput()*/

void Write_Tec_Headers(
    FILE* ff, FILE* fb, FILE* fg, SIM& svar, FLUID const& fvar, AERO const& avar,
    std::string const& prefix
)
{
    if (svar.out_encoding == 1)
    {
        Init_Binary_PLT(
            svar, fvar, avar, prefix, "_fluid.szplt", "Simulation Particles", svar.output_fluid_file
        );

        Init_Binary_PLT(
            svar, fvar, avar, prefix, "_boundary.szplt", "Boundary Particles", svar.output_bound_file
        );
    }
    else
    {
        /* Write first timestep */
        string mainput_file = prefix;
        mainput_file.append("_fluid.dat");
        ff = fopen(mainput_file.c_str(), "a");
        if (ff != NULL)
            Write_ASCII_header(ff, svar, "Simulation Particles");
        else
        {
            printf("Failed to open %s. Stopping.\n", mainput_file.c_str());
            exit(-1);
        }

        /*If the boundary exists, write it.*/
        string h5part_bound_file = prefix;
        h5part_bound_file.append("_boundary.dat");
        fb = fopen(h5part_bound_file.c_str(), "a");
        if (fb != NULL)
            Write_ASCII_header(fb, svar, "Boundary Particles");
        else
        {
            printf("Error opening %s file. Stopping\n", h5part_bound_file.c_str());
            exit(-1);
        }
    }
}

void Write_h5part_Headers(SIM& svar, FLUID const& fvar, AERO const& avar, std::string const& prefix)
{
    open_h5part_files(svar, fvar, avar, prefix, svar.h5part_fluid_file, svar.h5part_bound_file);
}

void Write_Timestep(
    FILE* ff, FILE* fb, FILE* fg, SIM& svar, FLUID const& fvar, AERO const& avar, LIMITS const& limits,
    SPHState const& pnp1
)
{
    if (svar.write_tecio)
    {
        std::ostringstream oss;
        oss << svar.current_time;
        if (!svar.single_file)
        {
            string prefix = svar.output_prefix + "_time_" + oss.str();
            Write_Tec_Headers(ff, fb, fg, svar, fvar, avar, prefix);
        }

        if (svar.out_encoding == 1)
        {
            for (size_t bound = 0; bound < svar.n_bound_blocks; bound++)
            {
                string title = "Boundary_" + std::to_string(bound) + "_" + limits[bound].name +
                               "_time_" + oss.str() + "s";
                Write_Binary_Timestep(
                    svar, fvar.rho_rest, pnp1, limits[bound], title.c_str(),
                    static_cast<int32_t>(bound + 1), svar.output_bound_file
                );
            }

            for (size_t block = svar.n_bound_blocks; block < svar.n_bound_blocks + svar.n_fluid_blocks;
                 block++)
            {
                string title = "Fluid_" + std::to_string(block) + "_" + limits[block].name + "_time_" +
                               oss.str() + "s";
                Write_Binary_Timestep(
                    svar, fvar.rho_rest, pnp1, limits[block], title.c_str(),
                    static_cast<int32_t>(block + 1), svar.output_fluid_file
                );
            }
        }
        else
        {
            for (size_t bound = 0; bound < svar.n_bound_blocks; bound++)
            {
                string title = "Boundary_" + std::to_string(bound) + "_" + limits[bound].name +
                               "_time_" + oss.str() + "s";
                Write_ASCII_Timestep(
                    svar, fvar.rho_rest, pnp1, limits[bound].index.first, limits[bound].index.second,
                    title.c_str(), bound + 1, fb
                );
            }

            for (size_t block = svar.n_bound_blocks; block < svar.n_bound_blocks + svar.n_fluid_blocks;
                 block++)
            {
                string title = "Fluid_" + std::to_string(block) + "_" + limits[block].name + "_time_" +
                               oss.str() + "s";

                Write_ASCII_Timestep(
                    svar, fvar.rho_rest, pnp1, limits[block].index.first, limits[block].index.second,
                    title.c_str(), block + 1, ff
                );
            }
        }

        if (!svar.single_file)
        { // Close files after use
            if (ff != NULL)
                fclose(ff);
            if (fb != NULL)
                fclose(fb);
            if (fg != NULL)
                fclose(fg);

            if (svar.out_encoding == 1)
            {
                /*Combine the szplt files*/
                close_file(&svar.output_fluid_file);
                close_file(&svar.output_bound_file);
            }
        }
    }

    if (svar.write_h5part)
    {
        string prefix = svar.output_prefix;
        if (!svar.single_file)
        {
            std::ostringstream oss;
            oss << svar.current_time;
            prefix += "_" + std::to_string(svar.current_frame);
        }

        Write_h5part_Headers(svar, fvar, avar, prefix);
        write_h5part_data(svar.h5part_fluid_file, svar.h5part_bound_file, svar, fvar, pnp1);
        // Files are closed in the write function
    }
}

void Remove_Old_Files(SIM const& svar)
{
    string dir;
    string prefix = svar.output_prefix;
    if (svar.output_prefix.find_last_of("/\\") != string::npos)
    {
        dir = svar.output_prefix.substr(0, svar.output_prefix.find_last_of("/\\"));
        // prefix = svar.output_prefix.substr(svar.output_prefix.find_last_of("/\\")+1);
    }
    else
    {
        dir = ".";
        prefix = dir + "/" + svar.output_prefix;
    }

    for (int ii = prefix.size() - 1; ii >= 0; ii--)
    {
        if (prefix[ii] == '.')
            prefix.insert(ii, "\\");

        if (prefix[ii] == '/')
            prefix.insert(ii, "\\");
    }

    int fault = 0;
    if (svar.write_tecio)
    {
        std::regex file_expr("(" + prefix + ")(.*)(\\.szplt\\.sz)(.*)$");
        // std::regex file_expr("(" + prefix + ")(.*)(\\.szplt)$");

        for (auto const& file : directory_iterator(dir))
        {
            // printf("%s: ",file.path().string().c_str());
            // cout << file.path().string() << endl;
            // Check it has the right prefix
            if (std::regex_match(file.path().string(), file_expr))
            {
                // printf("File matches the expression \n");
                // Remove the file
                try
                {
                    std::filesystem::remove(file.path());
                }
                catch (const std::filesystem::__cxx11::filesystem_error& e)
                {
                    printf("Could not remove tecplot file. Trying to remove other files.\n");
                    fault = 1;
                }
            }
            // else
            // 	printf("File doesn't match the expression\n");
        }
    }

    if (fault)
    {
        printf("Failed to remove some old simulation files, so cannot continue.\n");
        exit(-1);
    }

    if (svar.write_h5part)
    {
        std::regex fluid_expr("(" + prefix + ")(.*)(\\.h5part)");

        for (auto const& file : directory_iterator(dir))
        {
            // printf("%s: ",file.path().filename().string().c_str());
            // cout << file.path().string() << endl;
            // Check it has the right prefix

            if (std::regex_match(file.path().string(), fluid_expr))
            {
                // printf("File matches the expression \n");
                // Remove the file
                try
                {
                    std::filesystem::remove(file.path());
                }
                catch (const std::filesystem::__cxx11::filesystem_error& e)
                {
                    printf("Could not remove fluid h5part file. Trying to remove other files.\n");
                    fault = 1;
                }
            }
            // else
            // 	printf("File doesn't match the expression\n");
        }
    }

    if (fault)
    {
        printf("Failed to remove some files. Try closing paraview, as it doesn't correctly free the "
               "files.\n");
        exit(-1);
    }
}

void Check_Output_Variables(SIM& svar)
{
    OutputMap& outvars = svar.output_variables;
    // Mandatory output parameters.
    outvars.insert({"pos-vec", OutputVariable("X", realType, true, true)});
    outvars.insert({"vel-vec", OutputVariable("V", realType, true, true)});
    outvars.insert({"acc-vec", OutputVariable("A", realType, true, true)});
    outvars.insert({"press", OutputVariable("Pressure", realType, true, false)});
    outvars.insert({"dRho", OutputVariable("dRho", realType, true, false)});
    outvars.insert({"part_id", OutputVariable("part_id", int32Type, true, false)});
    outvars.insert({"cellID", OutputVariable("cellID", int32Type, true, false)});
    outvars.insert({"bound", OutputVariable("bound", uint8Type, true, false)});

    // Start of optional parameters.
    outvars.insert({"dens", OutputVariable("Density", realType, false, false)});
    outvars.insert({"densVar", OutputVariable("Density-Variation", realType, false, false)});
    outvars.insert({"vmag", OutputVariable("V-mag", realType, false, false)});
    outvars.insert({"surf", OutputVariable("Surface", realType, false, false)});
    outvars.insert({"surfZ", OutputVariable("Surface-Zone", realType, false, false)});
    outvars.insert({"aero-mag", OutputVariable("Aero-mag", realType, false, false)});
    outvars.insert({"aero-vec", OutputVariable("Aero-vector", realType, false, true)});
    outvars.insert({"curv", OutputVariable("curvature", realType, false, false)});
    outvars.insert({"occl", OutputVariable("occl", realType, false, false)});
    outvars.insert({"cellP", OutputVariable("cell-pressure", realType, false, false)});
    outvars.insert({"cellRho", OutputVariable("cell-density", realType, false, false)});
    outvars.insert({"cellV-mag", OutputVariable("cell-vel-mag", realType, false, false)});
    outvars.insert({"cellV-vec", OutputVariable("cell-vel-vector", realType, false, true)});
    outvars.insert({"dsphG-vec", OutputVariable("Dsph-grad-vector", realType, false, true)});
    outvars.insert({"lam", OutputVariable("lambda", realType, false, false)});
    outvars.insert({"lam-nb", OutputVariable("lambda-no-bound", realType, false, false)});
    outvars.insert({"colour", OutputVariable("colour", realType, false, false)});
    outvars.insert({"colour-G", OutputVariable("colour-grad-vector", realType, false, false)});
    outvars.insert({"norm-vec", OutputVariable("surf-normal-vector", realType, false, true)});
    outvars.insert({"shiftV-mag", OutputVariable("shift-vel-mag", realType, false, false)});
    outvars.insert({"shiftV-vec", OutputVariable("shift-vel-vector", realType, false, true)});

    if (!svar.output_names.empty())
    {
        // Split the names using underscore as a delimiter
        std::vector<std::string> vals;
        size_t start = 0;
        size_t end = svar.output_names.find_first_of("_");
        while (end != std::string::npos)
        {
            vals.emplace_back(svar.output_names.substr(start, end - start));
            start = end + 1;
            end = svar.output_names.find_first_of("_", start);
        }
        vals.emplace_back(svar.output_names.substr(start)); // catch the final value

        // Now check the variables for which ones want to be output
        for (std::string const& var : vals)
        {
            if (outvars.find(var) != outvars.end())
            {
                // Variable exists, so mark it to be output
                outvars.at(var).write = true;
            }
            else
            {
                std::printf("Unrecognised output variable \"%s\" defined.\n", var.c_str());
            }
        }
    }

    // Fill the var_names string and var_types vector
    bool write_x = true;
    bool write_y = true;
    bool write_z = true;
#if SIMDIM == 2
    switch (svar.offset_axis)
    {
    case x_axis:
        write_x = false;
        break;
    case y_axis:
        write_y = false;
        break;
    case z_axis:
        write_z = false;
        break;
    }
#endif

    svar.var_names.clear();
    svar.var_types.clear();
    for (auto const& [key, var] : outvars)
    {
        if (var.write)
        {
            if (var.is_vector)
            {
                // Add variables for a vector
                if (write_x)
                {
                    svar.var_names += var.output_name + "-x,";
                    svar.var_types.emplace_back(var.data_type);
                }
                if (write_y)
                {
                    svar.var_names += var.output_name + "-y,";
                    svar.var_types.emplace_back(var.data_type);
                }
                if (write_z)
                {
                    svar.var_names += var.output_name + "-z,";
                    svar.var_types.emplace_back(var.data_type);
                }
            }
            else
            {
                svar.var_names += var.output_name + ",";
                svar.var_types.emplace_back(var.data_type);
            }
        }
    }

    // Remove the trailing comma
    if (svar.var_names.back() == ',')
        svar.var_names.pop_back();
}
