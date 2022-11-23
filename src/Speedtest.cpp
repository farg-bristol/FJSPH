

/* Library to perform drag assessment for a range of particle resolutions */
#include "Speedtest.h"
#include "Add.h"
#include "Aero.h"
#include "CDFIO.h"
#include "Containment.h"
#include "FOAMIO.h"
#include "Init.h"
#include "Integration.h"
#include "IOFunctions.h"
#include "IPT.h"
#include "Neighbours.h"
#include "Shifting.h"
#include <chrono>

using namespace std::chrono;

void Make_Mesh(SIM& svar, AERO const& avar, MESH& cells, size_t const& ni)
{
    real length = 10.0;
    #if SIMDIM == 3
    // Create the new mesh
    size_t nj = 2;
    size_t nk = 2;
    real dx = length/(ni-1);
    real dy = 1.0/(nj-1);
    real dz = 1.0/(nk-1);
    
    size_t npnts = ni*nj*nk;
    size_t nface = nk*(ni-1)*(nj-1) + ni*(nj-1)*(nk-1)+ nj*(ni-1)*(nk-1);
    size_t nelem = (ni-1)*(nj-1)*(nk-1);

    vector<StateVecD> pVel;
    vector<real> pP;
    vector<real> pRho;

    cells.alloc(npnts,nelem,nface,0,1);
    cells.faces = vector<vector<size_t>>(cells.nFace,vector<size_t>(4));
    cells.leftright = vector<std::pair<int,int>>(cells.nFace);

    cells.maxlength = 2.0*std::max(dx,dy);
    cells.minlength = std::min(dx,dy);


    cout << "Building vertices..." << endl;
    #pragma omp parallel for
    for(size_t ii = 0; ii < cells.nElem; ++ii)
    {
        cells.cVel[ii] = avar.vInf;
        cells.cP[ii] = avar.pRef;
        cells.cRho[ii] = avar.rhog;
    }

    // cout << "Building vertices..." << endl;
    #pragma omp parallel for
    for(size_t kk = 0; kk < nk; kk++)
    {
        for(size_t jj = 0; jj < nj; jj++)
        {
            for(size_t ii = 0; ii < ni; ii++)
            {
                cells.verts[index(ii,jj,kk,ni,nj)][0] = ii*dx;
                cells.verts[index(ii,jj,kk,ni,nj)][1] = jj*dy;
                cells.verts[index(ii,jj,kk,ni,nj)][2] = kk*dz;
            }
        }
    } 

    // Create the faces
    cout << "Building faces..." << endl;
    size_t faceno = 0;

    // Do the front and back faces
    for(size_t jj = 0; jj < nj-1; jj++)
    {
        for(size_t ii = 0; ii < ni-1; ii++)
        {
            cells.faces[faceno][0] = index(ii,jj,0,ni,nj);
            cells.faces[faceno][1] = index(ii+1,jj,0,ni,nj);
            cells.faces[faceno][2] = index(ii+1,jj+1,0,ni,nj);
            cells.faces[faceno][3] = index(ii,jj+1,0,ni,nj);
            cells.leftright[faceno].first = static_cast<int>(index(ii,jj,0,ni-1,nj-1));
            cells.leftright[faceno].second = -2;
            cells.nSurf++;
            faceno++;

            cells.faces[faceno][0] = index(ii,jj,nk-1,ni,nj);
            cells.faces[faceno][1] = index(ii,jj+1,nk-1,ni,nj);
            cells.faces[faceno][2] = index(ii+1,jj+1,nk-1,ni,nj);
            cells.faces[faceno][3] = index(ii+1,jj,nk-1,ni,nj);
            cells.leftright[faceno].first = static_cast<int>(index(ii,jj,nk-2,ni-1,nj-1));
            cells.leftright[faceno].second = -2;
            cells.nSurf++;
            faceno++;
        }
    }

    // Do the bottom and top faces
    for(size_t kk = 0; kk < nk-1; kk++)
    {
        for(size_t ii = 0; ii < ni-1; ii++)
        {
            cells.faces[faceno][0] = index(ii,0,kk,ni,nj);
            cells.faces[faceno][1] = index(ii+1,0,kk,ni,nj);
            cells.faces[faceno][2] = index(ii+1,0,kk+1,ni,nj);
            cells.faces[faceno][3] = index(ii,0,kk+1,ni,nj);
            cells.leftright[faceno].first = static_cast<int>(index(ii,0,kk,ni-1,nj-1));
            cells.leftright[faceno].second = -2;
            cells.nSurf++;
            faceno++;

            cells.faces[faceno][0] = index(ii,nj-1,kk,ni,nj);
            cells.faces[faceno][1] = index(ii,nj-1,kk+1,ni,nj);
            cells.faces[faceno][2] = index(ii+1,nj-1,kk+1,ni,nj);
            cells.faces[faceno][3] = index(ii+1,nj-1,kk,ni,nj);
            cells.leftright[faceno].first = static_cast<int>(index(ii,nj-2,kk,ni-1,nj-1));
            cells.leftright[faceno].second = -2;
            cells.nSurf++;
            faceno++;
        }
    }


    // Do the left and right faces
    for(size_t jj = 0; jj < nj-1; jj++)
    {
        for(size_t kk = 0; kk < nk-1; kk++)
        {
            cells.faces[faceno][0] = index(0,jj,kk,ni,nj);
            cells.faces[faceno][1] = index(0,jj+1,kk,ni,nj);
            cells.faces[faceno][2] = index(0,jj+1,kk+1,ni,nj);
            cells.faces[faceno][3] = index(0,jj,kk+1,ni,nj);
            cells.leftright[faceno].first = static_cast<int>(index(0,jj,kk,ni-1,nj-1));
            cells.leftright[faceno].second = -2;
            cells.nSurf++;
            faceno++;

            cells.faces[faceno][0] = index(ni-1,jj,kk,ni,nj);
            cells.faces[faceno][1] = index(ni-1,jj,kk+1,ni,nj);
            cells.faces[faceno][2] = index(ni-1,jj+1,kk+1,ni,nj);
            cells.faces[faceno][3] = index(ni-1,jj+1,kk,ni,nj);
            cells.leftright[faceno].first = static_cast<int>(index(ni-2,jj,kk,ni-1,nj-1));
            cells.leftright[faceno].second = -2;
            cells.nSurf++;
            faceno++;
        }
    }

    

    // Do faces in i
    for(size_t ii = 1; ii < ni-1; ii++)
    {
        for(size_t jj = 0; jj < nj-1; jj++)
        {
            for(size_t kk = 0; kk < nk-1; kk++)
            {
                cells.faces[faceno][0] = index(ii,jj,kk,ni,nj);
                cells.faces[faceno][1] = index(ii,jj+1,kk,ni,nj);
                cells.faces[faceno][2] = index(ii,jj+1,kk+1,ni,nj);
                cells.faces[faceno][3] = index(ii,jj,kk+1,ni,nj);
                cells.leftright[faceno].first = static_cast<int>(index(ii-1,jj,kk,ni-1,nj-1));
                cells.leftright[faceno].second = static_cast<int>(index(ii,jj,kk,ni-1,nj-1));
                faceno++;
            }
        }
    }

    // Do faces in j
    for(size_t jj = 1; jj < nj-1; jj++)
    {
        for(size_t ii = 0; ii < ni-1; ii++)
        {
            for(size_t kk = 0; kk < nk-1; kk++)
            {
                cells.faces[faceno][0] = index(ii,jj,kk,ni,nj);
                cells.faces[faceno][1] = index(ii+1,jj,kk,ni,nj);
                cells.faces[faceno][2] = index(ii+1,jj,kk+1,ni,nj);
                cells.faces[faceno][3] = index(ii,jj,kk+1,ni,nj);
                cells.leftright[faceno].first = static_cast<int>(index(ii,jj-1,kk,ni-1,nj-1));
                cells.leftright[faceno].second = static_cast<int>(index(ii,jj,kk,ni-1,nj-1));
                faceno++;
            }
        }
    }

    // Do faces in k
    for(size_t kk = 1; kk < nk-1; kk++)
    {
        for(size_t ii = 0; ii < ni-1; ii++)
        {
            for(size_t jj = 0; jj < nj-1; jj++)
            {
                cells.faces[faceno][0] = index(ii,jj,kk,ni,nj);
                cells.faces[faceno][1] = index(ii+1,jj,kk,ni,nj);
                cells.faces[faceno][2] = index(ii+1,jj+1,kk,ni,nj);
                cells.faces[faceno][3] = index(ii,jj+1,kk,ni,nj);
                cells.leftright[faceno].first = static_cast<int>(index(ii,jj,kk-1,ni-1,nj-1));
                cells.leftright[faceno].second = static_cast<int>(index(ii,jj,kk,ni-1,nj-1));
                faceno++;
            }
        }
    }

    #endif

    #if SIMDIM == 2
    // Create the new mesh
    size_t nj = 2;
    real dx = length/(ni-1);
    real dz = 1.0/(nj-1);

    // 
    vector<StateVecD> pVel;
    vector<real> pP;
    vector<real> pRho;

    cells.alloc(ni*nj, (ni-1)*(nj-1),2*ni*nj - ni - nj,0,1);
    cells.faces = vector<vector<size_t>>(cells.nFace,vector<size_t>(2));
    cells.leftright = vector<std::pair<int,int>>(cells.nFace);

    cells.maxlength = 2.0*std::max(dx,dz);
    cells.minlength = std::min(dx,dz);

    cout << "Building vertices..." << endl;
    for(size_t jj = 0; jj < nj; jj++)
    {
        for(size_t ii = 0; ii < ni; ii++)
        {
            cells.verts[index(ii,jj,ni)][0] = ii*dx;
            cells.verts[index(ii,jj,ni)][1] = jj*dz;
        }
    }

    for(size_t ii = 0; ii < cells.nElem; ++ii)
    {
        cells.cVel[ii] = avar.vInf;
        cells.cP[ii] = avar.pRef;
        cells.cRho[ii] = avar.rhog;
    }

    // Create the faces
    
    size_t faceno = 0;
    cout << "Building faces..." << endl;
    // Do the top and bottom edges
    for(size_t ii = 0; ii < ni-1; ii++)
    {
        cells.faces[faceno][0] = index(ii,0,ni);
        cells.faces[faceno][1] = index(ii+1,0,ni);
        cells.leftright[faceno].first = static_cast<int>(index(ii,0,ni-1));
        cells.leftright[faceno].second = -2;
        cells.nSurf++;
        faceno++;

        
        
        cells.faces[faceno][0] = index(ii+1,nj-1,ni);
        cells.faces[faceno][1] = index(ii,nj-1,ni);
        cells.leftright[faceno].first = static_cast<int>(index(ii,nj-2,ni-1));
        cells.leftright[faceno].second = -2;
        cells.nSurf++;
        faceno++;				
    }

    // Do the left & right edges
    for(size_t jj = 0; jj < nj-1; jj++)
    {
        cells.faces[faceno][0] = index(0,jj+1,ni);
        cells.faces[faceno][1] = index(0,jj,ni);
        cells.leftright[faceno].first = static_cast<int>(index(0,jj,ni-1));
        cells.leftright[faceno].second = -2;
        cells.nSurf++;
        faceno++;

        cells.faces[faceno][0] = index(ni-1,jj,ni);
        cells.faces[faceno][1] = index(ni-1,jj+1,ni);
        cells.leftright[faceno].first = static_cast<int>(index(ni-2,jj,ni-1));
        cells.leftright[faceno].second = -2;
        cells.nSurf++;
        faceno++;
                    
    }

    // Do the horizontal faces
    for(size_t jj = 1; jj < nj-1; jj++)
        for(size_t ii = 0; ii < ni-1; ii++)
        {
            cells.faces[faceno][0] = index(ii,jj,ni);
            cells.faces[faceno][1] = index(ii+1,jj,ni);
            cells.leftright[faceno].first = static_cast<int>(index(ii,jj,ni-1));
            cells.leftright[faceno].second = static_cast<int>(index(ii,jj-1,ni-1));
            faceno++;
        }

    // Do the vertical faces
    for(size_t ii = 1; ii < ni-1; ii++)
        for(size_t jj = 0; jj < nj-1; jj++)
        {
            cells.faces[faceno][0] = index(ii,jj,ni);
            cells.faces[faceno][1] = index(ii,jj+1,ni);
            cells.leftright[faceno].first = static_cast<int>(index(ii-1,jj,ni-1));
            cells.leftright[faceno].second = static_cast<int>(index(ii,jj,ni-1));
            faceno++;
        }
    #endif


    // Build cell faces
    for(size_t ii = 0; ii < cells.nFace; ++ii)
    {
        int const& left = cells.leftright[ii].first;
        int const& right = cells.leftright[ii].second;
        if(left >=0)
            cells.cFaces[left].emplace_back(ii);

        if(right >= 0)
            cells.cFaces[right].emplace_back(ii);
    }

    // Average the cell centres
    Average_Point_to_Cell(cells.verts,cells.cCentre,cells.cFaces,cells.faces);
}

void copy_svar(SIM& svar, SIM& svar_)
{
    /* Input files */
    svar_.CDForFOAM = svar.CDForFOAM;
    svar_.isBinary = svar.isBinary; 
    svar_.buoyantSim = svar.buoyantSim; 
    svar_.incomp = svar.incomp;
    svar_.labelSize = svar.labelSize; 
    svar_.scalarSize = svar.scalarSize;
    svar_.offset_axis = svar.offset_axis;
    svar_.angle_alpha = svar.angle_alpha;
    svar_.out_encoding = svar.out_encoding;
    svar_.outvar = svar.outvar; 
    svar_.var_types = svar.var_types; 
    svar_.gout = svar.gout;

    svar_.restart = svar.restart;
    svar_.scale = svar.scale;

    svar_.partID = svar.partID;
    svar_.simPts = svar.simPts;
    svar_.bndPts = svar.bndPts; 
    svar_.totPts = svar.totPts;
    svar_.gstPts = svar.gstPts;  
    svar_.psnPts = svar.psnPts; 
    svar_.delNum = svar.delNum; 
    svar_.intNum = svar.intNum;	
    svar_.nrefresh = svar.nrefresh;
    svar_.addcount = svar.addcount;
    svar_.finPts = svar.finPts;
    
    svar_.Pstep = svar.Pstep; 
    svar_.Bstep = svar.Bstep;
    svar_.dx = svar.dx; 
    
    svar_.start_type = svar.start_type; 
    svar_.bound_type = svar.bound_type;
    svar_.Scase = svar.Scase; 
    svar_.Bcase = svar.Bcase;
    svar_.sim_start = svar.sim_start;
    svar_.bound_start = svar.bound_start;
    svar_.Bclosed = svar.Bclosed; 
    svar_.ghost = svar.ghost;
    svar_.Asource = svar.Asource;
    svar_.init_hydro_pressure = svar.init_hydro_pressure;
    svar_.hydro_height = svar.hydro_height ;

    svar_.jet_diam = svar.jet_diam;
    svar_.jet_depth = svar.jet_depth;
    svar_.Angle = svar.Angle;
    svar_.Rotate = svar.Rotate;
    svar_.Transp = svar.Transp;
    svar_.diam = svar.diam;
    svar_.nrad = svar.nrad;
    svar_.xyPART = svar.xyPART;
    svar_.sim_box = svar.sim_box;
    svar_.bound_box = svar.bound_box;

    /* Integration parameters */
    svar_.subits = svar.subits;
    svar_.Nframe = svar.Nframe;
    svar_.frame = svar.frame;
    svar_.bound_solver = svar.bound_solver;
    svar_.cfl = svar.cfl;
    svar_.t = svar.t;
    svar_.tframem1 = svar.tframem1;
    svar_.dt = svar.dt;
    svar_.framet = svar.framet;
    svar_.dt_max = svar.dt_max;
    svar_.dt_min = svar.dt_min;
    svar_.beta = svar.beta; 
    svar_.gamma = svar.gamma;
    svar_.maxmu = svar.maxmu;
    svar_.maxshift = svar.maxshift;
    svar_.grav = svar.grav;
    svar_.framecount = svar.framecount;
    svar_.restart_tol = svar.restart_tol;
    svar_.Force = svar.Force; 
    svar_.AForce = svar.AForce;
    svar_.mass = svar.mass;
    svar_.tMom = svar.tMom; 
    svar_.aMom = svar.aMom;

    /* SPHPart tracking settings */
    svar_.using_ipt = svar.using_ipt;
    svar_.eqOrder = svar.eqOrder;
    svar_.max_x_sph = svar.max_x_sph;
    svar_.max_x = svar.max_x;
    svar_.nSuccess = svar.nSuccess; 
    svar_.nFailed = svar.nFailed;
    svar_.IPT_diam = svar.IPT_diam; 
    svar_.IPT_area = svar.IPT_area;
    svar_.cellsout = svar.cellsout; 
    svar_.streakout = svar.streakout; 
    svar_.partout = svar.partout;

    svar_.nacross = svar.nacross;
    svar_.diameters = svar.diameters;
    svar_.velocities = svar.velocities;
    svar_.Reynolds = svar.Reynolds;
    svar_.dropDragSweep = svar.dropDragSweep;
    svar_.speedTest = svar.speedTest;
    svar_.nRuns = svar.nRuns;
    svar_.boundFile = svar.boundFile;
    svar_.fluidFile = svar.fluidFile;
    svar_.ghostFile = svar.ghostFile;
    svar_.surfaceHandle = svar.surfaceHandle;
    svar_.cellHandle = svar.cellHandle;
    svar_.streakHandle = svar.streakHandle;
    svar_.partHandle = svar.partHandle;

}

struct sph_group
{
    sph_group() : pn(SPHState(1,SPHPart(StateVecD::Zero(),StateVecD::Zero(),0.0,0.0,0.0,0,0))), 
            pnp1(SPHState(1,SPHPart(StateVecD::Zero(),StateVecD::Zero(),0.0,0.0,0.0,0,0))),
            NP1(SIMDIM,pnp1,20){};
    SIM svar;
    SPHState pn, pnp1;
    Sim_Tree NP1;
    OUTL outlist;
};


void Speed_Test_Sweep(SIM& svar, FLUID& fvar, AERO& avar)
{
    cout << "Speed testing sweep activated. Beginning sweep..." << endl;

    if(svar.Scase != SPHERE || svar.Bcase != NONE)
	{
        cout << "Only able to do test with a droplet at this stage..." << endl;
        exit(-1);
    }

    ///****** Initialise the particles memory *********/
    SPHState pn;	    /*Particles at n   */
    SPHState pnp1; 	/*Particles at n+1 */
    OUTL outlist;
    LIMITS limits;
    vector<IPTState> iptdata; /* Particle tracking data */

    svar.t = 0.0;				/*Total simulation time*/
    svar.simPts = 0;
    svar.totPts = 0;

    // InitSPH(svar,fvar,avar,pn,pnp1);
    // Only use a single particle for now
    
    pn.emplace_back(SPHPart(svar.sim_start,StateVecD::Zero(),fvar.rho0,fvar.simM,0.0,FREE,0));
    pnp1 = pn;
    limits.emplace_back(0,1);

    string file = svar.output_prefix + "_boundary.szplt.sz*";
    string cmd = "exec rm -f " + file;
    if(system(cmd.c_str()))
    {
        cout << "No prexisting boundary files deleted." << endl;
    }

    file = svar.output_prefix + "_fluid.szplt.sz*";
    cmd = "exec rm -f " + file;
    if(system(cmd.c_str()))
    {
        cout << "No prexisting fluid files deleted." << endl;
    }

    // #ifdef DEBUG
    //     ///*************** Open simulation files ***************/
    // std::fstream f1,f2,f3,fb,fg;
    // string framef = svar.output_prefix;
    // framef.append("_frame.info");
    // if(svar.restart == 1)
    //     f2.open(framef, std::ios::out | std::ios::app);
    // else
    //     f2.open(framef, std::ios::out);

    // f1 << std::scientific << std::setprecision(6);
    // f3 << std::scientific << std::setw(10);

    // // pertLog << std::scientific << std::setprecision(8);
    // Write_First_Step(f1,fb,fg,svar,pnp1);
    // #endif

    cout << "Starting counts: " << endl;
    cout << "Boundary: " << svar.bndPts << "  Sim: " << svar.simPts << endl << endl;
    cout << "dx: " << svar.dx << endl;
    cout << "Mass: " << pnp1[0].m << endl;

    SPHState pnp1_orig = pnp1;
    double dt_orig = svar.dt;


    // double error = 0.0;
    /* Data to write */
    vector<double> avg_times_sph(svar.nacross.size());
    vector<double> avg_times_ipt(svar.nacross.size());

    size_t outer_iter = 0;
    
    vector<sph_group> sph_vect(svar.nRuns);
    vector<IPTPart> IPT_nm1, IPT_n, IPT_np1;

    for(int run = 0; run < svar.nRuns; ++run)
        IPT_nm1.emplace_back(IPTPart(pnp1[0],svar.t,svar.IPT_diam,svar.IPT_area));
    

    for(int const res: svar.nacross)
    {
        /* Ingest the mesh */
        MESH cells;
        SURFS surf_marks;

        Make_Mesh(svar,avar,cells,res);
        
        cout << "Mesh size: " << cells.nElem << endl;

        ///********* Tree algorithm stuff ************/
        // KDTREE TREE(pnp1,cells);
        Vec_Tree CELL_TREE(SIMDIM,cells.cCentre,20);
        CELL_TREE.index->buildIndex();

        // Transform to do nRuns timed, and average, so that machine precision is less involved 
        // How to do efficiently for SPH however...        

        svar.t = 0.0;
        svar.dt = dt_orig;
        svar.frame = 0;
        svar.tframem1 = 0.0;
        svar.nSuccess = 0;
        svar.nFailed = 0;
        svar.totPts = pnp1.size();
        svar.simPts = pnp1.size();
        
        /* Need to find the cell it's in for IPT */
        pnp1 = pnp1_orig;
        uint to_del = 0;
        FirstCell(svar,CELL_TREE, cells, pnp1[0], to_del);

        for(int run = 0; run < svar.nRuns; ++run)
            IPT_nm1[run].reset(pnp1[0],0.0);

        IPT_n = IPT_nm1;
        IPT_np1 = IPT_nm1;
        
        // Perform the implicit particle tracking as well
        auto time1_ipt = high_resolution_clock::now();
        for(size_t ii = 0; ii < IPT_np1.size(); ++ii)
        {
            IPT::Integrate(svar, fvar,avar, cells,ii,IPT_nm1[ii],IPT_n[ii],IPT_np1[ii]
                        ,surf_marks,iptdata);
        }
        auto time2_ipt = high_resolution_clock::now();
        double time_ipt = duration_cast<microseconds>(time2_ipt-time1_ipt).count();
        
        cout <<  "\nIPT total time: " << time_ipt  << "  nSuccess: " << svar.nSuccess <<  "  nFailed: " << svar.nFailed << endl;

        svar.t = 0.0;
        svar.dt = dt_orig;
        svar.frame = 0;
        svar.tframem1 = 0.0;
        svar.nSuccess = 0;
        svar.nFailed = 0;
        svar.totPts = pnp1.size();
        svar.simPts = pnp1.size();

        for(int run = 0; run < svar.nRuns; ++run)
        {
            copy_svar(svar,sph_vect[run].svar);
            sph_vect[run].pn = pnp1_orig;
            sph_vect[run].pnp1 = pnp1_orig;
            sph_vect[run].NP1.index->buildIndex();
            FindNeighbours(sph_vect[run].NP1,fvar,sph_vect[run].pnp1,sph_vect[run].outlist);
        }

        auto time1_sph = high_resolution_clock::now();
        for(int run = 0; run < svar.nRuns; ++run)
        {
            First_Step(sph_vect[run].NP1,CELL_TREE,sph_vect[run].svar,fvar,avar,cells,limits,sph_vect[run].outlist,
                        sph_vect[run].pnp1,sph_vect[run].pn,iptdata);

            for (uint frame = 0; frame < sph_vect[run].svar.Nframe; ++frame)
            {
                //int stepits=0;
                real stept=0.0;
                
                while (stept + 0.1*sph_vect[run].svar.dt_min < sph_vect[run].svar.framet)
                {
                    Integrate(sph_vect[run].NP1,CELL_TREE,sph_vect[run].svar,fvar,avar,cells,surf_marks,
                        limits,sph_vect[run].outlist,sph_vect[run].pn,sph_vect[run].pnp1,iptdata);
                    stept+=sph_vect[run].svar.dt;
                    //++stepits;

                    if(sph_vect[run].svar.totPts == 0 || sph_vect[run].svar.simPts == 0)
                        break;
                }
                ++sph_vect[run].svar.frame;

                // #ifdef DEBUG
                // Write_Timestep(f1,fb,fg,sph_vect[run].svar,sph_vect[run].pnp1);
                // #endif
                // t2= high_resolution_clock::now();
                // duration = duration_cast<microseconds>(t2-t1).count()/1e6;

                sph_vect[run].svar.tframem1 += sph_vect[run].svar.framet; /* March frame time forward */
                if(sph_vect[run].svar.totPts == 0 || sph_vect[run].svar.simPts == 0)
                    break;
                
            }
        }
        auto time2_sph = high_resolution_clock::now();
        double time_sph = duration_cast<microseconds>(time2_sph-time1_sph).count();

        // cout << "nFailed: " << svar.nFailed << " nSuccess: " << svar.nSuccess << endl;

        /* Do some statistics of the timings */
        // double avg_time_sph = 0.0;
        // double avg_time_ipt = 0.0;
        // for(int run = 0; run < svar.nRuns; ++run)
        // {
        //     avg_time_sph += times_sph[run];
        //     avg_time_ipt += times_ipt[run];
        // }
        // avg_time_sph /= double(svar.nRuns);
        // avg_time_ipt /= double(svar.nRuns);
        // double stddev_time_sph = 0.0;
        // double stddev_time_ipt = 0.0;
        // for(int run = 0; run < svar.nRuns; ++run)
        // {
        //     stddev_time_sph += pow(times_sph[run] - avg_time_sph,2);
        //     stddev_time_ipt += pow(times_ipt[run] - avg_time_ipt,2);
        // }
        // stddev_time_sph = sqrt(stddev_time_sph/double(svar.nRuns));
        // stddev_time_ipt = sqrt(stddev_time_ipt/double(svar.nRuns));

        cout << "SPH total time: " << time_sph  << endl;
       
        avg_times_sph[outer_iter] = time_sph;
        // std_dev_times_sph[outer_iter] = stddev_time_sph;
        avg_times_ipt[outer_iter] = time_ipt;
        // std_dev_times_ipt[outer_iter] = stddev_time_ipt;
        outer_iter++;
    }
    

    /* Write data to file */
    string fileout = "Speed_Sweep_" + std::to_string(SIMDIM) + "D.dat"; 
    std::ofstream fout(fileout,std::ios::out);
    fout << "TITLE=\"SPH Speed Test\"" << endl;
    fout << "VARIABLES = \"Mesh resolution\", \"Average time (us)\"" << endl;

    fout << "ZONE T=\"SPH\"" << endl;

    for(size_t ii = 0; ii < svar.nacross.size(); ++ii)
    {
        fout << svar.nacross[ii]-1 << " " << avg_times_sph[ii] << endl;
    }
    fout << endl;

    fout << "ZONE T=\"IPT\"" << endl;
    for(size_t ii = 0; ii < svar.nacross.size(); ++ii)
    {
        fout << svar.nacross[ii]-1 << " " << avg_times_ipt[ii] << endl;
    }
    fout << endl;
    fout.close();

}
