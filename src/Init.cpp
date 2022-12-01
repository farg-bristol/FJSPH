/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "Init.h"
#include "IOFunctions.h"
#include "Kernel.h"
#include "Add.h"

#include "shapes/arc.h"
#include "shapes/circle.h"
#include "shapes/inlet.h"
#include "shapes/line.h"
#include "shapes/square.h"
#include "shapes/cylinder.h"

void get_boundary_velocity(shape_block& boundvar)
{
    size_t ntimes = boundvar.ntimes;
    if(ntimes != 0)
    {
        boundvar.vels.resize(ntimes-1);
        for(size_t jj = 0; jj < ntimes - 1; jj++)
        {
            boundvar.vels[jj] = 
                (boundvar.pos[jj+1] - boundvar.pos[jj])/
                (boundvar.times[jj+1] - boundvar.times[jj]);
        }
    }
    
}


void Read_Shapes(Shapes& var, real& globalspacing, SIM const& svar, FLUID const& fvar, std::string const& filename)
{
    /* Read shapes from file */
    std::ifstream inFile(filename);
    if (!inFile.is_open())
    {
        std::cerr << filename << " file missing\n";
        exit(1);
    }

    std::string line;
    size_t nblocks = 0;

    /* Find number of blocks to initialise */
    while(getline(inFile,line))
    {
        size_t end = line.find_first_of('#');
		if(end != std::string::npos)
			line = line.substr(0,end+1);

        if(line.find("block end")!=std::string::npos)
        {
            nblocks++;
        }
    }

	inFile.clear();
	inFile.seekg(0);
	
    std::vector<shape_block> shapes(nblocks);

    size_t block = 0;
    while(getline(inFile,line))
    {
        size_t end = line.find_first_of('#');
		if(end != std::string::npos)
			line = line.substr(0,end);

        Get_String(line, "Name",shapes[block].name);
        Get_String(line, "Shape", shapes[block].shape);
        Get_String(line, "Sub-shape", shapes[block].subshape);
		Get_String(line, "Boundary solver", shapes[block].solver_name);

		Get_Number(line, "Write surface data (0/1)", shapes[block].write_data);
		Get_Number(line, "Fixed velocity or dynamic inlet BC (0/1)",shapes[block].fixed_vel_or_dynamic);

		// Pipe exit plane.
		Get_Vector(line, "Aerodynamic entry normal", shapes[block].pipe_norm);
		Get_Vector(line, "Deletion normal", shapes[block].delete_norm);
		Get_Vector(line, "Insertion normal", shapes[block].insert_norm);
		Get_Number(line, "Aerodynamic entry plane constant", shapes[block].pipeconst);
		Get_Number(line, "Deletion plane constant", shapes[block].delconst);
		Get_Number(line, "Insertion plane constant", shapes[block].insconst);
        Get_Number(line, "Pipe depth", shapes[block].thickness);

        Get_Number(line, "i-direction count", shapes[block].ni);
		Get_Number(line, "j-direction count", shapes[block].nj);
		Get_Number(line, "k-direction count", shapes[block].nk);
        
       	Get_Vector(line, "Stretching factor", shapes[block].stretch);

		// Rotation definitions
		Get_Vector(line, "Normal vector", shapes[block].normal);
		Get_Vector(line, "Rotation angles", shapes[block].angles);
		Get_Number(line, "Rotation angle", shapes[block].angles[0]);

        // Square block definitions
		Get_Vector(line, "Start coordinate", shapes[block].start);
		Get_Vector(line, "End coordinate", shapes[block].end);
        Get_Vector(line, "Right coordinate", shapes[block].right);

        // Circle/arc definitions
        Get_Vector(line, "Midpoint coordinate", shapes[block].mid);
		Get_Vector(line, "Centre coordinate", shapes[block].centre);
        Get_Vector(line, "Arch normal", shapes[block].right);
        Get_Number(line, "Radius", shapes[block].radius);
        Get_Number(line, "Length", shapes[block].length);
		Get_Number(line, "Arc start (degree)", shapes[block].arc_start);
		Get_Number(line, "Arc end (degree)", shapes[block].arc_end);
		Get_Number(line, "Arc length (degree)", shapes[block].arclength);
        Get_Number(line, "Starting straight length", shapes[block].sstraight);
        Get_Number(line, "Ending straight length", shapes[block].estraight);

        Get_Number(line, "Particle spacing", shapes[block].dx);
		Get_Number(line, "Particle ordering (0=grid,1=HCP)", shapes[block].hcpl);
        Get_Number(line, "Wall thickness", shapes[block].thickness);
        Get_Number(line, "Wall radial particle count", shapes[block].nk);
        Get_Number(line, "Wall is no slip (0/1)", shapes[block].no_slip);

        Get_String(line, "Coordinate filename", shapes[block].filename);
        
        if(line.find("Coordinate data:") != std::string::npos)
        {
            std::string tmp;
            getline(inFile,tmp); /* Go to the next line and get data */
            size_t npts;
            inFile >> npts;
            shapes[block].npts = npts;
            shapes[block].coords = std::vector<StateVecD>(npts);

            for(size_t ii = 0; ii < shapes[block].npts; ++ii)
            {
                for(size_t dim = 0 ; dim < SIMDIM; ++dim)
                {
                    inFile >> shapes[block].coords[ii][dim];
                }
            }
        }

        Get_Vector(line, "Starting velocity", shapes[block].vel);
		Get_Number(line, "Starting pressure", shapes[block].press);
        Get_Number(line, "Starting density", shapes[block].dens);
        // Get_Number(line, "Starting mass", shapes[block].mass);

        Get_Number(line, "Cole EOS gamma", shapes[block].gamma);
		Get_Number(line, "Speed of sound", shapes[block].speedOfSound);
        Get_Number(line, "Resting density", shapes[block].rho0);
        Get_Number(line, "Volume to target", shapes[block].renorm_vol);

        Get_String(line, "Time data filename", shapes[block].position_filename);

        if(line.find("Time position data") != std::string::npos)
        {
            std::string tmp;
            getline(inFile,tmp); /* Go to the next line and get data */
            size_t ntimes;
            std::istringstream iss(tmp);
            iss >> ntimes;
            // inFile >> ntimes;
            shapes[block].ntimes = ntimes;
            shapes[block].times = std::vector<real>(ntimes);
            shapes[block].pos = std::vector<StateVecD>(ntimes);

            for(size_t ii = 0; ii < ntimes; ++ii)
            {
                getline(inFile,tmp); /* Go to the next line and get data */
                std::istringstream iss2(tmp);
                iss2 >> shapes[block].times[ii];
                for(size_t dim = 0 ; dim < SIMDIM; ++dim)
                {
                    iss2 >> shapes[block].pos[ii][dim];
                }
            }
        }

        if(line.find("Time velocity data") != std::string::npos)
        {
            std::string tmp;
            getline(inFile,tmp); /* Go to the next line and get data */
            size_t ntimes;
            std::istringstream iss(tmp);
            iss >> ntimes;
            // inFile >> ntimes;
            shapes[block].ntimes = ntimes;
            shapes[block].times = std::vector<real>(ntimes);
            shapes[block].vels = std::vector<StateVecD>(ntimes);

            for(size_t ii = 0; ii < ntimes; ++ii)
            {
                getline(inFile,tmp); /* Go to the next line and get data */
                std::istringstream iss2(tmp);
                iss2 >> shapes[block].times[ii];
                for(size_t dim = 0 ; dim < SIMDIM; ++dim)
                {
                    iss2 >> shapes[block].vels[ii][dim];
                }
            }
        }

        if (line.find("block end") != std::string::npos)
        {
            block++;
            if (block == nblocks)
                break;
        }
    }

    inFile.close();

    /* Check enough information has been provided */
    size_t count = 0;
    for(shape_block& bound:shapes)
    {
        #if SIMDIM == 2
        if (bound.shape == "Line")
        #else
        if (bound.shape == "Plane")
        #endif
        {
            bound.bound_type = linePlane;
            check_line_input(bound, globalspacing);
        }
        #if SIMDIM == 2
        if (bound.shape == "Square")
        #else
        if (bound.shape == "Cube")
        #endif
        {
            bound.bound_type = squareCube;
            check_square_input(bound,globalspacing);
        }
        #if SIMDIM == 2
        else if (bound.shape == "Circle")
        #else
        else if (bound.shape == "Sphere")
        #endif
        {
            bound.bound_type = circleSphere;
            check_circle_input(bound,globalspacing);
        }
        #if SIMDIM == 2
        else if (bound.shape == "Arc")
        #else 
        else if (bound.shape == "Arch")
        #endif
        {
            bound.bound_type = arcSection;
            check_arc_input(bound,globalspacing);
        }
        else if (bound.shape == "Cylinder")
        {
            bound.bound_type = cylinder;
            check_cylinder_input(bound,globalspacing);
        }
		else if (bound.shape == "Inlet")
		{
			bound.bound_type = inletZone;
			check_inlet_input(bound,globalspacing);
		}
        else if (bound.shape == "Coordinates")
        {
            bound.bound_type = coordDef;
            if(bound.filename.empty())
            {
                if(bound.npts == 0 || bound.coords.empty())
                {
                    std::cout << "ERROR: Block \"" << bound.name << "\" coordinates have not been ingested properly. Stopping." << std::endl;
                    exit(1);
                }
            }
        }

        if(bound.position_filename.empty())
        {
            if(bound.ntimes != 0) 
            {
                if(bound.times.empty())
                {
                    cout << "ERROR: Block \"" << bound.name << "\" has no time information" << endl;
                    exit(-1);
                } 
                if(bound.pos.empty() && bound.vels.empty())
                {
                    std::cout << "ERROR: Block \"" << bound.name << "\" position or velocity data has not been ingested properly. Stopping." << std::endl;
                    exit(1);
                }
            }
        }

        if(bound.bound_type == -1)
        {
            std::cout << "ERROR: Block \"" << bound.name << "\" shape has not been correctly defined. Stopping." << std::endl;
            std::cout << "File: " << filename << std::endl;
            std::cout << "Block: " << bound.name << " \tID: " << count << std::endl;

            exit(1);
        }

		if(bound.solver_name == "DBC")
		{
			bound.bound_solver = DBC;
		}
		else if (bound.solver_name == "Pressure-Gradient")
		{
			bound.bound_solver = pressure_G;
		}
        else if (bound.solver_name == "Ghost")
        {
            bound.bound_solver = ghost;
        }

        globalspacing = std::max(bound.dx,globalspacing);
            bound.npts = bound.npts > 1 ? bound.npts : 1;
        if(bound.bound_type != coordDef)
            bound.coords.reserve(bound.npts);

        if(bound.press != 0)
        {
            real Bconst = fvar.rho0 * (fvar.Cs * fvar.Cs) / fvar.gam;
            bound.dens = pow(((bound.press - fvar.pPress) / Bconst + 1.0), fvar.gam) * fvar.rho0;
        }
        else
        {
            bound.dens = fvar.rho0;
        }

        if(bound.nu < 0)
        {
            bound.nu = fvar.nu;
        }

        ++count;
    
    }/* End for bound:shapes */

    var.nblocks = nblocks;
    var.block = shapes;

    // Estimate the number of particles in each Block
    for (size_t ii = 0; ii < var.nblocks; ++ii)
    {
        var.totPts += var.block[ii].npts;
    }

    return;
}

std::vector<StateVecD> Read_Geom_File(std::string const& filename)
{
    std::ifstream inFile(filename,std::ios::in);

    int nPts;
    inFile >> nPts;
    std::vector<StateVecD> points(nPts);

    for(int ii = 0; ii < nPts; ++ii)
    {
        for(int dim = 0; dim < SIMDIM; ++dim)
        {
            inFile >> points[ii][dim];
        }
    }

    inFile.close();

    return points;
}

void Check_Intersection(SIM const& svar, Shapes& boundvar, Shapes& fluvar)
{
    // Search for overlap with particles
    double searchDist = 0.9 * svar.dx;
    searchDist = searchDist*searchDist;

	// Check boundary blocks first
    for(size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        boundvar.block[blockID].intersect.assign(boundvar.block[blockID].npts, 0);
    }

    for(size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
        nanoflann::SearchParams const params(0,0,false);
        Vec_Tree tree(SIMDIM,boundvar.block[blockID].coords,50);
        tree.index->buildIndex();

        // Search for overlap with boundary particles (do the upper diagonal only)
        for (size_t ii = blockID; ii < boundvar.nblocks; ++ii)
        {
            /*if (ii == blockID)
                continue;*/

            for(size_t jj = 0; jj < boundvar.block[ii].coords.size(); ++jj)
            {
                if(boundvar.block[ii].intersect[jj] == 0)
                {
                    // size_t no(Searchtree(tree, fluvar.block[blockID].coords, particle, searchDist, outlist_local));
                    std::vector<std::pair<size_t, double>> matches; /* Nearest Neighbour Search*/
                    #if SIMDIM == 3
                        matches.reserve(250);
                    #else
                        matches.reserve(47);
                    #endif

                    tree.index->radiusSearch(&boundvar.block[ii].coords[jj][0], searchDist, matches, params);

                    // std::cout << ll << "  " << blockID << "  " << matches.size() << std::endl; 
                    for (auto const& match : matches)
                    {
                        if (match.first != jj)
                            boundvar.block[blockID].intersect[match.first] = 1;
                    }
                }
            }
        }
    }
	

	// Now check fluid blocks
    for(size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        fluvar.block[blockID].intersect.assign(fluvar.block[blockID].npts, 0);
    }
    
    for(size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
        nanoflann::SearchParams const params(0,0,false);
        Vec_Tree tree(SIMDIM,fluvar.block[blockID].coords,10);
        tree.index->buildIndex();

        // Search for overlap with boundary particles
        for (size_t ii = 0; ii < boundvar.nblocks; ++ii)
        {
            for(size_t jj = 0; jj < boundvar.block[ii].coords.size(); ++jj)
            {
                if(boundvar.block[ii].intersect[jj] == 0)
                {
                    std::vector<std::pair<size_t, double>> matches; /* Nearest Neighbour Search*/
                    #if SIMDIM == 3
                        matches.reserve(250);
                    #else
                        matches.reserve(47);
                    #endif

                    tree.index->radiusSearch(&boundvar.block[ii].coords[jj][0], searchDist, matches, params);

                    // std::cout << ll << "  " << blockID << "  " << matches.size() << std::endl; 
                    for (auto const& match : matches)
                        fluvar.block[blockID].intersect[match.first] = 1;
                }
            }
        }

        // Check for intersection with other fluid blocks (do the upper diagonal only)
        for (size_t ii = blockID; ii < fluvar.nblocks; ++ii)
        {
            if(ii == blockID)
                continue; /* Ignore points within its own block. Should not self intersect */

            for(size_t jj = 0; jj < fluvar.block[ii].npts; ++jj)
            {
                if(fluvar.block[ii].intersect[jj] == 0)
                {
                    // size_t no(Searchtree(tree, fluvar.block[blockID].coords, 
                    //     fluvar.block[blockID].coords[ll], searchDist, outlist_local));
                    std::vector<std::pair<size_t, double>> matches; /* Nearest Neighbour Search*/
                    #if SIMDIM == 3
                        matches.reserve(50);
                    #else
                        matches.reserve(10);
                    #endif

                    tree.index->radiusSearch(&fluvar.block[ii].coords[jj][0], searchDist, matches, params);

                    // std::cout << ii << "  " << blockID << "  " << matches.size() << std::endl; 
                    for (size_t kk = 0; kk < matches.size(); kk++)
                    {          
                        fluvar.block[blockID].intersect[matches[kk].first] = 1;   
                    }
                }                
            }
        }
    }

    for(size_t ii = 0; ii < fluvar.nblocks; ++ii)
    {
        for(size_t bID = 0; bID < fluvar.block[ii].back.size(); bID++)
        {
            size_t backID = fluvar.block[ii].back[bID];
            int does_int = 0;
            if(fluvar.block[ii].intersect[backID])
            {   // If the back particle interstects, don't place any of the inlet
                does_int = 1;
            }
            else
            {   // Check buffer of that back particle now
                for(size_t const buffID : fluvar.block[ii].buffer[bID])
                {
                    if(fluvar.block[ii].intersect[buffID])
                        does_int = 1;
                }
            }

            if(does_int)
            {   // Set the intersection of all buffer and back particles to have correct numbers
                fluvar.block[ii].intersect[backID] = 1;

                for(size_t const buffID : fluvar.block[ii].buffer[bID])
                    fluvar.block[ii].intersect[buffID] = 1;
            }
        }
    }

    size_t nBoundPts = 0, nFluidPts = 0;
    for(size_t blockID = 0; blockID < boundvar.nblocks; blockID++)
    {
         // Update the points to be accurate to that retained.
        boundvar.block[blockID].npts = std::count(boundvar.block[blockID].intersect.begin(), boundvar.block[blockID].intersect.end(), 0);

        nBoundPts += boundvar.block[blockID].npts;
    }

    for(size_t blockID = 0; blockID < fluvar.nblocks; blockID++)
    {
         // Update the points to be accurate to that retained.
        fluvar.block[blockID].npts = std::count(fluvar.block[blockID].intersect.begin(), fluvar.block[blockID].intersect.end(), 0);

        nFluidPts += fluvar.block[blockID].npts;
    }

    // Get number of fluid particles to keep
    printf("Keeping %zu out of %zu boundary points. %.2f%% intersected.\n",
            nBoundPts,boundvar.totPts,100.0*(boundvar.totPts-nBoundPts)/boundvar.totPts);
            
    printf("Keeping %zu out of %zu fluid points. %.2f%% intersected.\n",
            nFluidPts,fluvar.totPts,100.0*(fluvar.totPts-nFluidPts)/fluvar.totPts);

    printf("Keeping %zu out of %zu total particles.\n",nBoundPts + nFluidPts,fluvar.totPts + boundvar.totPts );
	boundvar.totPts = nBoundPts;
	fluvar.totPts = nFluidPts;
}

size_t Generate_Points(SIM const& svar, FLUID const& fvar, double const& globalspacing, Shapes& var)
{
    size_t totPts = 0;
    size_t diff = 0;
    for(shape_block& bound:var.block)
    {
        std::cout << "Creating boundary block: " << bound.name << "\t...\t";
        if(bound.bound_type == linePlane)
        {
            #if SIMDIM == 2
            bound.coords = create_line(bound, globalspacing);
            #else
            bound.coords = create_plane(bound, globalspacing);
            #endif
        }
        else if (bound.bound_type == squareCube)
        {
            bound.coords = 
            create_square(bound.start,bound.end, globalspacing, bound.hcpl);
        }
        else if (bound.bound_type == circleSphere)
        {
            bound.coords = 
            create_circle(bound.centre,bound.radius, globalspacing, bound.hcpl);
        }
        else if (bound.bound_type == cylinder)
        {
            bound.coords = create_cylinder(bound,globalspacing);
        }
        else if (bound.bound_type == arcSection)
        {
            bound.coords = create_arc_segment(bound,globalspacing);
        }
        else if (bound.bound_type == inletZone)
        {
            bound.coords = create_inlet_zone(bound,globalspacing);
        }
        else if (bound.bound_type == coordDef)
        {
            if(!bound.filename.empty())
            {   /* Read the file */
                bound.coords = Read_Geom_File(bound.filename);
            }
        }

        if (bound.coords.size() != bound.npts)
		{
			std::cerr << std::endl << "Number of boundary points generated for boundary \"" << 
                bound.name << "\" differs from expected amount by "
			 << static_cast<int>(bound.npts) - static_cast<int>(bound.coords.size()) << std::endl;
             diff += bound.coords.size() - bound.npts;
		}

        bound.npts = bound.coords.size();
        totPts += bound.npts;
        std::cout << "npts: " << bound.npts << std::endl;

        if(svar.use_global_gas_law)
        {   // Set block gas law properties as the global properties
            bound.rho0 = fvar.rho0;
            bound.gamma = fvar.gam;
            bound.speedOfSound = fvar.Cs;
            bound.backgroundP = fvar.pPress;
        }
    }
    var.totPts = totPts;
    // var.totPts = boundary_intersection(globalspacing, var);
    
    // get_boundary_velocity(boundvar);

    return diff;
}



void Init_Particles(SIM& svar, FLUID& fvar, AERO& avar, SPHState& pn, SPHState& pnp1, LIMITS& limits)
{
	real dx = svar.dx;
	// Read boundary blocks
	Shapes boundvar;
	cout << "Reading boundary settings..." << endl;
	Read_Shapes(boundvar,dx,svar,fvar,svar.boundfile);

	// Read fluid
	Shapes fluvar;
	cout << "Reading fluid settings..." << endl;
	Read_Shapes(fluvar,dx,svar,fvar,svar.fluidfile);

	// Now generate points and add to indexes
	Generate_Points(svar, fvar, svar.dx, boundvar);
	Generate_Points(svar, fvar, svar.dx, fluvar);

	// Now check for intersections
	Check_Intersection(svar,boundvar,fluvar);

	// Insert points into the state vectors. Change to SoA?
	size_t pID = 0;
	for(size_t block = 0; block < boundvar.nblocks; block++)
	{
		limits.emplace_back(pID,0);
		
		for(size_t ii = 0; ii < boundvar.block[block].coords.size(); ii++)
		{
			if(!boundvar.block[block].intersect[ii])
			{
				pn.emplace_back(SPHPart(boundvar.block[block].coords[ii],
							boundvar.block[block].vel,
							boundvar.block[block].dens,fvar.bndM,
							boundvar.block[block].press,BOUND,pID));
				pID++;
			}
		}
		if(!boundvar.block[block].times.empty())
		{
            if(!boundvar.block[block].pos.empty() && boundvar.block[block].vels.empty())
                get_boundary_velocity(boundvar.block[block]);
            else if(boundvar.block[block].vels.empty())
            {
                cout << "No velocity or position data available for boundary block " << block <<
                    " even though times were defined." << endl;
                exit(-1);
            }

			limits.back().times = boundvar.block[block].times;
			limits.back().vels = boundvar.block[block].vels;
			limits.back().nTimes = boundvar.block[block].ntimes;
		}
		else
		{
			limits.back().nTimes = 0;
			limits.back().vels.emplace_back(boundvar.block[block].vel);
		}

		limits.back().index.second = pID;

        limits.back().name = boundvar.block[block].name;
		limits.back().fixed_vel_or_dynamic = boundvar.block[block].fixed_vel_or_dynamic;
		limits.back().hcpl = boundvar.block[block].hcpl;
		limits.back().no_slip = boundvar.block[block].no_slip;
		limits.back().bound_solver = boundvar.block[block].bound_solver;
		limits.back().block_type = boundvar.block[block].bound_type;
	}

	svar.nbound = boundvar.nblocks;
	svar.bndPts = pID;

    // Fluid points
	for(size_t block = 0; block < fluvar.nblocks; block++)
	{	
		limits.emplace_back(pID,0);
        if(fluvar.block[block].bound_type == inletZone)
        {
            size_t ii = 0;
            size_t nBuff = fluvar.block[block].hcpl == 1 ? 5 : 4;

            while(fluvar.block[block].bc[ii] == PIPE)
            {
                if(!fluvar.block[block].intersect[ii])
                {
                    pn.emplace_back(SPHPart(fluvar.block[block].coords[ii],
                                fluvar.block[block].vel,
                                fluvar.block[block].dens,fvar.bndM,
                                fluvar.block[block].press,PIPE,pID));
                    pID++;
                }
                ii++;
            }

            // Check if the back line or buffer ends up intersecting
            vector<size_t> back_intersect(fluvar.block[block].back.size(),0);
            
            for(size_t bID = 0; bID < fluvar.block[block].back.size(); bID++)
            {
                size_t does_int = 0;
                size_t backID = fluvar.block[block].back[bID];

                // If the back particle interstects, don't place any of the inlet
                if(fluvar.block[block].intersect[backID])
                    does_int = 1;

                back_intersect[bID] = does_int;
            }

            // Insert the back particles
            for(size_t bID = 0; bID < fluvar.block[block].back.size(); bID++)
            {
                if(!back_intersect[bID])
                {
                    size_t pointID = fluvar.block[block].back[bID];
                    pn.emplace_back(SPHPart(fluvar.block[block].coords[pointID],
                                fluvar.block[block].vel,
                                fluvar.block[block].dens,fvar.simM,
                                fluvar.block[block].press,BACK,pID));

                    limits.back().back.emplace_back(pID);
                    limits.back().buffer.emplace_back(nBuff,0);

                    pID++;
                }
            }
            // Insert the buffer particles
            for(size_t buffID = 0; buffID < nBuff; buffID++)
            {
                size_t bIndex = 0;
                for(size_t bID = 0; bID < fluvar.block[block].back.size(); bID++)
                {
                    if(!back_intersect[bID])
                    {
                        size_t pointID = fluvar.block[block].buffer[bID][buffID];
                        pn.emplace_back(SPHPart(fluvar.block[block].coords[pointID],
                                    fluvar.block[block].vel,
                                    fluvar.block[block].dens,fvar.simM,
                                    fluvar.block[block].press,BUFFER,pID));

                        limits.back().buffer[bIndex][buffID] = pID;

                        pID++;
                        bIndex++;
                    }
                }
            }
        }
        else
        {
            for(size_t ii = 0; ii < fluvar.block[block].coords.size(); ii++)
            {
                if(!fluvar.block[block].intersect[ii])
                {
                    pn.emplace_back(SPHPart(fluvar.block[block].coords[ii],
                                fluvar.block[block].vel,
                                fluvar.block[block].dens,fvar.simM,
                                fluvar.block[block].press,FREE,pID));
                    pID++;
                }
            }
        }

		limits.back().index.second = pID;
		if(!fluvar.block[block].times.empty())
		{
			limits.back().times = fluvar.block[block].times;
			limits.back().vels = fluvar.block[block].vels;
			limits.back().nTimes = fluvar.block[block].ntimes;
		}
        else
        {
            limits.back().nTimes = 0;
            limits.back().vels.emplace_back(StateVecD::Zero());
        }

        limits.back().name = fluvar.block[block].name;
		limits.back().insert_norm = fluvar.block[block].insert_norm;
		limits.back().delete_norm = fluvar.block[block].delete_norm;
		limits.back().pipe_norm = fluvar.block[block].pipe_norm;
		limits.back().insconst = fluvar.block[block].insconst;
		limits.back().delconst = fluvar.block[block].delconst;
		limits.back().pipeconst = fluvar.block[block].pipeconst;
		
		limits.back().fixed_vel_or_dynamic = fluvar.block[block].fixed_vel_or_dynamic;
		limits.back().hcpl = fluvar.block[block].hcpl;
		limits.back().no_slip = fluvar.block[block].no_slip;
		limits.back().bound_solver = fluvar.block[block].bound_solver;
		limits.back().block_type = fluvar.block[block].bound_type;

        if(fluvar.block[block].renorm_vol != -1)
        {
            /* Renormalise the volume */
            double nfill = limits.back().index.second - limits.back().index.first;

            double mass = fluvar.block[block].dens * fluvar.block[block].renorm_vol / nfill;
            #pragma omp for nowait 
			for (size_t ii=0; ii < pID; ++ii)
            {
                pn[ii].m = mass;
            }

            // Adjust the support radius and spacing. Only needed for lid driven cavity.
            fvar.simM = mass;
            fvar.bndM = mass;
            svar.Pstep = pow(mass/fluvar.block[block].dens, 1.0/SIMDIM);
            
            fvar.H = fvar.Hfac*svar.Pstep;
            fvar.HSQ = fvar.H*fvar.H; 
            fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/

            fvar.dCont = 2.0 * fvar.delta * fvar.H * fvar.Cs;
            
            #if SIMDIM == 2
                #ifdef CUBIC
                    fvar.correc = 10.0 / (7.0 * M_PI * fvar.H * fvar.H);
                #else
                    fvar.correc = 7.0 / (4.0 * M_PI * fvar.H * fvar.H);
                #endif
            #endif
            #if SIMDIM == 3
                #ifdef CUBIC
                    fvar.correc = (1.0/(M_PI*fvar.H*fvar.H*fvar.H));
                #else
                    fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
                #endif
            #endif
                
            fvar.Wdx = Kernel(svar.Pstep,fvar.H,fvar.correc);
            
            avar.GetYcoef(fvar, /*fvar.H*/ svar.Pstep);

            #if SIMDIM == 3
                avar.aPlate = svar.Pstep * svar.Pstep;
                // avar.aPlate = fvar.H*fvar.H;
            #else
                avar.aPlate = svar.Pstep /**svar.Pstep*/ /** pow(avar.L,0.5)*/;
                // avar.aPlate = fvar.H;
            #endif	
            
            /* Particle tracking values */
            svar.IPT_diam = pow((6.0*fvar.simM)/(M_PI*fvar.rho0),1.0/3.0);
            svar.IPT_area = M_PI * svar.IPT_diam*svar.IPT_diam/4.0; 
        }
	}
	svar.nfluid = fluvar.nblocks;
	svar.totPts = pID;
	svar.simPts = svar.totPts - svar.bndPts;
    svar.partID = pID;

    
	if(svar.init_hydro_pressure)
	{	/* If using hydrostatic initialisation, set pressure */
		cout << "Initialising hydrostatic pressure..." << endl;
        #pragma omp parallel for
		for(size_t ii = 0; ii < svar.totPts; ++ii)
		{
			real press = std::max(0.0,-fvar.rho0 * svar.grav[1] * (svar.hydro_height -  pn[ii].xi[1]));
			real dens = density_equation(press, fvar.B, fvar.gam, fvar.Cs, fvar.rho0, fvar.backP);
			pn[ii].p = press;
			pn[ii].rho = dens;
		}
	}

	pnp1 = pn;
}

void Init_Surface(SIM const& svar, MESH const& cells, vector<SURF>& surf_marks)
{
	surf_marks = vector<SURF>(svar.markers.size());

	vector<vector<size_t>> faceIDs(svar.markers.size());
	vector<vector<int>> markers(svar.markers.size());
	/* for each surface, find how many faces are in it. */
	for(std::pair<size_t,int> const& marker:cells.smarkers)
	{
		auto index = find(svar.markers.begin(),svar.markers.end(),marker.second);
		if(index != svar.markers.end())
		{
			size_t mark = index - svar.markers.begin();
			faceIDs[mark].emplace_back(marker.first);
			markers[mark].emplace_back(marker.second);

			// surf_faces[mark].back().faceID  = marker.first;
			// if()
			// cout << mark << "  " << markers.first << endl;
			// cout << surf_faces[mark].back().faceID << "  " << markers.first << endl;
			// surf_faces[mark].back().marker  = marker.second;
		}
		else
		{
			cout << "Couldn't find the marker in the index" << endl;
		}
	}

	for(size_t ii = 0; ii < svar.markers.size(); ii++)
	{
		surf_marks[ii].name = svar.bnames[ii];
		surf_marks[ii].marker = svar.markers[ii];
		surf_marks[ii].output = svar.bwrite[ii];
		
		size_t nFaces = faceIDs[ii].size();
		surf_marks[ii].faceIDs = faceIDs[ii];
		surf_marks[ii].face_count = vector<uint>(nFaces,0);
		surf_marks[ii].face_beta = vector<real>(nFaces,0.0);
		surf_marks[ii].face_area = vector<real>(nFaces,0.0);

	}
}
