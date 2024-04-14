/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#include "Add.h"
#include "IOFunctions.h"
#include "Init.h"
#include "Kernel.h"
#include "Neighbours.h"

#include "shapes/arc.h"
#include "shapes/circle.h"
#include "shapes/cylinder.h"
#include "shapes/inlet.h"
#include "shapes/line.h"
#include "shapes/square.h"

class fluid_block
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    fluid_block() {}

    /* Read the input for a block from the input file. */
    void read_block()
    {
        while (getline(inFile, line))
        {
            size_t end = line.find_first_of('#');
            if (end != std::string::npos)
                line = line.substr(0, end);

            Get_String(line, "Name", self.name);
            Get_String(line, "Shape", self.shape);
            Get_String(line, "Sub-shape", self.subshape);
            Get_String(line, "Boundary solver", self.solver_name);

            Get_Number(line, "Write surface data (0/1)", self.write_data);
            Get_Number(line, "Fixed velocity or dynamic inlet BC (0/1)", self.fixed_vel_or_dynamic);

            // Pipe exit plane.
            Get_Vector(line, "Aerodynamic entry normal", self.aero_norm);
            Get_Vector(line, "Deletion normal", self.delete_norm);
            Get_Vector(line, "Insertion normal", self.insert_norm);
            Get_Number(line, "Aerodynamic entry plane constant", self.aeroconst);
            Get_Number(line, "Deletion plane constant", self.delconst);
            Get_Number(line, "Insertion plane constant", self.insconst);
            Get_Number(line, "Pipe depth", self.thickness);

            Get_Number(line, "i-direction count", self.ni);
            Get_Number(line, "j-direction count", self.nj);
            Get_Number(line, "k-direction count", self.nk);

            Get_Vector(line, "Stretching factor", self.stretch);

            // Rotation definitions
            Get_Vector(line, "Normal vector", self.normal);
            Get_Vector(line, "Rotation angles", self.angles);
            Get_Number(line, "Rotation angle", self.angles[0]);

            // Square block definitions
            Get_Vector(line, "Start coordinate", self.start);
            Get_Vector(line, "End coordinate", self.end);
            Get_Vector(line, "Right coordinate", self.right);

            // Circle/arc definitions
            Get_Vector(line, "Midpoint coordinate", self.mid);
            Get_Vector(line, "Centre coordinate", self.centre);
            Get_Vector(line, "Arch normal", self.right);
            Get_Number(line, "Radius", self.radius);
            Get_Number(line, "Length", self.length);
            Get_Number(line, "Arc start (degree)", self.arc_start);
            Get_Number(line, "Arc end (degree)", self.arc_end);
            Get_Number(line, "Arc length (degree)", self.arclength);
            Get_Number(line, "Starting straight length", self.sstraight);
            Get_Number(line, "Ending straight length", self.estraight);

            Get_Number(line, "Particle spacing", self.dx);
            Get_Number(line, "Particle ordering (0=grid,1=HCP)", self.particle_order);
            Get_Number(line, "Wall thickness", self.thickness);
            Get_Number(line, "Wall radial particle count", self.nk);
            Get_Number(line, "Wall is no slip (0/1)", self.no_slip);

            Get_String(line, "Coordinate filename", self.filename);

            if (line.find("Coordinate data:") != std::string::npos)
            {
                std::string tmp;
                getline(inFile, tmp); /* Go to the next line and get data */
                size_t npts;
                inFile >> npts;
                self.npts = npts;
                self.coords = std::vector<StateVecD>(npts);

                for (size_t ii = 0; ii < self.npts; ++ii)
                {
                    for (size_t dim = 0; dim < SIMDIM; ++dim)
                    {
                        inFile >> self.coords[ii][dim];
                    }
                }
            }

            Get_Vector(line, "Starting velocity", self.vel);
            Get_Number(line, "Starting jet velocity", self.vmag);
            Get_Number(line, "Starting pressure", self.press);
            Get_Number(line, "Starting density", self.dens);
            // Get_Number(line, "Starting mass", self.mass);

            Get_Number(line, "Cole EOS gamma", self.gamma);
            Get_Number(line, "Speed of sound", self.speedOfSound);
            Get_Number(line, "Resting density", self.rho0);
            Get_Number(line, "Volume to target", self.renorm_vol);

            Get_String(line, "Time data filename", self.position_filename);

            if (line.find("Time position data") != std::string::npos)
            {
                std::string tmp;
                getline(inFile, tmp); /* Go to the next line and get data */
                size_t ntimes;
                std::istringstream iss(tmp);
                iss >> ntimes;
                // inFile >> ntimes;
                self.ntimes = ntimes;
                self.times = std::vector<real>(ntimes);
                self.pos = std::vector<StateVecD>(ntimes);

                for (size_t ii = 0; ii < ntimes; ++ii)
                {
                    getline(inFile, tmp); /* Go to the next line and get data */
                    std::istringstream iss2(tmp);
                    iss2 >> self.times[ii];
                    for (size_t dim = 0; dim < SIMDIM; ++dim)
                    {
                        iss2 >> self.pos[ii][dim];
                    }
                }
            }

            if (line.find("Time velocity data") != std::string::npos)
            {
                std::string tmp;
                getline(inFile, tmp); /* Go to the next line and get data */
                size_t ntimes;
                std::istringstream iss(tmp);
                iss >> ntimes;
                // inFile >> ntimes;
                self.ntimes = ntimes;
                self.times = std::vector<real>(ntimes);
                self.vels = std::vector<StateVecD>(ntimes);

                for (size_t ii = 0; ii < ntimes; ++ii)
                {
                    getline(inFile, tmp); /* Go to the next line and get data */
                    std::istringstream iss2(tmp);
                    iss2 >> self.times[ii];
                    for (size_t dim = 0; dim < SIMDIM; ++dim)
                    {
                        iss2 >> self.vels[ii][dim];
                    }
                }
            }

            if (line.find("block end") != std::string::npos)
            {
                // End if we enter a new block.
                break;
            }
        }
    }
}