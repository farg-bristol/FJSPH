/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "coordinates.h"

// Public functions

void CoordShape::check_input(SIM const& svar, int& fault)
{
    ShapeBlock::check_input(svar, fault);

    if (filename.empty())
    {
        if (coords.empty())
        {
            printf(
                "ERROR: Block \"%s\" coordinates have not been ingested "
                "properly. Stopping.\n",
                name.c_str()
            );
            fault = 1;
        }
        else
        {
            npts = coords.size();
        }
    }

    ShapeBlock::check_input_post();
}

void CoordShape::generate_points()
{
    // Nothing to do for coordinates. Should have already been read from the file.
}
