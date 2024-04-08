/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "shapes.h"

#include "arc.h"
#include "circle.h"
#include "cylinder.h"
#include "inlet.h"
#include "line.h"
#include "square.h"

void Read_Shapes(
    Shapes& var, real& globalspacing, SIM const& svar, FLUID const& fvar, std::string const& filename
)
{
}

void shape_block::check_input(real& globalspacing, int& fault)
{
    switch (bound_type)
    {
    case linePlane:
        LineShape::check_input(*this, globalspacing, fault);
        break;
    case squareCube:
        SquareShape::check_input(*this, globalspacing, fault);
        break;
    case circleSphere:
        CircleShape::check_input(*this, globalspacing, fault);
        break;
    case arcSection:
        ArcShape::check_input(*this, globalspacing, fault);
        break;
    case cylinder:
        CircleShape::check_input(*this, globalspacing, fault);
        break;
    case inletZone:
        InletShape::check_input(*this, globalspacing, fault);
        break;
    case coordDef:
        if (filename.empty())
        {
            if (npts == 0 || coords.empty())
            {
                printf(
                    "ERROR: Block \"%s\" coordinates have not been ingested "
                    "properly. Stopping.\n",
                    name.c_str()
                );
                fault = 1;
            }
        }
        break;
    default:
        printf("Unrecognised boundary type");
        break;
    }
}

void shape_block::generate_points(real const& globalspacing)
{
    switch (bound_type)
    {
    case linePlane:
        coords = LineShape::generate_points(*this, globalspacing);
        break;
    case squareCube:
        coords = SquareShape::generate_points(*this, globalspacing);
        break;
    case circleSphere:
        coords = CircleShape::generate_points(*this, globalspacing);
        break;
    case arcSection:
        coords = ArcShape::generate_points(*this, globalspacing);
        break;
    case cylinder:
        coords = CircleShape::generate_points(*this, globalspacing);
        break;
    case inletZone:
        coords = InletShape::generate_points(*this, globalspacing);
        break;
    case coordDef:
        break;
    default:
        printf("Unrecognised boundary type");
        break;
    }
}