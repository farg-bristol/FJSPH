#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Var.h"

#define PERTURB(i,j) pow(MEPSILON,pow(2,i*SIMDIM-j))

// #define X 0
// #define Y 1

inline StateMatD GetRotationMat(StateVecD& angles)
{
    #if SIMDIM == 3
        StateMatD rotx, roty, rotz;
        rotx << 1.0, 0.0            , 0.0           ,
                0.0, cos(angles(0)) , sin(angles(0)),
                0.0, -sin(angles(0)), cos(angles(0));

        roty << cos(angles(1)) , 0.0 , -sin(angles(1)),
                0.0            , 1.0 , 0.0            ,
                sin(angles(1)) , 0.0 , cos(angles(1));

        rotz << cos(angles(2)) , sin(angles(2)) , 0.0 ,
                -sin(angles(2)), cos(angles(2)) , 0.0 ,
                0.0            , 0.0            , 1.0 ;

        return rotx*roty*rotz;
    #else
        StateMatD rot;
        rot << cos(angles(0)), -sin(angles(0)),
               sin(angles(0)),  cos(angles(0));

        return rot;
    #endif
}

inline std::pair<StateVecD,StateVecD> Find_MinMax(SIM& svar, const SPHState& pnp1)
{
    /*Find the max and min positions*/
    auto xC = std::minmax_element(pnp1.begin(),pnp1.end(),
                [](SPHPart p1, SPHPart p2){return p1.xi(0)< p2.xi(0);});
    auto yC = std::minmax_element(pnp1.begin(),pnp1.end(),
                [](SPHPart p1, SPHPart p2){return p1.xi(1)< p2.xi(1);});
    #if SIMDIM == 3
        auto zC = std::minmax_element(pnp1.begin(),pnp1.end(),
                [](SPHPart p1, SPHPart p2){return p1.xi(2)< p2.xi(2);});

        StateVecD minC = StateVecD(xC.first->xi(0),yC.first->xi(1),zC.first->xi(2));
        StateVecD maxC = StateVecD(xC.second->xi(0),yC.second->xi(1),zC.second->xi(2));
    #else 
        StateVecD minC = StateVecD(xC.first->xi(0),yC.first->xi(1));
        StateVecD maxC = StateVecD(xC.second->xi(0),yC.second->xi(1));
    #endif  

    return std::pair<StateVecD,StateVecD>(minC,maxC);
}

/*Crossing test for 3 dimensions.*/
inline int LessThanREError(StateP1MatD const& A)
{
    real a1, a2, a3;

    /*Calculate components of the absolute*/
    a1 = fabs(A(0,2)-A(3,2))*(fabs((A(1,0)-A(3,0))*(A(2,1)-A(3,1)))+fabs((A(1,1)-A(3,1))*(A(2,0)-A(3,0))));
    a2 = fabs(A(1,2)-A(3,2))*(fabs((A(2,0)-A(3,0))*(A(0,1)-A(3,1)))+fabs((A(2,1)-A(3,1))*(A(0,0)-A(3,0))));
    a3 = fabs(A(2,2)-A(3,2))*(fabs((A(0,0)-A(3,0))*(A(1,1)-A(3,1)))+fabs((A(0,1)-A(3,1))*(A(1,0)-A(3,0))));

    if(fabs(A.determinant()) <=  MERROR*(a1+a2+a3))
    {
        // #pragma omp critical
        // {
        // cout << "Volume is inside the tolerance of round-off error. Need to do some form of tiebreaking." << endl;
        // cout << "Matrix: " << endl << A << endl;
        // cout << "MERROR: " << MERROR << endl;
        // cout << "Components: " << a1 << "  " << a2 << "  " << a3 << endl;
        // cout << A.determinant() <<  " <= " <<  MERROR*(a1+a2+a3) << endl;
        // }
        return 1;
    }

    return 0;
}

/*Surface detection as described by Marrone, Colagrossi, Le Touze, Graziani - (2010)*/
void Detect_Surface(SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t const& end,
                /*  DELTAP const& dp, */ OUTL const& outlist, MESH const& cells, SPHState& pnp1);

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
bool get_line_intersection(vector<StateVecD> const& verts, vector<size_t> const& edge, 
    StateVecD const& p1, StateVecD const& cellC);

/* ======= Crossings algorithm ============================================  */
/* By Eric Haines, 3D/Eye Inc, erich@eye.com                                 */
/* Shoot a test ray along +X axis.  The strategy, from MacMartin, is to      */
/* compare vertex Y values to the testing point's Y and quickly discard      */
/* edges which are entirely to one side of the test ray.                     */
/*                                                                           */
/* Input 2D polygon _pgon_ with _numverts_ number of vertices and test point */
/* _point_, returns 1 if inside, 0 if outside.  WINDING and CONVEX can be    */
/* defined for this test.                                                    */
#if SIMDIM == 2
    int Crossings2D(vector<StateVecD> const& verts, vector<size_t> const& edge, StateVecD const& point);

    /* Check for intersection with the infinite plane (I.e. just do two volumes) */
    /* Used for implicit particle tracking */
    int Cross_Plane(vector<StateVecD> const& verts, vector<size_t> const& face, 
        StateVecD const& point, StateVecD const& point2, bool& perturb);
#endif

#if SIMDIM == 3
    int Crossings3D(vector<StateVecD> const& verts, vector<size_t> const& face, 
        StateVecD const& point, StateVecD const& point2/* , uint& perturb */);

    /* Check for intersection with the infinite plane (I.e. just do two volumes) */
    /* Used for implicit particle tracking */
    int Cross_Plane(vector<StateVecD> const& verts, vector<size_t> const& face, 
        StateVecD const& point, StateVecD const& point2, bool& perturb);
#endif

/* Check for the distance from the face the point lies. It is already assumed that the ray crosses */
/* Used for implicit particle tracking */
void RayNormalIntersection(MESH const& cells, StateVecD const& rayOrigin, StateVecD const& rayVector,
                           vector<size_t> const& face, int const& cellID, real& dt, real& denom);

// Need cell elements to check size, and do the tet case. 
// Need cell centres for every other element
real Cell_Volume( vector<StateVecD> const& verts, vector<vector<size_t>> const& faces, vector<size_t> const& elems, 
				 vector<size_t> const& cell,  StateVecD const& cCentre);


#endif