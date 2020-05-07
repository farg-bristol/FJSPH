#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Var.h"

#define PERTURB(i,j) pow(MEPSILON,pow(2,i*SIMDIM-j))

#define X 0
#define Y 1

RotMat GetRotationMat(StateVecD& angles)
{
    if (SIMDIM == 3)
    {
        RotMat rotx, roty, rotz;
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
    }
    else if (SIMDIM == 2)
    {
        RotMat rot;
        rot << cos(angles(0)), -sin(angles(0)),
               sin(angles(0)),  cos(angles(0));

        return rot;
    }

    return RotMat::Zero();

}

std::pair<StateVecD,StateVecD> Find_MinMax(SIM& svar, const State& pnp1)
{
    /*Find the max and min positions*/
    auto xC = std::minmax_element(pnp1.begin(),pnp1.end(),
                [](Particle p1, Particle p2){return p1.xi(0)< p2.xi(0);});
    auto yC = std::minmax_element(pnp1.begin(),pnp1.end(),
                [](Particle p1, Particle p2){return p1.xi(1)< p2.xi(1);});
    #if SIMDIM == 3
        auto zC = std::minmax_element(pnp1.begin(),pnp1.end(),
                [](Particle p1, Particle p2){return p1.xi(2)< p2.xi(2);});

        StateVecD minC = StateVecD(xC.first->xi(0),yC.first->xi(1),zC.first->xi(2));
        StateVecD maxC = StateVecD(xC.second->xi(0),yC.second->xi(1),zC.second->xi(2));
    #else 
        StateVecD minC = StateVecD(xC.first->xi(0),yC.first->xi(1));
        StateVecD maxC = StateVecD(xC.second->xi(0),yC.second->xi(1));
    #endif  

    return std::pair<StateVecD,StateVecD>(minC,maxC);
}



/*Crossing test for 3 dimensions.*/
int LessThanREError(DensMatD const& A)
{
    real a1, a2, a3;

    /*Calculate components of the absolute*/
    a1 = fabs(A(0,2)-A(3,2))*(fabs((A(1,0)-A(3,0))*(A(2,1)-A(3,1)))+fabs((A(1,1)-A(3,1))*(A(2,0)-A(3,0))));
    a2 = fabs(A(1,2)-A(3,2))*(fabs((A(2,0)-A(3,0))*(A(0,1)-A(3,1)))+fabs((A(2,1)-A(3,1))*(A(0,0)-A(3,0))));
    a3 = fabs(A(2,2)-A(3,2))*(fabs((A(0,0)-A(3,0))*(A(1,1)-A(3,1)))+fabs((A(0,1)-A(3,1))*(A(1,0)-A(3,0))));

    if(fabs(A.determinant()) <=  MERROR*(a1+a2+a3))
    {
        #pragma omp critical
        {
        // cout << "Volume is inside the tolerance of round-off error. Need to do some form of tiebreaking." << endl;
        // cout << "Matrix: " << endl << A << endl;
        // cout << "MERROR: " << MERROR << endl;
        // cout << "Components: " << a1 << "  " << a2 << "  " << a3 << endl;
        // cout << A.determinant() <<  " <= " <<  MERROR*(a1+a2+a3) << endl;
        }
        return 1;
    }

    return 0;
}

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
bool get_line_intersection(vector<StateVecD> const& verts, vector<size_t> const& edge, 
    const StateVecD& p1, const StateVecD& cellC)
{
    const StateVecD& e1 = verts[edge[0]];
    const StateVecD& e2 = verts[edge[1]];
    StateVecD s,r;
    s = cellC - p1; r = e2 - e1;

    // Find the denominator of the calculation 
    real denom =  (-r(0) * s(1) + s(0) * r(1));

    // If the value of this is nearly 0, then 
    if(denom < MEPSILON)
    {   // Lines are colinear
        return 0;
    }

    real u, t;
    u = (-s(1) * (p1(0) - e1(0)) + s(0) * (p1(1) - e1(1))) / denom;
    t = ( r(0) * (p1(1) - e1(1)) - r(1) * (p1(0) - e1(0))) / denom;

    if (u > 0 && u < 1 && t > 0 && t < 1)
    {   // Collision detected
        return 1;
    }

    return 0; // No collision
}


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
    int Crossings2D(vector<StateVecD> const& verts, vector<size_t> const& edge, StateVecD const& point)
    {
        int  yflag0, yflag1, inside_flag;
        real  ty, tx;
        StateVecD vtx0, vtx1;

        tx = point[X];
        ty = point[Y];

        inside_flag = 0;

        vtx0 = verts[edge[0]];
        vtx1 = verts[edge[1]];
        /* Move to the next pair of vertices, retaining info as possible. */
        yflag0 = ( vtx0[Y] >= ty );
        yflag1 = ( vtx1[Y] >= ty );

        /* Check if endpoints straddle (are on opposite sides) of X axis
         * (i.e. the Y's differ); if so, +X ray could intersect this edge.
         * Credit to Joseph Samosky to try dropping
         * the "both left or both right" part of my code.
         */
        if ( yflag0 != yflag1 ) 
        {
            /* Check intersection of pgon segment with +X ray.
             * Note if >= point's X; if so, the ray hits it.
             * The division operation is avoided for the ">=" test by checking
             * the sign of the first vertex wrto the test point; idea inspired
             * by Joseph Samosky's and Mark Haigh-Hutchinson's different
             * polygon inclusion tests.
             */
            if ( ((vtx1[Y]-ty) * (vtx1[X]-vtx0[X]) >=
              (vtx1[X]-tx) * (vtx1[Y]-vtx0[Y])) == yflag1 )
            {
              inside_flag = !inside_flag;
            }

            /* For convex cells, further optimisation can be done: */
            /* A ray can only pass through a maximum of two faces.*/
            /* If this is second edge hit, then done testing. */ 
        }            
        return( inside_flag );
    }
#endif


#if SIMDIM == 3
int Crossings3D(vector<StateVecD> const& verts, vector<size_t> const& face, 
    StateVecD const& point, StateVecD const& point2, uint& perturb)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996, & 
    Robust Adaptive Floating-Point Geometric Predicates
    Michael Aftosmis, Cart3D Software*/
    /*https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/bool_intersection.html*/
    StateVecD const testp = point; 
    StateVecD const rayp = point2;
    DensMatD vol1;
    int flag1, flag2;
    vol1 << testp(0)         , testp(1)         , testp(2)         , 1.0,
            verts[face[0]](0), verts[face[0]](1), verts[face[0]](2), 1.0,
            verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
            verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0;
     

    if(LessThanREError(vol1))
    {   // Perturb the test point so that it doesn't go into roundoff error  
        perturb = TRUE;
        return 0;  
    }
    
    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0) << rayp(0), rayp(1), rayp(2),1.0;

    // ray point is very far so I don't see this falling into roundoff error 
    // if(LessThanREError(vol1))
    // {
    //     vol1.row(0) << rayp(0)+PERTURB(0,1), rayp(1)+PERTURB(0,2), rayp(2)+PERTURB(0,3),1.0;
    // }

    flag2 = (vol1.determinant() < 0.0); 

    /*If signs of the volumes alternate, then the points lie either side of the plane*/
    /*Now check if the line drawn by the two points intersects inside the bounds of the triangle plane*/
    if(flag1 != flag2)
    {   
        DensMatD vol;     
        int flag3, flag4;

        StateVecD vtx0, vtx1;
        vtx0 = verts[face.back()]; /*Start on the last - first point edge*/
        vtx1 = verts[face[0]];
        
        /*Find initial volume size*/
        
        vol.row(0) << testp(0), testp(1), testp(2), 1.0;
        vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
        vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
        vol.row(3) << rayp(0), rayp(1), rayp(2),1.0;

        if(LessThanREError(vol))
        {   /*Perturb the test point, since all volume calculations need to be done with the new location*/
            perturb = TRUE;
            return 0; 
        }
        
        flag3 = (vol.determinant() < 0.0);

        /*Check for each face, if the signs of all the tets are the same.*/
        for (size_t ii = 1; ii < face.size(); ++ii)
        {   /*Change the face vertices used*/
            // cout << face.size() << "  " << ii << endl;
            vtx0 = vtx1;
            vtx1 = verts[face[ii]];

            vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
            vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;

            if(LessThanREError(vol))
            {
                perturb = TRUE;
                return 0; 
            }
        
            flag4 = (vol.determinant() < 0.0);

            /*If the sign of the tet is different, this face isn't intersected.*/
            if (flag4 != flag3)
                return 0;   
        }  
        return 1;
    }  
    return 0;    
}
#endif

// Need cell elements to check size, and do the tet case. 
// Need cell centres for every other element
real Cell_Volume( vector<StateVecD> const& verts, vector<vector<size_t>> const& faces, vector<size_t> const& elems, 
				 vector<size_t> const& cell,  StateVecD const& cCentre)
{
#if SIMDIM == 3
	if(elems.size() == 4)
	{
		// Cell is a tetrahedron. Volume is trivial
		DensMatD vol;

		vol << verts[elems[0]](0), verts[elems[0]](1), verts[elems[0]](2), 1.0,
	           verts[elems[1]](0), verts[elems[1]](1), verts[elems[1]](2), 1.0,
	           verts[elems[2]](0), verts[elems[2]](1), verts[elems[2]](2), 1.0,
	           verts[elems[3]](0), verts[elems[3]](1), verts[elems[3]](2), 1.0; 

	    return abs(vol.determinant())/6.0;
	}
	
	// Otherwise, volume is not so trivial
	// Form tetrahedrons with the cell centre. 
	real sum = 0.0;
	
	for(auto const& faceID:cell)
	{
		DensMatD vol;

		const vector<size_t> face = faces[faceID];

		vol << cCentre(0), cCentre(1), cCentre(2), 1.0,
	           verts[face[0]](0), verts[face[0]](1), verts[face[0]](2), 1.0,
	           verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
	           verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0;

	    sum += abs(vol.determinant())/6.0;
	}
	
	return sum;
#else 

    if(elems.size() == 3)
    {
        DensMatD vol;

        vol << verts[elems[0]](0), verts[elems[0]](1), 1.0,
               verts[elems[1]](0), verts[elems[1]](1), 1.0,
               verts[elems[2]](0), verts[elems[2]](1), 1.0;

        return abs(vol.determinant())/2;
    }


    real sum = 0.0;
    for(auto const& faceID:cell)
    {
        DensMatD vol;
        
        vector<size_t> const& face = faces[faceID];

        vol << cCentre(0), cCentre(1), 1.0,
               verts[face[0]](0), verts[face[0]](1), 1.0,
               verts[face[1]](0), verts[face[1]](1), 1.0;

        sum += abs(vol.determinant())/2.0;
    }

    return sum;

#endif
}

void Make_Cell(FLUID const& fvar, AERO const& avar, MESH& cells)
{
	#if SIMDIM == 3

		cells.verts.emplace_back(StateVecD(-0.5,-0.5,-0.5));
		cells.verts.emplace_back(StateVecD(-0.5,-0.5,0.5));
		cells.verts.emplace_back(StateVecD(-0.5,0.5,0.5));
		cells.verts.emplace_back(StateVecD(-0.5,0.5,-0.5));
		cells.verts.emplace_back(StateVecD(0.5,-0.5,-0.5));
		cells.verts.emplace_back(StateVecD(0.5,-0.5,0.5));
		cells.verts.emplace_back(StateVecD(0.5,0.5,0.5));
		cells.verts.emplace_back(StateVecD(0.5,0.5,-0.5));
        cells.verts.emplace_back(StateVecD(1.5,-0.5,-0.5));
        cells.verts.emplace_back(StateVecD(1.5,-0.5,0.5));
        cells.verts.emplace_back(StateVecD(1.5,0.5,0.5));
        cells.verts.emplace_back(StateVecD(1.5,0.5,-0.5));

        // Cell 1
		cells.faces.emplace_back(vector<size_t>{0,1,2});
		cells.faces.emplace_back(vector<size_t>{0,2,3});
		cells.faces.emplace_back(vector<size_t>{0,3,7});
		cells.faces.emplace_back(vector<size_t>{0,7,4});
		cells.faces.emplace_back(vector<size_t>{0,5,1});
		cells.faces.emplace_back(vector<size_t>{0,4,5});
		cells.faces.emplace_back(vector<size_t>{1,5,6});
		cells.faces.emplace_back(vector<size_t>{1,6,2});
		cells.faces.emplace_back(vector<size_t>{3,2,6});
		cells.faces.emplace_back(vector<size_t>{3,6,7});

        // Interface
        cells.faces.emplace_back(vector<size_t>{4,6,5});
        cells.faces.emplace_back(vector<size_t>{4,7,6});

        // Cell 2
        cells.faces.emplace_back(vector<size_t>{4,7,11});
        cells.faces.emplace_back(vector<size_t>{4,11,8});
        cells.faces.emplace_back(vector<size_t>{4,9,5});
        cells.faces.emplace_back(vector<size_t>{4,8,9});
        cells.faces.emplace_back(vector<size_t>{5,9,10});
        cells.faces.emplace_back(vector<size_t>{5,10,6});
        cells.faces.emplace_back(vector<size_t>{8,10,9});
        cells.faces.emplace_back(vector<size_t>{8,11,10});
        cells.faces.emplace_back(vector<size_t>{7,6,10});
        cells.faces.emplace_back(vector<size_t>{7,10,11});

		cells.elems.emplace_back(vector<size_t>{0,1,2,3,4,5,6,7});
        cells.elems.emplace_back(vector<size_t>{4,5,6,7,8,9,10,11});


		cells.cFaces.emplace_back(vector<size_t>{0,1,2,3,4,5,6,7,8,9,10,11});
        cells.cFaces.emplace_back(vector<size_t>{10,11,12,13,14,15,16,17,18,19,20,21});

		cells.cCentre.emplace_back(StateVecD(0.0,0.0,0.0));
        cells.cCentre.emplace_back(StateVecD(1.0,0.0,0.0));

		for(uint ii = 0; ii < 10; ++ii)
			cells.leftright.emplace_back(std::pair<int,int>(0,-2));

        cells.leftright.emplace_back(std::pair<int,int>(0,1));
        cells.leftright.emplace_back(std::pair<int,int>(0,1));

        for(uint ii = 0; ii < 10; ++ii)
            cells.leftright.emplace_back(std::pair<int,int>(1,-2));

	#else 
		cells.verts.emplace_back(StateVecD(0.0,0.0));
		cells.verts.emplace_back(StateVecD(0.0,1.0));
		cells.verts.emplace_back(StateVecD(1.0,1.0));
		cells.verts.emplace_back(StateVecD(1.0,0.0));

		cells.faces.emplace_back(vector<size_t>{0,1});
		cells.faces.emplace_back(vector<size_t>{1,2});
		cells.faces.emplace_back(vector<size_t>{2,3});
		cells.faces.emplace_back(vector<size_t>{3,0});

		cells.elems.emplace_back(vector<size_t>{0,1,2,3});

		cells.cFaces.emplace_back(vector<size_t>{0,1,2,3});

		cells.cCentre.emplace_back(StateVecD(0.5,0.5));

		for(uint ii = 0; ii < 4; ++ii)
			cells.leftright.emplace_back(std::pair<int,int>(0,-2));
	#endif


	// Generate the solution vectors
	cells.cVel.emplace_back(avar.vInf);
	cells.cP.emplace_back(avar.pRef);
	// cells.SPHRho.emplace_back(fvar.rho0 * pow((fvar.gasPress/fvar.B + 1),1/fvar.gam));

	cells.cVol.emplace_back(1.0);
	cells.fNum.emplace_back(0.0);
	cells.fMass.emplace_back(fvar.simM);
	cells.vFn.emplace_back(StateVecD::Zero());
	cells.vFnp1.emplace_back(StateVecD::Zero());

	cells.cRho.emplace_back(avar.rhog);
	cells.cMass.emplace_back(avar.rhog); //Since volume = 1 for this cell

	cells.cPertn.emplace_back(StateVecD::Zero());
	cells.cPertnp1.emplace_back(StateVecD::Zero());

    // Second cell
    cells.cVel.emplace_back(avar.vInf);
    cells.cP.emplace_back(avar.pRef);
    // cells.SPHRho.emplace_back(fvar.rho0 * pow((fvar.gasPress/fvar.B + 1),1/fvar.gam));

    cells.cVol.emplace_back(1.0);
    cells.fNum.emplace_back(0.0);
    cells.fMass.emplace_back(fvar.simM);
    cells.vFn.emplace_back(StateVecD::Zero());
    cells.vFnp1.emplace_back(StateVecD::Zero());

    cells.cRho.emplace_back(avar.rhog);
    cells.cMass.emplace_back(avar.rhog); //Since volume = 1 for this cell

    cells.cPertn.emplace_back(StateVecD::Zero());
    cells.cPertnp1.emplace_back(StateVecD::Zero());

}


#endif