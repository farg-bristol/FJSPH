#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Var.h"
#include "IOFunctions.h"

#define PERTURB(i,j) pow(MEPSILON,pow(2,i*SIMDIM-j))

// #define X 0
// #define Y 1

/*Surface detection as described by Marrone, Colagrossi, Le Touze, Graziani - (2010)*/
void Detect_Surface(SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t const& end,
                 DELTAP const& dp, outl const& outlist, MESH const& cells, State& pnp1)
{
    #pragma omp parallel shared(pnp1)
    {
        real const h = 1.33 * svar.Pstep;

        vector<real> curves(end,0.0);
        vector<real> correcs(end,0.0);

        #pragma omp for schedule(static) nowait
        for(size_t ii = start; ii < end; ++ii)
        {
            if(dp.lam[ii] < 0.2)
            {   /*If less that 0.2, then its a surface particle by default*/
                pnp1[ii].surf = 1;
            }
            else if(dp.lam[ii] < 0.75)
            {   /*Particle could be a surface particle. Perform a test*/
            
                // Create point T
                StateVecD xi = pnp1[ii].xi;
                StateVecD pointT = xi + h * dp.norm[ii].normalized();

                uint surf = 1; /*Begin by assuming a surface, and proving otherwise*/
                for(size_t const& jj:outlist[ii])
                {
                    if(jj == ii)
                        continue;
                    
                    StateVecD x_jT = pnp1[jj].xi - pointT;
                    real r = (pnp1[jj].xi - xi).norm();

                    if(r >= sqrt(2)*h)
                    {
                        if (x_jT.norm() < h)
                        {
                            surf = 0;
                            break;
                        }
                    }
                    else  
                    {
#if SIMDIM == 2
                        StateVecD tau(dp.norm[ii](1),-dp.norm[ii](0));
                        if((abs(dp.norm[ii].normalized().dot(x_jT)) + abs(tau.normalized().dot(x_jT))) < h )
                        {
                            surf = 0;
                            break;
                        }
#else
                        StateVecD Rij = pnp1[jj].xi - xi; 
                        if(acos(dp.norm[ii].normalized().dot(Rij.normalized())) < M_PI/4.0)
                        {
                            surf = 0; 
                            break;
                        }
#endif
                    }
                    
                }
                pnp1[ii].surf = surf;
                
            }
            else
            {   /*If its eigenvalue is high, then by default it cannot be a surface*/
                pnp1[ii].surf = 0;
            }
        }

        #pragma omp for schedule(static) nowait
        for(size_t ii = start; ii < end; ++ii)
        {
            /*Find the curvature for the surface particles*/
            if(dp.lam[ii] < 0.75 )
            {
                real curve = 0.0;
                real correc = 0.0;

                real woccl_ = 0.0;
                StateVecD Vdiff = StateVecD::Zero();

                if(dp.lam[ii] < 0.75 && pnp1[ii].b == PartState.FREE_)
                {
                    if (svar.Asource == 1)
                    {
                        Vdiff = (pnp1[ii].cellV) - pnp1[ii].v/* dp.avgV[ii]*/;
                    }
                    else if (svar.Asource == 2)
                    {
                        Vdiff = (pnp1[ii].cellV+cells.cPertnp1[pnp1[ii].cellID]) - pnp1[ii].v /*dp.avgV[ii]*/;
                    }
                    else 
                    {
                        Vdiff = avar.vInf - pnp1[ii].v /*dp.avgV[ii]*/;
                    }

#if SIMDIM == 3
                    if(svar.Asource == 3)
                    {   
                        StateVecD Vel = svar.vortex.getVelocity(pnp1[ii].xi);
                        Vdiff = Vel - pnp1[ii].v /*dp.avgV[ii]*/;
                    }
#endif
                }

                for(size_t const& jj:outlist[ii])
                {
                    if (ii == jj)
                        continue;

                    StateVecD Rij = pnp1[jj].xi - pnp1[ii].xi;
                    real r = Rij.norm();
                    real volj = pnp1[jj].m/pnp1[jj].rho;
                    StateVecD diffK = volj *  GradK(Rij,r,fvar.H,fvar.correc);


                    /*Occlusion for Gissler Aero Method*/
                    if (pnp1[ii].b == PartState.FREE_ && (avar.acase == 4 || avar.acase == 1))
                    {
                        real const frac = -Rij.dot(Vdiff)/(Vdiff.norm()*r);
                        
                        if (frac > woccl_)
                        {
                            woccl_ = frac;
                        }
                    } 

                    if (pnp1[jj].surf == 1)
                    {
                        curve -= (dp.norm[jj].normalized()-dp.norm[ii].normalized()).dot(diffK);
                        correc += volj * Kernel(r,fvar.H,fvar.correc);
                        // curve += curves[jj] * pnp1[jj].m/pnp1[jj].rho* kern / dp.kernsum[ii];
                    }
                }

                curves[ii] = curve/correc;
                correcs[ii] = correc;
                pnp1[ii].woccl = woccl_; 
            }
        }

        #pragma omp for schedule(static) nowait
        for(size_t ii = start; ii < end; ++ii)
        {
            /*Find the curvature for the surface particles*/
            if(pnp1[ii].surf == 1)
            {
                real curve = 0.0;

                for(size_t const& jj:outlist[ii])
                {
                    if (pnp1[jj].surf != 1)
                        continue;

                    if (ii == jj)
                    {
                        curve += curves[ii] * fvar.correc/correcs[ii];
                        continue;
                    }

                    real r = (pnp1[jj].xi - pnp1[ii].xi).norm();
                    real kern = pnp1[jj].m/pnp1[jj].rho * Kernel(r,fvar.H,fvar.correc);

                    curve += curves[jj] * kern/correcs[ii];
                }

                pnp1[ii].curve = curve;
            }
        }   
    }
}


StateMatD GetRotationMat(StateVecD& angles)
{
    if (SIMDIM == 3)
    {
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
    }
    else if (SIMDIM == 2)
    {
        StateMatD rot;
        rot << cos(angles(0)), -sin(angles(0)),
               sin(angles(0)),  cos(angles(0));

        return rot;
    }

    return StateMatD::Zero();

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
int LessThanREError(StateP1MatD const& A)
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

        tx = point[0];
        ty = point[1];

        inside_flag = 0;

        vtx0 = verts[edge[0]];
        vtx1 = verts[edge[1]];
        /* Move to the next pair of vertices, retaining info as possible. */
        yflag0 = ( vtx0[1] >= ty );
        yflag1 = ( vtx1[1] >= ty );

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
            if ( ((vtx1[1]-ty) * (vtx1[0]-vtx0[0]) >=
              (vtx1[0]-tx) * (vtx1[1]-vtx0[1])) == yflag1 )
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
    StateP1MatD vol1;
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
        StateP1MatD vol;     
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
		StateP1MatD vol;

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
		StateP1MatD vol;

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
        StateP1MatD vol;

        vol << verts[elems[0]](0), verts[elems[0]](1), 1.0,
               verts[elems[1]](0), verts[elems[1]](1), 1.0,
               verts[elems[2]](0), verts[elems[2]](1), 1.0;

        return abs(vol.determinant())/2;
    }


    real sum = 0.0;
    for(auto const& faceID:cell)
    {
        StateP1MatD vol;
        
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


void Set_Mass(SIM& svar, FLUID& fvar, AERO& avar, State& pn, State& pnp1)
{

    if(svar.Bcase == 4)
    {
        real volume = pow(svar.Pstep,SIMDIM)*real(svar.totPts);
    #if SIMDIM == 3
        real dvol = (4.0*M_PI/3.0)*pow(0.5*svar.diam,3.0);
    #else
        real dvol = M_PI*pow(0.5*svar.diam,2.0);
    #endif
        
        cout << "SPH Volume: " << volume << "  Droplet expected volume: " << dvol << endl;
        cout << std::scientific;
        svar.mass = pnp1[0].m * svar.totPts;

        cout << "SPH Mass:    " << svar.mass << "     Droplet expected mass: " << dvol*fvar.rho0;   

        cout << " Error: " << std::fixed << 100.0*(dvol*fvar.rho0-svar.mass)/(dvol*fvar.rho0) << "%" << endl;   
        cout << std::scientific;

        cout << "Old SPH particle mass: " << fvar.simM << "  spacing: " << svar.Pstep << endl;
        // Adjust mass and spacing to create the correct mass/volume for the droplet.

        fvar.simM = dvol*fvar.rho0/real(svar.totPts);
        // svar.Pstep = pow(fvar.simM/fvar.rho0,1.0/SIMDIM);



    }
    else if (svar.Bcase == 3)
    {   /*Jet flow*/
        /*Take the height of the jet, and find the volume of the cylinder*/
#if SIMDIM == 3 
        real cVol = svar.Jet(1) * (M_PI * pow((svar.Jet(0)/2.0),2));
#else
        real cVol = svar.Jet(1) * svar.Jet(0);
#endif
        real volume = pow(svar.Pstep,SIMDIM)*real(svar.simPts);

        cout << "SPH Volume: " << volume << "  Starting jet expected volume: " << cVol << endl;
        cout << std::scientific;
        svar.mass = pnp1[0].m * real(svar.simPts);

        cout << "SPH Mass:    " << svar.mass << "     Starting jet expected mass: " << cVol*fvar.rho0;   

        cout << " Error: " << std::fixed << 100.0*(cVol*fvar.rho0-svar.mass)/(cVol*fvar.rho0) << "%" << endl;   
        cout << std::scientific;

        cout << "Old SPH particle mass: " << fvar.simM << "  spacing: " << svar.Pstep << endl;
        // Adjust mass and spacing to create the correct mass/volume for the droplet.

        fvar.simM = cVol*fvar.rho0/real(svar.simPts);

    }

#if SIMDIM == 3
    // svar.Pstep = 2*pow((3.0*fvar.simM)/(4.0*M_PI*fvar.rho0),1.0/3.0);
    svar.Pstep = pow(fvar.simM/fvar.rho0,1.0/3.0);
#else
    // svar.Pstep = 2*sqrt(fvar.simM/(fvar.rho0*M_PI));
    svar.Pstep = pow(fvar.simM/fvar.rho0,1.0/2.0);

#endif
    svar.dx = svar.Pstep * pow(fvar.rho0/fvar.rhoJ,1.0/SIMDIM);
    GetYcoef(avar, fvar, svar.Pstep);
    fvar.H = 2.0*svar.Pstep;
    fvar.HSQ = fvar.H*fvar.H; 

    fvar.sr = 4.0*fvar.HSQ;   /*KDtree search radius*/

    fvar.dCont = fvar.delta * fvar.H * fvar.Cs;
    fvar.dMom = fvar.dCont * fvar.rho0;


#if SIMDIM == 3
    fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
    // fvar.correc = (1.0/(M_PI*fvar.H*fvar.H*fvar.H));

    avar.pVol = 4.0/3.0 * M_PI * pow(avar.L,SIMDIM);
    // avar.aPlate = svar.Pstep*svar.Pstep;
    avar.aPlate = 4.0*avar.L*avar.L;
#else
    fvar.correc = 7.0/(4.0*M_PI*fvar.H*fvar.H);
    // fvar.correc = 10.0/(7.0*M_PI*fvar.H*fvar.H);

    avar.pVol = M_PI* avar.L*avar.L/4.0;
    avar.aPlate = svar.Pstep;
    // avar.aPlate = 2.0*avar.L;
#endif

    fvar.dCont = 2.0 * fvar.delta * fvar.H * fvar.Cs;
    fvar.artMu = std::max(fvar.mu, fvar.alpha * fvar.Cs * fvar.H * fvar.rho0);
    fvar.Wdx = Kernel(svar.Pstep, fvar.H, fvar.correc);

    for(size_t ii = 0; ii < svar.totPts; ii++)
    {
        pn[ii].m = fvar.simM;
        pnp1[ii].m = fvar.simM;
    }

    svar.mass = pnp1[0].m * svar.simPts;
    cout << "New SPH particle mass: " << fvar.simM << "  spacing: " << svar.Pstep << endl;
}

#endif