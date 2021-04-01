/*Cell containment test*/
#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include "Eigen/Core"
#include "Eigen/StdVector"
#include "Eigen/Geometry"

// #include "Var.h"
// #include "IO.h"
// #include "CDFIO.h"

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 3
#endif

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using std::vector;
using std::cout;
using std::endl;

/* Define data type. */
#ifndef FOD
#define FOD 1 /*0 = float, 1 = double*/
#endif

typedef unsigned int uint;

#if FOD == 1
typedef double real;
#else
typedef float real;
#endif

/*Get machine bit precision for Simulation of Simplicity*/
#ifndef MEPSILON
#define MEPSILON std::numeric_limits<real>::epsilon() /*For  float, power is -24*/
#endif

#ifndef MERROR
#define MERROR (7.0*MEPSILON + 56.0*MEPSILON*MEPSILON)
#endif

#define PERTURB(i,j) pow(MEPSILON,pow(2,i*SIMDIM-j))

#define X 0
#define Y 1

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif
/****** Eigen vector definitions ************/

typedef Eigen::Matrix<real,SIMDIM,1> StateVecD;
typedef Eigen::Matrix<int,SIMDIM,1> StateVecI;
typedef Eigen::Matrix<real,SIMDIM,SIMDIM> RotMat;

/*Vector definitions for Density reinitialisation*/
typedef Eigen::Matrix<real, SIMDIM+1,1> DensVecD;
typedef Eigen::Matrix<real, SIMDIM+1, SIMDIM+1> DensMatD;





/*Crossing test for 3 dimensions.*/
int LessThanREError(const DensMatD& A)
{
    real a1, a2, a3;

    /*Calculate components of the absolute*/
    a1 = fabs(A(0,2)-A(3,2))*(fabs((A(1,0)-A(3,0))*(A(2,1)-A(3,1)))+fabs((A(1,1)-A(3,1))*(A(2,0)-A(3,0))));
    a2 = fabs(A(1,2)-A(3,2))*(fabs((A(2,0)-A(3,0))*(A(0,1)-A(3,1)))+fabs((A(2,1)-A(3,1))*(A(0,0)-A(3,0))));
    a3 = fabs(A(2,2)-A(3,2))*(fabs((A(0,0)-A(3,0))*(A(1,1)-A(3,1)))+fabs((A(0,1)-A(3,1))*(A(1,0)-A(3,0))));

    cout << "Determinant: " << A.determinant() << " Matrix Permanant error: " <<  MERROR*(a1+a2+a3) << endl;
    if(fabs(A.determinant()) <=  MERROR*(a1+a2+a3))
    {
        cout << "Volume is inside the tolerance of round-off error. Need to do some form of tiebreaking." << endl;
        // cout << "Matrix: " << endl << A << endl;
        // cout << "MERROR: " << MERROR << endl;
        // cout << "Components: " << a1 << "  " << a2 << "  " << a3 << endl;
        // cout << A.determinant() <<  " <= " <<  MERROR*(a1+a2+a3) << endl;
        return 1;
    }

    return 0;
}

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
int Crossings3D(const vector<StateVecD>& verts, const vector<size_t>& face, 
    StateVecD const& point, StateVecD const& point2)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996, & 
    Robust Adaptive Floating-Point Geometric Predicates
    Michael Aftosmis, Cart3D Software*/
    /*https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/bool_intersection.html*/
    StateVecD testp = point; 
    StateVecD rayp = point2;
    DensMatD vol1;
    int flag1, flag2;
    vol1 << verts[face[0]](0), verts[face[0]](1), verts[face[0]](2), 1.0,
            verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
            verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0,
            testp(0)         , testp(1)         , testp(2)         , 1.0;
     
    cout << "Matrix 1:" << endl;
    cout << vol1 << endl;

    if(LessThanREError(vol1))
    {   // Perturb the test point so that it doesn't go into roundoff error
		// cout << "Peturbing point: " << endl;
		// testp += StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
        // cout << "New point: " << testp(0) << "  " << testp(1) << "  " << testp(2) << endl;
		// vol1.row(0) << testp(0),testp(1),testp(2),1.0;


		cout << "Peturbing Matrix: " << endl;
        vol1 << 
	verts[face[0]](0)+ PERTURB(0,1), verts[face[0]](1)+ PERTURB(0,2), verts[face[0]](2)+ PERTURB(0,3), 1.0,
	verts[face[1]](0)+ PERTURB(1,1), verts[face[1]](1)+ PERTURB(1,2), verts[face[1]](2)+ PERTURB(1,3), 1.0,
    verts[face[2]](0)+ PERTURB(2,1), verts[face[2]](1)+ PERTURB(2,2), verts[face[2]](2)+ PERTURB(2,3), 1.0,
	testp(0) + PERTURB(3,1),testp(1) + PERTURB(3,2), testp(2) + PERTURB(3,3), 1.0;
    	cout << "New matrix:" << endl << vol1 << endl;

        
    
	    /*Check if the roundoff has been fixed*/
	    if(LessThanREError(vol1))
	    {
	        cout << "This is still the perturbed volume." << endl;
	    }
    }

    flag1 = (vol1.determinant() < 0);

    vol1 << verts[face[0]](0), verts[face[0]](1), verts[face[0]](2), 1.0,
            verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
            verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0,
            rayp(0)          , rayp(1)          , rayp(2)          , 1.0;

    cout << "Matrix 2:" << endl;
    cout << vol1 << endl;

    if(LessThanREError(vol1))
    {
       /*To add Simulation of Simplicity*/ 
        rayp += StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
        cout << "Peturbing point: " << endl;
    }
    


    flag2 = (vol1.determinant() < 0); 

    /*If signs of the volumes alternate, then the points lie either side of the plane*/
    if(flag1 != flag2)
    {        
        int flag3, flag4;
        StateVecD vtx0, vtx1;
        vtx0 = verts[face.back()]; /*Start on the last - first point edge*/
        vtx1 = verts[face[0]];

        /*Find initial volume size*/
        DensMatD vol;
        vol.row(0) << point(0), point(1), point(2), 1.0;
        vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
        vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
        vol.row(3) << point2(0), point2(1), point2(2),1.0;

        cout << "Inside Matrix 0:" << endl;
	    cout << vol << endl;
        if(LessThanREError(vol))
        {
            /*To add Simulation of Simplicity*/
        }
        
        flag3 = (vol.determinant() < 0);

        /*Check for each face, if the signs of all the tets are the same.*/
        for (size_t ii = 1; ii < face.size(); ++ii)
        {   /*Change the face vertices used*/
            // cout << face.size() << "  " << ii << endl;
            vtx0 = vtx1;
            vtx1 = verts[face[ii]];

            vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
            vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
            cout << "Inside Matrix " << ii << ":" << endl;
		    cout << vol << endl;
            if(LessThanREError(vol))
            {
                /*To add Simulation of Simplicity*/
            }

            flag4 = (vol.determinant() < 0);

            /*If the sign of the tet is different, this face isn't intersected.*/
            if (flag4 != flag3)
                return 0;
            
        }
        
        return 1;
            // std::cout << "Crosses face!" << std::endl;
    }
    
    return 0;    
}
#endif



int main(int argc, char* argv[])
{
	cout << std::scientific << std::left << std::setprecision(8);
	// cout << std::numeric_limits<double>::epsilon() << endl;
	// cout << pow(2,-52)<< endl;
	// /*Get the test settings. */
	// string infolder = "Inputs/";
	// string mesh = "M6.taumesh.faces";

	// SIM svar;
	// FLUID fvar;
	// AERO avar;
	// outl outlist;
	// MESH cells;

	// GetInput(argc,argv,svar,fvar,avar);
	
	// /*Read in the tau mesh*/
	// if(svar.Bcase == 6)
	// {
	// 	#if SIMDIM == 3
	// 	Read_TAUMESH_FACE(svar,cells,fvar);
	// 	// Read_TAUMESH(svar,cells,fvar);
	// 	#else
	// 	Read_TAUMESH_EDGE(svar,cells,fvar);
	// 	#endif
	// }	

	// /*Which vertex do you want to play near?*/
	// size_t vertex = 972780;

	vector<StateVecD> verts;
	verts.emplace_back(StateVecD(0.0001*MEPSILON,-3*MEPSILON,-3*MEPSILON));
	verts.emplace_back(StateVecD(-0.000001*MEPSILON,1*MEPSILON,-2*MEPSILON));
	verts.emplace_back(StateVecD(-0.000001*MEPSILON,-2*MEPSILON,1*MEPSILON));

	vector<size_t> face = {0,1,2};

	StateVecD testp = StateVecD(-2.12404e-5*MEPSILON,0.001*MEPSILON,0.001*MEPSILON);
	// StateVecD testp = StateVecD(-2.124039999999999e-5*MEPSILON,0.001*MEPSILON,0.001*MEPSILON);
	/*Verify we can determine a cross of the face safely*/
	StateVecD rayp = testp;
	rayp(0) += 1e+5;

	cout << "Test Point: " << testp(0) << "  " << testp(1) << "  " << testp(2) << endl;
	cout << "Ray Point:  " << rayp(0) << "  " << rayp(1) << "  " << rayp(2) << endl;

	uint inside_flag = 0;
	     
#if SIMDIM == 3        
    if(Crossings3D(verts,face,testp,rayp)) 
#else
    if(Crossings2D(verts,face,testp))
#endif   
    {   
        inside_flag=!inside_flag;
        // if ( line_flag ) break; //Convex assumption

        // //  note that one edge has been hit by the ray's line 
        // line_flag = TRUE;
    }
    

    cout << "Result:  " << inside_flag << endl;


	return 0;
}