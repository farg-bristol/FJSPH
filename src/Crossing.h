#ifndef CROSSING_H
#define CROSSING_H

#include "Eigen/Core"
#include "Var.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmpxx.h>

#define X 0
#define Y 1

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

/* ======= Crossings algorithm ============================================  */
/* By Eric Haines, 3D/Eye Inc, erich@eye.com                                 */
/* Shoot a test ray along +X axis.  The strategy, from MacMartin, is to      */
/* compare vertex Y values to the testing point's Y and quickly discard      */
/* edges which are entirely to one side of the test ray.                     */
/*                                                                           */
/* Input 2D polygon _pgon_ with _numverts_ number of vertices and test point */
/* _point_, returns 1 if inside, 0 if outside.  WINDING and CONVEX can be    */
/* defined for this test.                                                    */
int Crossings2D(std::vector<StateVecD> &pgon, StateVecD &point)
{
  int numverts = pgon.size();
  int i, j, yflag0, yflag1, inside_flag, line_flag;
  double  ty, tx;
  StateVecD vtx0, vtx1;

    tx = point[X] ;
    ty = point[Y] ;

    vtx0 = pgon[numverts-1] ;
    /* get test bit for above/below X axis */
    yflag0 = ( vtx0[Y] >= ty ) ;
    i = 0;
    vtx1 = pgon[i] ;

    inside_flag = 0 ;
    line_flag = 0;
    for ( j = numverts+1 ; --j ; ) 
    {
    yflag1 = ( vtx1[Y] >= ty ) ;
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
      inside_flag = !inside_flag ;
        }

        /* For convex cells, further optimisation can be done: */
        /* A ray can only pass through a maximum of two faces.*/
        /* If this is second edge hit, then done testing. */
        if ( line_flag ) goto Exit ;

        /* note that one edge has been hit by the ray's line */
        line_flag = TRUE ;
    }

    /* Move to the next pair of vertices, retaining info as possible. */
    yflag0 = yflag1 ;
    vtx0 = vtx1 ;
    ++i;
    vtx1 = pgon[i] ;
    }
    Exit: ;
    return( inside_flag ) ;
}


int Crossing2DPrecise(std::vector<StateVecD> &pgon, StateVecD &point)
{
  int numverts = pgon.size();
  int i, j, yflag0, yflag1, inside_flag, line_flag;
  mpf_class  ty, tx;
  mpf_class vtx0[SIMDIM], vtx1[SIMDIM];
  i = 0;

  tx = point[X];
  ty = point[Y];


    vtx0[0] = pgon[numverts-1][0]; vtx0[1] =  pgon[numverts-1][1];
    /* get test bit for above/below X axis */
    yflag0 = ( vtx0[Y] > ty ) ;
    
    // vtx1 = pgon[i] ;

    inside_flag = 0 ;
    line_flag = 0;
    for ( j = numverts+1 ; --j ; ) 
    {
    yflag1 = (vtx1[Y] > ty) ;
    /* Check if endpoints straddle (are on opposite sides) of X axis
     * (i.e. the Y's differ); if so, +X ray could intersect this edge.
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
        if ( ((vtx1[Y]-ty) * (vtx1[X]-vtx0[X]) >
          (vtx1[X]-tx) * (vtx1[Y]-vtx0[Y])) == yflag1 )
      {
      inside_flag = !inside_flag ;
        }

        /* For convex cells, further optimisation can be done: */
        /* A ray can only pass through a maximum of two faces.*/
        /* If this is second edge hit, then done testing. */
        if ( line_flag ) goto Exit ;

        /* note that one edge has been hit by the ray's line */
        line_flag = TRUE ;
    }

    /* Move to the next pair of vertices, retaining info as possible. */
    yflag0 = yflag1 ;
    vtx0[0] = vtx1[0]; vtx0[1] = vtx1[1];
    ++i;
    vtx1[0] = pgon[i][0]; vtx1[1] = pgon[i][1];
    }
    Exit: ;
    return( inside_flag ) ;    
}

/*Crossing test for 3 dimensions.*/
/* Cast a ray across the X-dimension, and use signed volumes*/
/*  of tetrahedra to identify if the ray crosses the plane */
int Crossings3D(std::vector<std::vector<StateVecD>> &pfaces, StateVecD &point)
{
  double d;
  StateVecD p2f, norm, tp;

  tp = point;
  for (auto face: pfaces) 
    { /*Check that the test point is on the correct side of each face 
        (requires convex, anti-clockwise faces) */

    /*Find the face normal*/
    norm = (face[1]-face[0]).cross(face[2]-face[0]);
    norm = norm.normalized();
      p2f = face[0] - tp;
      d = p2f.dot(norm);
      d /= p2f.norm();

      if(d < -1e-15)
      {
        return 0;
      }
    }

  return 1;
}

#endif