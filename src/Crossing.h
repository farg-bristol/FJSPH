#ifndef CROSSING_H
#define CROSSING_H

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Var.h"
#ifdef DEBUG
#include <fstream>
#endif
// #include <gmpxx.h>

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
#if SIMDIM == 2
    int Crossings2D(vector<StateVecD> const& verts, vector<uint> const& pgon, StateVecD const& point)
    {
        
        int numverts = pgon.size();
        int  yflag0, yflag1, inside_flag, line_flag;
        double  ty, tx;
        StateVecD vtx0, vtx1;

        tx = point[X];
        ty = point[Y];

        vtx0 = verts[pgon[numverts-1]];
        /* get test bit for above/below X axis */
        yflag0 = ( vtx0[Y] >= ty );
        // i = 0;
        vtx1 = verts[pgon[0]];

        inside_flag = 0;
        line_flag = 0;

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

            /* note that one edge has been hit by the ray's line */
            line_flag = TRUE;
        }

        for (auto const& index:pgon) 
        {
            /* Move to the next pair of vertices, retaining info as possible. */
            yflag0 = yflag1;
            vtx0 = vtx1;
            vtx1 = verts[index];
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
                if ( line_flag ) return inside_flag;

                /* note that one edge has been hit by the ray's line */
                line_flag = TRUE;
            }            
        }
        
        return( inside_flag );
    }


    // int Crossing2DPrecise(std::vector<StateVecD> &pgon, StateVecD &point)
    // {
    //     int numverts = pgon.size();
    //     int i, j, yflag0, yflag1, inside_flag, line_flag;
    //     mpf_class  ty, tx;
    //     mpf_class vtx0[SIMDIM], vtx1[SIMDIM];
    //     i = 0;

    //     tx = point[X];
    //     ty = point[Y];


    //     vtx0[0] = pgon[numverts-1][0]; vtx0[1] =  pgon[numverts-1][1];
    //     /* get test bit for above/below X axis */
    //     yflag0 = ( vtx0[Y] > ty ) ;
        
    //     // vtx1 = pgon[i] ;

    //     inside_flag = 0 ;
    //     line_flag = 0;
    //     for ( j = numverts+1 ; --j ; ) 
    //     {
    //     yflag1 = (vtx1[Y] > ty) ;
    //     /* Check if endpoints straddle (are on opposite sides) of X axis
    //      * (i.e. the Y's differ); if so, +X ray could intersect this edge.
    //      */
    //     if ( yflag0 != yflag1 ) 
    //     {
    //          Check intersection of pgon segment with +X ray.
    //          * Note if >= point's X; if so, the ray hits it.
    //          * The division operation is avoided for the ">=" test by checking
    //          * the sign of the first vertex wrto the test point; idea inspired
    //          * by Joseph Samosky's and Mark Haigh-Hutchinson's different
    //          * polygon inclusion tests.
             
    //         if ( ((vtx1[Y]-ty) * (vtx1[X]-vtx0[X]) >
    //           (vtx1[X]-tx) * (vtx1[Y]-vtx0[Y])) == yflag1 )
    //       {
    //       inside_flag = !inside_flag ;
    //         }

    //         /* For convex cells, further optimisation can be done: */
    //         /* A ray can only pass through a maximum of two faces.*/
    //         /* If this is second edge hit, then done testing. */
    //         if ( line_flag ) goto Exit ;

    //         /* note that one edge has been hit by the ray's line */
    //         line_flag = TRUE ;
    //     }

    //     /* Move to the next pair of vertices, retaining info as possible. */
    //     yflag0 = yflag1 ;
    //     vtx0[0] = vtx1[0]; vtx0[1] = vtx1[1];
    //     ++i;
    //     vtx1[0] = pgon[i][0]; vtx1[1] = pgon[i][1];
    //     }
    //     Exit: ;
    //     return( inside_flag ) ;    
    // }
#endif


/*Crossing test for 3 dimensions.*/

ldouble LessThanREError(const DensMatD& A)
{
    ldouble a1, a2, a3;

    /*Calculate components of the absolute*/
    a1 = fabs(A(0,2)-A(3,2))*(fabs((A(1,0)-A(3,0))*(A(2,1)-A(3,1)))+fabs((A(1,1)-A(3,1))*(A(2,0)-A(3,0))));
    a2 = fabs(A(1,2)-A(3,2))*(fabs((A(2,0)-A(3,0))*(A(0,1)-A(3,1)))+fabs((A(2,1)-A(3,1))*(A(0,0)-A(3,0))));
    a3 = fabs(A(2,2)-A(3,2))*(fabs((A(0,0)-A(3,0))*(A(1,1)-A(3,1)))+fabs((A(0,1)-A(3,1))*(A(1,0)-A(3,0))));

    return (7*MEPSILON + 56*MEPSILON*MEPSILON)*(a1+a2+a3);
}

#if SIMDIM == 3
int Crossings3D(vector<StateVecD> const& verts, std::vector<std::vector<uint>> const& pfaces, 
                StateVecD const& point)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996, & 
    Robust Adaptive Floating-Point Geometric Predicates
    Michael Aftosmis, Cart3D Software*/
    /*https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/bool_intersection.html*/
    
    /*Take the unit vector of the point, and then add 50* it to make the second test point*/
    StateVecD point2 = point;
    point2(0) += 500;
    int inside_flag = 0;



    for (std::vector<uint> const& face: pfaces) 
    {
        DensMatD vol1;
        int flag1, flag2;
        vol1.row(0) << verts[face[0]](0),verts[face[0]](1),verts[face[0]](2),1.0;
        vol1.row(1) << verts[face[1]](0),verts[face[1]](1),verts[face[1]](2),1.0;
        vol1.row(2) << verts[face[2]](0),verts[face[2]](1),verts[face[2]](2),1.0;
        vol1.row(3) << point(0), point(1), point(2),1.0;

        if(fabs(vol1.determinant()) <= LessThanREError(vol1))
        {
            cout << "Volume is inside the tolerance of round-off error. Need to do some form of tiebreaking." << endl;
            cout << vol1.determinant() << endl;
        }

        flag1 = (vol1.determinant() < 0);

        vol1.row(3) << point2(0), point2(1), point2(2),1.0;
        if(fabs(vol1.determinant()) <= LessThanREError(vol1))
        {
            cout << "Volume is inside the tolerance of round-off error. Need to do some form of tiebreaking." << endl;
            cout << vol1.determinant() << endl;
        }


        flag2 = (vol1.determinant() < 0); 

        /*If signs of the volumes alternate, then the points lie either side of the plane*/
        if(flag1 != flag2)
        {
            uint numverts = face.size();
            
            int flag3, flag4, face_flag;
            face_flag = 1; /*Default to crossing the face first*/
            StateVecD vtx0, vtx1;
            vtx0 = verts[face[numverts-1]]; /*Start on the last - first point edge*/
            vtx1 = verts[face[0]];

            /*Find initial volume size*/
            DensMatD vol;
            vol.row(0) << point(0), point(1), point(2), 1.0;
            vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
            vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
            vol.row(3) << point2(0), point2(1), point2(2),1.0;
            if(fabs(vol.determinant()) <= LessThanREError(vol))
            {
            cout << "Volume is inside the tolerance of round-off error." << 
                    " Need to do some form of tiebreaking." << endl;
            cout << vol.determinant() << endl;
            }
            
            flag3 = (vol.determinant() < 0);

            /*Check for each face, if the signs of all the tets are the same.*/
            for (uint ii = 1; ii < face.size(); ++ii)
            {   /*Change the face vertices used*/
                // cout << face.size() << "  " << ii << endl;
                vtx0 = vtx1;
                vtx1 = verts[face[ii]];

                vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
                vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
                if(fabs(vol.determinant()) <= LessThanREError(vol))
                {
                cout << "Volume is inside the tolerance of round-off error." <<
                        " Need to do some form of tiebreaking." << endl;
                cout << vol.determinant() << endl;
                }

                flag4 = (vol.determinant() < 0);

                /*If the sign of the tet is different, this face isn't intersected.*/
                if (flag4 != flag3)
                {
                    face_flag = 0;
                    goto failedface;
                }
            }
            
            if(face_flag == 1)
            {
                // std::cout << "Crosses face!" << std::endl;
                inside_flag = !inside_flag;
            }
failedface: ;
        }
    }
    // std::cout << "inside_flag: " << inside_flag << std::endl;
    return inside_flag;
}
#endif

#ifdef DEBUG
    void Write_Containment(const vector<size_t>& ret_indexes, const MESH& cells, const StateVecD& testp)
    {
        cout << ret_indexes.size() << endl;
        vector<StateVecD> const& verts = cells.verts;
        cout << "First three cell volumes:" << endl;
        for (uint ii = 0; ii < 3; ii++)
        {
            auto index = ret_indexes[ii];
            auto cell = cells.cFaces[index];
            cout << "Cell " <<  static_cast<uint>(index) << " volumes:" << endl;

            DensMatD vol1;
            for (auto const& face:cell)
            {
                vol1.row(0) << verts[face[0]](0),verts[face[0]](1),verts[face[0]](2),1.0;
                vol1.row(1) << verts[face[1]](0),verts[face[1]](1),verts[face[1]](2),1.0;
                vol1.row(2) << verts[face[2]](0),verts[face[2]](1),verts[face[2]](2),1.0;
                vol1.row(3) << testp(0), testp(1), testp(2),1.0;
                cout << fabs(vol1.determinant()) << "  " << LessThanREError(vol1) << endl;
            }
        }

        std::ofstream f1("FailedCellContainment.dat",std::ios::out);
        // f1 << "VARIABLES = x, y, z" << endl;
        // f1 << "ZONE T=\"Test Point\", I=1, F=POINT" << endl;
        // f1 << std::scientific << std::setprecision(16);
        // f1 << testp[0] << "  " << testp[1] << "  " << testp[2] << endl;
        
        // for(auto const& index:ret_indexes)
        // {  
        //     uint size = cells.elems[index].size();
        //     f1 << "ZONE T=\"Cell " << index << "\", N=" << size << ", E=1, F=FEBLOCK," ;
        //     if(size ==4)
        //     {
        //         f1 << " ET=Tetrahedron" << endl;
        //     }
        //     else
        //     {
        //        f1 << " ET=Brick" << endl;
        //     }

        //     for(uint dim = 0; dim < SIMDIM; ++dim)
        //     {
        //        for(uint kk = 0; kk < cells.elems[index].size(); ++kk)
        //         {
        //             uint cellv = cells.elems[index][kk];
        //             f1 << cells.verts[cellv][dim] << " " ;
                    
        //         }
        //         f1 << endl;
        //     }
        //     if(size == 8)
        //     {   /*Hexahedron*/
        //         f1 << "1 2 3 4 5 7 8 6" << endl;
        //     }
        //     else if (size == 5)
        //     {   /*Pyramid*/
        //         f1 << "1 2 3 4 5 5 5 5" << endl;
        //     }
        //     else if (size == 6)
        //     {   /*Prism*/
        //         f1 << "1 2 3 3 4 5 6 6" << endl;
        //     }
        //     else if (size == 4)
        //     {   /*Tetrahedron*/
        //         f1 << "1 2 3 4" << endl;
        //     }
        // }

        

        f1 << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
        f1 << "ZONE T=\"Test Point\", I=1, F=POINT" << endl;
        uint w = 25;
        
        f1 << std::left << std::scientific << std::setprecision(16);
        f1 << std::setw(w) << testp[0] << std::setw(w) << testp[1] << std::setw(w) << testp[2] << endl;

        for(auto const& index:ret_indexes)
        {  
            /*Number of points*/
            uint size = cells.elems[index].size();
            uint numFaces = 0;
            uint TotalNumFaceNodes;
            if(size == 8)
            {   /*Hexahedron*/
                numFaces = 12;
            }
            else if (size == 6)
            {   /*Prism*/
                numFaces = 8;
            }
            else if (size == 5)
            {   /*Pyramid*/
                numFaces = 6;
            }
            else if (size == 4)
            {   /*Tetrahedron*/
                numFaces = 4;
            }
            
            TotalNumFaceNodes = numFaces*3;

            /*Write zone header information*/
            f1 << "ZONE T=\"Cell " << index << "\""<< endl;
            f1 << "ZONETYPE=FEPOLYHEDRON" << endl;
            f1 << "NODES=" << size << " ELEMENTS=" << 1 << " FACES=" << numFaces << endl;
            f1 << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
            f1 << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

            /*Write vertices in block format*/
            for(uint dim = 0; dim < SIMDIM; ++dim)
            {
               for(auto const& cellv:cells.elems[index])
                {
                    f1 << std::setw(w) << cells.verts[cellv][dim] ;   
                }
                f1 << endl;
            }

            /*Write how many vertices per face*/
            for(uint ii = 0; ii < numFaces; ++ii)
            {
                f1  << std::setw(5) << 3;
            }
            f1 << endl;
            /*Write the face vertex list*/
            for(auto const& faces:cells.cFaces[index])
            {
                for(auto const& vertex:faces)
                {   /*How to get the vertex list to start from 1...*/
                    /*State the position along the elems list*/
                    auto it = std::find(cells.elems[index].begin(),cells.elems[index].end(),vertex);
                    uint ind = std::distance(cells.elems[index].begin(),it);
                    f1 << std::setw(w) << ind+1;
                }
            }
            f1 << endl;
            /*Write left elements*/
            for(uint ii = 0; ii < numFaces; ++ii)
            {
                f1 << std::setw(w) << 1;
            }
            f1 << endl;
            /*Write right elements*/
            for(uint ii = 0; ii < numFaces; ++ii)
            {
                f1 << std::setw(w) << 0;
            }
            f1 << endl << endl;
        }

        f1.close();
    }
#endif

void FirstCell(const uint start, const uint end, const uint ii, Vec_Tree& CELL_INDEX,
     const MESH& cells, const outl& outlist, State& pnp1)
 {
        #if SIMDIM == 3
        uint found = 0;
        StateVecD testp = pnp1[ii].xi;
        const size_t num_results = 40;
        vector<size_t> ret_indexes(num_results);
        vector<double> out_dists_sqr(num_results);

        nanoflann::KNNResultSet<double> resultSet(num_results);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        
        CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

        // cout << "Test Point: " << testp(0) << "  " << testp(1)  << "  " << testp(2) << endl;
        // cout << cells.cFaces.size() << endl;
        uint count = 0;
        for(auto index:ret_indexes)
        {   
            const std::vector<std::vector<uint>> cell =
                 cells.cFaces[index];

            if(Crossings3D(cells.verts,cell,testp))
            {
                uint jj =  static_cast<uint>(index);
                pnp1[ii].cellID = jj;
                pnp1[ii].cellV = cells.cVel[jj];
                pnp1[ii].cellP = cells.cellP[jj];
                pnp1[ii].cellRho = cells.cellRho[jj];
                found = 1;
#ifdef DEBUG
                    dbout << "Cell found: " << jj << endl;
                    dbout << "Tries: " << count << endl;
                    dbout << "Test point: " << 
                    testp[0] << "  " << testp[1] << "  " << testp[2] << endl;
                    dbout << "Containing cell vertices: " << endl;
                    dbout << "Cell number of elements: " << cells.elems[jj].size() << endl;
                    for(uint kk = 0; kk < cells.elems[jj].size(); ++kk)
                    {
                        uint index = cells.elems[jj][kk];
                        dbout << kk << "  " << cells.verts[index][0] << " " 
                        << cells.verts[index][1] << " " << cells.verts[index][2] << endl;
                    }
#endif
                break;
            }
            count++;
        }

        if(found != 1)
        {
            cout << "Containing cell not found. Something is wrong." << endl;
#ifdef DEBUG
            Write_Containment(ret_indexes,cells,testp);
#endif
            exit(-1);
        }   
    #else
        uint found = 0;
        StateVecD testp = pnp1[ii].xi;
        /*Do a cell containment*/
        const size_t num_results = 20;
        vector<size_t> ret_indexes(num_results);
        vector<double> out_dists_sqr(num_results);

        nanoflann::KNNResultSet<double> resultSet(num_results);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        
        CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

        for(auto index:ret_indexes)
        {   
            // cout << "Cell: " << index << endl;
        const std::vector<uint> cell = 
            cells.elems[static_cast<uint>(index)];
            
            if(Crossings2D(cells.verts,cell,testp))
            {
                uint jj =  static_cast<uint>(index);
                pnp1[ii].cellID = jj;
                pnp1[ii].cellV = cells.cVel[jj];
                pnp1[ii].cellP = cells.cellP[jj];
                pnp1[ii].cellRho = cells.cellRho[jj];
                found = 1;
                // cout << "Cell found: " << jj << endl;
                // cout << "Found the containing cell!" << endl;
                // cout << jj << "  " << cells.cVel[jj][0] << endl;
                break;
                
            }
        }

        if(found != 1)
        {
        cout << "First containing cell not found. Something is wrong." << endl;
            exit(-1);
        }
    #endif
 }

void FindCell(const uint start, const uint end, const ldouble nfull, Vec_Tree& CELL_INDEX,
     const MESH& cells, const outl& outlist, State& pnp1)
{
    /*Find which cell the particle is in*/
    #pragma omp parallel for
    for (uint ii = start; ii < end; ++ii)
    {
        if (pnp1[ii].b == 2 && outlist[ii].size() < nfull )
        {   
            #if SIMDIM == 3
                uint found = 0;
                StateVecD testp = pnp1[ii].xi;
                /*Test the cell it is already in first*/
                if(Crossings3D(cells.verts,cells.cFaces[pnp1[ii].cellID],testp))
                {
                    found = 1;
                    #ifdef DEBUG
                    // cout << "In the same cell" << endl;
                    #endif
                    goto cellfound;
                }
                else
                {
                    for(const auto cell:cells.cNeighb[pnp1[ii].cellID])
                    {

                        if(Crossings3D(cells.verts,cells.cFaces[cell],testp))
                        {
                            pnp1[ii].cellID = cell;
                            pnp1[ii].cellV = cells.cVel[cell];
                            pnp1[ii].cellP = cells.cellP[cell];
                            pnp1[ii].cellRho = cells.cellRho[cell];
                            found = 1;
                            #ifdef DEBUG
                            // cout << "Moved to a neighbour cell." << endl;
                            #endif
                            goto cellfound;
                        }
                    }

                    if(found == 0)
                    {   /*The containing cell wasn't found in the neighbours.*/
                        /*Scan through the whole list again*/
                       
                        const size_t num_results = 40;
                        vector<size_t> ret_indexes(num_results);
                        vector<double> out_dists_sqr(num_results);
#ifdef DEBUG
                        cout << "Having to perform neighbour search again. size: " << num_results << endl;
#endif
                        

                        nanoflann::KNNResultSet<double> resultSet(num_results);
                        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
                        
        CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

                        

                        for(auto const& index:ret_indexes)
                        {   
                            auto cell = cells.cFaces[index];
                            // std::cout << "Cell: " << index << std::endl;
                            

                            if(Crossings3D(cells.verts,cell,testp))
                            {
                                uint jj =  static_cast<uint>(index);
                                pnp1[ii].cellID = jj;
                                pnp1[ii].cellV = cells.cVel[jj];
                                pnp1[ii].cellP = cells.cellP[jj];
                                pnp1[ii].cellRho = cells.cellRho[jj];
                                found = 1;
                                #ifdef DEBUG
                                // cout << "Cell found after a second search." << endl;
                                #endif
                                // cout << "Cell found: " << jj << endl;
                                goto cellfound;
                            }
                        }
                        if (found == 0)
                        {
                            std::cout << "Cell still wasn't found. Something needs changing." << std::endl;
                            #ifdef DEBUG
                            Write_Containment(ret_indexes,cells,testp);
                            #endif

                            exit(-1);
                        }
                    } 
                }  
            
            #else                   
                uint found = 0;
                StateVecD testp = pnp1[ii].xi;
                /*Do a cell containment*/
                if(Crossings2D(cells.verts,cells.elems[pnp1[ii].cellID],testp))
                {
                    found = 1;
                }
                else
                {
                    for(auto cell:cells.cNeighb[pnp1[ii].cellID])
                    {
                        if(Crossings2D(cells.verts,cells.elems[cell],testp))
                        {
                            pnp1[ii].cellID = cell;
                            pnp1[ii].cellV = cells.cVel[cell];
                            pnp1[ii].cellP = cells.cellP[cell];
                            pnp1[ii].cellRho = cells.cellRho[cell];
                            found = 1;
                            goto cellfound;
                        }
                    }

                    if(found != 1)
                    {   /*The containing cell wasn't found in the neighbours.*/
                        /*Redo the neighbour search around the test point*/
                        // std::cout << "Having to search tree again" << std::endl;
                        const size_t num_results = 40;
                        vector<size_t> ret_indexes(num_results);
                        vector<double> out_dists_sqr(num_results);

                        nanoflann::KNNResultSet<double> resultSet(num_results);
                        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
                        
        CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

                        for(auto index:ret_indexes)
                        {   
                            // std::cout << "Cell: " << index << std::endl;
                            auto cell = cells.elems[static_cast<uint>(index)];
                            if(Crossings2D(cells.verts,cell,testp))
                            {
                                uint jj =  static_cast<uint>(index);
                                pnp1[ii].cellID = jj;
                                pnp1[ii].cellV = cells.cVel[jj];
                                pnp1[ii].cellP = cells.cellP[jj];
                                pnp1[ii].cellRho = cells.cellRho[jj];
                                found = 1;
                                // cout << "Cell found: " << jj << endl;
                                goto cellfound;
                            }
                        }
                        
                    }

                    if(found != 1)
                    {   /*If it's still not found, do a huge search. Not pretty, but robust*/
                        /*The containing cell wasn't found in the neighbours.*/
                        /*Redo the neighbour search around the test point*/
                        // std::cout << "Expanding search..." << std::endl;
                        const size_t num_results = 50;
                        vector<size_t> ret_indexes(num_results);
                        vector<double> out_dists_sqr(num_results);

#ifdef DEBUG
                        cout << "Having to perform neighbour search again. size: " << num_results << endl;
#endif
                        
                        nanoflann::KNNResultSet<double> resultSet(num_results);
                        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
                        
        CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

                        for(auto index:ret_indexes)
                        {   
                            // std::cout << "Cell: " << index << std::endl;
                            auto cell = cells.elems[static_cast<uint>(index)];
                            if(Crossings2D(cells.verts,cell,testp))
                            {
                                uint jj =  static_cast<uint>(index);
                                pnp1[ii].cellID = jj;
                                pnp1[ii].cellV = cells.cVel[jj];
                                pnp1[ii].cellP = cells.cellP[jj];
                                pnp1[ii].cellRho = cells.cellRho[jj];
                                found = 1;
                                // cout << "Cell found: " << jj << endl;
                                goto cellfound;
                            }
                        }
                    }

                    if (found != 1)
                    {
        std::cout << "Cell still wasn't found. Something needs changing." << std::endl;
                        exit(-1);
                    }
                }
                
                
            #endif
        }
        cellfound: ;
    }

}
#endif