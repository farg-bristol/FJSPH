#ifndef CROSSING_H
#define CROSSING_H

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Var.h"
#include <iomanip>
// #include <gmpxx.h>

#define X 0
#define Y 1

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

/*Crossing test for 3 dimensions.*/
ldouble LessThanREError(const DensMatD& A)
{
    ldouble a1, a2, a3;

    /*Calculate components of the absolute*/
    a1 = fabs(A(0,2)-A(3,2))*(fabs((A(1,0)-A(3,0))*(A(2,1)-A(3,1)))+fabs((A(1,1)-A(3,1))*(A(2,0)-A(3,0))));
    a2 = fabs(A(1,2)-A(3,2))*(fabs((A(2,0)-A(3,0))*(A(0,1)-A(3,1)))+fabs((A(2,1)-A(3,1))*(A(0,0)-A(3,0))));
    a3 = fabs(A(2,2)-A(3,2))*(fabs((A(0,0)-A(3,0))*(A(1,1)-A(3,1)))+fabs((A(0,1)-A(3,1))*(A(1,0)-A(3,0))));

    return MERROR*(a1+a2+a3);
}

// #ifdef DEBUG
void Write_Containment(const vector<size_t>& ret_indexes, const MESH& cells, const StateVecD& testp)
{
    cout << ret_indexes.size() << endl;
    vector<StateVecD> const& verts = cells.verts;
    cout << "First three cell volumes:" << endl;
    for (size_t ii = 0; ii < 3; ii++)
    {
        auto index = ret_indexes[ii];
        auto cell = cells.cFaces[index];
        cout << "Cell " <<  index << " volumes:" << endl;

        DensMatD vol1;
        for (auto const& findex:cell)
        {
            const vector<size_t> face = cells.faces[findex];
            vol1.row(0) << verts[face[0]](0),verts[face[0]](1),verts[face[0]](2),1.0;
            vol1.row(1) << verts[face[1]](0),verts[face[1]](1),verts[face[1]](2),1.0;
            vol1.row(2) << verts[face[2]](0),verts[face[2]](1),verts[face[2]](2),1.0;
            vol1.row(3) << testp(0), testp(1), testp(2),1.0;
            cout << fabs(vol1.determinant()) << "  " << LessThanREError(vol1) << endl;
        }
    }

    std::ofstream f1("FailedCellContainment.dat",std::ios::out);
    f1 << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
    f1 << "ZONE T=\"Test Point\", I=1, F=POINT" << endl;
    uint w = 25;
    
    f1 << std::left << std::scientific << std::setprecision(16);
    f1 << std::setw(w) << testp[0] << std::setw(w) << testp[1] << std::setw(w) << testp[2] << endl;

    for(auto const& index:ret_indexes)
    {  
        /*Number of points*/
        size_t size = cells.elems[index].size();
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
        for(size_t dim = 0; dim < SIMDIM; ++dim)
        {
           for(auto const& cellv:cells.elems[index])
            {
                f1 << std::setw(w) << cells.verts[cellv][dim] ;   
            }
            f1 << endl;
        }

        /*Write how many vertices per face*/
        for(size_t ii = 0; ii < numFaces; ++ii)
        {
            f1  << std::setw(5) << 3;
        }
        f1 << endl;
        /*Write the face vertex list*/
        for(auto const& faces:cells.cFaces[index])
        {
            for(auto const& vertex:cells.faces[faces])
            {   /*How to get the vertex list to start from 1...*/
                /*State the position along the elems list*/
                auto it = std::find(cells.elems[index].begin(),cells.elems[index].end(),vertex);
                size_t ind = std::distance(cells.elems[index].begin(),it);
                f1 << std::setw(w) << ind+1;
            }
        }
        f1 << endl;
        /*Write left elements*/
        for(size_t ii = 0; ii < numFaces; ++ii)
        {
            f1 << std::setw(w) << 1;
        }
        f1 << endl;
        /*Write right elements*/
        for(size_t ii = 0; ii < numFaces; ++ii)
        {
            f1 << std::setw(w) << 0;
        }
        f1 << endl << endl;
    }

    f1.close();
}
// #endif


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
    ldouble denom =  (-r(0) * s(1) + s(0) * r(1));

    // If the value of this is nearly 0, then 
    if(denom < MEPSILON)
    {   // Lines are colinear
        return 0;
    }

    ldouble u, t;
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
        double  ty, tx;
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
int Crossings3D(const vector<StateVecD>& verts, const vector<size_t>& face, StateVecD const& point, StateVecD const& point2)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996, & 
    Robust Adaptive Floating-Point Geometric Predicates
    Michael Aftosmis, Cart3D Software*/
    /*https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/bool_intersection.html*/
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
        for (size_t ii = 1; ii < face.size(); ++ii)
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
                return 0;
            }
        }
        
        if(face_flag == 1)
        {
            // std::cout << "Crosses face!" << std::endl;
           return 1;
        }

    }
    
    return 0;    
}
#endif

void FirstCell(SIM& svar, size_t end, const uint ii, Vec_Tree& CELL_INDEX,
     const MESH& cells, const outl& outlist, State& pnp1, State& pn)
{

    uint found = 0;
    StateVecD testp = pnp1[ii].xi;
    StateVecD rayp;
#if SIMDIM == 3    
    rayp = testp;
    rayp(0) += 50000;
    const size_t num_results = 40;
#else
    const size_t num_results = 20;
#endif

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
        uint inside_flag = 0;
        uint line_flag = 0;

        for (uint const& findex:cells.cFaces[index] ) 
        {
            const vector<size_t>& face = cells.faces[findex];

#if SIMDIM == 3            
            if(Crossings3D(cells.verts,face,testp,rayp))
#else
            if(Crossings2D(cells.verts,face,testp))
#endif                
            {   
                inside_flag=!inside_flag;
                if ( line_flag ) goto wrongcell;

                /* note that one edge has been hit by the ray's line */
                line_flag = TRUE;
            }
        }

        if(inside_flag == TRUE)
        {
           
            pnp1[ii].cellID = index;
            pnp1[ii].cellV = cells.cVel[index];
            pnp1[ii].cellP = cells.cellP[index];
            pnp1[ii].cellRho = cells.cellRho[index];
            found = 1;
#ifdef DEBUG
                dbout << "Cell found: " << index << endl;
                dbout << "Tries: " << count << endl;
                dbout << "Test point: " << 
                testp[0] << "  " << testp[1];
#if SIMDIM == 3
                dbout << "  " << testp[2];
#endif
                dbout << endl;
                dbout << "Containing cell vertices: " << endl;
                dbout << "Cell number of elements: " << cells.elems[index].size() << endl;
                for(size_t kk = 0; kk < cells.elems[index].size(); ++kk)
                {
                    size_t elemI = cells.elems[index][kk];
                    dbout << kk << "  " << cells.verts[elemI][0] << " " 
                    << cells.verts[elemI][1];
#if SIMDIM == 3
                    dbout << " " << cells.verts[elemI][2];
#endif
                    dbout << endl;
                }
#endif
            break;
        }
wrongcell:
        count++;
    }

    if(found != 1)
    {
        /*Point could be across a boundary. Test if a ray crosses...*/
        cout << "Checking if particle has crossed a boundary." << endl;
        /*Check if the ray from the point to the cell centre crosses a boundary face.*/
        uint cross = 0;
        for(auto index:ret_indexes)
        {
            rayp = cells.cCentre[index];            
            for (size_t const& findex:cells.cFaces[index] ) 
            {   /*If its a boundary face, check if the point crosses it*/
                if(cells.leftright[findex].second < 0)
                {
                    const vector<size_t>& face = cells.faces[findex];
                    int ints;
#if SIMDIM == 3                   
                    ints = Crossings3D(cells.verts,face,testp,rayp);
#else               /*2D line intersection*/
                    ints = get_line_intersection(cells.verts,face,testp,rayp);
#endif
                    if(ints)
                    {
                        cross=!cross;
#ifdef DEBUG
                        cout << "Particle has crossed a boundary!" << endl;
#endif
                        if(cells.leftright[findex].second == -1)
                        {
                            cout << "Particle has crossed an inner boundary!" << endl;
                            StateVecD norm;
#if SIMDIM == 3
                            /*Get the face normal*/
                            StateVecD r1 = cells.verts[face[1]]- cells.verts[face[0]];
                            StateVecD r2 = cells.verts[face[2]]- cells.verts[face[0]];

                            norm = r1.cross(r2);
                            norm = norm.normalized();
#else
                            StateVecD r1 = cells.verts[face[1]]-cells.verts[face[0]];
                            norm = StateVecD(-r1(1),r1(0)); 
                            norm = norm.normalized();
#endif
                            /*Reflect the velocity away from the surface*/
                            pnp1[ii].v = pnp1[ii].v - 2*(pnp1[ii].v.dot(norm))*norm;

                        }
                        else if(cells.leftright[findex].second == -2)
                        {
                            #pragma omp critical
                            {
                            cout << "Particle has crossed an outer boundary!" << endl;
                            cout << "Particle will be deleted." << endl;
                            pnp1.erase(pnp1.begin()+ii);
                            pn.erase(pn.begin()+ii);
                            }
                            #pragma omp atomic
                                svar.totPts--;
                            #pragma omp atomic
                                end--;
                            
                        }
                    }  
                }
            }
        }

        if(cross == 0)
        {
            cout << "First containing cell not found. Something is wrong." << endl;
            Write_Containment(ret_indexes,cells,testp);
            exit(-1);
        }
    }   
}

void FindCell(SIM& svar, const ldouble nfull, Vec_Tree& CELL_INDEX,
     const MESH& cells, const outl& outlist, State& pnp1, State& pn)
{
    /*Find which cell the particle is in*/
    const size_t start = svar.bndPts;
    size_t end = svar.totPts;
    #pragma omp parallel for
    for (size_t ii = start; ii < end; ++ii)
    {
        if (pnp1[ii].b == 2 && outlist[ii].size() < nfull )
        {   
            StateVecD testp = pnp1[ii].xi;
            StateVecD rayp;
#if SIMDIM == 3
                
            rayp = testp;
            rayp(0) = 1e+100;
#endif
            uint inside_flag = 0;
            uint line_flag = 0;
            /*Test the cell it is already in first*/
            for (size_t const& findex:cells.cFaces[pnp1[ii].cellID] ) 
            {   /*Check each face of the cell*/
                const vector<size_t>& face = cells.faces[findex];

#if SIMDIM == 3            
                if(Crossings3D(cells.verts,face,testp,rayp))
#else
                if(Crossings2D(cells.verts,face,testp))
#endif                
                {   
                    inside_flag=!inside_flag;
                    if ( line_flag ) break;

                    /* note that one edge has been hit by the ray's line */
                    line_flag = TRUE;
                }
            
            }


            if(inside_flag == 1)
            {
#ifdef DEBUG
                // cout << "In the same cell" << endl;
#endif
                goto cellfound;
            }
            else
            {
                for(const auto cell:cells.cNeighb[pnp1[ii].cellID])
                {
                    uint inside_flag = 0;
                    uint line_flag = 0;
            
                    for (uint const& findex:cells.cFaces[cell] ) 
                    {
                        const vector<size_t>& face = cells.faces[findex];
#if SIMDIM == 3            
                        if(Crossings3D(cells.verts,face,testp,rayp))
#else
                        if(Crossings2D(cells.verts,face,testp))
#endif                
                        {   
                            inside_flag=!inside_flag;
                            if ( line_flag ) break;

                            /* note that one edge has been hit by the ray's line */
                            line_flag = TRUE;
                        }
                    }

                    if(inside_flag == 1)
                    {
                        pnp1[ii].cellID = cell;
                        pnp1[ii].cellV = cells.cVel[cell];
                        pnp1[ii].cellP = cells.cellP[cell];
                        pnp1[ii].cellRho = cells.cellRho[cell];
#ifdef DEBUG
                        // cout << "Moved to a neighbour cell." << endl;
#endif
                        goto cellfound;
                    }
                }

                /*The containing cell wasn't found in the neighbours.*/
                /*Scan through the whole list again*/
#if SIMDIM == 3
                const size_t num_results = 50;
#else
                const size_t num_results = 40;
#endif
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
                    uint inside_flag = 0;
                    uint line_flag = 0;

                    for (uint const& findex:cells.cFaces[index] ) 
                    {
                        const vector<size_t>& face = cells.faces[findex];
#if SIMDIM == 3            
                        if(Crossings3D(cells.verts,face,testp,rayp))
#else
                        if(Crossings2D(cells.verts,face,testp))
#endif                
                        {   
                            inside_flag=!inside_flag;
                            if ( line_flag ) break; /*Assumes convex cells*/

                            /* note that one edge has been hit by the ray's line */
                            line_flag = TRUE;
                        }
                    }

                    if(inside_flag == 1)
                    {
                        pnp1[ii].cellID = index;
                        pnp1[ii].cellV = cells.cVel[index];
                        pnp1[ii].cellP = cells.cellP[index];
                        pnp1[ii].cellRho = cells.cellRho[index];
#ifdef DEBUG
                        // cout << "Moved to a neighbour cell." << endl;
#endif
                        goto cellfound;
                    }
                }
                
                cout << "Checking if particle has crossed a boundary." << endl;
                // If still not found, then the point could be across a boundary.
                // Check if the ray from the point to the cell centre crosses a boundary face.
                uint cross = 0;
                for(auto index:ret_indexes)
                {
                    rayp = cells.cCentre[index];            
                    for (size_t const& findex:cells.cFaces[index] ) 
                    {   /*If its a boundary face, check if the point crosses it*/
                        if(cells.leftright[findex].second < 0)
                        {
                            const vector<size_t>& face = cells.faces[findex];
                            int ints;
#if SIMDIM == 3                   
                            ints = Crossings3D(cells.verts,face,testp,rayp);
#else                       /*2D line intersection*/
                            ints = get_line_intersection(cells.verts,face,testp,rayp);
#endif
                            if(ints)
                            {
                                cross=!cross;
#ifdef DEBUG
                                cout << "Particle has crossed a boundary!" << endl;
#endif
                                if(cells.leftright[findex].second == -1)
                                {
                                    cout << "Particle has crossed an inner boundary!" << endl;
                                    StateVecD norm;
#if SIMDIM == 3
                                    /*Get the face normal*/
                                    StateVecD r1 = cells.verts[face[1]]- cells.verts[face[0]];
                                    StateVecD r2 = cells.verts[face[2]]- cells.verts[face[0]];

                                    norm = r1.cross(r2);
                                    norm = norm.normalized();
#else
                                    StateVecD r1 = cells.verts[face[1]]-cells.verts[face[0]];
                                    norm = StateVecD(-r1(1),r1(0)); 
                                    norm = norm.normalized();
#endif
                                    /*Reflect the velocity away from the surface*/
                                    pnp1[ii].v = pnp1[ii].v - 2*(pnp1[ii].v.dot(norm))*norm;

                                }
                                else if(cells.leftright[findex].second == -2)
                                {
                                    #pragma omp critical
                                    {
                                    cout << "Particle has crossed an outer boundary!" << endl;
                                    cout << "Particle will be deleted." << endl;
                                    pnp1.erase(pnp1.begin()+ii);
                                    pn.erase(pn.begin()+ii);
                                    }
                                    #pragma omp atomic
                                        svar.totPts--;
                                    #pragma omp atomic
                                        end--;
                                    
                                }
                            }  
                        }
                    }
                }

                if(cross == 0)
                {
                    cout << "First containing cell not found. Something is wrong." << endl;
                    Write_Containment(ret_indexes,cells,testp);
                    exit(-1);
                }
                
                 
            }  
        }
        cellfound: ;
    }

}
#endif