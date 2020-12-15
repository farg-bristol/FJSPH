#ifndef CROSSING_H
#define CROSSING_H

#include "Var.h"
#include "Geometry.h"
#include <iomanip>
// #include <gmpxx.h>


// #ifdef DEBUG
void Write_Containment(vector<size_t> const& ret_indexes, MESH const& cells, StateVecD const& testp)
{
    cout << "Search size: " << ret_indexes.size() << endl;
    vector<StateVecD> const& verts = cells.verts;

#if SIMDIM == 3
    cout << "Test point: " << testp(0) << " " << testp(1) << " " << testp(2) << endl;    
    cout << "First cell containment:" << endl;
    auto const index = ret_indexes[0];
    auto const& cell = cells.cFaces[index];
    cout << "Cell " <<  index << " volumes:" << endl;
    
    StateVecD const rayp(testp(0)+1e+5,testp(1),testp(2));
    

    
    uint intersects = 0;
    StateP1MatD vol1, vol2;
    for (auto const& findex:cell)
    {
        cout << "Face: " << findex << " left: " << cells.leftright[findex].first << " right: "
         << cells.leftright[findex].second  << endl;

        vector<size_t> const face = cells.faces[findex];
        vol1 << testp(0)         , testp(1)         , testp(2)         , 1.0,
                verts[face[0]](0), verts[face[0]](1), verts[face[0]](2), 1.0,
                verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
                verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0;
                

        vol2 = vol1;

        vol2.row(0) << rayp(0), rayp(1), rayp(2), 1.0;


        cout << "Volume 1: " << vol1.determinant() << "  Volume 2: " << vol2.determinant() << endl;

        if(LessThanREError(vol1))
        {
            cout << "Volume 1 needs perturbing" << endl;

            vol1.row(0) << testp(0)+PERTURB(0,1),testp(1)+PERTURB(0,2),testp(2)+PERTURB(0,3),1.0;

            cout << "New volume:" << vol1.determinant() << endl;
        }

        if(LessThanREError(vol2))
        {
            cout << "Volume 2 needs perturbing" << endl;

            vol2.row(0) << rayp(0)+PERTURB(0,1),rayp(1)+PERTURB(0,2),rayp(2)+PERTURB(0,3),1.0;

            cout << "New volume:" << vol2.determinant() << endl;
        }


        if((vol1.determinant()<0.0)!=(vol2.determinant()<0.0))
        {
            cout << "Could cross this face..." << endl;
            
            
            uint flag3,flag4;
            uint face_cross = 1;
            StateVecD vtx0, vtx1;
            vtx0 = verts[face.back()]; /*Start on the last - first point edge*/
            vtx1 = verts[face[0]];

            /*Find initial volume size*/
            StateP1MatD vol;
            vol.row(0) << testp(0), testp(1), testp(2), 1.0;
            vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
            vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
            vol.row(3) << rayp(0), rayp(1), rayp(2),1.0;
                        
            cout << "Volume 0: " << vol.determinant() << endl;

            if(LessThanREError(vol))
            {
                cout << "Volume needs perturbing" << endl;
                StateP1MatD pert = vol;

                pert.row(0) << testp(0)+PERTURB(0,1),testp(1)+PERTURB(0,2),testp(2)+PERTURB(0,3),1.0;

                flag3 = (pert.determinant() < 0.0);
            }
            else
                flag3 = (vol.determinant() < 0.0);

            /*Check for each face, if the signs of all the tets are the same.*/
            for (size_t ii = 1; ii < face.size(); ++ii)
            {   /*Change the face vertices used*/
                // cout << face.size() << "  " << ii << endl;
                vtx0 = vtx1;
                vtx1 = verts[face[ii]];

                vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
                vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;

                cout << "Volume " << ii << ": " << vol.determinant() << endl;
                if(LessThanREError(vol))
                {
                    cout << "Volume needs perturbing" << endl;
                    StateP1MatD pert = vol;

                    pert.row(0) << testp(0)+PERTURB(0,1),testp(1)+PERTURB(0,2),testp(2)+PERTURB(0,3),1.0;

                    flag4 = (pert.determinant() < 0.0);
                }
                else
                    flag4 = (vol.determinant() < 0.0);
                

                /*If the sign of the tet is different, this face isn't intersected.*/
                if(flag4 != flag3)
                {
                    face_cross = 0;
                }
            }

            if(face_cross == 1)
            {
                cout << "This face is crossed!!!!" << endl;
                intersects++;    
            }
        }

    }    
    cout << "Face intersections: " << intersects << endl;

    std::ofstream f1("FailedCellContainment.dat",std::ios::out);
    f1 << "VARIABLES= \"X\", \"Y\", \"Z\"" << endl;
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
        // for(auto const& faces:cells.cFaces[index])
        //     f1 << faces;

        // f1 << endl;

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
#endif

#if SIMDIM == 2
    cout << "Test point: " << testp(0) << " " << testp(1) << endl; 
    cout << "First cell containment:" << endl;
   
    auto const index = ret_indexes[0];
    auto const cell = cells.cFaces[index];
    cout << "Cell " <<  index << " Evaluations:" << endl;


    int  yflag0, yflag1, inside_flag;
    real  ty, tx;
    StateVecD vtx0, vtx1;

    tx = testp[0];
    ty = testp[1];

    inside_flag = 0;

    for (auto const& findex:cell)
    {
        auto const& edge = cells.faces[findex];
        vtx0 = verts[edge[0]];
        vtx1 = verts[edge[1]];
        /* Move to the next pair of vertices, retaining info as possible. */
        yflag0 = ( vtx0[1] >= ty );
        yflag1 = ( vtx1[1] >= ty );

        cout << "Left vertex: " << vtx0(0) << "  " << vtx0(1) << endl;
        cout << "Right vertex: " << vtx1(0) << "  " << vtx1(1) << endl;

        cout << "yflag0:  " << yflag0 << "  yflag1: " << yflag1 << endl;
        if ( yflag0 != yflag1 ) 
        {
            cout << "Left side: " << (vtx1[1]-ty) * (vtx1[0]-vtx0[0]) << "  Right side: "
            << (vtx1[0]-tx) * (vtx1[1]-vtx0[1]) << endl;

            if ( ((vtx1[1]-ty) * (vtx1[0]-vtx0[0]) >=
              (vtx1[0]-tx) * (vtx1[1]-vtx0[1])) == yflag1 )
            {
                cout << "Crosses face!" << endl;
                inside_flag = !inside_flag;
            }
        }            
    }

    cout << "Intersect? " << inside_flag << endl;
    

    std::ofstream f1("FailedCellContainment.dat",std::ios::out);
    f1 << "VARIABLES= \"X\", \"Z\"" << endl;
    f1 << "ZONE T=\"Test Point\", I=1, F=POINT" << endl;
    uint w = 25;
    
    f1 << std::left << std::scientific << std::setprecision(16);
    f1 << std::setw(w) << testp[0] << std::setw(w) << testp[1] << endl;

    for(auto const& index:ret_indexes)
    {  
        /*Number of points*/
        size_t size = cells.elems[index].size();
        uint numFaces = 0;                
        numFaces = size;

        /*Write zone header information*/
        f1 << "ZONE T=\"Cell " << index << "\""<< endl;
        f1 << "ZONETYPE=FEPOLYGON" << endl;
        f1 << "NODES=" << size << " ELEMENTS=" << 1 << " FACES=" << numFaces << endl;
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



#endif

}

uint CheckCell(size_t const& cell, MESH const& cells, StateVecD const& testp, uint& pert)
{
    // StateVecD testp = xi;
#if SIMDIM == 3    
    const StateVecD rayp(testp(0)+1e+5,testp(1),testp(2));
#endif

    uint line_flag = 0;
    uint inside_flag = 0;
    uint perturb = FALSE;
    for (auto const& cFaces:cells.cFaces[cell]) 
    {   // Step through cell faces
        vector<size_t> const& face = cells.faces[cFaces];
#if SIMDIM == 3            
        if(Crossings3D(cells.verts,face,testp,rayp,perturb))
#else
        if(Crossings2D(cells.verts,face,testp))
#endif                
        {   

            inside_flag=!inside_flag;
            if ( line_flag ) break; //Convex assumption
            // // //  note that one edge has been hit by the ray's line 
            line_flag = TRUE;
        }

        if(perturb == TRUE)
        {
            cout << "point needs perturbing" << endl;
            pert++;
            return 0;
        }
    }



    return inside_flag;
}

uint Check_Pipe(Vec_Tree& CELL_INDEX, const MESH& cells, const StateVecD& xi)
{
    StateVecD testp = xi;
    const size_t num_results = 100;
    vector<size_t> ret_indexes(num_results);
    vector<real> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<real> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    
    CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));
    
    // cout << ret_indexes.size() << endl;

    uint inside_flag=0;
    uint pert = 0;
    for(auto cell:ret_indexes)
    {   
        inside_flag = CheckCell(cell, cells, testp, pert);

        if(pert == 1)
        {
            // cout << "Need to perturb point..." << endl;
#if SIMDIM == 3
            testp = xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
            testp = xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
            inside_flag = CheckCell(cell,cells,testp, pert);  
        }

        if(pert > 1)
        {
            cout << "Still can't identify which cell point is in even after perturbation" << endl;       
        }
        
        if(inside_flag == TRUE)
        {
            break;
        }
    }

    if(inside_flag == FALSE)
    {
        cout << "Containing cell not found..." << endl;
        Write_Containment(ret_indexes, cells, testp);
    }

    return inside_flag;
}




void FirstCell(SIM& svar, size_t end, uint const ii, Vec_Tree const& CELL_INDEX,
     MESH const& cells, State& pnp1, State& pn)
{

    uint found = 0;
    StateVecD testp = pnp1[ii].xi;
    StateVecD rayp;
#if SIMDIM == 3    
    rayp = testp;
    rayp(0) += 1e+10;
    size_t const num_results = 50;
#else
    size_t const num_results = 10;
#endif

    vector<size_t> ret_indexes(num_results);
    vector<real> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<real> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    
    CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

    // cout << "Test Point: " << testp(0) << "  " << testp(1)  << "  " << testp(2) << endl;
    // cout << cells.cFaces.size() << endl;
    uint count = 0;
    for(auto cell:ret_indexes)
    {   
        uint pert = 0;
        if(CheckCell(cell,cells,testp,pert))
        {
            found = 1; 
            pnp1[ii].cellID = cell;
            pnp1[ii].cellV = cells.cVel[cell];
            pnp1[ii].cellP = cells.cP[cell];

            break;
        }

        if(pert == 1)
        {
#if SIMDIM == 3
            testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
            testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
            if(CheckCell(cell,cells,testp,pert))
            {
                found = 1; 
                pnp1[ii].cellID = cell;
                pnp1[ii].cellV = cells.cVel[cell];
                pnp1[ii].cellP = cells.cP[cell];

                break;
            }

            if(pert == 2)
            {
                cout << "Still can't identify which cell point is in even after perturbation" << endl; 
            }
        }
        count++;
    }

    if(found != 1)
    {
        /*Point could be across a boundary. Test if a ray crosses...*/
        cout << "Checking if particle if outside the mesh." << endl;
        /*Check if the ray from the point to the cell centre crosses a boundary face.*/
        uint cross = 0;
        for(auto index:ret_indexes)
        {
            rayp = cells.cCentre[index];            
            for (size_t const& findex:cells.cFaces[index] ) 
            {   /*If its a boundary face, check if the point crosses it*/
                if(cells.leftright[findex].second < 0)
                {
                    vector<size_t> const& face = cells.faces[findex];
                    int ints;
                    
#if SIMDIM == 3     
                    uint perturb = FALSE;              
                    ints = Crossings3D(cells.verts,face,testp,rayp,perturb);
                    if(perturb == TRUE)
                    {
#if SIMDIM == 3
                        testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
                        testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
                        ints = Crossings3D(cells.verts,face,testp,rayp,perturb);
                    }
#else               /*2D line intersection*/
                    ints = get_line_intersection(cells.verts,face,testp,rayp);
#endif
                    if(ints)
                    {
                        cross=!cross;
#ifdef DEBUG
                        // cout << "Particle has crossed a boundary!" << endl;
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
                            real plane = norm.dot(cells.verts[face[1]]);
                            real dist =  (plane - pnp1[ii].xi.dot(norm))/(norm.dot(norm));
                            pnp1[ii].xi = pnp1[ii].xi + dist*norm;

                        }
                        else if(cells.leftright[findex].second == -2)
                        {
                            #pragma omp critical
                            {
                            // cout << "Particle has crossed an outer boundary!" << endl;
                            // cout << "Particle will be deleted." << endl;
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
            #pragma omp single
            {
                cout << "First containing cell not found. Something is wrong." << endl;
                Write_Containment(ret_indexes,cells,testp);
                exit(-1);
            }
        }
    }   
}

void FindCell(SIM& svar, real const sr, KDTREE& TREE, MESH& cells, State& pnp1, State& pn)
{
    /*Find which cell the particle is in*/
    vector<size_t> toDelete;
    size_t const start = svar.bndPts;
    size_t const end = svar.totPts;

#pragma omp parallel 
    {
    vector<size_t> localDel;
    #pragma omp for schedule(static) nowait 
    for (size_t ii = start; ii < end; ++ii)
    {
        if (pnp1[ii].b == PartState.FREE_)
        {   
            StateVecD testp = pnp1[ii].xi;  
            uint inside_flag = 0;
            uint pert = 0;
            if(CheckCell(pnp1[ii].cellID,cells,testp,pert))
            {
                inside_flag = 1;
                continue;
            }

            if(pert == 1)
            {
#if SIMDIM == 3
                testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
                testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
                if(CheckCell(pnp1[ii].cellID,cells,testp,pert))
                {
                    inside_flag = 1;
                    continue;
                }

                if(pert == 2)
                {
                    cout << "Still can't identify which cell point is in even after perturbation" << endl; 
                }
            }

            if(inside_flag == 0)
            {   
                /*Perform a small search, since most cells are found within 1 or 2 indexes.*/
#if SIMDIM == 3
                size_t const num_results = 5;
#else
                size_t const num_results = 5;
#endif
                vector<size_t> ret_indexes(num_results);
                vector<real> out_dists_sqr(num_results);
#ifdef DEBUG
                // cout << "Having to perform neighbour search again. size: " << num_results << endl;
#endif

                nanoflann::KNNResultSet<real> resultSet(num_results);
                resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
                
    TREE.CELL.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));
                uint count = 0;
                for(auto const& cell:ret_indexes)
                {   
                    if(CheckCell(cell,cells,testp,pert))
                    {
#ifdef DEBUG
                        // cout << "Moved to a neighbour cell." << endl;
#endif
                        pnp1[ii].cellID = cell;
                        pnp1[ii].cellV = cells.cVel[cell];
                        pnp1[ii].cellP = cells.cP[cell];
                        inside_flag = 1;

                        // cout << "Found new cell at " << count << " in list" << endl;
                    }

                    if(pert == 1)
                    {
#if SIMDIM == 3
                        testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
                        testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
                        if(CheckCell(pnp1[ii].cellID,cells,testp,pert))
                        {
                            inside_flag = 1;
                            pnp1[ii].cellID = cell;
                            pnp1[ii].cellV = cells.cVel[cell];
                            pnp1[ii].cellP = cells.cP[cell];
                            break;
                        }
                    }

                    if(pert > 1)
                    {
                        cout << "Still can't identify which cell point is in even after perturbation" << endl; 
                    }

                    count++;
                }
            }

            
            if(inside_flag == 0)
            {   /*If first search fails to find the cell, perform a large search*/
                // In tets, the cell centre can be very far away, making the index large
#if SIMDIM == 3
                size_t const num_results = 100;
#else
                size_t const num_results = 150;
#endif
                vector<size_t> ret_indexes(num_results);
                vector<real> out_dists_sqr(num_results);
#ifdef DEBUG
                // cout << "Having to perform neighbour search again. size: " << num_results << endl;
#endif

                nanoflann::KNNResultSet<real> resultSet(num_results);
                resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
                
    TREE.CELL.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));
                uint count = 0;
                for(auto const& cell:ret_indexes)
                {   
                    if(CheckCell(cell,cells,testp,pert))
                    {
#ifdef DEBUG
                        // cout << "Moved to a neighbour cell." << endl;
#endif
                        pnp1[ii].cellID = cell;
                        pnp1[ii].cellV = cells.cVel[cell];
                        pnp1[ii].cellP = cells.cP[cell];
                        inside_flag = 1;

                        // cout << "Found new cell at " << count << " in list" << endl;
                    }

                    if(pert == 1)
                    {
#if SIMDIM == 3
                        testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
                        testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
                        if(CheckCell(pnp1[ii].cellID,cells,testp,pert))
                        {
                            inside_flag = 1;
                            pnp1[ii].cellID = cell;
                            pnp1[ii].cellV = cells.cVel[cell];
                            pnp1[ii].cellP = cells.cP[cell];
                            break;
                        }
                    }

                    if(pert > 1)
                    {
                        cout << "Still can't identify which cell point is in even after perturbation" << endl; 
                    }

                    count++;
                }

                if(inside_flag == 0)
                {
                    // cout << "Checking if particle has crossed a boundary." << endl;
                    // If still not found, then the point could be across a boundary.
                    // Check if the ray from the point to the cell centre crosses a boundary face.
                    uint cross = 0;
                    for(auto const& index:ret_indexes)
                    {
                        StateVecD rayp = cells.cCentre[index];            
                        for (size_t const& findex:cells.cFaces[index] ) 
                        {   /*If its a boundary face, check if the point crosses it*/
                            if(cells.leftright[findex].second < 0)
                            {
                                vector<size_t> const& face = cells.faces[findex];
                                int ints;
#if SIMDIM == 3                   
                                uint perturb = FALSE;              
                                ints = Crossings3D(cells.verts,face,testp,rayp,perturb);
                                if(perturb == TRUE)
                                {

#if SIMDIM == 3
                                    testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
#else
                                    testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
#endif
                                    ints = Crossings3D(cells.verts,face,testp,rayp,perturb);
                                }
#else                           /*2D line intersection*/
                                ints = get_line_intersection(cells.verts,face,testp,rayp);
#endif
                                if(ints)
                                {
                                    cross=!cross;
                                    if(cells.leftright[findex].second == -1)
                                    {
                                        #pragma omp critical
                                        {
                                        cout << "Particle has crossed an inner boundary!" << endl;
                                        }
                                    }
                                    else if(cells.leftright[findex].second == -2)
                                    {
                                        // cout << "Particle has crossed an outer boundary!" << endl; 
                                        localDel.emplace_back(ii);         
                                    }
                                }  
                            }
                        }
                    }

                    if(cross == 0)
                    {
                        #pragma omp critical
                        {
                            cout << "Particle " << ii-start << " containing cell not found. Something is wrong." << endl;
                            Write_Containment(ret_indexes,cells,testp);
                            exit(-1);
                        }
                    }
                }
            }
                 
        } //end if particle is valid
    }   // End of particle loop

    #pragma omp for schedule(static) ordered
    for(int ii=0; ii<omp_get_num_threads(); ii++)
    {
        #pragma omp ordered
            toDelete.insert(toDelete.end(),localDel.begin(),localDel.end());
    }


    }/*End of parallel*/

    if(toDelete.size()!=0)
    {
        std::sort(toDelete.begin(),toDelete.end());
        for(vector<size_t>::reverse_iterator index = toDelete.rbegin(); 
            index!=toDelete.rend(); ++index)
        {
            pnp1.erase(pnp1.begin()+*index);
            pn.erase(pn.begin()+*index);
        }
        for(auto& back:svar.back)
            back-=toDelete.size();
    }
}
#endif