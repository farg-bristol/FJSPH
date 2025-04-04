/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Containment.h"
#include "Geometry.h"
#include <iomanip>
#include <set>

void Write_Containment(
    SIM const& svar, vector<size_t> const& ret_indexes, MESH const& cells, StateVecD const& testp
)
{
    cout << "Search size: " << ret_indexes.size() << endl;
    vector<StateVecD> const& verts = cells.verts;

#if SIMDIM == 3
    cout << "Test point: " << testp(0) << " " << testp(1) << " " << testp(2) << endl;
    cout << "First cell containment:" << endl;
    auto const index = ret_indexes[0];
    auto const& cell = cells.cFaces[index];
    cout << "Cell " << index << " volumes:" << endl;

    StateVecD const rayp(testp(0) + 1e+5, testp(1), testp(2));

    uint intersects = 0;
    StateP1MatD vol1, vol2;
    for (auto const& findex : cell)
    {
        cout << "Face: " << findex << " left: " << cells.leftright[findex].first
             << " right: " << cells.leftright[findex].second << endl;

        vector<size_t> const face = cells.faces[findex];
        vol1 << testp(0), testp(1), testp(2), 1.0, verts[face[0]](0), verts[face[0]](1),
            verts[face[0]](2), 1.0, verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
            verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0;

        vol2 = vol1;

        vol2.row(0) << rayp(0), rayp(1), rayp(2), 1.0;

        cout << "Volume 1: " << vol1.determinant() << "  Volume 2: " << vol2.determinant() << endl;

        if (LessThanREError(vol1))
        {
            cout << "Volume 1 needs perturbing" << endl;

            vol1.row(0) << testp(0) + PERTURB(0, 1), testp(1) + PERTURB(0, 2), testp(2) + PERTURB(0, 3),
                1.0;

            cout << "New volume:" << vol1.determinant() << endl;
        }

        if (LessThanREError(vol2))
        {
            cout << "Volume 2 needs perturbing" << endl;

            vol2.row(0) << rayp(0) + PERTURB(0, 1), rayp(1) + PERTURB(0, 2), rayp(2) + PERTURB(0, 3),
                1.0;

            cout << "New volume:" << vol2.determinant() << endl;
        }

        if ((vol1.determinant() < 0.0) != (vol2.determinant() < 0.0))
        {
            cout << "Could cross this face..." << endl;

            uint flag3, flag4;
            uint face_cross = 1;
            StateVecD vtx0, vtx1;
            vtx0 = verts[face.back()]; /*Start on the last - first point edge*/
            vtx1 = verts[face[0]];

            /*Find initial volume size*/
            StateP1MatD vol;
            vol.row(0) << testp(0), testp(1), testp(2), 1.0;
            vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
            vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;
            vol.row(3) << rayp(0), rayp(1), rayp(2), 1.0;

            cout << "Volume 0: " << vol.determinant() << endl;

            if (LessThanREError(vol))
            {
                cout << "Volume needs perturbing" << endl;
                StateP1MatD pert = vol;

                pert.row(0) << testp(0) + PERTURB(0, 1), testp(1) + PERTURB(0, 2),
                    testp(2) + PERTURB(0, 3), 1.0;

                flag3 = (pert.determinant() < 0.0);
            }
            else
                flag3 = (vol.determinant() < 0.0);

            /*Check for each face, if the signs of all the tets are the same.*/
            for (size_t ii = 1; ii < 3; ++ii)
            { /*Change the face vertices used*/
                // cout << face.size() << "  " << ii << endl;
                vtx0 = vtx1;
                vtx1 = verts[face[ii]];

                vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
                vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;

                cout << "Volume " << ii << ": " << vol.determinant() << endl;
                if (LessThanREError(vol))
                {
                    cout << "Volume needs perturbing" << endl;
                    StateP1MatD pert = vol;

                    pert.row(0) << testp(0) + PERTURB(0, 1), testp(1) + PERTURB(0, 2),
                        testp(2) + PERTURB(0, 3), 1.0;

                    flag4 = (pert.determinant() < 0.0);
                }
                else
                    flag4 = (vol.determinant() < 0.0);

                /*If the sign of the tet is different, this face isn't intersected.*/
                if (flag4 != flag3)
                {
                    face_cross = 0;
                }
            }

            if (face_cross == 1)
            {
                cout << "This face is crossed!!!!" << endl;
                intersects++;
            }
        }
    }
    cout << "Face intersections: " << intersects << endl;
    string fout = svar.io.output_prefix;
    fout.append("_FailedCellContainment.dat");

    cout << "Writing to file: " << fout << endl;

    std::ofstream f1(fout, std::ios::out);
    f1 << "VARIABLES= \"X\", \"Y\", \"Z\"" << endl;
    f1 << "ZONE T=\"Test Point\", I=1, F=POINT" << endl;
    uint w = 25;

    f1 << std::left << std::scientific << std::setprecision(16);
    f1 << std::setw(w) << testp[0] / svar.scale << std::setw(w) << testp[1] / svar.scale << std::setw(w)
       << testp[2] / svar.scale << endl;

    for (auto const& index : ret_indexes)
    {
        uint numFaces = 0;
        uint TotalNumFaceNodes = 0;
        uint nQuad = 0;
        uint nTrig = 0;
        std::set<size_t> elem;

        /* Build the element unique index list */
        for (auto const& faceID : cells.cFaces[index])
        {
            vector<size_t> const face = cells.faces[faceID];
            for (auto const& vert : face)
            { /* Use set to ensure unique set of vertices in elem */
                elem.insert(vert);
            }
        }

        /*Number of points*/
        size_t size = elem.size();

        if (size == 8)
        { /*Hexahedron*/
            numFaces = 6;
            TotalNumFaceNodes = numFaces * 4;
            nQuad = 6;
            nTrig = 0;
        }
        else if (size == 6)
        { /*Prism*/
            numFaces = 5;
            TotalNumFaceNodes = 3 * 4 + 2 * 3; /* 3 square faces, 2 triangle faces */
            nQuad = 3;
            nTrig = 2;
        }
        else if (size == 5)
        { /*Pyramid*/
            numFaces = 5;
            TotalNumFaceNodes = 4 + 4 * 3; /* 1 square faces, 4 triangle faces */
            nQuad = 1;
            nTrig = 4;
        }
        else if (size == 4)
        { /*Tetrahedron*/
            numFaces = 4;
            TotalNumFaceNodes = 12;
            nQuad = 0;
            nTrig = 4;
        }

        /*Write zone header information*/
        f1 << "ZONE T=\"Cell " << index << "\"" << endl;
        f1 << "ZONETYPE=FEPOLYHEDRON" << endl;
        f1 << "NODES=" << size << " ELEMENTS=" << 1 << " FACES=" << numFaces << endl;
        f1 << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
        f1 << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

        /*Write vertices in block format*/
        for (size_t dim = 0; dim < SIMDIM; ++dim)
        {
            for (auto const& cellv : elem)
            {
                f1 << std::setw(w) << cells.verts[cellv][dim] / cells.scale;
            }
            f1 << endl;
        }
        // for(auto const& faces:cells.cFaces[index])
        //     f1 << faces;

        // f1 << endl;

        /*Write how many vertices per face*/
        for (size_t ii = 0; ii < nQuad; ++ii)
        {
            f1 << std::setw(5) << 4;
        }
        for (size_t ii = 0; ii < nTrig; ++ii)
        {
            f1 << std::setw(5) << 3;
        }
        f1 << endl;

        /*Write the face vertex list*/
        for (auto const& faces : cells.cFaces[index])
        {
            for (auto const& vertex : cells.faces[faces])
            {
                /*state the position along the elems list*/
                auto it = std::find(elem.begin(), elem.end(), vertex);
                size_t ind = std::distance(elem.begin(), it);
                f1 << std::setw(w) << ind + 1;
            }
        }
        f1 << endl;
        /*Write left elements*/
        for (size_t ii = 0; ii < numFaces; ++ii)
        {
            f1 << std::setw(w) << 1;
        }
        f1 << endl;
        /*Write right elements*/
        for (size_t ii = 0; ii < numFaces; ++ii)
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
    cout << "Cell " << index << " Evaluations:" << endl;

    int yflag0, yflag1, inside_flag;
    real ty, tx;
    StateVecD vtx0, vtx1;

    tx = testp[0];
    ty = testp[1];

    inside_flag = 0;

    for (auto const& findex : cell)
    {
        auto const& edge = cells.faces[findex];
        vtx0 = verts[edge[0]];
        vtx1 = verts[edge[1]];
        /* Move to the next pair of vertices, retaining info as possible. */
        yflag0 = (vtx0[1] >= ty);
        yflag1 = (vtx1[1] >= ty);

        cout << "Left vertex: " << vtx0(0) << "  " << vtx0(1) << endl;
        cout << "Right vertex: " << vtx1(0) << "  " << vtx1(1) << endl;

        cout << "yflag0:  " << yflag0 << "  yflag1: " << yflag1 << endl;
        if (yflag0 != yflag1)
        {
            cout << "Left side: " << (vtx1[1] - ty) * (vtx1[0] - vtx0[0])
                 << "  Right side: " << (vtx1[0] - tx) * (vtx1[1] - vtx0[1]) << endl;

            if (((vtx1[1] - ty) * (vtx1[0] - vtx0[0]) >= (vtx1[0] - tx) * (vtx1[1] - vtx0[1])) == yflag1)
            {
                cout << "Crosses face!" << endl;
                inside_flag = !inside_flag;
            }
        }
    }

    cout << "Intersect? " << inside_flag << endl;

    string fout = svar.io.output_prefix;
    fout.append("_FailedCellContainment.dat");

    cout << "Writing to file: " << fout << endl;

    std::ofstream f1(fout, std::ios::out);
    f1 << "VARIABLES= \"X\", \"Z\"" << endl;
    f1 << "ZONE T=\"Test Point\", I=1, F=POINT" << endl;
    uint w = 25;

    f1 << std::left << std::scientific << std::setprecision(16);
    f1 << std::setw(w) << testp[0] << std::setw(w) << testp[1] << endl;

    for (size_t const& index : ret_indexes)
    {
        /* Build the element unique index list */
        vector<size_t> elem;
        for (size_t const& faceID : cells.cFaces[index])
        {
            vector<size_t> const face = cells.faces[faceID];
            for (size_t const& vert : face)
            {
                if (std::find(elem.begin(), elem.end(), vert) == elem.end())
                { /*Vertex doesn't exist in the elems vector yet.*/
                    elem.emplace_back(vert);
                }
            }
        }

        /*Number of points*/
        size_t size = elem.size();
        uint numFaces = size;

        /*Write zone header information*/
        f1 << "ZONE T=\"Cell " << index << "\"" << endl;
        f1 << "ZONETYPE=FEPOLYGON" << endl;
        f1 << "NODES=" << size << " ELEMENTS=" << 1 << " FACES=" << numFaces << endl;
        f1 << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

        /*Write vertices in block format*/
        for (size_t dim = 0; dim < SIMDIM; ++dim)
        {
            for (auto const& cellv : elem)
            {
                f1 << std::setw(w) << cells.verts[cellv][dim];
            }
            f1 << endl;
        }

        /*Write the face vertex list*/
        for (auto const& faces : cells.cFaces[index])
        {
            for (auto const& vertex : cells.faces[faces])
            { /* State the position along the elems list*/
                auto it = std::find(elem.begin(), elem.end(), vertex);
                size_t ind = std::distance(elem.begin(), it);
                f1 << std::setw(w) << ind + 1;
            }
        }
        f1 << endl;

        /*Write left elements*/
        for (size_t ii = 0; ii < numFaces; ++ii)
        {
            f1 << std::setw(w) << 1;
        }
        f1 << endl;

        /*Write right elements*/
        for (size_t ii = 0; ii < numFaces; ++ii)
        {
            f1 << std::setw(w) << 0;
        }
        f1 << endl << endl;
    }

    f1.close();

#endif
}

uint CheckCell(size_t const& cell, MESH const& cells, StateVecD const& testp)
{
// StateVecD testp = xi;
#if SIMDIM == 3
    const StateVecD rayp(testp(0) + 1e+5, testp(1), testp(2));
#endif

    uint line_flag = 0;
    uint inside_flag = 0;
    for (auto const& cFaces : cells.cFaces[cell])
    { // Step through cell faces
        vector<size_t> const& face = cells.faces[cFaces];
#if SIMDIM == 3
        if (Crossings3D(cells.verts, face, testp, rayp))
#else
        if (Crossings2D(cells.verts, face, testp))
#endif
        {

            inside_flag = !inside_flag;
            if (line_flag)
                break; // Convex assumption
            // // //  note that one edge has been hit by the ray's line
            line_flag = TRUE;
        }

        // if(perturb == TRUE)
        // {
        //     cout << "point needs perturbing" << endl;
        //     pert++;
        //     return 0;
        // }
    }

    return inside_flag;
}

/* <summary> Find the first cell for particles transitioning from the pipe to being in free air.  */
/* It is assumed that the particle does not have a previous cell defined, and so goes */
/* immediately to the KD Tree to find the nearest cells to iterate through. </summary> */
void FirstCell(SIM& svar, Vec_Tree const& CELL_INDEX, MESH const& cells, SPHPart& pi, uint& to_del)
{
    uint found = 0;
    StateVecD testp = pi.xi;
    StateVecD rayp;
#if SIMDIM == 3
    rayp = testp;
    rayp(0) += 1e+10;
    size_t const num_results = 150;
#else
    size_t const num_results = 20;
#endif

    vector<size_t> ret_indexes(num_results);
    vector<real> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<real> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    CELL_INDEX.index->findNeighbors(resultSet, &testp[0]);

    // cout << "Test Point: " << testp(0) << "  " << testp(1)  << "  " << testp(2) << endl;
    // cout << cells.cFaces.size() << endl;
    // uint count = 0;
    for (auto const& cell : ret_indexes)
    {
        // uint pert = 0;
        if (CheckCell(cell, cells, testp))
        {
            found = 1;
            pi.cellID = cell;
            pi.cellV = cells.cVel[cell];
            pi.cellP = cells.cP[cell];
            pi.cellRho = cells.cRho[cell];

            break;
        }

        // if(pert == 1)
        // {
        //     #if SIMDIM == 3
        //         testp = pi.xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
        //     #else
        //         testp = pi.xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
        //     #endif

        //     if(CheckCell(cell,cells,testp))
        //     {
        //         found = 1;
        //         pi.cellID = cell;
        //         pi.cellV = cells.cVel[cell];
        //         pi.cellP = cells.cP[cell];
        //         pi.cellRho = cells.cRho[cell];

        //         break;
        //     }

        //     if(pert == 2)
        //     {
        //         cout << "Still can't identify which cell point is in even after perturbation" << endl;
        //     }
        // }
        // count++;
    }

    if (found != 1)
    {
        /*Point could be across a boundary. Test if a ray crosses...*/
        cout << "Checking if particle if outside the mesh." << endl;
        /*Check if the ray from the point to the cell centre crosses a boundary face.*/
        uint cross = 0;
        for (auto index : ret_indexes)
        {
            rayp = cells.cCentre[index];
            for (size_t const& findex : cells.cFaces[index])
            { /*If its a boundary face, check if the point crosses it*/
                if (cells.leftright[findex].second < 0)
                {
                    vector<size_t> const& face = cells.faces[findex];
                    int ints;

#if SIMDIM == 3
                    // uint perturb = FALSE;
                    ints = Crossings3D(cells.verts, face, testp, rayp);
                    // if(perturb == TRUE)
                    // {
                    //     #if SIMDIM == 3
                    //     testp = pi.xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
                    //     #else
                    //     testp = pi.xi + StateVecD(PERTURB(0, 1), PERTURB(0, 2));
                    //     #endif
                    //     ints = Crossings3D(cells.verts,face,testp,rayp);
                    // }
#else /*2D line intersection*/
                    ints = get_line_intersection(cells.verts, face, testp, rayp);
#endif

                    if (ints)
                    {
                        cross = !cross;
#ifdef DEBUG
// cout << "Particle has crossed a boundary!" << endl;
#endif

                        if (cells.leftright[findex].second == -1)
                        {
                            cout << "Particle has crossed an inner boundary!" << endl;
                            StateVecD norm;

#if SIMDIM == 3
                            /*Get the face normal*/
                            StateVecD r1 = cells.verts[face[1]] - cells.verts[face[0]];
                            StateVecD r2 = cells.verts[face[2]] - cells.verts[face[0]];

                            norm = r1.cross(r2);
                            norm = norm.normalized();
#else
                            StateVecD r1 = cells.verts[face[1]] - cells.verts[face[0]];
                            norm = StateVecD(-r1(1), r1(0));
                            norm = norm.normalized();
#endif

                            /*Reflect the velocity away from the surface*/
                            pi.v = pi.v - 2 * (pi.v.dot(norm)) * norm;
                            real plane = norm.dot(cells.verts[face[1]]);
                            real dist = (plane - pi.xi.dot(norm)) / (norm.dot(norm));
                            pi.xi = pi.xi + dist * norm;
                        }
                        else if (cells.leftright[findex].second == -2)
                        {
                            to_del = 1;
                        }
                    }
                }
            }
        }

        if (cross == 0)
        {
#pragma omp single
            {
                cout << "First containing cell for particle " << pi.part_id
                     << " not found. Something is wrong." << endl;
                Write_Containment(svar, ret_indexes, cells, testp);
                exit(-1);
            }
        }
    }
}

/* <summary> Find the cells for all non-boundary particles. Checks for whether a paricle is free.  */
/* It is assumed that the particle does have a previous cell defined, and checks this cell */
/* before going to the KD Tree to find the nearest cells to iterate through. </summary> */
vector<size_t>
FindCell(SIM& svar, Vec_Tree const& CELL_TREE, MESH const& cells, SPHState& pn, SPHState& pnp1)
{
    /*Find which cell the particle is in*/
    vector<size_t> toDelete;
    size_t const& start = svar.bound_points;
    size_t const& end = svar.total_points;

#pragma omp parallel default(shared)
    {
        vector<size_t> localDel;
#pragma omp for schedule(static) nowait
        for (size_t ii = start; ii < end; ++ii)
        {
            if (pnp1[ii].b != FREE || pnp1[ii].lam_nb > svar.air.lam_cutoff)
            {
                pnp1[ii].cellID = c_no_cell;
                continue;
            }

            StateVecD testp = pnp1[ii].xi;
            uint inside_flag = 0;
            if (CheckCell(pnp1[ii].cellID, cells, testp))
            { /* Particle is still in the same cell, so do nothing */
                inside_flag = 1;
                pnp1[ii].cellV = cells.cVel[pnp1[ii].cellID];
                pnp1[ii].cellP = cells.cP[pnp1[ii].cellID];
                pnp1[ii].cellRho = cells.cRho[pnp1[ii].cellID];
                pnp1[ii].ipt_n_failed = 0;
                continue;
            }

            // if(pert == 1)
            // {
            //     #if SIMDIM == 3
            //         testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
            //     #else
            //         testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
            //     #endif

            //     if(CheckCell(pnp1[ii].cellID,cells,testp,pert))
            //     {
            //         inside_flag = 1;
            //         continue;
            //     }

            //     if(pert == 2)
            //     {
            //         cout << "Still can't identify which cell point is in even after perturbation" <<
            //         endl;
            //     }
            // }

            if (inside_flag == 0)
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

                CELL_TREE.index->findNeighbors(resultSet, &testp[0]);
                // uint count = 0;
                for (auto const& cell : ret_indexes)
                {
                    if (CheckCell(cell, cells, testp))
                    {
#ifdef DEBUG
// cout << "Moved to a neighbour cell." << endl;
#endif

                        pnp1[ii].cellID = cell;
                        pnp1[ii].cellV = cells.cVel[cell];
                        pnp1[ii].cellP = cells.cP[cell];
                        pnp1[ii].cellRho = cells.cRho[cell];
                        pnp1[ii].ipt_n_failed = 0;
                        pnp1[ii].internal = 0;
                        inside_flag = 1;
                        break;
                        // cout << "Found new cell at " << count << " in list" << endl;
                    }

                    // if(pert == 1)
                    // {
                    //     #if SIMDIM == 3
                    //         testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
                    //     #else
                    //         testp = pnp1[ii].xi + StateVecD(PERTURB(0,1),PERTURB(0,2));
                    //     #endif

                    //     if(CheckCell(pnp1[ii].cellID,cells,testp,pert))
                    //     {
                    //         inside_flag = 1;
                    //         pnp1[ii].cellID = cell;
                    //         pnp1[ii].cellV = cells.cVel[cell];
                    //         pnp1[ii].cellP = cells.cP[cell];
                    //         pnp1[ii].cellRho = cells.cRho[cell];
                    //         pnp1[ii].ipt_n_failed = 0;
                    //         break;
                    //     }
                    // }

                    // if(pert > 1)
                    // {
                    //     cout << "Still can't identify which cell point is in even after perturbation"
                    //     << endl;
                    // }

                    // count++;
                }
            }

            if (inside_flag == 0)
            { /*If first search fails to find the cell, perform a large search*/
                // In tets, the cell centre can be very far away, making the index large
                size_t const num_results = 500;
                vector<size_t> ret_indexes(num_results);
                vector<real> out_dists_sqr(num_results);

#ifdef DEBUG
// cout << "Having to perform neighbour search again. size: " << num_results << endl;
#endif

                nanoflann::KNNResultSet<real> resultSet(num_results);
                resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

                CELL_TREE.index->findNeighbors(resultSet, &testp[0]);
                // uint count = 0;
                for (auto const& cell : ret_indexes)
                {
                    if (CheckCell(cell, cells, testp))
                    {
#ifdef DEBUG
                        // cout << "Moved to a neighbour cell." << endl;
#endif

                        pnp1[ii].cellID = cell;
                        pnp1[ii].cellV = cells.cVel[cell];
                        pnp1[ii].cellP = cells.cP[cell];
                        pnp1[ii].cellRho = cells.cRho[cell];
                        pnp1[ii].ipt_n_failed = 0;
                        pnp1[ii].internal = 0;
                        inside_flag = 1;
                        break;
                        // cout << "Found new cell at " << count << " in list" << endl;
                    }
                }

                if (inside_flag == 0)
                {
                    // cout << "Checking if particle has crossed a boundary." << endl;
                    // If still not found, then the point could be across a boundary.
                    // Check if the ray from the point to the cell centre crosses a boundary face.
                    uint cross = 0;
                    for (auto const& index : ret_indexes)
                    {
                        StateVecD rayp = cells.cCentre[index];
                        for (size_t const& findex : cells.cFaces[index])
                        { /*If its a boundary face, check if the point crosses it*/
                            if (cells.leftright[findex].second < 0)
                            {
                                vector<size_t> const& face = cells.faces[findex];
                                int ints = FALSE;
#if SIMDIM == 3
                                // uint perturb = FALSE;
                                ints = Crossings3D(cells.verts, face, testp, rayp);
#else /*2D line intersection*/
                                ints = get_line_intersection(cells.verts, face, testp, rayp);
#endif
                                if (ints)
                                {
                                    cross = !cross;
                                    if (cells.leftright[findex].second == -1)
                                    {
#pragma omp critical
                                        {
                                            cout << "Particle has crossed an inner boundary!" << endl;
                                        }
                                        pnp1[ii].internal = 1;
                                        break;
                                    }
                                    else if (cells.leftright[findex].second == -2)
                                    {
                                        // cout << "Particle has crossed an outer boundary!" << endl;
                                        localDel.emplace_back(ii);
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    if (cross == 0)
                    {
                        if (pnp1[ii].ipt_n_failed > 10)
                        {
#pragma omp critical
                            {
                                cout << "Particle " << ii - start
                                     << " containing cell not found. Something is wrong." << endl;
                                Write_Containment(svar, ret_indexes, cells, testp);
                                // exit(-1);
                            }
                            // pnp1[ii].b = LOST;  /* Flag as lost to be deleted at the end of this
                            // integration loop */
                            localDel.emplace_back(ii
                            ); /* Just delete particle and allow simulation to continue. */
                               /* Allow user to examine in post whether simulation is bogus */
                        }
                        else
                        {
#pragma omp critical
                            cout << "Couldn't identify which cell particle " << ii - start
                                 << " is inside, but continuing..." << endl;

#pragma omp atomic
                            pnp1[ii].ipt_n_failed++;
                        }
                    }
                }
            }

        } // End of particle loop

#pragma omp for schedule(static) ordered
        for (int ii = 0; ii < omp_get_num_threads(); ii++)
        {
#pragma omp ordered
            toDelete.insert(toDelete.end(), localDel.begin(), localDel.end());
        }
    } /*End of parallel*/

    return toDelete;
}

void Check_Pipe_Outlet(
    Vec_Tree const& CELL_TREE, SIM& svar, MESH const& cells, LIMITS& limits, SPHState& pn,
    SPHState& pnp1, size_t& end
)
{
    vector<size_t> del;
    vector<size_t> ndel(limits.size(), 0);
    for (size_t block = svar.n_bound_blocks; block < svar.n_fluid_blocks + svar.n_bound_blocks; block++)
    {
#pragma omp for nowait
        for (size_t ii = limits[block].index.first; ii < limits[block].index.second; ++ii)
        {
            // Check if pipe vector is defined?
            if (pnp1[ii].b == PIPE)
            { /*Do a check to see if it needs to be given an aero force*/
                if (pnp1[ii].xi.dot(limits[block].aero_norm) > limits[block].aeroconst)
                {
                    pnp1[ii].b = FREE;
                    if (pnp1[ii].lam_nb < svar.air.lam_cutoff && svar.Asource == meshInfl)
                    { /* Retrieve the cell it's in */
                        uint to_del = 0;
                        FirstCell(svar, CELL_TREE, cells, pnp1[ii], to_del);

                        if (to_del)
                        {
                            ndel[block]++;
                            del.emplace_back(ii);
                        }
                    }
                }
            }
        }
    }

    if (!del.empty())
    {
        std::sort(del.begin(), del.end());
        for (vector<size_t>::reverse_iterator itr = del.rbegin(); itr != del.rend(); ++itr)
        {
            pnp1.erase(pnp1.begin() + *itr);
            svar.total_points--;
            svar.fluid_points--;
            end--;
            svar.delete_count++;
        }

        size_t delshift = 0;
        for (size_t block = svar.n_bound_blocks; block < svar.n_fluid_blocks + svar.n_bound_blocks;
             block++)
        {
            // Shift the blocks, which will happen cumulatively
            limits[block].index.first -= delshift;
            delshift += ndel[block];
            limits[block].index.second -= delshift;

            /* Need to shift back vector and buffer vector for correct inlet */
            for (size_t& back : limits[block].back)
            {
                back -= delshift;
            }

            for (vector<size_t>& buffer : limits[block].buffer)
                for (size_t& part : buffer)
                {
                    part -= delshift;
                }
        }
    }
}

/* IMPLICIT PARTICLE TRACKING FUNCTIONS */

/* Function to check if the cell contains the point being tested, and if it needs perturbing */
/* Change function to find the shortest intersection with the infinite planes of each face */
vector<uint> CheckCellFace(
    size_t const& cell, MESH const& cells, size_t const& curr_face, StateVecD const& testp,
    StateVecD const& rayp
)
{
    vector<uint> intersects;
    bool perturb = 0;
    // size_t face_intersect = std::numeric_limits<size_t>::max();
    for (auto const& cFace : cells.cFaces[cell])
    { // Step through cell faces

        vector<size_t> const& face = cells.faces[cFace];
        // if(Crossings2D(cells.verts,face,testp))

        if (Cross_Plane(cells.verts, face, testp, rayp, perturb))
        { /* Intersects a face */
            intersects.emplace_back(static_cast<uint>(cFace));
        }

        // if(perturb)
        // {
        //     perturb = FALSE;

        //     if(Cross_Plane_P(cells.verts,face,testp,rayp,perturb))
        //     {
        //         intersects.emplace_back(static_cast<int>(cFace));
        //     }

        //     if(perturb)
        //     {
        //         cout << "Still can't tell if face was crossed after perturbing" << endl;
        //     }
        //     // cout << "point needs perturbing" << endl;
        //     // pert++;
        //     // return FALSE;
        // }

        /* Assuming convex cells, and point is residing on a face, */
        /* stop at the first interecting face */
    }
    return intersects;
}

real FindFace(SIM const& svar, MESH const& cells, IPTPart const& pn, IPTPart& pnp1)
{
    // Want to find the intersection of the droplet vector with a cell face.

    // Velocity vector = average of start and end
    StateVecD testv = 0.5 * (pnp1.v + pn.v);

    // Test if the velocity is too small or zero, and provide an initial one if so
    if (testv.norm() < 1e-10)
    { // Use the containing cell velocity, since it's bascially guaranteed to follow it
        testv = pnp1.cellV;
    }

    // Test point: starting position
    StateVecD testp = pn.xi - 1e3 * testv.normalized();
    StateVecD rayp = pn.xi + 1e3 * testv.normalized();

    vector<uint> intersects = CheckCellFace(pn.cellID, cells, pn.faceID, testp, rayp);

    // {
    //     // if(pert == 1)
    //     // {
    //     //     testp = pn.xi + vec<real,3>(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));

    //     //     if(!CheckCellFace(pn.cellID,cells,pn.faceID,testp,rayp,pnp1.faceID,pert))
    //     //     {
    //     //         cout << "Still can't identify which cell point is in even after perturbation" <<
    //     endl;
    //     //         pnp1.going = 0;
    //     //     }
    //     // }
    //     // else
    //     // {
    //         cout << "Couldn't tell which face was crossed." << endl;
    //         pnp1.going = 0;
    //     // }
    // }

    /* Found the faces that were crossed. Now find the time spent in the cell */

    int nextface = pn.faceID; /* At least use a reasonable face */
    real mindist = 1e6;
    int hasintersect = 0;
    for (size_t ii = 0; ii < intersects.size(); ++ii)
    {
        // /* Ignore the current face. It's not an option. */
        if (intersects[ii] == pn.faceID)
            continue;

        vector<size_t> const& face = cells.faces[intersects[ii]];

        real dt = 0.0;
        real denom = 0.0;

        RayNormalIntersection(cells, pn.xi, testv, face, pn.cellID, dt, denom);

        /* Check for if the face normal points in or out of the cell */
        // if(cells.leftright[intersects[ii]].first != pn.cellID)
        // {
        //     /* face normal is pointing into the cell, so flip the sign*/
        //     dt = -dt;
        // }

        if (denom > 0)
        {
            /* face normal product with the velocity vector is positive */
            /* so it isn't behind the particle */

            if (dt < mindist)
            {
                nextface = intersects[ii];
                if (dt > MEPSILON)
                {
                    mindist = dt;
                    // cout << "Updated dt: " << dt  << " Crossing face: " << nextface << endl;
                }
                else
                {
                    /* Calculation has resulted in a roundoff suggesting a */
                    /* negative intersection distance, so set distance to zero */
                    /* so particle containment is updated, but not position */
                    // cout << "distance is near 0" << endl;
                    mindist = MEPSILON;
                }
            }

            hasintersect = 1;
        }
    }

    if (hasintersect == 0)
    {
        // cout << "Could not find a positive intersection for some reason..." << endl;
        // cout << "No of intersections with planes: " << intersects.size() << endl;
        // cout << "ParticleID: " << pnp1.part_id << " position: " << pn.xi[0] << "  " << pn.xi[1]
        // << "  " << pn.xi[2] << " velocity n: " << pn.v[0] << "  " << pn.v[1] << "  " << pn.v[2] <<
        //  " velocity np1: " << pnp1.v[0] << "  " << pnp1.v[1] << "  " << pnp1.v[2] << endl;
        pnp1.xi = pn.xi;
        pnp1.v = pn.v;

        pnp1.going = 0;
        return 0.0;
    }
    else
    {
        // auto min_t = std::min_element(dists.begin(),dists.end(),
        // [](std::pair<real,size_t> const& p1, std::pair<real,size_t> const& p2)
        // {return p1.first < p2.first;});

        /* This is the face that has actually been crossed */
        // pnp1.dt = min_t->first;
        // pnp1.faceID = min_t->second;

        // pnp1.dt = pnp1.relax * mindist + (1.0 - pnp1.relax) * pnp1.dt;
        pnp1.faceID = nextface;

        /* Update the face and cell data */
        int newcell = 0;

        if (cells.leftright[nextface].first == static_cast<int>(pn.cellID))
            newcell = cells.leftright[nextface].second;
        else
            newcell = cells.leftright[nextface].first;

        // cout << "first: " << cells.leftright[pnp1.faceID].first << " second: " <<
        //     cells.leftright[pnp1.faceID].second << " cellID: " << pn.cellID << " newcell: "
        //     << newcell << endl;

        if (newcell > -1)
        {
            pnp1.cellID = newcell;
            pnp1.cellV = cells.cVel[newcell];
            pnp1.cellRho = cells.cRho[newcell];
        }
        else
        {
            pnp1.cellID = newcell;
        }
        return mindist; // Return timestep to do something with.
    }
}
