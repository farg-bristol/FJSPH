/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#include "Geometry.h"
#include "Kernel.h"
#include "Neighbours.h"
#include "Var.h"

#include <Eigen/Geometry>

#define PERTURB(i, j) pow(MEPSILON, pow(2, i* SIMDIM - j))

// #define X 0
// #define Y 1

/*Surface detection as described by Marrone, Colagrossi, Le Touze, Graziani -
 * (2010)*/
void Detect_Surface(
    SIM& svar, FLUID const& fvar, AERO const& avar, size_t const& start, size_t const& end,
    OUTL const& outlist, MESH const& cells, VLM const& vortex, SPHState& pnp1
)
{

    vector<StateVecD> norms(end, StateVecD::Zero());

#pragma omp parallel default(shared) // shared(pnp1)
    {
        real const h = 1.33 * svar.Pstep;

#pragma omp for
        for (size_t ii = start; ii < end; ++ii)
        {

            StateVecD norm = StateVecD::Zero();

            SPHPart& pi = pnp1[ii];
            if (pi.b < PIPE)
            {
                pi.surf = 0;
                pi.woccl = 1;
            }
            else if (pi.lam_nb < 0.2)
            { /*If less that 0.2, then its a surface
                 particle by default*/
                pi.surf = 1;
            }
            else if (pi.lam_nb < 0.75)
            { /*Particle could be a surface particle.
                 Perform a test*/

                // Create point T
                StateVecD const xi = pi.xi;
                StateVecD const pointT = xi + h * pi.norm.normalized();

                uint surf = 1; /*Begin by assuming a surface, and proving otherwise*/
                for (neighbour_index const& jj : outlist[ii])
                {
                    if (jj.first == ii || pnp1[jj.first].b == GHOST)
                        continue;

                    SPHPart const& pj = pnp1[jj.first];

#if SIMDIM == 3
                    StateVecD const Rij = pi.xi - pj.xi;
#endif

                    StateVecD const x_jT = pj.xi - pointT;
                    real const r = sqrt(jj.second);

                    if (r >= sqrt(2) * h)
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
                        StateVecD tau(pi.norm(1), -pi.norm(0));
                        if ((abs(pi.norm.normalized().dot(x_jT)) + abs(tau.normalized().dot(x_jT))) < h)
                        {
                            surf = 0;
                            break;
                        }
#else
                        if (acos(pi.norm.normalized().dot(-Rij / r)) < M_PI / 4.0)
                        {
                            surf = 0;
                            break;
                        }
#endif
                    }

                    // #if defined(CSF) || defined(HESF)
                    // norm -= volj * (dp.colour[jj.first] -
                    // dp.colour[ii])*GradK(Rij,r,fvar.H,fvar.correc); #endif
                }

                pi.surf = surf;
            } /* end if lam < 0.75 */
            else
            { /*If its eigenvalue is high, then by default it cannot be a
                 surface*/
                pi.surf = 0;
            }

            // if(pi.lam_nb <  avar.cutoff || pi.surf == 1)
            // {
            /* Calculate normal using eigenvalue gradients */

            if (pi.lam > 0.7)
            {
                for (neighbour_index const& jj : outlist[ii])
                {
                    SPHPart const& pj = pnp1[jj.first];
                    if (jj.first == ii || pj.b == GHOST)
                        continue;
                    StateVecD const Rij = pi.xi - pj.xi;
                    real const r = sqrt(jj.second);
                    real const volj = pj.m / pj.rho;

                    norm += volj * (pj.lam - pi.lam) * (GradK(Rij, r, fvar.H, fvar.correc));
                }
            }
            else
            {
                for (neighbour_index const& jj : outlist[ii])
                {
                    SPHPart const& pj = pnp1[jj.first];
                    if (jj.first == ii || pj.b == GHOST)
                        continue;
                    StateVecD const Rij = pi.xi - pj.xi;
                    real const r = sqrt(jj.second);
                    real const volj = pj.m / pj.rho;

                    norm += volj * pj.lam * (GradK(Rij, r, fvar.H, fvar.correc));
                }
            }

            norm = pi.L * norm;
            if (norm.norm() > 0.1 * pi.lam / fvar.H)
                norms[ii] = norm.normalized();
            // }
        }

#pragma omp for
        for (size_t ii = start; ii < end; ++ii)
        { /* Now calculate curvature using the normals */
            real woccl_ = 0.0;
            StateVecD Vdiff = StateVecD::Zero();

            SPHPart& pi = pnp1[ii];
            if (pi.cellID != -1 && pi.b == FREE && (avar.acase == Gissler))
            {
                if (svar.Asource == meshInfl)
                {
                    Vdiff = (pi.cellV) - pi.v;
                }
#if SIMDIM == 3
                else if (svar.Asource == VLMInfl)
                {
                    StateVecD Vel = vortex.getVelocity(pi.xi);
                    Vdiff = Vel - pi.v;
                }
#endif
                else
                {
                    Vdiff = avar.vInf - pi.v;
                }
            }

            real curve = 0.0;
            for (neighbour_index const& jj : outlist[ii])
            {
                SPHPart const& pj = pnp1[jj.first];
                if (jj.first == ii || pj.b == GHOST)
                    continue;

                StateVecD const Rij = pj.xi - pi.xi;
                real const r = sqrt(jj.second);

                /* Needed for induced pressure model and CSF surface tension */
                real volj = pj.m / pj.rho;

// if(pi.lam > 0.7)
// {
//     curve += volj * (pi.L*(norms[jj.first] - norms[ii])).dot
//                 (GradK(Rij,r,fvar.H,fvar.correc));
//     // curve += volj * (dp.L[ii]*(dp.norm[jj.first].normalized() -
//     pi.norm.normalized())).dot
//     //         (GradK(Rij,r,fvar.H,fvar.correc));
// }
// else
// {
//     #if SIMDIM == 2
//     if(/* norms[ii].dot(norms[jj.first]) < 0.5 && */ norms[ii].norm() > 0 &&
//     norms[jj.first].norm() > 0) #else if(norms[ii].dot(norms[jj.first]) >
//     -0.333 && norms[ii].norm() > 0 && norms[jj.first].norm() > 0) #endif
//     // if(norms[ii].norm() > 0 && norms[jj.first].norm() > 0)
//     {
//         // curve += volj * (dp.L[ii]*(dp.norm[jj.first].normalized() -
//         pi.norm.normalized())).dot
//         //     (GradK(Rij,r,fvar.H,fvar.correc));
//         curve += volj * (pi.L*(norms[jj.first] - norms[ii])).dot
//                 (GradK(Rij,r,fvar.H,fvar.correc));
//         // curve += volj * ((norms[jj.first] - norms[ii])).dot
//         //         (pi.L*GradK(Rij,r,fvar.H,fvar.correc));
//     }
// }

/* Just do the check for all particles. Getting discontinuity in curvature */
#if SIMDIM == 2
                if (/* norms[ii].dot(norms[jj.first]) < 0.5 && */ norms[ii].norm() > 0 &&
                    norms[jj.first].norm() > 0)
#else
                if (/* norms[ii].dot(norms[jj.first]) > -0.333 && */ norms[ii].norm() > 0 &&
                    norms[jj.first].norm() > 0)
#endif
                // if(norms[ii].norm() > 0 && norms[jj.first].norm() > 0)
                {
                    // curve += volj * (dp.L[ii]*(dp.norm[jj.first].normalized() -
                    // pi.norm.normalized())).dot
                    //     (GradK(Rij,r,fvar.H,fvar.correc));
                    curve +=
                        volj *
                        (pi.L * (norms[jj.first] - norms[ii])).dot(GradK(Rij, r, fvar.H, fvar.correc));
                    // curve += volj * ((norms[jj.first] - norms[ii])).dot
                    //         (pi.L*GradK(Rij,r,fvar.H,fvar.correc));
                }

                // curve += volj * ((dp.norm[jj.first].normalized() -
                // pi.norm.normalized())).dot
                //             (/* dp.L[ii]* */GradK(Rij,r,fvar.H,fvar.correc));

                // curve -= volj * ((norms[jj.first] - norms[ii])).dot
                //         (/* dp.L[ii]* */GradK(Rij,r,fvar.H,fvar.correc));

                /*Occlusion for Gissler Aero Method*/
                if (pi.b == FREE && avar.acase == Gissler)
                {
                    real const frac = -Rij.dot(Vdiff) / (Vdiff.norm() * r);

                    if (frac > woccl_)
                    {
                        woccl_ = frac;
                    }
                }
            }

            if (pi.lam_nb < avar.cutoff)
            {
                pi.woccl = std::max(0.0, std::min(woccl_, 1.0));
            }
            else
            {
                pi.woccl = 1.0;
            }

            pi.norm = norms[ii];
            pi.curve = curve /* /correc */;
            pi.s = svar.dx * curve /* / correc */;
            pi.pDist = pi.lam;
        } /* end particle loop */

#if defined(ALE) || defined(TIC)
#pragma omp for nowait
        for (size_t ii = start; ii < end; ++ii)
        {
            pnp1[ii].surfzone = 0;
            for (neighbour_index const& jj : outlist[ii])
            {
                if (pnp1[jj.first].surf == 1)
                {
                    pnp1[ii].surfzone = 1;
                    break;
                }
            }
        }
#endif

    } /* end parallel section */
}

real get_n_full(real const& dx, real const& H)
{
    // Create a lattice of points from -2H to 2H
    std::vector<StateVecD> vec;
    for (real x = -2.0 * (H + dx); x <= 2.0 * (H + dx); x += dx)
    {
        for (real y = -2.0 * (H + dx); y <= 2.0 * (H + dx); y += dx)
        {
#if SIMDIM == 3
            for (real z = -2.0 * (H + dx); z <= 2.0 * (H + dx); z += dx)
                vec.emplace_back(StateVecD(x, y, z));
#else
            vec.emplace_back(StateVecD(x, y));
#endif
        }
    }

    // Create a tree of the vector?
    const real search_radius = 4.0 * H * H;
    Vec_Tree tree(SIMDIM, vec, 10);
    tree.index->buildIndex();
    StateVecD test = StateVecD::Zero();
    std::vector<neighbour_index> matches =
        radius_search(tree, test, search_radius); /* Nearest Neighbour Search*/

    return real(matches.size());
}

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
// intersect the intersection point may be stored in the floats i_x and i_y.
bool get_line_intersection(
    vector<StateVecD> const& verts, vector<size_t> const& edge, StateVecD const& p1,
    StateVecD const& cellC
)
{
    const StateVecD& e1 = verts[edge[0]];
    const StateVecD& e2 = verts[edge[1]];
    StateVecD s, r;
    s = cellC - p1;
    r = e2 - e1;

    // Find the denominator of the calculation
    real denom = (-r(0) * s(1) + s(0) * r(1));

    // If the value of this is nearly 0, then
    if (denom < MEPSILON)
    { // Lines are colinear
        return 0;
    }

    real u, t;
    u = (-s(1) * (p1(0) - e1(0)) + s(0) * (p1(1) - e1(1))) / denom;
    t = (r(0) * (p1(1) - e1(1)) - r(1) * (p1(0) - e1(0))) / denom;

    if (u > 0 && u < 1 && t > 0 && t < 1)
    { // Collision detected
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
    int yflag0, yflag1, inside_flag;
    real ty, tx;
    StateVecD vtx0, vtx1;

    tx = point[0];
    ty = point[1];

    inside_flag = 0;

    vtx0 = verts[edge[0]];
    vtx1 = verts[edge[1]];
    /* Move to the next pair of vertices, retaining info as possible. */
    yflag0 = (vtx0[1] >= ty);
    yflag1 = (vtx1[1] >= ty);

    /* Check if endpoints straddle (are on opposite sides) of X axis
     * (i.e. the Y's differ); if so, +X ray could intersect this edge.
     * Credit to Joseph Samosky to try dropping
     * the "both left or both right" part of my code.
     */
    if (yflag0 != yflag1)
    {
        /* Check intersection of pgon segment with +X ray.
         * Note if >= point's X; if so, the ray hits it.
         * The division operation is avoided for the ">=" test by checking
         * the sign of the first vertex wrto the test point; idea inspired
         * by Joseph Samosky's and Mark Haigh-Hutchinson's different
         * polygon inclusion tests.
         */
        if (((vtx1[1] - ty) * (vtx1[0] - vtx0[0]) >= (vtx1[0] - tx) * (vtx1[1] - vtx0[1])) == yflag1)
        {
            inside_flag = !inside_flag;
        }

        /* For convex cells, further optimisation can be done: */
        /* A ray can only pass through a maximum of two faces.*/
        /* If this is second edge hit, then done testing. */
    }
    return (inside_flag);
}

/* Check for intersection with the infinite plane (I.e. just do two volumes) */
/* Used for implicit particle tracking */
int Cross_Plane(
    vector<StateVecD> const& verts, vector<size_t> const& face, StateVecD const& point,
    StateVecD const& point2, bool& perturb
)
{ /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996 */
    StateVecD const testp = point;
    StateVecD const rayp = point2;
    real vol = testp[0] * (verts[face[0]][1] - verts[face[1]][1]) -
               testp[1] * (verts[face[0]][0] - verts[face[1]][0]) +
               (verts[face[0]][0] * verts[face[1]][1] - verts[face[1]][0] * verts[face[0]][1]);
    int flag1, flag2;

    // if(LessThanREError(vol1))
    // {   // Perturb the test point so that it doesn't go into roundoff error
    //     perturb = TRUE;
    //     return 0;
    // }

    flag1 = (vol < 0.0);

    vol = rayp[0] * (verts[face[0]][1] - verts[face[1]][1]) -
          rayp[1] * (verts[face[0]][0] - verts[face[1]][0]) +
          (verts[face[0]][0] * verts[face[1]][1] - verts[face[1]][0] * verts[face[0]][1]);
    // ray point is very far so I don't see this falling into roundoff error
    // if(LessThanREError(vol1))
    // {
    //     perturb = TRUE;
    //     return 0;
    // }

    flag2 = (vol < 0.0);

    /* If signs of the volumes alternate, */
    /* then the ray intersects the infinite plane*/
    if (flag1 != flag2)
    {
        return 1;
    }
    return 0;
}

/* Check for the distance from the face the point lies. It is already assumed
 * that the ray crosses */
/* Used for implicit particle tracking */
void RayNormalIntersection(
    MESH const& cells, StateVecD const& rayOrigin, StateVecD const& rayVector,
    vector<size_t> const& face, int const& cellID, real& dt, real& denom
)
{
    /* find the normal vector of the face */
    StateVecD norm(
        cells.verts[face[0]][1] - cells.verts[face[1]][1],
        cells.verts[face[1]][0] - cells.verts[face[0]][0]
    );

    /* Use the face centre as the second part for the direction */
    StateVecD face_c = cells.verts[face[0]] + cells.verts[face[1]];
    face_c /= 2.0;

    /* Check for normal orientation */
    StateVecD celldir = face_c - cells.cCentre[cellID];

    if (norm.dot(celldir) < 0)
    {
        /* Normal points inwards to the cell , so flip it*/
        norm = -1.0 * norm;
    }

    StateVecD temp_p = StateVecD::Zero();

    /* Find the most distant point from the current point */
    for (size_t ii = 0; ii < face.size(); ii++)
    {
        if ((cells.verts[face[ii]] - rayOrigin).squaredNorm() > temp_p.squaredNorm())
            temp_p = (cells.verts[face[ii]] - rayOrigin);
    }

    /* Find numerator */
    denom = rayVector.dot(norm);
    dt = temp_p.dot(norm) / denom;
}

#endif

#if SIMDIM == 3
int Crossings3D(
    vector<StateVecD> const& verts, vector<size_t> const& face, StateVecD const& point,
    StateVecD const& point2 /* , uint& perturb */
)
{ /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996, &
    Robust Adaptive Floating-Point Geometric Predicates
    Michael Aftosmis, Cart3D Software*/
    /*https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/bool_intersection.html*/
    StateVecD const testp = point;
    StateVecD const rayp = point2;
    StateP1MatD vol1;
    int flag1, flag2;
    vol1 << testp(0), testp(1), testp(2), 1.0, verts[face[0]](0), verts[face[0]](1), verts[face[0]](2),
        1.0, verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0, verts[face[2]](0),
        verts[face[2]](1), verts[face[2]](2), 1.0;

    // if(LessThanREError(vol1))
    // {   // Perturb the test point so that it doesn't go into roundoff error
    //     perturb = TRUE;
    //     return 0;
    // }

    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0) << rayp(0), rayp(1), rayp(2), 1.0;

    // ray point is very far so I don't see this falling into roundoff error
    // if(LessThanREError(vol1))
    // {
    //     vol1.row(0) << rayp(0)+PERTURB(0,1), rayp(1)+PERTURB(0,2),
    //     rayp(2)+PERTURB(0,3),1.0;
    // }

    flag2 = (vol1.determinant() < 0.0);

    /*If signs of the volumes alternate, then the points lie either side of the
     * plane*/
    /*Now check if the line drawn by the two points intersects inside the bounds
     * of the triangle plane*/
    if (flag1 != flag2)
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
        vol.row(3) << rayp(0), rayp(1), rayp(2), 1.0;

        // if(LessThanREError(vol))
        // {   /*Perturb the test point, since all volume calculations need to be
        // done with the new location*/
        //     perturb = TRUE;
        //     return 0;
        // }

        flag3 = (vol.determinant() < 0.0);

        /*Check for each face, if the signs of all the tets are the same.*/
        for (size_t ii = 1; ii < 3; ++ii)
        { /*Change the face vertices used*/
            // cout << face.size() << "  " << ii << endl;
            vtx0 = vtx1;
            vtx1 = verts[face[ii]];

            vol.row(1) << vtx0(0), vtx0(1), vtx0(2), 1.0;
            vol.row(2) << vtx1(0), vtx1(1), vtx1(2), 1.0;

            // if(LessThanREError(vol))
            // {
            //     perturb = TRUE;
            //     return 0;
            // }

            flag4 = (vol.determinant() < 0.0);

            /*If the sign of the tet is different, this face isn't intersected.*/
            if (flag4 != flag3)
                return 0;
        }
        return 1;
    }
    return 0;
}

/* Check for intersection with the infinite plane (I.e. just do two volumes) */
/* Used for implicit particle tracking */
int Cross_Plane(
    vector<StateVecD> const& verts, vector<size_t> const& face, StateVecD const& point,
    StateVecD const& point2, bool& perturb
)
{ /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996 */
    StateVecD const testp = point;
    StateVecD const rayp = point2;
    StateP1MatD vol1;
    int flag1, flag2;
    vol1 << testp[0], testp[1], testp[2], 1.0, verts[face[0]][0], verts[face[0]][1], verts[face[0]][2],
        1.0, verts[face[1]][0], verts[face[1]][1], verts[face[1]][2], 1.0, verts[face[2]][0],
        verts[face[2]][1], verts[face[2]][2], 1.0;

    // if(LessThanREError(vol1))
    // {   // Perturb the test point so that it doesn't go into roundoff error
    //     perturb = TRUE;
    //     return 0;
    // }

    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0) << rayp[0], rayp[1], rayp[2], 1.0;

    // ray point is very far so I don't see this falling into roundoff error
    // if(LessThanREError(vol1))
    // {
    //     perturb = TRUE;
    //     return 0;
    // }

    flag2 = (vol1.determinant() < 0.0);

    /* If signs of the volumes alternate, */
    /* then the ray intersects the infinite plane*/
    if (flag1 != flag2)
    {
        return 1;
    }
    return 0;
}

/* Check for intersection with the infinite plane (I.e. just do two volumes) */
/* Used for implicit particle tracking */
int Cross_Plane_P(
    vector<StateVecD> const& verts, vector<size_t> const& face, StateVecD const& point,
    StateVecD const& point2, bool& perturb
)
{ /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996 */
    StateVecD const testp = point;
    StateVecD const rayp = point2;
    StateP1MatD vol1;
    int flag1, flag2;
    vol1 << testp[0], testp[1], testp[2], 1.0, verts[face[0]][0] + PERTURB(0, 0),
        verts[face[0]][1] + PERTURB(0, 1), verts[face[0]][2] + PERTURB(0, 2), 1.0,
        verts[face[1]][0] + PERTURB(1, 0), verts[face[1]][1] + PERTURB(1, 1),
        verts[face[1]][2] + PERTURB(1, 2), 1.0, verts[face[2]][0] + PERTURB(2, 0),
        verts[face[2]][1] + PERTURB(2, 1), verts[face[2]][2] + PERTURB(2, 2), 1.0;

    if (LessThanREError(vol1))
    { // Perturb the test point so that it doesn't go
      // into roundoff error
        perturb = TRUE;
        return 0;
    }

    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0) << rayp[0], rayp[1], rayp[2], 1.0;

    // ray point is very far so I don't see this falling into roundoff error
    if (LessThanREError(vol1))
    {
        perturb = TRUE;
        return 0;
    }

    flag2 = (vol1.determinant() < 0.0);

    /* If signs of the volumes alternate, */
    /* then the ray intersects the infinite plane*/
    if (flag1 != flag2)
    {
        return 1;
    }
    return 0;
}

bool MollerTrumbore(
    StateVecD const& v0, StateVecD const& v1, StateVecD const& v2, StateVecD const& rayOrigin,
    StateVecD const& rayVector, real& dt, real& denom
)
{
    StateVecD const e1 = v1 - v0;
    StateVecD const e2 = v2 - v1;

    StateVecD h, s, q;
    real a, f, u, v;

    h = rayVector.cross(e2);
    a = e1.dot(h);
    if (a < MEPSILON && a > -MEPSILON)
    { // Ray is parallel to the triangle
        dt = MEPSILON;
        denom = MEPSILON;
        return false;
    }

    f = 1.0 / a;
    s = rayOrigin - v0;
    u = f * s.dot(h);

    if (u < 0.0 || u > 1.0)
    { // Does not intersect within the bounds of the triangle.
        dt = 9999999;
        denom = -1;
        return false;
    }

    q = s.cross(e1);
    v = f * rayVector.dot(q);
    if (v < 0.0 || u + v > 1.0)
    { // Does not intersect within the bounds of the triangle.
        dt = 9999999;
        denom = -1;
        return false;
    }

    // Does intersect, so find distance (dt) along the ray to intersection
    dt = f * e2.dot(q);
    if (dt > MEPSILON)
        denom = 1.0;
    else
        denom = -1.0;
    return true;
}

// Ray tracing intersection using the Moller Trubmore algorithm
void RayNormalIntersection(
    MESH const& cells, StateVecD const& rayOrigin, StateVecD const& rayVector,
    vector<size_t> const& face, int const& cellID, real& dt, real& denom
)
{

    if (face.size() == 3)
    { // Just need to do the algorithm, and not worry about
      // whether it does intersect
        MollerTrumbore(
            cells.verts[face[0]], cells.verts[face[1]], cells.verts[face[2]], rayOrigin, rayVector, dt,
            denom
        );
    }
    else
    { // Need to check which side of the face is intersected, so need to
      // know the result.
        if (!MollerTrumbore(
                cells.verts[face[0]], cells.verts[face[1]], cells.verts[face[2]], rayOrigin, rayVector,
                dt, denom
            ))
            MollerTrumbore(
                cells.verts[face[0]], cells.verts[face[3]], cells.verts[face[2]], rayOrigin, rayVector,
                dt, denom
            );
    }
}

/* Check for the distance from the face the point lies. It is already assumed
 * that the ray crosses */
/* Used for implicit particle tracking */
// void RayNormalIntersection(MESH const& cells, StateVecD const& rayOrigin,
// StateVecD const& rayVector,
//                            vector<size_t> const& face, int const& cellID,
//                            real& dt, real& denom)
// {
//     /* find the normal vector of the face */
//     StateVecD norm;
//     if(face.size() == 3)
//     {   /* Cross two edges of the triangle */
//         norm =
//         (cells.verts[face[1]]-cells.verts[face[0]]).cross(cells.verts[face[2]]-cells.verts[face[0]]);
//     }
//     else
//     {   /* Cross the two diagonals of the square */
//         norm =
//         (cells.verts[face[3]]-cells.verts[face[1]]).cross(cells.verts[face[2]]-cells.verts[face[0]]);
//     }

//     /* Check for normal orientation */
//     StateVecD celldir;

//     /* Use the face centre as the second part for the direction */
//     StateVecD face_c = StateVecD::Zero();
//     for(size_t ii = 0; ii < face.size(); ii++)
//     {
//         face_c += cells.verts[face[ii]];
//     }

//     face_c /= real(face.size());

//     celldir = face_c - cells.cCentre[cellID];

//     if(norm.dot(celldir) < 0)
//     {
//         /* Normal points inwards to the cell , so flip it*/
//         norm = -1.0 * norm;
//     }

//     StateVecD temp_p = StateVecD::Zero();

//     /* Find the most distant point from the current point */
//     for(size_t ii = 0; ii < face.size(); ii++)
//     {
//         if( (cells.verts[face[ii]] - rayOrigin).squaredNorm() >
//         temp_p.squaredNorm())
//             temp_p = (cells.verts[face[ii]] - rayOrigin);
//     }

//     /* Find numerator */
//     denom = rayVector.dot(norm);
//     dt = temp_p.dot(norm)/denom;
// }

#endif

// Need cell elements to check size, and do the tet case.
// Need cell centres for every other element
real Cell_Volume(
    vector<StateVecD> const& verts, vector<vector<size_t>> const& faces, vector<size_t> const& elems,
    vector<size_t> const& cell, StateVecD const& cCentre
)
{
#if SIMDIM == 3
    if (elems.size() == 4)
    {
        // Cell is a tetrahedron. Volume is trivial
        StateP1MatD vol;

        vol << verts[elems[0]](0), verts[elems[0]](1), verts[elems[0]](2), 1.0, verts[elems[1]](0),
            verts[elems[1]](1), verts[elems[1]](2), 1.0, verts[elems[2]](0), verts[elems[2]](1),
            verts[elems[2]](2), 1.0, verts[elems[3]](0), verts[elems[3]](1), verts[elems[3]](2), 1.0;

        return abs(vol.determinant()) / 6.0;
    }

    // Otherwise, volume is not so trivial
    // Form tetrahedrons with the cell centre.
    real sum = 0.0;

    for (auto const& faceID : cell)
    {
        StateP1MatD vol;

        const vector<size_t> face = faces[faceID];

        vol << cCentre(0), cCentre(1), cCentre(2), 1.0, verts[face[0]](0), verts[face[0]](1),
            verts[face[0]](2), 1.0, verts[face[1]](0), verts[face[1]](1), verts[face[1]](2), 1.0,
            verts[face[2]](0), verts[face[2]](1), verts[face[2]](2), 1.0;

        sum += abs(vol.determinant()) / 6.0;
    }

    return sum;
#else

    if (elems.size() == 3)
    {
        StateP1MatD vol;

        vol << verts[elems[0]](0), verts[elems[0]](1), 1.0, verts[elems[1]](0), verts[elems[1]](1), 1.0,
            verts[elems[2]](0), verts[elems[2]](1), 1.0;

        return abs(vol.determinant()) / 2;
    }

    real sum = 0.0;
    for (auto const& faceID : cell)
    {
        StateP1MatD vol;

        vector<size_t> const& face = faces[faceID];

        vol << cCentre(0), cCentre(1), 1.0, verts[face[0]](0), verts[face[0]](1), 1.0, verts[face[1]](0),
            verts[face[1]](1), 1.0;

        sum += abs(vol.determinant()) / 2.0;
    }

    return sum;

#endif
}

void Make_Cell(FLUID const& fvar, AERO const& avar, MESH& cells)
{
#if SIMDIM == 3

    cells.verts.emplace_back(StateVecD(-0.5, -0.5, -0.5));
    cells.verts.emplace_back(StateVecD(-0.5, -0.5, 0.5));
    cells.verts.emplace_back(StateVecD(-0.5, 0.5, 0.5));
    cells.verts.emplace_back(StateVecD(-0.5, 0.5, -0.5));
    cells.verts.emplace_back(StateVecD(0.5, -0.5, -0.5));
    cells.verts.emplace_back(StateVecD(0.5, -0.5, 0.5));
    cells.verts.emplace_back(StateVecD(0.5, 0.5, 0.5));
    cells.verts.emplace_back(StateVecD(0.5, 0.5, -0.5));
    cells.verts.emplace_back(StateVecD(1.5, -0.5, -0.5));
    cells.verts.emplace_back(StateVecD(1.5, -0.5, 0.5));
    cells.verts.emplace_back(StateVecD(1.5, 0.5, 0.5));
    cells.verts.emplace_back(StateVecD(1.5, 0.5, -0.5));

    // Cell 1
    cells.faces.emplace_back(vector<size_t>{0, 1, 2});
    cells.faces.emplace_back(vector<size_t>{0, 2, 3});
    cells.faces.emplace_back(vector<size_t>{0, 3, 7});
    cells.faces.emplace_back(vector<size_t>{0, 7, 4});
    cells.faces.emplace_back(vector<size_t>{0, 5, 1});
    cells.faces.emplace_back(vector<size_t>{0, 4, 5});
    cells.faces.emplace_back(vector<size_t>{1, 5, 6});
    cells.faces.emplace_back(vector<size_t>{1, 6, 2});
    cells.faces.emplace_back(vector<size_t>{3, 2, 6});
    cells.faces.emplace_back(vector<size_t>{3, 6, 7});

    // Interface
    cells.faces.emplace_back(vector<size_t>{4, 6, 5});
    cells.faces.emplace_back(vector<size_t>{4, 7, 6});

    // Cell 2
    cells.faces.emplace_back(vector<size_t>{4, 7, 11});
    cells.faces.emplace_back(vector<size_t>{4, 11, 8});
    cells.faces.emplace_back(vector<size_t>{4, 9, 5});
    cells.faces.emplace_back(vector<size_t>{4, 8, 9});
    cells.faces.emplace_back(vector<size_t>{5, 9, 10});
    cells.faces.emplace_back(vector<size_t>{5, 10, 6});
    cells.faces.emplace_back(vector<size_t>{8, 10, 9});
    cells.faces.emplace_back(vector<size_t>{8, 11, 10});
    cells.faces.emplace_back(vector<size_t>{7, 6, 10});
    cells.faces.emplace_back(vector<size_t>{7, 10, 11});

    // cells.elems.emplace_back(vector<size_t>{0,1,2,3,4,5,6,7});
    // cells.elems.emplace_back(vector<size_t>{4,5,6,7,8,9,10,11});

    cells.cFaces.emplace_back(vector<size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11});
    cells.cFaces.emplace_back(vector<size_t>{10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21});

    cells.cCentre.emplace_back(StateVecD(0.0, 0.0, 0.0));
    cells.cCentre.emplace_back(StateVecD(1.0, 0.0, 0.0));

    for (uint ii = 0; ii < 10; ++ii)
        cells.leftright.emplace_back(std::pair<int, int>(0, -2));

    cells.leftright.emplace_back(std::pair<int, int>(0, 1));
    cells.leftright.emplace_back(std::pair<int, int>(0, 1));

    for (uint ii = 0; ii < 10; ++ii)
        cells.leftright.emplace_back(std::pair<int, int>(1, -2));

#else
    cells.verts.emplace_back(StateVecD(0.0, 0.0));
    cells.verts.emplace_back(StateVecD(0.0, 1.0));
    cells.verts.emplace_back(StateVecD(1.0, 1.0));
    cells.verts.emplace_back(StateVecD(1.0, 0.0));

    cells.faces.emplace_back(vector<size_t>{0, 1});
    cells.faces.emplace_back(vector<size_t>{1, 2});
    cells.faces.emplace_back(vector<size_t>{2, 3});
    cells.faces.emplace_back(vector<size_t>{3, 0});

    // cells.elems.emplace_back(vector<size_t>{0,1,2,3});

    cells.cFaces.emplace_back(vector<size_t>{0, 1, 2, 3});

    cells.cCentre.emplace_back(StateVecD(0.5, 0.5));

    for (uint ii = 0; ii < 4; ++ii)
        cells.leftright.emplace_back(std::pair<int, int>(0, -2));
#endif

    // Generate the solution vectors
    cells.cVel.emplace_back(avar.vInf);
    cells.cP.emplace_back(avar.pRef);
    // cells.SPHRho.emplace_back(fvar.rho0 * pow((fvar.gasPress/fvar.B +
    // 1),1/fvar.gam));

    cells.fNum.emplace_back(0);
    cells.fMass.emplace_back(fvar.simM);
    cells.vFn.emplace_back(StateVecD::Zero());
    cells.vFnp1.emplace_back(StateVecD::Zero());

    cells.cRho.emplace_back(avar.rhog);

    cells.cPertn.emplace_back(StateVecD::Zero());
    cells.cPertnp1.emplace_back(StateVecD::Zero());

    // Second cell
    cells.cVel.emplace_back(avar.vInf);
    cells.cP.emplace_back(avar.pRef);
    // cells.SPHRho.emplace_back(fvar.rho0 * pow((fvar.gasPress/fvar.B +
    // 1),1/fvar.gam));

    cells.fNum.emplace_back(0);
    cells.fMass.emplace_back(fvar.simM);
    cells.vFn.emplace_back(StateVecD::Zero());
    cells.vFnp1.emplace_back(StateVecD::Zero());

    cells.cRho.emplace_back(avar.rhog);

    cells.cPertn.emplace_back(StateVecD::Zero());
    cells.cPertnp1.emplace_back(StateVecD::Zero());
}
