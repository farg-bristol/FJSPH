
#include "arc.h"

#include <Eigen/Dense>

void check_arc_input(shape_block& bound, real& globalspacing, int& fault)
{
    int has_config = 0;
    int arc_defined = 0;

    #if SIMDIM == 2
    // These checks are only valid for 2D
    if(check_vector(bound.centre))
    {   /* Default to these values first */
        if(bound.arc_start >= 0 && bound.arc_end >= 0 && bound.radius > 0)
        {
            has_config = 1;
        }
        else if(bound.arc_start >= 0 && bound.ni > 0 && bound.radius > 0)
        {
            has_config = 1;
        }
        else if(bound.arc_start >= 0 && bound.arclength != default_val && bound.radius > 0)
        {
            has_config = 1;
            arc_defined = 1;
        }
    }
    #endif
    
    if(!has_config)
    {
        if(check_vector(bound.centre) && check_vector(bound.start) && check_vector(bound.end))
        {   // No arclength, but just the end acquired
            has_config = 1;
            #if SIMDIM == 2
            get_arclength_centrepoint(bound.start,bound.end,bound.centre,
                                    bound.radius,bound.arc_start,bound.arc_end);
            #else
            get_arclength_centrepoint(bound.start,bound.end,bound.centre,
                        bound.right,bound.radius,bound.arc_start,bound.arc_end);
            #endif
        }
    }


    #if SIMDIM == 3
    if(!has_config)
    {
        // check arch geometry first. Want to use an arclength first
        // Midpoint is required to define the plane appropriately, but not for the circle
        if(check_vector(bound.centre) && check_vector(bound.start) && check_vector(bound.mid)
            && bound.arclength != default_val)
        {
            has_config = 1;
            arc_defined = 1;
        }
    }
    #else
    if(!has_config)
    {
        if(check_vector(bound.start) && check_vector(bound.end) && check_vector(bound.mid))
        {   // No arclength, but just the end acquired
            has_config = 1;
            #if SIMDIM == 2
            get_arclength_midpoint(bound.start,bound.end,bound.mid,
                        bound.centre,bound.radius,bound.arc_start,bound.arc_end);
            #else
            get_arclength_midpoint(bound.start,bound.end,bound.mid,bound.centre,
                        bound.right,bound.radius,bound.arc_start,bound.arc_end);
            #endif
        }
    }
    #endif

    if(!has_config)
    {   // Check if the arclength is defined
        if(check_vector(bound.centre) && check_vector(bound.start) && bound.arclength != default_val)
        {  
            #if SIMDIM == 3
            if(check_vector(bound.right))
            {   // Need to check if normal is defined to know the plane
                has_config = 1;
                arc_defined = 1;
            }
            #else
            has_config = 1;
            arc_defined = 1;
            #endif
        }
    }

    if(!has_config)
    {
        std::cout << "ERROR: Block \"" << bound.name << 
            "\" arc geometry has not been sufficiently defined. Stopping." << std::endl;
        fault = 1;
    }

    if(bound.dx < 0)
    {
        if(bound.ni < 0)
        {
            std::cout << "WARNING: Block globalspacing has not been defined. Using global globalspacing." << std::endl;
            bound.dx = globalspacing;
        }
        else
        {
            /* Use ni as a number along the diameter, so define globalspacing off that */
            real di = (2.0*bound.radius)/real(bound.ni);
            bound.dx = di;
        }
    }
    

    if(bound.thickness < 0)
    {
        if(bound.nk < 0)
        {
            std::cout << "ERROR: Block \"" << bound.name << "\" arc thickness has not been correctly defined. Stopping." << std::endl;
            fault = 1;
        }
    }
    else if (bound.nk < 0)
    {   /* Use a particle count over the real thickness value */
        if(bound.hcpl == 1)
            bound.nk = static_cast<int>(ceil(bound.thickness / globalspacing / sqrt(3.0) * 2.0));
        else
            bound.nk = static_cast<int>(ceil(bound.thickness / globalspacing));
    }

    // Check that start coordinate exists
    // #if SIMDIM == 2
    // // Need to find the start vector first if it's not defined
    // if(!check_vector(bound.start))
    // {
    //     bound.start(bound.radius*cos(bound.arc_start), bound.radius*sin(bound.arc_start));
    //     bound.start += bound.centre;
    // }
    // #endif

    real dtheta = globalspacing/bound.radius;
    if(arc_defined)
    {

        #if SIMDIM == 3
            get_arc_end(bound.start,bound.centre,bound.right,bound.arclength,bound.end,bound.radius);
        // #else
            
            // get_arc_end(bound.start,bound.centre,bound.arclength,bound.end,bound.radius,bound.arc_start);
        #endif
        // bound.arclength *= M_PI/180.0;
  
    }
    else
    {
        /* Estimate how many points on the Block */     
        if(bound.arc_end < 0 && bound.ni > 0)
        {   /* Can define arc using only starting angle and number of particles, to generate an end */
            bound.arc_end = 180.0/M_PI * (bound.arc_start*M_PI/180 + real(bound.ni) * dtheta);
        }
        bound.arclength = (bound.arc_end-bound.arc_start);
    }

    int ni = static_cast<int>(ceil((fabs(bound.arclength)*M_PI/180)/dtheta));
    ni = ni > 1 ? ni : 1;
    bound.ni = ni;
    // Need to add straights in to the count
    
    int smax = bound.sstraight > 0 ? ceil(bound.sstraight/globalspacing) : 0;
    int emax = bound.estraight > 0 ? ceil(bound.estraight/globalspacing) : 0;

    bound.npts = (ni + smax + emax) * bound.nk;
    #if SIMDIM == 3
    if(bound.length < 0)
    {
        if(bound.nj < 0)
        {
            std::cout << "ERROR: Block \"" << bound.name << "\" arch length has not been correctly defined. Stopping." << std::endl;
            fault = 1;
        }
    }
    else if (bound.nj < 0)
    {
        if(bound.hcpl == 1)
            bound.nj = static_cast<int>(ceil((bound.length) / globalspacing / sqrt(6.0) * 3.0));
        else
            bound.nj = static_cast<int>(ceil(bound.length / globalspacing));
        bound.nj = bound.nj > 1 ? bound.nj : 1;
    }
    bound.npts *= bound.nj;
    #endif
}

#if SIMDIM == 2
inline std::vector<StateVecD> make_arc(int const& ni, int const& nk, real const& globalspacing, 
                                real const& theta0, real const& dtheta, real const& radius, 
                                StateVecD const& centre, real const& slength, real const& elength,
                                 int const& hcpl)
{
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    // Find the starting point of the straight
    StateVecD svec(cos(theta0),sin(theta0));
    StateVecD snormal(svec[1],-svec[0]);
    // StateVecD startp = svec*radius - slength * snormal;
    int smax = ceil(slength/globalspacing);

    StateVecD evec(cos(theta0+static_cast<real>(ni-1)*dtheta),
                sin(theta0+static_cast<real>(ni-1)*dtheta));
    StateVecD enormal(-evec[1],evec[0]);
    int emax = ceil(elength/globalspacing);
    // StateVecD endp = evec - elength * enormal;

    std::vector<StateVecD> points;
    for(int jj = 0; jj < nk; jj++)
    {
        // Add start and ending straights if they're defined
        if(slength > 0)
        {
            for(int ii = smax; ii > 0; ii--)
            {
                StateVecD newPoint = svec*(radius - real(jj)*globalspacing) + snormal*real(ii)*globalspacing;
                newPoint += StateVecD(unif(re),unif(re)) + centre;
                points.push_back(newPoint);
            }
        }

        for(int ii = 0; ii < ni; ii++)
        {
            real theta = theta0 + real(ii)*dtheta;
            real r = radius;
            // newPoint = 0.5 * StateVecD(real(2 * i + (j % 2)), sqrt(3.0) * real(j));
            if(hcpl == 1)
            {   /* Shift by half a globalspacing if jj is even */
                theta += 0.5*dtheta*(jj%2);
                r -= 0.5*sqrt(3)*real(jj)*globalspacing;
            }
            else
                r -= real(jj)*globalspacing;
            
            StateVecD newPoint(cos(theta) * r, sin(theta) * r);
            newPoint += StateVecD(unif(re),unif(re));
            newPoint += centre;
            points.push_back(newPoint);
        }

        if(elength > 0)
        {
            for(int ii = 1; ii <= emax; ii++)
            {
                StateVecD newPoint = evec*(radius - real(jj)*globalspacing) + enormal*real(ii)*globalspacing;
                newPoint += StateVecD(unif(re),unif(re)) + centre;
                points.push_back(newPoint);
            }
        }
    }
    return points;
}

void get_arc_end(StateVecD const& start, StateVecD const& centre, real const& arclength, 
                StateVecD& end, real& radius, real& theta0)
{
    StateVecD u = start - centre;
    radius = u.norm();

    theta0 = atan2(u[0],u[1]);
    
    real alen = arclength * M_PI/180.0; 
    real theta1 = theta0 + alen;
    StateVecD d2(cos(theta1),sin(alen));
    d2 = d2 * radius + centre;
    end = d2;
}


void get_arclength_centrepoint(StateVecD const& start, StateVecD const& end, StateVecD const& centre,
                        real& radius, real& theta0, real& theta1)
{
    StateVecD d0 = start - centre;
    StateVecD d1 = end - centre;
    real r = d0.squaredNorm();

    if(fabs(d1.squaredNorm() - r) > 0.001)
    {
        std::cout << "Arc points are not correctly defined. The ending radius differs from the starting radius" << std::endl;
        exit(1);
    }

    radius = sqrt(r);
    
    real t0 = atan2(d0[1],d0[0]);
    real t1 = atan2(d1[1],d1[0]);

    // Want them stored in degrees
    // theta0 = std::min(t0,t1) * 180/M_PI;
    // theta1 = std::max(t0,t1) * 180/M_PI;
    theta0 = t0 * 180/M_PI;
    theta1 = t1 * 180/M_PI;

    /* Check arc lengths */
    // real arclength = theta1 - theta0;
}

void get_arclength_midpoint(StateVecD const& start, StateVecD const& end, StateVecD const& midpoint, 
                StateVecD& centre, real& radius, real& theta0, real& theta1)
{
    StateVecD r1(start[0] * start[0], start[1] * start[1]);
    StateVecD r2(midpoint[0] * midpoint[0], midpoint[1] * midpoint[1]);
    StateVecD r3(end[0] * end[0], end[1] * end[1]);

    /* Find the arc centre */
    Eigen::Matrix2d m;
    m << (start[0] - end[0]), (start[1] - end[1]), (start[0] - midpoint[0]), (start[1] - midpoint[1]);
    Eigen::Vector2d v(r1[0] - r3[0] + r1[1] - r3[1], r1[0] - r2[0] + r1[1] - r2[1]);
    m *= 2.0;
    Eigen::Vector2d result = m.colPivHouseholderQr().solve(v);

    centre = StateVecD(result[0],result[1]);

    get_arclength_centrepoint(start,end,centre,radius,theta0,theta1);
}

std::vector<StateVecD> create_arc_segment(shape_block const& block, real const& globalspacing)
{
    /* Convert to radians */
    real theta0_ = block.arc_start * M_PI/180.0; 
    real dtheta = globalspacing / block.radius;
    // int ni = static_cast<int>(ceil(arclength/dtheta));
    
    return make_arc(block.ni,block.nk,globalspacing,theta0_,dtheta,block.radius,
                block.centre,block.sstraight,block.estraight,block.hcpl);
}

#endif

#if SIMDIM == 3


inline std::vector<StateVecD> make_arch(StateVecD const& start, StateVecD const& centre, StateVecD const& right,
                                real const& slength, real const& elength, int const& nrad,
                                int const& nthick, int const& nlong, real const& globalspacing, 
                                real const& dtheta, real const& radius, int const& hcpl)
{
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    // Create local coordinate system (u,v,w)
    StateVecD u = (start - centre).normalized();
    StateVecD w = right.normalized();
    StateVecD v = u.cross(w).normalized();

    int smax = ceil(slength/globalspacing);
    int emax = ceil(elength/globalspacing);

    StateVecD evec = cos(real(nrad-1)*dtheta) * u + sin(real(nrad-1)*dtheta) * v;
    StateVecD enorm = evec.cross(w);
    
    std::vector<StateVecD> points;
    for(int kk = 0; kk < nlong; kk++)
    {
        for(int jj = 0; jj < nthick; jj++)
        {
            real r = radius;
            real l = real(kk)*globalspacing;
            real doffset = 0;
            if(hcpl == 1)
            {
                r += 1.0 / 3.0 * sqrt(6.0) * real(jj) * globalspacing;
                l = 0.5*sqrt(3)*(real(kk) + real(jj%2)/3.0)*globalspacing;
                doffset = 0.5*dtheta*((jj+kk)%2);
            }
            else
                r += real(jj)*globalspacing;

            if(slength > 0)
            {
                for(int ii = smax; ii > 0; ii--)
                {
                    real dist = real(ii)*globalspacing + doffset;
                    StateVecD newPoint = u * r - v*dist + l*w;
                    newPoint += StateVecD(unif(re),unif(re),unif(re)) + centre;
                    points.push_back(newPoint);
                }
            }

            for(int ii = 0; ii < nrad; ii++)
            {
                real theta = static_cast<real>(ii)*dtheta;
                real l = real(kk)*globalspacing;
                if(hcpl == 1)
                {   /* Shift by half a globalspacing if jj is even */
                    theta += doffset;
                }

                real a = cos(theta)*r; real b = sin(theta)*r;
                StateVecD newPoint = a * u + b * v + l*w;
                newPoint += StateVecD(unif(re),unif(re),unif(re));
                newPoint += centre;
                points.push_back(newPoint);
            }

            if(elength>0)
            {
                for(int ii = 1; ii <= emax; ii++)
                {
                    real dist = real(ii)*globalspacing + doffset;
                    StateVecD newPoint = evec * r + enorm*dist + l*w;
                    newPoint += StateVecD(unif(re),unif(re),unif(re)) + centre;
                    points.push_back(newPoint);
                }
            }
        }
    }
    return points;
}

void get_arc_end(StateVecD const& start, StateVecD const& centre, StateVecD const& normal, 
                 real const& arclength, StateVecD& end, real& radius)
{
    StateVecD u = start - centre;

    radius = u.norm();
          u /= radius;

    StateVecD w = normal.normalized();

    // Check that the start sits on the plane
    if(fabs(normal.dot(start) - normal.dot(centre)) > 1e-4)
    {
        std::cout << "Points do not exist upon the plane defined by the provided normal" << std::endl;
        exit(-1);
    }

    StateVecD v = u.cross(w).normalized();// Should be unit length but just to be sure
    // Have local plane coordinates now get the end point
    real alen = arclength * M_PI/180.0; 
    StateVecD d2 = radius*(cos(alen) * u + sin(alen) * v) + centre;
    end = d2;
}

void get_arclength_centrepoint(StateVecD const& start, StateVecD const& end, StateVecD const& centre, StateVecD& right,
                        real& radius, real& theta0, real& theta1)
{
    StateVecD d1 = start - centre;
    StateVecD d2 = end - centre;
    real r = d1.squaredNorm();

    if(fabs(d2.squaredNorm() - r) > 0.001)
    {
        std::cout << "Arc points are not correctly defined. The ending radius differs from the starting radius" << std::endl;
        exit(1);
    }

    radius = sqrt(r);
    // May need a check for colinearity? 
    // If it's a semicircle, normal will be poorly defined, unless I permit normal to be defined
    // Need to use projected plane coordinates (a,b)
    d1 = d1.normalized();
    d2 = d2.normalized();
    StateVecD w = d2.cross(d1).normalized();
    StateVecD v = w.cross(d1).normalized(); // the 'y' axis of the projected plane
    // So start should be (1,0) in projected coordinates
    // But end will be non trivial (end = a.U + b.V)
    // Check which axis has the best resolution, since the system is over constrained
    real maxx = fabs(d1[0]) + fabs(v[0]);
    real maxy = fabs(d1[1]) + fabs(v[1]);
    real maxz = fabs(d1[2]) + fabs(v[2]);


    Eigen::Matrix2d m; Eigen::Vector2d result;
    Eigen::Vector2d v_;
    if(maxx < maxy && maxx < maxz)
    {   //x is the worst dimension, use y and z
        m << d1[1], v[1], d1[2], v[2];
        v_(v[1],v[2]);
    }
    else if(maxy < maxx && maxy < maxz)
    {   // y is the worst dimension, use x and z
        m << d1[0], v[0], d1[2], v[2];
        v_(v[0],v[2]);
    }
    else
    {   // z is the worst dimension, use x and y
        m << d1[0], v[0], d1[1], v[1];
        v_(v[0],v[1]);
    }

    result = m.colPivHouseholderQr().solve(v_);
    real a = result[0]; real b = result[1];

    real t0 = atan2(1,0);
    real t1 = atan2(a,b);

    // Want them stored in degrees
    // theta0 = std::min(t0,t1) * 180/M_PI;
    // theta1 = std::max(t0,t1) * 180/M_PI;
    theta0 = t0 * 180/M_PI;
    theta1 = t1 * 180/M_PI;

    right = w; //Also retain the normal vector travelling along the archway.
    /* Check arc lengths */
    // real arclength = theta1 - theta0;
}

void get_arclength_midpoint(StateVecD const& start, StateVecD const& end, StateVecD const& midpoint, 
                StateVecD& centre, StateVecD& right, real& radius, real& theta0, real& theta1)
{
    StateVecD r1 = midpoint - start;
    StateVecD r3 = end - start;

    // Create local coordinate system (u,v,w)
    StateVecD u = r1.normalized();
    StateVecD w = r3.cross(r1).normalized();
    StateVecD v = w.cross(u);

    // Find the arc centre
    real bx = r1.dot(u);
    Eigen::Matrix<real,2,1> c(r3.dot(u),r3.dot(v));
    // Will exist upon the line x = bx/2, thus need to find y coordinate (bx/2,h).
    // Distance from point c must be the same as from the centre, thus
    real hdist = (pow(c[0]-bx/2.0,2.0) + c[1]*c[1] - bx*bx/4.0)/(2.0*c[1]);

    centre = start + bx*u/2.0 + hdist*v;

    get_arclength_centrepoint(start,end,centre,right,radius,theta0,theta1);
}

std::vector<StateVecD> create_arc_segment(shape_block const& block, real const& globalspacing)
{
    /* Convert to radians */
    // real theta0_ = block.arc_start * M_PI/180.0; 
    // real theta1_ = block.arc_end * M_PI/180.0;

    // real arclength = theta1_-theta0_;
    // if(arclength < 0)
    //     arclength += 2.0*M_PI;
    real arclength = block.arclength * M_PI/180.0;
    real dtheta = globalspacing / block.radius;
    int ni = static_cast<int>(ceil(arclength/dtheta));
    
    return make_arch(block.start,block.centre,block.right,block.sstraight,block.estraight,ni,block.nk,
                block.nj,globalspacing,dtheta,block.radius,block.hcpl);
}

std::vector<StateVecD> create_fibre_arch(shape_block& var, real const globalspacing)
{
    /* Create a plane with i and j indexing */
    std::vector<StateVecD> points = create_arc_segment(var,globalspacing);
    
    /* Add data about the index sizes */
    // var.ni = ni-1;
    // var.nj = nj-1;
    var.ni--;
    var.nj--;

    return points;
}

#endif // End 3D code