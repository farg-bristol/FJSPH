/*********     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code      *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/


#include "VLM.h"

#include <Eigen/LU>
#include <Eigen/Geometry>
// #include "IO.h"
#include "IOFunctions.h"

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

void VLM::Init(std::string input)
{
    // input.append("VLM.dat");
    std::ifstream fin(input, std::ios::in);

    if(fin.is_open())
    {
        string line;
        while (getline(fin,line))
        {
            Get_Vector(line,"VLM wing dimensions", coords);
            Get_Vector(line,"VLM panel counts",panels);
            Get_Number(line,"VLM angle alpha (degree)",AoA);
            Get_Number(line,"VLM sweep angle (degree)",sweep);
            Get_Number(line,"VLM taper ratio",taper);
            Get_Number(line,"VLM flap start (panels)",flap[0]);
            Get_Number(line,"VLM flap width (panels)",flap[1]);
            Get_Number(line,"VLM flap depth (panels)",flap[2]);
            Get_Number(line,"VLM flap angle (degree)",beta);
            Get_Number(line,"VLM write trajectories (0/1)",write_traj);
            Get_Number(line,"VLM trajectory maximum x coordinate",maxX);
            Get_Number(line,"VLM trajectory maximum iterations",maxIters);
            Get_Number(line,"VLM trajectory target dx",streamDx);
            Get_Number(line,"VLM trajectory minimum dx",streamDx);
        }
    }
    else 
    {
        std::cerr << "Couldn't open " << input << " to read VLM settings. Stopping." << std::endl;
        exit(-1);
    }
    fin.close();

    if(coords[0] == -1 || coords[1] == -1)
    {
        cout << "VLM coordinates not defined." << endl;
        exit(-1);
    }

    if(panels[0] == -1 || panels[1] == -1)
    {
        cout << "VLM panels numbers not defined." << endl;
        exit(-1);
    }

    npanels = 2*panels[0]*panels[1];

    AoA *= M_PI/180.0;
    sweep *= M_PI/180.0;
    beta *= M_PI/180.0;
    
    /*Initialise matrices*/
    panelData.reserve(npanels);
    aInf = Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic>(npanels,npanels);
    gamma = Eigen::Matrix<real,Eigen::Dynamic,1>(npanels);
    RHS = Eigen::Matrix<real,Eigen::Dynamic,1>(npanels);

    MakeMatrix();
}

void VLM::GetGamma(StateVecD inf)
{
    printf("Finding VLM influence...\n");
    Freestream = inf;
    /*Find the influence matrix, and invert it to find gamma*/

    /*Find the influence matrix aInf*/
    for(int i=0; i < npanels; ++i)
    {
        for(int j=0; j < npanels; ++j)
        {
            /*aInf[i,j] = the influence of vortex j on control point i*/
            StateVecD A, B, C;
            A = panelData[i].A;
            B = panelData[i].B;
            C = panelData[j].C;

            aInf(j,i) = FindInfluence(A,B,C).dot(panelData[i].norm);
        }

    RHS(i) = -Freestream.dot(panelData[i].norm);
    }
    /*Find Gamma */
    Eigen::FullPivLU<Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic>> lu(aInf);
    if (lu.isInvertible())
        gamma = lu.solve(RHS);
    else
        std::cerr << "Influence matrix is singular" << std::endl;

    // for (int i= 0; i < gamma.rows(); ++i)
    // {
    // 	std::cout << gamma[i] << std::endl;
    // }
}

const StateVecD VLM::getVelocity(StateVecD const& pos)
{	/*Find velocity for a particle at its position*/
    // Use a Kahan sum to avoid truncation errors in the accumulation
    StateVecD vel = StateVecD::Zero();
    StateVecD c = StateVecD::Zero();// A running compensation for lost low-order bits.

    for (int ii = 0; ii < npanels; ++ii)  
    {
        // StateVecD y = gamma(ii)*FindInfluence(panelData[ii].A,panelData[ii].B,pos) - c; 
        StateVecD y = gamma(ii)*FindInfluenceConditioned(panelData[ii].A,panelData[ii].B,pos) - c;         // c is zero the first time around.
        StateVecD t = vel + y;              // Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - vel) - y;            // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
        vel = t;                      // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
        // Next time around, the lost low part will be added to y in a fresh attempt.
    }

    // for(int i = 0; i<npanels; ++i)
    // 	vel += gamma(i)*FindInfluenceConditioned(panelData[i].A,panelData[i].B,pos);
    
    return vel+Freestream;
}

void VLM::write_VLM_Panels(string &prefix)
{	
    printf("Writing panels...\n");
    string file1 = prefix;
    file1.append("_Panels.dat");
    FILE* fp = fopen(file1.c_str(),"w");
    if(fp == NULL)
    {
        printf("Failed to open VLM Panel file: %s.\nStopping\n",file1.c_str());
        exit(-1);
    }

    
    fprintf(fp,"TITLE=\"VLM Panels\"\n");
    fprintf(fp,"VARIABLES = \"X\", \"Y\", \"Z\"\n");
    for (auto const& p:panelData)
    {
        fprintf(fp,"ZONE\n");
        // fp << "VARLOCATION=([1-3]=NODAL)" << endl;
        fprintf(fp,"%3.7e %3.7e %3.7e\n",p.p1(0), p.p1(1), p.p1(2));
        fprintf(fp,"%3.7e %3.7e %3.7e\n",p.p2(0), p.p2(1), p.p2(2));
        fprintf(fp,"%3.7e %3.7e %3.7e\n",p.p3(0), p.p3(1), p.p3(2));
        fprintf(fp,"%3.7e %3.7e %3.7e\n",p.p4(0), p.p4(1), p.p4(2));
        fprintf(fp,"%3.7e %3.7e %3.7e\n",p.p1(0), p.p1(1), p.p1(2));
    }
    fclose(fp);

    string file2 = prefix;
    file2.append("_Vortices.dat");
    FILE* fq = fopen(file2.c_str(),"w");
    if(fq == NULL)
    {
        printf("Failed to open VLM_Vortices.plt. Attempted path: %s\n",file2.c_str());
        exit(-1);
    }
    fprintf(fq,"TITLE=\"VLM Vortices and Control Points\"\n");
    fprintf(fq,"VARIABLES = \"X\", \"Y\", \"Z\"\n");
    for (auto const& p:panelData)
    {
        fprintf(fq,"ZONE\n");
        fprintf(fq,"%3.7e %3.7e %3.7e\n",p.A(0), p.A(1), p.A(2));
        fprintf(fq,"%3.7e %3.7e %3.7e\n",p.B(0), p.B(1), p.B(2));
        fprintf(fq,"ZONE\n");
        fprintf(fq,"%3.7e %3.7e %3.7e\n",p.C(0), p.C(1), p.C(2));
    }
    fclose(fq);
}

void VLM::MakeMatrix(void)
{
    printf("Contructing VLM Matrix...\n");
    /*Define the steps for each dimension*/
    real dr0 = coords[1]/real(panels[1]);
    real dy0 = (coords[0]/real(panels[0])); // /cos(sweep);
    real dx0 = dr0; // cos(AoA)*dr0;
    real dz0 = 0.0; // sin(AoA)*dr0;
    
    real x0 = 0.0;
    real y0 = 0.0; // cos(AoA)*coords[1];
    real z0 = 0.0; // sin(AoA)*coords[1];
    
    real yend = coords[0]; // *cos(sweep); /* End is now actually the span, and hypotenuse adjusted for the sweep */
    
    /*Find 1/4 panel points (assuming symmetry) */
    /*Find 3/4 Control Points*/
    for(int i = -panels[0]; i < panels[0]; ++i)
    {	/*i = count in z, (spanwise)*/
        /*X coordinates*/
        real y1 = real(i)*dy0;
        real y2 = real(i+1)*dy0;
        real yhalf = (real(i)+0.5)*dy0;
        
        /*x and y values at j = 0 to apply sweep*/
        real xi = tan(sweep)*(fabs(y1)-y0);
        real zi = 0.0; // sin(sweep)*(fabs(y1)-y0);
        real xip1 = tan(sweep)*(fabs(y2)-y0);
        real zip1 = 0.0; // sin(sweep)*(fabs(y2)-y0);

        /*Vertex 1 deltas (i)*/
        real frac1 = (fabs(y1)-y0)/(yend-y0);
        real fac1 = ((1-frac1) + taper*frac1);
        real dx1 = dx0*fac1;
        real dz1 = dz0*fac1;
        
        /*Vertex 2 deltas (i+1)*/
        real frac2 = (fabs(y2)-y0)/(yend-y0);
        real fac2 = ((1-frac2) + taper*frac2);
        real dx2 =  dx0*fac2;
        real dz2 =  dz0*fac2;

        /*Halfway values (i+1/2) for control point*/
        real frac3 = (fabs(yhalf)-y0)/(yend-y0);
        real fac3 = ((1-frac3) + taper*frac3);
        real xhalf = tan(sweep)*(fabs(yhalf)-y0);
        real zhalf = 0.0; // sin(sweep)*(fabs(yhalf)-y0);
        real dzhalf = dz0*fac3;
        real dxhalf = dx0*fac3;

        int jend;
        int doflap = 0; /*0 = no flap, 1 = flap*/

        
        if( -i > flap(0) && -i <= flap(1)+flap(0))
        {	/*Left wing*/
            jend = panels(1)-flap(2);
            doflap = 1;
        }
        else if ( i >= flap(0) && i < flap(1)+flap(0))
        {	/*Right wing*/
            jend = panels(1)-flap(2);
            doflap = 1;
        }
        else
        {	/*No flap*/
            jend = panels(1);
        }

        for (int j=0; j < jend; ++j)
        {	/*j = count in x and y, (chordwise)*/
            /*Panel Verticies (for visual)*/
            StateVecD p1(x0 + xi   + dx1*(real(j))  , y1, z0 + zi   + dz1*(real(j))  );
            StateVecD p2(x0 + xi   + dx1*(real(j+1)), y1, z0 + zi   + dz1*(real(j+1)));
            StateVecD p3(x0 + xip1 + dx2*(real(j+1)), y2, z0 + zip1 + dz2*(real(j+1)));
            StateVecD p4(x0 + xip1 + dx2*(real(j))  , y2, z0 + zip1 + dz2*(real(j))  );

            /*1/4 chord points*/
            StateVecD A(x0 + xi   + dx1*(0.25+real(j)), y1, z0 + zi   + dz1*(0.25+real(j)));
            StateVecD B(x0 + xip1 + dx2*(0.25+real(j)), y2, z0 + zip1 + dz2*(0.25+real(j)));
            /*3/4 control point*/
            StateVecD C(x0 + xhalf + dxhalf*(0.75+real(j)), yhalf, z0 + zhalf + dzhalf*(0.75+real(j)));
            panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
        }

        if(doflap == 1)
        {
            /*Flap delta parameters*/
            real dzflap1 = -(sin(beta))*dr0*fac1;
            real dxflap1 = (cos(beta))*dr0*fac1;
            real dzflap2 = -(sin(beta))*dr0*fac2;
            real dxflap2 = (cos(beta))*dr0*fac2;
            real dzflaphalf = -(sin(beta))*dr0*fac3;
            real dxflaphalf = (cos(beta))*dr0*fac3;

            /*Flap start position*/
            real fac = real(panels(1)-flap(2));
            real zflap1 = z0 + zi   + fac*dz1;
            real xflap1 = x0 + xi   + fac*dx1;
            real zflap2 = z0 + zip1 + fac*dz2;
            real xflap2 = x0 + xip1 + fac*dx2;
            real zflaphalf = z0 + zhalf + fac*dzhalf;
            real xflaphalf = x0 + xhalf + fac*dxhalf;

            for (int j = 0; j < flap(2); ++j)
            {
                /*j = count in x and y, (chordwise)*/
                /*Panel Verticies (for visual)*/
                StateVecD p1(xflap1+dxflap1*(real(j))  , y1, zflap1 + dzflap1*(real(j))  );
                StateVecD p2(xflap1+dxflap1*(real(j+1)), y1, zflap1 + dzflap1*(real(j+1)));
                StateVecD p3(xflap2+dxflap2*(real(j+1)), y2, zflap2 + dzflap2*(real(j+1)));
                StateVecD p4(xflap2+dxflap2*(real(j))  , y2, zflap2 + dzflap2*(real(j))  );

                /*1/4 chord points*/
                StateVecD A(xflap1+dxflap1*(0.25+real(j)), y1, zflap1 + dzflap1*(0.25+real(j)));
                StateVecD B(xflap2+dxflap2*(0.25+real(j)), y2, zflap2 + dzflap2*(0.25+real(j)));
                /*3/4 control point*/
                StateVecD C(xflaphalf+dxflaphalf*(0.75+real(j)), yhalf, zflaphalf + dzflaphalf*(0.75+real(j)));
                panelData.push_back(Panel(A,B,C,p1,p2,p3,p4));
            }
        }
    }

    /*End initialisation*/
}

const StateVecD VLM::FindInfluence(StateVecD const& A, StateVecD const& B, StateVecD const& P)
{	
    StateVecD r0, r1, r2, inf;

    real const static ipi4 = 1.0/(4.0*M_PI);

    /*Bounded Vortex*/
    StateVecD coefAB;

    r0 = B-A;
    r1 = P-A;
    r2 = P-B;

    coefAB = ipi4*((r1.cross(r2))/((r1.cross(r2)).squaredNorm()))*
            (r0.dot(r1.normalized()-r2.normalized()));


    /*Horseshoe vortex from point A*/
    StateVecD coefA;

    inf = A + Freestream;
    r2 = P - A;
    r1 = P - inf;
    r0 = A - inf;

    coefA = ipi4*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
            (1-r0.dot(r2)/(r0.norm()*r2.norm()));

    /*Horseshoe vortex from point B*/
    StateVecD coefB;

    /*Vector B to infinity*/
    inf = B + Freestream;
    r1 = P - B;
    r2 = P - inf;
    r0 = inf - B;

    coefB = ipi4*r0.norm()*((r1.cross(r2))/(r1.cross(r2)).squaredNorm())*
            (r0.dot(r1)/(r1.norm()*r0.norm())+1.0);

    return (coefAB+coefB+coefA);			
}

/* May need some conditioning for points that are nearly colinear */
/* Current test - If any vortex is ill conditioned, ignore the entire horseshoe */
/* Feels more appropriate to do this, but could be excessive. */
const StateVecD VLM::FindInfluenceConditioned(StateVecD const& A, StateVecD const& B, StateVecD const& P)
{	
    StateVecD r0, r1, r2, inf;

    real const static ipi4 = 1.0/(4.0*M_PI);
    StateVecD coefAB(StateVecD::Zero()), coefA(StateVecD::Zero()), coefB(StateVecD::Zero());

    real const cutoff = 1e-12;

    r0 = B-A;
    r1 = P-A;
    r2 = P-B;

    // Bounded Vortex
    // Psi term
    StateVecD numer = r1.cross(r2);
    double denom = numer.squaredNorm();
    // Cut off if denomenator is too small
    if((r1.normalized().cross(r2.normalized())).norm() > cutoff)
        coefAB = ipi4*(numer/denom)*
            (r0.dot(r1.normalized()-r2.normalized()));   //Omega term     
    else
        return StateVecD::Zero();

    // Horseshoe vortex from point A to infinity		
    inf = A + Freestream;
    r2 = P - A;
    r1 = P - inf;
    r0 = A - inf;

    numer = r1.cross(r2);
    denom = numer.squaredNorm();
    // Cut off if denomenator is too small
    if((r1.normalized().cross(r2.normalized())).norm() > cutoff)
        coefA = ipi4*r0.norm()*(numer/denom)*
                (1-r0.dot(r2)/(r0.norm()*r2.norm()));
    else
        return StateVecD::Zero();

    // Horseshoe vortex from point B to infinity
    inf = B + Freestream;
    r1 = P - B;
    r2 = P - inf;
    r0 = inf - B;

    numer = r1.cross(r2);
    denom = numer.squaredNorm();
    // Cut off if denomenator is too small
    if((r1.normalized().cross(r2.normalized())).norm() > cutoff)
        coefB = ipi4*r0.norm()*(numer/denom)*
                (r0.dot(r1)/(r1.norm()*r0.norm())+1.0);
    else
        return StateVecD::Zero();

    return (coefAB+coefB+coefA);			
}

// Create streamlines from the 
void VLM::Plot_Streamlines(std::string& prefix)
{
    printf("Calculating streamlines...\n");
    
    real const gamma = 0.5;
	real const a = 1 - gamma;
	real const b = gamma;
    size_t const max_subits = 15;

    std::vector<vector<StateVecD>> trajectories(npanels);
    
    #pragma omp parallel for
    for(int ii = 0; ii < npanels; ii++)
    {
        std::vector<StateVecD> traj;
        traj.reserve(maxIters);
        // Integrate the velocity
        StateVecD pos = panelData[ii].C;
        traj.emplace_back(pos);
        size_t iter = 0;
        while (pos[0] < maxX && iter < maxIters)
        {
            // Integrate the position using second order RK
            StateVecD veln = getVelocity(pos);
            /* Want to set dt to traverse a maximum distance, 
                and hopefully not lose too much info */
            real dt = streamDx/veln.norm();
            
            /* 2nd order Semi implicit Newmark Beta */
            StateVecD posp1 = pos + dt * veln;
            StateVecD velnp1 = getVelocity(posp1);
            StateVecD posh = posp1;
            real logbase = log10((posp1 - pos).norm());
            real error = 0;
            size_t subit = 0;
            while(error > -7.0 && subit < max_subits)
            {
                posp1 = pos + dt * (a * veln + b * velnp1);
                velnp1 = getVelocity(posp1);
                error = log10((posp1 - posh).norm()) - logbase;
                posh = posp1;
                subit++;
            }
            pos = posp1;

            /* 2nd order RK */
            // StateVecD posh = pos + 0.5 * dt * veln;
            // pos = pos + dt * getVelocity(posh);
            traj.emplace_back(pos);
            iter++;
        }

        trajectories[ii] = traj;
    }

    printf("Writing streamlines...\n");
    /* Need to decide a file */
    string file = prefix + "_Streams.dat";
    FILE* streams = fopen(file.c_str(),"w");
	if(streams == NULL)
	{
		printf("Failed to open ASCII VLM streamtraces file: %s\nStopping.\n", file.c_str());
		exit(-1);
	}

	fprintf(streams, "TITLE=\"VLM Streamtraces\"\n");
	fprintf(streams, "VARIABLES=\"X\",\"Y\",\"Z\"\n");

    for(int ii = 0; ii < npanels; ii++)
    {
        // Write the data
        fprintf(streams, "ZONE T=\"Streamtrace %d\"\n",ii);
        for(StateVecD const& x:trajectories[ii])
        {
            for(uint dim = 0; dim < SIMDIM; ++dim)
                fprintf(streams,"%3.7e ",x[dim]);
            fprintf(streams,"\n");  
        }
    }
    fclose(streams);
}