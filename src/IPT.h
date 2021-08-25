/*********   WCSPH (Weakly Compressible Smoothed Particle Hydrodynamics) Code   *************/
/*********        Created by Jamie MacLeod, University of Bristol               *************/

#ifndef IPT_H
#define IPT_H

#include "Var.h"
#include "Containment.h"
#include "IO.h"

namespace IPT
{
    void Terminate_Particle(SIM& svar, MESH const& cells, vector<IPTPart> const& time_record, size_t const& ii,
        IPTPart& pnp1, vector<SURF>& surface_data )
    {
        pnp1.going = 0;

        // cout << "Particle " << pnp1.partID << " stopped iterating"; 
        if(svar.Asource != 0)
        {
            if(pnp1.cellID < 0)
            {
                /* Find which surface it has impacted  */
                size_t const& face = pnp1.faceID;

                auto const index = std::find_if(cells.smarkers.begin(), cells.smarkers.end(), 
                [face](std::pair<size_t,int> const& p1){return p1.first == face;});

                if(index != cells.smarkers.end())
                {
                    /* Find which surface to add to, and which face of that*/
                    auto const surface = std::find(svar.markers.begin(), svar.markers.end(), index->second);

                    if(surface != svar.markers.end())
                    {
                        size_t const surf = surface - svar.markers.begin();

                        // cout << " after hitting surface: " << svar.bnames[surf] << endl;
                        #pragma omp critical
                        {
                            surface_data[surf].marker_count++;
                            surface_data[surf].marker_pIDs.emplace_back(ii);
                            surface_data[surf].end_pos.emplace_back(pnp1.xi);
                            surface_data[surf].start_pos.emplace_back(time_record[0].xi);
                        }
                        auto const index2 =
                            std::find(surface_data[surf].faceIDs.begin(), surface_data[surf].faceIDs.end(), face);

                        if(index2 != surface_data[surf].faceIDs.end())
                        {
                            size_t fIndex = index2 - surface_data[surf].faceIDs.begin();
                            #pragma omp critical
                            {
                                surface_data[surf].face_count[fIndex]++;
                                surface_data[surf].impacted_face.emplace_back(fIndex);
                            }
                        }
                        else
                        {
                            cout << "Couldn't find which face of the surface the particle impacted" << endl;
                        }
                    }
                    else
                    {
                        cout << "Couldn't find which surface index for the marker." << endl;
                    }
                    
                }
                else
                {
                    cout << "Couldn't find which surface the particle impacted" << endl;
                }

                
            }
            else
            {
                // cout << " prematurely" << endl;
            }
        }
        if(svar.streakout == 1)
        {
            Write_ASCII_Streaks(svar,time_record,pnp1);
        }

        if(svar.cellsout == 1)
        {
            Write_ASCII_Cells(svar,cells,time_record);
        }
    }


    inline real AeroForce(StateVecD const& Vdiff, FLUID const& fvar, AERO const& avar, IPTPart const& pi)
    {
        real const Re = 2.0*pi.faceRho*Vdiff.norm()*pi.d/avar.mug;
        real const Cd = GetCd(Re);


        return Vdiff.norm() * (3.0*Cd*pi.faceRho)/(4.0 * pi.d * fvar.rho0);
    }

    inline void BFD1(FLUID const& fvar, AERO const& avar, StateVecD const& g, IPTPart const& pn, IPTPart& pnp1)
    {
        StateVecD const Vdiff = pnp1.faceV - pnp1.v ;

        real res = AeroForce(Vdiff, fvar, avar, pnp1);
        
        pnp1.acc = res;
        pnp1.v = (pn.v + pnp1.dt*res*pnp1.faceV + pnp1.dt*g)/(1.0 + pnp1.dt*res);
        pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); 
    }

    inline void BFD2(FLUID const& fvar, AERO const& avar, StateVecD const& g, real const& relax, real const& dt, real const& dtm1, 
                            IPTPart const& pnm1, IPTPart const& pn, IPTPart& pnp1)
    {
        StateVecD const Vdiff = pnp1.faceV - pnp1.v ;

        real res = AeroForce(Vdiff, fvar, avar, pnp1);

        pnp1.acc = res;
        /* Second order velocity calculation */
        pnp1.v = ((dt+dtm1)*(dt*dtm1*(relax*res*pnp1.faceV+g) + (dt+dtm1)* pn.v) - dt*dt*pnm1.v)
                    /(dtm1*((2*dt+dtm1)+dt*(dt+dtm1)*res));

                    
        // pnp1.xi = pn.xi-pnp1.dt*(pnp1.v + pn.v)+0.5*pn.dt*(pn.v + pnm1.v); 
        pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); /* Just a first order space integration. */
    }

    void Integrate(SIM& svar, FLUID const& fvar, AERO const& avar, MESH const& cells, size_t const& ii, 
                            IPTPart& pnm1, IPTPart& pn, IPTPart& pnp1, vector<SURF>& marker_data)
    {
        /* Do the first step */
        IPTState time_record;
        time_record.emplace_back(pnp1);
        ofstream partfile;
        if(svar.partout == 1)
        {

        }

        /* Find intersecting face for np1 */
        // pnp1.faceV = pnp1.cellV; /* Start by guessing the face as the cell*/
        // pnp1.faceRho = pnp1.cellRho;
        FindFace(svar, cells, pn, pnp1);
        /* Now c+1 is known, use the central difference face properties*/
        pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
        pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

        // int iter = 0;
        // int min_iter = 3;
        // real error = 1.0;
        real relaxation = 1.0; /* Relaxation factor for 2nd order */
        for(size_t iter = 0; iter < 10; iter++)
        // while(error > -7.0 && iter < min_iter)
        {
            relaxation = (real(iter)+1.0)/6.0;
            if (relaxation > 1.0)
                relaxation = 1.0;
            // cout << "Iteration: " << iter << endl; 
            // vec<real,3> temp = pnp1.xi;

            if (svar.eqOrder == 2)
            {
                real const dt = pnp1.dt;
                real const dtm1 = pnp1.dt;
                BFD2(fvar, avar, svar.grav, relaxation, dt, dtm1, pnm1, pn, pnp1);
            }
            else
                BFD1(fvar, avar, svar.grav, pn, pnp1);

            FindFace(svar, cells, pn, pnp1);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            // error = log10((pnp1.xi-temp).norm());

            // if(iter < svar.max_iter)
            //     iter++;
            // else
            //     break;
        }

        if(svar.partout == 1)
            Write_ASCII_Point(partfile, svar.scale, pnp1);

        if (pnp1.going == 0 || pnp1.v.norm() > 1e3
            || (pnp1.xi - pn.xi).norm() > cells.maxlength )
        {
            pnp1.failed = 1;
            Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
            #pragma omp atomic
                svar.nFailed++;
        }
        else if (pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x)
        {      
            Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
            #pragma omp atomic
                svar.nSuccess++;
        }

        /* Move forward in time */
        // cout << "Time: " << pnp1.t << " dt: " << pnp1.dt  << " x: " <<
        // pnp1.xi[0] << " y: " << pnp1.xi[1] << " z: " << pnp1.xi[2] << 
        // " old cell: " << pn.cellID << " new cell: " << pnp1.cellID << endl;

        /* March forward in time */
        if(svar.eqOrder == 2)
            pnm1 = pn;

        pn = pnp1;
        pnp1.t += pnp1.dt;
        

        /* Only need the information if streaks or cells are output */
        if(svar.cellsout == 1 || svar.streakout == 1)
            time_record.emplace_back(pnp1);

        while(pnp1.going != 0)
        {
        
            /* Find intersecting face for np1 */
            // pnp1.faceV = pnp1.cellV; /* Start by guessing the face as the cell*/
            // pnp1.faceRho = pnp1.cellRho;
            
            if(pnp1.d < 1.0e-5)
            {   /* If particle is small, use cell velocity to start */
                pnp1.v = pnp1.cellV;
            }
            else if (pnp1.d < 100.0e-6)
            {   /* If particle is heavier, use sliding scale up to 100 microns */
                pnp1.v = (1.0 - (pnp1.d-1e-05)/(9e-5)) * pnp1.cellV + 
                (pnp1.d-1e-05)/(9e-5) * pnp1.v;
            }


            FindFace(svar, cells, pn, pnp1);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            // int iter = 0;
            // int min_iter = 3;
            // real error = 1.0;
            real relaxation = 1.0;
            for(size_t iter = 0; iter < 10; iter++)
            // while(error > -7.0 && iter < min_iter)
            {
                relaxation = (real(iter)+1.0)/6.0;
                if (relaxation > 1.0)
                    relaxation = 1.0;

                // cout << "Iteration: " << iter << endl; 
                // vec<real,3> temp = pnp1.xi;
                
                if(svar.eqOrder == 2)
                {
                    real const dt = pnp1.dt;
                    real dtm1 = pnp1.dt;
                    if(iter > 0)
                        dtm1 = pn.dt;
                    BFD2(fvar, avar, svar.grav, relaxation, dt, dtm1, pnm1, pn, pnp1);
                }
                else
                    BFD1(fvar, avar, svar.grav, pn, pnp1);


                FindFace(svar, cells, pn, pnp1);
                /* Now c+1 is known, use the central difference face properties*/
                pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
                pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

                // error = log10((pnp1.xi-temp).norm());

                // if(iter < svar.max_iter)
                //     iter++;
                // else
                //     break;
            }

            if(svar.partout == 1)
                Write_ASCII_Point(partfile, svar.scale, pnp1);

            if ( pnp1.going == 0 || pnp1.v.norm() > 1e3
                    || (pnp1.xi - pn.xi).norm() > cells.maxlength )
            {
                pnp1.failed = 1;
                Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
                #pragma omp atomic
                    svar.nFailed++;
                break;
            }
            else if ( pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x)
            {      
                Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
                #pragma omp atomic
                    svar.nSuccess++;
                break;
            }

            // /* Move forward in time */
            // cout << "Time: " << pnp1.t << " dt: " << pnp1.dt  << " x: " <<
            // pnp1.xi[0] << " y: " << pnp1.xi[1] << " z: " << pnp1.xi[2] << 
            // " old cell: " << pn.cellID << " new cell: " << pnp1.cellID << endl;

            /* March forward in time */
            if(svar.eqOrder == 2)
                pnm1 = pn;
            
            pn = pnp1;
            pnp1.t += pnp1.dt;
            
            /* Only need the information if streaks or cells are output */
            if(svar.cellsout == 1 || svar.streakout == 1)
                time_record.emplace_back(pnp1);

        }
        
    }
}
#endif