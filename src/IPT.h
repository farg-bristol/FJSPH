/******     FJSPH (Fuel Jettison Smoothed Particles Hydrodynamics) Code ***********/
/******          Created by Jamie MacLeod, University of Bristol        ***********/

#ifndef IPT_H
#define IPT_H

#include "Var.h"

/* Particle tracking functions */
namespace IPT
{
    void Init_IPT_Files(SIM& svar);

    void Write_Data(SIM& svar, MESH& cells, vector<IPTState>& time_record);

    void Integrate(
        SIM& svar, MESH const& cells, size_t const& ii, IPTPart& pnm1, IPTPart& pn, IPTPart& pnp1,
        vector<SURF>& marker_data, vector<IPTState>& iptdata
    );
} // namespace IPT
#endif