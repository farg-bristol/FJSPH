#ifndef DROPLET_H
#define DROPLET_H
/* Library to perform drag assessment for a range of particle resolutions */
#include "Var.h"
#include "Shifting.h"
#include "Neighbours.h"
#include "Aero.h"
#include "Add.h"


/*Sphere-Plate interpolation method - Gissler et al (2017)*/
inline void const Droplet_GisslerForce(AERO const& avar, StateVecD const& Vdiff, StateVecD const& norm, 
						real const& rho, real const& press, real const& mass, real const& lam, real const& woccl,
                        StateVecD& totAcc, StateVecD& sphereAcc, StateVecD& plateAcc)
{
	// real const nfull = avar.nfull;
	real const Re = 2*rho*Vdiff.norm()*avar.L/avar.mug;
	
	real const frac2 = std::min(1.5 * lam, 1.0);
	real const frac1 = (1.0 - frac2);

 	real const Cds  = GetCd(Re);

	real Cdl, Adrop;
	
	#if SIMDIM == 3 
		if(avar.useDef)
		{
			real ymax = Vdiff.squaredNorm()*avar.ycoef;
			if (ymax > 1.0)
				ymax = 1.0;
			Cdl = Cds*(1+2.632*ymax);
			Adrop = M_PI*pow((avar.L + avar.Cb*avar.L*ymax),2);
		}
		else
		{
			Cdl = Cds;
			Adrop = avar.aSphere;
		}
	#else
		if(avar.useDef)
		{
			real ymax = Vdiff.squaredNorm()*avar.ycoef;
			if (ymax > 1.0)
				ymax = 1.0;
			 Cdl = Cds*(1+2.632*ymax);
			Adrop = avar.aSphere + 2*(avar.Cb*avar.L*ymax) /** pow(avar.L,1)*/;
		}
		else
		{
			Cdl = Cds;
			Adrop = avar.aSphere;
		}
	#endif

	real const Cdi = frac1*Cdl + /* 0.5* */ /*1.37**/frac2;
	real const Aunocc = (frac1*Adrop + frac2*avar.aPlate);

	real const Ai = (1.0 - woccl)*Aunocc;

	// cout << "Fractions: " << frac1 << "  " << frac2 << endl;
	// cout << "Areas: " << Ai << "  " << Adrop << "  " << avar.aPlate << endl;
	// cout << "Re: " << Re << "  Cds: " << Cdi << "  "  << Cdl << "  " << Cds  << endl;
	// cout << "F: " << F(0) << "  " << F(1) << endl << endl;

	StateVecD const q = 0.5 * Vdiff.norm() * Vdiff * avar.isossqr;
	totAcc = q * avar.gamma * press * Cdi * Ai / mass;
			 
    sphereAcc = (1.0 - woccl) * frac1 * q *
			 avar.gamma * press * Cdl * Adrop / mass;
    
    plateAcc = (1.0 - woccl) * frac2 * q *
			 avar.gamma * press * avar.aPlate / mass; /* Cd of plate assumed to be 1 */

}

inline void const Droplet_InducedPressure(AERO const& avar, StateVecD const& Vdiff,
        StateVecD const& norm, real const& Pbasei, real const& lam, SPHPart const& pi,
        StateVecD& totAcc, StateVecD& sphereAcc, StateVecD& plateAcc)
{
	real theta = abs(acos(-norm.normalized().dot(Vdiff.normalized())));
		
	real Cp_p = 0.0;
	real Cp_s = 0.0;
	real Cp_tot = 0.0;

	/*Spherical Cp*/
	if(theta < 2.4455)
	{
		Cp_s = 1.0 - (2.25) * pow(sin(theta),2.0);
	}
	else
	{
		Cp_s = 0.075;
	}

	/*Plate Cp*/
	if(theta < 1.570797)
	{
		Cp_p = cos(theta);
	}
	else if (theta < 1.9918)
	{
		Cp_p = -pow(cos(6.0*theta+0.5*M_PI),1.5);
	}
	else if (theta < 2.0838)
	{
		Cp_p = 5.5836*theta-11.5601;
	}
	else
	{
		Cp_p = 0.075;
	}

	/* Overall Cp (Ideal to machine learn at some point) */
	if(pi.curve > 200.0)
	{
		if(theta < 1.570796)
			Cp_tot = 0.5*cos(2.0*theta)+0.5;
		else
			Cp_tot = 0.0;
	}
	else if (pi.curve > 0.0)
	{
		real const frac = (pi.curve)/(200.0);
		Cp_tot = frac + (1.0-frac)*Cp_p;
	}
	else if (pi.curve > -200.0)
	{
		real const frac = (pi.curve + 200.0)/(200.0);
		Cp_tot = frac*Cp_p + (1.0-frac)*Cp_s;
	}
	else
	{
		Cp_tot = Cp_s;
	}

	/* Compressible dynamic pressure */
	real const Plocali = 0.5*Vdiff.squaredNorm() * avar.isossqr * avar.gamma * pi.cellP * Cp_tot;
	// cout << Plocalj << endl;

	real const Pi = (Plocali/* +Pbasei */);
	
	// #pragma omp critical
	// cout << Cp_s << "  " << Cp_p << "  " << Cp_tot << "  " << theta << "  " << pi.curve << "  " << Pi <<  endl;

	real const Re = 2*pi.cellRho*Vdiff.norm()*avar.L/avar.mug;
	real const Cdi  = GetCd(Re);

	/* Pure droplet force */
	StateVecD const acc_drop = 0.5*Vdiff*Vdiff.norm() * avar.isossqr * avar.gamma * pi.cellP
					* avar.aSphere * Cdi/pi.m;
	// aeroD = -Pi * avar.aPlate * norm[ii].normalized();

	/* Induce pressure force */
	StateVecD const acc_kern =  -/*0.5 **/ Pi * avar.aPlate * norm.normalized()/pi.m;
	// StateVecD const F_kern = -(Pi/(pi.rho)) * norm;

	/*Next, consider a skin friction force acting parallel to the surface*/
	real const Vnorm = Vdiff.dot(norm.normalized());
	StateVecD const Vpar = Vdiff - Vnorm*norm.normalized();

	real const Re_par = 2.0*pi.cellRho*Vpar.norm()*avar.L/avar.mug;
	real const Cf = 0.027/pow(Re_par+1e-6,1.0/7.0); /*Prandtl seventh power law for turbulent BL*/
	// real const Cf = 0;

	// StateVecD const acc_skin = StateVecD::Zero();
	// StateVecD const acc_skin = 0.5*avar.rhog* Vpar.norm() * Cf * avar.aPlate * Vpar/pi.m;
	StateVecD const acc_skin = 0.5 * Vpar.norm() * 	Vpar / (avar.sos*avar.sos) *
			 avar.gamma * pi.cellP * Cf * avar.aPlate / pi.m;

	real const frac1 = std::min(1.5 * lam, 1.0);
	// real const frac2 = std::min(exp(pi.curve*0.001+200),1.0);


	totAcc = (frac1 * /*frac2**/ (acc_kern+acc_skin) + (1.0-frac1) * acc_drop);
    sphereAcc = (1.0-frac1)*acc_drop;
    plateAcc = (frac1 * /*frac2**/ (acc_kern/* +acc_skin */));
}


void Droplet_Drag_Sweep(SIM& svar, FLUID& fvar, AERO& avar)
{
	cout << "Droplet drag sweet selected. Beginning sweep..." << endl;
	vector<int> new_res;
	/* Check if velocities or Reynolds have been defined */
	if(svar.Reynolds.empty())
	{
		if(svar.velocities.empty())
		{
			cout << "ERROR: No Reynolds or velocities have been defined for the sweep. Stopping. " << endl;
			exit(-1);
		}
		else
		{
			for(real const& vel:svar.velocities)
			{
				real const Re = avar.rhog*vel*svar.diam/avar.mug;
				svar.Reynolds.emplace_back(Re);
			}
		}
	}

	if(svar.nacross.empty())
	{
		if(svar.diameters.empty())
		{
			cout << "ERROR: No resolutions have been defined for the sweep. Stopping. " << endl;
			exit(-1);
		}
		else
		{
			for(real const& diam : svar.diameters)
			{
				if(svar.diam/diam < 1.5)
				{
					new_res.emplace_back(1);
					svar.nacross.emplace_back(1);
				}
				else
				{
					real radius = 0.5*svar.diam;
					int nrad = ceil(radius/diam);
					svar.nacross.emplace_back(2*nrad+1);
					new_res.emplace_back(2*nrad + 1);
				}
			}
		}
	}
	else
	{
		for(int const& res : svar.nacross)
		{
			if(res == 1)
			{
				new_res.emplace_back(1);
			}
			else
			{
				real radius = 0.5*svar.diam;
				real dx = radius/real(0.5*(res-1));
				int nrad = ceil(radius/dx);
				new_res.emplace_back(2*nrad + 1);
			}
		}
	}

    /* Level 1 = resolution, level 2 = Reynolds, level 3 = different drag values */
    vector<vector<vector<real>>> drop_data(
		svar.nacross.size(),vector<vector<real>>(svar.Reynolds.size(),vector<real>(6,0.0)));

    MESH cells;

    if(cells.cCentre.size() == 0)
        cells.cCentre.emplace_back(StateVecD::Zero());

    if(cells.fNum.size() == 0)
    {
        cells.fNum.emplace_back();
        cells.fMass.emplace_back();
        cells.vFn.emplace_back();
        cells.vFnp1.emplace_back();
        cells.cVol.emplace_back();
        cells.cMass.emplace_back();
        cells.cRho.emplace_back();
        cells.cPertn.emplace_back();
        cells.cPertnp1.emplace_back();
    }

    for(size_t res_i = 0; res_i < svar.nacross.size(); ++res_i)
    {
		size_t const& res =  svar.nacross[res_i];

        drop_data.emplace_back();
        /* Create the droplet */
        SPHState pn;
        DELTAP dp;
        OUTL outlist;
        
		svar.simPts = 0;
		svar.totPts = 0;
        
        if(res == 1)
		{	/*spacing is too close to full size to use standard adjustment*/
			svar.nrad = 1;
			svar.dx = svar.diam;
			svar.Pstep = svar.diam;
			fvar.simM = fvar.rho0*pow(svar.Pstep,SIMDIM); 
			StateVecD xi = StateVecD::Zero();
			StateVecD v = StateVecD::Zero();
			real rho = fvar.rho0;
			real press = 0;
			size_t pID = 0;
			pn.emplace_back(SPHPart(xi,v,rho,fvar.simM,press,FREE,pID));
			svar.simPts=1;
			svar.totPts=1;
		}
		else
		{
			real radius = 0.5*svar.diam;
			svar.dx = radius/real(0.5*(res-1));
			int nrad = ceil(radius/svar.dx);
			svar.dx = radius / real(nrad);
			
			svar.Pstep = svar.dx;
			
			fvar.simM = fvar.rho0*pow(svar.Pstep,SIMDIM); 

			// #if SIMDIM == 2
			// CreateRDroplet(svar,fvar,pn);
			// #else
			CreateDroplet(svar,fvar,pn);
			// #endif
		}

		/* Set the SPH parameters for the new size */
		fvar.H = fvar.Hfac*svar.Pstep;
		fvar.HSQ = fvar.H*fvar.H; 
		fvar.sr = 4*fvar.HSQ; 	/*KDtree search radius*/

		#if SIMDIM == 2
			#ifdef CUBIC
				fvar.correc = 10.0 / (7.0 * M_PI * fvar.H * fvar.H);
			#else
				fvar.correc = 7.0 / (4.0 * M_PI * fvar.H * fvar.H);
			#endif
		#endif
		#if SIMDIM == 3
			#ifdef CUBIC
				fvar.correc = (1.0/(M_PI*fvar.H*fvar.H*fvar.H));
			#else
				fvar.correc = (21/(16*M_PI*fvar.H*fvar.H*fvar.H));
			#endif
		#endif

		avar.GetYcoef(fvar, /* fvar.H */ svar.Pstep);

		
        cout << endl << "Original resolution: " << res<< " Actual resolution: " << 
			new_res[res_i] << "  Particle count: " << svar.simPts << endl;
		
		// real tot_mass = fvar.simM * svar.simPts;

        ///********* Tree algorithm stuff ************/
        KDTREE TREE(pn,cells);
        
        TREE.NP1.index->buildIndex();
        FindNeighbours(TREE.NP1, fvar, pn, outlist);

        dSPH_PreStep(fvar,pn.size(),pn,outlist,dp);
        Detect_Surface(svar,fvar,avar,0,pn.size(),dp,outlist,cells,pn);

        /* Remove particles that aren't going to receive a force */
        SPHState to_test;
		vector<StateVecD> norm;
		vector<real> lam;
		#pragma omp parallel default(shared)/* shared(leftright) */
		{
			SPHState local;
			vector<StateVecD> norm_l;
			vector<real> lam_l;
			#pragma omp for schedule(static) nowait
			for(size_t ii = 0; ii < pn.size(); ++ii)
			{
				if( dp.lam_ng[ii] < avar.cutoff )
				{
					local.emplace_back(pn[ii]);
					norm_l.emplace_back(dp.norm[ii]);
					lam_l.emplace_back(dp.lam_ng[ii]);
				}
			}

			#pragma omp for schedule(static) ordered
			for(int ii = 0; ii < omp_get_num_threads(); ++ii)
			{
				#pragma omp ordered
				{
					to_test.insert(to_test.end(),local.begin(),local.end());
					norm.insert(norm.end(),norm_l.begin(),norm_l.end());
					lam.insert(lam.end(),lam_l.begin(),lam_l.end());
				}
			}
		}

		cout << "Test size: " << to_test.size() << endl;

        for(size_t Re_i = 0; Re_i < svar.Reynolds.size(); ++Re_i)
        {
			real const& Re =  svar.Reynolds[Re_i];
			real const vel = (avar.mug*Re)/(avar.rhog*svar.diam);
            StateVecD tAcc_g = StateVecD::Zero();
            StateVecD sAcc_g = StateVecD::Zero();
            StateVecD pAcc_g = StateVecD::Zero();
            StateVecD tAcc_ip = StateVecD::Zero();
            StateVecD sAcc_ip = StateVecD::Zero();
            StateVecD pAcc_ip = StateVecD::Zero();
			cout << "Calculating drag for Reynolds: " << Re;

            // #pragma omp parallel for
            for(size_t ii = 0; ii < to_test.size(); ++ii)
            {
                SPHPart& pi = to_test[ii];
                // pi.cellV[1] = vel;
                pi.cellP = avar.pRef;
				pi.cellRho = avar.rhog;

                StateVecD Vdiff = StateVecD::Zero();
                Vdiff[1] = vel;
				
                StateVecD totalAcc_g, sphereAcc_g, plateAcc_g, totalAcc_ip, sphereAcc_ip, plateAcc_ip;

                
				/* Original Gissler */
				Droplet_GisslerForce(avar,Vdiff,norm[ii],avar.rhog,avar.pRef,pi.m,lam[ii],pi.woccl,
					totalAcc_g,sphereAcc_g,plateAcc_g);
			
			
				/* Induced pressure based model */	
				Droplet_InducedPressure(avar,Vdiff,norm[ii],0.0,lam[ii],pi, 
						totalAcc_ip, sphereAcc_ip, plateAcc_ip);
                

                tAcc_g += totalAcc_g * pi.m;
                sAcc_g += sphereAcc_g * pi.m;
                pAcc_g += plateAcc_g * pi.m;
				tAcc_ip += totalAcc_ip * pi.m;
                sAcc_ip += sphereAcc_ip * pi.m;
                pAcc_ip += plateAcc_ip * pi.m;
            }

			// real dynp = 0.5*(tot_mass/(M_PI*0.25*svar.diam*svar.diam))*vel*vel*svar.diam;
			// real dynp = 1.0;
			#if SIMDIM == 3
			real dynp = 0.5 * (vel*vel/(avar.sos*avar.sos)) * avar.gamma * avar.pRef * 0.25 * M_PI * svar.diam * svar.diam;
			#else
			real dynp = 0.5 * (vel*vel/(avar.sos*avar.sos)) * avar.gamma * avar.pRef * svar.diam;
			#endif

            drop_data[res_i][Re_i][0] = tAcc_g[1]/dynp;
			drop_data[res_i][Re_i][1] = sAcc_g[1]/dynp;
			drop_data[res_i][Re_i][2] = pAcc_g[1]/dynp;
			drop_data[res_i][Re_i][3] = tAcc_ip[1]/dynp;
			drop_data[res_i][Re_i][4] = sAcc_ip[1]/dynp;
			drop_data[res_i][Re_i][5] = pAcc_ip[1]/dynp;

			cout << "    Gissler: " << drop_data[res_i][Re_i][0] << "  Induced Pressure: " << drop_data[res_i][Re_i][3] << endl;
        }
		cout << "Resolutions remaining: " << svar.nacross.size() - res_i - 1 << endl;
    }

    /* Write data to file */
    std::ofstream fout("Droplet_Sweep.dat",std::ios::out);

    fout << "VARIABLES = \"diameter resolution\", \"Gissler normalised drag\"," <<
        "\"Gissler total drag\", \"Gissler sphere drag\", \"Gissler plate drag\"" << 
		"\"Induced pressure normalised drag\", " << 
		"\"Induced pressure total drag\", \"Induced pressure sphere drag \", \"Induced pressure plate drag\"" << endl; 

    for(size_t ii = 0; ii < svar.Reynolds.size(); ++ii)
    {
		/* Drag for the full droplet */
        real const Re = svar.Reynolds[ii];
        real const big_cd = GetCd(Re);

        fout << "ZONE T=\"Re = " << std::to_string(Re) << "\"" << endl;
        for(size_t jj = 0; jj < svar.nacross.size(); ++jj)
        {
            real drag_g = (drop_data[jj][ii][0])/big_cd;
			real drag_ip = (drop_data[jj][ii][3])/big_cd  ;
			// real drag_g = (drop_data[jj][ii][0])/big_drag;
			// real drag_ip = (drop_data[jj][ii][3])/big_drag;
            fout << new_res[jj] << " " << drag_g << " " << drop_data[jj][ii][0] << " " << drop_data[jj][ii][1] << " " << drop_data[jj][ii][2] <<
				" " << drag_ip << " " << drop_data[jj][ii][3] << " " << drop_data[jj][ii][4] << " " << drop_data[jj][ii][5] << endl;
        }
		
		/* Ideal line curve */
		fout << "ZONE T=\"Re = " << std::to_string(Re) << " Idealised\"" << endl;
		for(size_t jj = 0; jj < svar.nacross.size(); ++jj)
		{
			fout << new_res[jj] << " 1 " << big_cd << " " << big_cd << " " << big_cd << 
						   " 1 " << big_cd << " " << big_cd << " " << big_cd << endl;
			// fout << res << " 1 " << big_drag << " " << big_drag << " " << big_drag << 
			// 			   " 1 " << big_drag << " " << big_drag << " " << big_drag << endl;
		}
    }

    fout.close();


    /* Write data to file */
    std::ofstream fout2("Droplet_Sweep_Transpose.dat",std::ios::out);

    fout2 << "VARIABLES = \"Reynolds\", \"Gissler normalised drag\"," <<
        "\"Gissler total drag\", \"Gissler sphere drag\", \"Gissler plate drag\"" << 
		"\"Induced pressure normalised drag\", " << 
		"\"Induced pressure total drag\", \"Induced pressure sphere drag \", \"Induced pressure plate drag\"" << endl; 

    for(size_t jj = 0; jj < svar.nacross.size(); ++jj)
    {
        fout2 << "ZONE T=\"Resolution = " << std::to_string(new_res[jj]) << "\"" << endl;
        for(size_t ii = 0; ii < svar.velocities.size(); ++ii)
        {
			real const vel = svar.velocities[ii];
			/* Drag for the full droplet */
			real const Re = avar.rhog*vel*svar.diam/avar.mug;
			real const big_cd = GetCd(Re);
			
            real drag_g = (drop_data[jj][ii][0] - big_cd)/big_cd + 1;
			real drag_ip = (drop_data[jj][ii][3] - big_cd)/big_cd + 1 ;
			// real drag_g = (drop_data[jj][ii][0] - big_drag)/big_drag;
			// real drag_ip = (drop_data[jj][ii][3] - big_drag)/big_drag;
            fout2 << Re << " " << drag_g << " " << drop_data[jj][ii][0] << " " << drop_data[jj][ii][1] << " " << drop_data[jj][ii][2] <<
				" " << drag_ip << " " << drop_data[jj][ii][3] << " " << drop_data[jj][ii][4] << " " << drop_data[jj][ii][5] << endl;
        }
		
		
    }

	/* Ideal line curve */
	fout2 << "ZONE T=\"Idealised\"" << endl;
	for(size_t ii = 0; ii < svar.velocities.size(); ++ii)
	{
		real const vel = svar.velocities[ii];
		/* Drag for the full droplet */
		real const Re = avar.rhog*vel*svar.diam/avar.mug;

		// real res = svar.diam/svar.dx;
		fout2 << Re << " 1 " << "1" << " " << "1" << " " << "1" << 
						" 1 " << "1" << " " << "1" << " " << "1" << endl;
	}

    fout2.close();
}


#endif