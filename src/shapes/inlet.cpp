#include "inlet.h"

void check_inlet_input(shape_block& block, real& globalspacing)
{
    int fault = 0;

    // if(block.subshape.empty())
    // {
    //     std::cout << "ERROR: Block \"" << block.name << 
    //         "\" sub shape type has not been defined." << std::endl;
    //     fault = 1;
    // }
    #if SIMDIM == 3
    if(block.subshape == "Square")
    {
        block.sub_bound_type = squareCube;
    }
    else if(block.subshape == "Circle")
    {
        block.sub_bound_type = circleSphere;
    }
    else
    {
        std::cout << "ERROR: Block \"" << block.name << 
            "\" sub shape type has not been correctly defined." << std::endl;
        cout << "Please choose from:" << endl;
        cout << "\t1. Square" << endl;
        cout << "\t2. Circle" << endl;    
        fault = 1;
    }
    #endif

    if(!check_vector(block.start))
    {
        std::cout << "ERROR: Block \"" << block.name << 
            "\" starting position has not been correctly defined." << std::endl;
        fault = 1;
    }    


    block.dx = globalspacing; // Potential to allow different size particles, but not right now.

    // if(block.dx < 0)
    // {
    //     block.dx = globalspacing;
    // }

    // int has_orientation = 0; // Use flag to ensure a viable configuration exists.
    int has_dimension = 0;
    if(check_vector(block.end) 
    #if SIMDIM == 3
        && check_vector(block.right)
    #endif    
        )
    {
        //Have three points to define the plane
        // has_orientation = 1;
        has_dimension = 1;
        StateVecD ab = (block.end - block.start).normalized();
        #if SIMDIM == 3
        StateVecD ac = (block.right - block.start).normalized();
        block.normal = ab.cross(ac).normalized();
        StateMatD rotmat;
        rotmat << ab[0]          , ab[1]          , ab[2]          ,
                  ac[0]          , ac[1]          , ac[2]          , 
                  block.normal[0], block.normal[1], block.normal[2];
        #else
        block.angles[0] = atan2(ab[1],ab[0]);
        StateMatD rotmat;
        rotmat << cos(block.angles(0)), sin(block.angles(0)),
            -sin(block.angles(0)),  cos(block.angles(0));
        block.normal = rotmat * StateVecD(0.0,1.0);
        #endif
        block.rotmat = rotmat;
        block.insert_norm = block.normal;
        block.insconst = block.insert_norm.dot(block.start)+0.01*globalspacing;

    }
    else if(block.angles.norm() != 0)
    {
        // has_orientation = 1;
        // Find the rotation matrix
        #if SIMDIM == 3
            StateMatD rotx, roty, rotz;
            rotx << 1.0, 0.0            , 0.0           ,
                    0.0, cos(block.angles(0)) , sin(block.angles(0)),
                    0.0, -sin(block.angles(0)), cos(block.angles(0));

            roty << cos(block.angles(1)) , 0.0 , -sin(block.angles(1)),
                    0.0            , 1.0 , 0.0            ,
                    sin(block.angles(1)) , 0.0 , cos(block.angles(1));

            rotz << cos(block.angles(2)) , sin(block.angles(2)) , 0.0 ,
                    -sin(block.angles(2)), cos(block.angles(2)) , 0.0 ,
                    0.0            , 0.0            , 1.0 ;

            block.rotmat = rotx*roty*rotz;
        #else
            StateMatD rotmat;
            rotmat << cos(block.angles(0)), sin(block.angles(0)),
                -sin(block.angles(0)),  cos(block.angles(0));

            block.rotmat = rotmat;
        #endif
    }   
    else if(!check_vector(block.normal))
    {
        // has_orientation = 1;
        // Need to find the rotation matrix now (will lose some information, as no knowledge of rotation around normal)
        block.normal.normalize();

        #if SIMDIM == 3
        StateVecD origin = StateVecD::Zero();
        origin[1] = 1.0;
        StateVecD fvec = (block.normal - origin).normalized(); // forward vector

        StateVecD axis = StateVecD::Zero();
        // check how close the forward vector is to the selected rotation axis
        if(fabs(fvec[0] - 1.0) > 0.1 && fabs(fvec[0]) < 0.1 && fabs(fvec[2]) < 0.1)
            axis[0] = 1.0;
        else
            axis[1] = 1.0;

        StateVecD right = fvec.cross(axis);

        StateMatD rotmat;
        rotmat << right[0], right[1], right[2],
                  axis[0] , axis[1] , axis[2] ,
                  fvec[0] , fvec[1] , fvec[2] ;
        
        block.rotmat = rotmat;
        #else
        // Find rotation matrix
        block.angles[0] = atan2(-block.normal[0],block.normal[1]);
        StateMatD rotmat;
        rotmat << cos(block.angles(0)), sin(block.angles(0)),
            -sin(block.angles(0)),  cos(block.angles(0));

        block.rotmat = rotmat;

        #endif
        StateVecD end = StateVecD::Zero();
        end[0] = 1.0;
        block.end = block.rotmat * end;

    } 
    else if(!check_vector(block.insert_norm))
    {
        // has_orientation = 1;
        // Need to find the rotation matrix now (will lose some information, as no knowledge of rotation around normal)
        block.normal = block.insert_norm;
        block.normal.normalize();

        #if SIMDIM == 3
        StateVecD origin = StateVecD::Zero();
        origin[1] = 1.0;
        StateVecD fvec = (block.normal - origin).normalized(); // forward vector

        StateVecD axis = StateVecD::Zero();
        // check how close the forward vector is to the selected rotation axis
        if(fabs(fvec[0] - 1.0) > 0.1 && fabs(fvec[0]) < 0.1 && fabs(fvec[2]) < 0.1)
            axis[0] = 1.0;
        else
            axis[1] = 1.0;

        StateVecD right = fvec.cross(axis);

        StateMatD rotmat;
        rotmat << right[0], right[1], right[2],
                  axis[0] , axis[1] , axis[2] ,
                  fvec[0] , fvec[1] , fvec[2] ;
        
        block.rotmat = rotmat;
        #else
        // Find rotation matrix
        block.angles[0] = atan2(block.normal[1],block.normal[0]) - M_PI/2.0;
        StateMatD rotmat;
        rotmat << cos(block.angles(0)), sin(block.angles(0)),
            -sin(block.angles(0)),  cos(block.angles(0));

        block.rotmat = rotmat;

        #endif
        StateVecD end = StateVecD::Zero();
        end[0] = 1.0;
        block.end = block.rotmat * end;

    } 
    else
    {
        // define a default normal vector
        block.normal = StateVecD::Zero();
        block.normal[1] = 1.0;
    }

    // if(!has_orientation)
    // {
    //     std::cout << "ERROR: Block \"" << block.name << "\" orientation has not been sufficiently defined." << endl;
    //     fault = 1;
    // }

    if(has_dimension)
    {
        /* Estimate how many points on the Block */
        // Need to have lengths along the axis to accurately have the counts
        real xlength = (block.end - block.start).norm();
        #if SIMDIM == 3
        real ylength = (block.right - block.start).norm();
        #endif

        int npoints;
        size_t ni, nk;
        #if SIMDIM == 3
        size_t nj;
        #endif
        if (block.hcpl == 1)
        {
            #if SIMDIM == 2
            ni = static_cast<int>(ceil(xlength / globalspacing));
            nk = static_cast<int>(ceil(block.thickness / globalspacing / sqrt(3.0) * 2.0));
            #else
            ni = static_cast<int>(ceil(xlength / globalspacing));
            nj = static_cast<int>(ceil(ylength / globalspacing / sqrt(3.0) * 2.0));
            nk = static_cast<int>(ceil(block.thickness / globalspacing / sqrt(6.0) * 3.0));
            #endif
        }
        else
        {
            ni = static_cast<int>(ceil(xlength / globalspacing));
            #if SIMDIM == 3
            nj = static_cast<int>(ceil(ylength/ globalspacing));
            #endif
            nk = static_cast<int>(ceil(block.thickness / globalspacing));
            
        }
        ni = ni > 1 ? ni : 1;
        nk = nk > 1 ? nk : 1;
        block.nk = nk;
        block.ni = ni;
        int nBuff = block.hcpl == 1 ? 5 : 4;
        #if SIMDIM == 3
        nj = nj > 1 ? nj : 1;
        block.nj = nj;
        npoints = ni * nk * nj;
        #else
        npoints = ni * (nk + nBuff);
        #endif
        npoints = npoints > 1 ? npoints : 1;

        block.npts = npoints;
    }
    else
    {
        // Check for counts
        if(block.ni > 0)
        {
            if(block.nj > 0)
            {
                #if SIMDIM == 3
                if(block.nk > 0)
                {
                #endif
                    has_dimension = 1;
                #if SIMDIM == 3
                }
                #endif
                // Find the end and right coordinates, then rotate
                block.end = StateVecD::Zero();
                block.right = StateVecD::Zero();
                block.end[0] = globalspacing * block.ni;
                block.end = block.rotmat * block.end;
                block.end += block.start;
                block.right[1] = globalspacing*block.nj;
                block.right = block.rotmat * block.right;
                block.right += block.start;
            }
        }
    }

    if(!has_dimension)
    {
        std::cout << "ERROR: Block \"" << block.name << "\" dimensions have not been sufficiently defined." << endl;
        fault = 1;
    }

    if(fault)
    {
        cout << "Geometry check finished with errors. Stopping" << endl;
        exit(-1);
    }
}

#if SIMDIM == 2
std::vector<StateVecD> create_inlet_zone(shape_block& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    StateVecD delta = (block.end - block.start) / static_cast<real>(block.ni);
    StateVecD norm(delta[1], -delta[0]);
    for (int jj = block.nk - 1; jj > 0; --jj)
    {
        for (int ii = 0; ii < block.ni; ++ii)
        {
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = delta * real(ii + 0.5*(jj % 2)) + norm * sqrt(3.0) * real(jj);
            else
                newPoint = delta * real(ii) + norm * real(jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.emplace_back(newPoint);
            block.bc.emplace_back(PIPE);

        } /* end i count */
    } /* end k count */

    // Create the back particles
    for (int ii = 0; ii < block.ni; ++ii)
    {
        StateVecD newPoint;
        if (block.hcpl == 1)
            newPoint = delta * real(ii);
        else
            newPoint = delta * real(ii);

        newPoint += StateVecD::Constant(unif(re));
        newPoint += block.start;
        points.emplace_back(newPoint);
        block.bc.emplace_back(BACK);
        block.back.emplace_back(points.size()-1);
    } /* end i count */
    
    // Buffer zone
    int nBuff = block.hcpl == 1 ? 5 : 4;
    block.buffer = vector<vector<size_t>>(block.ni, vector<size_t>(nBuff));
    size_t buff = 0;
    for (int jj = - 1; jj >= -nBuff; --jj)
    {
        for (int ii = 0; ii < block.ni; ++ii)
        {
            StateVecD newPoint;
            if (block.hcpl == 1)
                newPoint = delta * real(ii + 0.5*(jj % 2)) + norm * sqrt(3.0) * real(jj);
            else
                newPoint = delta * real(ii) + norm * real(jj);

            newPoint += StateVecD::Constant(unif(re));
            newPoint += block.start;
            points.emplace_back(newPoint);
            block.bc.emplace_back(BUFFER);
            block.buffer[ii][buff] = points.size()-1;
        } /* end i count */
        buff++;
    } /* end k count */
    return points;
}

#else

std::vector<StateVecD>  create_inlet_zone(shape_block& block, real const& globalspacing)
{
    std::vector<StateVecD> points;
    /* Perturb points so they aren't in a perfect grid */
    std::uniform_real_distribution<real> unif(0.0, MEPSILON*globalspacing);
    std::default_random_engine re; 

    // Which geometry to create?
    if(block.sub_bound_type == squareCube)
    {

        for (int kk =  block.nk - 1; kk > 0; --kk)
        {
            for (int jj = 0; jj < block.nj; ++jj)
            {
                for (int ii = 0; ii < block.ni; ++ii)
                {
                    StateVecD newPoint;
                    if (block.hcpl == 1)
                        newPoint = 0.5 * StateVecD(real(2 * ii + ((jj + kk) % 2)), 
                            sqrt(3.0) * (real(jj) + real((kk % 2)) / 3.0), 
                            2.0 / 3.0 * sqrt(6.0) * real(kk));
                    else
                        newPoint = StateVecD(real(ii), real(jj), real(kk));

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re),unif(re),unif(re));
                    newPoint = block.rotmat * newPoint;
                    newPoint += block.start;
                    points.emplace_back(newPoint);
                    block.bc.emplace_back(PIPE);
                }   /* end i count */
            }   /* end j count */
        }   /* end k count */

        // Create the back particles
        for (int jj = 0; jj < block.nj; ++jj)
        {
            for (int ii = 0; ii < block.ni; ++ii)
            {
                StateVecD newPoint;
                if (block.hcpl == 1)
                    newPoint = 0.5 * StateVecD(real(2 * ii + ((jj) % 2)), sqrt(3.0) * (real(jj)), 0.0);
                else
                    newPoint = StateVecD(real(ii), real(jj), 0.0);

                newPoint *= globalspacing;
                newPoint += StateVecD(unif(re),unif(re),unif(re));
                newPoint = block.rotmat * newPoint;
                newPoint += block.start;
                points.emplace_back(newPoint);
                block.bc.emplace_back(BACK);
                block.back.emplace_back(points.size()-1);

            }   /* end i count */
        }   /* end j count */

        // Create buffer zone
        size_t nBuff = block.hcpl == 1 ? 5 : 4;
        block.buffer = vector<vector<size_t>>(block.back.size(), vector<size_t>(nBuff));
        size_t buff = 0;
        for (int kk =  -1; kk >= -nBuff; --kk)
        {
            size_t partID = 0;
            for (int jj = 0; jj < block.nj; ++jj)
            {
                for (int ii = 0; ii < block.ni; ++ii)
                {
                    StateVecD newPoint;
                    if (block.hcpl == 1)
                        newPoint = 0.5 * StateVecD(real(2 * ii + ((jj + kk) % 2)), 
                            sqrt(3.0) * (real(jj) + real((kk % 2)) / 3.0), 
                            2.0 / 3.0 * sqrt(6.0) * real(kk));
                    else
                        newPoint = StateVecD(real(ii), real(jj), real(kk));

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re),unif(re),unif(re));
                    newPoint = block.rotmat * newPoint;
                    newPoint += block.start;
                    points.emplace_back(newPoint);
                    block.bc.emplace_back(BUFFER);
                    block.buffer[partID][buff] = points.size()-1;
                    partID++;
                }   /* end i count */
            }   /* end j count */
            buff++;
        }   /* end k count */
    }
    else if(block.sub_bound_type == circleSphere)
    {
        real const& radius = block.radius/globalspacing;
        real const rsq = radius*radius;
        StateVecD const& centre = block.centre;
        for (int kk =  block.nk - 1; kk > 0; --kk)
        {
            for (int jj = 0; jj < block.nj; ++jj)
            {
                for (int ii = 0; ii < block.ni; ++ii)
                {
                    StateVecD newPoint;
                    if (block.hcpl == 1)
                        newPoint = 0.5 * StateVecD(real(2 * ii + ((jj + kk) % 2)), 
                            sqrt(3.0) * (real(jj) + real((kk % 2)) / 3.0), 
                            2.0 / 3.0 * sqrt(6.0) * real(kk));
                    else
                        newPoint = StateVecD(real(ii), real(jj), real(kk));

                    real a = newPoint[0] - radius;
                    real b = newPoint[1] - radius;
                    if((a*a + b*b) > rsq)
                        continue;

                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re),unif(re),unif(re));
                    newPoint  = block.rotmat * newPoint;
                    newPoint += block.start;
                    points.emplace_back(newPoint);
                    block.bc.emplace_back(PIPE);

                }   /* end i count */
            }   /* end j count */
        }   /* end k count */

        // Create the back particles
        for (int jj = 0; jj < block.nj; ++jj)
        {
            for (int ii = 0; ii < block.ni; ++ii)
            {
                StateVecD newPoint;
                if (block.hcpl == 1)
                    newPoint = 0.5 * StateVecD(real(2 * ii + ((jj) % 2)), sqrt(3.0) * (real(jj)), 0.0);
                else
                    newPoint = StateVecD(real(ii), real(jj), 0.0);

                real a = newPoint[0] - radius;
                real b = newPoint[1] - radius;
                if((a*a + b*b) > rsq)
                    continue;

                newPoint *= globalspacing;
                newPoint += StateVecD(unif(re),unif(re),unif(re));
                newPoint = block.rotmat * newPoint;
                newPoint += block.start;
                points.emplace_back(newPoint);
                block.bc.emplace_back(BACK);
                block.back.emplace_back(points.size()-1);

            }   /* end i count */
        }   /* end j count */

        // Create buffer zone
        size_t nBuff = block.hcpl == 1 ? 5 : 4;
        block.buffer = vector<vector<size_t>>(block.back.size(), vector<size_t>(nBuff));
        size_t buff = 0;
        for (int kk =  -1; kk >= -nBuff; --kk)
        {
            size_t partID = 0;
            for (int jj = 0; jj < block.nj; ++jj)
            {
                for (int ii = 0; ii < block.ni; ++ii)
                {
                    StateVecD newPoint;
                    if (block.hcpl == 1)
                        newPoint = 0.5 * StateVecD(real(2 * ii + ((jj + kk) % 2)), 
                            sqrt(3.0) * (real(jj) + real((kk % 2)) / 3.0), 
                            2.0 / 3.0 * sqrt(6.0) * real(kk));
                    else
                        newPoint = StateVecD(real(ii), real(jj), real(kk));
                    
                    real a = newPoint[0] - radius;
                    real b = newPoint[1] - radius;
                    if((a*a + b*b) > rsq)
                        continue;
                        
                    newPoint *= globalspacing;
                    newPoint += StateVecD(unif(re),unif(re),unif(re));
                    newPoint = block.rotmat * newPoint;
                    newPoint += block.start;
                    points.emplace_back(newPoint);
                    block.bc.emplace_back(BUFFER);
                    block.buffer[partID][buff] = points.size()-1;
                    partID++;
                }   /* end i count */
            }   /* end j count */
            buff++;
        }   /* end k count */
    }
    return points;
}

#endif

uint update_buffer_region(SIM& svar, LIMITS& limits, SPHState& pnp1, size_t& end, size_t& end_ng)
{
    uint nAdd = 0;
    /*Check if more particles need to be created*/
	for(size_t block = svar.nbound; block < svar.nbound + svar.nfluid; block++)
	{
		if(limits[block].block_type == inletZone)
		{
			size_t& partID = svar.partID;
			/* Check if any buffer particles have become pipe particles */
			for (size_t ii = 0; ii < limits[block].back.size(); ++ii)
			{
				size_t const& pID = limits[block].back[ii];
                StateVecD const& pos = pnp1[pID].xi;
				/*Check that the starting area is clear first...*/
				if(pos.dot(limits[block].insert_norm) > limits[block].insconst)	
				{
					/* particle has cleared the back zone */
					pnp1[pID].b = PIPE;

					/* Update the back vector */
					pnp1[limits[block].buffer[ii][0]].b = BACK;
					limits[block].back[ii] = limits[block].buffer[ii][0];

					/* Update the buffer vector */
					for(size_t jj = 0; jj < limits[block].buffer[0].size()-1; ++jj)
						limits[block].buffer[ii][jj] = limits[block].buffer[ii][jj+1];
					
					/* Create a new particle */
					if(svar.totPts < svar.finPts)
					{
						StateVecD xi = pnp1[limits[block].buffer[ii].back()].xi-svar.dx * limits[block].insert_norm;

						pnp1.insert(pnp1.begin() + limits[block].index.second,
						SPHPart(xi,pnp1[limits[block].buffer[ii].back()],BUFFER,partID));
						limits[block].buffer[ii].back() = end_ng;

						limits[block].index.second++;
						// Also need to adjust all blocks after the current
						for(size_t jj = block+1; jj < svar.nbound + svar.nfluid; jj++)
						{
							limits[jj].index.first++;
							limits[jj].index.second++;
						}

						svar.simPts++;
						svar.totPts++;
						end_ng++;
						end++;
						partID++;
						nAdd++;
					}
				}
			}
		}
	}
    return nAdd;
}
