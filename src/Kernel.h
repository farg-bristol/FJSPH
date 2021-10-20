#ifndef KERNEL_H
#define KERNEL_H

#include "Var.h"

#ifdef CUBIC
	///******Cubic Spline Kernel*******///
	inline real const Kernel(real const dist, real const H, real const correc)
	{
		real const q = dist/H;

		if(q < 1.0)
		{
			return correc * (1.0 - 1.5 * q*q*(1-0.5*q));
		}
		else if (q < 2.0)
		{
			return correc * 0.25 * pow(2.0-q,3.0);
		}
		
		return 0;
	}

	inline StateVecD const GradK(StateVecD const& Rij, real const dist, real const H, real const correc)
	{
		real const q = dist/H;

		if(q < 1.0)
		{
			return Rij/(H*H) * correc * ( 3.0 * (0.75*q-1.0));
		}
		else if (q < 2.0)
		{
			return Rij/(dist*H) * correc * -0.75 * (2.0-q) * (2.0-q);
		}
		
		return StateVecD::Zero();
	}

#else
	///******Wendland's C2 Quintic Kernel*******///
	inline real const Kernel(real const& dist, real const& H, real const& correc) 
	{
		// if(dist/H > 2.0)
		// {
		// 	// cout << "Distance greater than 2H" << endl;
		// 	return 0.0;
		// }
		return (pow(1-0.5*dist/H,4))*(2*dist/H+1)*correc;
	}

	/*Gradient*/
	inline StateVecD const GradK(StateVecD const& Rij, real const& dist, real const& H, real const& correc)
	{
		if(dist/H < 1e-12)
		{
			cout << "Points are too close" << endl;
			return StateVecD::Zero();
		}
		// else if(dist/H > 2.0)
		// {
		// 	// cout << "Distance greater than 2H" << endl;
		// 	return StateVecD::Zero();
		// }
		return 5.0*(Rij/(H*H))*pow(1-0.5*dist/H,3)*correc;
	}
#endif


inline real const BoundaryKernel(real const dist, real const H, real const beta)
{
	real const q = dist/H;
	if (q < 2.0/3.0)
	{
		return beta*2.0/3.0;
	}
	else if (2.0/3.0 <= q && q < 1.0)
	{
		return beta*(2*q - 3.0/2.0*q*q);
	}
	else if (1<= q && q < 2)
	{
		return 0.5*beta*pow((2-q),2);
	}
	return 0;
}

#endif