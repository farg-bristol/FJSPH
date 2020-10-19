#ifndef KERNEL_H
#define KERNEL_H

#include "Var.h"

///******Wendland's C2 Quintic Kernel*******
real const W2Kernel(real const dist, real const H, real const correc) 
{
	return (pow(1-0.5*dist/H,4))*(2*dist/H+1)*correc;
}

/*Gradient*/
StateVecD const W2GradK(StateVecD const& Rij, real const dist, real const H, real const correc)
{
	return 5.0*(Rij/(H*H))*pow(1-0.5*dist/H,3)*correc;
}

/*2nd Gradient*/
real const W2Grad2(StateVecD const& Rij, real const dist, real const H, real const correc) 
{
	return Rij.dot(Rij)*(5.0*correc/(H*H))*(2*dist/H-1)*pow(1-0.5*dist/H,2);
}

real const BoundaryKernel(real const dist, real const H, real const beta)
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