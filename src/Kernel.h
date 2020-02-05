#ifndef KERNEL_H
#define KERNEL_H

#include "Var.h"

///******Wendland's C2 Quintic Kernel*******
inline real W2Kernel(const real dist, const real H, const real correc) 
{
	const real q = dist/H;
	return (pow(1-0.5*q,4))*(2*q+1)*correc;
}

/*Gradient*/
inline StateVecD W2GradK(const StateVecD& Rij, const real dist, const real H, const real correc)
{
	const real q = dist/H;
	return 5.0*(Rij/(H*H))*pow(1-0.5*q,3)*correc;
}

/*2nd Gradient*/
inline real W2Grad2(const StateVecD& Rij, const real dist, const real H, const real correc) 
{
	const real q = dist/H;
	return Rij.dot(Rij)*(5.0*correc/(H*H))*(2*q-1)*pow(1-0.5*q,2);
}

#endif