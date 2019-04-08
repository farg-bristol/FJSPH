#ifndef KERNEL_H
#define KERNEL_H

#include "Eigen/Core"
#include "Var.h"

///******Wendland's C2 Quintic Kernel*******
ldouble W2Kernel(ldouble dist,ldouble H, ldouble correc) 
{
	ldouble q = dist/H;
	return (pow(1-0.5*q,4))*(2*q+1)*correc;
}

/*Gradient*/
StateVecD W2GradK(StateVecD Rij, ldouble dist, ldouble H, ldouble correc)
{
	ldouble q = dist/H;
	return 5.0*(Rij/(H*H))*pow(1-0.5*q,3)*correc;
}

/*2nd Gradient*/
ldouble W2Grad2(StateVecD Rij, ldouble dist, ldouble H, ldouble correc) 
{
	ldouble q = dist/H;
	return Rij.dot(Rij)*(5.0*correc/(H*H))*(2*q-1)*pow(1-0.5*q,2);
}

#endif