/******************************************************************************
* FILE: hw1.c
* DESCRIPTION:
*
* Here are the functions
* i.) fn(x) - the function to be analyized
* ii.) dfn(x) - the true derivative of the function
* iii.) ifn(x) - the true integral of the function
*
* The function fn(x) should be smooth and continuous, and
* the derivative and integral should have analyitic expressions
* on the entire domain.
*
* AUTHOR: Christopher Mauney
* LAST REVISED: 08/18/12
******************************************************************************/
#include <math.h>

/* floating point precision type definitions */
typedef   double   FP_PREC;

//returns the function y(x) = fn
FP_PREC fn(FP_PREC x)
{
  return sqrt(x);
//  return x;
}

//returns the derivative d(fn)/dx = dy/dx
FP_PREC dfn(FP_PREC x)
{
  return 0.5*(1.0/sqrt(x));
//  return 1;
}

//returns the integral from a to b of y(x) = fn
FP_PREC ifn(FP_PREC a, FP_PREC b)
{
  return (2./3.) * (pow(sqrt(b), 3) - pow(sqrt(a),3));
//  return 0.5 * (b*b - a*a);
}

