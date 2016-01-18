/******************************************************************************
* FILE: hw1.c
* DESCRIPTION:
*
* Users will supply the functions
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* The number of grid points */
#define   NGRID           100
/* first grid point */
#define   XI              1.0
/* last grid point */
#define   XF              100.0

/* floating point precision type definitions */
typedef   double   FP_PREC;

/* function declarations */
FP_PREC     fn(FP_PREC);
FP_PREC     dfn(FP_PREC); 
FP_PREC     ifn(FP_PREC, FP_PREC);
void        print_function_data(int, FP_PREC*, FP_PREC*, FP_PREC*);
void        print_error_data(int np, FP_PREC, FP_PREC, FP_PREC*, FP_PREC*, FP_PREC);
int         main(int, char**);

int main (int argc, char *argv[])
{
  //loop index
  int         i;

  //domain array and step size
  FP_PREC       xc[NGRID + 2], dx;

  //function array and derivative
  //the size will be dependent on the
  //number of processors used
  //to the program
  FP_PREC     *yc, *dyc;

  //"real" grid indices
  int         imin, imax;  
  
  //integration values
  FP_PREC     intg;

  //error analysis array
  FP_PREC     *derr;

  //error analysis values
  FP_PREC     davg_err, dstd_dev, intg_err;

  imin = 1;
  imax = NGRID;

  //construct grid
  for (i = 1; i <= NGRID ; i++)
  {
    xc[i] = XI + (XF - XI) * (FP_PREC)(i - 1)/(FP_PREC)(NGRID - 1);
  }
  //step size and boundary points
  dx = xc[2] - xc[1];
  xc[0] = xc[1] - dx;
  xc[NGRID + 1] = xc[NGRID] + dx;

  //allocate function arrays
  yc  =   (FP_PREC*) malloc((NGRID + 2) * sizeof(FP_PREC));
  dyc =   (FP_PREC*) malloc((NGRID + 2) * sizeof(FP_PREC));

  //define the function
  for( i = imin; i <= imax; i++ )
  {
    yc[i] = fn(xc[i]);
  }

  //set boundary values
  yc[imin - 1] = 0.0;
  yc[imax + 1] = 0.0;

  //NB: boundary values of the whole domain
  //should be set
  yc[0] = fn(xc[0]);
  yc[imax + 1] = fn(xc[NGRID + 1]);

  //compute the derivative using first-order finite differencing
  //
  //  d           f(x + h) - f(x - h)
  // ---- f(x) ~ --------------------
  //  dx                 2 * dx
  //
  for (i = imin; i <= imax; i++)
  {
    dyc[i] = (yc[i + 1] - yc[i - 1])/(2.0 * dx);
  }

  //compute the integral using Trapazoidal rule
  //
  //    _b
  //   |
  //   | f(x) dx ~ (b - a) / 2 * (f(b) + f(a))
  //  _|
  //   a
  //
  intg = 0.0;
  for (i = imin; i <= imax; i++)
  {
    //there are NGRID points, so there are NGRID-1 integration zones.
    if(i - imin != NGRID - 1) intg += 0.5 * (xc[i + 1] - xc[i]) * (yc[i + 1] + yc[i]);
  }

  //compute the error, average error of the derivatives
  derr = (FP_PREC*)malloc(NGRID * sizeof(FP_PREC));

  //compute the errors
  for(i = imin; i <= imax; i++)
  {
    derr[i-imin] = fabs((dyc[i] - dfn(xc[i]))/dfn(xc[i]));
  }

  //find the average error
  davg_err = 0.0;
  for(i = 0; i < NGRID ; i++)
    davg_err += derr[i];
    
  davg_err /= (FP_PREC)NGRID;

  //find the standard deviation of the error
  //standard deviation is defined to be
  //
  //                   ____________________________
  //          __      /      _N_
  // \omega =   \    /  1    \   
  //             \  /  --- *  >  (x[i] - avg_x)^2 
  //              \/    N    /__  
  //                        i = 1 
  //
  dstd_dev = 0.0;
  for(i = 0; i< NGRID; i++)
  {
    dstd_dev += pow(derr[i] - davg_err, 2);
  }
  dstd_dev = sqrt(dstd_dev/(FP_PREC)NGRID);
  
  intg_err = fabs((ifn(XI, XF) - intg)/ifn(XI, XF));
        
  print_function_data(NGRID, &xc[1], &yc[1], &dyc[1]);
  print_error_data(NGRID, davg_err, dstd_dev, &xc[1], derr, intg_err);

 
  //free allocated memory 
  free(yc);
  free(dyc);
  free(derr);

  return 0;
}

//prints out the function and its derivative to a file
void print_function_data(int np, FP_PREC *x, FP_PREC *y, FP_PREC *dydx)
{
  int   i;

  FILE *fp = fopen("fn.dat", "w");

  for(i = 0; i < np; i++)
  {
    fprintf(fp, "%f %f %f\n", x[i], y[i], dydx[i]);
  }

  fclose(fp);
}

void print_error_data(int np, FP_PREC avgerr, FP_PREC stdd, FP_PREC *x, FP_PREC *err, FP_PREC ierr)
{
  int   i;
  FILE *fp = fopen("err.dat", "w");

  fprintf(fp, "%e\n%e\n%e\n", avgerr, stdd, ierr);
  for(i = 0; i < np; i++)
  {
    fprintf(fp, "%e %e \n", x[i], err[i]);
  }
  fclose(fp);
}
