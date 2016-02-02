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
/*
 Single Author info:
 avelayu Ashitha Velayudhan
 Group info:
 1. avelayu Ashitha Velayudhan
 2. prajago4 Priyadarshini Rajagopal
 3. smnatara Sekharan Muthusamy Natarajan
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* The number of grid points */
#define   NGRID           100
/* first grid point */
#define   XI              1.0
/* last grid point */
#define   XF              100.0

/* floating point precision type definitions */
typedef double FP_PREC;

#define LB(RANK,NUM_PROC,NUM_ELEMENTS) (((RANK)*(NUM_ELEMENTS))/(NUM_PROC))
#define UB(RANK,NUM_PROC,NUM_ELEMENTS) ((LB((RANK)+1,NUM_PROC,NUM_ELEMENTS))-1)
#define BLOCK_SIZE(RANK,NUM_PROC,NUM_ELEMENTS) ((UB(RANK,NUM_PROC,NUM_ELEMENTS) - LB(RANK,NUM_PROC,NUM_ELEMENTS))+1)
#define DEFAULT_TAG 0

/* function declarations */
FP_PREC fn(FP_PREC);
FP_PREC dfn(FP_PREC);
FP_PREC ifn(FP_PREC, FP_PREC);
void print_function_data(int, FP_PREC*, FP_PREC*, FP_PREC*);
void print_error_data(int np, FP_PREC, FP_PREC, FP_PREC*, FP_PREC*, FP_PREC);
int main(int, char**);

void MPI_Manual_Reduce(void* intg, void* intg_reduced, int num_procs) {
	int half = num_procs / 2;
	int i;
	int rank;
	int ROOT = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == ROOT) {
		*(FP_PREC*) intg_reduced = *(FP_PREC*) intg;
		for (i = 1; i < num_procs; i++) {
			/* Receive from everyone else */
			MPI_Recv(intg, 1, MPI_DOUBLE, i, DEFAULT_TAG, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
			*(FP_PREC*) intg_reduced = *(FP_PREC*) (intg_reduced)
					+ *(FP_PREC*) intg;
		}
	} else {
		/* Send to root for reduction */
		MPI_Send(intg, 1, MPI_DOUBLE, ROOT, DEFAULT_TAG, MPI_COMM_WORLD);
	}
}

int main(int argc, char *argv[]) {

	//loop index
	int i, j;

	//domain array and step size
	FP_PREC *xc, dx;

	//function array and derivative
	//the size will be dependent on the
	//number of processors used
	//to the program
	FP_PREC *yc, *dyc;

	//"real" grid indices
	int imin, imax;

	//integration values
	FP_PREC intg;
	FP_PREC intg_reduced = 0.0;

	//error analysis array
	FP_PREC *derr;
	FP_PREC *derr_root;

	//error analysis values
	FP_PREC davg_err, dstd_dev, intg_err = 0.0, global_avg_err = 0.0;
	FP_PREC davg_err_local;
	int rank, num_procs, len;
	int lower, upper, block_size;
	char name[MPI_MAX_PROCESSOR_NAME];
	double prg_exec_time = 0.0;

	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf("Error in MPI INIT");
		return -1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Get_processor_name(name, &len);

	int rcounts[num_procs];
	int *displs = malloc(num_procs * sizeof(int));

	int ROOT = 0;
	int LAST_NODE = num_procs - 1;

#ifdef DEBUG
	double time_comm = 0.0,
	time_comm_total = 0.0,
	time_diff = 0.0,
	time_diff_total = 0.0,
	time_integral = 0.0,
	time_integral_total = 0.0,
	time_reduction = 0.0,
	time_reduction_total = 0.0,
	time_err = 0.0,
	time_err_total = 0.0;
#endif

	MPI_Request reqs[4];
	MPI_Status stats[4];

	if (rank == ROOT) {
#if 1
#ifdef BLOCKING
		printf("Blocking using ");
#else
		printf("Non-blocking using ");
#endif
#ifdef MANUAL
		printf("manual reduction");
#else
		printf("single reduction");
#endif
		printf("\n");
#endif
	}

	prg_exec_time = -MPI_Wtime();

	// Each node will have a subset with size=block_size.
	block_size = BLOCK_SIZE(rank, num_procs, NGRID);
	imin = 1;
	imax = block_size;
	lower = LB(rank,num_procs,NGRID) + imin;
	upper = UB(rank,num_procs,NGRID) + imin;

	// Allocate function arrays
	xc = (FP_PREC*) malloc((block_size + 2) * sizeof(FP_PREC));
	yc = (FP_PREC*) malloc((block_size + 2) * sizeof(FP_PREC));
	dyc = (FP_PREC*) malloc((block_size + 2) * sizeof(FP_PREC));

	// Construct grid
	for (i = imin - 1; i <= imax + 1; i++) {
		xc[i] = XI
				+ (XF - XI) * (FP_PREC) (lower + i - 2) / (FP_PREC) (NGRID - 1);
	}

	// Step size and Boundary points
	if (rank == ROOT) {
		dx = xc[2] - xc[1];
	}
	/* dx calculation:
	 While testing with sin(x) function we noticed a slight difference in caclulation of dx value in different nodes
	 due to floating point roundoff errors. To avoid that dx is calculated in root node and broad catsed to all
	 other nodes.  */

	MPI_Bcast(&dx, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	//xc[0] = xc[1] - dx;
	//xc[imax + 1] = xc[imax] + dx;
	if (rank == ROOT) {
		xc[0] = xc[1] - dx;
	}
	if (rank == LAST_NODE) {
		xc[imax + 1] = xc[imax] + dx;
	}

	// Define the function
	for (i = imin; i <= imax; i++) {
		yc[i] = fn(xc[i]);
	}

	//------------------------------------------------//
	// set boundary values: 4 values needs to be sent //
	//------------------------------------------------//
	// 1: imax of even nodes to imin-1 of  odd nodes  //
	// 2: imin of  odd nodes to imax+1 of even nodes  //
	// 3: imin of even nodes to imax+1 of  odd nodes  //
	// 4: imax of  odd nodes to imin-1 of even nodes  //
	//------------------------------------------------//
	//  +------+          +------+          +------+  //
	//  |      |  imax    |      |  imax    |      |  //
	//  |      |--------->|      |--------->|      |  //
	//  | even |          | odd  |          | even |  //
	//  |      |  imin    |      |  imin    |      |  //
	//  |      |<---------|      |<---------|      |  //
	//  +------+          +------+          +------+  //
	//------------------------------------------------//

	MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
	time_comm = -MPI_Wtime();
#endif

#ifdef BLOCKING
	if(rank%2 == 0) {
		if (rank!= LAST_NODE)
		{
			/* Send and receive only if it is not the last node. */
			MPI_Send(&yc[imax] , 1, MPI_DOUBLE, rank+1, DEFAULT_TAG, MPI_COMM_WORLD);
			MPI_Recv(&yc[imax+1], 1, MPI_DOUBLE, rank+1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != ROOT) {
			MPI_Send(&yc[imin] , 1, MPI_DOUBLE, rank-1, DEFAULT_TAG, MPI_COMM_WORLD);
			MPI_Recv(&yc[imin-1], 1, MPI_DOUBLE, rank-1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Recv(&yc[imin-1], 1, MPI_DOUBLE, rank-1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(&yc[imin] , 1, MPI_DOUBLE, rank-1, DEFAULT_TAG, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != LAST_NODE) {
			MPI_Recv(&yc[imax+1], 1, MPI_DOUBLE, rank+1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&yc[imax] , 1, MPI_DOUBLE, rank+1, DEFAULT_TAG, MPI_COMM_WORLD);
		}
	}
#else
	if (rank % 2 == 0) {
		MPI_Irecv(&yc[imax + 1], 1, MPI_DOUBLE, rank + 1, DEFAULT_TAG,
				MPI_COMM_WORLD, &reqs[0]);
		MPI_Isend(&yc[imax], 1, MPI_DOUBLE, rank + 1, DEFAULT_TAG,
				MPI_COMM_WORLD, &reqs[1]);

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != ROOT) {
			MPI_Irecv(&yc[imin - 1], 1, MPI_DOUBLE, rank - 1, DEFAULT_TAG,
					MPI_COMM_WORLD, &reqs[2]);
			MPI_Isend(&yc[imin], 1, MPI_DOUBLE, rank - 1, DEFAULT_TAG,
					MPI_COMM_WORLD, &reqs[3]);
		}
	} else {
		MPI_Irecv(&yc[imin - 1], 1, MPI_DOUBLE, rank - 1, DEFAULT_TAG,
				MPI_COMM_WORLD, &reqs[1]);
		MPI_Isend(&yc[imin], 1, MPI_DOUBLE, rank - 1, DEFAULT_TAG,
				MPI_COMM_WORLD, &reqs[0]);

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != LAST_NODE) {
			MPI_Irecv(&yc[imax + 1], 1, MPI_DOUBLE, rank + 1, DEFAULT_TAG,
					MPI_COMM_WORLD, &reqs[3]);
			MPI_Isend(&yc[imax], 1, MPI_DOUBLE, rank + 1, DEFAULT_TAG,
					MPI_COMM_WORLD, &reqs[2]);
		}
	}
	if (rank == 0) {
		MPI_Waitall(2, reqs, stats);
	} else if (rank == LAST_NODE) {
		MPI_Waitall(2, reqs, stats);
	} else {
		MPI_Waitall(4, reqs, stats);
	}
#endif

	// NB: boundary values of the whole domain
	// should be set
	if (rank == ROOT)
		yc[0] = fn(xc[0]);
	if (rank == LAST_NODE)
		yc[imax + 1] = fn(xc[imax + 1]);

#ifdef DEBUG
	time_comm += MPI_Wtime();
	MPI_Reduce(&time_comm, &time_comm_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	time_comm_total/= num_procs;
#endif

#ifdef DEBUG
	printf("R:%d Lower=%d, Upper=%d, Block Size=%d\n", rank, lower, upper, block_size);
	for (i = imin-1; i <= imax+1; i++) {
		printf("R:%d, xc[%d]=%0.6lf, yc[%d]=%e\n",rank,i,xc[i],i,yc[i]);
	}
#endif // DEBUG

	//compute the derivative using first-order finite differencing
	//
	//  d           f(x + h) - f(x - h)
	// ---- f(x) ~ --------------------
	//  dx                 2 * dx
	//

#ifdef DEBUG
	time_diff = -MPI_Wtime();
#endif

	for (i = imin; i <= imax; i++) {
		dyc[i] = (yc[i + 1] - yc[i - 1]) / (2.0 * dx);
	}

#ifdef DEBUG
	/* find the average of all times from all nodes. */
	time_diff += MPI_Wtime();
	MPI_Reduce(&time_diff, &time_diff_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	time_diff_total/= num_procs;
#endif
	//compute the integral using Trapazoidal rule
	//
	//    _b
	//   |
	//   | f(x) dx ~ (b - a) / 2 * (f(b) + f(a))
	//  _|
	//   a
	//
#ifdef DEBUG
	time_integral = -MPI_Wtime();
#endif

	intg = 0.0;
	for (i = imin; i <= imax; i++) {
		//there are NGRID points, so there are NGRID-1 integration zones.
		//if(i - imin != NGRID - 1) intg += 0.5 * (xc[i + 1] - xc[i]) * (yc[i + 1] + yc[i]);
		if (rank == LAST_NODE && i == imax) {
#ifdef DEBUG
//            printf("skip :xc[%d]:%e\n", i,xc[i]);
#endif // DEBUG
			continue;
		}
		intg += 0.5 * (xc[i + 1] - xc[i]) * (yc[i + 1] + yc[i]);
	}

#ifdef DEBUG
	time_integral += MPI_Wtime();

	/* find the average of all times from all nodes. */
	MPI_Reduce(&time_integral, &time_integral_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	time_integral_total /= num_procs;

	time_reduction = -MPI_Wtime();
#endif

	// Reduction operation
	// Do MPI_Reduce() call here
#ifdef MANUAL
	/* Reduce using manual reduction implemented using MPI_Send */
	MPI_Manual_Reduce(&intg, &intg_reduced,num_procs);
#else
	/* Reduce using single reduction. */
	MPI_Reduce(&intg, &intg_reduced, 1, MPI_DOUBLE, MPI_SUM, ROOT,
			MPI_COMM_WORLD);
#endif

#ifdef DEBUG
	time_reduction += MPI_Wtime();
	MPI_Reduce(&time_reduction, &time_reduction_total , 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	time_reduction_total /= num_procs;
#endif

	// compute the error, average error of the derivatives
	derr = (FP_PREC*) malloc(block_size * sizeof(FP_PREC));

	if (rank == ROOT)
		derr_root = (FP_PREC*) malloc(NGRID * sizeof(FP_PREC));
	// compute the errors

#ifdef DEBUG
	time_err = -MPI_Wtime();
#endif
	for (i = imin; i <= imax; i++) {
		derr[i - imin] = fabs((dyc[i] - dfn(xc[i])) / dfn(xc[i]));
	}
#ifdef DEBUG
	time_err += MPI_Wtime();
	MPI_Reduce(&time_err, &time_err_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	time_err_total /= num_procs;
#endif

	if (rank == ROOT) {
#ifdef DEBUG
		time_err = -MPI_Wtime();
#endif
		displs[0] = 0;
		rcounts[0] = BLOCK_SIZE(0, num_procs, NGRID);
		for (j = 1; j < num_procs; j++) {
			rcounts[j] = BLOCK_SIZE(j, num_procs, NGRID);
			displs[j] = rcounts[j - 1] + displs[j - 1];
			//MPI_Recv(derr_root+, 1, MPI_DOUBLE, j, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	MPI_Gatherv(derr, block_size, MPI_DOUBLE, derr_root, rcounts, displs,
			MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

#if DEBUG
	if (rank == ROOT) {
		for(i = 0; i < NGRID; i++)
		printf("%e %e \n", xc[i+1], derr_root[i]);
	}
#endif

	//find the average error
	/* Here calculate the sum of errors and communicate that to the root, instead of average. */

#if 0
	//(i) Calculate the average error for the local domain
	davg_err_local = 0.0;
	for(i = 0; i < block_size; i++)
	davg_err_local += derr[i];
	//davg_err_local /= (FP_PREC)block_size;
	// (ii) Communicate the local average error to a single process, and calculate the average error for the entire domain
	if (rank != ROOT) {
		MPI_Send(&davg_err_local , 1, MPI_DOUBLE, ROOT, DEFAULT_TAG, MPI_COMM_WORLD);
	}
	else {
		global_avg_err += davg_err_local;
		for(j=1; j<num_procs;j++)
		{
			MPI_Recv(&davg_err, 1, MPI_DOUBLE, j, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			global_avg_err += davg_err;
		}
		global_avg_err /= (FP_PREC)NGRID;
	}
#else
	/* calculate the average error at the ROOT. This is more accurate */
	if (rank == ROOT) {
		for (i = 0; i < NGRID; i++) {
			global_avg_err += derr_root[i];
		}
		global_avg_err /= (FP_PREC) NGRID;
	}
#endif

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
#ifdef DEBUG
	//printf("davg_err:%e\n",global_avg_err);
#endif
	if (rank == ROOT) {
		// (iii) use the global average error to calculate the standard deviation
		for (i = 0; i < NGRID; i++) {
			dstd_dev += pow(derr_root[i] - global_avg_err, 2);
#ifdef DEBUG
			//printf("derr[%d] %.17g \n",i,derr_root[i]);
			//printf("davg_err[%d] %.17g \n",i,global_avg_err);
			//printf("dstdev:%d %.17g \n",i,derr_root[i] - global_avg_err);
#endif
		}
		dstd_dev = sqrt(dstd_dev / (FP_PREC) NGRID);
		intg_err = fabs((ifn(XI, XF) - intg_reduced) / ifn(XI, XF));
		//print_function_data(NGRID, &xc[1], &yc[1], &dyc[1]);
#ifdef DEBUG
		time_err += MPI_Wtime();
		time_err_total += time_err;
#endif

		print_error_data(NGRID, global_avg_err, dstd_dev, &xc[1], derr_root,
				intg_err);

	}

	//free allocated memory
	free(xc);
	free(yc);
	free(dyc);
	free(derr);
	if (rank == ROOT)
		free(derr_root);

	if (rank == ROOT) {
		prg_exec_time += MPI_Wtime();
#ifdef DEBUG
		printf("\nTime for finite difference    = %.17g\n", time_diff_total);
		printf("Time for integral calculation = %.17g\n", time_integral_total);
		printf("Time for error calculation    = %.17g\n", time_err_total);
		printf("Time for communication        = %.17g\n", time_comm_total);
		printf("Time for reduction            = %.17g\n", time_reduction_total);
#endif
		printf("End of Program: Total Exection Time = %.17g\n", prg_exec_time);

		// Printing Statistics
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

FP_PREC fn(FP_PREC x) {
#ifdef SIN
	return sin(x);
#elif POW
	return pow(x,2);
#else
	return sqrt(x);
#endif
	//return (sinh(x));
//  return x;

}

//returns the derivative d(fn)/dx = dy/dx
FP_PREC dfn(FP_PREC x) {
#ifdef SIN
	return cos(x);
#elif POW
	return (2*x);
#else
	return 0.5 * (1.0 / sqrt(x));
#endif

	//return (cosh(x));
//  return 1;

}

//returns the integral from a to b of y(x) = fn
FP_PREC ifn(FP_PREC a, FP_PREC b) {
#ifdef SIN
	return -cos(b) + cos(a);
#elif POW
	return (pow(b,3)/3.0) - (pow(a,3)/3.);
#else
	return (2. / 3.) * (pow(sqrt(b), 3) - pow(sqrt(a), 3));
#endif

	//return (cosh(b) - cosh(a));
//  return 0.5 * (b*b - a*a);
}

//prints out the function and its derivative to a file
void print_function_data(int np, FP_PREC *x, FP_PREC *y, FP_PREC *dydx) {
	int i;

	FILE *fp = fopen("fn3.dat", "w");

	for (i = 0; i < np; i++) {
		fprintf(fp, "%f %f %f\n", x[i], y[i], dydx[i]);
	}

	fclose(fp);
}

void print_error_data(int np, FP_PREC avgerr, FP_PREC stdd, FP_PREC *x,
		FP_PREC *err, FP_PREC ierr) {
	int i;
	FILE *fp = fopen("err3.dat", "w");
	double xc = 0.0;

	if (fp == NULL)
		printf("Error\n");
	fprintf(fp, "%e\n%e\n%e\n", avgerr, stdd, ierr);
	for (i = 0; i < np; i++) {
		/*TODO: How to get all xc values to root. */
		//fprintf(fp, "%e %e \n", x[i], err[i]);
		xc = XI + (XF - XI) * (FP_PREC) i / (FP_PREC) (NGRID - 1);
		fprintf(fp, "%e %e \n", xc, err[i]);
	}
	fclose(fp);
}
