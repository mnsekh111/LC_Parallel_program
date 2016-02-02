/*
 * MPI_Parallel_Int_dif_v1.c
 *
 *  Created on: Jan 18, 2016
 *      Author: mns
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* The number of grid points */
#define   NGRID           100
/* first grid point */
#define   XI              1.0
/* last grid point */
#define   XF              100.0

#define ROOT 0

/* floating point precision type definitions */
typedef double FP_PREC;

/* function declarations */
FP_PREC fn(FP_PREC);
FP_PREC dfn(FP_PREC);
FP_PREC ifn(FP_PREC, FP_PREC);
void print_function_data(int, FP_PREC*, FP_PREC*, FP_PREC*);
void print_error_data(int np, FP_PREC, FP_PREC, FP_PREC*, FP_PREC*, FP_PREC);
int main(int, char**);

int main(int argc, char *argv[]) {
	int taskId, totaltasks, i, j;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	MPI_Comm_size(MPI_COMM_WORLD, &totaltasks);

	FP_PREC xc[NGRID / totaltasks + 2];
	FP_PREC yc[NGRID / totaltasks + 2];
	FP_PREC dyc[NGRID / totaltasks + 2];
	FP_PREC derr[NGRID / totaltasks + 2];

	FP_PREC intg;
	FP_PREC dx;

	int prev_task = (taskId - 1) < 0 ? totaltasks - 1 : taskId - 1;
	int next_task = (taskId + 1) % totaltasks;

	MPI_Request reqs[4];
	MPI_Status stats[4];

	MPI_Irecv(&yc[0], 1, MPI_DOUBLE, prev_task, prev_task * 1000 + taskId,
	MPI_COMM_WORLD, &reqs[0]);

	MPI_Irecv(&yc[NGRID / totaltasks + 1], 1, MPI_DOUBLE, next_task,
			next_task * 1000 + taskId, MPI_COMM_WORLD, &reqs[1]);

	for (i = 1; i <= NGRID / totaltasks; i++) {
		xc[i] = (XI + (XF - XI) * (FP_PREC) (i - 1) / (FP_PREC) (NGRID - 1))
				+ taskId * NGRID / totaltasks;
	}

	//define the function
	for (i = 1; i <= NGRID / totaltasks; i++) {
		yc[i] = fn(xc[i]);
	}

	MPI_Isend(&yc[NGRID / totaltasks], 1, MPI_DOUBLE, next_task,
			taskId * 1000 + next_task, MPI_COMM_WORLD, &reqs[3]);

	MPI_Isend(&yc[1], 1, MPI_DOUBLE, prev_task, taskId * 1000 + prev_task,
	MPI_COMM_WORLD, &reqs[2]);



	MPI_Waitall(4, reqs, stats);
	dx = xc[2] - xc[1];
	if (taskId == ROOT) {
		xc[0] = xc[1] - dx;
		yc[0] = fn(xc[0]);
		printf("<<<<%f %f>>R>>\n", xc[0], yc[0]);
		printf("<<<<%f %f>>R>>\n", xc[NGRID/totaltasks+1], yc[NGRID/totaltasks+1]);
	}


	if (taskId == totaltasks - 1) {
		xc[NGRID / totaltasks + 1] = xc[NGRID / totaltasks] + dx;
		yc[NGRID / totaltasks + 1] = fn(xc[NGRID / totaltasks + 1]);
		printf("<<<<%f %f>>NR>>\n", xc[NGRID / totaltasks + 1],
				yc[NGRID / totaltasks + 1]);
		printf("<<<<%f %f>>NR>>\n", xc[0],
						yc[0]);
	}

	//compute the derivative using first-order finite differencing
	for (i = 1; i <= NGRID / totaltasks; i++) {
		dyc[i] = (yc[i + 1] - yc[i - 1]) / (2.0 * dx);
	}

	//compute the integral using Trapazoidal rule
	intg = 0.0;
	for (i = 1; i < NGRID / totaltasks; i++) {
		intg += 0.5 * (xc[i + 1] - xc[i]) * (yc[i + 1] + yc[i]);
	}

	//compute the errors
	for (i = 1; i <= NGRID / totaltasks; i++) {
		derr[i] = fabs((dyc[i] - dfn(xc[i])) / dfn(xc[i]));
	}

	if (taskId != ROOT) {

		MPI_Request nreqs[2];
		MPI_Status nstats[2];

		MPI_Isend(derr + 1, NGRID / totaltasks, MPI_DOUBLE, ROOT,
				taskId * 1000 + ROOT,
				MPI_COMM_WORLD, &nreqs[0]);
		MPI_Isend(&intg, 1, MPI_DOUBLE, ROOT, taskId * 1000 + ROOT,
		MPI_COMM_WORLD, &nreqs[1]);
		MPI_Waitall(2, nreqs, nstats);
		for (i = 1; i <= NGRID / totaltasks; i++)
			printf("%f----->\n", derr[i]);
	}

	else {

		for (i = 1; i <= NGRID / totaltasks; i++)
			printf("%f********>\n", derr[i]);
		FP_PREC allxc[NGRID];
		FP_PREC allderr[NGRID];
		FP_PREC allintg[totaltasks];
		FP_PREC davg_err = 0.0;
		FP_PREC dstd_dev = 0.0;
		FP_PREC intg_err = 0.0;
		MPI_Request nreqs[2 * (totaltasks - 1)];
		MPI_Status nstats[2 * (totaltasks - 1)];

		for (i = 1; i < totaltasks; i++) {
			MPI_Irecv(allderr + (i * NGRID / totaltasks), NGRID / totaltasks,
			MPI_DOUBLE, i, i * 1000 + ROOT, MPI_COMM_WORLD,
					&nreqs[2 * (i - 1)]);

			MPI_Irecv(allintg + i, 1, MPI_DOUBLE, i, i * 1000 + ROOT,
			MPI_COMM_WORLD, &nreqs[2 * (i - 1) + 1]);
		}

		for (i = 0; i < NGRID / totaltasks; i++) {
			allderr[i] = derr[i + 1];
		}

		MPI_Waitall(2 * (totaltasks - 1), nreqs, nstats);

		for (i = 0; i < NGRID; i++)
			printf("%d %f\n", i, allderr[i]);
		//find the average error
		for (i = 0; i < NGRID; i++)
			davg_err += allderr[i];

		for (i = 1; i < totaltasks; i++) {
			intg += allintg[i];
		}

		davg_err /= (FP_PREC) NGRID;

		dstd_dev = 0.0;
		for (i = 0; i < NGRID; i++) {
			dstd_dev += pow(allderr[i] - davg_err, 2);
		}
		dstd_dev = sqrt(dstd_dev / (FP_PREC) NGRID);

		intg_err = fabs((ifn(XI, XF) - intg) / ifn(XI, XF));

		for (i = 0; i < NGRID; i++) {
			allxc[i] = XI + (XF - XI) * (FP_PREC) i / (FP_PREC) (NGRID - 1);
		}

		//print_error_data(NGRID, davg_err, dstd_dev, &xc[1], derr, intg_err);
		print_error_data(NGRID, davg_err, dstd_dev, allxc, allderr, intg_err);
	}

	MPI_Finalize();
}

//prints out the function and its derivative to a file
void print_function_data(int np, FP_PREC *x, FP_PREC *y, FP_PREC *dydx) {
	int i;

	FILE *fp = fopen("fn2.dat", "w");

	for (i = 0; i < np; i++) {
		fprintf(fp, "%f %f %f\n", x[i], y[i], dydx[i]);
	}

	fclose(fp);
}

void print_error_data(int np, FP_PREC avgerr, FP_PREC stdd, FP_PREC *x,
		FP_PREC *err, FP_PREC ierr) {
	int i;
	FILE *fp = fopen("err2.dat", "w");

	fprintf(fp, "%e\n%e\n%e\n", avgerr, stdd, ierr);
	for (i = 0; i < np; i++) {
		fprintf(fp, "%e %e \n", x[i], err[i]);
	}
	fclose(fp);
}

FP_PREC fn(FP_PREC x) {
	return sqrt(x);
//  return x;
}

//returns the derivative d(fn)/dx = dy/dx
FP_PREC dfn(FP_PREC x) {
	return 0.5 * (1.0 / sqrt(x));
//  return 1;
}

//returns the integral from a to b of y(x) = fn
FP_PREC ifn(FP_PREC a, FP_PREC b) {
	return (2. / 3.) * (pow(sqrt(b), 3) - pow(sqrt(a), 3));
//  return 0.5 * (b*b - a*a);
}
