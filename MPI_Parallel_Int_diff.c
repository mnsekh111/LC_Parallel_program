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
	int CHUNK;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	MPI_Comm_size(MPI_COMM_WORLD, &totaltasks);

	CHUNK = NGRID / totaltasks;
	FP_PREC *xc = (FP_PREC*) malloc((CHUNK + 2) * sizeof(FP_PREC));
	FP_PREC *yc = (FP_PREC*) malloc((CHUNK + 2) * sizeof(FP_PREC));
	FP_PREC *dyc = (FP_PREC*) malloc((CHUNK + 2) * sizeof(FP_PREC));
	FP_PREC *derr = (FP_PREC*) malloc((CHUNK + 2) * sizeof(FP_PREC));

	FP_PREC intg;
	FP_PREC dx;

	int prev_task = (taskId - 1) < 0 ? totaltasks - 1 : taskId - 1;
	int next_task = (taskId + 1) % totaltasks;

	MPI_Request reqs[4];
	MPI_Status stats[4];

	MPI_Irecv(&yc[0], 1, MPI_DOUBLE, prev_task, prev_task * 1000 + taskId,
	MPI_COMM_WORLD, &reqs[0]);

	MPI_Irecv(&yc[CHUNK + 1], 1, MPI_DOUBLE, next_task,
			next_task * 1000 + taskId, MPI_COMM_WORLD, &reqs[1]);

	for (i = 1; i <= CHUNK; i++) {
		xc[i] = (XI + (XF - XI) * (FP_PREC) (i - 1) / (FP_PREC) (NGRID - 1))
				+ taskId * CHUNK;
	}

	//define the function
	for (i = 1; i <= CHUNK; i++) {
		yc[i] = fn(xc[i]);
	}

	MPI_Isend(&yc[1], 1, MPI_DOUBLE, prev_task, taskId * 1000 + prev_task,
	MPI_COMM_WORLD, &reqs[2]);

	MPI_Isend(&yc[CHUNK], 1, MPI_DOUBLE, next_task, taskId * 1000 + next_task,
	MPI_COMM_WORLD, &reqs[3]);

	MPI_Waitall(4, reqs, stats);

	printf("waiting complete ");
	dx = xc[2] - xc[1];
	if (taskId == ROOT) {
		xc[0] = xc[1] - dx;
		yc[0] = fn(xc[0]);
	}

	if (taskId == totaltasks - 1) {
		xc[CHUNK + 1] = xc[CHUNK] + dx;
		yc[CHUNK + 1] = fn(xc[CHUNK + 1]);
	}

	//compute the derivative using first-order finite differencing
	for (i = 1; i <= CHUNK; i++) {
		dyc[i] = (yc[i + 1] - yc[i - 1]) / (2.0 * dx);
	}

	//compute the integral using Trapazoidal rule
	intg = 0.0;
	for (i = 1; i < CHUNK; i++) {
		intg += 0.5 * (xc[i + 1] - xc[i]) * (yc[i + 1] + yc[i]);
	}

	//compute the errors
	for (i = 1; i <= CHUNK; i++) {
		derr[i] = fabs(dfn(xc[i]) - dyc[i]) / dfn(xc[i]);
	}

	if (taskId != ROOT) {

		MPI_Request nreqs[2];
		MPI_Status nstats[2];

		MPI_Isend(&derr[1], CHUNK, MPI_DOUBLE, ROOT, taskId * 1000 + ROOT,
		MPI_COMM_WORLD, &nreqs[0]);
		MPI_Isend(&intg, 1, MPI_DOUBLE, ROOT, taskId * 1000 + ROOT,
		MPI_COMM_WORLD, &nreqs[1]);
		MPI_Waitall(2, nreqs, nstats);

	} else {

		FP_PREC allxc[NGRID];
		FP_PREC davg_err = 0.0;
		FP_PREC dstd_dev = 0.0;
		FP_PREC intg_err = 0.0;

		FP_PREC allderr1D[NGRID];
		FP_PREC * allintg = (FP_PREC *) malloc(totaltasks * sizeof(FP_PREC));
		FP_PREC ** allderr = (FP_PREC **) malloc(totaltasks * sizeof(FP_PREC*));
		MPI_Request *nreqs = (MPI_Request*) malloc(
				2 * (totaltasks - 1) * sizeof(MPI_Request));
		MPI_Status *nstats = (MPI_Status*) malloc(
				2 * (totaltasks - 1) * sizeof(MPI_Status));

		for (i = 1; i < totaltasks; i++) {

			allderr[i] = (FP_PREC *) malloc(CHUNK * sizeof(FP_PREC*));

			MPI_Irecv(allderr[i], CHUNK,
			MPI_DOUBLE, i, i * 1000 + ROOT, MPI_COMM_WORLD,
					&nreqs[2 * (i - 1)]);

			MPI_Irecv(allintg + i, 1, MPI_DOUBLE, i, i * 1000 + ROOT,
			MPI_COMM_WORLD, &nreqs[2 * (i - 1) + 1]);
		}

		allderr[ROOT] = (FP_PREC *) malloc(CHUNK * sizeof(FP_PREC*));
		for (j = 0; j < CHUNK; j++) {
			allderr[ROOT][j] = derr[i + 1];
		}

		MPI_Waitall(2 * (totaltasks - 1), nreqs, nstats);

		for (i = 1; i < totaltasks; i++) {
			intg += allintg[i];
		}

		//find the average error
		for (i = 0; i < totaltasks; i++)
			for (j = 0; j < CHUNK; j++)
				davg_err += allderr[i][j];

		davg_err /= (FP_PREC) NGRID;

		dstd_dev = 0.0;

		for (i = 0; i < totaltasks; i++)
			for (j = 0; j < CHUNK; j++)
				dstd_dev += pow(allderr[i][j] - davg_err, 2);

		dstd_dev = sqrt(dstd_dev / (FP_PREC) NGRID);

		intg_err = fabs((ifn(XI, XF) - intg) / ifn(XI, XF));

		for (i = 0; i < NGRID; i++) {
			allxc[i] = XI + (XF - XI) * (FP_PREC) i / (FP_PREC) (NGRID - 1);
		}

		for (i = 0; i < totaltasks; i++)
			for (j = 0; j < NGRID / totaltasks; j++)
				allderr1D[i * NGRID + j] = allderr[i][j];
		//print_error_data(NGRID, davg_err, dstd_dev, &xc[1], derr, intg_err);
		print_error_data(NGRID, davg_err, dstd_dev, allxc, allderr1D, intg_err);

		free(allintg);
		free(nreqs);
		free(nstats);
		for (i = 0; i < totaltasks; i++)
			free(allderr[i]);
		free(allderr);
	}

	free(xc);
	free(yc);
	free(dyc);
	free(derr);

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
