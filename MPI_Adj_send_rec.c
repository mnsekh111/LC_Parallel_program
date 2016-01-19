/*
 * MPI_Adj_send_rec.c
 *
 *  Created on: Jan 18, 2016
 *      Author: mns
 */

#include "mpi.h"
#include <stdio.h>
#include <limits.h>
#include <math.h>

#define   NGRID           100
#define   XI              1.0
#define   XF              100.0

#define ROOT 0

typedef double FP_PREC;

/* function declarations */
FP_PREC fn(FP_PREC);
FP_PREC dfn(FP_PREC);
FP_PREC ifn(FP_PREC, FP_PREC);

int main(int argc, char* argv[]) {

	int taskId, totaltasks, i, j;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	MPI_Comm_size(MPI_COMM_WORLD, &totaltasks);

	FP_PREC xc[NGRID / totaltasks + 2], dx;
	FP_PREC yc[NGRID / totaltasks + 2], dyc;

	int prev_task = (taskId - 1) < 0 ? totaltasks - 1 : taskId - 1;
	int next_task = (taskId + 1) % totaltasks;

	MPI_Request reqs[4];
	MPI_Status stats[4];

	MPI_Irecv(&xc[0], 1,MPI_DOUBLE, prev_task, prev_task * 1000 + taskId,
			MPI_COMM_WORLD, &reqs[0]);
	MPI_Irecv(&xc[NGRID / totaltasks + 1], 1, MPI_DOUBLE, next_task,
			next_task * 1000 + taskId, MPI_COMM_WORLD, &reqs[1]);

	for (i = 1; i <= NGRID / totaltasks; i++) {
		xc[i] = XI + (XF - XI) * (FP_PREC) (i - 1) / (FP_PREC) (NGRID - 1);
	}

	MPI_Isend(&xc[1], 1, MPI_DOUBLE, prev_task, taskId * 1000 + prev_task,
	MPI_COMM_WORLD, &reqs[2]);
	MPI_Isend(&xc[NGRID / totaltasks], 1,MPI_DOUBLE, next_task,
			taskId * 1000 + next_task, MPI_COMM_WORLD, &reqs[3]);

	MPI_Waitall(4, reqs, stats);

	for (j = 0; j < totaltasks; j++) {
		if (taskId == j) {
			printf("Task Id %d\n", taskId);
			for (i = 0; i <= NGRID / totaltasks + 1; i++) {
				printf("%f ", xc[i]);
			}
			printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

//	dx = xc[2] - xc[1];
//	xc[0] = xc[1] - dx;
//	xc[NGRID/totaltasks + 1] = xc[NGRID/totaltasks] + dx;

	MPI_Finalize();
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

