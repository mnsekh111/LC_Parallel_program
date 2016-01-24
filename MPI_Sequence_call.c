/*
 * MPI_Scatter_Gather.c
 *
 *  Created on: Jan 18, 2016
 *      Author: mns
 */

#include "mpi.h"
#include <stdio.h>
#include <limits.h>
#define SIZE 4
#define ROOT 0

int min_element(int*, int);
int main(int argc, char* argv[]) {
	int taskId, totaltasks, i;
	int smallest_num;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	MPI_Comm_size(MPI_COMM_WORLD, &totaltasks);

	for (i = 0; i < totaltasks; i++) {
		if (i == taskId) {
			printf("Printing from %d\n", taskId);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Finalize();
}


