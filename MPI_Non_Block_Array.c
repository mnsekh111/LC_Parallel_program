#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 10
#define ROOT 0

int main(int argc, char *argv[]) {
	int taskId, totaltasks, i, j;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	MPI_Comm_size(MPI_COMM_WORLD, &totaltasks);

	int * sampleArray = (int *) malloc(SIZE * sizeof(int));
	for (i = 0; i < SIZE; i++)
		sampleArray[i] = taskId * 10 + i;

	if (taskId == ROOT) {
		MPI_Request *reqs = (MPI_Request*) malloc(
				(totaltasks - 1) * sizeof(MPI_Request));
		MPI_Status *stats = (MPI_Status*) malloc(
				(totaltasks - 1) * sizeof(MPI_Status));

		int ** aggrArray = (int **) malloc((totaltasks-1)*sizeof(int*));


		for (i = 1; i < totaltasks; i++) {
			aggrArray[i] = (int *) malloc(SIZE*sizeof(int));
			MPI_Irecv(aggrArray[i], SIZE, MPI_INT, i, i * 1000 + ROOT,
			MPI_COMM_WORLD, &reqs[i - 1]);
		}

		MPI_Waitall(totaltasks - 1, reqs, stats);

		aggrArray[ROOT] = (int *) malloc(SIZE*sizeof(int));
		for(j=0;j<SIZE;j++)
			aggrArray[ROOT][j] = sampleArray[j];
		for (i = 0; i < totaltasks; i++) {
			for(j=0;j<SIZE;j++)
				printf("%d ", aggrArray[i][j]);
			printf("\n");
		}

	}

	if (taskId > ROOT) {

		MPI_Request req;
		MPI_Status stat;
		MPI_Isend(sampleArray, SIZE, MPI_INT, ROOT, taskId * 1000 + ROOT,
		MPI_COMM_WORLD, &req);
		MPI_Wait(&req, &stat);
	}

	MPI_Finalize();
}
