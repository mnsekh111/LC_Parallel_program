#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ROOT 0
#define REPS 10

int min_element(int*, int);
int main(int argc, char* argv[]) {
	int taskId, totaltasks, i, j;
	int smallest_num;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	MPI_Comm_size(MPI_COMM_WORLD, &totaltasks);

	MPI_Status status;
	char * sendMessage = (char*) malloc(10);
	char * recMessage = (char*) malloc(10);
	memset(sendMessage, 'A', 9);
	sendMessage[9]='\0';

	for (i = 0; i < totaltasks / 2; i++) {
		if ((i * 2) == taskId || (i * 2 + 1) == taskId) {
			for (j = 0; j < REPS; j++) {
				if (taskId % 2 == 0) {
					MPI_Send(sendMessage, 10, MPI_CHAR, taskId + 1, 0,
					MPI_COMM_WORLD);
					printf("Task %d sent message %s to task %d\n", taskId,
							sendMessage, taskId + 1);
					MPI_Recv(recMessage, 10, MPI_CHAR, taskId + 1, 0,
					MPI_COMM_WORLD, &status);
					printf("Task %d received message %s from task %d\n", taskId,
							recMessage, taskId + 1);
				} else {
					MPI_Recv(recMessage, 10, MPI_CHAR, taskId - 1, 0,
					MPI_COMM_WORLD, &status);
					printf("Task %d received message %s from task %d\n", taskId,
							recMessage, taskId - 1);
					memset(recMessage, 'A'+j, 9);
					MPI_Send(recMessage, 10, MPI_CHAR, taskId - 1, 0,
					MPI_COMM_WORLD);
					printf("Task %d sent message %s to task %d\n", taskId,
							sendMessage, taskId - 1);

				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		MPI_Finalize();
	}
}
