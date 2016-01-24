#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ROOT 0
#define REPS 10
#define MSG_SIZES 17

int get_uniq_tag(int task_id, int rep);

int main(int argc, char * argv[]) {

	int num_tasks, task_id, rc, i, j, k, l;

	double start_time, end_time;
	double avg_rtt, std_dev;

	char * dummy_rec_data, *dummy_send_data;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

	double table[MSG_SIZES][num_tasks / 2][REPS];
//	printf("%d", num_tasks);
	MPI_Barrier(MPI_COMM_WORLD);

	int msg_size = 32;
	int curr_pair = 0;

	for (i = 0; i < MSG_SIZES; i++, msg_size *= 2) {
		dummy_send_data = (char *) malloc(msg_size);
		dummy_rec_data = (char *) malloc(msg_size);
		memset(dummy_send_data, 'A', msg_size - 1);

		for (j = 0; j < num_tasks / 2; j++) {
			if (task_id == j*2 || task_id == j*2 + 1) {
				avg_rtt = 0;
				std_dev = 0;
				for (k = 0; k < REPS; k++) {
					if (task_id % 2 == 0) {
						start_time = MPI_Wtime();
						rc = MPI_Send(dummy_send_data, msg_size, MPI_CHAR,
								task_id + 1, get_uniq_tag(task_id, k),
								MPI_COMM_WORLD);
						if (rc != MPI_SUCCESS) {
							fprintf(stderr, "Send to %d failed\n", task_id + 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

						rc = MPI_Recv(dummy_rec_data, msg_size, MPI_CHAR,
								task_id + 1, get_uniq_tag(task_id + 1, k),
								MPI_COMM_WORLD, &status);

						if (rc != MPI_SUCCESS) {
							fprintf(stderr, "Receive back from %d failed\n",
									task_id + 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

						end_time = MPI_Wtime();
						table[i][j][k] = end_time - start_time;
						avg_rtt += table[i][j][k];

					} else {
						rc = MPI_Recv(dummy_rec_data, msg_size, MPI_CHAR,
								task_id - 1, get_uniq_tag(task_id - 1, k),
								MPI_COMM_WORLD, &status);
						if (rc != MPI_SUCCESS) {
							fprintf(stderr, "Receive from %d failed\n",
									task_id - 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

						rc = MPI_Send(dummy_rec_data, msg_size, MPI_CHAR,
								task_id - 1, get_uniq_tag(task_id, k),
								MPI_COMM_WORLD);
						if (rc != MPI_SUCCESS) {
							fprintf(stderr, "Send back to %d failed\n",
									task_id - 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

					}
				}

				if (task_id % 2 == 0) {
					for (k = 0; k < REPS; k++) {
						std_dev += pow(fabs(avg_rtt - table[i][j][k]), 2);
					}
					std_dev /= REPS;
					std_dev = sqrt(std_dev);
					printf("Message size %d from %d to %d --> %f %f \n", i,
							task_id, task_id + 1, avg_rtt, std_dev);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		free(dummy_rec_data);
		free(dummy_send_data);
	}

//for (k = 0; k < num_tasks / 2; k++) {
//	if (task_id == curr_pair || task_id == curr_pair + 1) {
//		for (j = 0; j < MSG_SIZES; j++) {
//			avg_rtt = 0;
//			std_dev = 0;
//			for (i = 0; i < REPS; i++) {
//				if (task_id % 2 == 0) {
//					dummy_send_data = (char *) malloc(msg_size);
//					dummy_rec_data = (char *) malloc(msg_size);
//					memset(dummy_send_data, 'A', msg_size - 1);
//					start_time = MPI_Wtime();
//					rc = MPI_Send(dummy_send_data, msg_size, MPI_CHAR,
//							task_id + 1, get_uniq_tag(task_id, i),
//							MPI_COMM_WORLD);
//					if (rc != MPI_SUCCESS) {
//						fprintf(stderr, "Send to %d failed\n", task_id + 1);
//						MPI_Abort(MPI_COMM_WORLD, rc);
//						exit(1);
//					}
//
//					//printf("Sent data from task %d\n",task_id);
//					rc = MPI_Recv(dummy_rec_data, msg_size, MPI_CHAR,
//							task_id + 1, get_uniq_tag(task_id + 1, i),
//							MPI_COMM_WORLD, &status);
//
//					if (rc != MPI_SUCCESS) {
//						fprintf(stderr, "Receive back from %d failed\n",
//								task_id + 1);
//						MPI_Abort(MPI_COMM_WORLD, rc);
//						exit(1);
//					}
//
////						printf("Received back %d bytes of data from task %d\n",
////								strlen(dummy_rec_data), task_id + 1);
//					end_time = MPI_Wtime();
//					table[j][k][i] = end_time - start_time;
//					avg_rtt += table[j][k][i];
//					printf("Message size %d from %d to %d --> %f \n", j,
//							task_id, task_id + 1, end_time - start_time);
//					free(dummy_send_data);
//				} else {
//					dummy_rec_data = (char *) malloc(msg_size);
//					memset(dummy_rec_data, ' ', msg_size - 1);
//					rc = MPI_Recv(dummy_rec_data, msg_size, MPI_CHAR,
//							task_id - 1, get_uniq_tag(task_id - 1, i),
//							MPI_COMM_WORLD, &status);
//					if (rc != MPI_SUCCESS) {
//						fprintf(stderr, "Receive from %d failed\n",
//								task_id - 1);
//						MPI_Abort(MPI_COMM_WORLD, rc);
//						exit(1);
//					}
//
////						printf("Received data (tag:%d) from task %d\n",
////								status.MPI_TAG, task_id - 1);
//					rc = MPI_Send(dummy_rec_data, msg_size, MPI_CHAR,
//							task_id - 1, get_uniq_tag(task_id, i),
//							MPI_COMM_WORLD);
//					if (rc != MPI_SUCCESS) {
//						fprintf(stderr, "Send back to %d failed\n",
//								task_id - 1);
//						MPI_Abort(MPI_COMM_WORLD, rc);
//						exit(1);
//					}
//					//printf("Sent data back to task %d\n",task_id - 1);
//					free(dummy_rec_data);
//				}
//			}
//
//			for (l = 0; l < REPS; l++) {
//				std_dev += pow(fabs(avg_rtt - table[]))
//			}
//			msg_size *= 2;
//		}
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	curr_pair += 2;
//}
	MPI_Finalize();
	return 0;
}

int get_uniq_tag(int task_id, int rep) {
	return task_id * 10000 + rep;
}
