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

	int num_tasks, task_id, rc, len, i, j, k, l;

	double start_time, end_time;
	double avg_rtt, std_dev;

	char * dummy_rec_data, *dummy_send_data, *processor_name = (char*) malloc(
			100);
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
	MPI_Get_processor_name(processor_name, &len);

//	for (i = 0; i < num_tasks; i++) {
//		if (task_id == i)
//			printf("Node %d = %s", task_id, processor_name);
//		MPI_Barrier(MPI_COMM_WORLD);
//	}
	double rtt_table[MSG_SIZES][REPS];
	double avg_stddev_table[MSG_SIZES * 2];
	MPI_Barrier(MPI_COMM_WORLD);

	if (num_tasks % 2 != 0) {
		//Ignoring the last node which has no pair
		num_tasks--;
	}

	if (num_tasks == 0) {
		//This condition occurs when number of tasks is either 0 or 1
		//The program won't run
		MPI_Finalize();
		return 0;
	}

	int msg_size = 32;

	for (i = 0; i < MSG_SIZES; i++) {
		dummy_send_data = (char *) malloc(msg_size);
		dummy_rec_data = (char *) malloc(msg_size);
		memset(dummy_send_data, 'A', msg_size - 1);
		dummy_send_data[msg_size - 1] = '\0';

		for (j = 0; j < num_tasks / 2; j++) {
			if ((task_id == (j * 2)) || (task_id == (j * 2) + 1)) {
				avg_rtt = 0;
				std_dev = 0;
				for (k = 0; k < REPS; k++) {
					if (task_id % 2 == 0) {
						start_time = MPI_Wtime();
						rc = MPI_Send(dummy_send_data, msg_size, MPI_CHAR,
								task_id + 1, get_uniq_tag(task_id, k),
								MPI_COMM_WORLD);
						if (rc != MPI_SUCCESS) {
							fprintf(stdout, "Send to %d failed\n", task_id + 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

						rc = MPI_Recv(dummy_rec_data, msg_size, MPI_CHAR,
								task_id + 1, get_uniq_tag(task_id + 1, k),
								MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						if (rc != MPI_SUCCESS) {
							fprintf(stdout, "Receive back from %d failed\n",
									task_id + 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

						end_time = MPI_Wtime();
						rtt_table[i][k] = end_time - start_time;
						avg_rtt += rtt_table[i][k];

					} else {
						rc = MPI_Recv(dummy_rec_data, msg_size, MPI_CHAR,
								task_id - 1, get_uniq_tag(task_id - 1, k),
								MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						if (rc != MPI_SUCCESS) {
							fprintf(stdout, "Receive from %d failed\n",
									task_id - 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

						rc = MPI_Send(dummy_rec_data, msg_size, MPI_CHAR,
								task_id - 1, get_uniq_tag(task_id, k),
								MPI_COMM_WORLD);
						if (rc != MPI_SUCCESS) {
							fprintf(stdout, "Send back to %d failed\n",
									task_id - 1);
							MPI_Abort(MPI_COMM_WORLD, rc);
							exit(1);
						}

					}
				}

				if (task_id % 2 == 0) {
					avg_rtt -= rtt_table[i][0];
					avg_rtt /= (REPS - 1);
					for (k = 1; k < REPS; k++) {
						std_dev += pow(fabs(avg_rtt - rtt_table[i][k]), 2);
					}
					std_dev /= (REPS - 1);
					std_dev = sqrt(std_dev);
					avg_stddev_table[2 * i] = avg_rtt;
					avg_stddev_table[2 * i + 1] = std_dev;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		free(dummy_rec_data);
		free(dummy_send_data);
		msg_size *= 2;
	}

	if (task_id == ROOT) {
		double display_array[num_tasks / 2][MSG_SIZES * 2];
		//double*display_array = (double*) malloc(
		//num_tasks / 2 * sizeof(double) * 2 * MSG_SIZES);
		for (i = 1; i < num_tasks / 2; i++) {
			MPI_Recv(display_array[i], 2 * MSG_SIZES,
			MPI_DOUBLE, i * 2, get_uniq_tag(i * 2, ROOT),
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		for (j = 0; j < 2 * MSG_SIZES; j++)
			display_array[ROOT][j] = avg_stddev_table[j];

		msg_size = 32;

		for (i = 0; i < MSG_SIZES; i++, msg_size *= 2) {
			printf("%d ", msg_size);
			for (j = 0; j < num_tasks / 2; j++) {
				printf("%e %e ", display_array[j][2 * i],
						display_array[j][2 * i + 1]);
			}
			printf("\n");
		}
	}

	if (task_id % 2 == 0) {
//		printf("Task %d\n", task_id);
//		for (i = 0; i < 2 * MSG_SIZES; i += 2) {
//			printf("%f %f ", avg_stddev_table[i], avg_stddev_table[i + 1]);
//			printf("\n");
//		}

		if (task_id != ROOT)
			MPI_Send(avg_stddev_table, 2 * MSG_SIZES, MPI_DOUBLE, ROOT,
					get_uniq_tag(task_id, ROOT), MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}

int get_uniq_tag(int task_id, int rep) {
	return task_id * 10000 + rep;
}
