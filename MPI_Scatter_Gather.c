#include "mpi.h"
#include <stdio.h>
#include <limits.h>
#define SIZE 4
#define ROOT 0

int min_element(int*, int);
int main(int argc, char* argv[]){
	int num [SIZE][SIZE] = {4,3,2,1,9,8,7,6,14,13,12,11,19,18,17,16};
	int taskId,totaltasks,sendcount,receivecount,i,min;
	int recbuff[SIZE];
	int smallrecbuff[SIZE];

	for(i=0;i<SIZE;i++)
		smallrecbuff[i]=0;
	int smallest_num;
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&taskId);
	MPI_Comm_size(MPI_COMM_WORLD,&totaltasks);

	if(totaltasks == SIZE){
		MPI_Scatter(num,SIZE,MPI_INT,recbuff,SIZE,MPI_INT,ROOT,MPI_COMM_WORLD);
		printf("Task %d received {%d, %d, %d, %d}\n",taskId,recbuff[0],recbuff[1],recbuff[2],recbuff[3]);

		min = min_element(recbuff,SIZE);	

		MPI_Gather(&min,1,MPI_INT,smallrecbuff,1,MPI_INT,ROOT,MPI_COMM_WORLD);
		
		if(taskId == ROOT){
		  	printf("The minimum of all elements %d ", min_element(smallrecbuff,SIZE));
		}
	}
	else{
		printf("This needs %d  processes\n",SIZE);
	}	

	MPI_Finalize();
}

int min_element(int * arr , int size){
	int i,min=INT_MAX;
	for(i=0;i<size;i++){
		if(*(arr+i) < min)
			min = arr[i];
	}
	return min;
}
