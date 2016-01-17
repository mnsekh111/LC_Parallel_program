#include "mpi.h"
#include <stdio.h>

main(int argc, char *argv[])  {
int numtasks, rank, dest, source, rc, count, tag=1; 
int MSG_LEN = 10; 
char inmsg[10],*outmsg="x+x+x+x+x";
MPI_Status Stat;

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if (rank == 0) {
  dest = 1;
  source = 1;
  rc = MPI_Send(outmsg,MSG_LEN, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  rc = MPI_Recv(inmsg,MSG_LEN, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
  } 

else if (rank == 1) {
  dest = 0;
  source = 0;
  rc = MPI_Recv(inmsg,MSG_LEN, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
  rc = MPI_Send(outmsg,MSG_LEN, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  }

rc = MPI_Get_count(&Stat, MPI_CHAR, &count);
printf("Task %d: Received %d char(s) from task %d with tag %d %s\n",
       rank,MSG_LEN, Stat.MPI_SOURCE, Stat.MPI_TAG,inmsg);

MPI_Finalize();
}
