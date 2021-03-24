#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> 

#define MAX_NUMBERS 100
#define PI M_PI
double cos(double);

int main(int argc, char **argv) {  
  double list[MAX_NUMBERS];
  double r;
  int ncount, k, myid, pe, nproc;
  MPI_Status status;
  MPI_Request req;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (myid == 0) {
    for (pe = 1; pe < nproc; pe++) {
//-------------------- Pick a random count of doubles to send to process pe
      srandom(pe);
      r = ((double) random() / (double)RAND_MAX) ;
      ncount = (int) ( r * (double) MAX_NUMBERS);
/*-------------------- generate these numbers */
      for (k=0; k< ncount; k++)
	list[k] = cos((double) (2*k+1) *PI / ((double) MAX_NUMBERS));
//--------------------  Send these to process pe 
      MPI_Send(list, ncount, MPI_DOUBLE, pe, 0, MPI_COMM_WORLD);
      printf("0 sent %d numbers to %d\n", ncount,pe);
    }
  } else{
//--------------------Probe for an incoming message from process zero
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
//--------------------When probe returns, 'status' contains size + other
//                    info on incoming message. Get message size
    MPI_Get_count(&status, MPI_DOUBLE, &ncount);
//--------------------Allocate buffer big enough to hold incoming list
    double* list_buf = (double*)malloc(sizeof(double) * ncount);
//--------------------Now receive the message with the allocated buffer
    MPI_Recv(list_buf, ncount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("PE # %d -- dynamically received %d doubles from 0.\n",
           myid,ncount);
    free(list_buf);
  }
  MPI_Finalize();
  return(0);
}
