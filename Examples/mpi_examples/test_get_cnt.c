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
      // Pick a random count of doubles to send to process pe
      srandom(pe);
      r = ((double) random() / (double)RAND_MAX) ;
      ncount = (int) ( r * (double) MAX_NUMBERS);
/*-------------------- generate these numbers */
      for (k=0; k< ncount; k++)
	list[k] = cos((double) (2*k+1) *PI / ((double) MAX_NUMBERS));
      // Send these to process pe 
      MPI_Send(list, ncount, MPI_DOUBLE, pe, 0, MPI_COMM_WORLD);
      printf("0 sent %d numbers to %d\n", ncount,pe);
    }
  } else{
    // Receive at most MAX_NUMBERS from process zero
    MPI_Recv(list, MAX_NUMBERS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
	     &status);
    // After receiving the message, check  status to see
    // how many numbers were actually received
    MPI_Get_count(&status, MPI_DOUBLE, &ncount);
    // Print off the amount of numbers, and also print additional
    // information in the status object
    printf("PE %d received %d numbers from 0. Message source = %d, "
	   "tag = %d\n",myid,ncount, status.MPI_SOURCE, status.MPI_TAG);
  }
  MPI_Finalize();
  return(0);
}
