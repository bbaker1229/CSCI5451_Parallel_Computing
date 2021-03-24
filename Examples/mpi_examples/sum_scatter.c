#include "mpi.h"
#include <stdio.h>
#include <stdlib.h> 
#define MAXSIZE 1000

int main(int argc, char *argv[]){
  int myid, numprocs;
  int data[MAXSIZE], locsum, sum, sum0; 
  int j, n, sd, m; 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid); 
/*-------------------- generate n random numbers - m in each PE */
  m = 30;
  n = m*numprocs;
  sum0 = 0; 
  /*-------------------- a sort of seed */
  sd = 231;   
 
  if (myid == 0) {
    printf(" -- running partial sums -- output from proc 0\n");
    for (j=0; j<n; j++) {
      data[j] = ((j+1)*117 + sd)  %  100;
      sum0 += data[j] ;
    }
  } 
/*-------------------- scatter data to processors            
  Note: in PE0 data is overwritten  */
  MPI_Scatter(data, m, MPI_INT, data, m, MPI_INT,  0, MPI_COMM_WORLD); 

/*-------------------- calculate subsums                      */
  locsum = 0.0;
  for (j=0; j<m; j++)
    locsum += data[j] ;
  printf("I got a sum of %d in %d\n", locsum, myid);
/*-------------------- calculate sum of subsums (reduction)--- */
  MPI_Reduce(&locsum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid == 0){ 
    printf("  ***** The total parallel computed sum is  %d \n", sum);
    printf("  ***** The actual sum is                   %d \n", sum0);
  }
  MPI_Finalize(); 
  return(0);
}


