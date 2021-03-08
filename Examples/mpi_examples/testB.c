#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
/*----------------------------------------
 * this is to compute the exponential of a 
 * number in a pipelined fashion -- Broadcast version 
 * ---------------------------------------*/

int main(int argc, char *argv[]){
  int myPE, East, West, n;
  double x, val;
  MPI_Status stat;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n);
  MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
/*------------------- broadcast x from root and set val to 1   */
  if (myPE == 0) 
    x = 1.0;
  MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  East = myPE+1; West = myPE-1;
  /*-------------------- start */
  if (myPE == n-1) {
    val = 1.0 + x /( (double) myPE+1 ) ;
    MPI_Send(&val, 1, MPI_DOUBLE, West, 0, MPI_COMM_WORLD) ;
  }   else    {
    MPI_Recv (&val, 1, MPI_DOUBLE, East, 0, MPI_COMM_WORLD,&stat); 
    val = 1.0 + val *( x /( (double) myPE+1));
    if (myPE>0) 
      MPI_Send(&val, 1, MPI_DOUBLE, West, 0, MPI_COMM_WORLD) ;
    else
      printf("  ***** Value val at node 0 is  %e \n", val);
  }
  MPI_Finalize(); 
  return(0);
}

