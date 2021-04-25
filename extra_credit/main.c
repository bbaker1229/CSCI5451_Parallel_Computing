// Compile using: 
// mpicc -o main.ex main.c -lm

// Execute using:
// mpirun -np x main.ex
// where x is the number of processors and must be a square number

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int main (int argc, char **argv){
 // Variable Declaration
 int myid, procs;
 MPI_Comm comm2D, row_comm, col_comm;
 int dims[2], coords[2], periods[2];
 double ai, bj;  
 int i, j, ierr = 0;
 
 // Initialize MPI
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &procs);
 MPI_Comm_rank(MPI_COMM_WORLD, &myid);

 // Check that the number of processors is a square number
 if (myid == 0) {
   i = (int) sqrt((double)procs);
   if(i * i != procs) {
     printf(" ERROR: the number of procs is not a square.\n");
     ierr = 1;
   }
   j = i;
   dims[0] = i;
   dims[1] = j;
 }

 // Broadcast error and exit if there is
 MPI_Bcast(&ierr, 1, MPI_INT, 0, MPI_COMM_WORLD);
 if (ierr>0) {
   return(ierr);
   MPI_Finalize();
   exit(1);
 }

 // Broadcast the final matrix dimensions
 MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
 i = dims[0] ; 
 j = dims[1];
 if (procs != i*j) {
   if (myid == 0) 
     printf(" ERROR: i*j does not equal nprocs - aborting \n");
   exit(1);
 }

 // Create Cartesian communicator
 periods[0]=periods[1]=0;
 MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&comm2D);

 // Find my coordinates in the 2-D communicator
 MPI_Cart_coords(comm2D,myid,2,coords);
 i = coords[0];
 j = coords[1];

 // Create n x 1 matrix A
 if(j == 0)
   ai = (double) (i+1);

 // Create 1 x n matrix B
 if(i == 0)
   bj = (double) (j+1);

 // Create Cartesian row communicator
 dims[0] = 0;
 dims[1] = 1;
 MPI_Cart_sub(comm2D, dims, &row_comm);

 // Create Catesian column communicator
 dims[0] = 1;
 dims[1] = 0;
 MPI_Cart_sub(comm2D, dims, &col_comm);

 // Broadcast row value of A across rows
 MPI_Bcast(&ai, 1, MPI_DOUBLE, 0, row_comm);

 // Broadcast col value of B across columns
 MPI_Bcast(&bj, 1, MPI_DOUBLE, 0, col_comm);

 // Print results
 printf(" (i, j) (%d %d ) == myid %d: ai = %d, bj = %d, value == %d\n", i, j, myid, (int)ai, (int)bj, (int)(ai*bj));

 MPI_Finalize();
 return(0);
}   

