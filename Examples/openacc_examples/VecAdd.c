#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openacc.h>

int main( int argc, char* argv[] )
{

  // Size of vectors
  int n = 10000;

  // Input vectors
  double *a;
  double *b;
  // Output vector
  double *c;

  // Size, in bytes, of each vector
  size_t bytes = n*sizeof(double);

  // Allocate memory for each vector
  a = (double*)malloc(bytes);
  b = (double*)malloc(bytes);
  c = (double*)malloc(bytes);

  printf("Number of devices: %d\n", acc_get_num_devices(acc_device_nvidia));
// Initialize content of input vectors,
// vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
  int i;
  for(i=0; i<n; i++) {
    a[i] = sin(i)*sin(i);
    b[i] = cos(i)*cos(i);
  }
  // sum component wise and save result into vector c
#pragma acc kernels copyin(a[0:n],b[0:n]), copyout(c[0:n])
  for(i=0; i<n; i++) {
    c[i] = a[i] + b[i];
  }

// Sum up vector c and print result divided by n,
// this should equal 1 within error
  double sum = 0.0;
  for(i=0; i<n; i++) {
    sum += c[i];
  }
  sum = sum/n;
  printf("final result: %f\n", sum);

  // Release memory
  free(a);
  free(b);
  free(c);

  return 0;
}
