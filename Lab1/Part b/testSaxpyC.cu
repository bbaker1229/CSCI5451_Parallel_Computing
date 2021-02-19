#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>

#define NITER 100

__global__ void saxpy_par(unsigned long int n, float a, float *x, float *y) {
  unsigned long int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n)
    y[i] += a*x[i];
}

float saxpy_check(unsigned long int n, float a, float *x, float *y, float *z) {
  // a, x, y == original data for saxpy
  // z = result found -- with which to compare.
  float s=0.0, t=0.0;
  for (unsigned long int i=0; i<n; i++) {
    y[i] += a * x[i];
    s += (y[i] - z[i]) * (y[i] - z[i]);
    t += z[i] * z[i];
  }
  if (t == 0.0) return (-1);
  else
    return(sqrt(s/t));
}

double wctime()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

int main(int argc, char *argv[]) {
  unsigned long int n = 8 * 1024 * 1024;
  unsigned long int i;
  float a, value, valresult, nops;
  float *x, *y, *z, *xg, *yg;
  double t1[NITER];
  double avgt;
  a = 1.0;
  //a = a / (float) NITER;
  x = (float*) malloc(n*sizeof(float));
  y = (float*) malloc(n*sizeof(float));
  z = (float*) malloc(n*sizeof(float));
  for(i=0; i < n; i++) {
    value = (float)rand() / (float)RAND_MAX;
    x[i] = value;
    value = (float)rand() / (float)RAND_MAX;
    y[i] = value;
    value = (float)rand() / (float)RAND_MAX;
    z[i] = value;
  }
  for (unsigned long int vecLen = 2048; vecLen <= n; vecLen *= 2) {
    cudaMalloc(&xg, vecLen*sizeof(float));
    cudaMalloc(&yg, vecLen*sizeof(float));
    //cudaMemcpy(xg, x, vecLen*sizeof(float), cudaMemcpyHostToDevice);
    //cudaMemcpy(yg, y, vecLen*sizeof(float), cudaMemcpyHostToDevice);
    dim3 dimGrid(1024 * 8);
    dim3 dimBlock(1024);
    //t1 = wctime();
    for(int iter=0; iter<NITER; iter++){
      t1[iter] = wctime();
      cudaMemcpy(xg, x, vecLen*sizeof(float), cudaMemcpyHostToDevice);
      cudaMemcpy(yg, y, vecLen*sizeof(float), cudaMemcpyHostToDevice);
      saxpy_par<<<dimGrid, dimBlock>>>(vecLen, a, xg, yg);
      cudaMemcpy(z, yg, vecLen*sizeof(float), cudaMemcpyDeviceToHost);
      t1[iter] = (wctime() - t1[iter])*1.e+06;
      //valresult = saxpy_check(vecLen, a, x, y, z);
    }
    //t1 = (wctime() - t1)*1.e+06;
    for(i=0; i<NITER; i++)
      avgt += t1[i];
    avgt /= (float) NITER;
    valresult = saxpy_check(vecLen, a, x, y, z);
    nops = (float) vecLen * 2;
    printf("** vecLen = %7.0lu, Mflops = %10.2lf  err = %2.2e\n", vecLen, nops/avgt, valresult);
    cudaFree(xg);
    cudaFree(yg);
  }
  free(x);
  free(y);
  free(z);
  return(0);
}
