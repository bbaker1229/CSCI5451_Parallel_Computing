#include <stdio.h>

#define NITER 100

__global__ void saxpy_par(unsigned long int n, float a, float *x, float *y) {
  //
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

int main(int argc, char *argv[]) {
  unsigned long int n = 8 * 1024 * 1024;
  unsigned long int i;//, maxval;
  float a, value, valresult;
  float *x, *y, *z, *xg, *yg;//, *zg, *xtmp, *ytmp, *ztmp;
  a = 1.0;
  a = a / (float) NITER;
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
  //printf("Entering loop:\n");
  for (unsigned long int vecLen = 2048; vecLen <= n; vecLen *= 2) {
    //printf("vecLen = %lu\n", vecLen);
    /*xtmp = (float*) malloc(vecLen*sizeof(float));
    ytmp = (float*) malloc(vecLen*sizeof(float));
    ztmp = (float*) malloc(vecLen*sizeof(float));
    if(2*vecLen > n)
      maxval = n;
    else
      maxval = 2*vecLen;
    for(i=vecLen; i <= maxval; i++) {
      xtmp[i-vecLen] = x[i];
      ytmp[i-vecLen] = y[i];
      ztmp[i-vecLen] = z[i];
    }*/
    //printf("Copied data\n");
    cudaMalloc(&xg, vecLen*sizeof(float));
    cudaMalloc(&yg, vecLen*sizeof(float));
    //cudaMalloc(&zg, vecLen*sizeof(float));
    cudaMemcpy(xg, x, vecLen*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(yg, y, vecLen*sizeof(float), cudaMemcpyHostToDevice);
    //cudaMemcpy(zg, z, vecLen*sizeof(float), cudaMemcpyHostToDevice);
    dim3 dimGrid(1024 * 8);
    dim3 dimBlock(1024);
    for(int iter=0; iter<NITER; iter++){
      saxpy_par<<<dimGrid, dimBlock>>>(vecLen, a, xg, yg);
      cudaMemcpy(z, yg, vecLen*sizeof(float), cudaMemcpyDeviceToHost);
      valresult = saxpy_check(vecLen, a, x, y, z);
    }
    printf("** vecLen = %7.0lu, Mflops = xx.dd  err = %2.2e\n", vecLen, valresult);
  }
  //printf("working\n");
  return(0);
}
