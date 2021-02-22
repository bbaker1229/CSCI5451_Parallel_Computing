#include <stdio.h>
#include <cuda.h>
	   
   __global__ void vec_add(float *a, float *b)   {
//  int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = threadIdx.x;
      a[i] = a[i]+b[i];
    }
/*-------------------- main will execute on host */

int main(void){
/* -------------------- host & device arrays*/
  float *y_h, *y_d,  *x_h, *x_d;  
/* -------------------- size of arrays */
  const int N = 500;  
  size_t size = N * sizeof(float);
/* -------------------- Allocate array on host */
  y_h = (float *)malloc(size);        
  x_h = (float *)malloc(size);        
/* -------------------- Allocate array on device */
  cudaMalloc((void **) &y_d, size);   
  cudaMalloc((void **) &x_d, size);   
/* -------------------- Initialize host array & 
                        copy it to device */
  for (int i=0; i<N; i++) {
      y_h[i] = (float) i;
      x_h[i] = (float) (N-i);
  }     
  cudaMemcpy(y_d,y_h,size,cudaMemcpyHostToDevice);
  cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
/* -------------------- invoke kernel on device */
//   int block_size = 4;
//   int n_blocks = (N+block_size-1)/block_size; 
//   vec_add <<< n_blocks, block_size >>> (y_d, N);
    vec_add <<< 1, N >>> (y_d, x_d);    
/*-------------------- retrieve result from device*/
   cudaMemcpy(y_h,y_d,size,cudaMemcpyDeviceToHost);
/* --------------------/ Print 20 first results / */
   for (int i=0; i<20; i++) 
       printf("%d %f\n", i, y_h[i]);
/* -------------------- free memory */ 
   free(y_h); cudaFree(y_d);
 }

