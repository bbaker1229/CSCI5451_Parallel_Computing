#include <stdio.h>
  __global__ void helloFromGPU2(){
  int i = threadIdx.x;
  int j = threadIdx.y;
  printf("Hello World-Thread: %d, %d \n",i,j);                                                      
  }

  int main(void) {
  dim3 ThisBlock(4,4);
    helloFromGPU2<<<1,ThisBlock>>>();
    cudaDeviceSynchronize();
    return(0);
  }
