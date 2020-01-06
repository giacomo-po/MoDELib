// /usr/local/cuda/bin/nvcc main.cu -o add_cuda
// ln -s /usr/local/cuda/include/crt/math_functions.hpp /usr/local/cuda/include/math_functions.hpp
// https://devblogs.nvidia.com/easy-introduction-cuda-c-and-c/
// https://stackoverflow.com/questions/9985912/how-do-i-choose-grid-and-block-dimensions-for-cuda-kernels
// The number of threads per block should be a round multiple of the warp size, which is 32 on all current hardware.
// reference for cudaDeviceProp
// https://docs.nvidia.com/cuda/cuda-runtime-api/structcudaDeviceProp.html
#include <iostream>
#include <math.h>

#include <Eigen/Dense>

// Kernel function to add the elements of two arrays
__global__
void add(int n, float *x, float *y)
{
  for (int i = 0; i < n; i++)
    y[i] = x[i] + y[i];
}

int getSPcores(cudaDeviceProp devProp)
{
int cores = 0;
int mp = devProp.multiProcessorCount;
switch (devProp.major){
case 2: // Fermi
if (devProp.minor == 1) cores = mp * 48;
else cores = mp * 32;
break;
case 3: // Kepler
cores = mp * 192;
break;
case 5: // Maxwell
cores = mp * 128;
break;
case 6: // Pascal
if (devProp.minor == 1) cores = mp * 128;
else if (devProp.minor == 0) cores = mp * 64;
else printf("Unknown device type\n");
break;
case 7: // Volta
if (devProp.minor == 0) cores = mp * 64;
else printf("Unknown device type\n");
break;
default:
printf("Unknown device type\n");
break;
}
return cores;
}

int main(void)
{

    int nDevices;
    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
    prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
    prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n",
    2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);

    printf("  # Warp size: %d\n",prop.warpSize);

    printf("  # cores: %d\n",getSPcores(prop));

    printf("  # streaming processor (SP) units: %d\n",prop.multiProcessorCount);

    printf("  # max thread per block: %d\n",prop.maxThreadsPerBlock);

    printf("  # max thread per SP: %d\n",prop.maxThreadsPerMultiProcessor);


}



  int N = 1<<20;
  float *x, *y;

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&x, N*sizeof(float));
  cudaMallocManaged(&y, N*sizeof(float));

  // initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  // Run kernel on 1M elements on the GPU
  add<<<1, 1>>>(N, x, y);

  // Wait for GPU to finish before accessing on host
  cudaDeviceSynchronize();

  // Check for errors (all values should be 3.0f)
  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = fmax(maxError, fabs(y[i]-3.0f));
  std::cout << "Max error: " << maxError << std::endl;

  // Free memory
  cudaFree(x);
  cudaFree(y);
  
  return 0;
}
