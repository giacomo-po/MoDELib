// /usr/local/cuda/bin/nvcc main.cu -o ddcuda -I. -I/usr/local/include

// ln -s /usr/local/cuda/include/crt/math_functions.hpp /usr/local/cuda/include/math_functions.hpp
// https://devblogs.nvidia.com/easy-introduction-cuda-c-and-c/
// https://stackoverflow.com/questions/9985912/how-do-i-choose-grid-and-block-dimensions-for-cuda-kernels
// The number of threads per block should be a round multiple of the warp size, which is 32 on all current hardware.
// reference for cudaDeviceProp
// https://docs.nvidia.com/cuda/cuda-runtime-api/structcudaDeviceProp.html

// Speed problem, last post
// https://stackoverflow.com/questions/10710867/cuda-speedup-for-simple-calculations

#ifndef model_DDcuda_cpp_
#define model_DDcuda_cpp_

#include <DDcuda.h>
#include <iostream>
#include <Eigen/Dense>

namespace model
{
    
    __global__
    void saxpy(int n, float a, float *x, float *y)
    {
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if (i < n) y[i] = a*x[i] + y[i];
    }
    
    __global__
    void dot(int n, Eigen::Vector3d *x, double *y)
    {
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if (i < n) y[i] = x[i].dot(x[i]);
//        if (i < n) y[i] = i;

    }

//    __global__
//    void dot(int n, Eigen::Vector3d *x, double *y)
//    {
////        for(int i=0;i<n;++i)
////            y[i] = x[i].dot(x[i]);
//        if (i < n)
//            y[i] = x[i].dot(x[i]);
//        //        if (i < n) y[i] = i;
//        
//    }
    
    int getSPcores(const cudaDeviceProp& devProp)
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
    
    DDcuda::DDcuda()
    {
        std::cout<<"querying CUDA"<<std::endl;
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
        
        
    }
    
    template <typename T1,typename T2>
    void DDcuda::computeStressAtQuadrature(const std::vector<T1>& x,
                                    std::vector<T2>& y)
    {
    
        // Copy data do device
        T1* d_x; // x on device
        T2* d_y; // y on device
        
        cudaMalloc(&d_x, x.size()*sizeof(T1));
        cudaMalloc(&d_y, y.size()*sizeof(T2));
        
        cudaMemcpy(d_x, x.data(), x.size()*sizeof(T1), cudaMemcpyHostToDevice);
        cudaMemcpy(d_y, y.data(), y.size()*sizeof(T2), cudaMemcpyHostToDevice);
    
        // Call kernel
//        saxpy<<<(x.size()+255)/256, 256>>>(x.size(), 2.0f, d_x, d_y);
        dot<<<(x.size()+1023)/1024, 1024>>>(x.size(), d_x, d_y);

        cudaThreadSynchronize();

        // Copy back to host
        cudaMemcpy(y.data(), d_y, y.size()*sizeof(T2), cudaMemcpyDeviceToHost);

        // Free device memory
        cudaFree(d_x);
        cudaFree(d_y);

        // Check error
//        double maxError = 0.0f;
//        for (const auto& val : y)
//        {
//            std::cout<<val<<std::endl;
//            maxError = std::max(maxError, std::fabs(val-1.0f));
//        }
//        printf("Max error: %f\n", maxError);
        
    }

}
#endif


