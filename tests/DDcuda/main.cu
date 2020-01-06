// /usr/local/cuda/bin/nvcc main.cu -o DDcuda
// ln -s /usr/local/cuda/include/crt/math_functions.hpp /usr/local/cuda/include/math_functions.hpp
// https://devblogs.nvidia.com/easy-introduction-cuda-c-and-c/
// https://stackoverflow.com/questions/9985912/how-do-i-choose-grid-and-block-dimensions-for-cuda-kernels
// The number of threads per block should be a round multiple of the warp size, which is 32 on all current hardware.
// reference for cudaDeviceProp
// https://docs.nvidia.com/cuda/cuda-runtime-api/structcudaDeviceProp.html

#include <chrono>
#include <iomanip>
#include <omp.h>


#include <DDcuda.cpp>


int main(void)
{

    model::DDcuda ddc;

std::vector<Eigen::Vector3d> x;

 std::vector<double> y;
std::vector<double> z;

const int N=1e8;
x.reserve(N);
y.reserve(N);
z.reserve(N);

for (size_t n=0;n<N;++n)
{
x.push_back(Eigen::Vector3d::Random().normalized());
y.push_back(0);
z.push_back(0);

}

// CUDA code
std::cout<<"CUDA"<<std::endl;
const auto t0= std::chrono::system_clock::now();
ddc.computeStressAtQuadrature(x,y);
std::cout<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<std::endl;

// serial code
std::cout<<"serial"<<std::endl;
const auto t1= std::chrono::system_clock::now();
for (size_t n=0;n<N;++n)
{
    z[n]=x[n].dot(x[n]);
}
std::cout<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<std::endl;

// omp code
//std::cout<<"serial"<<std::endl;
//const auto t2= std::chrono::system_clock::now();
//#pragma omp parallel for
//for (int n=0;n<10000000;++n)
//{
//z[n]=x[n].dot(x[n]);
//}
//std::cout<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<std::endl;


  return 0;
}
