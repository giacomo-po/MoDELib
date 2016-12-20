#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <model/Math/SmithDecomposition.h>

using namespace model;



int main()
{

    const int N=4;
    Eigen::Matrix<long long int,N,N> A=Eigen::Matrix<long long int,N,N>::Random()/100000000;
    
    std::cout<<"matrix A="<<std::endl<<A<<std::endl<<std::endl;
    
    SmithDecomposition<N> sd(A);
    
    std::cout<<"Smith decomposition U*A*V=D"<<std::endl;
    std::cout<<"unimodular matrix U="<<std::endl<<sd.matrixU()<<std::endl;
    std::cout<<"unimodular matrix V="<<std::endl<<sd.matrixV()<<std::endl;
    std::cout<<"diagonal matrix D="<<std::endl<<sd.matrixD()<<std::endl;
      return 0;
}
