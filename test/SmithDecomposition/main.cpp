#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <model/Math/SmithDecomposition.h>

using namespace model;



int main()
{

    const int N=4;
    Eigen::Matrix<long long int,N,N> A=Eigen::Matrix<long long int,N,N>::Random()/100000000;

//    const int N=2;
//    Eigen::Matrix<long long int,N,N> A;
//    A<<4,-3,3,4;

    
    std::cout<<"matrix A="<<std::endl<<A<<std::endl<<std::endl;
    
    SmithDecomposition<N> sd(A);
    
    std::cout<<"Smith decomposition U*A*V=D"<<std::endl;
    std::cout<<"unimodular matrix U="<<std::endl<<sd.matrixU()<<std::endl;
    std::cout<<"unimodular matrix V="<<std::endl<<sd.matrixV()<<std::endl;
    std::cout<<"diagonal matrix D="<<std::endl<<sd.matrixD()<<std::endl;

    std::cout<<"Smith decomposition X*D*Y=A"<<std::endl;
    std::cout<<"unimodular matrix X="<<std::endl<<sd.matrixX()<<std::endl;
    std::cout<<"unimodular matrix Y="<<std::endl<<sd.matrixY()<<std::endl;
    std::cout<<"diagonal matrix D="<<std::endl<<sd.matrixD()<<std::endl;

    
    return 0;
}
