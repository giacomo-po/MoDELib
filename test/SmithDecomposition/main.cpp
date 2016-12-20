#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <model/Math/SmithDecomposition.h>

using namespace model;



int main()
{

    const int N=4;
    Eigen::Matrix<int,N,N> A=Eigen::Matrix<int,N,N>::Random()/10000000;
    
//    Eigen::Matrix<int,4,4> T1=SmithDecomposition::signChange<4>(3);
    std::cout<<A<<std::endl<<std::endl;

//    std::cout<<T1*B<<std::endl<<std::endl;
//    std::cout<<B*T1<<std::endl<<std::endl;
//    
//    Eigen::Matrix<int,4,4> T2=SmithDecomposition::swap<4>(2,3);
//        std::cout<<T2*B<<std::endl<<std::endl;
//    std::cout<<B*T2<<std::endl<<std::endl;
//
//    Eigen::Matrix<int,4,4> T3=SmithDecomposition::add<4>(2,3,2);
//    std::cout<<T3*B<<std::endl<<std::endl;
//    std::cout<<B*T3.transpose()<<std::endl<<std::endl;

    
    SmithDecomposition<N> sd(A);
//    std::cout<<sd.matrixD()<<std::endl<<std::endl;

      return 0;
}
