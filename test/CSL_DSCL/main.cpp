#include <iostream>
#include <chrono>
#include <deque>
#include <Eigen/Dense>
#include <model/Math/SmithDecomposition.h>
#include <model/LatticeMath/CSL.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>

using namespace model;


void fillMap(std::deque<Eigen::Matrix<double,3,3>>& gbMap)
{
// Rotation matrices from Grimmer, Acta Cryst. (1974). A30, 197, table 1

//    gbMap.push_back((Eigen::Matrix<double,3,3>()<<2,-1,2,2,2,-1,-1,2,2).finished()/3);
    gbMap.push_back((Eigen::Matrix<double,3,3>()<<5,0,0,0,4,-3,0,3,4).finished()/5);

}

int main()
{

    const int dim=3;
    Lattice<dim> A;
    Lattice<dim> B;
    
    
    std::deque<Eigen::Matrix<double,dim,dim>> gbMap;
    fillMap(gbMap);
    
    for (const auto& R : gbMap)
    {
    
        A.setLatticeBasis(FCC::getLatticeBasis<dim>());
        B.setLatticeBasis(R*FCC::getLatticeBasis<dim>());
        CSL<dim> csl(A,B);
        
        
    }
    
//    const int N=4;
//    Eigen::Matrix<long long int,N,N> A=Eigen::Matrix<long long int,N,N>::Random()/100000000;
//
////    const int N=2;
////    Eigen::Matrix<long long int,N,N> A;
////    A<<4,-3,3,4;
//
//    
//    std::cout<<"matrix A="<<std::endl<<A<<std::endl<<std::endl;
//    
//    SmithDecomposition<N> sd(A);
//    
//    std::cout<<"Smith decomposition U*A*V=D"<<std::endl;
//    std::cout<<"unimodular matrix U="<<std::endl<<sd.matrixU()<<std::endl;
//    std::cout<<"unimodular matrix V="<<std::endl<<sd.matrixV()<<std::endl;
//    std::cout<<"diagonal matrix D="<<std::endl<<sd.matrixD()<<std::endl;
//
//    std::cout<<"Smith decomposition X*D*Y=A"<<std::endl;
//    std::cout<<"unimodular matrix X="<<std::endl<<sd.matrixX()<<std::endl;
//    std::cout<<"unimodular matrix Y="<<std::endl<<sd.matrixY()<<std::endl;
//    std::cout<<"diagonal matrix D="<<std::endl<<sd.matrixD()<<std::endl;

    
    return 0;
}
