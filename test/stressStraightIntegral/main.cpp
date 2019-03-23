#include <iostream>
#include <chrono>
#define _MODEL_NON_SINGULAR_DD_ 1
#include <model/DislocationDynamics/ElasticFields/StressStraight.h>
#include <model/DislocationDynamics/ElasticFields/DislocationStress.h>

using namespace model;

template<typename ScalarType>
void timeStressStraight(const int& N)
{
    Eigen::Matrix<ScalarType,3,1> P0(Eigen::Matrix<ScalarType,3,1>::Random());
    Eigen::Matrix<ScalarType,3,1> P1(Eigen::Matrix<ScalarType,3,1>::Random());
    Eigen::Matrix<ScalarType,3,1> b(Eigen::Matrix<ScalarType,3,1>::Random());
    Eigen::Matrix<ScalarType,3,1> x(Eigen::Matrix<ScalarType,3,1>::Random());

    auto t0 = std::chrono::system_clock::now();
    StressStraight<3,ScalarType> ss(P0,P1,b);
    for(int n=0;n<N;++n)
    {
        ss.stress(x);
    }
    std::cout<<std::chrono::duration<double>(std::chrono::system_clock::now()-t0).count()<<std::endl;
//    std::cout<<StressStraight<3,ScalarType>(P0,P1,b).stress(x)<<std::endl;
}


int main(int argc, char** argv)
{
    
    // Take meshID as a user input
    int N(1000);
    if (argc>1)
    {
        N=atoi(argv[1]);
    }
    

    timeStressStraight<double>(N);
//    timeStressStraight<float>(N);


    return 0;
}
