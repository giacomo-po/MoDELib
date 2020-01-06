#include <iostream>
#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadratureDynamic.h>
#include <model/Quadrature/QuadPowDynamic.h>

using namespace model;

int main (int argc, char * const argv[])
{
    
    QuadratureDynamic<1,GaussLegendre,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> qd;
    
    const int N=8;
    
    for (int k=0;k<N;++k)
    {
        std::cout<<qd.abscissa(N,k)<<std::endl;//<<" vs "<<Quadrature<1,N,GaussLegendre>::abscissa(k)<<std::endl;
    }
    
    for(const auto& order : qd.orderSet)
    {
        std::cout<<order<<std::endl;
    }
    
    std::cout<<qd.lowerOrder(-7)<<std::endl;
    std::cout<<qd.lowerOrder(6)<<std::endl;
    std::cout<<qd.lowerOrder(150)<<std::endl;
    std::cout<<qd.lowerOrder(3.5)<<std::endl;
    std::cout<<qd.lowerOrder(2048)<<std::endl;
    
    
    QuadPowDynamic<3,GaussLegendre,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> qpd;
    
    std::cout<<qpd.uPow(7)<<std::endl;

    
    std::cout<<qpd.duPow(7)<<std::endl;

    
    std::cout<<qpd.dduPow(7)<<std::endl;

    return 0;
}
