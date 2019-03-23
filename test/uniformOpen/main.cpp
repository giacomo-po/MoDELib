#include <iostream>
#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadratureDynamic.h>

using namespace model;
int main()
{

    const size_t qOrder=64;

//    Quadrature<2,qOrder,UniformOpen> quad;
//
//    std::cout<<quad.abscissas.transpose()<<std::endl<<std::endl;
//    std::cout<<quad.weights.transpose()<<std::endl<<std::endl;

    //
//    std::cout<<UniformOpen<2,qOrder>::abcsissasAndWeights().transpose()<<std::endl;
    
    QuadratureDynamic<2,UniformOpen,1,4,9,16,25,36,49,64,81,100> qd;
//    QuadratureDynamic<2,UniformOpen,4,49,16> qd;

//    std::cout<<qd.abscissa(49,3).transpose()<<std::endl;//<<" vs "<<Quadrature<1,N,GaussLegendre>::abscissa(k)<<std::endl;

    
//    std::cout<<    Quadrature<2,49,UniformOpen>::abscissa(3)<<std::endl;
    
    for (int N=1;N<=10;++N)
    {
        for (int k=0;k<N*N;++k)
        {
            std::cout<<qd.abscissa(N*N,k).transpose()<<"\n";//<<" vs "<<Quadrature<1,N,GaussLegendre>::abscissa(k)<<std::endl;
        }
            std::cout<<"\n";//<<" vs "<<Quadrature<1,N,GaussLegendre>::abscissa(k)<<std::endl;
    }
    
    
    return 0;
}