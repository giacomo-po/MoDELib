#include <iostream>
#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadratureDynamic.h>

using namespace model;

int main (int argc, char * const argv[])
{
    
    QuadratureDynamic<1,GaussLegendre,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024,2048> qd;
	
    const int N=9;
    
    for (int k=0;k<N;++k)
    {
        std::cout<<qd.abscissa(N,k)<<" vs "<<Quadrature<1,N,GaussLegendre>::abscissa(k)<<std::endl;
    }
    
    return 0;
}
