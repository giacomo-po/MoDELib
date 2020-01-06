#include <iostream>
#include <chrono>
#include <iomanip>
#include <model/Geometry/SegmentSegmentDistance.h>

using namespace model;
int main()
{
    
    const int dim=3;
    typedef typename SegmentSegmentDistance<dim>::VectorDim VectorDim;

    VectorDim A;
    VectorDim B;
    VectorDim C;
    VectorDim D;

    if(false)
    {// Test cases from V. LUMELSKY, Information Processing Letters 21 (1985) 55-61.
        A<<0, 0, 0;
        B<<1, 2, 1;
        C<<1, 0, 0;
        D<<2, 1, 0;
        
        SegmentSegmentDistance<dim> ssi(A,B,C,D);
        std::cout<<"t="<<ssi.t<<" (expected 0.167)"<<std::endl;
        std::cout<<"u="<<ssi.u<<" (expected 0)"<<std::endl;
        std::cout<<"dMin="<<ssi.dMin<<" (expected 0.9129)"<<std::endl;
    }
    
    if(false)
    {// Speed test
        const auto t0= std::chrono::system_clock::now();
        const int nRep=1000000;
        for(int i=0;i<nRep;++i)
        {
            SegmentSegmentDistance<dim> ssi(A,B,C,D);
        }
        const auto dt=(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count();
        std::cout<<"average computation time="<<std::setprecision(3)<<std::scientific<<dt/nRep<<" sec"<<std::endl;
        
    }
    
    if(false)
    {// Coincident test
        A<<0, 0, 0;
        B<<1, 0, 0;
        C<<-0.5, 0, 0;
        D<<0.5, 0, 0;
        
        SegmentSegmentDistance<dim> ssi(A,B,C,D);
        std::cout<<"t="<<ssi.t<<std::endl;
        std::cout<<"u="<<ssi.u<<std::endl;
        std::cout<<"dMin="<<ssi.dMin<<std::endl;

    }
    
    if(true)
    {// Degenerate test
        A<<4.508746034732568e+0,2 -5.545084859017030e+02,  5.926782361066371e+02;
        B<<9.301489269216704e+01, -6.182861388340581e+02,  8.867602939553722e+02;
        C=B;
        D<<1.231196733319394e+02, -3.735588426060161e+02,  1.101382809543642e+03;
        
        SegmentSegmentDistance<dim> ssi(A,B,D,D,0.0);
        std::cout<<"t="<<ssi.t<<std::endl;
        std::cout<<"u="<<ssi.u<<std::endl;
        std::cout<<"dMin="<<ssi.dMin<<std::endl;
        
    }
    
    return 0;
}
