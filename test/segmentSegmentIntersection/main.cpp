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

    if(true)
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
    
    if(true)
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
    
    if(true)
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
    
    return 0;
}
