#include <iostream>
#include <chrono>
#include <iomanip>
#include <model/Geometry/SegmentSegmentIntersection.h>

using namespace model;
int main()
{

    const int dim=3;
    typedef typename SegmentSegmentIntersection<dim>::VectorDimD VectorDimD;
    
// Test cases from V. LUMELSKY, Information Processing Letters 21 (1985) 55-61.
    
    VectorDimD A;
    A<<0, 0, 0;
    VectorDimD B;
    B<<1, 2, 1;
    VectorDimD C;
    C<<1, 0, 0;
    VectorDimD D;
    D<<2, 1, 0;
    
    // Results test
    SegmentSegmentIntersection<dim> ssi(A,B,C,D);
    std::cout<<"t="<<ssi.t<<" (expected 0.167)"<<std::endl;
    std::cout<<"u="<<ssi.u<<" (expected 0)"<<std::endl;
    std::cout<<"dMin="<<ssi.dMin<<" (expected 0.9129)"<<std::endl;

    // Speed test
    const auto t0= std::chrono::system_clock::now();
    const int nRep=1000000;
    for(int i=0;i<nRep;++i)
    {
        SegmentSegmentIntersection<dim> ssi1(A,B,C,D);
    }
    const auto dt=(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count();
    std::cout<<"average computation time="<<std::setprecision(3)<<std::scientific<<dt/nRep<<" sec"<<std::endl;
    
    return 0;
}
