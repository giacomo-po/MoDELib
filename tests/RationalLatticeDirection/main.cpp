#include <iostream>
#include <FCClattice.h>

using namespace model;

int main()
{
    
    FCClattice<3> lat(Eigen::Matrix<double,3,3>::Identity());
    
    Eigen::Matrix<double,3,1> v;
    v<<sqrt(2.0)/2.0,sqrt(2.0)/2.0,0.0;
    const auto lv(lat.latticeVector(v));
    
    v<<sqrt(2.0)/4.0,-sqrt(2.0)/4.0,0.0;
    const auto rld(lat.rationalLatticeDirection(v));
    std::cout<<rld.rat<<std::endl;
    std::cout<<rld.dir.cartesian().transpose()<<std::endl;

    
    return 0;
}
