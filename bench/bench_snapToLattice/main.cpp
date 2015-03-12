#include <iostream>
#include <chrono>
#include <model/DislocationDynamics/SnapToLattice.h>
#include <LatticeVector.h>

using namespace model;


int main()
{

    SnapToLattice<3> s;
    
    Eigen::Matrix<double,3,1> d1;
    d1<<2.0,2,3;
    
    Eigen::Matrix<long int,3,1> i2;
    i2<<3,5,1;
    
    LatticeVector<3> v1(d1);
    
    LatticeVector<3> v2(i2);
//    v2<<3,5,1;
    
    LatticeVector<3> v3(v1+v2);

    LatticeVector<3> v4(v1.cross(v2));

    
    
    std::cout<<v3.transpose()<<std::endl;

    std::cout<<v4.transpose()<<std::endl;

    
//    std::cout<<s.gcd(-2,-4,6)<<std::endl;
    
    return 0;
}