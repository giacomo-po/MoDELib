#include <iostream>
#include <chrono>
#include <model/DislocationDynamics/SnapToLattice.h>
using namespace model;


int main()
{

    SnapToLattice<3> s;
    
    
    
    std::cout<<s.gcd(-2,-4,6)<<std::endl;
    
    return 0;
}