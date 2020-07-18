#include <iostream>
#include <CompareVectorsByComponent.h>

using namespace model;


int main()
{
    Eigen::Matrix<double,3,1> a(1.714870967880501e+03,4.339526089740382e-13,2.690404278911689e+03);
    Eigen::Matrix<double,3,1> b(1.714870967880502e+03,-2.273736754432321e-12,2.690404278911691e+03);
                                
    CompareVectorsByComponent<double,3,float> comp;
    
    std::cout<<comp(a,b)<<std::endl;

    return 0;
}
