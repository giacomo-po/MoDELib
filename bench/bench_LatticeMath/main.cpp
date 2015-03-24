#include <iostream>
#include <chrono>
//#include <model/DislocationDynamics/SnapToLattice.h>
#include <model/DislocationDynamics/LatticeMath/Rational.h>

#include <model/DislocationDynamics/LatticeMath/LatticeBase.h>
#include <model/DislocationDynamics/LatticeMath/LatticeVector.h>
//#include <model/DislocationDynamics/LatticeMath/LatticeDirection.h>
//#include <model/DislocationDynamics/LatticeMath/LatticeLine.h>

using namespace model;




int main()
{

    const double a=1.0;
    const double b=0.5*sqrt(2.0)*a;
    
//    SnapToLattice<3> s;
    Eigen::Matrix<double,3,3> A;
    A << 0.0, 1.0, 1.0,
    /**/ 1.0, 0.0, 1.0,
    /**/ 1.0, 1.0, 0.0;
    A*=0.5*a;
    
    LatticeBase<3>::setLatticeMatrix(A);
    
    Eigen::Matrix<double,3,1> d1;
    d1<<1.0,1.0,-1.0;
    LatticeVector<3> v1(d1);

    
    Eigen::Matrix<double,3,1> d2;
    d2<<1.0,0.0,1.0;
    LatticeVector<3> v2(d2);

//
//    Eigen::Matrix<long int,3,1> i2;
//    i2<<-1,0,-5;
//    
    std::cout<<"v1.contra="<<v1.contra().transpose()<<std::endl;
    std::cout<<"v1.cov="<<v1.cov().transpose()<<std::endl;

    Rational r1(2,3); // 2/3
    Rational r2(3,-4); // 2/3
    
    std::cout<<r1*r2<<std::endl;
    std::cout<<r1/r2<<std::endl;
    std::cout<<r1+r2<<std::endl;
    std::cout<<r1-r2<<std::endl;

    std::cout<<r1-2<<std::endl;

    
    //    LatticeVector<3> v1(d1);
//    std::cout<<"v1.contra="<<v1.contra().transpose()<<std::endl;

    
    Eigen::Matrix<Rational,3,1> v3;
    v3<< Rational(2,3) , Rational(4,-5) , Rational(1,2);

    std::cout<<v3<<std::endl;
    
    //
//    LatticeVector<3> v2(i2);
////    v2<<3,5,1;
//    
//    LatticeVector<3> v3(v1+v2);
//
//    LatticeDirection<3> dir1(v1);
//    
////    LatticePoint<3> v4(v1,v2);
//
//    
//    
//    std::cout<<v3.transpose()<<std::endl;
//
//    std::cout<<dir1.transpose()<<std::endl;
//
//    LatticeLine l1(v2,dir1);
    
    
//    std::cout<<s.gcd(-2,-4,6)<<std::endl;
    
    return 0;
}