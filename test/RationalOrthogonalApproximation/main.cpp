#include <iostream>
#include <model/LatticeMath/RationalOrthogonalApproximation.h>



int main()
{
    srand((unsigned int) time(0));
    
    for(int k=0;k<1000;++k)
    {
        Eigen::Matrix<double,3,1> a=Eigen::Matrix<double,3,1>::Random();
        double theta=Eigen::Matrix<double,1,1>::Random()(0,0)*M_PI;
        
        Eigen::Matrix<double,3,3> Q(Eigen::AngleAxisd(theta, a.normalized()).toRotationMatrix());
        
        model::RationalOrthogonalApproximation rma(theta,a,10000);
        
        Eigen::Matrix<double,3,3> R(rma.m.template cast<double>()/rma.den);
        
        std::cout<<(R*R.transpose()-Eigen::Matrix<double,3,3>::Identity()).norm()/Eigen::Matrix<double,3,3>::Identity().norm()<<" ";
        std::cout<<(R-Q).norm()/Q.norm()<<std::endl;
    }
    
    return 0;
}
