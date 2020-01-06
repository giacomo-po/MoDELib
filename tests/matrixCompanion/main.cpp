#include <iostream>
#include <MatrixCompanion.h>


int main()
{

    Eigen::Matrix<double,1,5> coeffs;
//    coeffs<<0.0,0.0,2.3,4.2,9.0;
    coeffs<<1.406e-11,-6.245e-7,-9.964e-1,-9.097e3,1.143e5;

    model::MatrixCompanion<4> mc(coeffs);
    for(int k=0;k<4;++k)
    {
        std::cout<<mc.root(k)<<std::endl;
    }
    
    return 0;
}
