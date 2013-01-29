

// /opt/local/bin/g++-mp-4.7 bench_quadratureDynamic.cpp -o bench_quadratureDynamic -H -O2 -std=c++11 -I../ -I/usr/local/include

#include <iostream>
#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadratureDynamic.h>

using namespace model;

int main (int argc, char * const argv[]) {
    
    const int qMax=200;
    QuadratureDynamic<1,qMax> qd;
	
    for (int q=1;q<=qMax;++q)
    {
        std::cout<<"QUADRATURE ORDER="<<q<<std::endl;
        for (int k=0;k<q;++k)
        {
            std::cout<<qd.abscissa(q,k)<<std::endl;
        }
        std::cout<<std::endl;
    }
    
    
//    Quadrature<1,1> q1;
//    Quadrature<1,2> q2;
//    //Quadrature<1,3> q3;
//    Quadrature<1,4> q4;
//    //Quadrature<1,5> q5;
//    //Quadrature<1,6> q6;
//    //Quadrature<1,7> q7;
//    Quadrature<1,8> q8;
//    //Quadrature<1,9> q9;
//    //Quadrature<1,10> q10;
//    //Quadrature<1,11> q11;
//    //Quadrature<1,12> q12;
//    //Quadrature<1,13> q13;
//    //Quadrature<1,14> q14;
//    //Quadrature<1,15> q15;
//    Quadrature<1,16> q16;
//    //Quadrature<1,17> q17;
//    //Quadrature<1,18> q18;
//    //Quadrature<1,19> q19;
//    //Quadrature<1,20> q20;
//
//    std::cout<<q1.abscissa(0)<<std::endl;
//    std::cout<<q2.abscissa(0)<<std::endl;
//  //  std::cout<<q3.abscissa(0)<<std::endl;
//    std::cout<<q4.abscissa(0)<<std::endl;
//  //  std::cout<<q5.abscissa(0)<<std::endl;
//  //  std::cout<<q6.abscissa(0)<<std::endl;
//  //  std::cout<<q7.abscissa(0)<<std::endl;
//        std::cout<<q8.abscissa(0)<<std::endl;
//   //     std::cout<<q9.abscissa(0)<<std::endl;
//   //     std::cout<<q10.abscissa(0)<<std::endl;
//   //     std::cout<<q11.abscissa(0)<<std::endl;
//   // std::cout<<q12.abscissa(0)<<std::endl;
//   // std::cout<<q13.abscissa(0)<<std::endl;
//   // std::cout<<q14.abscissa(0)<<std::endl;
//   // std::cout<<q15.abscissa(0)<<std::endl;
//    std::cout<<q16.abscissa(0)<<std::endl;
//   // std::cout<<q17.abscissa(0)<<std::endl;
//   // std::cout<<q18.abscissa(0)<<std::endl;
//   // std::cout<<q19.abscissa(0)<<std::endl;
//   // std::cout<<q20.abscissa(0)<<std::endl;
    
    
    return 0;
}
