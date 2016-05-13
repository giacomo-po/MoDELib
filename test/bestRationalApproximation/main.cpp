/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
//#include <fstream>
//#include <deque>
#include <model/Math/BestRationalApproximation.h>


using namespace model;

int main()
{
    
    const double f=0.333333;
    BestRationalApproximation bra(f,100);
    
    std::cout<<"num="<<bra.num<<std::endl;
    std::cout<<"den="<<bra.den<<std::endl;
    const double r=double(bra.num)/bra.den;
    std::cout<<"f="<<f<<std::endl;
    std::cout<<"r="<<r<<std::endl;

    //    int i;
//    int64_t d, n;
//    double f;
//    
//    printf("f = %16.14f\n", f = 1.0/7);
//    for (i = 1; i <= 20000000; i *= 16) {
//        printf("denom <= %d: ", i);
//        rat_approx(f, i, &n, &d);
//        printf("%lld/%lld\n", n, d);
//    }
//    
//    printf("\nf = %16.14f\n", f = atan2(1,1) * 4);
//    for (i = 1; i <= 20000000; i *= 16) {
//        printf("denom <= %d: ", i);
//        rat_approx(f, i, &n, &d);
//        printf("%lld/%lld\n", n, d);
//    }
    
    return 0;
}


