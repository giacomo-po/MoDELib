/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BESTRATIONALAPPROXIMATION_H_
#define model_BESTRATIONALAPPROXIMATION_H_

#include <math.h>

// https://rosettacode.org/wiki/Convert_decimal_number_to_rational#C

namespace model
{
    
    class BestRationalApproximation
    {
    public:
        typedef int64_t LongIntType;
        
    private:
        static std::pair<LongIntType,LongIntType> rat_approx(double f, LongIntType md)
        {
            if(fabs(f)<1.0/md)
            {
                return std::pair<LongIntType,LongIntType>(0,1);
            }

            
            /*  a: continued fraction coefficients. */
            LongIntType a, h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
            LongIntType x, d, n = 1;
            long int i, neg = 0;
            
            if (md <= 1)
            { // *denom = 1; *num = (LongIntType) f; return;
                return std::pair<LongIntType,LongIntType>(f,1);
            }
            
            if (f < 0.0)
            {
                neg = 1;
                f = -f;
            }
            
            
            while (f != floor(f))
            {
                n <<= 1;
                f *= 2;
            }

            d = f;
            
            
            /* continued fraction and check denominator each step */
            for (i = 0; i < 64; i++)
            {
                
                a = n ? d / n : 0;
                if (i && !a) break;
                
                x = d;
                d = n;
                

                n = x % n;
                
                x = a;
                if (k[1] * a + k[0] >= md)
                {
                    x = (md - k[0]) / k[1];
                    if (x * 2 >= a || k[1] >= md)
                        i = 65;
                    else
                        break;
                }
                
                h[2] = x * h[1] + h[0]; h[0] = h[1]; h[1] = h[2];
                k[2] = x * k[1] + k[0]; k[0] = k[1]; k[1] = k[2];
            }
            return std::pair<LongIntType,LongIntType>(neg ? -h[1] : h[1],k[1]);
        }
        
        const std::pair<LongIntType,LongIntType> intPair;
        
    public:
        
        const LongIntType& num;
        const LongIntType& den;
        
        BestRationalApproximation(const double& f, const LongIntType& md) :
        /* init list */ intPair(rat_approx(f,md)),
        /* init list */ num(intPair.first),
        /* init list */ den(intPair.second)
        {
            
        }
        
    };
    
}
#endif
