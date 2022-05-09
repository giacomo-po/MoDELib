/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_IntegerMath_h_
#define model_IntegerMath_h_


namespace model
{
    template <typename IntScalarType>
    struct IntegerMath
    {
        
        static IntScalarType sgn(const IntScalarType& a)
        {
            return a<0? -1 : 1;
        }
        
        static IntScalarType gcd(const IntScalarType &a, const IntScalarType &b)
        {
            const IntScalarType absA(abs(a));
            const IntScalarType absB(abs(b));
            return absB > 0 ? gcd(absB, absA % absB) : (absA > 0 ? absA : 1);
        }
        
        template<typename ArrayType>
        static IntScalarType gcd(const ArrayType& a)
        {
            
            const size_t aSize(a.size());
            switch (aSize)
            {
                case 0:
                {
                    throw std::runtime_error("gcd: array size is zero\n");
                    return 0;
                    break;
                }
                    
                case 1:
                {
                    return a[0];
                    break;
                }
                    
                case 2:
                {
                    return gcd(a[0],a[1]);
                    break;
                }
                    
                default:
                {
                    IntScalarType temp(a[0]);
                    for(size_t k=1;k<aSize;++k)
                    {
                        temp=gcd(temp,a[k]);
                    }
                    return temp;
                    break;
                }
            }
        }
        
        static IntScalarType lcm(const IntScalarType &a, const IntScalarType &b)
        {
            return a * b / gcd(a, b);
        }
        
    };
}
#endif
