/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_Rational_h_
#define model_Rational_h_

#include <Eigen/Core>
#include <IntegerMath.h>

namespace model
{
    template<typename IntScalarType>
    struct Rational
    {
        
        IntScalarType n; // numerator
        IntScalarType d; // denominator

        Rational() :
        /* init */ n(0),
        /* init */ d(1)
        {

            
        }
        
        Rational(const IntScalarType& n_in) :
        /* init */ n(n_in),
        /* init */ d(1)
        {
            
            
        }
        
        Rational(const IntScalarType& n_in, const IntScalarType& d_in) :
        /* init */ n(IntegerMath<IntScalarType>::sgn(d_in)*n_in/IntegerMath<IntScalarType>::gcd(abs(n_in),abs(d_in))),
        /* init */ d(IntegerMath<IntScalarType>::sgn(d_in)*d_in/IntegerMath<IntScalarType>::gcd(abs(n_in),abs(d_in)))
        {

            
        }
        
        double asDouble() const
        {
            return double(n)/double(d);
        }
        
        bool operator==(const IntScalarType& other) const
        {
            switch (d)
            {
                case 0:
                    return false;
                    break;
                    
                case 1:
                    return n==other;
                    break;
                    
                default:
                    return n==0 && other==0;
                    break;
            }
        }
        
        bool operator!=(const IntScalarType& other) const
        {
            return !(*this==other);
        }
        
        Rational operator*(const Rational& r2) const
        {
            return Rational(n*r2.n,d*r2.d);
        }
        
        Rational operator/(const Rational& r2) const
        {
            return Rational(n*r2.d,d*r2.n);
        }
        
        Rational operator+(const Rational& r2) const
        {
            return Rational(n*r2.d+r2.n*d,d*r2.d);
        }
        
        Rational operator-(const Rational& r2) const
        {
            return Rational(n*r2.d-r2.n*d,d*r2.d);
        }
        
        Rational operator*(const IntScalarType& i) const
        {
            return Rational(n*i,d);
        }
        
        Rational operator/(const IntScalarType& i) const
        {
            return Rational(n,d*i);
        }
        
        Rational operator+(const IntScalarType& i) const
        {
            return *this+Rational(i,1);
        }
        
        Rational operator-(const IntScalarType& i) const
        {
            return *this-Rational(i,1);
        }
        
        friend std::ostream& operator << (std::ostream& os, const Rational& r)
        {
            os<<r.n<<"/"<<r.d;
            return os;
        }
        
    };
    

    
    
} // end namespace


namespace Eigen
{
    template<typename IntScalarType> struct NumTraits<model::Rational<IntScalarType>>
    : NumTraits<IntScalarType> // permits to get the epsilon, dummy_precision, lowest, highest functions
    {
        typedef long int Real;
        typedef long int NonInteger;
        typedef long int Nested;
        enum
        {
            IsComplex = 0,
            IsInteger = 1,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = int(NumTraits<IntScalarType>::ReadCost),
            AddCost = int(NumTraits<IntScalarType>::AddCost),
            MulCost = int(NumTraits<IntScalarType>::MulCost)
        };
    };
}


#endif
