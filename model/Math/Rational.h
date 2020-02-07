/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Rational_h_
#define model_Rational_h_

#include <Eigen/Core>


namespace model
{
    class Rational
    {
        
        /**********************************************************************/
        static int gcd(const size_t& a,const size_t& b)
        {
            return b>0? gcd(b, a % b) : a;
        }
        
        /**********************************************************************/
        static int sgn(const long int& a)
        {
            return a<0? -1 : 1;
        }
        
    public:
        
        long int n; // numerator
        long int d; // denominator

        Rational() :
        /* init */ n(0),
        /* init */ d(1)
        {

            
        }
        
        Rational(const long int& n_in) :
        /* init */ n(n_in),
        /* init */ d(1)
        {
            
            
        }
        
        Rational(const long int& n_in, const long int& d_in) :
        /* init */ n(sgn(d_in)*n_in/gcd(abs(n_in),abs(d_in))),
        /* init */ d(sgn(d_in)*d_in/gcd(abs(n_in),abs(d_in)))
        {

            
        }
        
        double asDouble() const
        {
            return double(n)/double(d);
        }
        
        bool operator==(const long int& other) const
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
        
        bool operator!=(const long int& other) const
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
        
        Rational operator*(const long int& i) const
        {
            return Rational(n*i,d);
        }
        
        Rational operator/(const long int& i) const
        {
            return Rational(n,d*i);
        }
        
        Rational operator+(const long int& i) const
        {
            return *this+Rational(i,1);
        }
        
        Rational operator-(const long int& i) const
        {
            return *this-Rational(i,1);
        }
        
        /*************************************************************/
        // friend T& operator <<
        //template <class T>
        friend std::ostream& operator << (std::ostream& os, const Rational& r)
        {
            os<<r.n<<"/"<<r.d;
            return os;
        }
        
    };
    

    
//    Rational rat(const long int& n,const long int& d)
//    {
//        return Rational(n,d);
//    }
    
} // end namespace


namespace Eigen
{
    template<> struct NumTraits<model::Rational>
    : NumTraits<long int> // permits to get the epsilon, dummy_precision, lowest, highest functions
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
            ReadCost = 1,
            AddCost = 3,
            MulCost = 3
        };
    };
}


#endif
