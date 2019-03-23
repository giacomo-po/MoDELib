/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeGCD_h_
#define model_LatticeGCD_h_

#include <Eigen/Dense>

namespace model
{
    template <int dim>
    struct LatticeGCD
    {

        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        
        /**********************************************************************/
        static long int gcd(const long int& a,const long int& b)
        {
            return abs(b)>0? gcd(abs(b), abs(a) % abs(b)) : abs(a);
        }
        
        /**********************************************************************/
        static long int gcd(const long int& a,const long int& b,const long int& c)
        {
            return gcd(a,gcd(b,c));
        }
        
        /**********************************************************************/
        static long int gcd(const VectorDimI& v)
        {
            //            assert(v.squaredNorm()>0 && "Cannot extract LatticeDirection from zero-vector.");
            //            int g=gcd(abs(v(0)),abs(v(1)),abs(v(2)));
            //            return v/g;
            
            return gcd(v(0),v(1),v(2));
        }
        
        
        const long int gCD;
        
        LatticeGCD(const VectorDimI& v) :
        /* base init */ gCD(gcd(v))
        {}
        
    };
    
} // end namespace
#endif
