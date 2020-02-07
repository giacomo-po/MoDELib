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
            const long int absA(abs(a));
            const long int absB(abs(b));
            return absB>0? gcd(absB, absA % absB) : (absA>0? absA : 1);
        }
        
        /**********************************************************************/
        static long int gcd(const long int& a,const long int& b,const long int& c)
        {
            return gcd(a,gcd(b,c));
        }
        
        /**********************************************************************/
        static long int gcd(const VectorDimI& v)
        {
            return gcd(v(0),v(1),v(2));
        }
    };
    
} // end namespace
#endif
