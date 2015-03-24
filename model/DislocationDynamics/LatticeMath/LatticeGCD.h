/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeGCD_h_
#define model_LatticeGCD_h_

//#include <model/DislocationDynamics/LatticeMath/LatticeVector.h>
//#include <model/DislocationDynamics/LatticeMath/LatticePoint.h>

namespace model
{
    template <int dim>
    struct LatticeGCD
    {
        
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        /**********************************************************************/
        static int gcd(const size_t& a,const size_t& b)
        {
            return b>0? gcd(b, a % b) : a;
        }
        
        /**********************************************************************/
        static int gcd(const size_t& a,const size_t& b,const size_t& c)
        {
            return gcd(a,gcd(b,c));
        }
        
        static VectorDimI gcd(const VectorDimI& v)
        {
            int g=gcd(abs(v(0)),abs(v(1)),abs(v(2)));
            return v/g;
        }
        
    };
    
} // end namespace
#endif
