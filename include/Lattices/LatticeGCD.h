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
        static long int gcd(const long int& a,const long int& b);
        
        /**********************************************************************/
        static long int gcd(const long int& a,const long int& b,const long int& c);
        
        /**********************************************************************/
        static long int gcd(const VectorDimI& v);
    };
    
} // end namespace
#endif
