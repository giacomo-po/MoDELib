/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RationalOrthogonalApproximation_h_
#define model_RationalOrthogonalApproximation_h_

#include <iostream>
#include <iomanip>
#include <cfloat> // FLT_EPSILON
#include <assert.h> // FLT_EPSILON
#include <utility>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// #include <RoundEigen.h>


namespace model
{
    /*!A quaternion-based implementation of the rational approximations of an 
     * orthogonal matrix representing a rotation.
     *
     * Implementation from
     * [1] Rational orthogonal approximations to orthogonal matrices
     * Victor J. Milenkovic, Veljko Milenkovic
     * Computational Geometry 7 (1997) 25-35
     */
    class RationalOrthogonalApproximation
    {
        
        /**********************************************************************/
        static std::pair<Eigen::Matrix<long int,3,3>,long int> compute(const double& angle,
                                                                       const Eigen::Matrix<double,3,1>& axis,
                                                                       const double& maxDen);
        
        const std::pair<Eigen::Matrix<long int,3,3>,long int> resultPair;
        
        
    public:
        
        const Eigen::Matrix<long int,3,3>& m;
        const long int& den;

        
        /**********************************************************************/
        RationalOrthogonalApproximation(const double& angle,
                                    const Eigen::Matrix<double,3,1>& axis,
                                    const long int& maxDen) ;
    };
    
}
#endif


