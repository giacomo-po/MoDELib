/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RationalOrthogonalApproximation_cpp_
#define model_RationalOrthogonalApproximation_cpp_

#include <RationalOrthogonalApproximation.h>

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

    /**********************************************************************/
    std::pair<Eigen::Matrix<long int, 3, 3>, long int> RationalOrthogonalApproximation::compute(const double &angle,
                                                                                                const Eigen::Matrix<double, 3, 1> &axis,
                                                                                                const double &maxDen)
    {
        const Eigen::AngleAxisd aa(angle, axis.normalized());
        const Eigen::Quaterniond q(aa);

        const double root(sqrt(fabs(maxDen)));

        Eigen::Quaterniond q1(std::round(q.w() * root),
                              std::round(q.x() * root),
                              std::round(q.y() * root),
                              std::round(q.z() * root));

        const double n2(q1.squaredNorm());
        q1.normalize();

        Eigen::Matrix<double, 3, 3> R(q1.toRotationMatrix());
        Eigen::Matrix<long int, 3, 3> Q;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (std::fabs(R(i, j) * n2 - std::round(R(i, j) * n2)) > FLT_EPSILON)
                {
                    std::cout << i << " " << j << std::endl;
                    std::cout << R(i, j) * n2 << std::endl;
                    std::cout << std::round(R(i, j) * n2) << std::endl;
                    std::cout << R * n2 << std::endl;
                    assert(0 && "RationalOrthogonalApproximation failed.");
                }
                else
                {
                    Q(i, j) = std::round(R(i, j) * n2);
                }
            }
        }

        return std::make_pair(Q, n2);
    }


  /**********************************************************************/
    RationalOrthogonalApproximation::RationalOrthogonalApproximation(const double& angle,
                                                                    const Eigen::Matrix<double,3,1>& axis,
                                                                    const long int& maxDen) :
                                                                    /* init */ resultPair(compute(angle,axis,maxDen)),
                                                                    /* init */ m(resultPair.first),
                                                                    /* init */ den(resultPair.second)
    {
    
    }
}
#endif
