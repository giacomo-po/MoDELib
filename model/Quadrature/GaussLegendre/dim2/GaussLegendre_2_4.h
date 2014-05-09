/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GAUSSLEGENDRE_2_4_H_
#define model_GAUSSLEGENDRE_2_4_H_

namespace model
{
    /**************************************************************************/
	template <>
	struct GaussLegendre<2,4>
    {
		static Eigen::Matrix<double,3,4> abcsissasAndWeights()
        {
			Eigen::Matrix<double,4,3> aw;
			aw<< 3.333333333333333e-01, 3.333333333333333e-01, -2.812500000000000e-01, // yes, negative sign in weight
            /**/ 2.000000000000000e-01, 2.000000000000000e-01,  2.604166666666667e-01,
            /**/ 6.000000000000000e-01, 2.000000000000000e-01,  2.604166666666667e-01,
            /**/ 2.000000000000000e-01, 6.000000000000000e-01,  2.604166666666667e-01;
			return aw.transpose();
		}
	};
}
#endif

