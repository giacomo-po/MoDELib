/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GAUSSLEGENDRE_3_1_H_
#define model_GAUSSLEGENDRE_3_1_H_

namespace model{
    
    /**************************************************************************/
	template <>
	struct GaussLegendre<3,1>
    {
		static Eigen::Matrix<double,4,1> abcsissasAndWeights()
        {
			Eigen::Matrix<double,1,4> U;
			U<< 2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01, 1.666666666666667e-01;
			return U.transpose();
		}
	};
}
#endif

