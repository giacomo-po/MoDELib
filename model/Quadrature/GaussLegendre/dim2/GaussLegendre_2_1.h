 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_2_1_H_
#define model_GAUSSLEGENDRE_2_1_H_

namespace model{

	//////////////////////////////////////////////////////////////
	template <>
	struct GaussLegendre<2,1>{
		static Eigen::Matrix<double,3,1> abcsissasAndWeights(){
			Eigen::Matrix<double,1,3> aw;
			aw<< 3.333333333333333e-01, 3.333333333333333e-01, 5.000000000000000e-01;
			return aw.transpose();
		}
	};
/*************************************************/
} 
#endif 

