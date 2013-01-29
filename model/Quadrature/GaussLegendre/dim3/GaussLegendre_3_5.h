 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_3_5_H_
#define model_GAUSSLEGENDRE_3_5_H_

namespace model{

	//////////////////////////////////////////////////////////////
	template <>
	struct GaussLegendre<3,5> {
		static Eigen::Matrix<double,4,5> abcsissasAndWeights(){
			Eigen::Matrix<double,4,5> U;
			U(0,0)=  2.500000000000000e-01;		U(3,0)= -1.333333333333333e-01;
			U(1,0)=  2.500000000000000e-01;
			U(2,0)=  2.500000000000000e-01;
			
			U(0,1)=	 1.666666666666667e-01;		U(3,1)= 7.500000000000000e-02;
			U(1,1)=  1.666666666666667e-01;
			U(2,1)=  1.666666666666667e-01;
			
			U(0,2)=  1.666666666666667e-01;		U(3,2)= 7.500000000000000e-02;
			U(1,2)=  5.000000000000000e-01;
			U(2,2)=  1.666666666666667e-01;
			
			U(0,3)=  1.666666666666667e-01;		U(3,3)= 7.500000000000000e-02;
			U(1,3)=  1.666666666666667e-01;
			U(2,3)=  5.000000000000000e-01;
			
			U(0,4)=  5.000000000000000e-01;		U(3,4)= 7.500000000000000e-02;
			U(1,4)=  1.666666666666667e-01;
			U(2,4)=  1.666666666666667e-01;
			return U;
		}
	};
/*************************************************/
} 
#endif 

