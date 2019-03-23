 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_3_4_H_
#define model_GAUSSLEGENDRE_3_4_H_

namespace model{

	//////////////////////////////////////////////////////////////
	template <>
	struct GaussLegendre<3,4>{
		//a=(5-sqrt(5))/20=1.381966011250105e-01
		//b=(5+3*sqrt(5))/20=5.854101966249685e-01
		static Eigen::Matrix<double,4,4> abcsissasAndWeights(){
			Eigen::Matrix<double,4,4> U;
			U(0,0)=  1.381966011250105e-01;		U(3,0)= 4.166666666666667e-02;
			U(1,0)=  1.381966011250105e-01;
			U(2,0)=  1.381966011250105e-01;
			
			U(0,1)=	 1.381966011250105e-01;		U(3,1)= 4.166666666666667e-02;
			U(1,1)=  5.854101966249685e-01;
			U(2,1)=  1.381966011250105e-01;
			
			U(0,2)=  1.381966011250105e-01;		U(3,2)= 4.166666666666667e-02;
			U(1,2)=  1.381966011250105e-01;
			U(2,2)=  5.854101966249685e-01;
			
			U(0,3)=  5.854101966249685e-01;		U(3,3)= 4.166666666666667e-02;
			U(1,3)=  1.381966011250105e-01;
			U(2,3)=  1.381966011250105e-01;
			return U;
		}
	};
/*************************************************/
} 
#endif 

