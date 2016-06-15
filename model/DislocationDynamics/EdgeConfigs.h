/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EDGECONFIGS_H_
#define model_EDGECONFIGS_H_

#include <Eigen/Core>
#include <model/Math/CompileTimeMath/Pow.h>

namespace model {
	
	/**********************************************************************/
	/**********************************************************************/
	template <unsigned int N>
	struct EdgeConfigs
    {
		
		static Eigen::Matrix<int,Pow<2,N-1>::value,N> get_Ci()
        {
			Eigen::Matrix<int,Pow<2,N-1>::value,N> temp;
			temp.template block<Pow<2,N-2>::value,N-1>(0,1)=EdgeConfigs<N-1>::get_Ci();
			temp.template block<Pow<2,N-2>::value,N-1>(Pow<2,N-2>::value,1)=-EdgeConfigs<N-1>::get_Ci();
			temp.col(0).setOnes();
			return temp;
		}
		
	};	
	
	template <>
	struct EdgeConfigs<1>
    {

		enum {N=1};
		
		static Eigen::Matrix<int,Pow<2,N-1>::value,N> get_Ci()
        {
			return Eigen::Matrix<int,Pow<2,N-1>::value,N>::Ones();
		}
		
	};


	/**********************************************************************/
	/**********************************************************************/
	struct EdgeDynamicConfigs
    {
        
        static const size_t maxEdge;
		
		static Eigen::MatrixXi getCi(const unsigned int& n)
        {
			Eigen::MatrixXi Ci;
			switch (n)
            {
				case 1:
					Ci=EdgeConfigs<1>::get_Ci();
					break;
				case 2:
					Ci=EdgeConfigs<2>::get_Ci();
					break;
				case 3:
					Ci=EdgeConfigs<3>::get_Ci();
					break;
				case 4:
					Ci=EdgeConfigs<4>::get_Ci();
					break;
				case 5:
					Ci=EdgeConfigs<5>::get_Ci();
					break;
				case 6:
					Ci=EdgeConfigs<6>::get_Ci();
					break;
				case 7:
					Ci=EdgeConfigs<7>::get_Ci();
					break;
				case 8:
					Ci=EdgeConfigs<8>::get_Ci();
					break;
//				case 9:
//					Ci=EdgeConfigs<9>::get_Ci();
//					break;
//				case 10:
//					Ci=EdgeConfigs<10>::get_Ci();
//					break;
//				case 11:
//					Ci=EdgeConfigs<11>::get_Ci();
//					break;
//				case 12:
//					Ci=EdgeConfigs<12>::get_Ci();
//					break;
//                case 13:
//					Ci=EdgeConfigs<13>::get_Ci();
//					break;
//                case 14:
//					Ci=EdgeConfigs<14>::get_Ci();
//					break;
//                case 15:
//					Ci=EdgeConfigs<15>::get_Ci();
//					break;
					
					
				default:
					assert(0 && "EDGE CONFIGURATION NOT IMPLEMENTED.");				
					break;
			}
			return Ci;
		}
		
	};
    
    const size_t EdgeDynamicConfigs::maxEdge=8;
	
	/*********************************************************************/
	/*********************************************************************/
} // end namespace
#endif

