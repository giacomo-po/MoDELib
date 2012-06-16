/* This file is part of finite element solution of BVP attached with model "the Mechanics of Material Defects Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_point_H_
#define bvpfe_point_H_

#include <Eigen/Dense>

namespace bvpfe{


	template<short unsigned int dim = 3>
	struct Point {


		typedef Eigen::Matrix<double,dim,1> VectorDim;
		const VectorDim P;
		
		
		Point(const VectorDim& Pin) : P(Pin) {}
		
//			Vector P;     // point position

//		public :
//			//---------------------------------------
//			Point()  {P.resize(dim);}
//			//---------------------------------------
//			void setP (const Vector& P_in)
//			{
//				for (int i=0; i<dim ; i++) {
//    		
//					P[i] = double(P_in[i]);
//				}
//			}
//			//--------------------------------------
//			Vector &  getP ()	
//			{
//				return P;
//			}
//			//--------------------------------------
	};	

}  //  namespace bvpfe
#endif
