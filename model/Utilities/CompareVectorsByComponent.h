/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_COMPAREVECTORSBYCOMPONENT_H_
#define model_COMPAREVECTORSBYCOMPONENT_H_

#include <Eigen/Dense>

namespace model {

    template <short unsigned int dim>
	Eigen::Matrix<int,dim,1> floorEigen(const Eigen::Matrix<double,dim,1>& P) {
		return (Eigen::Matrix<int,dim,1>()<< (int)std::floor(P(0)), floorEigen<dim-1>(P.template segment<dim-1>(1))).finished();
	}
	
	template <>
	Eigen::Matrix<int,1,1> floorEigen<1>(const Eigen::Matrix<double,1,1>& P) {
		return (Eigen::Matrix<int,1,1>()<< (int)std::floor(P(0)) ).finished();
	}
	
	template<typename T, unsigned int N, typename CastType=T>
	struct CompareVectorsByComponent {
		bool operator() (const Eigen::Matrix<T,N,1>& lhs, const Eigen::Matrix<T,N,1>& rhs) const{
			return ( ((CastType) lhs(0))==((CastType) rhs(0)) )? CompareVectorsByComponent<T,N-1,CastType>()(lhs.template segment<N-1>(1),rhs.template segment<N-1>(1)) : (lhs(0)<rhs(0));
		}
	};
	
	template<typename T,typename CastType>
	struct CompareVectorsByComponent<T,1,CastType> {
		bool operator() (const Eigen::Matrix<T,1,1>& lhs, const Eigen::Matrix<T,1,1>& rhs) const{
			return ((CastType) lhs(0))<((CastType) rhs(0));
		}
	};
		
} // namespace model
#endif
