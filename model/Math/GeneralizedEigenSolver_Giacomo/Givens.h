/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_GIVENS_H_
#define  model_GIVENS_H_

#include <Eigen/Dense>
#include <math.h>  // for std::abs
#include <utility> // for std::pair

namespace model {
	
	struct Givens{

		template <typename scalar>
		static std::pair<scalar,scalar> get_cs(const scalar&a, const scalar&b) {
			/*! The pair (c,s) where c=cos(theta) and s=sin(theta)
			 */
			std::pair<scalar,scalar> temp(std::make_pair(1.0,0.0)); // initialize to the case b=0
			if (b!=0.0){
				if (std::fabs(b)>std::fabs(a)){
					scalar tau=-a/b;
					temp.second=1.0*std::pow(1.0+std::pow(tau,2),-0.5);
					temp.first=temp.second*tau;
				}
				else {
					scalar tau=-b/a;
					temp.first=1.0*std::pow(1.0+std::pow(tau,2),-0.5);
					temp.second=temp.first*tau;
				}
			}
			return temp;
		}

		template <typename scalar>
		static Eigen::Matrix<scalar,2,2> Q(const scalar& a,const scalar& b){
			const std::pair<scalar,scalar> cs (get_cs<scalar>(a,b));
			return (Eigen::Matrix<scalar,2,2>() << cs.first , -cs.second,
			/*                                  */ cs.second,  cs.first).finished();
		 }
	};

} // close namespace
#endif


