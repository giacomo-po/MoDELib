/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LINEARELASTICGREENSFUNCTION_H_
#define model_LINEARELASTICGREENSFUNCTION_H_


//#include <vector>
#include <Eigen/Dense>


namespace model {
	
	struct Isotropic{};
	
	/****************************************************/
	/* LinearElasticGreensFunction: general case ********/
	/****************************************************/
	template <short unsigned int dim, typename Type>
	struct LinearElasticGreensFunction {
		
		
	};
	
	
	/****************************************************/
	/* LinearElasticGreensFunction: isotropic case ******/
	/****************************************************/
	template <short unsigned int dim>
	struct LinearElasticGreensFunction<dim,Isotropic> {
		
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		

		static MatrixDim stressKernel(const VectorDim& R, const VectorDim& t, const VectorDim& t){
			return   (C1*(1.0+1.5*coreLsquared/RaSquared)*rugauss.col(k)*(Burgers.cross(R)).transpose()
					  + 	R*(rugauss.col(k).cross(Burgers)).transpose() 
					  +  0.5* R.cross(Burgers).dot(rugauss.col(k)) * (I*(1.0+3.0*coreLsquared/RaSquared) + 3.0/RaSquared*R*R.transpose())
					  )/std::pow(RaSquared,1.5);		
		}
		
		
		
	};
	
	
}
#endif

