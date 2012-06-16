/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLCATIONSHAREDOBJECTS_H_
#define model_DISLCATIONSHAREDOBJECTS_H_

#include <Eigen/Dense>
#include <model/BVP/Domain.h>

namespace model {
	
	template <short unsigned int dim, typename MaterialType>
	struct DislocationSharedObjects {

		static unsigned int boundary_type;
		static unsigned int use_bvp;
		static bvpfe::Domain domain;					// temporary
		static MaterialType material;
		static Eigen::Matrix<double,dim,dim> externalStress;
		
	};
	
	// Instantiate Static data members
	template <short unsigned int dim,typename MaterialType>
	unsigned int DislocationSharedObjects<dim,MaterialType>::boundary_type=0;
	
	template <short unsigned int dim,typename MaterialType>
	unsigned int DislocationSharedObjects<dim,MaterialType>::use_bvp=0;
	
	template <short unsigned int dim,typename MaterialType>
	bvpfe::Domain DislocationSharedObjects<dim,MaterialType>::domain;
		
	template <short unsigned int dim,typename MaterialType>
	MaterialType DislocationSharedObjects<dim,MaterialType>::material;

	template <short unsigned int dim,typename MaterialType>
	Eigen::Matrix<double,dim,dim> DislocationSharedObjects<dim,MaterialType>::externalStress=Eigen::Matrix<double,dim,dim>::Zero();
	
} // namespace model
#endif


//#include <model/Dislocations/ExternalStresses/LoadController.h>
//		static model::LoadController<dim> loadController;

//	template <short unsigned int dim,typename MaterialType>
//	LoadController<dim> DislocationSharedObjects<dim,MaterialType>::loadController;
