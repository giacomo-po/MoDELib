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
#include <model/Dislocations/DislocationNetworkTraits.h>
#include <model/BVP/VirtualBoundarySlipContainer.h>
#include <model/BVP/Domain.h>

namespace model {
	
//	template <short unsigned int dim, typename _MaterialType, typename LinkType>
	template <typename LinkType>
	struct DislocationSharedObjects {


		static unsigned int boundary_type;
		static unsigned int use_bvp;
		static bvpfe::Domain domain;					// temporary
		static typename TypeTraits<LinkType>::MaterialType material;
		static Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> externalStress;
		static VirtualBoundarySlipContainer<LinkType> vbsc;
	};
	
	// Instantiate Static data members
	template <typename LinkType>
	unsigned int DislocationSharedObjects<LinkType>::boundary_type=0;
	
	template <typename LinkType>
	unsigned int DislocationSharedObjects<LinkType>::use_bvp=0;
	
	template <typename LinkType>
	bvpfe::Domain DislocationSharedObjects<LinkType>::domain;
		
	template <typename LinkType>
	typename TypeTraits<LinkType>::MaterialType DislocationSharedObjects<LinkType>::material;

	template <typename LinkType>
	Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::externalStress=Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim>::Zero();

	template <typename LinkType>
	VirtualBoundarySlipContainer<LinkType> DislocationSharedObjects<LinkType>::vbsc;


	
} // namespace model
#endif


