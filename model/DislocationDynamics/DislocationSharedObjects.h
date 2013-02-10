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
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/BVP/VirtualBoundarySlipContainer.h>
#include <model/BVP/Domain.h>

namespace model {
	
	template <typename LinkType>
	struct DislocationSharedObjects {

        static size_t minSNorderForSolve;
		static unsigned int boundary_type;
		static unsigned int use_bvp;
		static Domain domain;					// temporary
		static Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> externalStress;
		static VirtualBoundarySlipContainer<LinkType> vbsc;
	};
	
	// Instantiate Static data members
	template <typename LinkType>
	unsigned int DislocationSharedObjects<LinkType>::boundary_type=0;
	
	template <typename LinkType>
	unsigned int DislocationSharedObjects<LinkType>::use_bvp=0;
	
	template <typename LinkType>
	Domain DislocationSharedObjects<LinkType>::domain;
    
    template <typename LinkType>
	size_t DislocationSharedObjects<LinkType>::minSNorderForSolve=0;
		
	template <typename LinkType>
	Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::externalStress=Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim>::Zero();

	template <typename LinkType>
	VirtualBoundarySlipContainer<LinkType> DislocationSharedObjects<LinkType>::vbsc;
	
} // namespace model
#endif


