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
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/DislocationDynamics/VirtualBoundarySlipContainer.h>
#include <model/BVP/Domain.h>

#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/ParticleInteraction/ParticleSystem.h>

#include <model/Mesh/SimplicialMesh.h> // defines mode::cout
#include <model/DislocationDynamics/BVPsolver.h>


namespace model {
	
	template <typename LinkType>
	struct DislocationSharedObjects
    {
        //enum{dim=3};

        static size_t minSNorderForSolve;
		static unsigned int boundary_type;
//		static enum BoundaryType boundary_type;
		static unsigned int use_bvp;
		static Domain domain;					// temporary
		static Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> externalStress;
		static VirtualBoundarySlipContainer<LinkType> vbsc;

//        static SimplicialMesh<TypeTraits<LinkType>::dim> mesh;

//        static BVPsolver<TypeTraits<LinkType>::dim,2> bvpSolver;

//        static ParticleSystem<DislocationParticle<TypeTraits<LinkType>::dim> > particleSystem;
//         ParticleSystem<DislocationParticle<TypeTraits<LinkType>::dim> > particleSystem;

//        static SpatialCellObserver<DislocationParticle<TypeTraits<LinkType>::dim>,TypeTraits<LinkType>::dim> spatialCellObserver;

        
	};
	
	// Instantiate Static data members
	template <typename LinkType>
	unsigned int DislocationSharedObjects<LinkType>::boundary_type=0;
//    enum BoundaryType DislocationSharedObjects<LinkType>::boundary_type=noBoundary;

	
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

//    template <typename LinkType>
//	SimplicialMesh<TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::mesh(0); // this reads mesh files N/N_0.txt and T/T_0.txt
    
//    template <typename LinkType>
//	BVPsolver<TypeTraits<LinkType>::dim,2> DislocationSharedObjects<LinkType>::bvpSolver(DislocationSharedObjects<LinkType>::mesh);

    
    //0x7f8732c04dc0
    
//    template <typename LinkType>
//    ParticleSystem<DislocationParticle<TypeTraits<LinkType>::dim> > DislocationSharedObjects<LinkType>::particleSystem;
    
//	template <typename LinkType>
//    SpatialCellObserver<DislocationParticle<TypeTraits<LinkType>::dim>,TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::spatialCellObserver;

} // namespace model
#endif


