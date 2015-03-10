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
//#include <model/DislocationDynamics/BVP/VirtualBoundarySlipContainer.h>

#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/ParticleInteraction/ParticleSystem.h>

#include <model/Mesh/SimplicialMesh.h> // defines model::cout
#include <model/DislocationDynamics/BVP/BVPsolver.h>
//#include <model/DislocationDynamics/BVP/BoundaryDislocationNetwork.h>


namespace model {
	
	template <typename LinkType>
	struct DislocationSharedObjects
    {
        
        typedef BVPsolver<TypeTraits<LinkType>::dim,2> BvpSolverType;
        
		static bool use_boundary;
        static bool use_meshRegions;

        static unsigned int use_bvp;
		static bool use_virtualSegments;
		static Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> externalStress;

        static SimplicialMesh<TypeTraits<LinkType>::dim> mesh;
        /*!\todo Fix the Segmentation fault that takes place when mesh is destroyed.
         * When mesh is destroyed, the first Simplex<dim,dim> is destroyed, and 
         * the call to SimplexObserver<dim,dim>::removeSimplex() fails. This only
         * happens if mesh is declared static. There must be some conflict in the 
         * order of destruction of the static mesh and the static map in SimplexObserver.
         */
        
//        static BoundaryDislocationNetwork<TypeTraits<LinkType>::dim> bdn;
        
        static BvpSolverType bvpSolver;
    
	};
	
	// Static data members
	template <typename LinkType>
	bool DislocationSharedObjects<LinkType>::use_boundary=false;

    template <typename LinkType>
    bool DislocationSharedObjects<LinkType>::use_meshRegions=false;

    
	template <typename LinkType>
	unsigned int DislocationSharedObjects<LinkType>::use_bvp=0;

    template <typename LinkType>
	bool DislocationSharedObjects<LinkType>::use_virtualSegments=true;

    
	template <typename LinkType>
	Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::externalStress=Eigen::Matrix<double,TypeTraits<LinkType>::dim,TypeTraits<LinkType>::dim>::Zero();

    template <typename LinkType>
	SimplicialMesh<TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::mesh;
    
//    template <typename LinkType>
//	BoundaryDislocationNetwork<TypeTraits<LinkType>::dim> DislocationSharedObjects<LinkType>::bdn;
    
    template <typename LinkType>
	BVPsolver<TypeTraits<LinkType>::dim,2> DislocationSharedObjects<LinkType>::bvpSolver(DislocationSharedObjects<LinkType>::mesh);
    
} // namespace model
#endif


