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

#include <model/DislocationDynamics/NearestNeighbor/DislocationParticle.h>
#include <model/ParticleInteraction/ParticleSystem.h>

#include <model/Mesh/SimplicialMesh.h> // defines model::cout
#include <model/DislocationDynamics/BVP/BVPsolver.h>
#include <model/DislocationDynamics/ExternalStressFieldController.h>

namespace model
{
    
    template <int dim>
    struct DislocationSharedObjects
    {
        
        typedef BVPsolver<dim,2> BvpSolverType;
        typedef ExternalStressFieldController<dim> ExternalStressFieldControllerType;
        
        static bool use_boundary;
        static bool use_externalStress;
        static bool use_externaldislocationstressfield;
        //        static bool use_meshRegions;
        
        static unsigned int use_bvp;
        static bool use_virtualSegments;
        //		static Eigen::Matrix<double,dim,dim> externalStress;
        
        static SimplicialMesh<dim> mesh;
        /*!\todo Fix the Segmentation fault that takes place when mesh is destroyed.
         * When mesh is destroyed, the first Simplex<dim,dim> is destroyed, and
         * the call to SimplexObserver<dim,dim>::removeSimplex() fails. This only
         * happens if mesh is declared static. There must be some conflict in the
         * order of destruction of the static mesh and the static map in SimplexObserver.
         */
        
        
        
        //        static BoundaryDislocationNetwork<TypeTraits<LinkType>::dim> bdn;
        
        static BvpSolverType bvpSolver;
        static ExternalStressFieldControllerType extStressController;
        static std::vector<StressStraight<dim>> ssdeq;
        
    };
    
    // Static data members
    template <int dim>
    bool DislocationSharedObjects<dim>::use_boundary=false;
    
    template <int dim>
    bool DislocationSharedObjects<dim>::use_externalStress=false;
    
    template <int dim>
    bool DislocationSharedObjects<dim>::use_externaldislocationstressfield=false;
    
    //    template <int dim>
    //    bool DislocationSharedObjects<dim>::use_meshRegions=false;
    
    
    template <int dim>
    unsigned int DislocationSharedObjects<dim>::use_bvp=0;
    
    
    
    template <int dim>
    bool DislocationSharedObjects<dim>::use_virtualSegments=true;
    
    
    //    template <int dim>
    //	Eigen::Matrix<double,dim,dim> DislocationSharedObjects<dim>::externalStress=Eigen::Matrix<double,dim,dim>::Zero();
    
    template <int dim>
    SimplicialMesh<dim> DislocationSharedObjects<dim>::mesh;
    
    template <int dim>
    BVPsolver<dim,2> DislocationSharedObjects<dim>::bvpSolver(DislocationSharedObjects<dim>::mesh);
    
    template <int dim>
    ExternalStressFieldController<dim> DislocationSharedObjects<dim>::extStressController;
    
    template <int dim>
    std::vector<StressStraight<dim>> DislocationSharedObjects<dim>::ssdeq;    
} // namespace model
#endif


