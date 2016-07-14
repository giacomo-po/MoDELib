/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_H_
#define model_Polycrystal_H_

#include <utility>
#include <tuple>
#include <map>
#include <Eigen/Core>
#include <model/MPI/MPIcout.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundary.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>


namespace model
{
    
    
    
    template <int dim>
    struct Polycrystal :
    /* base */ public std::map<size_t,Grain<dim>>,
    /* base */ public std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>
    {
        
        typedef SimplicialMesh<dim> SimplicialMeshType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;


        const SimplicialMeshType& mesh;
        
        /**********************************************************************/
        Polycrystal(const SimplicialMeshType& mesh_in) :
        mesh(mesh_in)
        {
            model::cout<<"Creating Polycrystal"<<std::endl;
            for(const auto& rIter : MeshRegionObserverType::regions())
            {
                grains().emplace(rIter.second->regionID,*(rIter.second));
                //model::cout<<"mesh region "<<rIter.second->regionID<<" contains "<<rIter.second->size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            }
            
            for(const auto& rgnBnd : mesh.regionBoundaries())
            {
                grainBoundaries().emplace(std::piecewise_construct,
                                          std::forward_as_tuple(rgnBnd.first),
                                          std::forward_as_tuple(rgnBnd.second,grain(rgnBnd.first.first),grain(rgnBnd.first.second)));
                //model::cout<<"mesh region "<<rIter.second->regionID<<" contains "<<rIter.second->size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            }
            
//            for(const auto& gb : grainBoundaries())
//            {
//                std::cout<<"("<<gb.first.first<<","<<gb.first.second<<")"
//                <<"["<<gb.second.regionBoundary.regionBndID.first<<","<<gb.second.regionBoundary.regionBndID.second<<"] "
//                <<gb.second.regionBoundary.size()<<std::endl;
//            }
            
//            m.emplace(std::piecewise_construct,
//                      std::forward_as_tuple("c"),
//                      std::forward_as_tuple(10, 'c'));
            
        }
        
        /**********************************************************************/
        Grain<dim>& grain(const size_t& k)
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const Grain<dim>& grain(const size_t& k) const
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const std::map<size_t,Grain<dim>>& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t,Grain<dim>>& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundary<dim>& grainBoundary(const size_t& i,
                                                const size_t& j) const
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(1,2)) : grainBoundaries().at(std::make_pair(2,1));
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDimD& p,
                                                    const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).latticeVectorFromPosition(p);
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDimD& p) const
        {
            return latticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDimD& p,
                                                    const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).reciprocalLatticeVectorFromPosition(p);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDimD& p) const
        {
            return reciprocalLatticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        
    };
    
    
} // end namespace
#endif

