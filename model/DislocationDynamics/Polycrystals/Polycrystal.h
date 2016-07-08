/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_H_
#define model_Polycrystal_H_

#include <map>
#include <Eigen/Core>
#include <model/MPI/MPIcout.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundary.h>
#include <model/Mesh/SimplicialMesh.h>

namespace model
{
    
    
    
    template <int dim>
    struct Polycrystal :
    /* base */ public std::map<size_t,Grain<dim>>
//    /* base */ public std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>
    {
        
        typedef SimplicialMesh<dim> SimplicialMeshType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        const SimplicialMeshType& mesh;
        
        /**********************************************************************/
        Polycrystal(const SimplicialMeshType& mesh_in) :
        mesh(mesh_in)
        {
            model::cout<<"Creating Polycrystal"<<std::endl;
            for(auto rIter : MeshRegionObserverType::regions())
            {
                this->emplace(rIter.second->regionID,*(rIter.second));
                //model::cout<<"mesh region "<<rIter.second->regionID<<" contains "<<rIter.second->size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            }
            
        }
        
	};	
	

} // end namespace
#endif

