/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegionBoundary_cpp_
#define model_MeshRegionBoundary_cpp_


#include <MeshRegionBoundary.h>

namespace model
{

    template<int dim>
    void MeshRegionBoundary<dim>::buildSingleFace(const Simplex<dim,dim-1>* newStart,
                                                  std::shared_ptr<PlanarMeshFace<dim>> newFace,
                                                  std::set<const Simplex<dim,dim-1>*>& allSimplices)
    {
        allSimplices.erase(newStart);
        const auto regionBoundaryNeighbors(newStart->regionBoundaryNeighbors());
        for(const auto& neighbor : regionBoundaryNeighbors)
        {
            const auto neighborRgnIDs(neighbor->regionIDs());
            if(   neighborRgnIDs.find(regionBndID.first)!=neighborRgnIDs.end()  // neighbor in region
               && neighborRgnIDs.find(regionBndID.second)!=neighborRgnIDs.end()
               && allSimplices.find(neighbor)!=allSimplices.end())   // neighbor found in allSimplices
            {
                if(   (newStart->outNormal(regionBndID.first)-neighbor->outNormal(regionBndID.first)).norm()<FLT_EPSILON
                   && (newStart->outNormal(regionBndID.second)-neighbor->outNormal(regionBndID.second)).norm()<FLT_EPSILON
                   )
                {// same plane
                    allSimplices.erase(neighbor);
                    newFace->insert(neighbor);
                    buildSingleFace(neighbor,newFace,allSimplices);
                }
            }
        }
        faces().emplace(newFace->sID,newFace);
    }


    template<int dim>
    void MeshRegionBoundary<dim>::buildFaces()
    {/*!Constructs and stores the PlanarMeshFace(s) of this MeshRegionBoundary
      * This function supports non-convex regions.
      */
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"MeshRegionBoundary ("<<regionBndID.first<<","<<regionBndID.second<<") buiding faces "<<std::flush;
        
        faces().clear();
        std::set<const Simplex<dim,dim-1>*> allSimplices(simplices()); //copy simplices
        while(allSimplices.size())
        {
            std::shared_ptr<PlanarMeshFace<dim>> newFace(new PlanarMeshFace<dim>(*allSimplices.begin()));
            buildSingleFace(*allSimplices.begin(),newFace,allSimplices);
        }
        
        while(allSimplices.size())
        {
            std::shared_ptr<PlanarMeshFace<dim>> newFace(new PlanarMeshFace<dim>(*allSimplices.begin()));
            buildSingleFace(*allSimplices.begin(),newFace,allSimplices);
        }
        
        for(auto& face : faces())
        {
            face.second->finalize();
        }
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
    }



    template<int dim>
    MeshRegionBoundary<dim>::MeshRegionBoundary(const std::pair<size_t,size_t>& rbndID) :
    /* init */ regionBndID(rbndID)
    {
    }

    template<int dim>
    void MeshRegionBoundary<dim>::update()
    {
        buildFaces();
    }

    template<int dim>
    const std::set<const typename MeshRegionBoundary<dim>::SimplexType*>& MeshRegionBoundary<dim>::simplices() const
    {
        return *this;
    }

    template<int dim>
    std::set<const typename MeshRegionBoundary<dim>::SimplexType*>& MeshRegionBoundary<dim>::simplices()
    {
        return *this;
    }

    template<int dim>
    const typename MeshRegionBoundary<dim>::MeshFacesContainerType& MeshRegionBoundary<dim>::faces() const
    {
        return *this;
    }

    template<int dim>
    typename MeshRegionBoundary<dim>::MeshFacesContainerType& MeshRegionBoundary<dim>::faces()
    {
        return *this;
    }

    /**********************************************************************/
    //        std::set<const Simplex<dim,dim-2>*> unsortedBoundary() const  __attribute__ ((deprecated))
    //        {
    //            std::set<const Simplex<dim,dim-2>*> temp;
    //            for(const auto& simplex : this->simplices())
    //            {
    //                for(const auto& child : simplex->children())
    //                {
    //                    if(child->regionIDs().size()>2 || child->isBoundarySimplex())
    //                    {
    //                        temp.insert(child.get());
    //                    }
    //                }
    //            }
    //            return temp;
    //        }



//    template struct MeshRegionBoundary<2>;
    template struct MeshRegionBoundary<3>;

}	// close namespace
#endif
