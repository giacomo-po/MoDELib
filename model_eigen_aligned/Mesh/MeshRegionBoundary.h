/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegionBoundary_H_
#define model_MeshRegionBoundary_H_

#include <set>
#include <deque>
#include <assert.h>
#include <MPIcout.h>
#include <MeshRegionObserver.h>
#include <PlanarMeshFace.h>

namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename _SimplexType>
    struct MeshRegionBoundary : public std::set<const _SimplexType*>
    /*                       */,public std::map<int,std::shared_ptr<PlanarMeshFace<_SimplexType::dim>>>
    {
        typedef _SimplexType SimplexType;
        static constexpr int dim=SimplexType::dim;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::map<int,std::shared_ptr<PlanarMeshFaceType>> MeshFacesContainerType;
        
    private:
        
        
        /**********************************************************************/
        void buildSingleFace(const Simplex<dim,dim-1>* newStart,
                             std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace,
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
        
        
        /**********************************************************************/
        void buildFaces()
        {/*!Constructs and stores the PlanarMeshFace(s) of this MeshRegionBoundary
          * This function supports non-convex regions.
          */
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"MeshRegionBoundary ("<<regionBndID.first<<","<<regionBndID.second<<") buiding faces "<<std::flush;

            faces().clear();
            std::set<const Simplex<dim,dim-1>*> allSimplices(simplices()); //copy simplices
            while(allSimplices.size())
            {
                std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace(new PlanarMeshFace<_SimplexType::dim>(*allSimplices.begin()));
                buildSingleFace(*allSimplices.begin(),newFace,allSimplices);
            }
            
            while(allSimplices.size())
            {
                std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace(new PlanarMeshFace<_SimplexType::dim>(*allSimplices.begin()));
                buildSingleFace(*allSimplices.begin(),newFace,allSimplices);
            }
            
            for(auto& face : faces())
            {
                face.second->finalize();
            }
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

        }
        
        
    public:
        
        const std::pair<int,int> regionBndID;
        
        /**********************************************************************/
        MeshRegionBoundary(const std::pair<size_t,size_t>& rbndID) :
        /* init */ regionBndID(rbndID)
        {
        }
        
        /**********************************************************************/
        void update()
        {
            buildFaces();
        }
        
        /**********************************************************************/
        const std::set<const _SimplexType*>& simplices() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::set<const _SimplexType*>& simplices()
        {
            return *this;
        }
        
        /**********************************************************************/
        const MeshFacesContainerType& faces() const
        {
            return *this;
        }
        
        /**********************************************************************/
        MeshFacesContainerType& faces()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::set<const Simplex<_SimplexType::dim,_SimplexType::dim-2>*> unsortedBoundary() const  __attribute__ ((deprecated))
        {
            std::set<const Simplex<_SimplexType::dim,_SimplexType::dim-2>*> temp;
            for(const auto& simplex : this->simplices())
            {
                for(const auto& child : simplex->children())
                {
                    if(child->regionIDs().size()>2 || child->isBoundarySimplex())
                    {
                        temp.insert(child.get());
                    }
                }
            }
            return temp;
        }
        
    };
    
}	// close namespace
#endif
