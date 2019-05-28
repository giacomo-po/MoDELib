/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegion_H_
#define model_MeshRegion_H_

#include <set>
#include <deque>
#include <memory>
#include <assert.h>
#include <MPIcout.h>
#include <MeshRegionObserver.h>
#include <PlanarMeshFace.h>
#include <Simplex.h>
#include <TerminalColors.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename _SimplexType>
    struct MeshRegion : public std::set<const _SimplexType*>
    /*               */,public std::map<int,std::shared_ptr<PlanarMeshFace<_SimplexType::dim>>> // MeshRegionBoundary container
    {
        typedef _SimplexType SimplexType;
        static constexpr int dim=SimplexType::dim;
        typedef MeshRegion<SimplexType> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef std::map<int,std::shared_ptr<PlanarMeshFace<dim>>> MeshFacesContainerType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
    private:
        
        /**********************************************************************/
        void buildSingleFace(const Simplex<dim,dim-1>* newStart,
                        std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace,
                        std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*>& allSimplices)
        {
            allSimplices.erase(newStart->xID);
            const auto boundaryNeighbors(newStart->boundaryNeighbors());
            for(const auto& neighbor : boundaryNeighbors)
            {
                const auto neighborRgnIDs(neighbor->regionIDs());
                if(   neighborRgnIDs.find(regionID)!=neighborRgnIDs.end()  // neighbor in region
                   && allSimplices.find(neighbor->xID)!=allSimplices.end())   // neighbor found in allSimplices
                {
                    if((newStart->outNormal()-neighbor->outNormal()).norm()<FLT_EPSILON)
                    {// same plane
                        allSimplices.erase(neighbor->xID);
                        newFace->insert(neighbor);
                        buildSingleFace(neighbor,newFace,allSimplices);
                    }
                }
            }
            faces().emplace(newFace->sID,newFace);
        }
        
        /**********************************************************************/
        void buildFaces()
        {/*!Constructs and stores the external PlanarMeshFace(s) of this MeshRegion
          * This function supports non-convex regions.
          */
            faces().clear();
            std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*> allSimplices;
            for(const auto& simplex : simplices())
            {// collect each possible boundary/region_boundary simplices
                for(const auto& child : simplex->children())
                {
                    if(child->isBoundarySimplex())
                    {
                        allSimplices.emplace(child->xID,child.get());
                    }
                }
            }

            while(allSimplices.size())
            {
                std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace(new PlanarMeshFace<_SimplexType::dim>(allSimplices.begin()->second));
                buildSingleFace(allSimplices.begin()->second,newFace,allSimplices);
            }

            for(auto& face : faces())
            {
                face.second->finalize();
            }
        }
        
        MeshRegionObserverType& regionObserver;
        
    public:
        const int regionID;
        
        /**********************************************************************/
        MeshRegion(MeshRegionObserverType& ro,
                   const int& rID) :
        /* init */ regionObserver(ro),
        /* init */ regionID(rID)
        {
            const bool success=regionObserver.emplace(regionID,this).second;
            assert(success && "COULD NOT INSERT MeshRegion in MeshRegionObserver.");
        }
        
        /**********************************************************************/
        ~MeshRegion()
        {
            const size_t n=regionObserver.erase(regionID);
            assert(n==1 && "COULD NOT ERASE MeshRegion in MeshRegionObserver.");
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
        
    };
    
}	// close namespace
#endif
