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
#include <PlanarMeshFacet.h>

namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename _SimplexType>
    struct MeshRegionBoundary :
    /* base */ public std::set<const _SimplexType*>
//    /* base */ private std::deque<PlanarMeshFacet<_SimplexType::dim>>
    {
        typedef _SimplexType SimplexType;
//        typedef MeshRegion<SimplexType> MeshRegionType;
//        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        const std::pair<int,int> regionBndID;
        
        /**********************************************************************/
        MeshRegionBoundary(const std::pair<size_t,size_t>& rbndID) : regionBndID(rbndID)
        {
//            const bool success=MeshRegionObserverType::emplace(regionID,this).second;
//            assert(success && "COULD NOT INSERT MeshRegion in MeshRegionObserver.");
        }
//
//        /**********************************************************************/
//        MeshRegionBoundary()
//        {
//            
//            const size_t n=MeshRegionObserverType::erase(regionID);
//            assert(n==1 && "COULD NOT ERASE MeshRegion in MeshRegionObserver.");
//            
//        }
        
        /**********************************************************************/
        const std::set<const _SimplexType*>& simplices() const
        {
            return *this;
        }
        
        std::set<const _SimplexType*>& simplices()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::set<const Simplex<_SimplexType::dim,_SimplexType::dim-2>*> unsortedBoundary() const
        {
            std::set<const Simplex<_SimplexType::dim,_SimplexType::dim-2>*> temp;
            for(const auto& simplex : simplices())
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
