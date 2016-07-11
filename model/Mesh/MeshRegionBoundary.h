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
#include <model/MPI/MPIcout.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/Mesh/PlanarMeshFacet.h>

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
        
        const std::pair<size_t,size_t> regionIDs;
        
//        /**********************************************************************/
//        MeshRegionBoundary(const int& rID) : regionID(rID)
//        {
//            const bool success=MeshRegionObserverType::emplace(regionID,this).second;
//            assert(success && "COULD NOT INSERT MeshRegion in MeshRegionObserver.");
//        }
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

    };
    
}	// close namespace
#endif