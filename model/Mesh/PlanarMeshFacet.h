/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarMeshFacet_H_
#define model_PlanarMeshFacet_H_

#include <set>
#include <assert.h>
#include <Simplex.h>
#include <MPIcout.h>

namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct PlanarMeshFacet:
    /* base */ std::set<const Simplex<dim,dim-1>*>
    {
//        typedef _SimplexType SimplexType;
//        typedef MeshRegion<SimplexType> MeshRegionType;
//        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
//        
//        const int regionID;
//        
//        /**********************************************************************/
//        MeshRegion(const int& rID) : regionID(rID)
//        {
//            const bool success=MeshRegionObserverType::emplace(regionID,this).second;
//            assert(success && "COULD NOT INSERT MeshRegion in MeshRegionObserver.");
//        }
//        
//        /**********************************************************************/
//        ~MeshRegion()
//        {
//            
//            const size_t n=MeshRegionObserverType::erase(regionID);
//            assert(n==1 && "COULD NOT ERASE MeshRegion in MeshRegionObserver.");
//            
//        }

    };
    
}	// close namespace
#endif
