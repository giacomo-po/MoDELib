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
#include <assert.h>
#include <model/MPI/MPIcout.h>
#include <model/Mesh/MeshRegionObserver.h>
//#include <model/Mesh/PlanarMeshFacet.h>

namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename _SimplexType>
    struct MeshRegion : public std::set<const _SimplexType*>
//    /* base */ private std::deque<PlanarMeshFacet<_SimplexType::dim>>
    {
        typedef _SimplexType SimplexType;
        typedef MeshRegion<SimplexType> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        MeshRegionObserverType& regionObserver;
        const int regionID;
        
        /**********************************************************************/
        MeshRegion(MeshRegionObserverType& ro,
                   const int& rID) :
        /* init */ regionObserver(ro),
        /* init */ regionID(rID)
        {
//            model::cout<<"Creating Mesh Region "<<rID<<std::flush;
            const bool success=regionObserver.emplace(regionID,this).second;
            assert(success && "COULD NOT INSERT MeshRegion in MeshRegionObserver.");
//            model::cout<<" done"<<std::endl;
        }
        
        /**********************************************************************/
        ~MeshRegion()
        {
            std::cout<<"Destroying Region..."<<std::flush;
            const size_t n=regionObserver.erase(regionID);
            assert(n==1 && "COULD NOT ERASE MeshRegion in MeshRegionObserver.");
            std::cout<<"done"<<std::endl;
        }
        
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
