/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegionObserver_H_
#define model_MeshRegionObserver_H_

#include <memory> // shared_ptr
#include <map>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename RegionType>
    class MeshRegionObserver
    {
        typedef typename RegionType::SimplexType SimplexType;
        typedef std::map<int,RegionType* const> RegionMapType;
        typedef std::shared_ptr<RegionType> SharedPtrType;
        
        static RegionMapType regionMap;
        
    public:
        
        /**********************************************************************/
        static const RegionMapType& regions()
        {
            return regionMap;
        }
        
        /**********************************************************************/
        static SharedPtrType getRegion(const int& k)
        {
            typename RegionMapType::const_iterator iter(regionMap.find(k));
            return (iter!=regionMap.end())? (*(iter->second->simplices().begin()))->region : SharedPtrType(new RegionType(k));
        }
        
        /**********************************************************************/
        static size_t erase(const int& k)
        {
            return regionMap.erase(k);
        }
        
        /**********************************************************************/
        static std::pair<typename RegionMapType::iterator,bool> emplace(const int& k, RegionType* const reg)
        {
            return regionMap.emplace(k,reg);
        }
        
    };
    
    
    template<typename RegionType>
    std::map<int,RegionType* const> MeshRegionObserver<RegionType>::regionMap;
}	// close namespace
#endif