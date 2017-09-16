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
    struct MeshRegionObserver : public std::map<int,RegionType* const>
    {
        typedef typename RegionType::SimplexType SimplexType;
        typedef std::map<int,RegionType* const> RegionMapType;
        typedef std::shared_ptr<RegionType> SharedPtrType;
        
////        static RegionMapType regionMap;
//        
//    public:
//        
        /**********************************************************************/
        const RegionMapType& regions() const
        {
            return *this;
        }
        
        /**********************************************************************/
        RegionMapType& regions()
        {
            return *this;
        }

        /**********************************************************************/
        RegionType* region(const int& k)
        {
            typename RegionMapType::const_iterator iter(this->find(k));
            assert(iter!=this->end());
            return iter.second;
//            return (iter!=this->end())? (*(iter->second->simplices().begin()))->region : SharedPtrType(new RegionType(*this,k));
        }
        
        /**********************************************************************/
        const RegionType* region(const int& k) const
        {
            typename RegionMapType::const_iterator iter(this->find(k));
            assert(iter!=this->end());
            return iter->second;
            //            return (iter!=this->end())? (*(iter->second->simplices().begin()))->region : SharedPtrType(new RegionType(*this,k));
        }
        
        /**********************************************************************/
        SharedPtrType getSharedRegion(const int& k)
        {
            typename RegionMapType::const_iterator iter(this->find(k));
            return (iter!=this->end())? (*(iter->second->simplices().begin()))->region : SharedPtrType(new RegionType(*this,k));
        }
        
//        /**********************************************************************/
//        static size_t erase(const int& k)
//        {
//            return regionMap.erase(k);
//        }
        
//        /**********************************************************************/
//        static std::pair<typename RegionMapType::iterator,bool> emplace(const int& k, RegionType* const reg)
//        {
//            return regionMap.emplace(k,reg);
//        }
        
    };
    
    
//    template<typename RegionType>
//    std::map<int,RegionType* const> MeshRegionObserver<RegionType>::regionMap;
}	// close namespace
#endif
