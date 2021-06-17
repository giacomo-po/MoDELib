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
    struct MeshRegionObserver : public std::map<size_t,RegionType* const>
    {
        typedef typename RegionType::SimplexType SimplexType;
        typedef std::map<size_t,RegionType* const> RegionMapType;
        typedef std::shared_ptr<RegionType> SharedPtrType;
        
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
        RegionType* region(const size_t& k)
        {
            typename RegionMapType::const_iterator iter(this->find(k));
            assert(iter!=this->end());
            return iter->second;
        }
        
        /**********************************************************************/
        const RegionType* region(const size_t& k) const
        {
            typename RegionMapType::const_iterator iter(this->find(k));
            assert(iter!=this->end());
            return iter->second;
        }
        
        /**********************************************************************/
        SharedPtrType getSharedRegion(const size_t& k)
        {
            typename RegionMapType::const_iterator iter(this->find(k));
            return (iter!=this->end())? (*(iter->second->simplices().begin()))->region : SharedPtrType(new RegionType(*this,k));
        }
    };

}	// close namespace
#endif
