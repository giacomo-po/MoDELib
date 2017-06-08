/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopObserver_H_
#define model_LoopObserver_H_

#include <map>
#include <utility>

namespace model
{
    template<typename LoopType>
    class LoopObserver
    {

    public:
        
        typedef std::map<size_t,const LoopType* const> LoopContainerType;
        
        typedef std::pair<bool, const LoopType* const> IsConstLoopType;
        
    private:
        
        static LoopContainerType loopMap;
        
    public:
        
        static IsConstLoopType loop(const size_t& i)
        {
            typename LoopContainerType::const_iterator loopIter(loopMap.find(i));
            return (loopIter==loopMap.end())?  std::make_pair(false,static_cast<const LoopType* const>(nullptr)) :
            /*                              */ std::make_pair(true,loopIter->second);
            
        }

        
        /**********************************************************************/
        static LoopContainerType& loops()
        {
            return loopMap;
        }
        
        /**********************************************************************/
        static void addLoop(const LoopType* const pL)
        {
            const bool success=loopMap.insert(std::make_pair(pL->sID,pL)).second;
            assert(success && "Could not insert in loopMap");
        }
        
        /**********************************************************************/
        static void removeLoop(const LoopType* const pL)
        {
            const size_t erased=loopMap.erase(pL->sID);
            assert(erased==1 && "Could not erase from loopMap");
        }
        
    };
    
    template<typename LoopType>
    std::map<size_t,const LoopType* const> LoopObserver<LoopType>::loopMap;
    
}
#endif
