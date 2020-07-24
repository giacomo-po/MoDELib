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
    class LoopObserver : public std::map<size_t,LoopType* const>
    {

    public:
        
        typedef std::map<size_t,LoopType* const> LoopContainerType;
        
        typedef std::pair<bool, const LoopType* const> IsConstLoopType;
        
//    private:
//        
//        static LoopContainerType loopMap;
        
    public:
        
        IsConstLoopType loop(const size_t& i)
        {
            typename LoopContainerType::const_iterator loopIter(this->find(i));
            return (loopIter==this->end())?  std::make_pair(false,static_cast<const LoopType* const>(nullptr)) :
            /*                              */ std::make_pair(true,loopIter->second);
            
        }

        
        /**********************************************************************/
        LoopContainerType& loops()
        {
            return *this;
        }
        
        /**********************************************************************/
        const LoopContainerType& loops() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addLoop(LoopType* const pL)
        {
            const bool success=this->insert(std::make_pair(pL->sID,pL)).second;
            assert(success && "Could not insert in loop");
        }
        
        /**********************************************************************/
        void removeLoop(const LoopType* const pL)
        {
            const size_t erased=this->erase(pL->sID);
            assert(erased==1 && "Could not erase loop");
        }
        
    };
    
//    template<typename LoopType>
//    std::map<size_t,const LoopType* const> LoopObserver<LoopType>::loopMap;
    
}
#endif
