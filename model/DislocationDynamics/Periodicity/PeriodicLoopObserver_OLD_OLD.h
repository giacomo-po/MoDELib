/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLoopObserver_H_
#define model_PeriodicLoopObserver_H_

#include <map>
#include <utility>

namespace model
{
    template<typename PeriodicLoopType>
    struct PeriodicLoopObserver : public std::map<size_t,PeriodicLoopType* const>
    {

        
        typedef std::map<size_t,PeriodicLoopType* const> LoopContainerType;
        
        
        /**********************************************************************/
        LoopContainerType& periodicLoops()
        {
            return *this;
        }
        
        /**********************************************************************/
        const LoopContainerType& periodicLoops() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addPeriodicLoop(PeriodicLoopType* const pL)
        {
            const bool success=this->insert(std::make_pair(pL->sID,pL)).second;
            assert(success && "Could not insert in loop");
        }
        
        /**********************************************************************/
        void removePeriodicLoop(const PeriodicLoopType* const pL)
        {
            const size_t erased=this->erase(pL->sID);
            assert(erased==1 && "Could not erase loop");
        }
        
    };
    
}
#endif
