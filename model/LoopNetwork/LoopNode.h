/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopNode_H_
#define model_LoopNode_H_

#include <iostream>
//#include <list>
#include <set>
#include <memory>
#include <assert.h>

//#include <iterator>

#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
//#include <iomanip>
//#include <vector>
//#include <Eigen/Dense>
//#include <Eigen/StdVector>
#include <model/LoopNetwork/NodeObserver.h>
#include <model/LoopNetwork/LoopLink.h>

namespace model
{
    template<typename Derived>
    class LoopNode : public StaticID<Derived>,
    /*            */ public CRTP<Derived>,
    /*            */ private std::set<LoopLink<typename TypeTraits<Derived>::LinkType>*>
    {
        
        
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        
    public:
        
        /**********************************************************************/

                LoopNode(const LoopNode&) =delete;
//        const std::shared_ptr<LoopType>& pLoop;
    
        /**********************************************************************/
        LoopNode()
        {
            std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
            NodeObserver<Derived>::addNode(this->p_derived());
        }
        
        /**********************************************************************/
        ~LoopNode()
        {
//            std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
            NodeObserver<Derived>::removeNode(this->p_derived());
        }
    
        /**********************************************************************/
        const LoopLinkContainerType& loopLinks() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
            // Store the connecting link
            const bool inserted=this->insert(pL).second;
            assert(inserted);

            // form the prev/next structure
            for(const auto& other : loopLinks())
            {
                if(pL->pLoop.get()==other->pLoop.get() && pL!=other) //two links in same loop
                {
                    
                    if(pL->source.get()==other->sink.get())
                    {
                        pL->prev=other;
                        other->next=pL;
                    }
                    if (pL->sink.get()==other->source.get())
                    {
                        pL->next=other;
                        other->prev=pL;
                    }
                    if(pL->source.get()==other->source.get() || pL->sink.get()==other->sink.get())
                    {
                        assert(0 && "Not a Loop");
                    }
//                    else
//                    {
//                        assert(pL->source.get()!=other->source.get());
//                        assert(pL->sink.get()!=other->sink.get());
//                    }
                }

            }

            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            const int erased=this->erase(pL);
            assert(erased==1);
            
        }
        
        
//        void addNetworkLink(LinkType* const pL
        
    };
    
    
}
#endif
