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
#include <set>
#include <memory>
#include <assert.h>
#include <tuple>
#include <limits.h>

#include <StaticID.h>
#include <CRTP.h>

#include <LoopLink.h>
#include <NetworkNode.h>
#include <NetworkBase.h>



#ifndef NDEBUG
#define VerboseLoopNode(N,x) if(verboseLevel>=N){std::cout<<greenColor<<x<<defaultColor;}
#else
#define VerboseLoopNode(N,x)
#endif

namespace model
{
    
    template<typename Derived>
    class LoopNode : public StaticID<Derived>
    /*            */,public CRTP<Derived>
    /*            */,public NetworkBase<Derived,size_t>
    {
        
    public:
        
        typedef Derived LoopNodeType;
        typedef typename TypeTraits<Derived>::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::NetworkNodeType NetworkNodeType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef NetworkBase<Derived,size_t> NetworkBaseType;

        
    private:
        

        
        std::shared_ptr<LoopType> _loop;

        
    public:

        std::shared_ptr<NetworkNodeType> networkNode;
        std::pair<LoopNodeType*,LoopLinkType*> prev;
        std::pair<LoopNodeType*,LoopLinkType*> next;
        
        
        static int verboseLevel;
        
        
        
    
        /**********************************************************************/
        LoopNode(LoopNetworkType* const loopNetwork_in,
                 const std::shared_ptr<LoopType>& loop_in,
                 const std::shared_ptr<NetworkNodeType>& networkNode_in) :
        /* init */ NetworkBaseType(loopNetwork_in,&loopNetwork_in->loopNodes(),this->sID)
        /* init */,_loop(loop_in)
        /* init */,networkNode(networkNode_in)
        /* init */,prev(nullptr,nullptr)
        /* init */,next(nullptr,nullptr)
        {
            VerboseLoopNode(1,"Constructing LoopNode "<<tag()<<std::endl);

            networkNode->addLoopNode(this->p_derived());
        }
        
        /**********************************************************************/
        ~LoopNode()
        {
            VerboseLoopNode(1,"Destroying LoopNode "<<tag()<<std::endl);

            assert(prev.first==nullptr && "~LoopNode prev.first!=nullptr");
            assert(prev.second==nullptr && "~LoopNode prev.second!=nullptr");
            assert(next.first==nullptr && "~LoopNode next.first!=nullptr");
            assert(next.second==nullptr && "~LoopNode next.second!=nullptr");
            
            networkNode->removeLoopNode(this->p_derived());
        }
        
//        size_t networkID() const
//        {
//            return std::distance(this->network().loopNodes().begin(),this->network().loopNodes().find(this->key));
//        }
        
        /**********************************************************************/
        void resetLoop(const std::shared_ptr<LoopType>& newLoop)
        {
            VerboseLoopNode(1,"LoopNode "<<tag()<<" resetLoop to "<<newLoop->tag()<<std::endl);
            assert((isIsolated() || _loop==newLoop) && "Trying to change loop of LoopNode connected to LoopLinks");
            if(isIsolated() && _loop!=newLoop)
            {
                _loop=newLoop;
            }
        }
        
        /**********************************************************************/
        bool isIsolated() const
        {
            return prev.first==nullptr && next.first==nullptr;
        }
        
        /**********************************************************************/
        const std::shared_ptr<LoopType>& loop() const
        {
            return _loop;
        }
        
        /**********************************************************************/
        bool isContractableTo(const Derived* const other) const
        {
            return loop()==other->loop();
        }
        
        /**********************************************************************/
        size_t gID() const
        {/*!\returns The NetworkComponent::snID() of the component
          * containing this.
          */
            return this->network().globalNodeID(this->sID);
        }
    
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
            VerboseLoopNode(2,"LoopNode "<<tag()<<" addLoopLink "<<pL->tag()<<std::endl);
            VerboseLoopNode(3,"LoopNode "<<tag()<<" loop()= "<<loop()->tag()<<std::endl);
            VerboseLoopNode(3,"LoopLink "<<pL->tag()<<" loop= "<<pL->loop->tag()<<std::endl);

            assert(pL->loop==_loop && "LoopNode and LoopLink in different loops");

            if(pL->source.get()==this)
            {// outlink (next)
                assert((next.first==nullptr || next.first==pL->sink.get()) && "LoopNode::addLoopLink wrong next node ");
                assert((next.second==nullptr || next.second==pL) && "LoopNode::addLoopLink wrong next link");
                next.first=pL->sink.get();
                next.second=pL;
                
                assert((pL->prev==nullptr || pL->prev==prev.second));
                pL->prev=prev.second;
                if(prev.second)
                {
                    assert(prev.second->loop==pL->loop && "links in different loops");
                    prev.second->next=pL;
                }
            }
            else if(pL->sink.get()==this)
            {// inlink (prev)
                assert((prev.first==nullptr || prev.first==pL->source.get()) && "LoopNode::addLoopLink wrong prev node ");
                assert((prev.second==nullptr || prev.second==pL) && "LoopNode::addLoopLink wrong prev link");
                prev.first=pL->source.get();
                prev.second=pL;
                
                assert((pL->next==nullptr || pL->next==next.second));
                pL->next=next.second;
                if(next.second)
                {
                    assert(next.second->loop==pL->loop && "links in different loops");
                    next.second->prev=pL;
                }
            }
            else
            {
                assert(false && "LoopLink connected to wrong node.");
            }
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            VerboseLoopNode(2,"LoopNode "<<tag()<<" removeLoopLink "<<pL->tag()<<std::endl);
            assert(pL->loop==_loop && "LoopNode and LoopLink in different loops");

            
            if(pL->source.get()==this)
            {// outlink (next)
                VerboseLoopNode(3,"LoopNode "<<tag()<<" removeLoopLink (next) "<<pL->tag()<<std::endl);
                
                assert((next.first==nullptr || next.first==pL->sink.get()) && "LoopNode::addLoopLink wrong next node ");
                assert((next.second==nullptr || next.second==pL) && "LoopNode::addLoopLink wrong next link");
                next.first=nullptr;
                next.second=nullptr;
                
                assert((pL->prev==nullptr || pL->prev==prev.second));
                if(prev.second)
                {
                    prev.second->next=nullptr;
                }
            }
            else if(pL->sink.get()==this)
            {// inlink (prev)
                VerboseLoopNode(3,"LoopNode "<<tag()<<" removeLoopLink (prev) "<<pL->tag()<<std::endl);

                assert((prev.first==nullptr || prev.first==pL->source.get()) && "LoopNode::addLoopLink wrong prev node ");
                assert((prev.second==nullptr || prev.second==pL) && "LoopNode::addLoopLink wrong prev link");
                prev.first=nullptr;
                prev.second=nullptr;

                assert((pL->next==nullptr || pL->next==next.second));
                if(next.second)
                {
                    next.second->prev=nullptr;
                }
            }
            else
            {
                assert(false && "LoopLink connected to wrong node.");
            }

        }
        
        /**********************************************************************/
        // void updateConfinedGeometry()
        // {
        //     //This will change the confinement for the loop node in the derived layer
        //     //Added by Yash
        // }
        /**********************************************************************/
        std::string tag() const
        {/*!\returns the string "i" where i is this->sID
          */
            return std::to_string(this->sID) + "[" + (networkNode? networkNode->tag() : "?") + "]";
        }
        
    };
    
    template<typename Derived>
    int LoopNode<Derived>::verboseLevel=0;

    
    
}
#endif
