/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopLink_H_
#define model_LoopLink_H_

#include <algorithm>
#include <memory>
#include <string>
#include <iterator>
#include <TypeTraits.h>
#include <MPIcout.h>
//#include <NodeObserver.h>
#include <NetworkBase.h>
#include <NetworkLink.h>

#define VerboseLoopLink(N,x) if(verboseLevel>=N){model::cout<<blueColor<<x<<defaultColor;}

namespace model
{
    template<typename Derived>
    class LoopLink : public  CRTP<Derived>
    /*            */,public NetworkBase<typename TypeTraits<Derived>::LoopNetworkType,std::array<size_t,3>>

    {
        
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef NetworkBase<LoopNetworkType,std::array<size_t,3>> NetworkBaseType;
        typedef typename TypeTraits<Derived>::LoopNodeType LoopNodeType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::NetworkLinkType NetworkLinkType;
        typedef typename TypeTraits<Derived>::FlowType FlowType;
        
        
        
    public:
        /**********************************************************************/
        static typename NetworkBaseType::KeyType loopLinkKey(const std::shared_ptr<LoopType>& loop,
                                      const std::shared_ptr<LoopNodeType>& Ni,
                                      const std::shared_ptr<LoopNodeType>& Nj)
        {
            return typename NetworkBaseType::KeyType{loop->sID,std::min(Ni->sID,Nj->sID),std::max(Ni->sID,Nj->sID)};
        }
        
        static int verboseLevel;
        
        
        const std::shared_ptr<LoopNodeType> source;
        const std::shared_ptr<LoopNodeType> sink;
        const std::shared_ptr<LoopType> loop; // THIS SHOULD BE CONST, THAT WAY WE COULD USE A MAP WITH 3 IDs TO STORE LOOPLINKS

    private:
        std::shared_ptr<NetworkLinkType> _networkLink;
        

    public:

        Derived* prev;
        Derived* next;
        
//        friend void LoopNode<LoopNodeType>::addLoopLink(Derived* const);
//        friend void LoopNode<LoopNodeType>::removeLoopLink(Derived* const);

        
        /**********************************************************************/
        LoopLink(const std::shared_ptr<LoopNodeType>& so,
                 const std::shared_ptr<LoopNodeType>& si,
                 const std::shared_ptr<LoopType>& pL) :
        /* init */ NetworkBaseType(pL->p_network(),loopLinkKey(pL,so,si))
        /* init */,source(so)
        /* init */,sink(si)
        /* init */,loop(pL)
        /* init */,_networkLink(getNetworkLink())
        /* init */,prev(nullptr)
        /* init */,next(nullptr)
        {
            VerboseLoopLink(1,"Constructing LoopLink "<<tag()<<std::endl);
            source->addLoopLink(this->p_derived());
            sink->addLoopLink(this->p_derived());
            loop->addLoopLink(this->p_derived());
            if(_networkLink)
            {
                _networkLink->addLoopLink(this->p_derived());
            }
        }

        /**********************************************************************/
        ~LoopLink()
        {// call removeLoopLink in inverse order compared to addLoopLink
            VerboseLoopLink(1,"Destroying LoopLink "<<tag()<<std::endl);
            if(_networkLink)
            {
                _networkLink->removeLoopLink(this->p_derived());
            }
            loop->removeLoopLink(this->p_derived());
            sink->removeLoopLink(this->p_derived());
            source->removeLoopLink(this->p_derived());
        }
        
        const  std::shared_ptr<NetworkLinkType>& networkLink() const
        {
            return _networkLink;
        }
        
        std::shared_ptr<NetworkLinkType> getNetworkLink()
        {
            return this->derived().hasNetworkLink()? this->network().networkLinks().get(source->networkNode,sink->networkNode) : nullptr;
        }
        
        void resetNetworkLink()
        {
            if(_networkLink)
            {
                _networkLink->removeLoopLink(this->p_derived());
            }
            _networkLink=getNetworkLink();
            if(_networkLink)
            {
                _networkLink->addLoopLink(this->p_derived());
            }
        }
        
        /**********************************************************************/
        const FlowType& flow() const
        {/*!\returns a reference to the loop flow
          */
            return loop->flow();
        }
        
        /**********************************************************************/
        std::string tag() const
        {/*!\returns the string "i->j" where i is source->sID and j=sink->sID
          */
            return source->tag() + "->" + sink->tag() + " ("+std::to_string(loop->sID)+")";
        }
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const LoopLink<Derived>& ll)
        {
            os  << ll.loop->sID<<"\t"
            /**/<< ll.source->sID<<"\t"
            /**/<< ll.sink->sID<<"\t";
            
            return os;
        }
        
    };
    
    template<typename LinkType>
    int LoopLink<LinkType>::verboseLevel=0;
    
}
#endif
