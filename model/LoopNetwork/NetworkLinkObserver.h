/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLinkObserver_H_
#define model_NetworkLinkObserver_H_

//#include <iostream>
//#include <list>
#include <map>
#include <utility>

//#include <memory>
//#include <iterator>



namespace model
{
    template<typename LinkType>
    struct NetworkLinkObserver : public std::map<std::pair<size_t,size_t>,LinkType* const>
    {
        
        typedef typename LinkType::NodeType NodeType;
        
        typedef std::pair<bool,LinkType* const>			IsNetworkLinkType;
        typedef std::pair<bool,const LinkType* const>	IsConstNetworkLinkType;
        
        
        typedef std::map<std::pair<size_t,size_t>,LinkType* const> LinkContainerType;
        
        typedef std::pair<bool,LinkType* const>			IsNetworkEdgeType;
        typedef std::pair<bool,const LinkType* const>	IsConstNetworkEdgeType;
        
        /**********************************************************************/
        IsNetworkEdgeType link(const size_t & i, const size_t & j)
        {
            typename LinkContainerType::iterator edgeIter(this->find(std::make_pair(i,j)));
            return (edgeIter==this->end())?  std::make_pair(false,static_cast<LinkType* const>(nullptr)) :
            /*                              */ std::make_pair(true,edgeIter->second);
        }
        
        /**********************************************************************/
        IsConstNetworkEdgeType link(const size_t & i, const size_t & j) const
        {
            typename LinkContainerType::const_iterator edgeIter(this->find(std::make_pair(i,j)));
            return (edgeIter==this->end())?  std::make_pair(false,static_cast<const LinkType* const>(nullptr)) :
            /*                              */ std::make_pair(true,edgeIter->second);
        }
        
        /**********************************************************************/
        const LinkContainerType& links() const
        {
            return *this;
        }
        
        /**********************************************************************/
        LinkContainerType& links()
        {
            return *this;
        }
        
        /**********************************************************************/
        void addLink(LinkType* const pL)
        {
            const bool success=this->insert(std::make_pair(std::make_pair(pL->source->sID,pL->sink->sID),pL)).second;
            assert(success && "Could not insert in linkMap");
        }
        
        /**********************************************************************/
        void removeLink(LinkType* const pL)
        {
            const size_t erased=this->erase(std::make_pair(pL->source->sID,pL->sink->sID));
            assert(erased==1 && "Could not erase from linkMap");
        }
        
        /**********************************************************************/
        std::shared_ptr<LinkType> sharedLink(const std::shared_ptr<NodeType>& nI, const std::shared_ptr<NodeType>& nJ) const
        {
        
            const std::pair<size_t,size_t> key=std::make_pair(std::min(nI->sID,nJ->sID),std::max(nI->sID,nJ->sID));
            const auto iterIJ(this->find(key));
            
            return iterIJ==this->end()? std::shared_ptr<LinkType>(new LinkType(nI,nJ)) : (*iterIJ->second->loopLinks().begin())->pLink;

        }

        
    };
        
}
#endif


//        /**********************************************************************/
//        static IsNetworkEdgeType link(const size_t & i, const size_t & j)
//        {
//            typename LinkContainerType::iterator edgeIter(this->find(std::make_pair(i,j)));
//            return (edgeIter==this->end())?  std::make_pair(false,static_cast<LinkType* const>(nullptr)) :
//            /*                              */ std::make_pair(true,edgeIter->second);
//        }

//        /**********************************************************************/
//        static IsConstNetworkEdgeType link(const size_t & i, const size_t & j)
//        {
//            typename LinkContainerType::const_iterator edgeIter(this->find(std::make_pair(i,j)));
//            return (edgeIter==this->end())?  std::make_pair(false,static_cast<const LinkType* const>(nullptr)) :
//            /*                              */ std::make_pair(true,edgeIter->second);
//        }


