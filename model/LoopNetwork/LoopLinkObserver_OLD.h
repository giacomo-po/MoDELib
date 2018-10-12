/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopLinkObserver_H_
#define model_LoopLinkObserver_H_

//#include <iostream>
//#include <list>
#include <map>
#include <utility>

//#include <memory>
//#include <iterator>



namespace model
{
    template<typename LinkType>
    struct LoopLinkObserver : public std::set<const LoopLink<LinkType>*>
    {
        
//        typedef std::pair<bool,LinkType* const>			IsNetworkLinkType;
//        typedef std::pair<bool,const LinkType* const>	IsConstNetworkLinkType;
//        
//        
//        typedef std::map<std::pair<size_t,size_t>,LinkType* const> LinkContainerType;
//        
//        typedef std::pair<bool,LinkType* const>			IsNetworkEdgeType;
//        typedef std::pair<bool,const LinkType* const>	IsConstNetworkEdgeType;
//        
//        /**********************************************************************/
//        IsNetworkEdgeType link(const size_t & i, const size_t & j)
//        {
//            typename LinkContainerType::iterator edgeIter(this->find(std::make_pair(i,j)));
//            return (edgeIter==this->end())?  std::make_pair(false,static_cast<LinkType* const>(nullptr)) :
//            /*                              */ std::make_pair(true,edgeIter->second);
//        }
//        
//        /**********************************************************************/
//        IsConstNetworkEdgeType link(const size_t & i, const size_t & j) const
//        {
//            typename LinkContainerType::const_iterator edgeIter(this->find(std::make_pair(i,j)));
//            return (edgeIter==this->end())?  std::make_pair(false,static_cast<const LinkType* const>(nullptr)) :
//            /*                              */ std::make_pair(true,edgeIter->second);
//        }
        
        /**********************************************************************/
        const LoopLinkObserver<LinkType>& loopLinkObserver() const
        {
            return *this;
        }
        
//        /**********************************************************************/
//        const LinkContainerType& loopLinks() const
//        {
//            return *this;
//        }
//        
//        /**********************************************************************/
//        LinkContainerType& links()
//        {
//            return *this;
//        }
        
        /**********************************************************************/
        void addLink(const LoopLink<LinkType>* const pL)
        {
            const bool success=this->insert(pL).second;
            assert(success && "Could not insert in LoopLinkSet");
        }
        
        /**********************************************************************/
        void removeLink(LoopLink<LinkType>* const pL)
        {
            const size_t erased=this->erase(pL);
            assert(erased==1 && "Could not erase from LoopLinkSet");
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


