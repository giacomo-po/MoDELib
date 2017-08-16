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
    class NetworkLinkObserver
    {
        
    public:
        typedef std::pair<bool,LinkType* const>			IsNetworkLinkType;
        typedef std::pair<bool,const LinkType* const>	IsConstNetworkLinkType;
        
        
        typedef std::map<std::pair<size_t,size_t>,LinkType* const> LinkContainerType;
        
        typedef std::pair<bool,LinkType* const>			IsNetworkEdgeType;
        typedef std::pair<bool,const LinkType* const>	IsConstNetworkEdgeType;
        
    private:
        
        static LinkContainerType linkMap;
        
    public:
        
        /**********************************************************************/
        static IsConstNetworkEdgeType link(const size_t & i, const size_t & j)
        {
            typename LinkContainerType::const_iterator edgeIter(linkMap.find(std::make_pair(i,j)));
            return (edgeIter==linkMap.end())?  std::make_pair(false,static_cast<const LinkType* const>(nullptr)) :
            /*                              */ std::make_pair(true,edgeIter->second);
        }
        
//        /**********************************************************************/
//        static IsConstNetworkEdgeType link(const size_t & i, const size_t & j)
//        {
//            typename LinkContainerType::const_iterator edgeIter(linkMap.find(std::make_pair(i,j)));
//            return (edgeIter==linkMap.end())?  std::make_pair(false,static_cast<const LinkType* const>(nullptr)) :
//            /*                              */ std::make_pair(true,edgeIter->second);
//        }
        
        /**********************************************************************/
        static LinkContainerType& links()
        {
            return linkMap;
        }
        
        /**********************************************************************/
        static void addLink(LinkType* const pL)
        {
            const bool success=linkMap.insert(std::make_pair(std::make_pair(pL->source->sID,pL->sink->sID),pL)).second;
            assert(success && "Could not insert in linkMap");
        }
        
        /**********************************************************************/
        static void removeLink(LinkType* const pL)
        {
            const size_t erased=linkMap.erase(std::make_pair(pL->source->sID,pL->sink->sID));
            assert(erased==1 && "Could not erase from linkMap");
        }
        
    };
    
    template<typename LinkType>
    std::map<std::pair<size_t,size_t>,LinkType* const> NetworkLinkObserver<LinkType>::linkMap;
    
}
#endif
