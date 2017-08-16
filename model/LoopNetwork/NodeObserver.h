/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NodeObserver_H_
#define model_NodeObserver_H_

#include <map>
#include <utility>

namespace model
{
    template<typename NodeType>
    class NodeObserver
    {

    public:
        
        typedef std::map<size_t,NodeType* const> NodeContainerType;
        
        typedef std::pair<bool,NodeType* const> IsNodeType;
        
//        typedef std::map<size_t,std::shared_ptr<NodeType> > SharedNodePtrMapType;
//
//
        typedef std::shared_ptr<NodeType> SharedNodePtrType;
        typedef std::pair<bool,SharedNodePtrType> IsSharedNodeType;

        
    private:
        
        static NodeContainerType nodeMap;

//        static SharedNodePtrMapType sharedNodePtrMap;
        
        
    public:
        
        
        /**********************************************************************/
        static IsNodeType node(const size_t& i)
        {
            typename NodeContainerType::const_iterator nodeIter(nodeMap.find(i));
            return (nodeIter==nodeMap.end())?  std::make_pair(false,static_cast<NodeType* const>(nullptr)) :
            /*                              */ std::make_pair(true,nodeIter->second);

        }
        
        /**********************************************************************/
        static IsSharedNodeType sharedNode(const size_t& i)
        {
            IsSharedNodeType temp=std::make_pair(false,std::shared_ptr<NodeType>(nullptr));

            
            const IsNodeType ni=node(i);
            if(ni.first)
            {
                if(ni.second->loopLinks().size())
                {
                    const auto pL=*ni.second->loopLinks().begin();
                    if(pL->source()->sID==i)
                    {
                        temp=std::make_pair(true,pL->source());
                    }
                    else if(pL->sink()->sID==i)
                    {
                        temp=std::make_pair(true,pL->sink());

                    }
                    else
                    {
                        assert(0 && "source or sink must be i");
                    }
                }
                
            }
            
            return temp;
            
//            typename SharedNodePtrMapType::const_iterator nodeIter(sharedNodePtrMap.find(i));
//            return (nodeIter==sharedNodePtrMap.end())?  std::make_pair(false,std::shared_ptr<NodeType>(nullptr)) :
//            /*                              */ std::make_pair(true,nodeIter->second);
//            
        }
        
        /**********************************************************************/
        static NodeContainerType& nodes()
        {
            return nodeMap;
        }
        
//        /**********************************************************************/
//        static SharedNodePtrMapType& sharedNodes()
//        {
//            return sharedNodePtrMap;
//        }
        
        /**********************************************************************/
        static void addNode(NodeType* const pL)
        {
            const bool success=nodeMap.insert(std::make_pair(pL->sID,pL)).second;
            assert(success && "Could not insert in NodeMap");
        }
        
        /**********************************************************************/
        static void removeNode(const NodeType* const pL)
        {
            const size_t erased=nodeMap.erase(pL->sID);
            assert(erased==1 && "Could not erase from NodeMap");

//            const size_t erased1=sharedNodePtrMap.erase(pL->sID);
//            assert(erased1==1 && "Could not erase from sharedNodePtrMap");

        }
        
//        /**********************************************************************/
//        static void addNode(const std::shared_ptr<NodeType>& pL)
//        {
//            sharedNodePtrMap.insert(std::make_pair(pL->sID,pL)).second;
////            assert(success && "Could not insert in NodeMap");
//        }
        
//        /**********************************************************************/
//        static void removeNode(const NodeType& pL)
//        {
//                const size_t erased=sharedNodePtrMap.erase(pL->sID);
//                assert(erased==1 && "Could not erase from NodeMap");
//        }
        
    };
    
    template<typename NodeType>
    std::map<size_t,NodeType* const> NodeObserver<NodeType>::nodeMap;

//    template<typename NodeType>
//    std::map<size_t,std::shared_ptr<NodeType> > NodeObserver<NodeType>::sharedNodePtrMap;

}
#endif
