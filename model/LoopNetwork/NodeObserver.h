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
    struct NodeObserver : public std::map<size_t,NodeType* const>
    {
        
        typedef std::map<size_t,NodeType* const> NodeContainerType;
        typedef std::pair<bool,NodeType* const> IsNodeType;
        typedef std::shared_ptr<NodeType> SharedNodePtrType;
        typedef std::pair<bool,SharedNodePtrType> IsSharedNodeType;

        /**********************************************************************/
        IsNodeType node(const size_t& i)
        {
            typename NodeContainerType::const_iterator nodeIter(this->find(i));
            return (nodeIter==this->end())?  std::make_pair(false,static_cast<NodeType* const>(nullptr)) :
            /*                              */ std::make_pair(true,nodeIter->second);

        }
        
        /**********************************************************************/
        IsSharedNodeType sharedNode(const size_t& i)
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
        }
        
        /**********************************************************************/
        const NodeContainerType& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addNode(NodeType* const pL)
        {
            const bool success=this->insert(std::make_pair(pL->sID,pL)).second;
            assert(success && "Could not insert in NodeMap");
        }
        
        /**********************************************************************/
        void removeNode(const NodeType* const pL)
        {
            const size_t erased=this->erase(pL->sID);
            assert(erased==1 && "Could not erase from NodeMap");

        }
        
    };
    
}
#endif
