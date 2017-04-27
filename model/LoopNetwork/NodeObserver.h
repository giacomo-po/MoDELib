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
        
        typedef std::map<size_t,const NodeType* const> NodeContainerType;
        
    private:
        
        static NodeContainerType nodeMap;
        
    public:
        
        /**********************************************************************/
        static NodeContainerType& nodes()
        {
            return nodeMap;
        }
        
        /**********************************************************************/
        static void addNode(const NodeType* const pL)
        {
            const bool success=nodeMap.insert(std::make_pair(pL->sID,pL)).second;
            assert(success && "Could not insert in NodeMap");
        }
        
        /**********************************************************************/
        static void removeNode(const NodeType* const pL)
        {
            const size_t erased=nodeMap.erase(pL->sID);
            assert(erased==1 && "Could not erase from NodeMap");
        }
        
    };
    
    template<typename NodeType>
    std::map<size_t,const NodeType* const> NodeObserver<NodeType>::nodeMap;
    
}
#endif
