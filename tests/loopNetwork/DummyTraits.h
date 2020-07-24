/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <TypeTraits.h>

namespace model
{
    struct DummyNetwork;
        struct DummyNetworkNode;
        struct DummyNetworkLink;
    struct DummyLoop;
    struct DummyLoopNode;
    struct DummyLoopLink;
    
    template<>
    struct TypeTraits<DummyNetwork>
    {
        typedef DummyNetwork LoopNetworkType;
        //        typedef DummyNetworkNode  NetworkNodeType;
        typedef DummyLoop LoopType;
        typedef DummyLoopNode LoopNodeType;
        typedef DummyLoopLink LoopLinkType;
        
        typedef DummyNetworkNode NetworkNodeType;
        typedef DummyNetworkLink NetworkLinkType;

        
        typedef double FlowType;
        
        static constexpr FlowType zeroFlow=0.0;
    };
    
    template<>
    struct TypeTraits<DummyLoopNode> : public TypeTraits<DummyNetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<DummyLoopLink> : public TypeTraits<DummyNetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<DummyLoop> : public TypeTraits<DummyNetwork>
    {
        typedef DummyNetwork LoopNetworkType;
        
    };
    
    template<>
    struct TypeTraits<DummyNetworkNode> : public TypeTraits<DummyNetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<DummyNetworkLink> : public TypeTraits<DummyNetwork>
    {
        
    };
    
    //
    //    template<>
    //    struct TypeTraits<DummyLoopNode> : public TypeTraits<DummyNetwork>
    //    {
    //
    //    };
    
    //    template<>
    //    struct TypeTraits<DummyLoopLink> : public TypeTraits<DummyNetwork>
    //    {
    //
    //    };
    
    //    template<>
    //    struct TypeTraits<DLnode> : public TypeTraits<Dnetwork>
    //    {
    //
    //    };
}
