/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>

//#include <Eigen/Dense>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/StaticID.h>


namespace model
{
    struct Dnetwork;
    struct Dnode;
//    struct DLnode;
    struct Dlink;
    struct Dloop;
    
    template<>
    struct TypeTraits<Dnetwork>
    {
        typedef Dnetwork LoopNetworkType;
        typedef Dnode  NodeType;
//        typedef DLnode LoopNodeType;
        typedef Dlink LinkType;
        typedef Dloop LoopType;
        typedef double FlowType;
        
        static constexpr FlowType zeroFlow=0.0;
    };
    
    template<>
    struct TypeTraits<Dnode> : public TypeTraits<Dnetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<Dlink> : public TypeTraits<Dnetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<Dloop> : public TypeTraits<Dnetwork>
    {
        
    };
    
//    template<>
//    struct TypeTraits<DLnode> : public TypeTraits<Dnetwork>
//    {
//        
//    };
}


#include <model/LoopNetwork/LoopNetwork.h>
//#include <Loop.h>



namespace model
{
    struct Dnetwork : public LoopNetwork<Dnetwork>{};
    
    struct Dnode : public LoopNode<Dnode>
    {
    
//    Node(const int& a)
//        {
//        }
        
    };
    
//    struct DLnode : public LoopNode<DLnode>
//    {
//        
//        //    Node(const int& a)
//        //        {
//        //        }
//        
//    };
    
    struct Dlink : public NetworkLink<Dlink>
    {
    
        Dlink(const std::shared_ptr<Dnode>& Ni,
              const std::shared_ptr<Dnode>& Nj) : NetworkLink<Dlink>(Ni,Nj){}
        
    };
    
    
    struct Dloop : public Loop<Dloop>
    {
    
        Dloop(const Dnetwork& net, const double& flow) : Loop<Dloop>(net,flow){}
    };


}

using namespace model;


int main()
{

    Dnetwork DN;
    
    int nNodes=8;
    for(int i=0;i<nNodes;++i)
    {
        DN.insertDanglingNode();
//        nodeIDs.push_back(id);
    }
    
    
//    std::vector<size_t> loop0={0,1,2,3,4,5,6,7};
    std::vector<size_t> loop0={0,1,2,3,4,};
        std::vector<size_t> loop1={5,1,2,7,6};
    //  std::vector<size_t> loop1={};

    DN.insertLoop(loop0,0.2);
    DN.insertLoop(loop1,0.2);
//    DN.insertLoop(loop2,0.2);
   
    DN.clearDanglingNodes();

    
//    DN.merge(0,1,4,5);
    DN.merge(0,1,5,1);
    
//    std::vector<size_t> loop1={0,3,4,5};
//    std::vector<size_t> loop2={6,7};
//    
//    
//    DN.insertLoop(loop0,0.2);
//    DN.insertLoop(loop1,4.5);
//    DN.insertLoop(loop2,2.3);
//
//
    
    std::cout<<"# loops="<<DN.loops().size()<<std::endl;
    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
//
//    std::cout<<"Expanding"<<std::endl;
//    DN.expand(3,4);
//    
//    DN.contract(0,1);
//    DN.contract(0,2);
//    
//    DN.contract(0,5);
//    DN.contract(0,4);
//    DN.contract(0,8);
//    
//    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
//    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
//    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
//    
//    DN.checkLoops();
    
    
    for(const auto& node : DN.nodes())
    {
        std::cout<<"node "<<node.second->sID<<" neighbor-size="<<node.second->loopLinks().size()<<std::endl;
    }

    DN.nodes();

    DN.printLoops();

    

    return 0;
    
}
