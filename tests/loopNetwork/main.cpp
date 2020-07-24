/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include <Eigen/Dense>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/StaticID.h>


#include <iostream>

#include <Eigen/Dense>
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
    
    struct Dlink : public NetworkLink<Dlink>
    {
        
        Dlink(const std::shared_ptr<Dnode>& Ni,
              const std::shared_ptr<Dnode>& Nj) : NetworkLink<Dlink>(Ni,Nj){}
        
    };
    
    struct Dnetwork : public LoopNetwork<Dnetwork>{};
    
    struct Dnode : public LoopNode<Dnode>
    {
        
        Dnode(Dnetwork* const net) :
        /* init */LoopNode<Dnode>(net)
        {
            
        }
        
    };
    
    
    
    
    
    struct Dloop : public Loop<Dloop>
    {
        
        Dloop(Dnetwork* const net, const double& flow) : Loop<Dloop>(net,flow){}
    };
    
    
}

using namespace model;


int main()
{

    // Create the network
    Dnetwork DN;
    
    // Create the nodes
    int nNodes=20;
    for(int i=0;i<nNodes;++i)
    {
        DN.insertDanglingNode();
    }
    
    
//    std::vector<size_t> loop0={0,1,2,3,4,5,6,7};
    std::vector<size_t> loop0={0,1,2,3,4,5,6,7,8,9,10,11};
        std::vector<size_t> loop1={1,2,12,13,7,8,14,15};
//    std::vector<size_t> loop2={8,7,16,17,18,19};

    //  std::vector<size_t> loop1={};

    DN.insertLoop(loop0,0.2);
    DN.insertLoop(loop1,0.3);
   // DN.insertLoop(loop2,0.5);

    //    DN.insertLoop(loop2,0.2);
   
    DN.clearDanglingNodes();

    DN.printLoops();

    
//    DN.merge(0,1,4,5);
//    DN.annihilateTransfer(1,2,8,7); // this creates a mess with loops
//    DN.transfer(1,2,8,7); // this creates a mess with loops

//    DN.cutLoop(0,1,8);
    
        DN.contractSecond(1,8);
    DN.contractSecond(2,7);
    
    DN.printLoops();

 //   DN.contract(2,12);
//    DN.printLoops();

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
    
    
//    for(const auto& node : DN.nodes())
//    {
//        std::cout<<"node "<<node.second->sID<<" neighbor-size="<<node.second->loopLinks().size()<<std::endl;
//    }
//
//    DN.nodes();


    

    return 0;
    
}
