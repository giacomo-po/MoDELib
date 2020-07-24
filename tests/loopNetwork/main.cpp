/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include <Eigen/Dense>
#include <StaticID.h>
#include <DummyTraits.h>
#include <LoopNetwork.h>



namespace model
{
    
    // .h
    
    struct DummyLoop : public Loop<DummyLoop>
    {
        
        DummyLoop(DummyNetwork* const net, const double& flow);
        std::shared_ptr<DummyLoop> clone() const;
    };
    
    struct DummyLoopNode : public LoopNode<DummyLoopNode>
    {
        
        DummyLoopNode(DummyNetwork* const net,
                      const std::shared_ptr<DummyLoop>&,
                      const std::shared_ptr<DummyNetworkNode>&);
        
        std::shared_ptr<DummyLoopNode> clone(const std::shared_ptr<DummyLoop>&,const std::shared_ptr<DummyNetworkNode>&) const;
    };

    
    struct DummyLoopLink: public LoopLink<DummyLoopLink>
    {

        DummyLoopLink(const std::shared_ptr<DummyLoopNode>& so,
                      const std::shared_ptr<DummyLoopNode>& si,
                      const std::shared_ptr<DummyLoop>& pL);
        
        bool hasNetworkLink() const;

    };
    
    
    
    struct DummyNetworkNode : public NetworkNode<DummyNetworkNode>
    {
        
        DummyNetworkNode(DummyNetwork* const net);
        
        std::shared_ptr<DummyNetworkNode> clone() const;
    };
    
    struct DummyNetworkLink : public NetworkLink<DummyNetworkLink>
    {
        
        DummyNetworkLink(DummyNetwork* const net,
                         const std::shared_ptr<DummyNetworkNode>&,
                         const std::shared_ptr<DummyNetworkNode>&
                         );
        
//        std::shared_ptr<DummyNetworkLink> clone() const;
    };
    
    struct DummyNetwork : public LoopNetwork<DummyNetwork>{};

    
    // .cpp

    DummyLoop::DummyLoop(DummyNetwork* const net, const double& flow) : Loop<DummyLoop>(net,flow)
    {
    }
    
    std::shared_ptr<DummyLoop> DummyLoop::clone() const
    {
        return std::shared_ptr<DummyLoop>(new DummyLoop(this->p_network(),this->flow()));
    }
    
    DummyLoopNode::DummyLoopNode(DummyNetwork* const net,
                                 const std::shared_ptr<DummyLoop>& loop,
                                 const std::shared_ptr<DummyNetworkNode>& networkNode) :
    /* init */ LoopNode<DummyLoopNode>(net,loop,networkNode)
    {
        //            std::cout<<"Creating DummyLoopNode "<<this->tag()<<std::endl;
    }
    
    std::shared_ptr<DummyLoopNode> DummyLoopNode::clone(const std::shared_ptr<DummyLoop>& otherLoop,
                                                        const std::shared_ptr<DummyNetworkNode>& otherNetworkNode) const
    {
//        return std::shared_ptr<DummyLoopNode>(new DummyLoopNode(this->p_network(),this->network().networkNodes().clone(networkNode->key)));
        return std::shared_ptr<DummyLoopNode>(new DummyLoopNode(this->p_network(),otherLoop,otherNetworkNode));
    }
    
    DummyLoopLink::DummyLoopLink(const std::shared_ptr<DummyLoopNode>& so,
                  const std::shared_ptr<DummyLoopNode>& si,
                  const std::shared_ptr<DummyLoop>& pL) :
    /* init */ LoopLink(so,si,pL)
    {// intersect existing patches, add as many networkLinks as # of intersected patches
        // fist netwrokNode being source->networkNode, last netwrokNode being sink->networkNode
        // and all nodes in the middle being new netwrokNodes

    }
    
    bool DummyLoopLink::hasNetworkLink() const
    {
        return true;
    }

    
    DummyNetworkNode::DummyNetworkNode(DummyNetwork* const net) :
    /* init */ NetworkNode<DummyNetworkNode>(net)
    {
        //            std::cout<<"Creating DummyLoopNode "<<this->tag()<<std::endl;
    }
    
    std::shared_ptr<DummyNetworkNode> DummyNetworkNode::clone() const
    {
        return std::shared_ptr<DummyNetworkNode>(new DummyNetworkNode(this->p_network()));
    }
    
    
    DummyNetworkLink::DummyNetworkLink(DummyNetwork* const net,
                                       const std::shared_ptr<DummyNetworkNode>& nI,
                                       const std::shared_ptr<DummyNetworkNode>& nJ) :
    /* init */ NetworkLink<DummyNetworkLink>(net,nI,nJ)
    {
        //            std::cout<<"Creating DummyLoopNode "<<this->tag()<<std::endl;
    }
    
//    std::shared_ptr<DummyNetworkLink> DummyNetworkLink::clone() const
//    {
//        assert(false && "undefined clone operation for DummyNetworkLink");
//        return std::shared_ptr<DummyNetworkLink>(nullptr);
//    }
    
}

using namespace model;


std::vector<std::shared_ptr<DummyLoopNode>> createNewLoop(DummyNetwork& DN, const size_t& nNodes, const double& flow)
{
    auto loop(DN.loops().create(flow));
    std::vector<std::shared_ptr<DummyLoopNode>> loopNodes;
    for(size_t k=0;k<nNodes;++k)
    {
        const auto networkNode(DN.networkNodes().create());
        //            networkNodes.push_back(networkNode);
        loopNodes.push_back(DN.loopNodes().create(loop,networkNode));
    }
    DN.insertLoop(loop,loopNodes);
    return loopNodes;
}

int main()
{
    DummyNetwork::verboseLevel=0;
    DummyLoopLink::verboseLevel=0;
    DummyLoopNode::verboseLevel=0;
    DummyLoop::verboseLevel=0;
    DummyNetworkNode::verboseLevel=0;
    DummyNetworkLink::verboseLevel=0;

    // Create the network
    DummyNetwork DN;
    

    
    auto loopNodes0(createNewLoop(DN,8,0.4));
    auto loopNodes1(createNewLoop(DN,8,0.4));

//    std::vector<std::shared_ptr<DummyNetworkNode>> networkNodes;

//    auto newNetworkNode(DN.networkNodes().create());
//    auto newLoopNode(DN.loopNodes().create(loop0,newNetworkNode));
    
    

    
    std::cout<<"BEFORE CONTRACTION"<<std::endl;
    DN.printLoops();
    
//    for(const auto& node : nodes)
//    {
//        std::cout<<"node "<<node->sID<<" count="<<node.use_count()<<std::endl;
//    }
    
    std::cout<<"AFTER CONTRACTION"<<std::endl;
//    DN.contractLoopNodes(loopNodes[0],loopNodes[4]);
    
    DN.contractNetworkNodes(loopNodes0[0]->networkNode,loopNodes0[1]->networkNode);

    
//        DN.cutLoop(nodes[0],nodes[4]);
//        DN.removeLoopNode(nodes[0]);
    
//    const auto newNetworkNode(DN.networkNodes().create());
//    const auto newLoopNode(DN.loopNodes().create(loop0,newNetworkNode));
//    DN.expandLoopLink(DN.loopLinks().begin()->second,newLoopNode);

//        const auto newNetworkNode(DN.networkNodes().create());
//        DN.expandNetworkLink(DN.networkLinks().begin()->second.lock(),newNetworkNode);

//    const auto networkLinks(DN.networkLinks());
//    for(const auto& networkLink : networkLinks)
//    {
//        if(!networkLink.second.expired())
//        {
//            auto curNetNode(DN.networkNodes().create());
//            DN.expandNetworkLink(networkLink.second.lock(),curNetNode);
//        }
//    }

//    DN.printLoopNodes();
//    DN.printNetworkNodes();

    
    
    DN.printLoops();

    
//    for(const auto& node : nodes)
//    {
//        std::cout<<"node "<<node->sID<<" count="<<node.use_count()<<std::endl;
//    }
    
//    DN.deleteLoop(0);
    std::cout<<"DONE"<<std::endl;

//    DN.connect(node0,node1,loop0);
//    DN.connect(node1,node0,loop0);
//    DN.disconnect(node1,node0,loop0);

//    DN.connect(node0,node1);
    
    //    loop0->printLoop();
//
//    // Create the nodes
//    int nNodes=20;
//    for(int i=0;i<nNodes;++i)
//    {
//        DN.insertDanglingNode();
//    }
//
//
////    std::vector<size_t> loop0={0,1,2,3,4,5,6,7};
//    std::vector<size_t> loop0={0,1,2,3,4,5,6,7,8,9,10,11};
//        std::vector<size_t> loop1={1,2,12,13,7,8,14,15};
////    std::vector<size_t> loop2={8,7,16,17,18,19};
//
//    //  std::vector<size_t> loop1={};
//
//    DN.insertLoop(loop0,0.2);
//    DN.insertLoop(loop1,0.3);
//   // DN.insertLoop(loop2,0.5);
//
//    //    DN.insertLoop(loop2,0.2);
//
//    DN.clearDanglingNodes();
//
//    DN.printLoops();
//
//
////    DN.merge(0,1,4,5);
////    DN.annihilateTransfer(1,2,8,7); // this creates a mess with loops
////    DN.transfer(1,2,8,7); // this creates a mess with loops
//
////    DN.cutLoop(0,1,8);
//
//        DN.contractSecond(1,8);
//    DN.contractSecond(2,7);
//
//    DN.printLoops();
//
// //   DN.contract(2,12);
////    DN.printLoops();
//
////    std::vector<size_t> loop1={0,3,4,5};
////    std::vector<size_t> loop2={6,7};
////
////
////    DN.insertLoop(loop0,0.2);
////    DN.insertLoop(loop1,4.5);
////    DN.insertLoop(loop2,2.3);
////
////
//
//    std::cout<<"# loops="<<DN.loops().size()<<std::endl;
//    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
//    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
//    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
//
//
//
//
////    std::cout<<"Expanding"<<std::endl;
////    DN.expand(3,4);
////
////    DN.contract(0,1);
////    DN.contract(0,2);
////
////    DN.contract(0,5);
////    DN.contract(0,4);
////    DN.contract(0,8);
////
////    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
////    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
////    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
////
////    DN.checkLoops();
//
//
////    for(const auto& node : DN.nodes())
////    {
////        std::cout<<"node "<<node.second->sID<<" neighbor-size="<<node.second->loopLinks().size()<<std::endl;
////    }
////
////    DN.nodes();


    

    return 0;
    
}
