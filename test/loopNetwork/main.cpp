/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>

//#include <TypeTraits.h>
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
    
        Dlink(const Dnode* const Ni,const Dnode* const Nj) : NetworkLink<Dlink>(Ni,Nj){}
        
    };
    
    
    struct Dloop : public Loop<Dloop>
    {
    
        Dloop(const Dnetwork& net) : Loop<Dloop>(net){}
    };


}

using namespace model;


int main()
{

    Dnetwork DN;
    
    
    for(int i=0;i<6;++i)
    {
        DN.insertDanglingNode();
//        nodeIDs.push_back(id);
    }
    
    std::vector<size_t> loop0={0,1,2,3};
    std::vector<size_t> loop1={0,3,4,5};
    
    
    DN.insertLoop(loop0);
    DN.insertLoop(loop1);
    DN.clearDanglingNodes();
    
    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
    
    std::cout<<"Expanding"<<std::endl;
    DN.expand(0,3);
    
    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
    
    DN.checkLoops();
    

    return 0;
    
}
