/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetworkRemesh_hpp_
#define model_DislocationNetworkRemesh_hpp_

#include <DislocationNetworkRemesh.h>
#include <TextFileParser.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    DislocationNetworkRemesh<DislocationNetworkType>::DislocationNetworkRemesh(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    /* init */,Lmax(TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmax",true)*minMeshSize(DN.mesh))
    /* init */,Lmin(TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmin",true)*minMeshSize(DN.mesh))
    /* init */,relativeAreaThreshold(TextFileParser("inputFiles/DD.txt").readScalar<double>("relativeAreaThreshold",true))
    /* init */,remeshFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("remeshFrequency",true))
    {
        
        assert(Lmin<=Lmax);
        assert(Lmax>3.0*Lmin);
        assert(Lmin>=0.0);
        assert(relativeAreaThreshold>=0.0);
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::minMeshSize(const SimplicialMesh<dim>& mesh)
    {
        return std::min(mesh.xMax(0)-mesh.xMin(0),std::min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remesh(const long int& runID)
    {/*! Performs remeshByContraction and then remeshByExpansion.
      * This order guarantees that 2-vertex NetworkComponents are expanded.
      */
        if (remeshFrequency)
        {
            if(!(runID%remeshFrequency))
            {
                remeshByRemoval();
                //                    remeshByContraction();
                remeshByExpansion();
                contract0chordSegments();
                remove0AreaLoopAcrossBnd();
                contractBoundaryNodes();//Added by Yash
            }
        }
    }
    
    /**************************************************************************/
    //Working Version of remesh by removal (In accordance with the commented version of DislocationLoopNode is Removable)
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remeshByRemoval()
    {
        const auto t0 = std::chrono::system_clock::now();
        std::cout << "        Remeshing network: removing... " << std::flush;
        DN.danglingBoundaryLoopNodes.clear();
        
        std::vector<std::shared_ptr<typename DislocationNetworkType::LoopNodeType>> tempLoopNodes;
        for (const auto &loopNode : DN.loopNodes())
        {
            tempLoopNodes.push_back(loopNode.second.lock());
        }
        std::set<size_t> bndLoopNodes;
        size_t Nremoved = 0;
        std::set<const typename DislocationNetworkType::LoopNodeType*> removedLoopNodes;
        
        // std::deque<size_t> toBeRemoved;
        for (const auto &loopNode : DN.loopNodes())
        {
            //                std::cout<<"node "<<node.second->sID<<" "<<node.second->isSimpleBoundaryNode()<<" "<<node.second->isSimpleGrainBoundaryNode()<<std::endl;
            if (!loopNode.second.lock()->networkNode->masterNode && removedLoopNodes.find(loopNode.second.lock().get())==removedLoopNodes.end())
            {
                
                const auto isRemovableLoopNode(loopNode.second.lock()->isRemovable(Lmin, relativeAreaThreshold));
                // std::cout<<" Trying to remove "<<loopNode.second.lock()->tag()<<" => "<<isRemovableLoopNode.first<<std::endl;
                if (isRemovableLoopNode.first)
                {
                    //Do this for all the loopNode of the networkNode
                    if (isRemovableLoopNode.second==loopNode.second.lock()->sID)
                    {
                        //Only the current node is needed to be removed
                        if (!loopNode.second.lock()->periodicPlaneEdge) //No need to update the boundary nodes for the nodes at the boundary
                        {
                            for (const auto &bndNode : loopNode.second.lock()->boundaryPrev())
                            {
                                // std::cout << " For " << loopNode.second.lock()->tag() << " bndPrev is " << bndNode->tag() << std::endl;
                                bndLoopNodes.insert(bndNode->sID);
                            }
                            for (const auto &bndNode : loopNode.second.lock()->boundaryNext())
                            {
                                // std::cout << " For " << loopNode.second.lock()->tag() << " bndNext is " << bndNode->tag() << std::endl;
                                bndLoopNodes.insert(bndNode->sID);
                            }
                        }
                        
                        // std::cout<<" Starting to remove the loop Node "<<std::endl;
                        
                        DN.removeLoopNode(loopNode.second.lock()->sID);
                        // std::cout<<" Removed "<<loopNode.second.lock()->tag()<<std::endl;
                        
                        removedLoopNodes.insert(loopNode.second.lock().get());
                        Nremoved++;
                    }
                    else
                    {
                        //Both the current and the twinned loopNode are needed to be removed
                        
                        //Collect corresponding to the twinned loopNode and current loopNode
                        for (const auto& loopN : loopNode.second.lock()->networkNode->loopNodes())
                        {
                            if (loopN->sID==loopNode.second.lock()->sID || loopN->sID==isRemovableLoopNode.second)
                            {
                                if (!loopN->periodicPlaneEdge)
                                {
                                    for (const auto &bndNode : loopN->boundaryPrev())
                                    {
                                        bndLoopNodes.insert(bndNode->sID);
                                    }
                                    for (const auto &bndNode : loopN->boundaryNext())
                                    {
                                        bndLoopNodes.insert(bndNode->sID);
                                    }
                                }
                                DN.removeLoopNode(loopN->sID);
                                removedLoopNodes.insert(loopN);
                                Nremoved++;
                            }
                        }
                    }
                }
                //                else
                //                {
                //                    const auto pPrev(loopNode.second.lock()->periodicPrev());
                //                    const auto pNext(loopNode.second.lock()->periodicNext());
                //                    if (pPrev && pNext)
                //                    {
                //                        if (pPrev->networkNode==loopNode.second.lock()->networkNode && pNext->networkNode==loopNode.second.lock()->networkNode)
                //                        {
                //                            DN.removeLoopNode(loopNode.second.lock()->sID);
                //                            removedLoopNodes.insert(loopNode.second.lock().get());
                //                            Nremoved++;
                //
                //                        }
                //                    }
                //                }
                
            }
        }
        //
        // for (const auto &nodeID : toBeRemoved)
        // {
        //     DN.removeLoopNode(nodeID);
        //     Nremoved++;
        // }
        tempLoopNodes.clear(); //Very important to clear the loop node here
        std::cout << " (" << Nremoved << " LoopNodes removed) " << std::flush;
        Nremoved=0;
        for (const size_t &nodeID : bndLoopNodes)
        {
            auto nodeIter(DN.loopNodes().find(nodeID));
            if (nodeIter != DN.loopNodes().end())
            {
                auto bndNode(nodeIter->second.lock());
                
                const auto pPrev(bndNode->periodicPrev());
                const auto pNext(bndNode->periodicNext());
                // std::cout<< "pPrev remesh= " << pPrev->tag() << " @ " << pPrev->get_P().transpose() <<"with size "<<pPrev->networkNode->loopNodes().size() << std::endl;
                // std::cout<< "pNext remesh= " << pNext->tag() << " @ " << pNext->get_P().transpose()  <<"with size "<<pNext->networkNode->loopNodes().size() << std::endl;
                if (pPrev && pNext)
                {
                    const auto pPrevLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                    const auto pNextLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                    //
                    
                    SegmentSegmentDistance<dim - 1> ssd(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge->source, *bndNode->periodicPlaneEdge->sink);
                    if (ssd.dMin < FLT_EPSILON)
                    {
                        //                    VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                        const VectorDim tempShift(bndNode->loop()->periodicGlidePlane->referencePlane->globalPosition(0.5 * (ssd.x0 + ssd.x1)));
                        // std::cout<<std::scientific<<std::setprecision(15)<<" bndNode tag "<<bndNode->tag()<<" ==> "<<bndNode->get_P().transpose()<<" Updated Position "<<tempShift.transpose()<<std::endl;
                        // std::cout<<"From 3D intersection "<<(0.5 * (ssd3.x0 + ssd3.x1)-bndNode->periodicPlanePatch()->shift).transpose()<<std::endl;
                        bndNode->set_P(VectorLowerDim(0.5 * (ssd.x0 + ssd.x1)));
                    }
                    
                    else
                    {
                        bndNode->set_P(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                        DN.danglingBoundaryLoopNodes.insert(bndNode.get());
                    }
                }
                else
                {
                    //Sole boundary node as the survivor..Remove them
                    DN.removeLoopNode(bndNode->sID);
                    Nremoved++;
                }
                
                
            }
        }
        
        std::cout << " (" << Nremoved << " bndLoopNodes removed). Updating boundary nodes " << std::flush;
        
        if (DN.danglingBoundaryLoopNodes.size())
        {
            DN.updateBoundaryNodes();
        }
        
        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]." << defaultColor << std::endl;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remeshByExpansion()
    {
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"        Remeshing network: expanding... "<<std::flush;
        std::set<std::pair<size_t,size_t> > toBeExpanded;
        for (const auto& linkIter : DN.networkLinks())
        {
            const auto sharedLink(linkIter.second.lock());
            
            if( !sharedLink->hasZeroBurgers()
               && !sharedLink->isSessile()
               && !sharedLink->isBoundarySegment()
               && !sharedLink->isGrainBoundarySegment()
               && !sharedLink->isVirtualBoundarySegment()
               )
            {
                const VectorDim chord(sharedLink->chord()); // this is sink->get_P() - source->get_P()
                const double chordLength(chord.norm());
                
                if (sharedLink->source->loopNodes().size()>1 && sharedLink->sink->loopNodes().size()>1
                    /*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/ && chordLength>3.0*Lmin)
                { // also expands a straight line to generate glissile segment
                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID,sharedLink->sink->sID));
                    
                }
                
                // Expand segments shorter than Lmax
                if (chordLength>Lmax)
                {
                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID,sharedLink->sink->sID));
                }
            }
            
        }
        
        
        // Call Network::expand
        unsigned int Nexpanded(0);
        const double expand_at(0.5);
        for (const auto& expIter : toBeExpanded)
        {
            //            const size_t i(expIter->first);
            //            const size_t j(expIter->second);
            const auto source(DN.networkNodes().get(expIter.first ));
            const auto sink  (DN.networkNodes().get(expIter.second));
            
            const auto Lij(DN.networkLinks().get(expIter));
            if(Lij)
            {
                VectorDim expandPoint(Lij->get_r(expand_at));
                
                
                if(  (expandPoint-source->get_P()).squaredNorm() > FLT_EPSILON
                   &&(expandPoint-  sink->get_P()).squaredNorm() > FLT_EPSILON)
                {
                    
                    const auto newNetNode(DN.networkNodes().create(expandPoint,0.5*(source->get_V()+sink->get_V()),0.5*(source->velocityReduction()+sink->velocityReduction())));
                    DN.expandNetworkLink(Lij,newNetNode);
                    Nexpanded++;
                }
            }
        }
        std::cout<<" ("<<Nexpanded<<" expanded)"<<std::flush;
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::contractBoundaryNodes()
    {
        std::cout<<"Contracting BoundaryNodes... "<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        typedef std::pair<std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>,std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>> NetNodePairType;
        typedef std::set<NetNodePairType> BndSetType;
        
        std::map<NetNodePairType,BndSetType> bndNodesContractionMap;
        
        for (const auto& netNode : DN.networkNodes())
        {
            if (netNode.second.lock()->isBoundaryNode())
            {
                //                std::cout<<"NetNode is "<<netNode.second.lock()->tag()<<std::endl;
                
                std::set<std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>> internalNodes; // used to be regular ptrs
                for(const auto& loopNode : netNode.second.lock()->loopNodes())
                {
                    internalNodes.insert(loopNode->periodicPrev()->networkNode);
                    internalNodes.insert(loopNode->periodicNext()->networkNode);
                }
                
                if(internalNodes.size()==2)
                {
                    
                    //                    std::cout<<"internal nodes are "<<(*internalNodes.begin())->tag()<<","<<(*internalNodes.rbegin())->tag()<<std::endl;
                    
                    const NetNodePairType key((*internalNodes.begin())->sID<(*internalNodes.rbegin())->sID? std::make_pair(*internalNodes.begin(),*internalNodes.rbegin()) : std::make_pair(*internalNodes.rbegin(),*internalNodes.begin())); // Line that was causing problem
                    
                    if(bndNodesContractionMap.find(key)==bndNodesContractionMap.end())
                    {// key not found
                        std::map<std::set<const PlanarMeshFace<dim>*>,std::set<typename DislocationNetworkType::NetworkNodeType*>> bndNodeMap;
                        
                        for(const auto& loopNode : (*internalNodes.begin())->loopNodes())
                        {
                            if(loopNode->periodicPrev()->networkNode==*internalNodes.rbegin())
                            {
                                for(const auto& bndNode : loopNode->boundaryPrev())
                                {
                                    //                                                        std::cout<<"adding bndNode "<<bndNode->networkNode->tag()<<std::endl;
                                    bndNodeMap[bndNode->networkNode->meshFaces()].insert(bndNode->networkNode.get());
                                }
                            }
                            else if(loopNode->periodicNext()->networkNode==*internalNodes.rbegin())
                            {
                                for(const auto& bndNode : loopNode->boundaryNext())
                                {
                                    //                                                        std::cout<<"adding bndNode "<<bndNode->networkNode->tag()<<std::endl;
                                    bndNodeMap[bndNode->networkNode->meshFaces()].insert(bndNode->networkNode.get());
                                }
                            }
                        }
                        
                        for(const auto& pair : bndNodeMap)
                        {
                            //                                                std::cout<<"Nodes of face"<<std::endl;
                            //                                                for(const auto& face : pair.first)
                            //                                                {
                            //                                                    std::cout<<face->outNormal().transpose()<<std::endl;
                            //                                                }
                            
                            //                                                std::cout<<"nodes size="<<pair.second.size()<<std::endl;
                            if(pair.second.size()>1)
                            {
                                for(const auto& tempNode : pair.second)
                                {
                                    //                                                        std::cout<<"node size="<<tempNode->sID<<std::endl;
                                    if(tempNode!=*pair.second.begin())
                                    {
                                        if((tempNode->get_P()-(*pair.second.begin())->get_P()).norm()<FLT_EPSILON)
                                        {
                                            const NetNodePairType contractPair(tempNode->sID<(*pair.second.begin())->sID? std::make_pair(tempNode,*pair.second.begin()) : std::make_pair(*pair.second.begin(),tempNode));
                                            //                                                                std::cout<<"contract pair="<<contractPair.first->sID<<" "<<contractPair.second->sID<<std::endl;
                                            bndNodesContractionMap[key].insert(contractPair);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
//        std::cout<<"Finished selection"<<std::endl;
        
        for(const auto& pair : bndNodesContractionMap)
        {
            std::cout<<"Contracting bnd nodes of intenal nodes "<<pair.first.first->tag()<<" "<<pair.first.second->tag()<<std::endl;
            for(const auto& bndPair : pair.second)
            {
                std::cout<<"Contracting bnd nodes  "<<bndPair.first->tag()<<" "<<bndPair.second->tag()<<std::endl;
                DN.contract(bndPair.first,bndPair.second);
            }
        }
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }
    
    //     template <typename DislocationNetworkType>
    //     void DislocationNetworkRemesh<DislocationNetworkType>::contractBoundaryNodes()
    //     {
    //     const auto t0= std::chrono::system_clock::now();
    //     std::cout<<"        Remeshing network: Contracting BoundaryNodes... "<<std::flush;
    //     std::set<std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>> bndNodes;
    //     size_t nContracted(0);
    //         for (const auto& netNode : DN.networkNodes())
    //         {
    //             if (netNode.second.lock()->isBoundaryNode())
    //             {
    //                 // std::cout<<" insertryy "<<netNode.second.lock()->sID<<std::endl;
    //
    //                 bndNodes.insert(netNode.second.lock());
    //             }
    //
    //         }
    //     // std::cout<<" bnd node size "<<bndNodes.size()<<std::endl;
    //     std::deque<std::pair<std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>,std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>>> toBeRemoved;
    
    //     for (size_t i = 0; i < bndNodes.size(); ++i)
    //     {
    //         auto linkIterI(bndNodes.begin());
    //         std::advance(linkIterI, i);
    //         for (size_t j=i+1; j<bndNodes.size(); ++j)
    //         {
    //             auto linkIterJ(bndNodes.begin());
    //             std::advance(linkIterJ, j);
    //             if (((*linkIterI)->get_P()-(*linkIterJ)->get_P()).squaredNorm()<FLT_EPSILON)
    //             {
    //                 //Get the neighbors of I
    //                 bool temp=false;
    //                 for (const auto& iNeigh : (*linkIterI)->neighbors())
    //                 {
    //                     for (const auto& jNeigh : (*linkIterJ)->neighbors())
    //                     {
    //                         if (std::get<0>(iNeigh.second)==std::get<0>(jNeigh.second) &&
    //                         (!std::get<1>(iNeigh.second)->hasZeroBurgers() && !std::get<1>(jNeigh.second)->hasZeroBurgers()))
    //                         {
    //                             temp=true;
    //                             break;
    //                         }
    //                     }
    //                     if (temp)
    //                     {
    //                         break;
    //                     }
    
    //                 }
    //                 if (temp)
    //                 {
    //                     toBeRemoved.push_back(std::make_pair(*linkIterI,*linkIterJ));
    //                 }
    //             }
    //         }
    //     }
    
    //     // std::cout<<" to be removed bnd node size "<<bndNodes.size()<<std::endl;
    
    
    //     for (const auto& nodeCont : toBeRemoved)
    //     {
    //         const auto success(DN.contract(nodeCont.first,nodeCont.second));
    //         nContracted+=success;
    //     }
    
    //     std::cout << "(" << nContracted << " contracted)" << std::flush;
    //     std::cout << magentaColor << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]" << defaultColor << std::endl;
    // }
    
    
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::contract0chordSegments()
    {
        std::cout<<"        Contracting zero-chord segments... "<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
        
        //            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
        for (const auto& linkIter : DN.networkLinks())
        {
            VectorDim chord(linkIter.second.lock()->chord()); // this is sink->get_P() - source->get_P()
            const double chordLength(chord.norm());
            //    if (!linkIter.second.lock()->source->isBoundaryNode() && !linkIter.second.lock()->sink->isBoundaryNode() && !linkIter.second.lock()->hasZeroBurgers())
            if (!linkIter.second.lock()->source->isBoundaryNode() && !linkIter.second.lock()->sink->isBoundaryNode() )
            {
                //    if (chordLength <= FLT_EPSILON && (linkIter.second.lock()->source->get_V() - linkIter.second.lock()->sink->get_V()).squaredNorm() < FLT_EPSILON)
                //Removing the second condition is important (both the nodes may not have the same velocity due to complex network formation)
                //If the two nodes are not connected to any other links then this may be true otherwise it is not
                if (chordLength <= FLT_EPSILON)
                { // toBeContracted part
                    toBeContracted.insert(std::make_pair(chordLength, std::make_pair(linkIter.second.lock()->source->sID, linkIter.second.lock()->sink->sID)));
                }
                else
                {
                    const VectorDim velChange(linkIter.second.lock()->sink->get_V() - linkIter.second.lock()->source->get_V());
                    const double velChangeNorm(velChange.norm());
                    if (velChangeNorm > FLT_EPSILON)
                    {
                        const double dotPTemp((chord / chordLength).dot(velChange / velChangeNorm));
                        // std::cout<<"For "<<linkIter.second.lock()->tag()<<"dotP Temp is "<<dotPTemp<<"=>"<<((velChangeNorm * DN.simulationParameters.dt) >= chordLength)<<std::endl;
                        // if (fabs(dotPTemp + 1) < FLT_EPSILON)
                        if (dotPTemp  > FLT_EPSILON)
                        {
                            if ((velChangeNorm * DN.simulationParameters.dt) >= chordLength)
                            {
                                toBeContracted.insert(std::make_pair(chordLength, std::make_pair(linkIter.second.lock()->source->sID, linkIter.second.lock()->sink->sID)));
                            }
                        }
                    }
                }
            }
        }
        
        // Call Network::contract
        unsigned int Ncontracted(0);
        for (const auto& smallIter : toBeContracted)
        {
            const size_t i(smallIter.second.first);
            const size_t j(smallIter.second.second);
            const auto Lij(DN.networkLinks().get(std::make_pair(i,j)));
            
            if (Lij)
            {
                std::set<size_t> bndLoopNodes;
                for (const auto &loopN : Lij->source->loopNodes())
                {
                    for (const auto &bndNode : loopN->boundaryPrev())
                    {
                        bndLoopNodes.insert(bndNode->sID);
                    }
                    for (const auto &bndNode : loopN->boundaryNext())
                    {
                        bndLoopNodes.insert(bndNode->sID);
                    }
                }
                for (const auto &loopN : Lij->sink->loopNodes())
                {
                    for (const auto &bndNode : loopN->boundaryPrev())
                    {
                        bndLoopNodes.insert(bndNode->sID);
                    }
                    for (const auto &bndNode : loopN->boundaryNext())
                    {
                        bndLoopNodes.insert(bndNode->sID);
                    }
                }
                
                Ncontracted+=DN.contract(Lij->source,Lij->sink);
                
                for (const size_t &nodeID : bndLoopNodes)
                {
                    auto nodeIter(DN.loopNodes().find(nodeID));
                    if (nodeIter != DN.loopNodes().end())
                    {
                        auto bndNode(nodeIter->second.lock());
                        
                        const auto pPrev(bndNode->periodicPrev());
                        const auto pNext(bndNode->periodicNext());
                        if (pPrev && pNext)
                        {
                            const auto pPrevLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                            const auto pNextLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                            
                            SegmentSegmentDistance<dim - 1> ssd(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge->source, *bndNode->periodicPlaneEdge->sink);
                            // /*  Giacomo's Version
                            if (ssd.dMin < FLT_EPSILON)
                            {
                                //                    VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                                bndNode->set_P(VectorLowerDim(0.5 * (ssd.x0 + ssd.x1)));
                            }
                            
                            else
                            {
                                bndNode->set_P(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                                DN.danglingBoundaryLoopNodes.insert(bndNode.get());
                            }
                        }
                    }
                }
            }
        }
        std::cout<<"("<<Ncontracted<<" contracted)"<<std::flush;
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
        std::cout << "Updating Boundary Nodes after 0chord contraction" << std::endl;
        if (DN.danglingBoundaryLoopNodes.size())
        {
            DN.updateBoundaryNodes();
        }
    }
    
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remove0AreaLoopAcrossBnd()
    {
        std::cout << "        Removing zero-area loops across boundary... " << std::flush;
        const auto t0 = std::chrono::system_clock::now();
        size_t nRemoved(0);
        
        //Snippet to remove zero area loops across boundary
        if (DN.simulationParameters.isPeriodicSimulation())
        {
            std::vector<size_t> loopToBeRemoved;
            for (const auto &weakLoop : DN.loops())
            {
                const auto loop(weakLoop.second.lock());
                if (loop->glidePlane)
                {
                    const double slippedArea(loop->slippedArea());
                    if (slippedArea < FLT_EPSILON)
                    {
                        std::vector<const typename DislocationNetworkType::LoopNodeType *> loopNodesPos;
                        std::vector<const typename DislocationNetworkType::LoopNodeType *> bndNodesMap;
                        for (const auto &loopLink : loop->linkSequence())
                        {
                            if (!loopLink->source->periodicPlaneEdge)
                            {
                                loopNodesPos.emplace_back(loopLink->source.get());
                            }
                            else
                            {
                                bndNodesMap.emplace_back(loopLink->source.get());
                            }
                        }
                        if (loopNodesPos.size() < 3 && bndNodesMap.size()>0)
                        {
                            loopToBeRemoved.push_back(loop->sID);
                        }
                    }
                }
            }
            
            //Remove the loop
            for (const auto &loopID : loopToBeRemoved)
            {
                DN.deleteLoop(loopID);
                nRemoved++;
            }
        }
        std::cout << "(" << nRemoved << " Loops Removed)" << std::flush;
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
    }
    
    
    
} // namespace model
#endif

