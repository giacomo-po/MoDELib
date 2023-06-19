/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetworkRemesh_cpp_
#define model_DislocationNetworkRemesh_cpp_

#include <DislocationNetworkRemesh.h>

namespace model
{
    /**************************************************************************/
    template <typename DislocationNetworkType>
    DislocationNetworkRemesh<DislocationNetworkType>::DislocationNetworkRemesh(DislocationNetworkType& DN_in):
    /* init */ DN(DN_in)
    /* init */,Lmax(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<double>("Lmax",true)*minMeshSize(DN.mesh))
    /* init */,Lmin(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<double>("Lmin",true)*minMeshSize(DN.mesh))
    /* init */,relativeAreaThreshold(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<double>("relativeAreaThreshold",true))
    /* init */,remeshFrequency(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("remeshFrequency",true))
    {
        
        assert(Lmin<=Lmax);
        assert(Lmax>3.0*Lmin);
        assert(Lmin>=0.0);
        assert(relativeAreaThreshold>=0.0);
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remesh(const long int &runID)
    { /*! Performs remeshByContraction and then remeshByExpansion.
       * This order guarantees that 2-vertex NetworkComponents are expanded.
       */
        if (remeshFrequency)
        {
            if (!(runID % remeshFrequency))
            {
                remeshByRemoval();
                //                    remeshByContraction();
                remeshByExpansion();
                contract0chordSegments();
                remove0AreaLoopAcrossBnd();
//                contractBoundaryNodes(); //Added by Yash
            }
        }
    }
    
    /**************************************************************************/
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

        std::set<std::tuple<size_t,size_t,size_t>> bndLoopNodes; //First is the ID of the loop node, second is the periodicPrev and third is the periodic Next
        //Definition of bndLoopNode as a set of tuple is very important as remesh may remove the loop nodes hence changing the connectivity of the loopnodes at the boundary
        //IF that happens the corresponding boundary node should be classified as a dangling loop node
        //Also believe that this is only required here as nothing else can change the configuration of the loops near the boundary
        size_t Nremoved = 0;
        std::set<const typename DislocationNetworkType::LoopNodeType*> removedLoopNodes;

        // std::deque<size_t> toBeRemoved;
        for (const auto &loopNode : DN.loopNodes())
        {
            //                std::cout<<"node "<<node.second->sID<<" "<<node.second->isSimpleBoundaryNode()<<" "<<node.second->isSimpleGrainBoundaryNode()<<std::endl;
//            if (!loopNode.second.lock()->networkNode->masterNode && removedLoopNodes.find(loopNode.second.lock().get())==removedLoopNodes.end())
            if (removedLoopNodes.find(loopNode.second.lock().get())==removedLoopNodes.end())
            {
                
                const auto isRemovableLoopNode(loopNode.second.lock()->isRemovable(Lmin, relativeAreaThreshold));
                // std::cout<<" Trying to remove "<<loopNode.second.lock()->tag()<<" => "<<isRemovableLoopNode.first<<std::endl;
                if (isRemovableLoopNode.first)
                {
                    //Do this for all the loopNode of the networkNode 
                    if (isRemovableLoopNode.second==loopNode.second.lock()->sID)
                    {
                        //Only the current node is needed to be removed
                        if (!loopNode.second.lock()->periodicPlaneEdge.first) //No need to update the boundary nodes for the nodes at the boundary
                        {
                            for (const auto &bndNode : loopNode.second.lock()->boundaryPrev())
                            {
                                // std::cout << " For " << loopNode.second.lock()->tag() << " bndPrev is " << bndNode->tag() << std::endl;
                                bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);
                            }
                            for (const auto &bndNode : loopNode.second.lock()->boundaryNext())
                            {
                                // std::cout << " For " << loopNode.second.lock()->tag() << " bndNext is " << bndNode->tag() << std::endl;
                                bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);
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
                                if (!loopN->periodicPlaneEdge.first)
                                {
                                    for (const auto &bndNode : loopN->boundaryPrev())
                                    {
                                        bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);

                                    }
                                    for (const auto &bndNode : loopN->boundaryNext())
                                    {
                                        bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);
                                    }
                                }
                                DN.removeLoopNode(loopN->sID);
                                removedLoopNodes.insert(loopN);
                                Nremoved++;
                            }
                        }
                    }
                }
            }
        }
        
        tempLoopNodes.clear(); //Very important to clear the loop node here
        std::cout << " (" << Nremoved << " LoopNodes removed) " << std::flush;
        Nremoved=0;
        for (const auto &nodeID : bndLoopNodes)
        {
            auto nodeIter(DN.loopNodes().find(std::get<0>(nodeID)));
            if (nodeIter != DN.loopNodes().end())
            {
                auto bndNode(nodeIter->second.lock());

                const auto pPrev(bndNode->periodicPrev());
                const auto pNext(bndNode->periodicNext());
                if (pPrev && pNext)
                {
                    if (pPrev->sID == std::get<1>(nodeID) && pNext->sID== std::get<2>(nodeID))
                    {
                        //If the periodic Prev and periodic next are the same as the remesh may change this
                        const auto pPrevLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                        const auto pNextLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                        //

                        SegmentSegmentDistance<dim - 1> ssd1(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                        if (ssd1.dMin < FLT_EPSILON)
                        {
                            //                    VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                            if (bndNode->periodicPlaneEdge.second)
                            {
                                SegmentSegmentDistance<dim - 1> ssd2(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.second->source, *bndNode->periodicPlaneEdge.second->sink);
                                const VectorLowerDim ssd1Pos(0.5 * (ssd1.x0 + ssd1.x1));
                                const VectorLowerDim ssd2Pos(0.5 * (ssd2.x0 + ssd2.x1));
                                assert((ssd1Pos - ssd2Pos).norm() < FLT_EPSILON && "The two positions must match ");
                            }
                            bndNode->set_P(VectorLowerDim(0.5 * (ssd1.x0 + ssd1.x1)));
                        }
                        else
                        {
                            bndNode->set_P(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                            DN.danglingBoundaryLoopNodes.insert(bndNode.get());
                        }
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
        const auto t0 = std::chrono::system_clock::now();
        std::cout << "        Remeshing network: expanding... " << std::flush;
        std::set<std::pair<size_t, size_t>> toBeExpanded;
        for (const auto &linkIter : DN.networkLinks())
        {
            const auto sharedLink(linkIter.second.lock());
            
            if (!sharedLink->hasZeroBurgers() && !sharedLink->isSessile() && !sharedLink->isBoundarySegment() && !sharedLink->isGrainBoundarySegment())
            {
                const VectorDim chord(sharedLink->chord()); // this is sink->get_P() - source->get_P()
                const double chordLength(chord.norm());
                
                if (sharedLink->source->loopNodes().size() > 1 && sharedLink->sink->loopNodes().size() > 1
                    /*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/
                    && chordLength > 3.0 * Lmin)
                { // also expands a straight line to generate glissile segment
                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID, sharedLink->sink->sID));
                }
                
                // Expand segments shorter than Lmax
                if (chordLength > Lmax)
                {
                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID, sharedLink->sink->sID));
                }
            }
        }
        
        // Call Network::expand
        unsigned int Nexpanded(0);
        const double expand_at(0.5);
        for (const auto &expIter : toBeExpanded)
        {
            //            const size_t i(expIter->first);
            //            const size_t j(expIter->second);
            const auto source(DN.networkNodes().get(expIter.first));
            const auto sink(DN.networkNodes().get(expIter.second));
            
            const auto Lij(DN.networkLinks().get(expIter));
            if (Lij)
            {
                VectorDim expandPoint(Lij->get_r(expand_at));
                
                if ((expandPoint - source->get_P()).squaredNorm() > FLT_EPSILON && (expandPoint - sink->get_P()).squaredNorm() > FLT_EPSILON)
                {
                    
                    const auto newNetNode(DN.networkNodes().create(expandPoint, 0.5 * (source->get_V() + sink->get_V()), 0.5 * (source->velocityReduction() + sink->velocityReduction())));
                    DN.expandNetworkLink(Lij, newNetNode);
                    Nexpanded++;
                }
            }
        }
        std::cout << " (" << Nexpanded << " expanded)" << std::flush;
        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]." << defaultColor << std::endl;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::contractBoundaryNodes() //This function contracts boundary nodes if they are at the same position and share atleast one same neightbors
    {
        std::cout << "Contracting BoundaryNodes... " << std::flush;
        const auto t0 = std::chrono::system_clock::now();
        
        typedef std::pair<std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>, std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>> NetNodePairType;
        typedef std::set<NetNodePairType> BndSetType;
        
        std::map<NetNodePairType, BndSetType> bndNodesContractionMap;
        
        for (const auto &netNode : DN.networkNodes())
        {
            if (netNode.second.lock()->isBoundaryNode())
            {
                //                std::cout<<"NetNode is "<<netNode.second.lock()->tag()<<std::endl;
                
                std::set<std::shared_ptr<typename DislocationNetworkType::NetworkNodeType>> internalNodes; // used to be regular ptrs
                for (const auto &loopNode : netNode.second.lock()->loopNodes())
                {
                    internalNodes.insert(loopNode->periodicPrev()->networkNode);
                    internalNodes.insert(loopNode->periodicNext()->networkNode);
                }
                
                if (internalNodes.size() == 2)
                {
                    
                    //                    std::cout<<"internal nodes are "<<(*internalNodes.begin())->tag()<<","<<(*internalNodes.rbegin())->tag()<<std::endl;
                    
                    const NetNodePairType key((*internalNodes.begin())->sID < (*internalNodes.rbegin())->sID ? std::make_pair(*internalNodes.begin(), *internalNodes.rbegin()) : std::make_pair(*internalNodes.rbegin(), *internalNodes.begin())); // Line that was causing problem
                    
                    if (bndNodesContractionMap.find(key) == bndNodesContractionMap.end())
                    { // key not found
                        std::map<std::set<const PlanarMeshFace<dim> *>, std::set<typename DislocationNetworkType::NetworkNodeType *>> bndNodeMap;
                        
                        for (const auto &loopNode : (*internalNodes.begin())->loopNodes())
                        {
                            if (loopNode->periodicPrev()->networkNode == *internalNodes.rbegin())
                            {
                                for (const auto &bndNode : loopNode->boundaryPrev())
                                {
                                    //                                                        std::cout<<"adding bndNode "<<bndNode->networkNode->tag()<<std::endl;
                                    bndNodeMap[bndNode->networkNode->meshFaces()].insert(bndNode->networkNode.get());
                                }
                            }
                            else if (loopNode->periodicNext()->networkNode == *internalNodes.rbegin())
                            {
                                for (const auto &bndNode : loopNode->boundaryNext())
                                {
                                    //                                                        std::cout<<"adding bndNode "<<bndNode->networkNode->tag()<<std::endl;
                                    bndNodeMap[bndNode->networkNode->meshFaces()].insert(bndNode->networkNode.get());
                                }
                            }
                        }
                        
                        for (const auto &pair : bndNodeMap)
                        {
                            //                                                std::cout<<"Nodes of face"<<std::endl;
                            //                                                for(const auto& face : pair.first)
                            //                                                {
                            //                                                    std::cout<<face->outNormal().transpose()<<std::endl;
                            //                                                }
                            
                            //                                                std::cout<<"nodes size="<<pair.second.size()<<std::endl;
                            if (pair.second.size() > 1)
                            {
                                for (const auto &tempNode : pair.second)
                                {
                                    //                                                        std::cout<<"node size="<<tempNode->sID<<std::endl;
                                    if (tempNode != *pair.second.begin())
                                    {
                                        if ((tempNode->get_P() - (*pair.second.begin())->get_P()).norm() < FLT_EPSILON)
                                        {
                                            const NetNodePairType contractPair(tempNode->sID < (*pair.second.begin())->sID ? std::make_pair(tempNode, *pair.second.begin()) : std::make_pair(*pair.second.begin(), tempNode));
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
        
        for (const auto &pair : bndNodesContractionMap)
        {
            std::cout << "Contracting bnd nodes of intenal nodes " << pair.first.first->tag() << " " << pair.first.second->tag() << std::endl;
            for (const auto &bndPair : pair.second)
            {
                std::cout << "Contracting bnd nodes  " << bndPair.first->tag() << " " << bndPair.second->tag() << std::endl;
                DN.contract(bndPair.first, bndPair.second);
            }
        }
        std::cout << magentaColor << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]" << defaultColor << std::endl;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::contract0chordSegments()
    {
        std::cout<<"        Contracting zero-chord segments... "<<std::flush;
       const auto t0= std::chrono::system_clock::now();

       std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length

       //            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
       for (const auto& linkIter : DN.networkLinks())
       {
           const auto link(linkIter.second.lock());
           if(!linkIter.second.lock()->hasZeroBurgers())
           {
            //    if (!linkIter.second.lock()->source->isBoundaryNode() && !linkIter.second.lock()->sink->isBoundaryNode() && !linkIter.second.lock()->hasZeroBurgers())
               if (!link->source->isBoundaryNode() && !link->sink->isBoundaryNode() )
               {
//                   VectorDim chord(linkIter.second.lock()->chord()); // this is sink->get_P() - source->get_P()
//                   const double chordLength(chord.norm());

                //    if (chordLength <= FLT_EPSILON && (linkIter.second.lock()->source->get_V() - linkIter.second.lock()->sink->get_V()).squaredNorm() < FLT_EPSILON)
                //Removing the second condition is important (both the nodes may not have the same velocity due to complex network formation)
                //If the two nodes are not connected to any other links then this may be true otherwise it is not
                if (link->chordLength() <= FLT_EPSILON)
                { // toBeContracted part
                    toBeContracted.insert(std::make_pair(link->chordLength(), std::make_pair(link->source->sID, link->sink->sID)));
                }
                else
                {
                    const auto sourcePlanes(link->source->glidePlanes());
                    const auto   sinkPlanes(link->sink->glidePlanes());

                    if(sourcePlanes.size()==dim-1 && sourcePlanes==sinkPlanes)
                    {
                        if((link->chordLengthSquared()+link->chord().dot(link->sink->get_V() - link->source->get_V())* DN.simulationParameters.dt)<FLT_EPSILON)
                        {
                            toBeContracted.insert(std::make_pair(link->chordLength(), std::make_pair(link->source->sID, link->sink->sID)));
                        }
                    }
                    
                    
//                    const VectorDim velChange(linkIter.second.lock()->sink->get_V() - linkIter.second.lock()->source->get_V());
//                    const double velChangeNorm(velChange.norm());
//                    if (velChangeNorm > FLT_EPSILON)
//                    {
//                        const double dotPTemp((chord / chordLength).dot(velChange / velChangeNorm));
//                        // std::cout<<"For "<<linkIter.second.lock()->tag()<<"dotP Temp is "<<dotPTemp<<"=>"<<((velChangeNorm * DN.simulationParameters.dt) >= chordLength)<<std::endl;
//                        // if (fabs(dotPTemp + 1) < FLT_EPSILON)
//                        if (dotPTemp  > FLT_EPSILON)
//                        {
//                            if ((velChangeNorm * DN.simulationParameters.dt) >= chordLength)
//                            {
//                                toBeContracted.insert(std::make_pair(chordLength, std::make_pair(linkIter.second.lock()->source->sID, linkIter.second.lock()->sink->sID)));
//                            }
//                        }
//                    }
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
               
               std::cout<<"Contracting "<<Lij->tag()<<std::endl;
               Ncontracted+=DN.contract(Lij->source,Lij->sink);

               std::cout<<"Managing bnd nodes "<<std::endl;

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

                           SegmentSegmentDistance<dim - 1> ssd1(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                           // /*  Giacomo's Version
                           if (ssd1.dMin < FLT_EPSILON)
                           {
                               if (bndNode->periodicPlaneEdge.second)
                               {
                                   SegmentSegmentDistance<dim - 1> ssd2(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.second->source, *bndNode->periodicPlaneEdge.second->sink);
                                   const VectorLowerDim ssd1Pos(0.5 * (ssd1.x0 + ssd1.x1));
                                   const VectorLowerDim ssd2Pos(0.5 * (ssd2.x0 + ssd2.x1));
                               }
                               bndNode->set_P(VectorLowerDim(0.5 * (ssd1.x0 + ssd1.x1)));
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
    
    /**************************************************************************/
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
                            if (!loopLink->source->periodicPlaneEdge.first)
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
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::minMeshSize(const SimplicialMesh<dim> &mesh)
    {
        return std::min(mesh.xMax(0) - mesh.xMin(0), std::min(mesh.xMax(1) - mesh.xMin(1), mesh.xMax(2) - mesh.xMin(2)));
    }
    
    template class DislocationNetworkRemesh<DislocationNetwork<3,0>>;

    
} // namespace model
#endif

