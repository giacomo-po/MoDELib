/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_cpp_
#define model_DislocationCrossSlip_cpp_

#include <numbers>
#include <memory>

#include <DislocationCrossSlip.h>
#include <CrossSlipModels.h>



namespace model
{

    template <typename DislocationNetworkType>
    std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>> DislocationCrossSlip<DislocationNetworkType>::getModel(const PolycrystallineMaterialBase& material,const DDtraitsIO& traitsIO)
    {
        const int crossSlipModel(TextFileParser(traitsIO.ddFile).readScalar<int>("crossSlipModel",true));
        switch (crossSlipModel)
        {
            case 0:
            {// no cross-slip
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                break;
            }
                
            case 1:
            {// deterministic cross-slip model based on max glide PK force
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(new AthermalCrossSlipModel<DislocationNetworkType>(traitsIO));
                break;
            }
                
            case 2:
            {
                if(material.crystalStructure=="HEX")
                {
                    return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(new EscaigCrossSlipModelHEX<DislocationNetworkType>(material,traitsIO));
                }
                else
                {
                    throw std::runtime_error("Unknown cross-slip model "+std::to_string(crossSlipModel)+" for crystal structure "+material.crystalStructure);
                    return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                }
                break;
            }
                
            default:
            {
                throw std::runtime_error("Unknown cross-slip model "+std::to_string(crossSlipModel));
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                break;
            }
        }
    }

    /**********************************************************************/
    template <typename DislocationNetworkType>
    DislocationCrossSlip<DislocationNetworkType>::DislocationCrossSlip(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    /* init */,verboseCrossSlip(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseCrossSlip",true))
    {
    }

    /**********************************************************************/
    template <typename DislocationNetworkType>
    void  DislocationCrossSlip<DislocationNetworkType>::findCrossSlipSegments()
    {
        
        if(DN.crossSlipModel)
        {
            std::cout<<"CrossSlip: "<<std::flush;
            const auto t0= std::chrono::system_clock::now();

            crossSlipDeq.clear();
            for(const auto& linkIter : DN.networkLinks())
            {
                const auto link(linkIter.second.lock());
                if(   !link->isBoundarySegment()
                   && !link->source->isBoundaryNode()
                   && !link->sink->isBoundaryNode()
                   && !link->isGrainBoundarySegment()
                   && !link->source->isGrainBoundaryNode()
                   && !link->sink->isGrainBoundaryNode()
                   && !link->hasZeroBurgers()
                   && link->isGlissile()
                   && link->chord().normalized().cross(link->burgers().normalized()).norm()<=DN.crossSlipModel->sinCrossSlip
                   && link->chord().norm()>2.0*DN.networkRemesher.Lmin
                   )
                {
                    DN.crossSlipModel->addToCrossSlip(*link,crossSlipDeq);
                }
            }
            std::cout<<crossSlipDeq.size()<<" found "<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
            std::cout<<"Collecting branches"<<std::endl;
//            std::deque<std::deque<std::pair<const LoopLinkType*,int>>> csBranches;
            csNodes.clear();
            for(const auto& loopIter : DN.loops())
            {
                const auto loop(loopIter.second.lock());
                loop->crossSlipBranches(csNodes);
            }
            
        
        }
    }

    /**********************************************************************/
    template <typename DislocationNetworkType>
    void  DislocationCrossSlip<DislocationNetworkType>::execute()
    {
        
        if(DN.crossSlipModel && false)
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Cross slip: "<<std::flush;
            std::cout<<"csNodes.size()="<<csNodes.size()<<std::endl;
            size_t executed(0);
            for(const auto& branch : csNodes)
            {// Align all nodes in a brach on the same cross-slip plane
                if(branch.first.size())
                {
                    std::cout<<"branch: "<<branch.second<<std::endl;
                    const auto& loop(branch.first.back()->loop());
                    if(loop->glidePlane)
                    {
                        VectorDim midPoint(VectorDim::Zero());
                        for(const auto& node : branch.first)
                        {
                            midPoint+=node->get_P();
                            std::cout<<"    "<<node->tag()<<std::endl;
                        }
                        midPoint/=branch.first.size();

                        const auto& grain(loop->grain);
                        const auto& crosSlipSystem(grain.slipSystems()[branch.second]);
                        const long int height(crosSlipSystem->n.closestPlaneIndexOfPoint(midPoint));
                        const VectorDim planePoint(height*crosSlipSystem->n.planeSpacing()*crosSlipSystem->unitNormal);
                        PlanePlaneIntersection<dim> ppi(loop->glidePlane->P,
                                                        loop->glidePlane->unitNormal,
                                                        planePoint,
                                                        crosSlipSystem->unitNormal);

                        for(auto& node : branch.first)
                        {// align all nodes in the branch to perfect screw direction
                            const VectorDim deltaP(ppi.P-node->get_P()+(node->get_P()-ppi.P).dot(ppi.d)*ppi.d);
                            node->networkNode->trySet_P(node->networkNode->get_P()+deltaP);
                        }
                    }
                }
            }
            DN.updateBoundaryNodes(); // After aligning bnd nodes must be updated
            
            for(const auto& branch : csNodes)
            {
                if(branch.first.size())
                {
                    std::cout<<"branch: "<<branch.second<<std::endl;
                    const auto& loop(branch.first.back()->loop());
                    if(loop->glidePlane)
                    {
                        const auto& grain(loop->grain);
                        const auto& crosSlipSystem(grain.slipSystems()[branch.second]);
                        GlidePlaneKey<dim> crossSlipPlaneKey(branch.first.back()->get_P(), crosSlipSystem->n);
                        const auto crossSlipPlane(DN.glidePlaneFactory.getFromKey(crossSlipPlaneKey));
                        auto crossSlipLoop(DN.loops().create(crosSlipSystem->s.cartesian(), crossSlipPlane));
                        const auto periodicCrossSlipPlane(DN.periodicGlidePlaneFactory->get(crossSlipPlane->key));
                        
                        VectorDim shift(VectorDim::Zero());
                        std::vector<std::shared_ptr<LoopNodeType>> crossSlipLoopNodes;
                        LoopNodeType* currentLoopNode(branch.first.back().get());
                        while(true)
                        {
                            
//                            if()
                            
                            if(currentLoopNode==branch.first.front().get())
                            {
                                break;
                            }
                            else
                            {
                                currentLoopNode=currentLoopNode->prev.first;
                            }
                        }
                        
//                        for(typename std::deque<std::shared_ptr<LoopNodeType>>::const_reverse_iterator iter= branch.first.rbegin();
//                            iter!=branch.first.rend();++iter)
//                        {// align all nodes in the branch to perfect screw direction
//                            shift=periodicCrossSlipPlane->findPatch(periodicCrossSlipPlane->referencePlane->localPosition((*iter)->get_P()),shift);
//                            crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, (*iter)->networkNode, (*iter)->networkNode->get_P()-shift, periodicCrossSlipPlane->getPatch(shift), std::make_pair(nullptr,nullptr)));
//                        }
//                        for(typename std::deque<std::shared_ptr<LoopNodeType>>::const_iterator iter= branch.first.begin();
//                            iter!=branch.first.end();++iter)
//                        {// align all nodes in the branch to perfect screw direction
//                            shift=periodicCrossSlipPlane->findPatch(periodicCrossSlipPlane->referencePlane->localPosition((*iter)->get_P()),shift);
//                            std::shared_ptr<NetworkNodeType> newNode(DN.networkNodes().create((*iter)->networkNode->get_P(), VectorDim::Zero(), 1.0));
//                            crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, newNode, newNode->get_P()-shift, periodicCrossSlipPlane->getPatch(shift), std::make_pair(nullptr,nullptr)));
//                        }
//                        DN.insertLoop(crossSlipLoop, crossSlipLoopNodes);
                    }
                }
            }

            std::cout<<executed<< "executed "<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
        
        
        
        
        
        if(DN.crossSlipModel)
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Cross slip: "<<std::flush;
            size_t executed(0);
            for(const auto& tup : crossSlipDeq)
            {
                const std::shared_ptr<NetworkNodeType>& source(std::get<0>(tup));
                const std::shared_ptr<NetworkNodeType>& sink(std::get<1>(tup));
                const size_t& sourceID(source->sID);
                const size_t& sinkID(sink->sID);
                const size_t& grainID(std::get<2>(tup));
                const size_t& slipID(std::get<3>(tup));
                
                const std::shared_ptr<NetworkNodeType> isSource(DN.networkNodes().get(sourceID));
                const std::shared_ptr<NetworkNodeType> isSink(DN.networkNodes().get(sinkID));
                const auto isLink(DN.networkLinks().get(std::make_pair(sourceID,sinkID)));
                
                const auto& crosSlipSystem(DN.poly.grain(grainID).slipSystems()[slipID]); // last element in map has highest pkGlide
                
                if(isSource && isSink && isLink)
                {
                    
                    // Align source and sink to perfect screw orientation
                    const VectorDim midPoint(0.5*(isSource->get_P()+isSink->get_P()));
                    const long int height(crosSlipSystem->n.closestPlaneIndexOfPoint(midPoint));
                    //
                    //                        const int height=LatticePlane::computeHeight(crosSlipSystem->n,midPoint).second;
                    const VectorDim planePoint(height*crosSlipSystem->n.planeSpacing()*crosSlipSystem->unitNormal);
                    
                    //const VectorDim planePoint2=midPoint-(midPoint-planePoint).dot(crosSlipSystem->unitNormal)*crosSlipSystem->unitNormal; // closest point to midPoint on the crossSlip plane
                    
                    //                        PlanePlaneIntersection<dim> ppi(midPoint,isLink.second->glidePlaneNormal(),
                    //                                                        planePoint2,crosSlipSystem->unitNormal);
                    
                    
                    
                    
                    PlanePlaneIntersection<dim> ppi((*isLink->loopLinks().begin())->loop->glidePlane->P,
                                                    (*isLink->loopLinks().begin())->loop->glidePlane->unitNormal,
                                                    planePoint,
                                                    crosSlipSystem->unitNormal);
                    
                    
                    const VectorDim newSourceP(ppi.P+(isSource->get_P()-ppi.P).dot(ppi.d)*ppi.d);
                    const VectorDim newSinkP(ppi.P+(isSink->get_P()-ppi.P).dot(ppi.d)*ppi.d);
                    
                    //Only one loopNode
                    //New position do not belong to the patch boundary based on the glide plane intersection
                    bool sourcePositiononBoundary(false);
                    bool sinkPositiononBoundary(false);
                    
                    for (const auto& gp : isSource->glidePlanes())
                    {
                        if (gp->meshIntersections.contains(newSourceP))
                        {
                            sourcePositiononBoundary=true;
                            break;
                        }
                    }
                    
                    for (const auto& gp : isSink->glidePlanes())
                    {
                        if (gp->meshIntersections.contains(newSourceP))
                        {
                            sinkPositiononBoundary=true;
                            break;
                        }
                    }
                    
                    if (!sourcePositiononBoundary && !sinkPositiononBoundary)
                    {
                        if (isSource->isMovableTo(newSourceP) && isSink->isMovableTo(newSinkP))
                        {
                            
                            VerboseCrossSlip(1, "cross-slip " << sourceID << "->" << sinkID << std::endl;);
                            
                            // Re-align source and sink
                            isSource->trySet_P(newSourceP);
                            isSink->trySet_P(newSinkP);
                            
                            DN.updateBoundaryNodes();
                            if ((isSource->get_P() - newSourceP).norm() < FLT_EPSILON && (isSink->get_P() - newSinkP).norm() < FLT_EPSILON)
                            {
                                
                                // Check if source and sink are already part of loops on the conjugate plane
                                
                                // Construct and insert new loop in conjugate plane
                                const VectorDim newNodeP(0.5 * (isSource->get_P() + isSink->get_P()));
                                //                                const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
                                std::shared_ptr<NetworkNodeType> newNode(DN.networkNodes().create(newNodeP, VectorDim::Zero(), 1.0));
                                
                                //                                std::vector<size_t> nodeIDs;
                                std::vector<std::shared_ptr<NetworkNodeType>> networkNodes;
                                
                                networkNodes.push_back(isSink);   // insert in reverse order, sink first, source second
                                networkNodes.push_back(isSource); // insert in reverse order, sink first, source second
                                networkNodes.push_back(newNode);
                                
                                GlidePlaneKey<dim> loopPlaneKey(newNodeP, DN.poly.grain(grainID).slipSystems()[slipID]->n);
                                const auto glidePlane(DN.glidePlaneFactory.getFromKey(loopPlaneKey));
                                auto glissileLoop(DN.loops().create(DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(), glidePlane));
                                
                                std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
                                
                                const auto periodicGlidePlane(DN.periodicGlidePlaneFactory->get(glidePlane->key));
                                const auto periodicPatch(periodicGlidePlane->getPatch(VectorDim::Zero()));
                                
                                if (isSink->isBoundaryNode())
                                {
                                    assert(!isSource->isBoundaryNode() && "Cross-slip cannot happen at the boundary");
                                    //Get the loopNodes of Sink
                                    std::set<short int> edgeIDs;
                                    for (const auto &edge : periodicPatch->edges())
                                    {
                                        if (((isSink->get_P() - edge->meshIntersection->P0).cross(isSink->get_P() - edge->meshIntersection->P1)).squaredNorm() < FLT_EPSILON)
                                        {
                                            edgeIDs.insert(edge->edgeID);
                                        }
                                    }
                                    assert(edgeIDs.size() == 1 && "Cross-Slip  at corner");
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSink, isSink->get_P(), periodicPatch, std::make_pair(periodicPatch->edges()[*edgeIDs.begin()],nullptr)));
                                }
                                else
                                {
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSink, isSink->get_P(), periodicPatch, std::make_pair(nullptr,nullptr)));
                                }
                                
                                if (isSource->isBoundaryNode())
                                {
                                    assert(!isSink->isBoundaryNode() && "Cross-slip cannot happen at the boundary");
                                    //Get the loopNodes of Sink
                                    std::set<short int> edgeIDs;
                                    for (const auto &edge : periodicPatch->edges())
                                    {
                                        if (((isSource->get_P() - edge->meshIntersection->P0).cross(isSource->get_P() - edge->meshIntersection->P1)).squaredNorm() < FLT_EPSILON)
                                        {
                                            edgeIDs.insert(edge->edgeID);
                                        }
                                    }
                                    assert(edgeIDs.size() == 1 && "Glissile Junction Intersection at corner");
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSource, isSource->get_P(), periodicPatch, std::make_pair(periodicPatch->edges()[*edgeIDs.begin()],nullptr)));
                                }
                                else
                                {
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSource, isSource->get_P(), periodicPatch, std::make_pair(nullptr,nullptr)));
                                }
                                
                                //New node cannot be a boundary node
                                loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, newNode, newNode->get_P(), periodicPatch, std::make_pair(nullptr,nullptr)));
                                //                                nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
                                //                                nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
                                //                                nodeIDs.push_back(newNodeID);
                                
                                //                                LatticePlane loopPlane(newNodeP,DN.poly.grain(grainID).slipSystems()[slipID]->n);
                                //                                GlidePlaneKey<dim> loopPlaneKey(grainID,loopPlane);
                                
                                //                                DN.insertLoop(nodeIDs,
                                //                                              DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
                                //                                              DN.glidePlaneFactory.get(loopPlaneKey));
                                
                                DN.insertLoop(glissileLoop, loopNodes);
                                
                                executed++;
                            }
                        }
                    }
                }
            }
            std::cout<<executed<< "executed "<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
    }

    template class DislocationCrossSlip<DislocationNetwork<3,0>>;
}
#endif
