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
            std::cout<<"CrossSlip: finding branches "<<std::flush;
            const auto t0= std::chrono::system_clock::now();

//            crossSlipDeq.clear();
//            for(const auto& linkIter : DN.networkLinks())
//            {
//                const auto link(linkIter.second.lock());
//                if(   !link->isBoundarySegment()
//                   && !link->source->isBoundaryNode()
//                   && !link->sink->isBoundaryNode()
//                   && !link->isGrainBoundarySegment()
//                   && !link->source->isGrainBoundaryNode()
//                   && !link->sink->isGrainBoundaryNode()
//                   && !link->hasZeroBurgers()
//                   && link->isGlissile()
//                   && link->chord().normalized().cross(link->burgers().normalized()).norm()<=DN.crossSlipModel->sinCrossSlip
//                   && link->chord().norm()>2.0*DN.networkRemesher.Lmin
//                   )
//                {
//                    DN.crossSlipModel->addToCrossSlip(*link,crossSlipDeq);
//                }
//            }
        
//            std::cout<<"Collecting branches"<<std::endl;
//            std::deque<std::deque<std::pair<const LoopLinkType*,int>>> csBranches;
            csNodes.clear();
            for(const auto& loopIter : DN.loops())
            {
                const auto loop(loopIter.second.lock());
                loop->crossSlipBranches(csNodes);
            }
            std::cout<<csNodes.size()<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

        
        }
    }

//    /**********************************************************************/
//    template <typename DislocationNetworkType>
//    void  DislocationCrossSlip<DislocationNetworkType>::executeKernel()
//    {
//
//    }


    template <typename DislocationNetworkType>
    void  DislocationCrossSlip<DislocationNetworkType>::execute()
    {
        
        if(DN.crossSlipModel)
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"CrossSlip: aligning branches"<<std::flush;
//            std::cout<<"csNodes.size()="<<csNodes.size()<<std::endl;
//            std::cout<<"aligning nodes"<<std::endl;
            
            const double dcsMin(20.0);
            
            size_t executed(0);
            for(const auto& branch : csNodes)
            {// Align all nodes in a brach on the same cross-slip plane
                if(branch.first.size())
                {
//                    std::cout<<"branch: "<<branch.second<<std::endl;
                    const auto& loop(branch.first.back()->loop());
                    if(loop->glidePlane)
                    {
                        
                        // Compute the average position of the brach nodes in the superspace
                        VectorDim midPoint(VectorDim::Zero());
                        for(const auto& node : branch.first)
                        {
                            midPoint+=node->get_P();
//                            std::cout<<"    "<<node->tag()<<std::endl;
                        }
                        midPoint/=branch.first.size();
                        
                        const auto& grain(loop->grain);
                        const auto& crosSlipSystem(grain.singleCrystal->slipSystems()[branch.second]);
                        // Compute the heigth of the crossSlip plane closest to midPoint
                        
                        const int dhMin(dcsMin/crosSlipSystem->n.planeSpacing());
                        const double height0(crosSlipSystem->n.closestPlaneIndexOfPoint(midPoint));
                        const long int height(std::round(height0/dhMin)*dhMin);
                        std::cout<<"dhMin="<<dhMin<<", height="<<height0<<" -> "<<height<<std::endl;
                        
                        const VectorDim planePoint(height*crosSlipSystem->n.planeSpacing()*crosSlipSystem->unitNormal);
                        PlanePlaneIntersection<dim> ppi(loop->glidePlane->P,
                                                        loop->glidePlane->unitNormal,
                                                        planePoint,
                                                        crosSlipSystem->unitNormal);

                        for(auto& node : branch.first)
                        {// align all nodes in the branch to perfect screw direction
                            const VectorDim oldP(node->networkNode->get_P());
                            const VectorDim deltaP(ppi.P-node->get_P()+(node->get_P()-ppi.P).dot(ppi.d)*ppi.d);
//                            std::cout<<"node "<<node->tag()<<", P="<<node->networkNode->get_P().transpose()<<std::endl;
//                            std::cout<<"node "<<node->tag()<<", dP"<<deltaP.transpose()<<std::endl;
                            node->networkNode->trySet_P(node->networkNode->get_P()+deltaP);
                            if((node->networkNode->get_P()-oldP-deltaP).norm()>FLT_EPSILON)
                            {
                                std::cout<<"node "<<node->tag()<<", cannot move by"<<deltaP.transpose()<<std::endl;
                            }
//
                        }
                    }
                }
            }
            DN.updateBoundaryNodes(); // After aligning, bnd nodes must be updated
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

            const auto t1= std::chrono::system_clock::now();
            std::cout<<"CrossSlip: creating loops "<<std::flush;
            std::vector<std::shared_ptr<LoopType>> newLoops;
            for(const auto& branch : csNodes)
            {
                if(branch.first.size())
                {
//                    std::cout<<"branch: "<<branch.second<<std::endl;
                    const auto& loop(branch.first.back()->loop());
                    if(loop->glidePlane)
                    {
                        const auto& grain(loop->grain);
                        const auto& crosSlipSystem(grain.singleCrystal->slipSystems()[branch.second]); //THERE COULD BE A PROBLEM HERE WITH SIGN OF B.
//                        std::cout<<"Point="<<branch.first.back()->get_P().transpose()<<std::endl;
//                        std::cout<<"normal="<<crosSlipSystem->unitNormal.transpose()<<std::endl;

                        // Define the cross-slip plane through last NetworkNode in branch
                        GlidePlaneKey<dim> crossSlipPlaneKey(branch.first.back()->networkNode->get_P(), crosSlipSystem->n);
                        const auto crossSlipPlane(DN.glidePlaneFactory.getFromKey(crossSlipPlaneKey));
                        auto crossSlipLoop(DN.loops().create(crosSlipSystem->s.cartesian(), crossSlipPlane));
                        const auto periodicCrossSlipPlane(DN.periodicGlidePlaneFactory->get(crossSlipPlane->key));
//                        std::cout<<"Constructed crossSlipPlane="<<std::endl;
                        
                        std::vector<std::shared_ptr<LoopNodeType>> crossSlipLoopNodes;
                        LoopNodeType* currentLoopNode(branch.first.back().get());
//                        std::cout<<"currentLoopNode="<<currentLoopNode->get_P().transpose()<<std::endl;
//                        std::cout<<"Computing shift"<<std::endl;

//                        VectorDim shift(periodicCrossSlipPlane->findPatch(periodicCrossSlipPlane->referencePlane->localPosition(currentLoopNode->get_P()),VectorDim::Zero()));
                        VectorDim shift(VectorDim::Zero()); // shift of first crossSlipLoopNode is zero, since crossSlipPlane is defined throught this point

//                        std::cout<<"Getting patch"<<std::endl;

                        std::shared_ptr<PeriodicPlanePatch<dim>> crossSlipPatch(periodicCrossSlipPlane->getPatch(shift));
                        
                        // add a new cross-slip node at the end of the branch
                        std::shared_ptr<NetworkNodeType> newNetNode(DN.networkNodes().create(currentLoopNode->networkNode->get_P(), VectorDim::Zero(), 1.0));
                        
                        bool allContained(crossSlipPlane->contains(newNetNode->get_P()));
                        std::cout<<"crossSlipPlane->contains("<<newNetNode->tag()<<")? "<<crossSlipPlane->contains(newNetNode->get_P())<<std::endl;
                        for(auto& node : branch.first)
                        {// align all nodes in the branch to perfect screw direction
                            allContained=(allContained && crossSlipPlane->contains(node->networkNode->get_P()));
                            std::cout<<"crossSlipPlane->contains("<<node->networkNode->tag()<<")? "<<crossSlipPlane->contains(node->networkNode->get_P())<<std::endl;
                        }
                        
                        if(allContained)
                        {
                            crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, newNetNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(nullptr,nullptr)));
                            
                            int crossID(0);
                            while(true)
                            {// kepping adding new cross-slip nodes traversing the branch in reverse order
    //                            std::cout<<"While loop"<<std::endl;
                                
//                                std::cout<<"branchNode "<<currentLoopNode->networkNode->sID<<std::endl;

                                std::shared_ptr<PeriodicPlaneEdge<dim>> edge1(nullptr);
                                std::shared_ptr<PeriodicPlaneEdge<dim>> edge2(nullptr);
                                if(currentLoopNode->periodicPlaneEdge.first)
                                {// a boundary node on one edge
                                    edge1=crossSlipPatch->edgeOnFaces(currentLoopNode->periodicPlaneEdge.first->meshIntersection->faces);
                                    if(!edge1)
                                    {
                                        throw std::runtime_error("crossSlipPatch patch does not have edge1.");
                                    }
                                    if(currentLoopNode->periodicPlaneEdge.second)
                                    {// a boundary node on two edges
                                        edge2=crossSlipPatch->edgeOnFaces(currentLoopNode->periodicPlaneEdge.second->meshIntersection->faces);
                                        if(!edge2)
                                        {
                                            throw std::runtime_error("crossSlipPatch patch does not have edge2.");
                                        }
                                    }
                                }
                                else
                                {// not a boundary node
                                    edge1=nullptr;
                                    edge2=nullptr;
                                }
                                
    //                            std::cout<<"          currentLoopNode->get_P()="<<currentLoopNode->get_P().transpose()<<std::endl;
    //                            std::cout<<"currentLoopNode->networkNode->get_P()="<<currentLoopNode->networkNode->get_P().transpose()<<std::endl;
    //                            std::cout<<"shift="<<shift.transpose()<<std::endl;

    //                            crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
                                crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
    //                            std::cout<<"crossSlipLoopNodes.back()->get_P()="<<crossSlipLoopNodes.back()->get_P().transpose()<<std::endl;

                                if(edge1)
                                {
                                    if(crossID%2 == 0)
                                    {
                                        shift+=edge1->deltaShift;
                                        if(edge2)
                                        {
                                            shift+=edge2->deltaShift;
                                        }
                                        crossSlipPatch=periodicCrossSlipPlane->getPatch(shift);
                                    }
                                    crossID++;
                                }
                                
                                
                                if(currentLoopNode==branch.first.front().get())
                                {
                                    // add a new cross-slip node at the beginning of the branch and break
                                    newNetNode= (currentLoopNode->networkNode==branch.first.back()->networkNode? crossSlipLoopNodes.front()->networkNode : DN.networkNodes().create(currentLoopNode->networkNode->get_P(), VectorDim::Zero(), 1.0));
                                    crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, newNetNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(nullptr,nullptr)));
                                    break;
                                }
                                else
                                {
                                    currentLoopNode=currentLoopNode->prev.first;
                                }
                            }
                            DN.insertLoop(crossSlipLoop, crossSlipLoopNodes);
                            newLoops.push_back(crossSlipLoop);
                            executed++;
                        }
                        else
                        {
                            std::cout<<" not all nodes contained in crossSlipPlane "<<std::endl;
                        }
                    }
                }
            }
            DN.updateBoundaryNodes(); // After creating loops, bnd nodes must be updated
            
            for(auto& loop :  newLoops)
            {
                for(const auto& link : loop->loopLinks())
                {
                    if(link->networkLink())
                    {
                        link->networkLink()->updateQuadraturePointsSeg();
                    }
                }
                loop->computeStackingFaultForces();
                for(const auto& link : loop->loopLinks())
                {
                    if(link->networkLink())
                    {
                        link->networkLink()->assembleGlide(true);
                    }
                }

            }
            
            std::cout<<executed<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;


//            std::cout<<executed<< "executed "<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
        
        
        
        
        
//        if(DN.crossSlipModel && false)
//        {
//            const auto t0= std::chrono::system_clock::now();
//            std::cout<<"Cross slip: "<<std::flush;
//            size_t executed(0);
//            for(const auto& tup : crossSlipDeq)
//            {
//                const std::shared_ptr<NetworkNodeType>& source(std::get<0>(tup));
//                const std::shared_ptr<NetworkNodeType>& sink(std::get<1>(tup));
//                const size_t& sourceID(source->sID);
//                const size_t& sinkID(sink->sID);
//                const size_t& grainID(std::get<2>(tup));
//                const size_t& slipID(std::get<3>(tup));
//                
//                const std::shared_ptr<NetworkNodeType> isSource(DN.networkNodes().get(sourceID));
//                const std::shared_ptr<NetworkNodeType> isSink(DN.networkNodes().get(sinkID));
//                const auto isLink(DN.networkLinks().get(std::make_pair(sourceID,sinkID)));
//                
//                const auto& crosSlipSystem(DN.poly.grain(grainID).singleCrystal->slipSystems()[slipID]); // last element in map has highest pkGlide
//                
//                if(isSource && isSink && isLink)
//                {
//                    
//                    // Align source and sink to perfect screw orientation
//                    const VectorDim midPoint(0.5*(isSource->get_P()+isSink->get_P()));
//                    const long int height(crosSlipSystem->n.closestPlaneIndexOfPoint(midPoint));
//                    //
//                    //                        const int height=LatticePlane::computeHeight(crosSlipSystem->n,midPoint).second;
//                    const VectorDim planePoint(height*crosSlipSystem->n.planeSpacing()*crosSlipSystem->unitNormal);
//                    
//                    //const VectorDim planePoint2=midPoint-(midPoint-planePoint).dot(crosSlipSystem->unitNormal)*crosSlipSystem->unitNormal; // closest point to midPoint on the crossSlip plane
//                    
//                    //                        PlanePlaneIntersection<dim> ppi(midPoint,isLink.second->glidePlaneNormal(),
//                    //                                                        planePoint2,crosSlipSystem->unitNormal);
//                    
//                    
//                    
//                    
//                    PlanePlaneIntersection<dim> ppi((*isLink->loopLinks().begin())->loop->glidePlane->P,
//                                                    (*isLink->loopLinks().begin())->loop->glidePlane->unitNormal,
//                                                    planePoint,
//                                                    crosSlipSystem->unitNormal);
//                    
//                    
//                    const VectorDim newSourceP(ppi.P+(isSource->get_P()-ppi.P).dot(ppi.d)*ppi.d);
//                    const VectorDim newSinkP(ppi.P+(isSink->get_P()-ppi.P).dot(ppi.d)*ppi.d);
//                    
//                    //Only one loopNode
//                    //New position do not belong to the patch boundary based on the glide plane intersection
//                    bool sourcePositiononBoundary(false);
//                    bool sinkPositiononBoundary(false);
//                    
//                    for (const auto& gp : isSource->glidePlanes())
//                    {
//                        if (gp->meshIntersections.contains(newSourceP))
//                        {
//                            sourcePositiononBoundary=true;
//                            break;
//                        }
//                    }
//                    
//                    for (const auto& gp : isSink->glidePlanes())
//                    {
//                        if (gp->meshIntersections.contains(newSourceP))
//                        {
//                            sinkPositiononBoundary=true;
//                            break;
//                        }
//                    }
//                    
//                    if (!sourcePositiononBoundary && !sinkPositiononBoundary)
//                    {
//                        if (isSource->isMovableTo(newSourceP) && isSink->isMovableTo(newSinkP))
//                        {
//                            
//                            VerboseCrossSlip(1, "cross-slip " << sourceID << "->" << sinkID << std::endl;);
//                            
//                            // Re-align source and sink
//                            isSource->trySet_P(newSourceP);
//                            isSink->trySet_P(newSinkP);
//                            
//                            DN.updateBoundaryNodes();
//                            if ((isSource->get_P() - newSourceP).norm() < FLT_EPSILON && (isSink->get_P() - newSinkP).norm() < FLT_EPSILON)
//                            {
//                                
//                                // Check if source and sink are already part of loops on the conjugate plane
//                                
//                                // Construct and insert new loop in conjugate plane
//                                const VectorDim newNodeP(0.5 * (isSource->get_P() + isSink->get_P()));
//                                //                                const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
//                                std::shared_ptr<NetworkNodeType> newNode(DN.networkNodes().create(newNodeP, VectorDim::Zero(), 1.0));
//                                
//                                //                                std::vector<size_t> nodeIDs;
//                                std::vector<std::shared_ptr<NetworkNodeType>> networkNodes;
//                                
//                                networkNodes.push_back(isSink);   // insert in reverse order, sink first, source second
//                                networkNodes.push_back(isSource); // insert in reverse order, sink first, source second
//                                networkNodes.push_back(newNode);
//                                
//                                GlidePlaneKey<dim> loopPlaneKey(newNodeP, DN.poly.grain(grainID).singleCrystal->slipSystems()[slipID]->n);
//                                const auto glidePlane(DN.glidePlaneFactory.getFromKey(loopPlaneKey));
//                                auto glissileLoop(DN.loops().create(DN.poly.grain(grainID).singleCrystal->slipSystems()[slipID]->s.cartesian(), glidePlane));
//                                
//                                std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
//                                
//                                const auto periodicGlidePlane(DN.periodicGlidePlaneFactory->get(glidePlane->key));
//                                const auto periodicPatch(periodicGlidePlane->getPatch(VectorDim::Zero()));
//                                
//                                if (isSink->isBoundaryNode())
//                                {
//                                    assert(!isSource->isBoundaryNode() && "Cross-slip cannot happen at the boundary");
//                                    //Get the loopNodes of Sink
//                                    std::set<short int> edgeIDs;
//                                    for (const auto &edge : periodicPatch->edges())
//                                    {
//                                        if (((isSink->get_P() - edge->meshIntersection->P0).cross(isSink->get_P() - edge->meshIntersection->P1)).squaredNorm() < FLT_EPSILON)
//                                        {
//                                            edgeIDs.insert(edge->edgeID);
//                                        }
//                                    }
//                                    assert(edgeIDs.size() == 1 && "Cross-Slip  at corner");
//                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSink, isSink->get_P(), periodicPatch, std::make_pair(periodicPatch->edges()[*edgeIDs.begin()],nullptr)));
//                                }
//                                else
//                                {
//                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSink, isSink->get_P(), periodicPatch, std::make_pair(nullptr,nullptr)));
//                                }
//                                
//                                if (isSource->isBoundaryNode())
//                                {
//                                    assert(!isSink->isBoundaryNode() && "Cross-slip cannot happen at the boundary");
//                                    //Get the loopNodes of Sink
//                                    std::set<short int> edgeIDs;
//                                    for (const auto &edge : periodicPatch->edges())
//                                    {
//                                        if (((isSource->get_P() - edge->meshIntersection->P0).cross(isSource->get_P() - edge->meshIntersection->P1)).squaredNorm() < FLT_EPSILON)
//                                        {
//                                            edgeIDs.insert(edge->edgeID);
//                                        }
//                                    }
//                                    assert(edgeIDs.size() == 1 && "Glissile Junction Intersection at corner");
//                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSource, isSource->get_P(), periodicPatch, std::make_pair(periodicPatch->edges()[*edgeIDs.begin()],nullptr)));
//                                }
//                                else
//                                {
//                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSource, isSource->get_P(), periodicPatch, std::make_pair(nullptr,nullptr)));
//                                }
//                                
//                                //New node cannot be a boundary node
//                                loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, newNode, newNode->get_P(), periodicPatch, std::make_pair(nullptr,nullptr)));
//                                //                                nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
//                                //                                nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
//                                //                                nodeIDs.push_back(newNodeID);
//                                
//                                //                                LatticePlane loopPlane(newNodeP,DN.poly.grain(grainID).slipSystems()[slipID]->n);
//                                //                                GlidePlaneKey<dim> loopPlaneKey(grainID,loopPlane);
//                                
//                                //                                DN.insertLoop(nodeIDs,
//                                //                                              DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                                //                                              DN.glidePlaneFactory.get(loopPlaneKey));
//                                
//                                DN.insertLoop(glissileLoop, loopNodes);
//                                
//                                executed++;
//                            }
//                        }
//                    }
//                }
//            }
//            std::cout<<executed<< "executed "<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//        }
    }

    template class DislocationCrossSlip<DislocationNetwork<3,0>>;
}
#endif
