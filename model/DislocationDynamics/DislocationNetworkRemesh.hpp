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
    
    template <typename DislocationNetworkType>
    double DislocationNetworkRemesh<DislocationNetworkType>::minMeshSize(const SimplicialMesh<dim>& mesh)
    {
        return std::min(mesh.xMax(0)-mesh.xMin(0),std::min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
    }
    
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
//                contract0chordSegments();
            }
        }
    }
    
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remeshByRemoval()
    {
        const auto t0= std::chrono::system_clock::now();
        model::cout<<"        Remeshing network: removing... "<<std::flush;
        DN.danglingBoundaryLoopNodes.clear();
        
        std::set<size_t> bndLoopNodes;
        
        std::deque<size_t> toBeRemoved;
        for(const auto& node : DN.networkNodes())
        {
            //                std::cout<<"node "<<node.second->sID<<" "<<node.second->isSimpleBoundaryNode()<<" "<<node.second->isSimpleGrainBoundaryNode()<<std::endl;
            if(!node.second.lock()->masterNode)
            {
                if(node.second.lock()->isRemovable(Lmin,relativeAreaThreshold))
                {
                    toBeRemoved.push_back(node.second.lock()->sID); // insert image IDs BEFORE master ID, so that images are removed FIRST
                
                    for(const auto& loopNode : node.second.lock()->loopNodes())
                    {
                        for(const auto& bndNode : loopNode->boundaryPrev())
                        {
                            bndLoopNodes.insert(bndNode->sID);
                        }
                        for(const auto& bndNode : loopNode->boundaryNext())
                        {
                            bndLoopNodes.insert(bndNode->sID);
                        }
                    }
                }
            }
        }
//
        size_t Nremoved=0;
        for(const auto& nodeID : toBeRemoved)
        {
            DN.removeNetworkNode(nodeID);
            Nremoved++;
        }
        model::cout<<" ("<<Nremoved<<" removed). Updating boundary nodes "<<std::flush;

        for(const size_t& nodeID : bndLoopNodes)
        {
            auto nodeIter(DN.loopNodes().find(nodeID));
            if(nodeIter!=DN.loopNodes().end())
            {
                auto bndNode(nodeIter->second.lock());
                
                const auto pPrev(bndNode->periodicPrev());
                const auto pNext(bndNode->periodicNext());
//                VerboseDislocationLoopNode(4,"pPrev= "<<pPrev->tag()<<" @ "<<pPrev->get_P().transpose()<<std::endl;);
//                VerboseDislocationLoopNode(4,"pNext= "<<pNext->tag()<<" @ "<<pNext->get_P().transpose()<<std::endl;);
                
                const auto pPrevLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
//                VerboseDislocationLoopNode(4,"pPrevLocal= "<<pPrevLocal.transpose()<<std::endl;);
                const auto pNextLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
//                VerboseDislocationLoopNode(4,"pNextLocal= "<<pNextLocal.transpose()<<std::endl;);
//                VerboseDislocationLoopNode(4,"periodicPlaneEdge->source= "<<bndNode->periodicPlaneEdge->source->transpose()<<std::endl;);
//                VerboseDislocationLoopNode(4,"bndNode->periodicPlaneEdge->sink= "<<bndNode->periodicPlaneEdge->sink->transpose()<<std::endl;);
//
                SegmentSegmentDistance<dim> ssd3(pPrev->get_P(),pNext->get_P(),bndNode->periodicPlaneEdge->meshIntersection->P0,bndNode->periodicPlaneEdge->meshIntersection->P1);
//                VerboseDislocationLoopNode(4,"ssd3.dMin= "<<ssd3.dMin<<std::endl;);
                
                
                SegmentSegmentDistance<dim-1> ssd(pPrevLocal,pNextLocal,*bndNode->periodicPlaneEdge->source,*bndNode->periodicPlaneEdge->sink);
                if(ssd.dMin<FLT_EPSILON)
                {
//                    VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                    bndNode->set_P(VectorLowerDim(0.5*(ssd.x0+ssd.x1)));
                }
                else
                {
                    bndNode->set_P(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                    DN.danglingBoundaryLoopNodes.insert(bndNode.get());
                }
                
            }
        }

        DN.updateBoundaryNodes();

        model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
    }
    
    template <typename DislocationNetworkType>
    void DislocationNetworkRemesh<DislocationNetworkType>::remeshByExpansion()
    {
        const auto t0= std::chrono::system_clock::now();
        model::cout<<"        Remeshing network: expanding... "<<std::flush;
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
                //                const VectorDimD dv(linkIter->second.sink->get_V()-linkIter->second.source->get_V());
                
                
                // Always expand single FR source segment
                //                    if (sharedLink->source->openOrder()==1 && sharedLink->sink->openOrder()==1)
                //                    {
                //                        toBeExpanded.insert(linkIter->second.nodeIDPair);
                //                    }
                
                // Expand pin points
//                if (   sharedLink->source->constraintNormals().size()>2
//                    && sharedLink->  sink->constraintNormals().size()>2
//                    && chordLength>3.0*Lmin)
//                {
//                    toBeExpanded.insert(sharedLink->nodeIDPair);
//                }
                
//                if (!sharedLink->source->isSimple() && !sharedLink->sink->isSimple()
                    if (sharedLink->source->loopNodes().size()>1 && sharedLink->sink->loopNodes().size()>1
                    /*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/ && chordLength>3.0*Lmin)
                { // also expands a straight line to generate glissile segment
//                    toBeExpanded.insert(sharedLink->nodeIDPair);
                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID,sharedLink->sink->sID));

                }
                
                // Expand segments shorter than Lmax
                if (chordLength>Lmax)
                {
//                    toBeExpanded.insert(sharedLink->nodeIDPair);
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
        model::cout<<" ("<<Nexpanded<<" expanded)"<<std::flush;
        model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        
        }

//    template <typename DislocationNetworkType>
//    void DislocationNetworkRemesh<DislocationNetworkType>::contract0chordSegments()
//    {
//        model::cout<<"        Contracting zero-chord segments... "<<std::flush;
//        const auto t0= std::chrono::system_clock::now();
//
//        std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
//
//        //            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
//        for (const auto& linkIter : DN.networkLinks())
//        {
//            VectorDimD chord(linkIter.second.lock()->chord()); // this is sink->get_P() - source->get_P()
//            const double chordLength(chord.norm());
//            if (chordLength<=FLT_EPSILON
//                && (linkIter.second.lock()->source->get_V()-linkIter.second.lock()->sink->get_V()).squaredNorm()<FLT_EPSILON)
//            {// toBeContracted part
//                toBeContracted.insert(std::make_pair(chordLength,std::make_pair(linkIter.second.lock()->source->sID,linkIter.second.lock()->sink->sID));
//            }
//        }
//
//        // Call Network::contract
//        unsigned int Ncontracted(0);
//        for (const auto& smallIter : toBeContracted)
//        {
//            const size_t i(smallIter.second.first);
//            const size_t j(smallIter.second.second);
//            const IsConstNetworkLinkType Lij(DN.link(i,j));
//
//            if (Lij.first )
//            {
//                Ncontracted+=DN.contract(Lij.second->source,Lij.second->sink);
//            }
//        }
//        model::cout<<"("<<Ncontracted<<" contracted)"<<std::flush;
//        model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//
//    }


    
} // namespace model
#endif

