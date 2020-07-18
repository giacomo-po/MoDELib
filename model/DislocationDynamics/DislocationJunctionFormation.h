/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONJUNCTIONFORMATION_H_
#define model_DISLOCATIONJUNCTIONFORMATION_H_

#include <utility> // for std::pair
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <SegmentSegmentDistance.h>
#include <SweepPlane.h>

#include <DislocationNetworkRemesh.h>
#include <MPIcout.h>
#include <EqualIteratorRange.h>
#include <N2IteratorRange.h>

#ifndef NDEBUG
#define VerboseJunctions(N,x) if(verboseJunctions>=N){model::cout<<x;}
#else
#define VerboseJunctions(N,x)
#endif


namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationJunctionFormation
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename DislocationNetworkType::IsNetworkEdgeType IsNetworkLinkType;
        typedef typename DislocationNetworkType::IsNodeType IsNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef std::pair<size_t,size_t> EdgeIDType;
        
        typedef std::tuple<EdgeIDType,EdgeIDType,SegmentSegmentDistance<dim>> IntersectionType;
        typedef std::deque<IntersectionType> IntersectionTypeContainerType;
        
        /**********************************************************************/
        void insertIntersection(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                                const LinkType* const linkA,
                                const LinkType* const linkB,
                                const SegmentSegmentDistance<dim>& ssd,
                                const double& currentcCollisionTOL,
                                const bool& bndJunction,const bool& gbndJunction)
        {
            VerboseJunctions(2,"insertIntersection "<<linkA->tag()<<","<<linkB->tag()<<std::endl;);
            VerboseJunctions(2,"insertIntersection "<<"ssd.dMin="<<ssd.dMin<<std::endl;);
            VerboseJunctions(2,"insertIntersection "<<"currentcCollisionTOL="<<currentcCollisionTOL<<std::endl;);
            
            if(ssd.dMin<currentcCollisionTOL)
            {
                bool isValidJunction((bndJunction || gbndJunction) && DN.simulationParameters.simulationType!=2);
                if(   !isValidJunction
                   && !linkA->isBoundarySegment()
                   && !linkB->isBoundarySegment()
                   && !linkA->isGrainBoundarySegment()
                   && !linkB->isGrainBoundarySegment())
                {// Check force condition for internal segments
                    
                    const VectorDim chordA(linkA->sink->get_P()-linkA->source->get_P());
                    const double LA(chordA.norm());
                    const VectorDim chordB(linkB->sink->get_P()-linkB->source->get_P());
                    const double LB(chordB.norm());
                    
                    if(LA>FLT_EPSILON && LB>FLT_EPSILON)
                    {
                        StressStraight<dim> stressA(ssd.x0-infiniteLineLength/LA*chordA,
                                                    ssd.x0+infiniteLineLength/LA*chordA,
                                                    linkA->burgers());
                        
                        StressStraight<dim> stressB(ssd.x1-infiniteLineLength/LB*chordB,
                                                    ssd.x1+infiniteLineLength/LB*chordB,
                                                    linkB->burgers());
                        
                        const VectorDim forceOnA=(stressB.stress(ssd.x0)*linkA->burgers()).cross(chordA);
                        const VectorDim forceOnB=(stressA.stress(ssd.x1)*linkB->burgers()).cross(chordB);
                        
                        if(forceOnA.dot(ssd.x1-ssd.x0)>0.0 && forceOnB.dot(ssd.x1-ssd.x0)<0.0)
                        {
                            VerboseJunctions(3,"attractive pair"<<std::endl;);
                            isValidJunction=true; // for non-parallel lines this neglects the energy of rotation
                        }
                        else
                        {
                            VerboseJunctions(3,"non-attractive pair"<<std::endl;);
                            
                        }
                    }
                }
                
                if(isValidJunction)
                {
#ifdef _OPENMP
                    intersectionContainer[omp_get_thread_num()].emplace_back(linkA->nodeIDPair,
                                                                             linkB->nodeIDPair,
                                                                             ssd);
#else
                    intersectionContainer[0].emplace_back(linkA->nodeIDPair,
                                                          linkB->nodeIDPair,
                                                          ssd);
#endif
                }
            }
            else
            {
                VerboseJunctions(3,"dMin="<<ssd.dMin<<", collisionTol="<<currentcCollisionTOL<<std::endl;);
            }
        }
        
        /**********************************************************************/
        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                               const size_t& nThreads)
        
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Finding collisions "<<std::flush;
            
            // Use SweepPlane to compute possible intersections
            SweepPlane<LinkType,dim> swp;
            for(const auto& link : DN.links())
            {
                if(   (!link.second->hasZeroBurgers() && !link.second->isVirtualBoundarySegment())
                   //                   || (link.second->isBoundarySegment() && DN.useVirtualExternalLoops))
                   || link.second->isBoundarySegment() )
                    
                {
                    //                    swp.addSegment(link.second->source->get_P()(0),link.second->source->get_P()(1),*link.second);
                    swp.addSegment(link.second->source->get_P()(0),link.second->sink->get_P()(0),*link.second);
                }
            }
            swp.computeIntersectionPairs();
            model::cout<<"("<<swp.potentialIntersectionPairs().size()<<" sweep-line pairs) "<<defaultColor<<std::flush;
            
            
            std::deque<std::pair<const LinkType*,const LinkType*>> reducedIntersectionPairs;
            for(size_t k=0;k<swp.potentialIntersectionPairs().size();++k)
            {
                const auto& linkA(swp.potentialIntersectionPair(k).first);
                const auto& linkB(swp.potentialIntersectionPair(k).second);
                
                const VectorDim& sourceA(linkA->source->get_P());
                const VectorDim&   sinkA(linkA->sink->get_P());
                const VectorDim& sourceB(linkB->source->get_P());
                const VectorDim&   sinkB(linkB->sink->get_P());
                
                // SweepPlane sweeps along x direction, so now check possible overlap in y and z direction
                bool yDontOverlap((std::min(sourceA(1),sinkA(1))>std::max(sourceB(1),sinkB(1))+collisionTol) || (std::max(sourceA(1),sinkA(1))<std::min(sourceB(1),sinkB(1))-collisionTol));
                bool zDontOverlap((std::min(sourceA(2),sinkA(2))>std::max(sourceB(2),sinkB(2))+collisionTol) || (std::max(sourceA(2),sinkA(2))<std::min(sourceB(2),sinkB(2))-collisionTol));
                
                if(!yDontOverlap && !zDontOverlap)
                {
                    reducedIntersectionPairs.emplace_back(linkA,linkB);
                }
            }
            
            
            model::cout<<" ("<<reducedIntersectionPairs.size()<<" reduced pairs) "<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            const auto t1= std::chrono::system_clock::now();
            model::cout<<"        Selecting junctions ("<<nThreads<<" threads): "<<std::flush;
            
            //! 2- loop over all links and determine their intersections
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<reducedIntersectionPairs.size();++k)
            {
                
                
                
                const auto& linkA(reducedIntersectionPairs[k].first);
                const auto& linkB(reducedIntersectionPairs[k].second);
                
                VerboseJunctions(2,"checking intersection "<<linkA->tag()<<","<<linkB->tag()<<std::endl;);
                
                
                const bool linkAisBnd(linkA->isBoundarySegment());
                const bool linkBisBnd(linkB->isBoundarySegment());
                const bool linkAisGBnd(linkA->isGrainBoundarySegment());
                const bool linkBisGBnd(linkB->isGrainBoundarySegment());
                
                
                const bool bndJunction(   (linkAisBnd || linkBisBnd)
                                       &&  linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).norm()<FLT_EPSILON // same plane
                                       &&  linkA->chord().cross(linkB->chord()).norm()<FLT_EPSILON // colinear on boundary
                                       );
                
                const bool gbndJunction(   (linkAisGBnd || linkBisGBnd)
                                        &&  linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).norm()<FLT_EPSILON
                                        &&  linkA->chord().cross(linkB->chord()).norm()<FLT_EPSILON // colinear on grain-boundary
                                        );
                
                VerboseJunctions(2,"links are bnd: "<<linkAisBnd<<" "<<linkBisBnd<<std::endl;);
                VerboseJunctions(2,"links are grainBnd: "<<linkAisGBnd<<" "<<linkBisGBnd<<std::endl;);
                VerboseJunctions(2,"bndJunction: "<<bndJunction<<std::endl;);
                VerboseJunctions(2,"gbndJunction: "<<gbndJunction<<std::endl;);
                
                const bool intersectionIsSourceSource(linkA->source->sID==linkB->source->sID);
                const bool intersectionIsSourceSink(  linkA->source->sID==linkB->  sink->sID);
                const bool intersectionIsSinkSource(  linkA->  sink->sID==linkB->source->sID);
                const bool intersectionIsSinkSink(    linkA->  sink->sID==linkB->  sink->sID);
                
                
                
                double currentcCollisionTOL=collisionTol;
                
                //                std::set<const GlidePlane<dim>*> commonGlidePlanes;
                //                std::set_intersection(linkA->glidePlanes().begin(),linkA->glidePlanes().end(),linkB->glidePlanes().begin(),linkB->glidePlanes().end(),std::inserter(commonGlidePlanes,commonGlidePlanes.begin()));
                //                if(commonGlidePlanes.size())
                //                {// links have a GlidePlane in common
                //                    const double cosTheta=(intersectionIsSourceSource||intersectionIsSinkSink)? linkA->chord().normalized().dot(linkB->chord().normalized()) : ((intersectionIsSourceSink ||intersectionIsSinkSource)? -linkA->chord().normalized().dot(linkB->chord().normalized()) : -1.0);
                //                    if(cosTheta<0.7)
                //                    {// links are either disconnected, or connected to a node and forming a large angle at that node. Reduce tolerance
                //                        currentcCollisionTOL=FLT_EPSILON;
                //                    }
                //                }
                
                const auto pcPlanes(linkA->parallelAndCoincidentGlidePlanes(linkB->glidePlanes()));
                VerboseJunctions(2,"pcPlanes.size()= "<<pcPlanes.size()<<std::endl;);
                for(const auto& pair : pcPlanes)
                {
                    if(pair.first==pair.second)
                    {// a common coincident plane found
                        const double cosTheta=(intersectionIsSourceSource||intersectionIsSinkSink)? linkA->chord().normalized().dot(linkB->chord().normalized()) : ((intersectionIsSourceSink ||intersectionIsSinkSource)? -linkA->chord().normalized().dot(linkB->chord().normalized()) : -1.0);
                        if(cosTheta<0.7)
                        {// links are either disconnected, or connected to a node and forming a large angle at that node. Reduce tolerance
                            currentcCollisionTOL=FLT_EPSILON;
                        }
                        break;
                    }
                    else
                    {// segments on parallel (non-coincident) planes
                        currentcCollisionTOL=FLT_EPSILON;
                    }
                }
                
                
                
                
                //                if(   linkA->glidePlaneNormal().squaredNorm()>FLT_EPSILON
                //                   && linkB->glidePlaneNormal().squaredNorm()>FLT_EPSILON
                //                   && linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).squaredNorm()<FLT_EPSILON)
                //                {// segments on parallel or coincident planes, reduce tolerance
                //
                //                    const double cosTheta=(intersectionIsSourceSource||intersectionIsSinkSink)? linkA->chord().normalized().dot(linkB->chord().normalized()) : ((intersectionIsSourceSink ||intersectionIsSinkSource)? -linkA->chord().normalized().dot(linkB->chord().normalized()) : -1.0);
                //
                //                    if(cosTheta<0.7)
                //                    {
                //                        currentcCollisionTOL=FLT_EPSILON;
                //                    }
                //
                //                }
                
                
                if(intersectionIsSourceSource)
                {
                    if(!linkA->source->isSimple())
                    {// non-simple common node, increase tolerance
                        currentcCollisionTOL=0.5*collisionTol;
                    }
                    
                    // intersect sink of A with link B
                    SegmentSegmentDistance<dim> ssdA(linkA->sink->get_P(), // fake degenerete segment at sink of A
                                                     linkA->sink->get_P(),  // fake degenerete segment at sink of A
                                                     linkB->source->get_P(),
                                                     linkB->sink->get_P(),
                                                     1.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                    // intersect sink of B with link A
                    SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                     linkA->sink->get_P(),
                                                     linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                     linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                     1.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                }
                else if(intersectionIsSourceSink)
                {
                    if(!linkA->source->isSimple())
                    {// non-simple common node, increase tolerance
                        currentcCollisionTOL=0.5*collisionTol;
                    }
                    
                    // intersect sink of A with link B
                    SegmentSegmentDistance<dim> ssdA(linkA->sink->get_P(), // fake degenerete segment at sink of A
                                                     linkA->sink->get_P(),  // fake degenerete segment at sink of A
                                                     linkB->source->get_P(),
                                                     linkB->sink->get_P(),
                                                     1.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                    // intersect source of B with link A
                    SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                     linkA->sink->get_P(),
                                                     linkB->source->get_P(),    // fake degenerete segment at source of B
                                                     linkB->source->get_P(),    // fake degenerete segment at source of B
                                                     0.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                    
                }
                else if(intersectionIsSinkSource)
                {
                    if(!linkA->sink->isSimple())
                    {// non-simple common node, increase tolerance
                        currentcCollisionTOL=0.5*collisionTol;
                    }
                    
                    // intersect source of A with link B
                    SegmentSegmentDistance<dim> ssdA(linkA->source->get_P(), // fake degenerete segment at source of A
                                                     linkA->source->get_P(),  // fake degenerete segment at source of A
                                                     linkB->source->get_P(),
                                                     linkB->sink->get_P(),
                                                     0.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                    
                    // intersect sink of B with link A
                    SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                     linkA->sink->get_P(),
                                                     linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                     linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                     1.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction);
                }
                else if(intersectionIsSinkSink)
                {
                    
                    if(!linkA->sink->isSimple())
                    {// non-simple common node, increase tolerance
                        currentcCollisionTOL=0.5*collisionTol;
                    }
                    
                    // intersect source of A with link B
                    SegmentSegmentDistance<dim> ssdA(linkA->source->get_P(), // fake degenerete segment at source of A
                                                     linkA->source->get_P(),  // fake degenerete segment at source of A
                                                     linkB->source->get_P(),
                                                     linkB->sink->get_P(),
                                                     0.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                    
                    // intersect source of B with link A
                    SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                     linkA->sink->get_P(),
                                                     linkB->source->get_P(),    // fake degenerete segment at source of B
                                                     linkB->source->get_P(),    // fake degenerete segment at source of B
                                                     0.0);
                    insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction);
                    
                }
                else
                {// segments don't share vertices
                    
                    SegmentSegmentDistance<dim> ssd(linkA->source->get_P(),
                                                    linkA->sink->get_P(),
                                                    linkB->source->get_P(),
                                                    linkB->sink->get_P()
                                                    );
                    insertIntersection(intersectionContainer,linkA,linkB,ssd,currentcCollisionTOL,bndJunction,gbndJunction);
                }
                //                }
                
            }
            
            int nIntersections=0;
            for (const auto& intersectionByThreadContainer : intersectionContainer)
            {
                nIntersections+=intersectionByThreadContainer.size();
            }
            model::cout<<nIntersections<<" physical junctions "<<std::flush;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        std::pair<std::shared_ptr<NodeType>,bool> junctionNode(const double& t,
                                                               const VectorDim& x,
                                                               const IsNetworkLinkType& L,
                                                               const EdgeIDType& key)
        {
            VerboseJunctions(4,"JunctionNode for segment: "<<key.first<<"->"<<key.second<<" @ "<<t<<std::endl;);
            
            auto  Nclose= t <0.5? L.second->source : L.second->sink;
            auto  Nfar  = t>=0.5? L.second->source : L.second->sink;
            VerboseJunctions(4,"Nclose="<<Nclose->sID<<std::endl;);
            VerboseJunctions(4,"Nfar="<<Nfar->sID<<std::endl;);
            if(t>FLT_EPSILON && t<1.0-FLT_EPSILON)
            {// intersection point is not an end node
                if((Nclose->get_P()-x).norm()>DN.networkRemesher.Lmin)
                {
                    VerboseJunctions(4,"JunctionNode case a"<<std::endl;);
                    return std::make_pair(DN.expand(key.first,key.second,t),true);
                }
                else
                {
                    if(Nclose->isMovableTo(x))
                    {
                        VerboseJunctions(4,"JunctionNode case b"<<std::endl;);
                        return std::make_pair(Nclose,false);
                    }
                    else
                    {
                        if((Nfar->get_P()-x).norm()>DN.networkRemesher.Lmin)
                        {
                            VerboseJunctions(4,"JunctionNode case c"<<std::endl;);
                            return std::make_pair(DN.expand(key.first,key.second,t),true);
                        }
                        else
                        {
                            if(Nfar->isMovableTo(x))
                            {
                                VerboseJunctions(4,"JunctionNode case d"<<std::endl;);
                                return std::make_pair(Nfar,false);
                            }
                            else
                            {
                                VerboseJunctions(4,"JunctionNode case e"<<std::endl;);
                                return std::make_pair(DN.expand(key.first,key.second,t),true);
                            }
                        }
                    }
                }
            }
            else
            {
                VerboseJunctions(4,"JunctionNode case f"<<std::endl;);
                return std::make_pair(Nclose,false);
            }
        }
        
        /**********************************************************************/
        size_t junctionStep()
        {
            
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            std::deque<IntersectionTypeContainerType> intersectionContainer;
            intersectionContainer.resize(nThreads);
            findIntersections(intersectionContainer,nThreads);
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Forming Junctions: "<<std::flush;
            
            size_t nContracted=0;
            for (const auto& intersectionByThreadContainer : intersectionContainer)
            {
                for (const auto& intersection : intersectionByThreadContainer)
                {
                    const EdgeIDType& key1(std::get<0>(intersection));
                    const EdgeIDType& key2(std::get<1>(intersection));
                    const SegmentSegmentDistance<dim>& ssd(std::get<2>(intersection));
                    const double& t(ssd.t);
                    const double& u(ssd.u);
                    
                    const IsNetworkLinkType L1(DN.link(key1.first,key1.second));
                    const IsNetworkLinkType L2(DN.link(key2.first,key2.second));
                    
                    if(L1.first && L2.first) // Links exist
                    {
                        
                        
                        
                        
                        VerboseJunctions(1,"forming junction "<<key1.first<<"->"<<key1.second<<", "
                                         /*                   */ <<key2.first<<"->"<<key2.second<<", "
                                         /*                   */ <<"dMin="<<ssd.dMin<<", "
                                         /*                   */ <<"@ ("<<t<<","<<u<<"), "<<std::flush);
                        std::pair<std::shared_ptr<NodeType>,bool> Ni(junctionNode(t,ssd.x0,L1,key1));
                        std::pair<std::shared_ptr<NodeType>,bool> Nj(junctionNode(u,ssd.x1,L2,key2));
                        
                        VerboseJunctions(1,"contracting "<<Ni.first->sID<<" "<<Nj.first->sID<<std::endl;);
                        
                        if(Ni.first->sID!=Nj.first->sID)
                        {
                            const bool success=DN.contract(Ni.first,Nj.first);
                            nContracted+=success;
                            if(!success)
                            {
                                if(Ni.second)
                                {
                                    DN.remove(Ni.first->sID);
                                }
                                if(Nj.second)
                                {
                                    DN.remove(Nj.first->sID);
                                }
                            }
                            else
                            {// first contraction happended. We want to generate a finite-length junction
                                
                                
                            }
                        }
                    }
                }
            }
            model::cout<<" ("<<nContracted<<" contracted)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            return nContracted;
        }
        
        /**********************************************************************/
        void glissileJunctions(const double &dx)
        {
            const auto t0 = std::chrono::system_clock::now();
            model::cout << "        Forming Glissile Junctions: " << std::flush;
            
            std::deque<std::tuple<std::shared_ptr<NodeType>, std::shared_ptr<NodeType>, size_t, size_t>> glissDeq;
            
            std::deque<std::tuple<std::shared_ptr<NodeType>, std::shared_ptr<NodeType>, std::shared_ptr<NodeType>>> expDeq;
            
            for (const auto &link : DN.links())
            {
                
                if (link.second->isSessile() && link.second->loopLinks().size() > 1) // a junction
                {
                    const VectorDim chord(link.second->sink->get_P() - link.second->source->get_P());
                    const double chordNorm(chord.norm());
                    
                    
                    if (fabs(link.second->burgers().norm() - 1.0) < FLT_EPSILON // a non-zero link with minimum Burgers
                        && chordNorm > dx)
                    {
                        
                        const VectorDim unitChord(chord / chordNorm);
                        
                        
                        if (!link.second->isGrainBoundarySegment() && !link.second->isBoundarySegment())
                        {
                            
                            VerboseJunctions(2,"glissele junction, segment "<<link.second->tag()<<std::endl;);

                            for (const auto &gr : link.second->grains())
                            {
                                for (size_t s = 0; s < gr->slipSystems().size(); ++s)
                                {
                                    const auto &slipSystem(gr->slipSystems()[s]);
                                    if ((slipSystem->s.cartesian() - link.second->burgers()).norm() < FLT_EPSILON && fabs(slipSystem->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
                                    {
                                        VerboseJunctions(3,"glissDeq, emplacing"<<std::endl;);
                                        
                                        glissDeq.emplace_back(link.second->source, link.second->sink, gr->grainID, s);
                                    }
                                }
                            }
                            
                            //                            std::set<const LoopType*> sourceLoops;
//                            int inserted(0);
//                            for (const auto &nodelink : link.second->source->outLoopLinks())
//                            {
//                                if (nodelink->pLink->loopLinks().size() == 1)
//                                {
//                                    const auto &loop(nodelink->loop());
//                                    if(loop->slipSystem())
//                                    {
//                                        if ((loop->slipSystem()->s.cartesian() + link.second->burgers()).norm() < FLT_EPSILON && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
//                                        {
//                                            expDeq.emplace_back(nodelink->pLink->source, nodelink->pLink->sink, link.second->sink);
//                                            inserted++;
//                                        }
//                                    }
//                                }
//                            }
//                            VerboseJunctions(3,"source->outLoopLinks(), inserted= "<<inserted<<std::endl;);
//
//
//                            for (const auto &nodelink : link.second->source->inLoopLinks())
//                            {
//                                if (nodelink->pLink->loopLinks().size() == 1)
//                                {
//                                    const auto &loop(nodelink->loop());
//                                    if (loop->slipSystem())
//                                    {
//                                        if ((loop->slipSystem()->s.cartesian() - link.second->burgers()).norm() < FLT_EPSILON && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
//                                        {
//                                            expDeq.emplace_back(nodelink->pLink->source, nodelink->pLink->sink, link.second->sink);
//                                            inserted++;
//                                        }
//                                    }
//                                }
//                            }
//                            VerboseJunctions(3,"source->inLoopLinks(), inserted= "<<inserted<<std::endl;);
//
//                            assert(inserted <= 1);
//
//                            if (inserted == 0)
//                            {
//                                for (const auto &nodelink : link.second->sink->outLoopLinks())
//                                {
//                                    if (nodelink->pLink->loopLinks().size() == 1)
//                                    {
//                                        const auto &loop(nodelink->loop());
//                                        if (loop->slipSystem())
//                                        {
//                                            if ((loop->slipSystem()->s.cartesian() + link.second->burgers()).norm() < FLT_EPSILON && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
//                                            {
//                                                expDeq.emplace_back(nodelink->pLink->source, nodelink->pLink->sink, link.second->source);
//                                                inserted++;
//                                            }
//                                        }
//                                    }
//                                }
//                                VerboseJunctions(3,"sink->inLoopLinks(), inserted= "<<inserted<<std::endl;);
//
//
//                                for (const auto &nodelink : link.second->sink->inLoopLinks())
//                                {
//                                    if (nodelink->pLink->loopLinks().size() == 1)
//                                    {
//                                        const auto &loop(nodelink->loop());
//                                        if (loop->slipSystem())
//                                        {
//                                            if ((loop->slipSystem()->s.cartesian() - link.second->burgers()).norm() < FLT_EPSILON && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
//                                            {
//                                                expDeq.emplace_back(nodelink->pLink->source, nodelink->pLink->sink, link.second->source);
//                                                inserted++;
//                                            }
//                                        }
//                                    }
//                                }
//                                VerboseJunctions(3,"sink->inLoopLinks(), inserted= "<<inserted<<std::endl;);
//
//                                assert(inserted <= 1);
//                            }
//
//                            if (inserted == 0)
//                            {
//                                for (const auto &gr : link.second->grains())
//                                {
//                                    for (size_t s = 0; s < gr->slipSystems().size(); ++s)
//                                    {
//                                        const auto &slipSystem(gr->slipSystems()[s]);
//                                        if ((slipSystem->s.cartesian() - link.second->burgers()).norm() < FLT_EPSILON && fabs(slipSystem->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
//                                        {
//                                            VerboseJunctions(3,"glissDeq, emplacing"<<std::endl;);
//
//                                            glissDeq.emplace_back(link.second->source, link.second->sink, gr->grainID, s);
//                                        }
//                                    }
//                                }
//                            }
                        }
                    }
                }
            }
            
            size_t formedJunctions = 0;
            
            for (const auto &tup : expDeq)
            {
                const std::shared_ptr<NodeType> &source(std::get<0>(tup));
                const std::shared_ptr<NodeType> &sink(std::get<1>(tup));
                const std::shared_ptr<NodeType> &exp(std::get<2>(tup));
                const size_t &sourceID(source->sID);
                const size_t &sinkID(sink->sID);
                const auto isLink(DN.link(sourceID, sinkID));
                
                if (isLink.first)
                {
                    VerboseJunctions(3,"expanding junction "<<isLink.second->tag()<<" @node "<<exp->sID<<std::endl;);
                    DN.expand(isLink.second, exp);
                    formedJunctions++;
                }
            }
            
            for (const auto &tup : glissDeq)
            {
                const std::shared_ptr<NodeType> &source(std::get<0>(tup));
                const std::shared_ptr<NodeType> &sink(std::get<1>(tup));
                const size_t &sourceID(source->sID);
                const size_t &sinkID(sink->sID);
                const size_t &grainID(std::get<2>(tup));
                const size_t &slipID(std::get<3>(tup));
                
                const auto isLink(DN.link(sourceID, sinkID));
                if (isLink.first)
                {
                    
                    const VectorDim newNodeP(0.5 * (isLink.second->source->get_P() + isLink.second->sink->get_P()));
                    const long int planeIndex(DN.poly.grain(grainID).slipSystems()[slipID]->n.closestPlaneIndexOfPoint(newNodeP));
                    const GlidePlaneKey<dim> glissilePlaneKey(planeIndex, DN.poly.grain(grainID).slipSystems()[slipID]->n);
                    const auto glidePlane(DN.glidePlaneFactory.get(glissilePlaneKey));
                    
                    std::shared_ptr<NodeType> newNode(new NodeType(&DN, glidePlane->snapToPlane(newNodeP), VectorDim::Zero(), 1.0));
                    
                    std::vector<std::shared_ptr<NodeType>> loopNodes;
                    
                    loopNodes.push_back(sink);   // insert in reverse order, sink first, source second
                    loopNodes.push_back(source); // insert in reverse order, sink first, source second
                    loopNodes.push_back(newNode);
                    
                    DN.insertLoop(loopNodes,
                                  DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
                                  glidePlane);
                    
                    formedJunctions++;
                }
            }
            model::cout << "(" << formedJunctions << " junctions)" << magentaColor << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]" << defaultColor << std::endl;
        }
        
//        /**********************************************************************/
//        void glissileJunctions(const double& dx)
//        {
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"        Forming Glissile Junctions: "<<std::flush;
//
//            std::deque<std::tuple<std::shared_ptr<NodeType>,std::shared_ptr<NodeType>,size_t,size_t>> glissDeq;
//
//            std::deque<std::tuple<std::shared_ptr<NodeType>,std::shared_ptr<NodeType>,std::shared_ptr<NodeType>>> expDeq;
//
//
//            for(const auto& link : DN.links())
//            {
//
//                if(   link.second->isSessile()
//                   && link.second->loopLinks().size()>1      // a junction
//                   )
//                {
//                    const VectorDim chord(link.second->sink->get_P()-link.second->source->get_P());
//                    const double chordNorm(chord.norm());
//
//
//                    if(   fabs(link.second->burgers().norm()-1.0)<FLT_EPSILON // a non-zero link with minimum Burgers
//                       && chordNorm>dx)
//                    {
//
//                        const VectorDim unitChord(chord/chordNorm);
//
//
//                        //                        //std::cout<<"link "<<link.second->source->sID<<"->"<<link.second->sink->sID<<std::endl;
//                        //                        //std::cout<<"burgers="<<link.second->burgers().transpose()<<std::endl;
//                        //                        //std::cout<<"chord="<<unitChord.transpose()<<std::endl;
//                        //
//                        //                        for(const auto& loopLink : link.second->loopLinks())
//                        //                        {
//                        //                            //std::cout<<"loopLink "<<loopLink->source()->sID<<"->"<<loopLink->sink()->sID<<", flow="<<loopLink->flow().cartesian().transpose()<<std::endl;
//                        //
//                        //                        }
//
//
//
//                        if(   !link.second->isGrainBoundarySegment()
//                           && !link.second->isBoundarySegment() )
//                        {
//
//
//                            //                            std::set<const LoopType*> sourceLoops;
//
//                            int inserted(0);
//
////                            for(const auto& link : link.second->source->outLoopLinks())
////                            {
////                                if(link->pLink->loopLinks().size()==1)
////                                {
////                                    const auto& loop(link->loop());
////                                    if(  (loop->slipSystem()->s.cartesian()+link.second->burgers()).norm()<FLT_EPSILON
////                                       && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
////                                    {
////                                        expDeq.emplace_back(link->pLink->source,link->pLink->sink,link.second->sink);
////                                        inserted++;
////                                    }
////                                }
////                            }
////
////                            for(const auto& link : source->inLoopLinks())
////                            {
////                                if(link->pLink->loopLinks().size()==1)
////                                {
////                                    const auto& loop(link->loop());
////                                    if(  (loop->slipSystem()->s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON
////                                       && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
////                                    {
////                                        expDeq.emplace_back(link->pLink->source,link->pLink->sink,sink);
////                                        inserted++;
////                                    }
////                                }
////                            }
////                            assert(inserted<1);
////
////                            if(inserted==0)
////                            {
////                                for(const auto& link : sink->outLoopLinks())
////                                {
////                                    if(link->pLink->loopLinks().size()==1)
////                                    {
////                                        const auto& loop(link->loop());
////                                        if(  (loop->slipSystem()->s.cartesian()+link.second->burgers()).norm()<FLT_EPSILON
////                                           && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
////                                        {
////                                            expDeq.emplace_back(link->pLink->source,link->pLink->sink,source);
////                                            inserted++;
////                                        }
////                                    }
////                                }
////
////                                for(const auto& link : sink->inLoopLinks())
////                                {
////                                    if(link->pLink->loopLinks().size()==1)
////                                    {
////                                        const auto& loop(link->loop());
////                                        if(  (loop->slipSystem()->s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON
////                                           && fabs(loop->slipSystem()->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
////                                        {
////                                            expDeq.emplace_back(link->pLink->source,link->pLink->sink,source);
////                                            inserted++;
////                                        }
////                                    }
////                                }
////                                assert(inserted<1);
////                            }
//
//                            if(inserted==0)
//                            {
//                                for(const auto& gr : link.second->grains())
//                                {
//                                    for(size_t s=0;s<gr->slipSystems().size();++s)
//                                    {
//                                        const auto& slipSystem(gr->slipSystems()[s]);
//                                        if(  (slipSystem->s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON
//                                           && fabs(slipSystem->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
//                                        {
//                                            glissDeq.emplace_back(link.second->source,link.second->sink,gr->grainID,s);
//                                        }
//                                    }
//                                }
//                            }
//
//                        }
//                    }
//                }
//            }
//
//            size_t formedJunctions=0;
//
//
//            for(const auto& tup : expDeq)
//            {
//                const std::shared_ptr<NodeType>& source(std::get<0>(tup));
//                const std::shared_ptr<NodeType>& sink(std::get<1>(tup));
//                const std::shared_ptr<NodeType>& exp(std::get<2>(tup));
//                const size_t& sourceID(source->sID);
//                const size_t& sinkID(sink->sID);
//                const auto isLink(DN.link(sourceID,sinkID));
//
//
//                if(isLink.first)
//                {
//                    DN.expand(isLink.second,exp);
//                }
//            }
//
//
//            for(const auto& tup : glissDeq)
//            {
//                const std::shared_ptr<NodeType>& source(std::get<0>(tup));
//                const std::shared_ptr<NodeType>& sink(std::get<1>(tup));
//                const size_t& sourceID(source->sID);
//                const size_t& sinkID(sink->sID);
//                const size_t& grainID(std::get<2>(tup));
//                const size_t& slipID(std::get<3>(tup));
//
//                const auto isLink(DN.link(sourceID,sinkID));
//                if(isLink.first)
//                {
//
//                    const VectorDim newNodeP(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()));
//                    const long int planeIndex(DN.poly.grain(grainID).slipSystems()[slipID]->n.closestPlaneIndexOfPoint(newNodeP));
//                    const GlidePlaneKey<dim> glissilePlaneKey(planeIndex,DN.poly.grain(grainID).slipSystems()[slipID]->n);
//                    const auto glidePlane(DN.glidePlaneFactory.get(glissilePlaneKey));
//
//                    std::shared_ptr<NodeType> newNode(new NodeType(&DN,glidePlane->snapToPlane(newNodeP),VectorDim::Zero(),1.0));
//
//                    std::vector<std::shared_ptr<NodeType>> loopNodes;
//
//                    loopNodes.push_back(sink);      // insert in reverse order, sink first, source second
//                    loopNodes.push_back(source);    // insert in reverse order, sink first, source second
//                    loopNodes.push_back(newNode);
//
//                    DN.insertLoop(loopNodes,
//                                  DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                                  glidePlane);
//
//                    //                    std::set<const LoopType*> sourceLoops;
//                    //                    for(const auto& loop : source->loops())
//                    //                    {
//                    //                        if(loop->slipSystem())
//                    //                        {
//                    //                            if(  (loop->slipSystem()->s.cartesian()-DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian()).norm()<FLT_EPSILON
//                    //                               && fabs(slipSystem->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
//                    //                            if(loop->slipSystem()==DN.poly.grain(grainID).slipSystems()[slipID])
//                    //                            {
//                    //                                sourceLoops.insert(loop);
//                    //                            }
//                    //                        }
//                    //                    }
//                    //
//                    //                    std::set<const LoopType*> sinkLoops;
//                    //                    for(const auto& loop : sink->loops())
//                    //                    {
//                    //                        if(loop->slipSystem())
//                    //                        {
//                    //                            if(loop->slipSystem()==DN.poly.grain(grainID).slipSystems()[slipID])
//                    //                            {
//                    //                                sinkLoops.insert(loop);
//                    //                            }
//                    //                        }
//                    //                    }
//                    //
//                    //                    if(sourceLoops.size()==1 && sinkLoops.size()==1)
//                    //                    {
//                    //                        DN.insertLoop(loopNodes,
//                    //                                      DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                    //                                      glidePlane);
//                    //                    }
//                    //                    else if(sourceLoops.size()==1 && sinkLoops.size()!=1)
//                    //                    {// source meets expand condition
//                    //                        const auto linksSet(source->linksByLoopID()[(*sourceLoops.begin())->sID]);
//                    //                        if(   (*linksSet. begin())->pLink->loopLinks().size()==1
//                    //                           && (*linksSet.rbegin())->pLink->loopLinks().size()!=1)
//                    //                        {
//                    //                            DN.expand((*linksSet. begin())->pLink.get(),sink);
//                    //                        }
//                    //                        else if(   (*linksSet. begin())->pLink->loopLinks().size()!=1
//                    //                                && (*linksSet.rbegin())->pLink->loopLinks().size()==1)
//                    //                        {
//                    //                            DN.expand((*linksSet.rbegin())->pLink.get(),sink);
//                    //                        }
//                    //                        else
//                    //                        {
//                    //                            DN.insertLoop(loopNodes,
//                    //                                          DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                    //                                          glidePlane);
//                    //                        }
//                    //                    }
//                    //                    else if(sourceLoops.size()!=1 && sinkLoops.size()==1)
//                    //                    {
//                    //                        const auto linksSet(sink->linksByLoopID()[(*sinkLoops.begin())->sID]);
//                    //                        if(   (*linksSet. begin())->pLink->loopLinks().size()==1
//                    //                           && (*linksSet.rbegin())->pLink->loopLinks().size()!=1)
//                    //                        {
//                    //                            DN.expand((*linksSet. begin())->pLink.get(),source);
//                    //                        }
//                    //                        else if(   (*linksSet. begin())->pLink->loopLinks().size()!=1
//                    //                                && (*linksSet.rbegin())->pLink->loopLinks().size()==1)
//                    //                        {
//                    //                            DN.expand((*linksSet.rbegin())->pLink.get(),source);
//                    //                        }
//                    //                        else
//                    //                        {
//                    //                            DN.insertLoop(loopNodes,
//                    //                                          DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                    //                                          glidePlane);
//                    //                        }
//                    //                    }
//                    //                    else
//                    //                    {
//                    //                        DN.insertLoop(loopNodes,
//                    //                                      DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                    //                                      glidePlane);
//                    //                    }
//
//
//
//                    formedJunctions++;
//                }
//
//
//
//            }
//            model::cout<<"("<<formedJunctions<<" junctions)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//        }
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        
        static double collisionTol;     //! The tolerance (in units of distance) used for collision detection
        const size_t maxJunctionIterations;
        const int verboseJunctions;
        const double infiniteLineLength;
        
        /**********************************************************************/
        DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
        /* init */,verboseJunctions(TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseJunctions",true))
        /* init */,infiniteLineLength(10000.0)
        {
            
        }
        
        /**********************************************************************/
        void formJunctions(const double& dx)
        {
            size_t nContracted=1;
            size_t iterations=0;
            while(nContracted && iterations<maxJunctionIterations)
            {
                nContracted=junctionStep();
                glissileJunctions(dx);
                iterations++;
            }
        }
        
    };
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
}
#endif
