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
                            isValidJunction=true;
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
                
                VerboseJunctions(2,"checking intersection "<<linkA->source->sID<<"->"<<linkA->sink->sID
                                 /*                   */ <<" "<<linkB->source->sID<<"->"<<linkB->sink->sID<<std::endl;);
                
                
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
                
                //                const bool frankRule(linkA->burgers().dot(linkB->burgers())*linkA->chord().dot(linkB->chord())<=0.0);
                //
                //
                //                const bool isValidJunction(   (frankRule   && !linkAisBnd && !linkBisBnd && !linkAisGBnd && !linkBisGBnd)  // energy rule is satisfied for internal segments
                //                                           || bndJunction   // junciton between parallel boundary segments
                //                                           || gbndJunction  // junciton between parallel grain-boundary segments
                //                                           );
                
                //                bool isValidJunction(bndJunction || gbndJunction);
                //                if(!isValidJunction)
                //                {
                //
                //                    SegmentSegmentDistance<dim> ssd();
                //
                //
                //
                //                }
                
                VerboseJunctions(2,"links are bnd: "<<linkAisBnd<<" "<<linkBisBnd<<std::endl;);
                VerboseJunctions(2,"links are grainBnd: "<<linkAisGBnd<<" "<<linkBisGBnd<<std::endl;);
                VerboseJunctions(2,"bndJunction: "<<bndJunction<<std::endl;);
                VerboseJunctions(2,"gbndJunction: "<<gbndJunction<<std::endl;);
                //                VerboseJunctions(2,"isValidJunction: "<<isValidJunction<<std::endl;);
                
                
                
                //                if(isValidJunction)
                //                {
                
                const bool intersectionIsSourceSource(linkA->source->sID==linkB->source->sID);
                const bool intersectionIsSourceSink(  linkA->source->sID==linkB->  sink->sID);
                const bool intersectionIsSinkSource(  linkA->  sink->sID==linkB->source->sID);
                const bool intersectionIsSinkSink(    linkA->  sink->sID==linkB->  sink->sID);
                
                
                
                double currentcCollisionTOL=collisionTol;
                if(   linkA->glidePlaneNormal().squaredNorm()>FLT_EPSILON
                   && linkB->glidePlaneNormal().squaredNorm()>FLT_EPSILON
                   && linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).squaredNorm()<FLT_EPSILON)
                {// segments on parallel or coincident planes, reduce tolerance
                    
                    
                    
                    
                    const double cosTheta=(intersectionIsSourceSource||intersectionIsSinkSink)? linkA->chord().normalized().dot(linkB->chord().normalized()) : ((intersectionIsSourceSink ||intersectionIsSinkSource)? -linkA->chord().normalized().dot(linkB->chord().normalized()) : -1.0);
                    
                    if(cosTheta<0.7)
                    {
                        currentcCollisionTOL=FLT_EPSILON;
                    }
                    
                }
                
                
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
        void glissileJunctions(const double& dx)
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Forming Glissile Junctions: "<<std::flush;
            
            std::deque<std::tuple<std::shared_ptr<NodeType>,std::shared_ptr<NodeType>,size_t,size_t>> glissDeq;
            
            for(const auto& link : DN.links())
            {
                
                if(   link.second->isSessile()
                   && link.second->loopLinks().size()>1      // a junction
                   )
                {
                    const VectorDim chord(link.second->sink->get_P()-link.second->source->get_P());
                    const double chordNorm(chord.norm());
                    
                    
                    if(   fabs(link.second->burgers().norm()-1.0)<FLT_EPSILON // a non-zero link with minimum Burgers
                       && chordNorm>dx)
                    {
                        
                        const VectorDim unitChord(chord/chordNorm);
                        
                        
                        //                        std::cout<<"link "<<link.second->source->sID<<"->"<<link.second->sink->sID<<std::endl;
                        //                        std::cout<<"burgers="<<link.second->burgers().transpose()<<std::endl;
                        //                        std::cout<<"chord="<<unitChord.transpose()<<std::endl;
                        //
                        //                        for(const auto& loopLink : link.second->loopLinks())
                        //                        {
                        //                            std::cout<<"loopLink "<<loopLink->source()->sID<<"->"<<loopLink->sink()->sID<<", flow="<<loopLink->flow().cartesian().transpose()<<std::endl;
                        //
                        //                        }
                        
                        
                        
                        if(   !link.second->isGrainBoundarySegment()
                           && !link.second->isBoundarySegment() )
                        {
                            //                            std::cout<<"here 1"<<std::endl;
                            for(const auto& gr : link.second->grains())
                            {
                                //                                std::cout<<"here 2"<<std::endl;
                                
                                for(size_t s=0;s<gr->slipSystems().size();++s)
                                {
                                    const auto& slipSystem(gr->slipSystems()[s]);
                                    
                                    //                                    std::cout<<"here 3 "<<"\n"<<slipSystem->s.cartesian().transpose()<<"\n"<<link.second->burgers().transpose()<<std::endl;
                                    //                                    std::cout<<((slipSystem->s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON)<<std::endl;
                                    //                                    std::cout<<(fabs(slipSystem->n.cartesian().normalized().dot(unitChord)))<<std::endl;
                                    if(  (slipSystem->s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON
                                       && fabs(slipSystem->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
                                    {
                                        //                                        std::cout<<"here 4"<<std::endl;
                                        
                                        glissDeq.emplace_back(link.second->source,link.second->sink,gr->grainID,s);
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }
            
            size_t formedJunctions=0;
            for(const auto& tup : glissDeq)
            {
                const std::shared_ptr<NodeType>& source(std::get<0>(tup));
                const std::shared_ptr<NodeType>& sink(std::get<1>(tup));
                const size_t& sourceID(source->sID);
                const size_t& sinkID(sink->sID);
                const size_t& grainID(std::get<2>(tup));
                const size_t& slipID(std::get<3>(tup));
                
                const auto isLink(DN.link(sourceID,sinkID));
                if(isLink.first)
                {
                    
                    const VectorDim newNodeP(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()));
//                    const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
                    std::shared_ptr<NodeType> newNode(new NodeType(&DN,newNodeP,VectorDim::Zero(),1.0));
                    
                    std::vector<std::shared_ptr<NodeType>> loopNodes;
//                    std::vector<size_t> nodeIDs;
  
                    loopNodes.push_back(sink);      // insert in reverse order, sink first, source second
                    loopNodes.push_back(source);    // insert in reverse order, sink first, source second
                    loopNodes.push_back(newNode);

                    
//                    nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
//                    nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
//                    nodeIDs.push_back(newNodeID);
                    
                    LatticePlane glissilePlane(newNodeP,DN.poly.grain(grainID).slipSystems()[slipID]->n);
                    GlidePlaneKey<dim> glissilePlaneKey(grainID,glissilePlane);

//                    DN.insertLoop(nodeIDs,
//                                  DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                                  DN.glidePlaneFactory.get(glissilePlaneKey));

                    DN.insertLoop(loopNodes,
                                  DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
                                  DN.glidePlaneFactory.get(glissilePlaneKey));

                    
                    formedJunctions++;
                }
                
                
                
            }
            
//            DN.clearDanglingNodes();
            
            
//            std::cout<<"glissDeq.size="<<glissDeq.size()<<std::endl;
            model::cout<<"("<<formedJunctions<<" junctions)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            //            static_assert(0,"FINISH HERE");
        }
        
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
