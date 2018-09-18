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
#include <model/Geometry/SegmentSegmentDistance.h>
#include <model/Geometry/SweepPlane.h>

//#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>
#include <model/DislocationDynamics/DislocationNetworkRemesh.h>
#include <model/MPI/MPIcout.h>
#include <model/Threads/EqualIteratorRange.h>
#include <model/Threads/N2IteratorRange.h>

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
        
        
        //        typedef std::pair<std::pair<size_t,size_t>, double> EdgeIntersectionType;
        
        //        typedef std::pair<EdgeIntersectionType,EdgeIntersectionType> EdgeIntersectionPairType;
        //		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
        typedef std::deque<IntersectionType,Eigen::aligned_allocator<IntersectionType>> IntersectionTypeContainerType;
        
        
        
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
                bool isValidJunction(bndJunction || gbndJunction);
                if(!isValidJunction)
                {
                    
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
            model::cout<<"		Finding Junctions: sweep-line "<<std::flush;
            
            // Use SweepPlane to compute possible intersections
            SweepPlane<LinkType,dim> swp;
            for(const auto& link : DN.links())
            {
                if(!link.second->hasZeroBurgers())
                {
                    //                    swp.addSegment(link.second->source->get_P()(0),link.second->source->get_P()(1),*link.second);
                    swp.addSegment(link.second->source->get_P()(0),link.second->sink->get_P()(0),*link.second);
                }
            }
            swp.computeIntersectionPairs();
            model::cout<<"("<<swp.potentialIntersectionPairs().size()<<" pairs)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::flush;
            
            
            const auto t1= std::chrono::system_clock::now();
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
            
            
            model::cout<<" ("<<reducedIntersectionPairs.size()<<" reduced pairs). Pair-intersections  ("<<nThreads<<" threads) "<<std::flush;
            
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
                const bool linkAisGBnd(linkA->grainBoundaries().size());
                const bool linkBisGBnd(linkB->grainBoundaries().size());
                
                
                const bool bndJunction(   (linkAisBnd || linkBisBnd)
                                       &&  linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).norm()<FLT_EPSILON);
                
                const bool gbndJunction(   (linkAisGBnd || linkBisGBnd)
                                        &&  linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).norm()<FLT_EPSILON);
                
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
            model::cout<<nIntersections<<" physical intersections. "<<std::flush;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        
        
        
        /**********************************************************************/
        size_t junctionStep(const double& dx)
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
                        
                        auto  Ni= t<0.5? L1.second->source : L1.second->sink;
                        bool NiIsNew=false;
                        if((ssd.x0-(t<0.5? L1.second->source->get_P() : L1.second->sink->get_P())).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                        {
                            if((ssd.x0-(t<0.5? L1.second->sink->get_P() : L1.second->source->get_P())).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                            {
                                Ni=DN.expand(key1.first,key1.second,ssd.x0);
                                NiIsNew=true;
                            }
                            else
                            {
                                Ni=t<0.5? L1.second->sink : L1.second->source;
                            }
                        }
                        
                        auto Nj=u<0.5? L2.second->source : L2.second->sink;
                        bool NjIsNew=false;
                        if((ssd.x1-(u<0.5? L2.second->source->get_P() : L2.second->sink->get_P())).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                        {
                            if((ssd.x1-(u<0.5? L2.second->sink->get_P() : L2.second->source->get_P())).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                            {
                                Nj=DN.expand(key2.first,key2.second,ssd.x1);
                                NjIsNew=true;
                            }
                            else
                            {
                                Nj=u<0.5? L2.second->sink : L2.second->source;;
                            }
                        }
                        
                        VerboseJunctions(1,"forming junction "<<key1.first<<"->"<<key1.second<<", "
                                         /*                   */ <<key2.first<<"->"<<key2.second<<", "
                                         /*                   */ <<"dMin="<<ssd.dMin<<", "
                                         /*                   */ <<"@ ("<<t<<","<<u<<"), "
                                         /*                   */ <<"contracting "<<Ni->sID<<" "<<Nj->sID<<std::endl;);
                        
                        if(Ni->sID!=Nj->sID)
                        {
                            const bool success=DN.contract(Ni,Nj);
                            nContracted+=success;
                            if(!success)
                            {
                                if(NiIsNew)
                                {
                                    DN.remove(Ni->sID);
                                }
                                if(NjIsNew)
                                {
                                    DN.remove(Nj->sID);
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
            
            std::deque<std::tuple<size_t,size_t,size_t,size_t>> glissDeq;
            
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
                                    
                                    //                                    std::cout<<"here 3 "<<"\n"<<slipSystem.s.cartesian().transpose()<<"\n"<<link.second->burgers().transpose()<<std::endl;
                                    //                                    std::cout<<((slipSystem.s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON)<<std::endl;
                                    //                                    std::cout<<(fabs(slipSystem.n.cartesian().normalized().dot(unitChord)))<<std::endl;
                                    if(  (slipSystem.s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON
                                       && fabs(slipSystem.n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
                                    {
                                        //                                        std::cout<<"here 4"<<std::endl;
                                        
                                        glissDeq.emplace_back(link.second->source->sID,link.second->sink->sID,gr->grainID,s);
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }
            
            for(const auto& tup : glissDeq)
            {
                const size_t& sourceID(std::get<0>(tup));
                const size_t& sinkID(std::get<1>(tup));
                const size_t& grainID(std::get<2>(tup));
                const size_t& slipID(std::get<3>(tup));
                
                const auto isLink(DN.link(sourceID,sinkID));
                if(isLink.first)
                {
                    
                    const VectorDim newNodeP(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()));
                    const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
                    
                    std::vector<size_t> nodeIDs;
                    
                    nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
                    nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
                    nodeIDs.push_back(newNodeID);
                    
                    DN.insertLoop(nodeIDs,
                                  DN.poly.grain(grainID).slipSystems()[slipID].s.cartesian(),
                                  DN.poly.grain(grainID).slipSystems()[slipID].n.cartesian(),
                                  newNodeP,
                                  grainID);
                }
                
                
                
            }
            
            DN.clearDanglingNodes();
            
            
            std::cout<<"glissDeq.size="<<glissDeq.size()<<std::endl;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            //            static_assert(0,"FINISH HERE");
        }
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        
        static double collisionTol;     //! The tolerance (in units of distance) used for collision detection
        static int verboseJunctions;
        const double infiniteLineLength;
        //        static size_t maxJunctionIterations;
        
        /**********************************************************************/
        DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        /* init list */,infiniteLineLength(10000.0)
        {
            
        }
        
        /**********************************************************************/
        void formJunctions(const bool& maxJunctionIterations, const double& dx)
        {
            size_t nContracted=1;
            size_t iterations=0;
            while(nContracted && iterations<maxJunctionIterations)
            {
                nContracted=junctionStep(dx);
                glissileJunctions(dx);
                iterations++;
            }
        }
        
    };
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
    
    template <typename DislocationNetworkType>
    int DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions=0;
    
}
#endif





//        /**********************************************************************/
//        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
//                               const size_t& nThreads)
//
//        {/*! @param[in]  avoidNodeIntersection
//          *  Computes all the intersections between the edges of the DislocationNetwork
//          */
//
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"		Finding Junctions ("<<nThreads<<" threads)... "<<std::flush;
//
//            N2IteratorRange<typename NetworkLinkContainerType::const_iterator> eir(DN.links().begin(),DN.links().end(),nThreads);
//
//            //! 2- loop over all links and determine their intersections
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//            for (size_t thread=0;thread<eir.size();thread++)
//            {
//                for (typename NetworkLinkContainerType::const_iterator linkIterA=eir[thread].first;linkIterA!=eir[thread].second;linkIterA++)
//                {
//                    if(!linkIterA->second->hasZeroBurgers()) // don't loop over zero-Burgers segments
//                    {
//                        for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.links().end();linkIterB++)
//                        {
//                            if (   linkIterA->second->sID!=linkIterB->second->sID           // don't intersect with itself
//                                && !linkIterB->second->hasZeroBurgers()   // don't intersect with zero-Burgers segments
//                                )
//                            {
//
//                                const bool linkAisBnd(linkIterA->second->isBoundarySegment());
//                                const bool linkBisBnd(linkIterB->second->isBoundarySegment());
//                                const bool linkAisGBnd(linkIterA->second->grainBoundaries().size());
//                                const bool linkBisGBnd(linkIterB->second->grainBoundaries().size());
//
//                                const bool frankRule(linkIterA->second->burgers().dot(linkIterB->second->burgers())*linkIterA->second->chord().dot(linkIterB->second->chord())<=0.0);
//
//                                const bool bndJunction(   (linkAisBnd || linkBisBnd)
//                                                       &&  linkIterA->second->glidePlaneNormal().cross(linkIterB->second->glidePlaneNormal()).norm()<FLT_EPSILON);
//
//                                const bool gbndJunction(   (linkAisGBnd || linkBisGBnd)
//                                                       &&  linkIterA->second->glidePlaneNormal().cross(linkIterB->second->glidePlaneNormal()).norm()<FLT_EPSILON);
//
//
//                                const bool isValidJunction(   (frankRule   && !linkAisBnd && !linkBisBnd && !linkAisGBnd && !linkBisGBnd)  // energy rule is satisfied for internal segments
//                                                           || bndJunction   // junciton between parallel boundary segments
//                                                           || gbndJunction  // junciton between parallel grain-boundary segments
//                                                           );
//
//                                if(isValidJunction)
//                                {
//
//
//                                    double currentcCollisionTOL=collisionTol;
//                                    if(   linkIterA->second->glidePlaneNormal().squaredNorm()>FLT_EPSILON
//                                       && linkIterB->second->glidePlaneNormal().squaredNorm()>FLT_EPSILON
//                                       && linkIterA->second->glidePlaneNormal().cross(linkIterB->second->glidePlaneNormal()).squaredNorm()<FLT_EPSILON)
//                                    {// segments on parallel or coincident planes, reduce tolerance
//                                        currentcCollisionTOL=FLT_EPSILON;
//                                    }
//
//                                    const bool intersectionIsSourceSource(linkIterA->second->source->sID==linkIterB->second->source->sID);
//                                    const bool intersectionIsSourceSink(linkIterA->second->source->sID==linkIterB->second->sink->sID);
//                                    const bool intersectionIsSinkSource(linkIterA->second->sink->sID==linkIterB->second->source->sID);
//                                    const bool intersectionIsSinkSink(linkIterA->second->sink->sID==linkIterB->second->sink->sID);
//
//                                    if(intersectionIsSourceSource)
//                                    {
//                                        if(!linkIterA->second->source->isSimple())
//                                        {// non-simple common node, increase tolerance
//                                            currentcCollisionTOL=0.5*collisionTol;
//                                        }
//
//                                        // intersect sink of A with link B
//                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->sink->get_P(), // fake degenerete segment at sink of A
//                                                                         linkIterA->second->sink->get_P(),  // fake degenerete segment at sink of A
//                                                                         linkIterB->second->source->get_P(),
//                                                                         linkIterB->second->sink->get_P(),
//                                                                         1.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
//
//                                        // intersect sink of B with link A
//                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
//                                                                         linkIterA->second->sink->get_P(),
//                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
//                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
//                                                                         1.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
//
//                                    }
//                                    else if(intersectionIsSourceSink)
//                                    {
//                                        if(!linkIterA->second->source->isSimple())
//                                        {// non-simple common node, increase tolerance
//                                            currentcCollisionTOL=0.5*collisionTol;
//                                        }
//
//                                        // intersect sink of A with link B
//                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->sink->get_P(), // fake degenerete segment at sink of A
//                                                                         linkIterA->second->sink->get_P(),  // fake degenerete segment at sink of A
//                                                                         linkIterB->second->source->get_P(),
//                                                                         linkIterB->second->sink->get_P(),
//                                                                         1.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
//
//                                        // intersect source of B with link A
//                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
//                                                                         linkIterA->second->sink->get_P(),
//                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
//                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
//                                                                         0.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
//
//
//                                    }
//                                    else if(intersectionIsSinkSource)
//                                    {
//                                        if(!linkIterA->second->sink->isSimple())
//                                        {// non-simple common node, increase tolerance
//                                            currentcCollisionTOL=0.5*collisionTol;
//                                        }
//
//                                        // intersect source of A with link B
//                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->source->get_P(), // fake degenerete segment at source of A
//                                                                         linkIterA->second->source->get_P(),  // fake degenerete segment at source of A
//                                                                         linkIterB->second->source->get_P(),
//                                                                         linkIterB->second->sink->get_P(),
//                                                                         0.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
//
//
//                                        // intersect sink of B with link A
//                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
//                                                                         linkIterA->second->sink->get_P(),
//                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
//                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
//                                                                         1.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
//
//
//                                    }
//                                    else if(intersectionIsSinkSink)
//                                    {
//
//                                        if(!linkIterA->second->sink->isSimple())
//                                        {// non-simple common node, increase tolerance
//                                            currentcCollisionTOL=0.5*collisionTol;
//                                        }
//
//                                        // intersect source of A with link B
//                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->source->get_P(), // fake degenerete segment at source of A
//                                                                         linkIterA->second->source->get_P(),  // fake degenerete segment at source of A
//                                                                         linkIterB->second->source->get_P(),
//                                                                         linkIterB->second->sink->get_P(),
//                                                                         0.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
//
//
//                                        // intersect source of B with link A
//                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
//                                                                         linkIterA->second->sink->get_P(),
//                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
//                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
//                                                                         0.0);
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
//
//                                    }
//                                    else
//                                    {// segments don't share vertices
//
//                                        SegmentSegmentDistance<dim> ssd(linkIterA->second->source->get_P(),
//                                                                        linkIterA->second->sink->get_P(),
//                                                                        linkIterB->second->source->get_P(),
//                                                                        linkIterB->second->sink->get_P()
//                                                                        );
//                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssd,currentcCollisionTOL);
//                                    }
//                                }
//                            } // end don't self-intersect
//                        }
//                    }
//                } // end loop over first segment
//            }// end loop ever threads
//
//            int nIntersections=0;
//            for (const auto& intersectionByThreadContainer : intersectionContainer)
//            {
//                nIntersections+=intersectionByThreadContainer.size();
//            }
//            model::cout<<nIntersections<<" physical intersections. "<<std::flush;
//            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//
//        }

//            /**********************************************************************/
//            void breakZeroLengthJunctions()
//            {
//                static_assert(0,"THIS IS DDEPRECATED, USE NEW NODE-REMOVE FUNCTIONALITY INSTEAD");
//                const auto t0= std::chrono::system_clock::now();
//                model::cout<<"		Breaking zero-length Junctions... "<<std::flush;
//
//                std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecompFirst;
//                std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecompSecond;
//
//
//                // TO DO: PARALLELIZE THIS LOOP
//                for (auto& nodePair : DN.nodes())
//                {
//                    const auto temp=nodePair.second.edgeDecomposition();
//                    //std::deque<std::pair<size_t,size_t> > temp=nodePair.second.edgeDecomposition();
//
//                    if(temp.first.size() && temp.second.size())
//                    {
//                        nodeDecompFirst.emplace_back(nodePair.second.sID,temp.first);
//                        nodeDecompSecond.emplace_back(nodePair.second.sID,temp.second);
//                    }
//                }
//
//                int broken=0;
//
//                for(size_t n=0;n<nodeDecompFirst.size();++n)
//                {
//                    const size_t& i=nodeDecompFirst[n].first;
//                    auto Ni=DN.node(i);
//                    assert(Ni.first);
//                    //                size_t m=DN.insertVertex(Ni.second->get_P());
//
//                    // Check that Links still exist
//                    //std::deque<VectorDimD> linkDirs;
//
//                    bool linksFirstexist=true;
//                    VectorDimD avrFirst=VectorDimD::Zero();
//                    for(size_t d=0;d<nodeDecompFirst[n].second.size();++d)
//                    {
//                        const size_t& j=nodeDecompFirst[n].second[d].first;
//                        const size_t& k=nodeDecompFirst[n].second[d].second;
//                        if(i==j)
//                        {
//                            auto Lik=DN.link(i,k);
//                            if(Lik.first)
//                            {
//                                avrFirst+=Lik.second->chord().normalized();
//                            }
//                            linksFirstexist*=Lik.first;
//                        }
//                        else if(i==k)
//                        {
//                            auto Lji=DN.link(j,i);
//                            if(Lji.first)
//                            {
//                                avrFirst-=Lji.second->chord().normalized();
//                            }
//                            linksFirstexist*=Lji.first;
//                        }
//                        else
//                        {
//                            assert(0 && "i must be equal to either j or k.");
//                        }
//                    }
//                    const double avrFirstNorm=avrFirst.norm();
//                    if(avrFirstNorm>FLT_EPSILON)
//                    {
//                        avrFirst/=avrFirstNorm;
//                    }
//
//                    bool linksSecondexist=true;
//                    VectorDimD avrSecond=VectorDimD::Zero();
//                    for(size_t d=0;d<nodeDecompSecond[n].second.size();++d)
//                    {
//                        const size_t& j=nodeDecompSecond[n].second[d].first;
//                        const size_t& k=nodeDecompSecond[n].second[d].second;
//                        if(i==j)
//                        {
//                            auto Lik=DN.link(i,k);
//                            if(Lik.first)
//                            {
//                                avrSecond+=Lik.second->chord().normalized();
//                            }
//                            linksSecondexist*=Lik.first;
//                        }
//                        else if(i==k)
//                        {
//                            auto Lji=DN.link(j,i);
//                            if(Lji.first)
//                            {
//                                avrSecond-=Lji.second->chord().normalized();
//                            }
//                            linksSecondexist*=Lji.first;
//                        }
//                        else
//                        {
//                            assert(0 && "i must be equal to either j or k.");
//                        }
//                    }
//                    const double avrSecondNorm=avrSecond.norm();
//                    if(avrSecondNorm>FLT_EPSILON)
//                    {
//                        avrSecond/=avrSecondNorm;
//                    }
//
//
//                    if(linksFirstexist && linksSecondexist && avrSecond.dot(avrFirst)<-0.0)
//                    {
//                        std::cout<<"NodeBreaking "<<Ni.second->sID<<" "<<avrSecond.dot(avrFirst)<<std::endl;
//
//                        size_t m=DN.insertVertex(Ni.second->get_P(),Ni.second->grain.grainID).first->first;
//
//                        for(size_t d=0;d<nodeDecompFirst[n].second.size();++d)
//                        {
//                            const size_t& j=nodeDecompFirst[n].second[d].first;
//                            const size_t& k=nodeDecompFirst[n].second[d].second;
//                            if(i==j)
//                            {
//                                auto Lik=DN.link(i,k);
//                                assert(Lik.first);
//                                DN.connect(m,k,Lik.second->flow);
//                                DN.template disconnect<0>(i,k);
//                            }
//                            else if(i==k)
//                            {
//                                auto Lji=DN.link(j,i);
//                                assert(Lji.first);
//                                DN.connect(j,m,Lji.second->flow);
//                                DN.template disconnect<0>(j,i);
//                            }
//                            else
//                            {
//                                assert(0 && "i must be equal to either j or k.");
//                            }
//                        }
//
//
//                        broken++;
//                    }
//                }
//                model::cout<<broken<<" broken."<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }

