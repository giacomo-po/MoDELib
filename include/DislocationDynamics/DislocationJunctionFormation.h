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
#include <TypeTraits.h>
#include <StressStraight.h>

#include <DislocationNetworkRemesh.h>
//#include <MPIcout.h>
// #include <EqualIteratorRange.h>
// #include <N2IteratorRange.h>

#ifndef NDEBUG
#define VerboseJunctions(N,x) if(verboseJunctions>=N){std::cout<<x;}
#else
#define VerboseJunctions(N,x)
#endif


namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationJunctionFormation
    {
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        typedef typename DislocationNetworkType::NetworkLinkType NetworkLinkType;
        typedef typename DislocationNetworkType::NetworkNodeType NetworkNodeType;
        typedef typename DislocationNetworkType::LoopNodeType LoopNodeType;

        typedef typename DislocationNetworkType::LoopType LoopType;
//        typedef typename DislocationNetworkType::IsNetworkEdgeType IsNetworkNetworkLinkType;
//        typedef typename DislocationNetworkType::IsNetworkNodeType IsNetworkNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename NetworkLinkType::KeyType KeyType;
        
        typedef std::tuple<KeyType,KeyType,SegmentSegmentDistance<dim>> IntersectionType;
        typedef std::deque<IntersectionType> IntersectionTypeContainerType;
        
//Original Working Version
//         void insertIntersection(std::deque<IntersectionTypeContainerType> &intersectionContainer,
//                                 const NetworkLinkType *const linkA,
//                                 const NetworkLinkType *const linkB,
//                                 const SegmentSegmentDistance<dim> &ssd,
//                                 const double &currentcCollisionTOL,
//                                 const bool &bndJunction, const bool &gbndJunction)
//         {
//             VerboseJunctions(2, "insertIntersection " << linkA->tag() <<" [ "<<linkA->source->get_P().transpose()<<", "<<linkA->sink->get_P().transpose()
//             <<" ]" << "," << linkB->tag() <<" [ "<<linkB->source->get_P().transpose()<<", "<<linkB->sink->get_P().transpose()
//             <<" ]" << std::endl;);
//             VerboseJunctions(2, "insertIntersection "
//                                     << "ssd.dMin=" << ssd.dMin << std::endl;);
//             VerboseJunctions(2, "insertIntersection "
//                                     << "currentcCollisionTOL=" << currentcCollisionTOL << std::endl;);

//             if (ssd.dMin < currentcCollisionTOL)
//             {
//                 bool isValidJunction((bndJunction || gbndJunction) && DN.simulationParameters.simulationType != 2);
//                 if (!isValidJunction && !linkA->isBoundarySegment() && !linkB->isBoundarySegment() && !linkA->isGrainBoundarySegment() && !linkB->isGrainBoundarySegment())
//                 { // Check force condition for internal segments

//                     const VectorDim chordA(linkA->sink->get_P() - linkA->source->get_P());
//                     const double LA(chordA.norm());
//                     const VectorDim chordB(linkB->sink->get_P() - linkB->source->get_P());
//                     const double LB(chordB.norm());

//                     if (LA > FLT_EPSILON && LB > FLT_EPSILON)
//                     {

//                         StressStraight<dim> stressA(ssd.x0 - infiniteLineLength / LA * chordA,
//                                                     ssd.x0 + infiniteLineLength / LA * chordA,
//                                                     linkA->burgers());

//                         StressStraight<dim> stressB(ssd.x1 - infiniteLineLength / LB * chordB,
//                                                     ssd.x1 + infiniteLineLength / LB * chordB,
//                                                     linkB->burgers());

//                         if (ssd.dMin > FLT_EPSILON)
//                         {
//                             VerboseJunctions(3, "Non-intersecting pair" << std::endl;);
//                             const VectorDim forceOnA = (stressB.stress(ssd.x0) * linkA->burgers()).cross(chordA);
//                             const VectorDim forceOnB = (stressA.stress(ssd.x1) * linkB->burgers()).cross(chordB);
//                             const VectorDim dxShift(ssd.x1 - ssd.x0);
//                             if (forceOnA.dot(dxShift) > FLT_EPSILON && forceOnB.dot(dxShift) < -FLT_EPSILON)
//                             {
//                                 VerboseJunctions(3, "attractive pair 1" << std::endl;);
//                                 isValidJunction = true; // for non-parallel lines this neglects the energy of rotation
//                             }
//                             else
//                             {
//                                 VerboseJunctions(3, "non-attractive pair" << std::endl;);
//                             }
//                         }
//                         else
//                         {
//                             VerboseJunctions(3, "Intersecting pair" << std::endl;);

//                             const VectorDim dXA(linkA->source->get_V() + linkA->sink->get_V());
//                             const double dXAnorm(dXA.norm());
//                             const VectorDim dXB(linkB->source->get_V() + linkB->sink->get_V());
//                             const double dXBnorm(dXB.norm());

//                             if (dXAnorm > FLT_EPSILON && dXBnorm > FLT_EPSILON)
//                             {
//                                 const VectorDim x0Shift(ssd.x0 - dXA / dXAnorm);
//                                 const VectorDim x1Shift(ssd.x1 - dXB / dXBnorm);
//                                 const VectorDim forceOnA = (stressB.stress(x0Shift) * linkA->burgers()).cross(chordA);
//                                 const VectorDim forceOnB = (stressA.stress(x1Shift) * linkB->burgers()).cross(chordB);
//                                 const VectorDim dxShift(x1Shift - x0Shift);
//                                 if (forceOnA.dot(dxShift) > FLT_EPSILON && forceOnB.dot(dxShift) < -FLT_EPSILON)
//                                 {
//                                     VerboseJunctions(3, "attractive pair 2" << std::endl;);
//                                     isValidJunction = true; // for non-parallel lines this neglects the energy of rotation
//                                 }
//                             }
//                             else if (dXAnorm > FLT_EPSILON && dXBnorm <= FLT_EPSILON)
//                             {
//                                 const VectorDim x0Shift(ssd.x0 - dXA / dXAnorm);
//                                 const VectorDim forceOnA = (stressB.stress(x0Shift) * linkA->burgers()).cross(chordA);
//                                 const VectorDim forceOnB = (stressA.stress(ssd.x1) * linkB->burgers()).cross(chordB);
//                                 const VectorDim dxShift(ssd.x1 - x0Shift);
//                                 if (forceOnA.dot(dxShift) > FLT_EPSILON && forceOnB.dot(dxShift) < -FLT_EPSILON)
//                                 {
//                                     VerboseJunctions(3, "attractive pair 3" << std::endl;);
//                                     isValidJunction = true; // for non-parallel lines this neglects the energy of rotation
//                                 }
//                             }
//                             else if (dXAnorm <= FLT_EPSILON && dXBnorm > FLT_EPSILON)
//                             {
//                                 const VectorDim x1Shift(ssd.x1 - dXB / dXBnorm);
//                                 const VectorDim forceOnA = (stressB.stress(ssd.x0) * linkA->burgers()).cross(chordA);
//                                 const VectorDim forceOnB = (stressA.stress(x1Shift) * linkB->burgers()).cross(chordB);
//                                 const VectorDim dxShift(x1Shift - ssd.x0);
//                                 if (forceOnA.dot(dxShift) > FLT_EPSILON && forceOnB.dot(dxShift) < -FLT_EPSILON)
//                                 {
//                                     VerboseJunctions(3, "attractive pair 4" << std::endl;);
//                                     isValidJunction = true; // for non-parallel lines this neglects the energy of rotation
//                                 }
//                             }
//                             else
//                             { // cannot determine seprate points
//                                 VerboseJunctions(3, "cannot determine seprate points" << std::endl;);
//                             }

//                             // //Coding the frank rule for this criteria
//                             // const double b1((linkA->burgers()).squaredNorm());
//                             // const double b2((linkB->burgers()).squaredNorm());
//                             // const double bJunction((linkA->burgers()+linkB->burgers()).squaredNorm());

//                             // isValidJunction= ((bJunction)<=(b1+b2+FLT_EPSILON));
//                             // VerboseJunctions(3, " From Frank's Rule "<<isValidJunction << std::endl;);
//                         }
//                     }
//                 }

//                 if (isValidJunction)
//                 {
// #ifdef _OPENMP
//                     intersectionContainer[omp_get_thread_num()].emplace_back(linkA->key,
//                                                                              linkB->key,
//                                                                              ssd);
// #else
//                     intersectionContainer[0].emplace_back(linkA->key,
//                                                           linkB->key,
//                                                           ssd);
// #endif
//                 }
//             }
//             else
//             {
//                 VerboseJunctions(3, "dMin=" << ssd.dMin << ", collisionTol=" << currentcCollisionTOL << std::endl;);
//             }
//         }

//Testing frank rule for boundary coontraction
        void insertIntersection(std::deque<IntersectionTypeContainerType> &intersectionContainer,
                                const NetworkLinkType *const linkA,
                                const NetworkLinkType *const linkB,
                                const SegmentSegmentDistance<dim> &ssd,
                                const double &currentcCollisionTOL,
                                const bool &bndJunction, const bool &gbndJunction)
        {
            VerboseJunctions(2, "insertIntersection " << linkA->tag() << " [ " << linkA->source->get_P().transpose() << ", " << linkA->sink->get_P().transpose()
                                                      << " ]"
                                                      << "," << linkB->tag() << " [ " << linkB->source->get_P().transpose() << ", " << linkB->sink->get_P().transpose()
                                                      << " ]" << std::endl;);
            VerboseJunctions(2, "insertIntersection "
                                    << "ssd.dMin=" << ssd.dMin << std::endl;);
            VerboseJunctions(2, "insertIntersection "
                                    << "currentcCollisionTOL=" << currentcCollisionTOL << std::endl;);

            if (ssd.dMin < currentcCollisionTOL)
            {
                bool isValidJunction((bndJunction || gbndJunction) && DN.simulationParameters.simulationType != 2);
                if (!isValidJunction && !linkA->isBoundarySegment() && !linkB->isBoundarySegment() && !linkA->isGrainBoundarySegment() && !linkB->isGrainBoundarySegment())
                { // Check force condition for internal segments

                    const VectorDim chordA(linkA->sink->get_P() - linkA->source->get_P());
                    const double LA(chordA.norm());
                    const VectorDim chordB(linkB->sink->get_P() - linkB->source->get_P());
                    const double LB(chordB.norm());

                    if (LA > FLT_EPSILON && LB > FLT_EPSILON)
                    {

                        StressStraight<dim> stressA(DN.poly,ssd.x0 - infiniteLineLength / LA * chordA,
                                                    ssd.x0 + infiniteLineLength / LA * chordA,
                                                    linkA->burgers());

                        StressStraight<dim> stressB(DN.poly,ssd.x1 - infiniteLineLength / LB * chordB,
                                                    ssd.x1 + infiniteLineLength / LB * chordB,
                                                    linkB->burgers());

                        if (ssd.dMin > FLT_EPSILON)
                        {
                            VerboseJunctions(3, "Non-intersecting pair" << std::endl;);
                            const VectorDim forceOnA = (stressB.stress(ssd.x0) * linkA->burgers()).cross(chordA);
                            const VectorDim forceOnB = (stressA.stress(ssd.x1) * linkB->burgers()).cross(chordB);
                            const VectorDim dxShift(ssd.x1 - ssd.x0);
                            if (forceOnA.dot(dxShift) > FLT_EPSILON && forceOnB.dot(dxShift) < -FLT_EPSILON)
                            {
                                VerboseJunctions(3, "attractive pair 1" << std::endl;);
                                isValidJunction = true; // for non-parallel lines this neglects the energy of rotation
                            }
                            else
                            {
                                VerboseJunctions(3, "non-attractive pair" << std::endl;);
                            }
                        }
                        else
                        {
                            VerboseJunctions(3, "Intersecting pair .. Determining via frank rule" << std::endl;);
                            // const bool frankRule(linkIterA->second->burgers().dot(linkIterB->second->burgers())*linkIterA->second->chord().dot(linkIterB->second->chord())<=0.0);
                            isValidJunction=(linkA->burgers().dot(linkB->burgers())*linkA->chord().dot(linkB->chord())<=0.0);
                            // //Coding the frank rule for this criteria
                            // const double b1((linkA->burgers()).squaredNorm());
                            // const double b2((linkB->burgers()).squaredNorm());
                            // const double bJunction((linkA->burgers()+linkB->burgers()).squaredNorm());

                            // isValidJunction= ((bJunction)<=(b1+b2+FLT_EPSILON));
                            // VerboseJunctions(3, " From Frank's Rule "<<isValidJunction << std::endl;);
                        }
                    }
                }

                if (isValidJunction)
                {
#ifdef _OPENMP
                    intersectionContainer[omp_get_thread_num()].emplace_back(linkA->key,
                                                                             linkB->key,
                                                                             ssd);
#else
                    intersectionContainer[0].emplace_back(linkA->key,
                                                          linkB->key,
                                                          ssd);
#endif
                }
            }
            else
            {
                VerboseJunctions(3, "dMin=" << ssd.dMin << ", collisionTol=" << currentcCollisionTOL << std::endl;);
            }
        }

        /**********************************************************************/
        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                               const size_t& nThreads)
        
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
            
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"		Finding collisions "<<std::flush;
            
            // Use SweepPlane to compute possible intersections
            SweepPlane<NetworkLinkType,dim> swp;
            for(const auto& link : DN.networkLinks())
            {
                const auto sharedLink(link.second.lock());
                
                if(   (!sharedLink->hasZeroBurgers() && !sharedLink->isVirtualBoundarySegment())
                   //                   || (link.second->isBoundarySegment() && DN.useVirtualExternalLoops))
                   || sharedLink->isBoundarySegment() )
                    
                {
                    //                    swp.addSegment(link.second->source->get_P()(0),link.second->source->get_P()(1),*link.second);
                    swp.addSegment(sharedLink->source->get_P()(0),sharedLink->sink->get_P()(0),*sharedLink);
                }
            }
            swp.computeIntersectionPairs();
            std::cout<<"("<<swp.potentialIntersectionPairs().size()<<" sweep-line pairs) "<<defaultColor<<std::flush;
            
            
            std::deque<std::pair<const NetworkLinkType*,const NetworkLinkType*>> reducedIntersectionPairs;
            for(size_t k=0;k<swp.potentialIntersectionPairs().size();++k)
            {
                const auto& linkA(swp.potentialIntersectionPair(k).first);
                const auto& linkB(swp.potentialIntersectionPair(k).second);
                


                const VectorDim& sourceA(linkA->source->get_P());
                const VectorDim&   sinkA(linkA->sink->get_P());
                const VectorDim& sourceB(linkB->source->get_P());
                const VectorDim&   sinkB(linkB->sink->get_P());

                // std::cout << "LinkA->Tag() " << linkA->tag() << " LinkB->Tag() " << linkB->tag() << std::endl;
                // std::cout << " SourceA=>sinkA " << sourceA.transpose() << " => " << sinkA.transpose() << std::endl;
                // std::cout << " SourceB=>sinkB " << sourceB.transpose() << " => " << sinkB.transpose() << std::endl;
                // std::cout << " Collision Tolerance " << collisionTol << std::endl;

                // SweepPlane sweeps along x direction, so now check possible overlap in y and z direction
                bool yDontOverlap((std::min(sourceA(1),sinkA(1))>std::max(sourceB(1),sinkB(1))+collisionTol) || (std::max(sourceA(1),sinkA(1))<std::min(sourceB(1),sinkB(1))-collisionTol));
                bool zDontOverlap((std::min(sourceA(2),sinkA(2))>std::max(sourceB(2),sinkB(2))+collisionTol) || (std::max(sourceA(2),sinkA(2))<std::min(sourceB(2),sinkB(2))-collisionTol));
                
                // std::cout<<"YdontoverLap ==> Zdontoverlap"<<yDontOverlap<<" ==> "<<zDontOverlap<<std::endl;

                if(!yDontOverlap && !zDontOverlap)
                {
                    reducedIntersectionPairs.emplace_back(linkA,linkB);
                }
            }
            
            
            std::cout<<" ("<<reducedIntersectionPairs.size()<<" reduced pairs) "<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            const auto t1= std::chrono::system_clock::now();
            std::cout<<"        Selecting junctions ("<<nThreads<<" threads): "<<std::flush;
            
            //! 2- loop over all links and determine their intersections
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
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
//                    if(!linkA->source->isSimple())
//                    {// non-simple common node, increase tolerance
//                        currentcCollisionTOL=0.5*collisionTol;
//                    }
                    
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
//                    if(!linkA->source->isSimple())
//                    {// non-simple common node, increase tolerance
//                        currentcCollisionTOL=0.5*collisionTol;
//                    }
                    
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
//                    if(!linkA->sink->isSimple())
//                    {// non-simple common node, increase tolerance
//                        currentcCollisionTOL=0.5*collisionTol;
//                    }
                    
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
                    
//                    if(!linkA->sink->isSimple())
//                    {// non-simple common node, increase tolerance
//                        currentcCollisionTOL=0.5*collisionTol;
//                    }
                    
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
            std::cout<<nIntersections<<" physical junctions "<<std::flush;
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        
    /**********************************************************************/
        //     Yash's version(Working)
        std::pair<std::shared_ptr<NetworkNodeType>, bool> junctionNode(const double &t,
                                                                       const VectorDim &x,
                                                                       const std::shared_ptr<NetworkLinkType> &L)
        {
//            VerboseJunctions(4, "JunctionNode for segment: " << key.first << "->" << key.second << " @ " << t << std::endl;);

            auto Nclose = t < 0.5 ? L->source : L->sink;
            auto Nfar = t >= 0.5 ? L->source : L->sink;
            VerboseJunctions(4, "Nclose=" << Nclose->sID << " (GP Size)=> " << Nclose->glidePlanes().size() << std::endl;);
            VerboseJunctions(4, "Nfar=" << Nfar->sID << " (GP Size)=> " << Nfar->glidePlanes().size() << std::endl;);
            //            if(t>FLT_EPSILON && t<1.0-FLT_EPSILON && !Nclose->isBoundaryNode() && !Nfar->isBoundaryNode()) //This condition is the problem It needs to be enabled
            if (t > FLT_EPSILON && t < 1.0 - FLT_EPSILON)
            { // intersection point is not an end node
                if ((Nclose->get_P() - x).norm() > DN.networkRemesher.Lmin)
                { // far enough from Nclose
                    VerboseJunctions(4, "JunctionNode case a" << std::endl;);

                    const VectorDim snappedPosition(L->glidePlanes().size() > 0 ? L->snapToGlidePlanes(x) : x);
                    const auto newNetNode(DN.networkNodes().create(snappedPosition, L->source->get_V() * (1.0 - t) + L->sink->get_V() * t, L->source->velocityReduction() * (1.0 - t) + L->sink->velocityReduction() * t));
                    DN.expandNetworkLink(L, newNetNode);
                    return std::make_pair(newNetNode, true);
                }
                else
                { // close to from Nclose
                    if (Nclose->isMovableTo(x) && !Nclose->isBoundaryNode())
                    {
                        VerboseJunctions(4, "JunctionNode case b" << std::endl;);
                        return std::make_pair(Nclose, false);
                    }
                    else
                    {
                        if ((Nfar->get_P() - x).norm() > DN.networkRemesher.Lmin)
                        {
                            VerboseJunctions(4, "JunctionNode case c" << std::endl;);

                            const VectorDim snappedPosition(L->glidePlanes().size() > 0 ? L->snapToGlidePlanes(x) : x);

                            const auto newNetNode(DN.networkNodes().create(snappedPosition, L->source->get_V() * (1.0 - t) + L->sink->get_V() * t, L->source->velocityReduction() * (1.0 - t) + L->sink->velocityReduction() * t));
                            // std::cout << "x is " << x.transpose() << std::endl;

                            // std::cout<<"Expanding network Links c "<<L->tag()<<" [ "<<L->source->get_P().transpose()<<" -> "<<L->sink->get_P().transpose()<<" ] GlidePlane size "
                            // <<" [ "<<L->source->glidePlanes().size()<<" --> "<<L->sink->glidePlanes().size() <<" ] "
                            // <<L->glidePlanes().size()<<" loopLinks size C"<<L->loopLinks().size()<<std::endl;

                            DN.expandNetworkLink(L, newNetNode);
                            // std::cout<<"Expanded network Links"<<std::endl;

                            return std::make_pair(newNetNode, true);
                        }
                        else
                        {
                            if (Nfar->isMovableTo(x) && !Nfar->isBoundaryNode())
                            {
                                VerboseJunctions(4, "JunctionNode case d" << std::endl;);
                                return std::make_pair(Nfar, false);
                            }
                            else
                            {
                                VerboseJunctions(4, "JunctionNode case e" << std::endl;);
                                // std::cout<<"x is "<<x.transpose()<<std::endl;
                                //   std::cout<<"Expanding network Links E "<<L->tag()<<" [ "<<L->source->get_P().transpose()<<" -> "<<L->sink->get_P().transpose()<<" ] GlidePlane size "
                                // <<" [ "<<L->source->glidePlanes().size()<<" --> "<<L->sink->glidePlanes().size() <<" ] "
                                // <<L->glidePlanes().size()<<" loopLinks size E"<<L->loopLinks().size()<<std::endl;

                                const VectorDim snappedPosition(L->glidePlanes().size() > 0 ? L->snapToGlidePlanes(x) : x);

                                const auto newNetNode(DN.networkNodes().create(snappedPosition, L->source->get_V() * (1.0 - t) + L->sink->get_V() * t, L->source->velocityReduction() * (1.0 - t) + L->sink->velocityReduction() * t));
                                DN.expandNetworkLink(L, newNetNode);
                                return std::make_pair(newNetNode, true);
                            }
                        }
                    }
                }
            }
            else
            {// junction node is an existing node
                if (!Nclose->isBoundaryNode())
                {
                    VerboseJunctions(4, "JunctionNode case f" << std::endl;);
                    return std::make_pair(Nclose, false);
                }
                else
                {
                    VerboseJunctions(4, "JunctionNode case g (returning nullptr)" << std::endl;);
                    return std::make_pair(nullptr, false);
                }
            }
        }

        //New Implemnetation with update of the boundary nodes
        size_t contractJunctions(const std::deque<IntersectionTypeContainerType>& intersectionContainer)
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"        : "<<std::flush;

            DN.danglingBoundaryLoopNodes.clear();
            size_t nContracted=0;
            for (const auto& intersectionByThreadContainer : intersectionContainer)
            {
                for (const auto& intersection : intersectionByThreadContainer)
                {
                    const KeyType& key1(std::get<0>(intersection));
                    const KeyType& key2(std::get<1>(intersection));
                    const SegmentSegmentDistance<dim>& ssd(std::get<2>(intersection));
                    const double& t(ssd.t);
                    const double& u(ssd.u);

                    const auto L1(DN.networkLinks().get(key1));
                    const auto L2(DN.networkLinks().get(key2));

                    if(L1 && L2) // Links exist
                    {




                        VerboseJunctions(1,"forming junction "<<key1.first<<"->"<<key1.second<<", "
                                         /*                   */ <<key2.first<<"->"<<key2.second<<", "
                                         /*                   */ <<"dMin="<<ssd.dMin<<", "
                                         /*                   */ <<"@ ("<<t<<","<<u<<"), "<<std::flush);
                        VerboseJunctions(4, "JunctionNode for segment: " << key1.first << "->" << key1.second << " @ " << t << std::endl;);
                        std::pair<std::shared_ptr<NetworkNodeType>,bool> Ni(junctionNode(t,ssd.x0,L1));
                        VerboseJunctions(4, "JunctionNode for segment: " << key2.first << "->" << key2.second << " @ " << u << std::endl;);
                        std::pair<std::shared_ptr<NetworkNodeType>,bool> Nj(junctionNode(u,ssd.x1,L2));

                        if (Ni.first && Nj.first)
                        {

                            VerboseJunctions(1, "contracting " << Ni.first->sID << " " << Nj.first->sID << std::endl;);

                            if (Ni.first->sID != Nj.first->sID)
                            {
                                //Collect the boundary nodes for Nj for (const auto &loopN : loopNode.second.lock()->networkNode->loopNodes())
                                std::set<size_t> bndLoopNodes;

                                for (const auto &loopN : Ni.first->loopNodes())
                                {
                                    if (!loopN->periodicPlaneEdge.first) //if we are allowing for the boundary node contraction
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
                                }
                                for (const auto & loopN : Nj.first->loopNodes())
                                {
                                    if (!loopN->periodicPlaneEdge.first)
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
                                }

                                const bool success = DN.contract(Ni.first, Nj.first);
                                nContracted += success;
                                if (!success)
                                {
                                    if (Ni.second)
                                    {
                                        DN.removeNetworkNode(Ni.first->sID);
                                    }
                                    if (Nj.second)
                                    {
                                        DN.removeNetworkNode(Nj.first->sID);
                                    }
                                }
                                else
                                { //Update the boundary nodes here
                                            for (const size_t &nodeID : bndLoopNodes)
                                            {
                                                auto nodeIter(DN.loopNodes().find(nodeID));
                                                if (nodeIter != DN.loopNodes().end())
                                                {
                                                    auto bndNode(nodeIter->second.lock());

                                                    const auto pPrev(bndNode->periodicPrev());
                                                    const auto pNext(bndNode->periodicNext());
                                                    //                VerboseDislocationLoopNode(4,"pPrev= "<<pPrev->tag()<<" @ "<<pPrev->get_P().transpose()<<std::endl;);
                                                    //                VerboseDislocationLoopNode(4,"pNext= "<<pNext->tag()<<" @ "<<pNext->get_P().transpose()<<std::endl;);
                                                    if (pPrev && pNext)
                                                    {
                                                        const auto pPrevLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                                                        //                VerboseDislocationLoopNode(4,"pPrevLocal= "<<pPrevLocal.transpose()<<std::endl;);
                                                        const auto pNextLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                                                        //                VerboseDislocationLoopNode(4,"pNextLocal= "<<pNextLocal.transpose()<<std::endl;);
                                                        //                VerboseDislocationLoopNode(4,"periodicPlaneEdge->source= "<<bndNode->periodicPlaneEdge->source->transpose()<<std::endl;);
                                                        //                VerboseDislocationLoopNode(4,"bndNode->periodicPlaneEdge->sink= "<<bndNode->periodicPlaneEdge->sink->transpose()<<std::endl;);
                                                        //
                                                        // SegmentSegmentDistance<dim> ssd3(pPrev->get_P(), pNext->get_P(), bndNode->periodicPlaneEdge->meshIntersection->P0, bndNode->periodicPlaneEdge->meshIntersection->P1);
                                                        //                VerboseDislocationLoopNode(4,"ssd3.dMin= "<<ssd3.dMin<<std::endl;);

                                                        SegmentSegmentDistance<dim - 1> ssd1(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                                                        
                                                        if (ssd1.dMin < FLT_EPSILON)
                                                        {
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
                                                }
                                            }
                                }
                            }
                        }
                        else
                        {
                            //Either Ni or Nj might have been created...
                            //Get rid of that
                             if (Ni.first)
                             {
                                 if (Ni.second)
                                 {
                                     DN.removeNetworkNode(Ni.first->sID);
                                 }
                             }

                             if (Nj.first)
                             {
                                 if (Nj.second)
                                 {
                                     DN.removeNetworkNode(Nj.first->sID);
                                 }
                             }

                        }
                        
                    }
                }
            }

            std::cout<<" ("<<nContracted<<" contracted)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;


            std::cout << "Updating Boundary Nodes after junction contraction" << std::endl;
            if (DN.danglingBoundaryLoopNodes.size())
            {
                DN.updateBoundaryNodes();
            }

            return nContracted;

        }

 
        

        //This version also compares the network links burgers across the boundary
        void glissileJunctions(const double &dx)
        {
            const auto t0 = std::chrono::system_clock::now();
            std::cout << "        Forming Glissile Junctions: " << std::flush;

            std::deque<std::tuple<std::shared_ptr<NetworkNodeType>, std::shared_ptr<NetworkNodeType>, size_t, size_t>> glissDeq;

            //    std::deque<std::tuple<std::shared_ptr<NetworkNodeType>, std::shared_ptr<NetworkNodeType>, std::shared_ptr<NetworkNodeType>>> expDeq;

            for (const auto &linkIter : DN.networkLinks())
            {
                const auto link(linkIter.second.lock());

                //    if (link->isSessile() && link->loopLinks().size() > 1) // a junction
                if (link->loopLinks().size() > 1) // a junction
                {
                    //    const VectorDim chord(link->sink->get_P() - link->source->get_P());
                    //    const double chordNorm(chord.norm());
                    // const double dx_updated((DN.simulationParameters.isPeriodicSimulation() && link.second->isConnectedtoBoundaryNodes()) ? 10 * dx : dx);

                    if (fabs(link->burgers().norm() - 1.0) < FLT_EPSILON // a non-zero link with minimum Burgers
                        && link->chordLength() > dx)
                    {

                        //    const VectorDim unitChord(chord / chordNorm);
                        const VectorDim unitChord(link->chord() / link->chordLength());

                        if (!link->isGrainBoundarySegment() && !link->isBoundarySegment())
                        {

                            VerboseJunctions(2, "glissele junction, segment " << link->tag() << std::endl;);

                            for (const auto &gr : link->grains())
                            {
                                for (size_t s = 0; s < gr->singleCrystal->slipSystems().size(); ++s)
                                {
                                    const auto &slipSystem(gr->singleCrystal->slipSystems()[s]);
                                    if ((slipSystem->s.cartesian() - link->burgers()).norm() < FLT_EPSILON && fabs(slipSystem->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
                                    {
                                        VerboseJunctions(3, "glissDeq, emplacing" << std::endl;);

                                        glissDeq.emplace_back(link->source, link->sink, gr->grainID, s);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            size_t formedJunctions = 0;

            for (const auto &tup : glissDeq)
            {
                const std::shared_ptr<NetworkNodeType> &source(std::get<0>(tup));
                const std::shared_ptr<NetworkNodeType> &sink(std::get<1>(tup));
                const size_t &sourceID(source->sID);
                const size_t &sinkID(sink->sID);
                const size_t &grainID(std::get<2>(tup));
                const size_t &slipID(std::get<3>(tup));

                const auto isLink(DN.networkLinks().get(std::make_pair(sourceID, sinkID)));
                if (isLink)
                {
                    if (!isLink->hasZeroBurgers())
                    {
                        const VectorDim newNodeP(0.5 * (isLink->source->get_P() + isLink->sink->get_P())); //This new node is only on one side

                        const long int planeIndex(DN.poly.grain(grainID).singleCrystal->slipSystems()[slipID]->n.closestPlaneIndexOfPoint(newNodeP));
                        const GlidePlaneKey<dim> glissilePlaneKey(planeIndex, DN.poly.grain(grainID).singleCrystal->slipSystems()[slipID]->n);
                        const auto glidePlane(DN.glidePlaneFactory.getFromKey(glissilePlaneKey));
                        auto glissileLoop(DN.loops().create(DN.poly.grain(grainID).singleCrystal->slipSystems()[slipID]->s.cartesian(), glidePlane));

                        VerboseJunctions(3, "Glissile Junction from Link" << isLink->tag() << std::endl;);

                        if (!source->isBoundaryNode() && !sink->isBoundaryNode())
                        {
                            VerboseJunctions(3, "Case (a) internal node" << std::endl;);

                            std::shared_ptr<NetworkNodeType> newNode(DN.networkNodes().create(glidePlane->snapToPlane(newNodeP), VectorDim::Zero(), 1.0));

                            std::vector<std::shared_ptr<NetworkNodeType>> networkNodes;

                            networkNodes.push_back(sink);   // insert in reverse order, sink first, source second
                            networkNodes.push_back(source); // insert in reverse order, sink first, source second
                            networkNodes.push_back(newNode);

                            std::vector<std::shared_ptr<LoopNodeType>> loopNodes;

                            const auto periodicGlidePlane(DN.periodicGlidePlaneFactory->get(glidePlane->key));
                            const auto periodicPatch(periodicGlidePlane->getPatch(VectorDim::Zero()));

                            loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, sink, sink->get_P(), periodicPatch, std::make_pair (nullptr,nullptr)));
                            loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, source, source->get_P(), periodicPatch, std::make_pair (nullptr,nullptr)));

                            //New node cannot be a boundary node
                            loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, newNode, newNode->get_P(), periodicPatch, std::make_pair (nullptr,nullptr)));

                            DN.insertLoop(glissileLoop, loopNodes);

                            formedJunctions++;
                        }
                        else if (source->isBoundaryNode() && sink->isBoundaryNode())
                        {
                            //Skip this case
                            VerboseJunctions(3, "Case (b) both nodes at boundary" << std::endl;);
                        }
                        else
                        {
                            //Either of the node is boundary
                            VerboseJunctions(3, "Case (c) one of the nodes boundary " << std::endl;);

                            if (source->isBoundaryNode())
                            {
                                VerboseJunctions(3, "Case (c.A) Source is Boundary " << std::endl;);

                                //Get the loop Links of the networkLink
                                std::set<std::shared_ptr<NetworkNodeType>> networkNodeOtherEnd;
                                std::vector<std::pair<std::shared_ptr<NetworkNodeType>, VectorDim>> networkNodeOtherEndBnd;
                                std::set<std::shared_ptr<NetworkNodeType>> endNodesInserted;
                                std::set<std::shared_ptr<NetworkLinkType>> equivalentNetworkLinks;

                                for (const auto &sourceLN : source->loopNodes())
                                {
                                    const auto pPrev(sourceLN->periodicPrev());
                                    const auto pNext(sourceLN->periodicNext());

                                    if (pPrev && pNext)
                                    {
                                        if (pPrev->networkNode == sink)
                                        {
                                            networkNodeOtherEnd.insert(pNext->networkNode);
                                            const LoopNodeType *LNtemp(sourceLN);

                                            while ((LNtemp->periodicPlanePatch()->shift - pNext->periodicPlanePatch()->shift).squaredNorm() > FLT_EPSILON)
                                            {
                                                if (LNtemp->next.second->networkLink())
                                                {
                                                    equivalentNetworkLinks.insert(LNtemp->next.second->networkLink());
                                                }
                                                LNtemp = LNtemp->next.first;
                                                assert(LNtemp->periodicPlaneEdge.first != nullptr && "Next loop node must be on the edge");
                                                const VectorDim patchShift(LNtemp->periodicPlanePatch()->shift - sourceLN->periodicPlanePatch()->shift);
                                                if (endNodesInserted.find(LNtemp->networkNode) == endNodesInserted.end()) //A new networknode is being inserted
                                                {
                                                    networkNodeOtherEndBnd.emplace_back(LNtemp->networkNode, patchShift);
                                                    endNodesInserted.insert(LNtemp->networkNode);
                                                }
                                            }
                                        }
                                        else if (pNext->networkNode == sink)
                                        {
                                            networkNodeOtherEnd.insert(pPrev->networkNode);
                                            const LoopNodeType *LNtemp(sourceLN);

                                            while ((LNtemp->periodicPlanePatch()->shift - pPrev->periodicPlanePatch()->shift).squaredNorm() > FLT_EPSILON)
                                            {
                                                if (LNtemp->prev.second->networkLink())
                                                {
                                                    equivalentNetworkLinks.insert(LNtemp->prev.second->networkLink());
                                                }
                                                LNtemp = LNtemp->prev.first;
                                                assert(LNtemp->periodicPlaneEdge.first != nullptr && "Next loop node must be on the edge");
                                                const VectorDim patchShift(LNtemp->periodicPlanePatch()->shift - sourceLN->periodicPlanePatch()->shift);
                                                if (endNodesInserted.find(LNtemp->networkNode) == endNodesInserted.end()) //A new networknode is being inserted
                                                {
                                                    networkNodeOtherEndBnd.emplace_back(LNtemp->networkNode, patchShift);
                                                    endNodesInserted.insert(LNtemp->networkNode);
                                                }
                                            }
                                        }
                                        else
                                        {
                                            assert(false && "Network connectivity ill-defined for the boundary glissile junction");
                                        }
                                    }
                                    else
                                    {
                                        assert(false && "Periodic Previous and Periodic Next must exist");
                                    }
                                }

                                //Check all the network links have the same burgers vector
                                bool tempLink(true);
                                for (const auto &tnLink : equivalentNetworkLinks)
                                {
                                    tempLink = (tempLink && ((tnLink->burgers() - isLink->burgers()).squaredNorm() < FLT_EPSILON ||
                                                 (tnLink->burgers() + isLink->burgers()).squaredNorm() < FLT_EPSILON));
                                    if (!tempLink)
                                    {
                                        break;
                                    }
                                }

                                if (networkNodeOtherEnd.size() == 1 && tempLink)
                                {
                                    //If the boundary nodes are uniquely mapped, then form the glissile junctions
                                    //If the boundary node contraction is enabled this condition can be relaxed
                                    std::vector<std::shared_ptr<LoopNodeType>> loopNodes;

                                    const auto periodicGlidePlane(DN.periodicGlidePlaneFactory->get(glidePlane->key));
                                    const auto periodicPatch1(periodicGlidePlane->getPatch(VectorDim::Zero()));

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, sink, sink->get_P(), periodicPatch1, std::make_pair(nullptr,nullptr)));

                                    std::set<short int> edgeIDsFirstPatch;

                                    for (const auto &edge : periodicPatch1->edges())
                                    {
                                        if (edge->meshIntersection->contains(source->get_P()))
                                        {
                                            edgeIDsFirstPatch.insert(edge->edgeID);
                                        }
                                    }
                                    assert(edgeIDsFirstPatch.size() == 1 && "Glissile Junction Intersection at corner");

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, source, source->get_P(), periodicPatch1, std::make_pair(periodicPatch1->edges()[*edgeIDsFirstPatch.begin()],nullptr)));

                                    for (const auto &bndGNodes : networkNodeOtherEndBnd)
                                    {
                                        const auto periodicPatchI(periodicGlidePlane->getPatch(bndGNodes.second));
                                        std::set<short int> edgeIDsSecondPatch;

                                        for (const auto &edge : periodicPatchI->edges())
                                        {
                                            if (edge->meshIntersection->contains(bndGNodes.first->get_P()))
                                            {
                                                edgeIDsSecondPatch.insert(edge->edgeID);
                                            }
                                        }

                                        assert(edgeIDsSecondPatch.size() == 1 && "Glissile Junction Intersection at corner");
                                        loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, bndGNodes.first, bndGNodes.first->get_P() - periodicPatchI->shift, periodicPatchI, 
                                        std::make_pair(periodicPatchI->edges()[*edgeIDsSecondPatch.begin()],nullptr)));
                                    }

                                    const auto periodicPatch2(periodicGlidePlane->getPatch(networkNodeOtherEndBnd.rbegin()->second));

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, *networkNodeOtherEnd.begin(), (*networkNodeOtherEnd.begin())->get_P() - periodicPatch2->shift, periodicPatch2, std::make_pair(nullptr,nullptr)));

                                    for (auto iter = networkNodeOtherEndBnd.rbegin(); iter != networkNodeOtherEndBnd.rend(); ++iter)
                                    {
                                        const auto periodicPatchI(periodicGlidePlane->getPatch(iter->second));
                                        std::set<short int> edgeIDsSecondPatch;

                                        for (const auto &edge : periodicPatchI->edges())
                                        {
                                            if (edge->meshIntersection->contains(iter->first->get_P()))
                                            {
                                                edgeIDsSecondPatch.insert(edge->edgeID);
                                            }
                                        }

                                        assert(edgeIDsSecondPatch.size() == 1 && "Glissile Junction Intersection at corner");
                                        const auto newNode(DN.networkNodes().create(iter->first->get_P(), VectorDim::Zero(), 1.0));
                                        loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, newNode, newNode->get_P() - periodicPatchI->shift, periodicPatchI, std::make_pair(periodicPatchI->edges()[*edgeIDsSecondPatch.begin()],nullptr)));
                                    }
                                    //Add the node corresponding to the source loop node

                                    std::shared_ptr<NetworkNodeType> sourceEquivalentNode(DN.networkNodes().create(source->get_P(), VectorDim::Zero(), 1.0));
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, sourceEquivalentNode, sourceEquivalentNode->get_P(), periodicPatch1, std::make_pair(periodicPatch1->edges()[*edgeIDsFirstPatch.begin()],nullptr)));

                                    //Add a middle node
                                    const VectorDim newNodeP(0.5 * (source->get_P() + sink->get_P())); //This new node is only on one side
                                    std::shared_ptr<NetworkNodeType> middleNode(DN.networkNodes().create(glidePlane->snapToPlane(newNodeP), VectorDim::Zero(), 1.0));

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, middleNode, middleNode->get_P(), periodicPatch1, std::make_pair(nullptr,nullptr)));

                                    DN.insertLoop(glissileLoop, loopNodes);

                                    formedJunctions++;
                                }
                            }
                            else if (sink->isBoundaryNode())
                            {
                                VerboseJunctions(3, "Case (c.B) Sink is Boundary " << std::endl;);

                                //Get the loop Links of the networkLink
                                std::set<std::shared_ptr<NetworkNodeType>> networkNodeOtherEnd;
                                std::vector<std::pair<std::shared_ptr<NetworkNodeType>, VectorDim>> networkNodeOtherEndBnd;
                                std::set<std::shared_ptr<NetworkNodeType>> endNodesInserted;
                                std::set<std::shared_ptr<NetworkLinkType>> equivalentNetworkLinks;

                                for (const auto &sinkLN : sink->loopNodes())
                                {
                                    const auto pPrev(sinkLN->periodicPrev());
                                    const auto pNext(sinkLN->periodicNext());

                                    if (pPrev && pNext)
                                    {
                                        if (pPrev->networkNode == source)
                                        {
                                            networkNodeOtherEnd.insert(pNext->networkNode);
                                            const LoopNodeType *LNtemp(sinkLN);

                                            while ((LNtemp->periodicPlanePatch()->shift - pNext->periodicPlanePatch()->shift).squaredNorm() > FLT_EPSILON)
                                            {
                                                if (LNtemp->next.second->networkLink())
                                                {
                                                    equivalentNetworkLinks.insert(LNtemp->next.second->networkLink());
                                                }
                                                LNtemp = LNtemp->next.first;
                                                assert(LNtemp->periodicPlaneEdge.first != nullptr && "Next loop node must be on the edge");
                                                const VectorDim patchShift(LNtemp->periodicPlanePatch()->shift - sinkLN->periodicPlanePatch()->shift);
                                                if (endNodesInserted.find(LNtemp->networkNode) == endNodesInserted.end()) //A new networknode is being inserted
                                                {
                                                    networkNodeOtherEndBnd.emplace_back(LNtemp->networkNode, patchShift);
                                                    endNodesInserted.insert(LNtemp->networkNode);
                                                }
                                            }
                                        }
                                        else if (pNext->networkNode == source)
                                        {

                                            networkNodeOtherEnd.insert(pPrev->networkNode);
                                            const LoopNodeType *LNtemp(sinkLN);

                                            while ((LNtemp->periodicPlanePatch()->shift - pPrev->periodicPlanePatch()->shift).squaredNorm() > FLT_EPSILON)
                                            {
                                                if (LNtemp->prev.second->networkLink())
                                                {
                                                    equivalentNetworkLinks.insert(LNtemp->prev.second->networkLink());
                                                }
                                                LNtemp = LNtemp->prev.first;
                                                assert(LNtemp->periodicPlaneEdge.first != nullptr && "Next loop node must be on the edge");
                                                const VectorDim patchShift(LNtemp->periodicPlanePatch()->shift - sinkLN->periodicPlanePatch()->shift);
                                                if (endNodesInserted.find(LNtemp->networkNode) == endNodesInserted.end()) //A new networknode is being inserted
                                                {
                                                    networkNodeOtherEndBnd.emplace_back(LNtemp->networkNode, patchShift);
                                                    endNodesInserted.insert(LNtemp->networkNode);
                                                }
                                            }
                                        }
                                        else
                                        {
                                            std::cout << " Node under consideration " << sink->tag() << std::endl;
                                            assert(false && "Network connectivity ill-defined for the boundary glissile junction");
                                        }
                                    }
                                    else
                                    {
                                        assert(false && "Periodic Previous and Periodic Next must exist");
                                    }
                                }
                                //Check all the network links have the same burgers vector
                                bool tempLink(true);
                                for (const auto &tnLink : equivalentNetworkLinks)
                                {
                                    tempLink =(tempLink && ((tnLink->burgers() - isLink->burgers()).squaredNorm() < FLT_EPSILON ||
                                                 (tnLink->burgers() + isLink->burgers()).squaredNorm() < FLT_EPSILON));
                                    if (!tempLink)
                                    {
                                        break;
                                    }
                                }

                                if (networkNodeOtherEnd.size() == 1 && tempLink)
                                {
                                    //If the boundary nodes are uniquely mapped, then form the glissile junctions
                                    //If the boundary node contraction is enabled this condition can be relaxed
                                    std::vector<std::shared_ptr<LoopNodeType>> loopNodes;

                                    const auto periodicGlidePlane(DN.periodicGlidePlaneFactory->get(glidePlane->key));
                                    const auto periodicPatch1(periodicGlidePlane->getPatch(VectorDim::Zero()));

                                    std::set<short int> edgeIDsFirstPatch;

                                    for (const auto &edge : periodicPatch1->edges())
                                    {
                                        if (edge->meshIntersection->contains(sink->get_P()))
                                        {
                                            edgeIDsFirstPatch.insert(edge->edgeID);
                                        }
                                    }
                                    assert(edgeIDsFirstPatch.size() == 1 && "Glissile Junction Intersection at corner");

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, sink, sink->get_P(), periodicPatch1, std::make_pair(periodicPatch1->edges()[*edgeIDsFirstPatch.begin()],nullptr)));

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, source, source->get_P(), periodicPatch1, std::make_pair(nullptr,nullptr)));

                                    //Add a middle node
                                    const VectorDim newNodeP(0.5 * (source->get_P() + sink->get_P())); //This new node is only on one side
                                    std::shared_ptr<NetworkNodeType> middleNode(DN.networkNodes().create(glidePlane->snapToPlane(newNodeP), VectorDim::Zero(), 1.0));
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, middleNode, middleNode->get_P(), periodicPatch1, std::make_pair(nullptr,nullptr)));

                                    //Add sink equivalent node
                                    std::shared_ptr<NetworkNodeType> sinkEquivalentNode(DN.networkNodes().create(sink->get_P(), VectorDim::Zero(), 1.0));
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, sinkEquivalentNode, sinkEquivalentNode->get_P(), periodicPatch1, std::make_pair(periodicPatch1->edges()[*edgeIDsFirstPatch.begin()],nullptr)));

                                    for (const auto &bndGNodes : networkNodeOtherEndBnd)
                                    {
                                        const auto periodicPatchI(periodicGlidePlane->getPatch(bndGNodes.second));
                                        std::set<short int> edgeIDsSecondPatch;

                                        for (const auto &edge : periodicPatchI->edges())
                                        {
                                            if (edge->meshIntersection->contains(bndGNodes.first->get_P()))
                                            {
                                                edgeIDsSecondPatch.insert(edge->edgeID);
                                            }
                                        }

                                        assert(edgeIDsSecondPatch.size() == 1 && "Glissile Junction Intersection at corner");
                                        const auto newNode(DN.networkNodes().create(bndGNodes.first->get_P(), VectorDim::Zero(), 1.0));
                                        loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, newNode, newNode->get_P() - periodicPatchI->shift, periodicPatchI, std::make_pair(periodicPatchI->edges()[*edgeIDsSecondPatch.begin()],nullptr)));
                                    }

                                    const auto periodicPatch2(periodicGlidePlane->getPatch(networkNodeOtherEndBnd.rbegin()->second));

                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, *networkNodeOtherEnd.begin(), (*networkNodeOtherEnd.begin())->get_P() - periodicPatch2->shift, periodicPatch2, std::make_pair(nullptr,nullptr)));

                                    for (auto iter = networkNodeOtherEndBnd.rbegin(); iter != networkNodeOtherEndBnd.rend(); ++iter)
                                    {
                                        const auto periodicPatchI(periodicGlidePlane->getPatch(iter->second));
                                        std::set<short int> edgeIDsSecondPatch;

                                        for (const auto &edge : periodicPatchI->edges())
                                        {
                                            if (edge->meshIntersection->contains(iter->first->get_P()))
                                            {
                                                edgeIDsSecondPatch.insert(edge->edgeID);
                                            }
                                        }

                                        assert(edgeIDsSecondPatch.size() == 1 && "Glissile Junction Intersection at corner");
                                        loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, iter->first, iter->first->get_P() - periodicPatchI->shift, periodicPatchI, std::make_pair(periodicPatchI->edges()[*edgeIDsSecondPatch.begin()],nullptr)));
                                    }

                                    DN.insertLoop(glissileLoop, loopNodes);

                                    formedJunctions++;
                                }
                            }
                            else
                            {
                                assert(false && "Error (Source and Sink must not be boundary for the glissile junctions)");
                            }
                        }
                    }
                }
            }
            std::cout << "(" << formedJunctions << " junctions)" << magentaColor << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]" << defaultColor << std::endl;
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
        /* init */,maxJunctionIterations(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("maxJunctionIterations",true))
        /* init */,verboseJunctions(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseJunctions",true))
        /* init */,infiniteLineLength(10000.0)
        {
            
        }

        static void initFromFile(const std::string& fileName)
        {
            collisionTol=TextFileParser(fileName).readScalar<double>("collisionTol",true);
        }
        
        /**********************************************************************/
        void formJunctions(const double& dx)
        {
#ifdef _OPENMP
            // const size_t nThreads = omp_get_max_threads();
            const size_t nThreads = 1;
#else
            const size_t nThreads = 1;
#endif
            

            
            size_t nContracted=1;
            size_t iterations=0;
            while(nContracted && iterations<maxJunctionIterations)
            {
                std::deque<IntersectionTypeContainerType> intersectionContainer;
                intersectionContainer.resize(nThreads);
                findIntersections(intersectionContainer,nThreads);

                nContracted=contractJunctions(intersectionContainer);
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
