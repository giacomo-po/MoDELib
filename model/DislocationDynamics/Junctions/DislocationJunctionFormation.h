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
//#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>
#include <model/DislocationDynamics/DislocationNetworkRemesh.h>
#include <model/MPI/MPIcout.h>
#include <model/Threads/EqualIteratorRange.h>
#include <model/Threads/N2IteratorRange.h>
#define VerboseJunctions(N,x) if(verboseJunctions>=N){model::cout<<x;}

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
//        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef std::pair<size_t,size_t> EdgeIDType;
        
        typedef std::tuple<EdgeIDType,EdgeIDType,SegmentSegmentDistance<dim>> IntersectionType;
        
        
        //        typedef std::pair<std::pair<size_t,size_t>, double> EdgeIntersectionType;
        
        //        typedef std::pair<EdgeIntersectionType,EdgeIntersectionType> EdgeIntersectionPairType;
        //		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
        typedef std::deque<IntersectionType,Eigen::aligned_allocator<IntersectionType>> IntersectionTypeContainerType;
        
        /**********************************************************************/
        void insertIntersection(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                                const std::pair<size_t,size_t>& linkApair,
                                const std::pair<size_t,size_t>& linkBpair,
                                const SegmentSegmentDistance<dim>& ssd,
                                const double& currentcCollisionTOL)
        {
            if(ssd.dMin<currentcCollisionTOL)
            {
#ifdef _OPENMP
                intersectionContainer[omp_get_thread_num()].emplace_back(linkApair,
                                                                         linkBpair,
                                                                         ssd);
#else
                intersectionContainer[0].emplace_back(linkApair,
                                                      linkBpair,
                                                      ssd);
#endif
            }
        }
        
        /**********************************************************************/
        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                               const size_t& nThreads)
        
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Finding Junctions ("<<nThreads<<" threads)... "<<std::flush;
            
            N2IteratorRange<typename NetworkLinkContainerType::const_iterator> eir(DN.links().begin(),DN.links().end(),nThreads);
            
            //! 2- loop over all links and determine their intersections
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t thread=0;thread<eir.size();thread++)
            {
                for (typename NetworkLinkContainerType::const_iterator linkIterA=eir[thread].first;linkIterA!=eir[thread].second;linkIterA++)
                {
                    if(!linkIterA->second->hasZeroBurgers()) // don't loop over zero-Burgers segments
                    {
                        for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.links().end();linkIterB++)
                        {
                            if (   linkIterA->second->sID!=linkIterB->second->sID           // don't intersect with itself
                                && !linkIterB->second->hasZeroBurgers()   // don't intersect with zero-Burgers segments
                                )
                            {
                                
                                const bool frankRule(linkIterA->second->burgers().dot(linkIterB->second->burgers())*linkIterA->second->chord().dot(linkIterB->second->chord())<=0.0);
                                const bool isValidJunction(frankRule ||
                                                           linkIterA->second->isBoundarySegment() || linkIterB->second->isBoundarySegment() ||
                                                           linkIterA->second->grainBoundaries().size() || linkIterB->second->grainBoundaries().size());
                                
                                if(isValidJunction)
                                {
                                    
                                    
                                    double currentcCollisionTOL=collisionTol;
                                    if(   linkIterA->second->glidePlaneNormal().squaredNorm()>FLT_EPSILON
                                       && linkIterB->second->glidePlaneNormal().squaredNorm()>FLT_EPSILON
                                       && linkIterA->second->glidePlaneNormal().cross(linkIterB->second->glidePlaneNormal()).squaredNorm()<FLT_EPSILON)
                                    {// segments on parallel or coincident planes
                                        currentcCollisionTOL=FLT_EPSILON;
                                    }
                                    
                                    const bool intersectionIsSourceSource(linkIterA->second->source->sID==linkIterB->second->source->sID);
                                    
                                    const bool intersectionIsSourceSink(linkIterA->second->source->sID==linkIterB->second->sink->sID);
                                    
                                    const bool intersectionIsSinkSource(linkIterA->second->sink->sID==linkIterB->second->source->sID);
                                    
                                    const bool intersectionIsSinkSink(linkIterA->second->sink->sID==linkIterB->second->sink->sID);
                                    
                                    if(intersectionIsSourceSource)
                                    {
                                        // intersect sink of A with link B
                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->sink->get_P(), // fake degenerete segment at sink of A
                                                                         linkIterA->second->sink->get_P(),  // fake degenerete segment at sink of A
                                                                         linkIterB->second->source->get_P(),
                                                                         linkIterB->second->sink->get_P(),
                                                                         1.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
                                        
                                        // intersect sink of B with link A
                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
                                                                         linkIterA->second->sink->get_P(),
                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
                                                                         1.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
                                        
                                    }
                                    else if(intersectionIsSourceSink)
                                    {
                                        // intersect sink of A with link B
                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->sink->get_P(), // fake degenerete segment at sink of A
                                                                         linkIterA->second->sink->get_P(),  // fake degenerete segment at sink of A
                                                                         linkIterB->second->source->get_P(),
                                                                         linkIterB->second->sink->get_P(),
                                                                         1.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
                                        
                                        // intersect source of B with link A
                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
                                                                         linkIterA->second->sink->get_P(),
                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
                                                                         0.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
                                        
                                        
                                    }
                                    else if(intersectionIsSinkSource)
                                    {
                                        // intersect source of A with link B
                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->source->get_P(), // fake degenerete segment at source of A
                                                                         linkIterA->second->source->get_P(),  // fake degenerete segment at source of A
                                                                         linkIterB->second->source->get_P(),
                                                                         linkIterB->second->sink->get_P(),
                                                                         0.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
                                        
                                        
                                        // intersect sink of B with link A
                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
                                                                         linkIterA->second->sink->get_P(),
                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
                                                                         linkIterB->second->sink->get_P(),    // fake degenerete segment at sink of B
                                                                         1.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
                                        
                                        
                                    }
                                    else if(intersectionIsSinkSink)
                                    {
                                        // intersect source of A with link B
                                        SegmentSegmentDistance<dim> ssdA(linkIterA->second->source->get_P(), // fake degenerete segment at source of A
                                                                         linkIterA->second->source->get_P(),  // fake degenerete segment at source of A
                                                                         linkIterB->second->source->get_P(),
                                                                         linkIterB->second->sink->get_P(),
                                                                         0.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdA,currentcCollisionTOL);
                                        
                                        
                                        // intersect source of B with link A
                                        SegmentSegmentDistance<dim> ssdB(linkIterA->second->source->get_P(),
                                                                         linkIterA->second->sink->get_P(),
                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
                                                                         linkIterB->second->source->get_P(),    // fake degenerete segment at source of B
                                                                         0.0);
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssdB,currentcCollisionTOL);
                                        
                                    }
                                    else
                                    {// segments don't share vertices
                                        
                                        SegmentSegmentDistance<dim> ssd(linkIterA->second->source->get_P(),
                                                                        linkIterA->second->sink->get_P(),
                                                                        linkIterB->second->source->get_P(),
                                                                        linkIterB->second->sink->get_P()
                                                                        );
                                        insertIntersection(intersectionContainer,linkIterA->second->nodeIDPair,linkIterB->second->nodeIDPair,ssd,currentcCollisionTOL);
                                    }
                                }
                            } // end don't self-intersect
                        }
                    }
                } // end loop over first segment
            }// end loop ever threads
            
            int nIntersections=0;
            for (const auto& intersectionByThreadContainer : intersectionContainer)
            {
                nIntersections+=intersectionByThreadContainer.size();
            }
            model::cout<<nIntersections<<" physical intersections. "<<std::flush;
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        void junctionStep(const double& dx)
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
                        
                        auto  Ni=L1.second->source;
                        if((ssd.x0-L1.second->source->get_P()).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                        {
                            if((ssd.x0-L1.second->sink->get_P()).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                            {
                                Ni=DN.expand(key1.first,key1.second,ssd.x0);
                            }
                            else
                            {
                                Ni=L1.second->sink;
                            }
                        }
                        
                        auto Nj=L2.second->source;
                        if((ssd.x1-L2.second->source->get_P()).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                        {
                            if((ssd.x1-L2.second->sink->get_P()).norm()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                            {
                                Nj=DN.expand(key2.first,key2.second,ssd.x1);
                            }
                            else
                            {
                                Nj=L2.second->sink;
                            }
                        }
                        
                        VerboseJunctions(1,"forming junction "<<key1.first<<"->"<<key1.second<<", "
                                         /*                   */ <<key2.first<<"->"<<key2.second<<", "
                                         /*                   */ <<"@ ("<<t<<","<<u<<"), "
                                         /*                   */ <<"contracting "<<Ni->sID<<" "<<Nj->sID<<std::endl;);
                        
                        DN.contract(Ni,Nj);
                        
                    }
                }
            }
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        //! The tolerance (in units of distance) used for collision detection
        static double collisionTol;
        //        static bool useVertexEdgeJunctions;
        static int verboseJunctions;
        
        
        /**********************************************************************/
        DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {
            
        }
        
        /**********************************************************************/
        void formJunctions(const bool& use_junctions, const double& dx)
        {
            if (use_junctions)
            {
                junctionStep(dx);
            }
        }
        
    };
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
    
    template <typename DislocationNetworkType>
    int DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions=0;
    
    
} // namespace model
#endif

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

