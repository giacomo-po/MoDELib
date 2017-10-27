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

namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationJunctionFormation
    {
        
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        
        typedef typename DislocationNetworkType::IsNetworkEdgeType IsNetworkLinkType;
        typedef typename DislocationNetworkType::IsNodeType IsNodeType;
        
        enum {dim=3};
        enum {pOrder=3};
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        //        typedef typename DislocationNetworkType::NetworkNodeContainerType NetworkNodeContainerType;
        
        //		typedef std::pair<const LinkType*, double> EdgeIntersectionType;
        
        typedef std::pair<size_t,size_t> EdgeIDType;
        
        typedef std::tuple<EdgeIDType,EdgeIDType,SegmentSegmentDistance<dim>> IntersectionType;
        
        
        //        typedef std::pair<std::pair<size_t,size_t>, double> EdgeIntersectionType;
        
        //        typedef std::pair<EdgeIntersectionType,EdgeIntersectionType> EdgeIntersectionPairType;
        //		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
        typedef std::deque<IntersectionType,Eigen::aligned_allocator<IntersectionType>> IntersectionTypeContainerType;
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        //! The tolerance (in units of distance) used for collision detection
        static double collisionTol;
        //        static bool useVertexEdgeJunctions;
        
        
        /* Constructor ********************************************************/
        DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {
            
        }
        
        /* findIntersections **************************************************/
        //		EdgeIntersectionPairContainerType findIntersections(const double& avoidNodeIntersection) const
        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                               //                               std::deque<std::deque<int>>& dirVector,
                               const size_t& nThreads) const __attribute__ ((deprecated)) // SESSILE PLANE NORMAL DOES NOT EXIST ANYMORE
        
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
            
            //                                        const double avoidNodeIntersection=0.2;
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Finding Junctions ("<<nThreads<<" threads)... "<<std::flush;
            
            // Create an EqualConstIteratorRange over links
            //            EqualConstIteratorRange<NetworkLinkContainerType> eir(DN.linkBegin(),DN.linkEnd(),nThreads);
            N2IteratorRange<typename NetworkLinkContainerType::const_iterator> eir(DN.links().begin(),DN.links().end(),nThreads);
            
            
            //            std::cout<<"#links="<<DN.linkOrder()<<std::endl;
            //            for(auto pair : eir)
            //            {
            //                std::cout<<std::distance(pair.first,pair.second)<<std::endl;
            //            }
            //
            //            std::vector<int> threadVector(nThreads,0);
            
            
            const double avoidNodeIntersection=0.1;
            
            
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
                                SegmentSegmentDistance<dim> ssd(linkIterA->second->source->get_P(),
                                                                linkIterA->second->sink->get_P(),
                                                                linkIterB->second->source->get_P(),
                                                                linkIterB->second->sink->get_P());
                                
                                if(ssd.dMin<collisionTol)
                                {
                                    const bool intersectionIsSourceSource(   ssd.t  <avoidNodeIntersection
                                                                          && ssd.u <avoidNodeIntersection
                                                                          && linkIterA->second->source->sID==linkIterB->second->source->sID);
                                    
                                    const bool   intersectionIsSourceSink(   ssd.t  <avoidNodeIntersection
                                                                          && ssd.u >1.0-avoidNodeIntersection
                                                                          && linkIterA->second->source->sID==linkIterB->second->sink->sID);
                                    
                                    const bool   intersectionIsSinkSource(   ssd.t  > 1.0-avoidNodeIntersection
                                                                          && ssd.u <avoidNodeIntersection
                                                                          && linkIterA->second->sink->sID==linkIterB->second->source->sID);
                                    
                                    const bool     intersectionIsSinkSink(   ssd.t  > 1.0-avoidNodeIntersection
                                                                          && ssd.u > 1.0-avoidNodeIntersection
                                                                          && linkIterA->second->sink->sID==linkIterB->second->sink->sID);
                                    
                                    const bool frankRule(linkIterA->second->burgers().dot(linkIterB->second->burgers())*linkIterA->second->chord().dot(linkIterB->second->chord())<=0.0);
                                    const bool isValidJunction(frankRule ||
                                                               linkIterA->second->is_boundarySegment() || linkIterB->second->is_boundarySegment() ||
                                                               linkIterA->second->grainBoundaries().size() || linkIterB->second->grainBoundaries().size());
                                    
                                    
                                    if(   isValidJunction
                                       && !intersectionIsSourceSource
                                       && !intersectionIsSourceSink
                                       && !intersectionIsSinkSource
                                       && !intersectionIsSinkSink)
                                    {
#ifdef _OPENMP
                                        intersectionContainer[omp_get_thread_num()].emplace_back(linkIterA->second->nodeIDPair,
                                                                                                 linkIterB->second->nodeIDPair,
                                                                                                 ssd);
#else
                                        intersectionContainer[0].emplace_back(linkIterA->second->nodeIDPair,
                                                                              linkIterB->second->nodeIDPair,
                                                                              ssd);
#endif
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
                //                assert(intersectionByThreadContainer.size()==dirVector[tt].size());
                nIntersections+=intersectionByThreadContainer.size();
            }
            model::cout<<nIntersections<<" physical intersections. ";
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            }
            
            /**********************************************************************/
            void formJunctions(const bool& use_junctions, const double& dx)
            {
                
                if (use_junctions)
                {
                    //! 1- Initialize intersectionContainer calling findIntersections deque<Pair<Pair<Link*double>,Pair<Link*,double>>>
                    std::deque<IntersectionTypeContainerType> intersectionContainer;
                    //                std::deque<std::deque<int>> dirVector;
                    
#ifdef _OPENMP
                    const size_t nThreads = omp_get_max_threads();
#else
                    const size_t nThreads = 1;
#endif
                    
                    intersectionContainer.resize(nThreads);
                    //                dirVector.resize(nThreads);
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
                                
                                DN.contract(Ni,Nj);
                                
                                
                                
                                
                            }
                            //std::cout<<"done forming Junction "<<std::endl;
                            
                        }
                    } // loop over threads
                    model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                }
                
                
                
                
                
            }
            
            /**********************************************************************/
            void breakZeroLengthJunctions()
            {
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"		Breaking zero-length Junctions... "<<std::flush;
                
                std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecompFirst;
                std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecompSecond;
                
                
                // TO DO: PARALLELIZE THIS LOOP
                for (auto& nodePair : DN.nodes())
                {
                    const auto temp=nodePair.second.edgeDecomposition();
                    //std::deque<std::pair<size_t,size_t> > temp=nodePair.second.edgeDecomposition();
                    
                    if(temp.first.size() && temp.second.size())
                    {
                        nodeDecompFirst.emplace_back(nodePair.second.sID,temp.first);
                        nodeDecompSecond.emplace_back(nodePair.second.sID,temp.second);
                    }
                }
                
                int broken=0;
                
                for(size_t n=0;n<nodeDecompFirst.size();++n)
                {
                    const size_t& i=nodeDecompFirst[n].first;
                    auto Ni=DN.node(i);
                    assert(Ni.first);
                    //                size_t m=DN.insertVertex(Ni.second->get_P());
                    
                    // Check that Links still exist
                    //std::deque<VectorDimD> linkDirs;
                    
                    bool linksFirstexist=true;
                    VectorDimD avrFirst=VectorDimD::Zero();
                    for(size_t d=0;d<nodeDecompFirst[n].second.size();++d)
                    {
                        const size_t& j=nodeDecompFirst[n].second[d].first;
                        const size_t& k=nodeDecompFirst[n].second[d].second;
                        if(i==j)
                        {
                            auto Lik=DN.link(i,k);
                            if(Lik.first)
                            {
                                avrFirst+=Lik.second->chord().normalized();
                            }
                            linksFirstexist*=Lik.first;
                        }
                        else if(i==k)
                        {
                            auto Lji=DN.link(j,i);
                            if(Lji.first)
                            {
                                avrFirst-=Lji.second->chord().normalized();
                            }
                            linksFirstexist*=Lji.first;
                        }
                        else
                        {
                            assert(0 && "i must be equal to either j or k.");
                        }
                    }
                    const double avrFirstNorm=avrFirst.norm();
                    if(avrFirstNorm>FLT_EPSILON)
                    {
                        avrFirst/=avrFirstNorm;
                    }
                    
                    bool linksSecondexist=true;
                    VectorDimD avrSecond=VectorDimD::Zero();
                    for(size_t d=0;d<nodeDecompSecond[n].second.size();++d)
                    {
                        const size_t& j=nodeDecompSecond[n].second[d].first;
                        const size_t& k=nodeDecompSecond[n].second[d].second;
                        if(i==j)
                        {
                            auto Lik=DN.link(i,k);
                            if(Lik.first)
                            {
                                avrSecond+=Lik.second->chord().normalized();
                            }
                            linksSecondexist*=Lik.first;
                        }
                        else if(i==k)
                        {
                            auto Lji=DN.link(j,i);
                            if(Lji.first)
                            {
                                avrSecond-=Lji.second->chord().normalized();
                            }
                            linksSecondexist*=Lji.first;
                        }
                        else
                        {
                            assert(0 && "i must be equal to either j or k.");
                        }
                    }
                    const double avrSecondNorm=avrSecond.norm();
                    if(avrSecondNorm>FLT_EPSILON)
                    {
                        avrSecond/=avrSecondNorm;
                    }
                    
                    
                    if(linksFirstexist && linksSecondexist && avrSecond.dot(avrFirst)<-0.0)
                    {
                        std::cout<<"NodeBreaking "<<Ni.second->sID<<" "<<avrSecond.dot(avrFirst)<<std::endl;
                        
                        size_t m=DN.insertVertex(Ni.second->get_P(),Ni.second->grain.grainID).first->first;
                        
                        for(size_t d=0;d<nodeDecompFirst[n].second.size();++d)
                        {
                            const size_t& j=nodeDecompFirst[n].second[d].first;
                            const size_t& k=nodeDecompFirst[n].second[d].second;
                            if(i==j)
                            {
                                auto Lik=DN.link(i,k);
                                assert(Lik.first);
                                DN.connect(m,k,Lik.second->flow);
                                DN.template disconnect<0>(i,k);
                            }
                            else if(i==k)
                            {
                                auto Lji=DN.link(j,i);
                                assert(Lji.first);
                                DN.connect(j,m,Lji.second->flow);
                                DN.template disconnect<0>(j,i);
                            }
                            else
                            {
                                assert(0 && "i must be equal to either j or k.");
                            }
                        }
                        
                        
                        broken++;
                    }
                }
                model::cout<<broken<<" broken."<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            
            
            };
            
            // Declare Static Data
            template <typename DislocationNetworkType>
            double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
            
            } // namespace model
#endif
            
