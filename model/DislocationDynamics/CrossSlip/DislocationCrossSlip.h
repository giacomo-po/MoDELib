/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_H_
#define model_DislocationCrossSlip_H_

#include <utility> // for std::pair
#include <vector>
#include <Eigen/Dense>
#include <model/Geometry/SegmentSegmentDistance.h>
//#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>
#include <model/DislocationDynamics/DislocationNetworkRemesh.h>
#include <model/Geometry/PlanePlaneIntersection.h>

#include <model/MPI/MPIcout.h>
#include <model/Threads/EqualIteratorRange.h>
#include <model/Threads/N2IteratorRange.h>

#ifndef NDEBUG
#define VerboseCrossSlip(N,x) if(verboseCrossSlip>=N){model::cout<<x;}
#else
#define VerboseCrossSlip(N,x)
#endif


namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationCrossSlip
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename DislocationNetworkType::IsNetworkEdgeType IsNetworkLinkType;
        typedef typename DislocationNetworkType::IsNodeType IsNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        
        typedef std::tuple<size_t,size_t,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
        
        /**********************************************************************/
        CrossSlipContainerType findCrossSlipSegments() const
        {
            
            const double sinCrossSlipRad(std::sin(crossSlipDeg*M_PI/180.0));
            CrossSlipContainerType crossSlipDeq;
            
            for(const auto& link : DN.links())
            {
                if(   !link.second->isBoundarySegment()
                   && !link.second->isGrainBoundarySegment()
                   && !link.second->isGrainBoundarySegment()
                   && !link.second->hasZeroBurgers()
                   && link.second->isGlissile()
                   && link.second->chord().normalized().cross(link.second->burgers().normalized()).norm()<=sinCrossSlipRad
                   )
                {
                    
                    assert(link.second->grains().size()==1 && "THERE SHOULD BE ONLY ONE GRAIN FOR AN INTERNAL SEGMENT");
                    
                    // Find slip system on which screw segment has highest glide force
                    const Grain<dim>& grain(**link.second->grains().begin());
                    const VectorDim pki=link.second->pkIntegral();
                    std::map<double,int> ssMap;
                    for(size_t s=0;s<grain.slipSystems().size();++s)
                    {
                        const auto& slipSystem(grain.slipSystems()[s]);
                        //const VectorDim ncs=slipSystem.unitNormal.normalized();
                        if((slipSystem.s.cartesian()-link.second->burgers()).norm()<FLT_EPSILON)
                        {
                            const double pkGlide=(pki-pki.dot(slipSystem.unitNormal)*slipSystem.unitNormal).norm();
                            ssMap.emplace(pkGlide,s);
                        }
                    }
                    
                    if(ssMap.size())
                    {
                        const auto& crosSlipSystem(grain.slipSystems()[ssMap.rbegin()->second]); // last element in map has highest pkGlide
                        if(crosSlipSystem.unitNormal.cross(link.second->glidePlaneNormal()).squaredNorm()>FLT_EPSILON)
                        {// cross slip system is different than original system
                            crossSlipDeq.emplace_back(link.second->source->sID,link.second->sink->sID,grain.grainID,ssMap.rbegin()->second);
                        }
                    }
                }
            }
            
            return crossSlipDeq;
        }
        
    public:
        
        
        static int verboseCrossSlip;
        static double crossSlipDeg;
        
        /**********************************************************************/
        DislocationCrossSlip(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {
            if(DN.use_crossSlip)
            {
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"		CrossSlip "<<std::flush;
                
                
                const CrossSlipContainerType crossSlipDeq=findCrossSlipSegments();
                VerboseCrossSlip(1,"crossSlipDeq.size()="<<crossSlipDeq.size()<<std::endl;);

                
                
                for(const auto& tup : crossSlipDeq)
                {
                    const size_t& sourceID(std::get<0>(tup));
                    const size_t& sinkID(std::get<1>(tup));
                    const size_t& grainID(std::get<2>(tup));
                    const size_t& slipID(std::get<3>(tup));
                    
                    const IsNodeType isSource(DN.node(sourceID));
                    const IsNodeType isSink(DN.node(sinkID));
                    const auto isLink(DN.link(sourceID,sinkID));
                    
                    const auto& crosSlipSystem(DN.poly.grain(grainID).slipSystems()[slipID]); // last element in map has highest pkGlide
                    
                    if(isSource.first && isSink.first && isLink.first)
                    {
                        
                        // Align source and sink to perfect screw orientation
                        const VectorDim midPoint(0.5*(isSource.second->get_P()+isSink.second->get_P()));
                        const int height=LatticePlane::computeHeight(crosSlipSystem.n,midPoint).second;
                        const VectorDim planePoint(height*crosSlipSystem.n.planeSpacing()*crosSlipSystem.unitNormal);
                        
                        PlanePlaneIntersection<dim> ppi(midPoint,isLink.second->glidePlaneNormal(),
                                                        planePoint,crosSlipSystem.unitNormal);
                        
                        const VectorDim newSourceP(ppi.P+(isSource.second->get_P()-ppi.P).dot(ppi.d)*ppi.d);
                        const VectorDim newSinkP(ppi.P+(isSink.second->get_P()-ppi.P).dot(ppi.d)*ppi.d);
                        
                        if(   isSource.second->isMovableTo(newSourceP)
                           &&   isSink.second->isMovableTo(newSinkP))
                        {
                            
                            VerboseCrossSlip(1,"cross-slip "<<sourceID<<"->"<<sinkID<<std::endl;);
                            
                            // Re-align source and sink
                            isSource.second->set_P(newSourceP);
                            isSink.second->set_P(newSinkP);
                            
                            // Construct and insert new loop in conjugate plane
                            const VectorDim newNodeP(0.5*(isSource.second->get_P()+isSink.second->get_P()));
                            const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
                            
                            std::vector<size_t> nodeIDs;
                            
                            nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
                            nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
                            nodeIDs.push_back(newNodeID);
                            
                            DN.insertLoop(nodeIDs,
                                          DN.poly.grain(grainID).slipSystems()[slipID].s.cartesian(),
                                          DN.poly.grain(grainID).slipSystems()[slipID].unitNormal,
                                          newNodeP,
                                          grainID);
                        }
                    }
                }
                
                DN.clearDanglingNodes();
                
                std::cout<<"crossSlipDeq.size="<<crossSlipDeq.size()<<std::endl;
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                
            }
            
        }
        
    };
    
    template <typename DislocationNetworkType>
    int DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip=0;
    
    template <typename DislocationNetworkType>
    double DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg=2.0;
    
}
#endif


