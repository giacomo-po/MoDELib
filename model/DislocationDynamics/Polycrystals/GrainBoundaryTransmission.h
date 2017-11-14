/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2016 by Nathaniel Burbery <nbur049@aucklanduni.ac.nz>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundaryTransmission_H_
#define model_GrainBoundaryTransmission_H_

#include <deque>
#include <utility>
#include <math.h>
#include <model/Utilities/TerminalColors.h>
//#include <model/DislocationDynamics/Polycrystals/TransmitSegment.h>
#include <model/MPI/MPIcout.h> // defines mode::cout


namespace model
{
    
    
    
    template <int modelID>
    struct GrainBoundaryTransmissionEnergyModel
    {
        
    };

    template <>
    struct GrainBoundaryTransmissionEnergyModel<1>
    {
        
        static constexpr double alpha=1.0;
        static constexpr double lambda=5.0;
        
        template<typename SegmentType>
        static double dG(const SegmentType& seg,
                         const Eigen::Matrix<double,SegmentType::dim,1>& grain2OutNormal,
                         const Eigen::Matrix<double,SegmentType::dim,1>& b2,
                         const Eigen::Matrix<double,SegmentType::dim,1>& n2)
        {
        
            int q=seg.qOrder/2;
            
            double temp=0.0;
            std::cout<<"here 0"<<std::endl;
            std::cout<<"seg.quadratureParticleContainer.size()="<<seg.quadratureParticleContainer.size()<<std::endl;
//            if((seg.stressAtQuadrature(q)*seg.burgers()).cross(seg.chord().normalized()).dot(grain2OutNormal)<0.0)
            const Eigen::Matrix<double,SegmentType::dim,1> pk((seg.stressAtQuadrature(q)*b2).cross(seg.chord().normalized()));
            
                        std::cout<<"here 1"<<std::endl;
            
            const Eigen::Matrix<double,SegmentType::dim,1> pkg(pk-pk.dot(n2)*n2);
            
                if(pkg.dot(grain2OutNormal)<0.0)
            {// pk force on segment must be pushing into new grain
                const Eigen::Matrix<double,SegmentType::dim,1>& b1(seg.burgers());
                const double tau=fabs((seg.stressAtQuadrature(q)*b2.normalized()).dot(n2));

                temp=0.5*alpha*((b1-b2).squaredNorm()+b2.squaredNorm()-b1.squaredNorm())-tau*b2.norm()*lambda;
            }
            return temp;
        }
        
    };
    
    template <typename DislocationNetworkType>
    class GrainBoundaryTransmission
    //: private std::deque<TransmitSegment<typename DislocationNetworkType::LinkType> >
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef GrainBoundary<DislocationNetworkType> GrainBoundaryType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        DislocationNetworkType& DN;
        
        static constexpr double chordTol=1.0;
        
        
        typedef std::tuple<size_t,  // sourceID
        /*              */ size_t,  // sinkID
        /*              */ size_t,  // grainID
        /*              */ size_t,   // slipSystemID
        /*              */ std::pair<int,int>   // grainBoundaryID
        /*                     */ > TransmissionDataType;

        typedef std::map<double,TransmissionDataType> LinkTransmissionDataContainerType;
        typedef std::deque<TransmissionDataType> TransmissionDataContainerType;
    public:
        
        static size_t grainBoundaryTransmissionModel;
        
        /**********************************************************************/
        GrainBoundaryTransmission(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        {

        }
        
        /**********************************************************************/
        size_t directTransmit()
        {
            size_t nTransmitted=0;
            if(grainBoundaryTransmissionModel)
            {
                model::cout<<"		GrainBoundaryTransmission... "<<std::flush;
                const auto t0=std::chrono::system_clock::now();
                
                
                TransmissionDataContainerType transmissionDeq;
                
                for (const auto& link : DN.links() )
                {
                    LinkTransmissionDataContainerType transmissionMap;

                    
                    const VectorDim chord(link.second->chord());
                    
                    if(   link.second->isGrainBoundarySegment()
                       && chord.norm()>chordTol*DislocationNetworkRemesh<DislocationNetworkType>::Lmin
                       )
                    {
                        
                        const VectorDim unitChord(chord.normalized());
                        
                        
                        
                        for(const auto& grainBoundary : link.second->grainBoundaries())
                        {
                            for(const auto& grain : grainBoundary->grains())
                            {
                                for(size_t ssID=0;ssID<grain.second->slipSystems().size();++ssID)
                                {
                                    const auto& slipSystem(grain.second->slipSystems()[ssID]);
                                    if(fabs(slipSystem.n.cartesian().normalized().dot(unitChord))<FLT_EPSILON)
                                    {
                                        std::cout<<"GBsegment "<<link.second->source->sID<<"->"<<link.second->sink->sID<<" "<<grain.second->grainID<<" "<<ssID<<std::flush;
                                        double dG=0.0;
                                        switch (grainBoundaryTransmissionModel)
                                        {
                                            case 1:
                                            {
                                                dG=GrainBoundaryTransmissionEnergyModel<1>::dG(*link.second,
                                                                                               grainBoundary->glidePlane(grain.second->grainID).unitNormal,
                                                                                               slipSystem.s.cartesian(),
                                                                                               slipSystem.n.cartesian().normalized()
                                                                                               );
                                                
                                                break;
                                            }
                                                
                                            default:
                                                //                                                std::cout<"GrainBonudaryTransmissionModel not implemented."<<std::endl;
                                                break;
                                        }
                                        
                                        std::cout<<" dG="<<dG<<std::endl;

                                        
                                        if(dG<0.0)
                                        {
                                            const TransmissionDataType data(link.second->source->sID,
                                                                            link.second->  sink->sID,
                                                                            grain.second->grainID,
                                                                            ssID,
                                                                            grainBoundary->grainBndID);

                                            transmissionMap.emplace(dG,data);

                                        }
                                    }
                                }
                            }
                        }
                        
                    }
                    
                    if(transmissionMap.size())
                    {
                        transmissionDeq.push_back(transmissionMap.begin()->second); // best transmission choice for the segment
                    }
                    
                }
                
                
                
                // Insert new loops
                for(const auto& tup : transmissionDeq)
                {
                    const size_t& sourceID(std::get<0>(tup));
                    const size_t& sinkID(std::get<1>(tup));
                    const size_t& grainID(std::get<2>(tup));
                    const size_t& slipID(std::get<3>(tup));
                    const std::pair<int,int>& gbID(std::get<4>(tup));
                    
                    const auto isLink(DN.link(sourceID,sinkID));
                    if(isLink.first)
                    {
                        
                        std::cout<<"Transmitting"<<std::endl;
                        
                        
                        const VectorDim gbInNormal(-DN.poly.grainBoundary(gbID.first,gbID.second).glidePlane(grainID).unitNormal);
                        const VectorDim dir=(gbInNormal-gbInNormal.dot(DN.poly.grain(grainID).slipSystems()[slipID].n.cartesian().normalized())*DN.poly.grain(grainID).slipSystems()[slipID].n.cartesian().normalized()).normalized();
                        
                        const VectorDim newNodeP(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P())+dir*10.0);
                        const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
                        
                        std::vector<size_t> nodeIDs;
                        
                        nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
                        nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
                        nodeIDs.push_back(newNodeID);
                        
                        DN.insertLoop(nodeIDs,
                                      DN.poly.grain(grainID).slipSystems()[slipID].s.cartesian(),
                                      DN.poly.grain(grainID).slipSystems()[slipID].n.cartesian(),
                                      0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()),
                                      grainID);

                        nTransmitted++;
                    }
                    
                    
                    
                }
                
                DN.clearDanglingNodes();
                
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            return nTransmitted;
        }
        
//        /**********************************************************************/
//        void transmit()
//        {
//            if(grainBoundaryTransmission)
//            {
//                model::cout<<"		grain-boundary transmission..."<<std::flush;
//                
//                const auto t0= std::chrono::system_clock::now();
//                size_t nTransmitted=0;
//                for(const auto& gbSegment : *this)
//                {
//                    if(DN.link(gbSegment.sourceID,gbSegment.sinkID).first)
//                    {
//                        if(gbSegment.isValidTransmission)
//                        {
//                            const int originalGrainID=DN.link(gbSegment.sourceID,gbSegment.sinkID).second->grain.grainID;
//                            const int      newGrainID=gbSegment.transmitGrain->grainID;
//                            
//                            std::cout<<"segment GB ID="<<DN.link(gbSegment.sourceID,gbSegment.sinkID).second->grain.grainID<<std::endl;
//                            std::cout<<"transmit GB ID="<<gbSegment.transmitGrain->grainID<<std::endl;
//                            
//                            size_t newSourceID=gbSegment.sourceID;
//                            size_t newSinkID=gbSegment.sinkID;
//                            
//                            if(originalGrainID==newGrainID)
//                            {
//                                DN.template disconnect<false>(newSourceID,newSinkID);
//                                DN.connect(newSourceID,newSinkID,gbSegment.originalBurgers+*gbSegment.transmitBurgers);
//                            }
//                            else
//                            {
//                                newSourceID=DN.insertVertex(gbSegment.transmitSourceP,gbSegment.transmitGrain->grainID).first->first;
//                                newSinkID=DN.insertVertex(gbSegment.transmitSinkP,gbSegment.transmitGrain->grainID).first->first;
//                                DN.connect(newSourceID,newSinkID,*gbSegment.transmitBurgers);
//                            }
//                            
//                            const size_t newMidpointID(DN.insertVertex(gbSegment.transmitMidpoint->cartesian(),gbSegment.transmitGrain->grainID).first->first);
//                            DN.connect(newSinkID,newMidpointID,*gbSegment.transmitBurgers);
//                            DN.connect(newMidpointID,newSourceID,*gbSegment.transmitBurgers);
//                            
//                            nTransmitted++;
//                        }
//                    }
//                }
//                
//                
//                
//                model::cout<<"("<<nTransmitted<<" transmissions)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//                
//            }
//        }
        
    };
    
    template <typename DislocationNetworkType>
    size_t GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel=0;
    
    
} // end namespace
#endif

