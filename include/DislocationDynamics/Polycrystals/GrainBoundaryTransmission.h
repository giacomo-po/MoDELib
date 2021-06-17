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
#include <TerminalColors.h>
//#include <TransmitSegment.h>
 // defines mode::cout


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
        static double dGdirect(const SegmentType& seg,
                               const Eigen::Matrix<double,SegmentType::dim,1>& grain2OutNormal,
                               const Eigen::Matrix<double,SegmentType::dim,1>& b2,
                               const Eigen::Matrix<double,SegmentType::dim,1>& n2)
        {
            
            double temp=0.0;
            
            const size_t q=seg.quadraturePoints().size()/2;
            if(seg.quadraturePoints().size()>q)
            {
                const Eigen::Matrix<double,SegmentType::dim,1> pk((seg.quadraturePoint(q).stress*b2).cross(seg.chord().normalized()));
                
                //		 std::cout<<"here 1"<<std::endl;
                
                const Eigen::Matrix<double,SegmentType::dim,1> pkg(pk-pk.dot(n2)*n2);
                
                if(   pkg.dot(grain2OutNormal)<0.0
                   )
                {// pk force on segment must be pushing into new grain
                    const Eigen::Matrix<double,SegmentType::dim,1>& b1(seg.burgers());
                    const double tau=fabs((seg.quadraturePoint(q).stress*b2.normalized()).dot(n2));
                    
                    temp=0.5*alpha*((b1-b2).squaredNorm()+b2.squaredNorm()-b1.squaredNorm())-tau*b2.norm()*lambda;
                }
            }
            
            //		 std::cout<<"here 0"<<std::endl;
            //		 std::cout<<"seg.quadratureParticleContainer.size()="<<seg.quadratureParticleContainer.size()<<std::endl;
            //            if((seg.stressAtQuadrature(q)*seg.burgers()).cross(seg.chord().normalized()).dot(grain2OutNormal)<0.0)
            
            return temp;
        }
        
        
        template<typename SegmentType>
        static double dGindirect(const SegmentType& seg,
                                 const Eigen::Matrix<double,SegmentType::dim,1>& grain2OutNormal,
                                 const Eigen::Matrix<double,SegmentType::dim,1>& b2,
                                 const Eigen::Matrix<double,SegmentType::dim,1>& n2)
        {
            
            return dGdirect(seg,grain2OutNormal,b2,n2);
        }
        
    };
    
    template <typename DislocationNetworkType>
    class GrainBoundaryTransmission
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef Grain<dim> GrainType;
        typedef GrainBoundary<dim> GrainBoundaryType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        DislocationNetworkType& DN;
        
        static constexpr double chordTol=1.0;
        
        typedef std::tuple<size_t,  // sourceID
        /*              */ size_t,  // sinkID
        /*              */ size_t,  // grainID
        /*              */ size_t,   // slipSystemID
        /*              */ std::pair<int,int>,   // grainBoundaryID
        /*              */ int     // transmission type
        /*                     */ > TransmissionDataType;
        
        typedef std::map<double,TransmissionDataType> LinkTransmissionDataContainerType;
        typedef std::deque<TransmissionDataType> TransmissionDataContainerType;
        
        
        /**********************************************************************/
        void emplaceData(const LinkType& seg,
                                const GrainBoundaryType& grainBoundary,
                                const GrainType& grain,
                                const std::shared_ptr<SlipSystem>& slipSystem,
                                const size_t& ssID,
                                const int& transmissionType,
                                LinkTransmissionDataContainerType& transmissionMap) const
        {
            
            std::cout<<"GBsegment "<<seg.source->sID<<"->"<<seg.sink->sID<<" "<<grain.grainID<<" "<<ssID<<std::flush;

            
            double dG=0.0;
            switch (grainBoundaryTransmissionModel)
            {
                case 1:
                {
                    switch (transmissionType)
                    {
                        case 0: // direct transmit
                        {
                            std::cout<<" directTransmit "<<std::flush;
                            dG=GrainBoundaryTransmissionEnergyModel<1>::dGdirect(seg,
                                                                                 grainBoundary.outNormal(grain.grainID),
                                                                                 slipSystem->s.cartesian(),
                                                                                 slipSystem->n.cartesian().normalized()
                                                                                 );
                            break;
                        }
                            
                        case 1: // offset transmit
                        {
                            std::cout<<" offsetTransmit "<<std::flush;
                            //                    dG=GrainBoundaryTransmissionEnergyModel<1>::dGdirect(*link.second,
                            //                                                                         grainBoundary->glidePlane(grain.second->grainID).unitNormal,
                            //                                                                         slipSystem->s.cartesian(),
                            //                                                                         slipSystem->n.cartesian().normalized()
                            //                                                                         );
                            break;
                        }
                            
                        default:
                        {
                            std::cout<<" indirectTransmit "<<std::flush;
                            dG=GrainBoundaryTransmissionEnergyModel<1>::dGindirect(seg,
                                                                                 grainBoundary.outNormal(grain.grainID),
                                                                                 slipSystem->s.cartesian(),
                                                                                 slipSystem->n.cartesian().normalized()
                                                                                 );
                            break;
                        }
                    }
                    break;
                }
                    
                default:
                    std::cout<<" Transmission Model not implemented. "<<std::flush;
                
                    break;
            }
            
            std::cout<<"dG= "<<dG<<std::endl;
            
            
            if(dG<0.0)
            {
                const TransmissionDataType data(seg.source->sID,
                                                seg.sink->sID,
                                                grain.grainID,
                                                ssID,
                                                grainBoundary.grainBndID,
                                                transmissionType);
                
                transmissionMap.emplace(dG,data);
                
            }
            
        }
        
    public:
        
        const size_t grainBoundaryTransmissionModel;
        
        /**********************************************************************/
        GrainBoundaryTransmission(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,grainBoundaryTransmissionModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("grainBoundaryTransmissionModel",true))
        {
            
        }
        
        
        
        /**********************************************************************/
        size_t directTransmit()
        {
            assert(false && "REWORK THIS");
            return 0;
        }
//        {
//            size_t nTransmitted=0;
//            if(grainBoundaryTransmissionModel)
//            {
//                std::cout<<"        GrainBounTransmission... "<<std::flush;
//                const auto t0=std::chrono::system_clock::now();
//
//
//                TransmissionDataContainerType transmissionDeq;
//
//                for (const auto& link : DN.links() )
//                {
//                    LinkTransmissionDataContainerType transmissionMap;
//
//
//                    const VectorDim chord(link.second->chord());
//                    const VectorDim glidePlaneNormal(link.second->glidePlaneNormal());
//
//                    if(   link.second->isGrainBoundarySegment()
//                       && chord.norm()>chordTol*DislocationNetworkRemesh<DislocationNetworkType>::Lmin
//                       )
//                    {
//                        const VectorDim unitChord(chord.normalized());
//
//                        for(const auto& grainBoundary : link.second->grainBoundaries())
//                        {
//                            for(const auto& grain : grainBoundary->grains())
//                            {
//                                for(size_t ssID=0;ssID<grain.second->slipSystems().size();++ssID)
//                                {
//                                    const auto& slipSystem(grain.second->slipSystems()[ssID]);
////                                    if(grainBoundary->glidePlane(grain.second->grainID).unitNormal.cross(slipSystem->n.cartesian().normalized()).norm()>FLT_EPSILON) // slip plane is not GB plane
//                                        if(grainBoundary->outNormal(grain.second->grainID).cross(slipSystem->n.cartesian().normalized()).norm()>FLT_EPSILON // slip plane is not GB plane
//                                           //&& slipSystem->n.cartesian().normalized().cross(glidePlaneNormal).norm()>FLT_EPSILON
//                                           ) // NOT SAME PLANE, REMOVE THIS
//
//                                    {
//                                        if(fabs(slipSystem->n.cartesian().normalized().dot(unitChord))<FLT_EPSILON )
//                                        {// slip plane cointains chord -> direct or indirect (offset) transmit
//                                            const std::pair<bool,long int> temp=LatticePlane::computeHeight(slipSystem->n,0.5*(link.second->source->get_P()+link.second->sink->get_P()));
//                                            if(temp.first)
//                                            {// direct transmit possible
//                                                emplaceData(*link.second,*grainBoundary,*grain.second,slipSystem,ssID,0,transmissionMap);
//                                            }
//                                            else
//                                            {// direct transmit not possible, use offset transmit
//                                                emplaceData(*link.second,*grainBoundary,*grain.second,slipSystem,ssID,1,transmissionMap);
//                                            }
//                                        }
//                                        else
//                                        {// slip plane doesn't cointains chord -> indirect (incident) transmit
//                                            emplaceData(*link.second,*grainBoundary,*grain.second,slipSystem,ssID,2,transmissionMap);
//                                        }
//                                    }
//                                }
//                            }
//                        }
//
//                    }
//
//                    if(transmissionMap.size())
//                    {
//                        transmissionDeq.push_back(transmissionMap.begin()->second); // best transmission choice for the segment
//                    }
//
//                }
//
//
//
//                // Insert new loops
//                for(const auto& tup : transmissionDeq)
//                {
//                    const size_t& sourceID(std::get<0>(tup));
//                    const size_t& sinkID(std::get<1>(tup));
//                    const size_t& grainID(std::get<2>(tup));
//                    const size_t& slipID(std::get<3>(tup));
//                    const std::pair<int,int>& gbID(std::get<4>(tup));
//                    const int& transmissionType(std::get<5>(tup));
//
//                    const auto isLink(DN.link(sourceID,sinkID));
//                    if(isLink.first)
//                    {
//                        switch(transmissionType)
//                        {
//                            case 0:
//                            {
//                                std::cout<<"direct transmission"<<std::endl;
//
//
//                                const VectorDim gbInNormal(-DN.poly.grainBoundary(gbID.first,gbID.second).outNormal(grainID));
//                                const VectorDim dir=(gbInNormal-gbInNormal.dot(DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().normalized())*DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().normalized()).normalized();
//
//                                const VectorDim newNodeP(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P())+dir*10.0);
//                                const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
//
//                                std::vector<size_t> nodeIDs;
//
//                                nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
//                                nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
//                                nodeIDs.push_back(newNodeID);
//
//                                DN.insertLoop(nodeIDs,
//                                              DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
//                                              DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian(),
//                                              0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()),
//                                              grainID);
//
//                                nTransmitted++;
//                                break;
//                            }
//                            case 1: // offset transmit
//                            {
//                                std::cout<<"indirect transmission"<<std::endl;
//
//
//                                const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(DN.poly.grain(grainID).slipSystems()[slipID]->n,
//                                                                                                      0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()));
//
//                                //                                const double planeSpacing=1.0/DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().norm();
//
//
//                                PlanePlaneIntersection<dim> ppi(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()),
//                                                                DN.poly.grainBoundary(gbID.first,gbID.second).outNormal(grainID),
//                                                                heightPair.second*DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian()/DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().squaredNorm(),
//                                                                DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian());
//
//                                assert(0 && "FINISH HERE");
//
//
//
//
//                                break;
//                            }
//                            default: // indirect (incident) transmit
//                            {
////                                std::cout<<"indirect transmission"<<std::endl;
////
////
////                                const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(DN.poly.grain(grainID).slipSystems()[slipID]->n,
////                                                                                                      0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()));
////
////                                const VectorDim newLoopP(heightPair.second*DN.poly.grain(grainID).slipSystems()[slipID]->n.interplaneVector());
////                                const VectorDim newLoopN(DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().normalized());
////
//////                                PlanePlaneIntersection<dim> ppi(0.5*(isLink.second->source->get_P()+isLink.second->sink->get_P()),
//////                                                                DN.poly.grainBoundary(gbID.first,gbID.second).outNormal(grainID),
//////                                                                newLoopP,
//////                                                                newLoopN);
////
////
////                                MeshPlane<dim> mp(DN.mesh,grainID,newLoopP,newLoopN);
////                                BoundingLineSegments<dim> bls(mp);
////
////                                const VectorDim newSourceP(std::get<0>(bls.snap(isLink.second->source->get_P())));
////                                const size_t newSourceID=DN.insertDanglingNode(newSourceP,VectorDim::Zero(),1.0).first->first;
////
////                                const VectorDim newSinkP(std::get<0>(bls.snap(isLink.second->sink->get_P())));
////                                const size_t newSinkID=DN.insertDanglingNode(newSinkP,VectorDim::Zero(),1.0).first->first;
////
////
//////                                assert(0 && "FINISH HERE");
////
////                                const VectorDim gbInNormal(-DN.poly.grainBoundary(gbID.first,gbID.second).outNormal(grainID));
////                                const VectorDim dir=(gbInNormal-gbInNormal.dot(DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().normalized())*DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian().normalized()).normalized();
////
////                                const VectorDim newNodeP(0.5*(newSourceP+newSinkP)+dir*10.0);
////                                const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
//////
////                                std::vector<size_t> nodeIDs;
////
////                                nodeIDs.push_back(newSinkID);      // insert in reverse order, sink first, source second
////                                nodeIDs.push_back(newSourceID);    // insert in reverse order, sink first, source second
////                                nodeIDs.push_back(newNodeID);
//////
////                                DN.insertLoop(nodeIDs,
////                                              DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
////                                              DN.poly.grain(grainID).slipSystems()[slipID]->n.cartesian(),
////                                              newLoopP,
////                                              grainID);
////
////                                nTransmitted++;
//
//
//
//                                break;
//                            }
//                        }
//
//
//                    }
//                }
//
//                DN.clearDanglingNodes();
//
//                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
//            return nTransmitted;
//        }
        
        
    };
    
//    template <typename DislocationNetworkType>
//    size_t GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel=0;
    
    
} // end namespace
#endif



//        /**********************************************************************/
//        void transmit()
//        {
//            if(grainBoundaryTransmission)
//            {
//                std::cout<<"		grain-boundary transmission..."<<std::flush;
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
//                std::cout<<"("<<nTransmitted<<" transmissions)"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//
//            }
//        }

