/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DISLOCATIONCROSSSLIP_H
#define model_DISLOCATIONCROSSSLIP_H

#include <deque>
#include <utility>
#include <math.h>

#include <model/DislocationDynamics/CrossSlip/CrossSlipModels.h>
#include <model/DislocationDynamics/CrossSlip/CrossSlipSegment.h>

//#include <model/BVP/SearchData.h>


namespace model
{
    
    
    /*!\brief Class template which performs topological operations on a
     * DislocationNetwork as a consequence of cross-slip events.
     *
     * \tparam DislocationNetworkType the type of DislocationNetwork
     */
    template <typename DislocationNetworkType>
    class DislocationCrossSlip : private std::deque<CrossSlipSegment<typename DislocationNetworkType::LinkType> >
    {
        
        typedef typename DislocationNetworkType::NetworkNodeContainerType NetworkNodeContainerType;
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename NodeType::LatticePlaneContainerType LatticePlaneContainerType;
        
        typedef CrossSlipSegment<LinkType> CrossSlipSegmentType;
        typedef std::deque<CrossSlipSegmentType> CrossSlipContainerType;
        typedef typename DislocationNetworkType::GlidePlaneObserverType GlidePlaneObserverType;
        typedef typename DislocationNetworkType::VectorDimD VectorDimD;
        typedef std::pair<VectorDimD,VectorDimD> CrossSlipPairType;
        
        // A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        enum{dim=3};
        typedef LatticeVector<dim> LatticeVectorType;
        
        
        static double crossSlipDeg;
        static double crossSlipLength;
        
        /* Constructor *******************************************************/
        DislocationCrossSlip(DislocationNetworkType& dislocationNetwork_in) :
        /* init list */ DN(dislocationNetwork_in)
        {/* @param[in] dislocationNetwork_in A reference to the DislocationNetwork
          * Constructor finds and stores DislocationSegments that must cross-slip
          */
            const double sinCrossSlipRad(std::sin(crossSlipDeg*M_PI/180.0));
            
            //! 1-Loop over DislocationSegment(s), check cross-slip criterion and store CrossSlipSegment(s)
            for (typename NetworkLinkContainerType::iterator linkIter =DN.linkBegin();
                 /*                                       */ linkIter!=DN.linkEnd();
                 /*                                       */ linkIter++)
            {
                
                if ( !linkIter->second.source->isBoundaryNode() && !linkIter->second.sink->isBoundaryNode()
                    && linkIter->second.chord().normalized().cross(linkIter->second.Burgers.normalized()).norm()<=sinCrossSlipRad
                    && !linkIter->second.isSessile
                    && linkIter->second.chord().norm()>2.5*DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
                {
                    const LatticePlaneBase& conjugatePlaneBase=CrossSlipModels::maxRSS(linkIter->second);
                    if(conjugatePlaneBase != linkIter->second.glidePlane.n)
                    {
                        this->emplace_back(linkIter->second,conjugatePlaneBase);
                    }
                }
            }
        }
        
        /* crossSlip *******************************************************/
        size_t crossSlip()
        {
            
            size_t  n_crossSlips(0);
            
            //! 2- Loop over container of CrossSlipSegment(s) and perform cross-slip
            for (const auto& css : *this)
            {

                const size_t sourceCP=css.source.confiningPlanes().size();
                const size_t   sinkCP=css.  sink.confiningPlanes().size();
                
                LatticeVectorType sourceL=css.source.get_L();
                LatticeVectorType   sinkL=css.  sink.get_L();
                
                bool nodesOk=true;
                
                // Re-align nodes to perfect screw direction
                if(sourceCP==1 && sinkCP==1)
                {
                    //std::cout<<"CrossSlip case 1"<<std::endl;
                    // 2.1- correct midPoint if there is an existing conjugate plane close to it
                    
                    //                    const double hP(iterCS->midPoint.dot(normalConjugate)); // heigth of midPoint along normalConjugate
                    //                    for (typename GlidePlaneObserverType::const_iterator gpIter=GlidePlaneObserverType::begin();
                    //                         gpIter!=GlidePlaneObserverType::end();
                    //                         ++gpIter)
                    //                    {
                    //                        if((gpIter->second->reciprocalNormal-iterCS->reciprocalConjugateNormal).squaredNorm()==0)
                    //                        {
                    //                            const double num(gpIter->second->height-hP);
                    //
                    //                            if(std::fabs(num)<planeTol)
                    //                            {
                    //                                const double den(1.0-std::pow(iterCS->normalPrimary.dot(normalConjugate),2));
                    //                                const double u(num/den);
                    //                                if(std::fabs(u)<planeTol)
                    //                                {
                    //                                    midPoint += u*(normalConjugate-(normalConjugate.dot(iterCS->normalPrimary)*iterCS->normalPrimary));
                    //                                    break;
                    //                                }
                    //                            }
                    //                        }
                    //                    }

                    
                    const LatticeLine line(css.midPoint,css.Burgers);
                    sourceL=LatticeVectorType(line.snapToLattice(css.source.get_P()));
                    sinkL=LatticeVectorType(line.snapToLattice(css.  sink.get_P()));
                }
                else if(sourceCP>1 && sinkCP==1)
                {
                    //std::cout<<"CrossSlip case 2"<<std::endl;
                    const LatticeLine line(css.source.get_L(),css.Burgers);
                    sinkL=LatticeVectorType(line.snapToLattice(css.sink.get_P()));
                }
                else if(sourceCP==1 && sinkCP>1)
                {
                    //std::cout<<"CrossSlip case 3"<<std::endl;
                    const LatticeLine line(css.sink.get_L(),css.Burgers);
                    sourceL=LatticeVectorType(line.snapToLattice(css.source.get_P()));
                    
                }
                else // if nodes cannot be moved, make sure that they are in perfect screw direction
                {
                                        //std::cout<<"CrossSlip case 4"<<std::endl;
                    if((css.source.get_L()-css.sink.get_L()).cross(css.Burgers).squaredNorm()>0)
                    {// points perfectly aligned
                        nodesOk=false;
                    }
                }
                
                // Compute position of the new point on the conjugate plane
                const VectorDimD midPoint=sourceL.cartesian()*0.5+sinkL.cartesian()*0.5;
                const VectorDimD bc(css.Burgers.cartesian());
                const VectorDimD dir=bc.cross(css.conjugatePlaneBase.cartesian()).normalized();
                const double dirDotPK(dir.dot(css.pkForce));
                const double sgnDir((dirDotPK > 0.0) ? 1.0 : ((dirDotPK < 0.0) ? -1.0 : 0.0));
                const LatticePlane conjugatePlane(sourceL,css.conjugatePlaneBase);
                const VectorDimD conjugatePoint=conjugatePlane.snapToLattice(midPoint+sgnDir*dir*(css.source.get_V()*0.5+css.sink.get_V()*0.5).norm()*DN.get_dt());
                const LatticeVectorType conjugateL(conjugatePoint);
                const VectorDimD crossSlipVelocity((conjugatePoint-midPoint)/DN.get_dt());

                
                //std::cout<<(conjugateL-sourceL).dot(css.conjugatePlaneBase)<<std::endl;
                //std::cout<<(sinkL-conjugateL).dot(css.conjugatePlaneBase)<<std::endl;
                
                // Make sure that in the new position nodes are not coincident
                nodesOk*=((sourceL-sinkL).squaredNorm()>0);
                nodesOk*=((sourceL-conjugateL).squaredNorm()>0);
                nodesOk*=((conjugateL-sinkL).squaredNorm()>0);
                
                // Make sure that in new positions are inside mesh
                if (DN.shared.use_boundary)
                {
                    nodesOk*=DN.shared.mesh.isStrictlyInsideMesh(sourceL.cartesian(), css.source.includingSimplex(),FLT_EPSILON).first;
                    if (nodesOk)
                    {
                        nodesOk*=DN.shared.mesh.isStrictlyInsideMesh(sinkL.cartesian(),css.sink.includingSimplex(),FLT_EPSILON).first;
                        if (nodesOk)
                        {
                            nodesOk*=DN.shared.mesh.isStrictlyInsideMesh(conjugatePoint,css.source.includingSimplex(),FLT_EPSILON).first;
                        }
                    }
                    
                }
                
                if(nodesOk)
                {
                    // actually move nodes
                    css.source.set(sourceL);
                    css.  sink.set(sinkL);
                    
                    // expand
                    std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(css.source.sID,css.sink.sID,conjugateL,crossSlipVelocity);
                    assert(temp.second && "COULD NOT DO THIRD EXPANSION IN CROSS SLIP");
                     n_crossSlips++;
                }
                
            } // end for loop
            return  n_crossSlips;
        }
        
    };
    
    // Static data
    template <typename DislocationNetworkType>
    double DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg=2.0;
    
    //    template <typename DislocationNetworkType>
    //    double DislocationCrossSlip<DislocationNetworkType>::crossSlipLength=100.0;
    
} // namespace model
#endif



//                    const VectorDimD bc(iterCS->Burgers.cartesian());
//                    const double dL=round(0.5*crossSlipLength/bc.norm());
//
//                    // 2.2- Compute position of cross-slip points on the common line
//                    CrossSlipPairType crossPoints((bc.dot(iterCS->chord)>=0.0) ? CrossSlipPairType(midPoint-bc*dL,midPoint+bc*dL)
//                                                  /*                                                      */ : CrossSlipPairType(midPoint+bc*dL,midPoint-bc*dL));
//
//
//                    // 2.3- Compute position of the new point on the conjugate plane
//                    VectorDimD dir(bc.cross(iterCS->conjugatePlane.n.cartesian()).normalized());
//                    double dirDotPK(dir.dot(iterCS->pkForce));
//                    double sgnDir((dirDotPK > 0.0) ? 1.0 : ((dirDotPK < 0.0) ? -1.0 : 0.0));
//                    VectorDimD crossSlipDisplacement(iterCS->conjugatePlane.n.snapToLattice(sgnDir*dir*conjugatePointDistance));
//                    VectorDimD crossSlipVelocity(crossSlipDisplacement/DN.get_dt());
//                    //                    VectorDimD conjugatePoint(0.5*(crossPoints.first+crossPoints.second)+crossSlipDisplacement);
//                    VectorDimD conjugatePoint(midPoint+crossSlipDisplacement);
//
//                    bool crossSlipPointsInsideMesh(true);
//                    if (DN.shared.use_boundary)
//                    {
//                        crossSlipPointsInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(conjugatePoint,
//                                                                                       DN.node(iterCS->source.sID).second->includingSimplex(),
//                                                                                       FLT_EPSILON).first;
//                        if (crossSlipPointsInsideMesh)
//                        {
//                            crossSlipPointsInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(crossPoints.first,
//                                                                                           DN.node(iterCS->source.sID).second->includingSimplex(),
//                                                                                           FLT_EPSILON).first;
//                            if(crossSlipPointsInsideMesh)
//                            {
//                                crossSlipPointsInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(crossPoints.second,
//                                                                                               DN.node(iterCS->source.sID).second->includingSimplex(),
//                                                                                               FLT_EPSILON).first;
//                            }
//                        }
//                    }
//
//                    if (crossSlipPointsInsideMesh)
//                    {
//
//                        //                        //std::cout<<"dL="<<dL<<std::endl;
//                        //                        //std::cout<<"crossPoints.first"<<crossPoints.first.transpose()<<std::endl;
//                        //                        //std::cout<<"crossPoints.second"<<crossPoints.second.transpose()<<std::endl;
//                        //                        //std::cout<<"conjugate point"<<conjugatePoint.transpose()<<std::endl;
//
//                        // 2.4- call expand
//
//                        std::pair<typename NetworkNodeContainerType::iterator,bool> expand1(DN.expand(iterCS->source.sID,iterCS->sink.sID,LatticeVectorType(crossPoints.first))); // place first point on common line
//                        assert(expand1.second && "COULD NOT DO FIRST EXPANSION IN CROSS SLIP");
//
//                        std::pair<typename NetworkNodeContainerType::iterator,bool> expand2(DN.expand(expand1.first->first,iterCS->sink.sID,LatticeVectorType(crossPoints.second)));  // place second point on common line
//                        assert(expand2.second && "COULD NOT DO SECOND EXPANSION IN CROSS SLIP");
//
//                        std::pair<typename NetworkNodeContainerType::iterator,bool> expand3(DN.expand(expand1.first->first,expand2.first->first,LatticeVectorType(conjugatePoint),crossSlipVelocity));
//                        assert(expand3.second && "COULD NOT DO THIRD EXPANSION IN CROSS SLIP");
//
//                         n_crossSlips++;
//                    }


//! 1-Loop over DislocationSegment(s), check cross-slip criterion and store CrossSlipSegment(s)
//			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
//            {
////				CrossSlipSegmentType css(linkIter->second.isCrossSlipSegment(sinCrossSlipRad,crossSlipLength));
////				if(css.isCrossSlipSegment)
////                {
////					CSC.push_back(css);
////				}
//                CSC.emplace_back(linkIter->second.isCrossSlipSegment(sinCrossSlipRad,crossSlipLength));
//
//			}

