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

#include <vector>
#include <utility>
#include <math.h>

#include <model/DislocationDynamics/CrossSlip/CrossSlipSegment.h>
//#include <model/BVP/SearchData.h>


namespace model {
    
    
	
	template <typename DislocationNetworkType>
	class DislocationCrossSlip
    {
		
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef CrossSlipSegment<LinkType> CrossSlipSegmentType;
        typedef std::vector<CrossSlipSegmentType> CrossSlipContainerType;
        typedef typename DislocationNetworkType::GlidePlaneObserverType GlidePlaneObserverType;	
        typedef typename DislocationNetworkType::VectorDimD VectorDimD;
        
        
        CrossSlipContainerType CSC;
        
        GlidePlaneObserverType gpo;
        
        // A reference to the DislocationNetwork
        DislocationNetworkType& dislocationNetwork; 
        
    public:
		
		enum{dim=3};
		
		
		/* Constructor *******************************************************/		
		DislocationCrossSlip(DislocationNetworkType& dislocationNetwork_in) :
        /* init list */ dislocationNetwork(dislocationNetwork_in)
        {/* @param[in] dislocationNetwork_in A reference to the DislocationNetwork
          * Constructor Initializes
          */
        }
		
		/* crossSlip *******************************************************/		
		size_t crossSlip(const double& crossSlipDeg,const double& crossSlipLength)
        {
            
            const double sinCrossSlipRad(std::sin(crossSlipDeg*M_PI/180.0));
            const double planeTol(5.0);
            const double conjugatePointDistance(0.01*crossSlipLength);
			
			//! 1-Loop over DislocationSegment(s), check cross-slip criterion and store CrossSlipSegment(s)
			for (typename NetworkLinkContainerType::const_iterator linkIter=dislocationNetwork.linkBegin();linkIter!=dislocationNetwork.linkEnd();++linkIter)
            {
				CrossSlipSegmentType css(linkIter->second->isCrossSlipSegment(sinCrossSlipRad,crossSlipLength));
				if(css.isCrossSlipSegment)
                {
					CSC.push_back(css);
				}
			}
			
			//! 2- Loop over container of CrossSlipSegment(s) and perform cross-slip
			for (typename CrossSlipContainerType::const_iterator iterCS=CSC.begin();iterCS!=CSC.end();++iterCS)
            {
				VectorDimD midPoint(iterCS->midPoint);
				const double hP(iterCS->midPoint.dot(iterCS->normalConjugate)); // heigth of midPoint along normalConjugate
				
				// 2.1- correct midPoint if there is an existing conjugate plane close to it
				for (typename GlidePlaneObserverType::const_iterator gpIter=gpo.begin(); gpIter!=gpo.end();++gpIter)
                {
					if((gpIter->second->planeNormal-iterCS->normalConjugate).norm()<FLT_EPSILON)
                    {
                        const double num(gpIter->second->height-hP);

                        if(std::fabs(num)<planeTol)
                        {
							const double den(1.0-std::pow(iterCS->normalPrimary.dot(iterCS->normalConjugate),2));
							const double u(num/den);
                            if(std::fabs(u)<planeTol)
                            {
                                midPoint += u*(iterCS->normalConjugate-(iterCS->normalConjugate.dot(iterCS->normalPrimary)*iterCS->normalPrimary));
                                break;
                            }
						}
					}
				}
				
				// 2.2- Compute position of cross-slip points on the common line
				typedef std::pair<VectorDimD,VectorDimD> CrossSlipPairType;
				CrossSlipPairType crossPoints((iterCS->Burgers.dot(iterCS->chord)>=0.0) ? CrossSlipPairType(midPoint-iterCS->Burgers.normalized()*crossSlipLength*0.5,midPoint+iterCS->Burgers.normalized()*crossSlipLength*0.5)
                /*                                                                   */ : CrossSlipPairType(midPoint+iterCS->Burgers.normalized()*crossSlipLength*0.5,midPoint-iterCS->Burgers.normalized()*crossSlipLength*0.5));
				
 				// 2.3- Compute position of the new point on the conjugate plane
				VectorDimD dir(iterCS->Burgers.cross(iterCS->normalConjugate).normalized());
				double dirDotPK(dir.dot(iterCS->pkForce));
				double sgnDir((dirDotPK > 0.0) ? 1.0 : ((dirDotPK < 0.0) ? -1.0 : 0.0));
				VectorDimD crossSlipDisplacement(sgnDir*dir*conjugatePointDistance);		
				VectorDimD crossSlipVelocity(crossSlipDisplacement/dislocationNetwork.get_dt());		
				VectorDimD conjugatePoint(0.5*(crossPoints.first+crossPoints.second)+crossSlipDisplacement);
				
				bool crossSlipPointsInsideMesh(true);
				if (dislocationNetwork.shared.use_boundary)
                {
//					SearchData<dim> SD(conjugatePoint);
//                    //                    dislocationNetwork.shared.domain.findIncludingTet(SD);
//                    dislocationNetwork.shared.domain.findIncludingTet(SD,dislocationNetwork.node(iterCS->sourceID).second->meshID());
//					crossSlipPointsInsideMesh*=(SD.meshLocation()==insideMesh);
					
                    crossSlipPointsInsideMesh*=dislocationNetwork.shared.mesh.isStrictlyInsideMesh(conjugatePoint,
                                                                                                   dislocationNetwork.node(iterCS->sourceID).second->includingSimplex(),
                                                                                                   FLT_EPSILON).first;
                    
                    if (crossSlipPointsInsideMesh)
                    {
//						SearchData<dim> SD1(crossPoints.first);
//						//dislocationNetwork.shared.domain.findIncludingTet(SD1);
//                        dislocationNetwork.shared.domain.findIncludingTet(SD1,dislocationNetwork.node(iterCS->sourceID).second->meshID());
//						crossSlipPointsInsideMesh*=(SD1.meshLocation()==insideMesh);
                        
                        crossSlipPointsInsideMesh*=dislocationNetwork.shared.mesh.isStrictlyInsideMesh(crossPoints.first,
                                                                                                       dislocationNetwork.node(iterCS->sourceID).second->includingSimplex(),
                                                                                                       FLT_EPSILON).first;

                        
						if(crossSlipPointsInsideMesh)
                        {
//							SearchData<dim> SD2(crossPoints.second);
//                            //							dislocationNetwork.shared.domain.findIncludingTet(SD2);
//                            dislocationNetwork.shared.domain.findIncludingTet(SD2,dislocationNetwork.node(iterCS->sourceID).second->meshID());
//							crossSlipPointsInsideMesh*=(SD2.meshLocation()==insideMesh);
                            
                            crossSlipPointsInsideMesh*=dislocationNetwork.shared.mesh.isStrictlyInsideMesh(crossPoints.second,
                                                                                                           dislocationNetwork.node(iterCS->sourceID).second->includingSimplex(),
                                                                                                           FLT_EPSILON).first;

						}
					}
				}
                
				if (crossSlipPointsInsideMesh)
                {
                    std::cout<<" cross-slipping "<<std::endl;
					// 2.4- call expand
					std::pair<bool,size_t> expand1(dislocationNetwork.expand(iterCS->sourceID,iterCS->sinkID,crossPoints.first)); // place first point on common line
                    assert(expand1.first && "COULD NOT DO FIRST EXPANSION IN CROSS SLIP");
                    
					std::pair<bool,size_t> expand2(dislocationNetwork.expand(expand1.second,iterCS->sinkID,crossPoints.second));  // place second point on common line
					assert(expand2.first && "COULD NOT DO SECOND EXPANSION IN CROSS SLIP");
                    
					std::pair<bool,size_t> expand3(dislocationNetwork.expand(expand1.second,expand2.second,conjugatePoint,crossSlipVelocity)); 
					assert(expand3.first && "COULD NOT DO THIRD EXPANSION IN CROSS SLIP");
				}
                
			}
			return CSC.size();
		}
        
	};
	
	
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

