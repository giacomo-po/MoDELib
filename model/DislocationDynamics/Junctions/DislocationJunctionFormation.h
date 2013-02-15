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
#include <model/Network/Operations/EdgeFinder.h>

#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>
#include <model/BVP/SearchData.h>


//#include <touple.h>


namespace model {
	
	template <typename DislocationNetworkType>
	class DislocationJunctionFormation{
		
		typedef typename DislocationNetworkType::LinkType LinkType;
		typedef typename DislocationNetworkType::NodeType NodeType;

        typedef typename EdgeFinder<LinkType>::isNetworkEdgeType isNetworkLinkType;
		typedef typename DislocationNetworkType::isNetworkNodeType isNetworkNodeType;
		
		enum {dim=3};
        enum {pOrder=3};
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		
		typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
		typedef typename DislocationNetworkType::NetworkNodeContainerType NetworkNodeContainerType;
		
		typedef std::pair<const LinkType*, double> EdgeIntersectionType;
		typedef std::pair<EdgeIntersectionType,EdgeIntersectionType> EdgeIntersectionPairType;
		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
		
		
		
		// a reference to the dislocation network
		DislocationNetworkType& DN;
		
		//	std::touple<int,double> tp;
        
        
        
        
        
		
	public:
		
        
        static double collisionTol;

		
		/* Constructor ******************************/
		DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {}
		
		/* findIntersections **************************************************/
		EdgeIntersectionPairContainerType findIntersections(const double& avoidNodeIntersection) const
        {
			
			EdgeIntersectionPairContainerType intersectionContainer;
			
			
			//! 2- loop over all links and determine their intersections
			for (typename NetworkLinkContainerType::const_iterator linkIterA=DN.linkBegin();linkIterA!=DN.linkEnd();linkIterA++){
                
                
                
                const DislocationSegmentIntersection<dim,pOrder> dsi(linkIterA->second->hermiteCoefficients(),linkIterA->second->glidePlaneNormal);
                
				for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.linkEnd();linkIterB++){
					if (linkIterA->second->sID!=linkIterB->second->sID){	// don't intersect with itself
                        
        //                std::cout<<"intersecting "<<linkIterA->second->source->sID<<"->"<<linkIterA->second->sink->sID<<" and "<<linkIterB->second->source->sID<<"->"<<linkIterB->second->sink->sID<<std::endl;
                        
                        //                        const bool areIncidentAtNodes(   linkIterA->second->source->sID==linkIterB->second->source->sID 
                        //                                                      || linkIterA->second->source->sID==linkIterB->second->sink->sID 
                        //                                                      || linkIterA->second->  sink->sID==linkIterB->second->source->sID 
                        //                                                      || linkIterA->second->  sink->sID==linkIterB->second->sink->sID); 
                        //                        
                        //                       const bool areOnDifferentPlanes((linkIterA->second->glidePlaneNormal-linkIterB->second->glidePlaneNormal).squaredNorm()>FLT_EPSILON);
                        const bool L1isSessile(linkIterA->second->sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
                        const bool L2isSessile(linkIterB->second->sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
                        //                       const bool eitherOneIsSessile(L1isSessile || L2isSessile);            
                        //                       const bool areIncidentAndOnDifferentPlanes(areIncidentAtNodes && areOnDifferentPlanes && eitherOneIsSessile);
                        
                        
                        //                        const bool dontIntersect(areIncidentAndOnDifferentPlanes || bothAreSessile || L2NotOnL1Planes || L1NotOnL2Planes);
                        //                        const bool dontIntersect(areIncidentAndOnDifferentPlanes);
                        
                        
                        //                      if(!dontIntersect){
                        
                        
                        
                        
                        std::set<std::pair<double,double> > temp;
                        
                        
                        
                        if (!L1isSessile && !L2isSessile) { // both are glissile
                            temp = dsi.intersectWith(linkIterB->second->hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                        }
                        
                        else if (!L1isSessile && L2isSessile){ // L1 is glissile and L2 is sessile
                            const bool gnAgnB((linkIterA->second->glidePlaneNormal-linkIterB->second->glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                            const bool gnAsnB((linkIterA->second->glidePlaneNormal-linkIterB->second->sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                            
                            if (!gnAgnB && !gnAsnB){
                                // cannot intersect 
                            }
                            else if(gnAgnB && !gnAsnB){ // use planeNormal of A and planeNormal of B
                                temp = dsi.intersectWith(linkIterB->second->hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else if (!gnAgnB && gnAsnB){ // use planeNormal of A and sessileNormal of B
                                temp = dsi.intersectWith(linkIterB->second->hermiteCoefficients(),linkIterB->second->sessilePlaneNormal,collisionTol);
                            }
                            else{
                                assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                            }
                        }
                        else if (L1isSessile && !L2isSessile){ // L1 is sessile and L2 is glissile
                            const bool gnBgnA((linkIterB->second->glidePlaneNormal-linkIterA->second->glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                            const bool gnBsnA((linkIterB->second->glidePlaneNormal-linkIterA->second->sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                            
                            if (!gnBgnA && !gnBsnA){
                                // cannot intersect 
                            }
                            else if(gnBgnA && !gnBsnA){ // use planeNormal of A and planeNormal of B
                                temp = dsi.intersectWith(linkIterB->second->hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else if (!gnBgnA && gnBsnA){ // use sessileNormal of A and use planeNormal of B
                                const DislocationSegmentIntersection<dim,pOrder> dsi2(linkIterA->second->hermiteCoefficients(),linkIterA->second->sessilePlaneNormal);
                                temp = dsi2.intersectWith(linkIterB->second->hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else{
                                assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                            }
                        }
                        else { // both are sessile
                            // cannot intersect
                        }
                        
                        
                        
                        // std::set<std::pair<double,double> > temp ( dsi.intersectWith(linkIterB->second->hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol));
                        
                        
                        for (std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paramIter!=temp.end();++paramIter){					
                            if (   paramIter->first >avoidNodeIntersection && paramIter-> first<1.0-avoidNodeIntersection
                                && paramIter->second>avoidNodeIntersection && paramIter->second<1.0-avoidNodeIntersection){		// avoid node intersection						
                                EdgeIntersectionType intersectionOnA(std::make_pair(&(*linkIterA->second),paramIter->first ));
                                EdgeIntersectionType intersectionOnB(std::make_pair(&(*linkIterB->second),paramIter->second));
                                EdgeIntersectionPairType intersectionPair(std::make_pair(intersectionOnA,intersectionOnB));
                                intersectionContainer.push_back(intersectionPair);
                            }
                        }
                        //                       }
                    }   
				}
			}	
			return intersectionContainer;
		}
		
		
		
		
		
		
		
		/* formJunctions **************************************************/
		void formJunctions(const double& dx, const double& avoidNodeIntersection) {
			
			
			//! 1- Initialize intersectionContainer calling findIntersections
			EdgeIntersectionPairContainerType intersectionContainer(findIntersections(avoidNodeIntersection));
			std::cout<<intersectionContainer.size()<<" geometric) ";
			
			//! 2- Remove from intersectionContainer all intersections that don't satisfy Frank's rule
			std::vector<int> dirVector;			
			for (typename EdgeIntersectionPairContainerType::iterator iter=intersectionContainer.begin();iter!=intersectionContainer.end();){
				const double u1(iter-> first.second);
				const double u2(iter->second.second);
				const VectorDimD b1(iter-> first.first->flow);
				const VectorDimD b2(iter->second.first->flow);
				const VectorDimD rl1(iter-> first.first->get_rl(u1));
				const VectorDimD rl2(iter->second.first->get_rl(u2));
                
                
                const bool isIsessile(iter-> first.first->sessilePlaneNormal.norm()>FLT_EPSILON);
                const bool isJsessile(iter->second.first->sessilePlaneNormal.norm()>FLT_EPSILON);
                
                VectorDimD prjDir(VectorDimD::Zero());
                if (!isIsessile && !isJsessile){
                    const VectorDimD commonLine(iter-> first.first->glidePlaneNormal.cross(iter->second.first->glidePlaneNormal));
                    if(commonLine.norm()>FLT_EPSILON){ // planes are not parallel, intersection will be on common line
                        prjDir=commonLine.normalized();
                    }
                    else{
                        const double rl1Norm(rl1.norm()); 
                        assert(rl1Norm>FLT_EPSILON && "TANGENT HAS ZERO NORM");
                        prjDir=rl1/rl1Norm;
                    }
                }
                else if(isIsessile && !isJsessile){ // use chord of I
                    prjDir=iter-> first.first->chord().normalized();
                }
                else if(!isIsessile && isJsessile){ // use chord of J
                    prjDir=iter->second.first->chord().normalized();
                }
                else{
                    assert(0 && "CANNOT DETERMINE COMMON LINE BETWEEN TWO SESSILE SEGMENTS.");
                }
                
                
                const double sgnrl1rl2(rl1.dot(prjDir)*rl2.dot(prjDir));
                
                
				//const bool frankRule(b1.dot(b2)*rl1.dot(rl2)<=0.0);
				const bool frankRule(b1.dot(b2)*sgnrl1rl2<=0.0);
                
                
                //				const bool frankRule(true);
				if (!frankRule){
					iter=intersectionContainer.erase(iter);
				}
				else{ // determine the relative orientation of the links
                    //					const int dir((rl1.dot(rl2) > 0.0) ? 1 : ((rl1.dot(rl2) < 0.0) ? -1 : 0)); 
					const int dir(( sgnrl1rl2 > 0.0) ? 1 : ((sgnrl1rl2 < 0.0) ? -1 : 0)); 
					dirVector.push_back(dir);
					++iter; // increment iterator only if didn't erase
				}
			}
			std::cout<<intersectionContainer.size()<<" physical intersections. ";
			assert(intersectionContainer.size()==dirVector.size());
            
			//! 3- Organize intersections by segments to use multiexpand
			typedef size_t IntersectionIDType;
			typedef std::map<double, IntersectionIDType> MultiExpandInputType;
			typedef std::pair<size_t,size_t> EdgeIDType;
			typedef std::map<EdgeIDType,MultiExpandInputType> EdgeIntersectionContainerType;
			
			EdgeIntersectionContainerType edgeIntersectionContainer;
			for (IntersectionIDType interID=0;interID!=intersectionContainer.size();++interID){
				
				const EdgeIDType key1(intersectionContainer[interID]. first.first->nodeIDPair);
				const EdgeIDType key2(intersectionContainer[interID].second.first->nodeIDPair);
				
                
                const isNetworkLinkType L1(DN.link(key1.first,key1.second));
                const isNetworkLinkType L2(DN.link(key2.first,key2.second));
                
//                const double cR1(dx/L1.second->chordLength());
//                const double cR2(dx/L2.second->chordLength());
//                const double cRmax(0.1);
//                
//                const double du1((cR1<cRmax)? cR1 : cRmax);
//				const double du2((cR2<cRmax)? cR2 : cRmax);
                
				const double du1(dx/L1.second->chordLength());
				const double du2(dx/L2.second->chordLength());
				
				if (intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection &&
					intersectionContainer[interID].second.second-du2 > avoidNodeIntersection && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection &&
					dirVector[interID]!=0){
					// Limit to 1 intersection per segment
					if(edgeIntersectionContainer.find(key1)==edgeIntersectionContainer.end() && edgeIntersectionContainer.find(key2)==edgeIntersectionContainer.end()){
						std::cout<<"key1 is "<<key1.first<<" "<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
						std::cout<<"key2 is "<<key2.first<<" "<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
						
                        
                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        const double u2m(intersectionContainer[interID].second.second-dirVector[interID]*du2);

                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.boundary_type){
                            SearchData<dim> SD1(L1.second->get_r(u1m));
                            DN.shared.domain.findIncludingTet(SD1,L1.second->source->meshID());
                            firstIntersectionInsideMesh*=(SD1.nodeMeshLocation==insideMesh);
                            if(firstIntersectionInsideMesh){
                                SearchData<dim> SD2(L2.second->get_r(u2m));
                                DN.shared.domain.findIncludingTet(SD2,L2.second->source->meshID());
                                firstIntersectionInsideMesh*=(SD2.nodeMeshLocation==insideMesh);
                            }
                        }
                        if(firstIntersectionInsideMesh){
                            assert(edgeIntersectionContainer[key1].insert(std::make_pair(u1m,2*interID)).second);
                            assert(edgeIntersectionContainer[key2].insert(std::make_pair(u2m,2*interID)).second);
                        }
                        
//						assert(edgeIntersectionContainer[key1].insert(std::make_pair(intersectionContainer[interID]. first.second-du1,2*interID)).second);
//						assert(edgeIntersectionContainer[key2].insert(std::make_pair(intersectionContainer[interID].second.second-dirVector[interID]*du2,2*interID)).second);


						
                        const double u1p(intersectionContainer[interID]. first.second+du1);
                        const double u2p(intersectionContainer[interID].second.second+dirVector[interID]*du2);
                        
                        bool secondIntersectionInsideMesh(true);
                        if (DN.shared.boundary_type){
                            SearchData<dim> SD1(L1.second->get_r(u1p));
                            DN.shared.domain.findIncludingTet(SD1,L1.second->source->meshID());
                            secondIntersectionInsideMesh*=(SD1.nodeMeshLocation==insideMesh);
                            if(secondIntersectionInsideMesh){
                                SearchData<dim> SD2(L2.second->get_r(u2p));
                                DN.shared.domain.findIncludingTet(SD2,L2.second->source->meshID());
                                secondIntersectionInsideMesh*=(SD2.nodeMeshLocation==insideMesh);
                            }
                        }
                        if(secondIntersectionInsideMesh){
                            assert(edgeIntersectionContainer[key1].insert(std::make_pair(u1p,2*interID+1)).second);
                            assert(edgeIntersectionContainer[key2].insert(std::make_pair(u2p,2*interID+1)).second);	
                        }
                        
                        
//						assert(edgeIntersectionContainer[key1].insert(std::make_pair(intersectionContainer[interID]. first.second+du1,2*interID+1)).second);
//						assert(edgeIntersectionContainer[key2].insert(std::make_pair(intersectionContainer[interID].second.second+dirVector[interID]*du2,2*interID+1)).second);						
					

					}
				}
			}
			
			std::cout<<"Using "<< edgeIntersectionContainer.size()/2<<" intersections. ";
			
			// Call Multiexpand
			std::map<IntersectionIDType, std::set<size_t> > nodeIntersectionMap;
			for (typename EdgeIntersectionContainerType::const_iterator edgeIter=edgeIntersectionContainer.begin();edgeIter!=edgeIntersectionContainer.end();++edgeIter){
				const size_t i(edgeIter->first. first);
				const size_t j(edgeIter->first.second);
				std::map<IntersectionIDType,size_t> multiExpandOut=DN.multiExpand(i,j,edgeIter->second);
				for (std::map<IntersectionIDType,size_t>::const_iterator iter=multiExpandOut.begin();iter!=multiExpandOut.end();++iter){
					nodeIntersectionMap[iter->first].insert(iter->second);
				}
			}
			
            // Call Contract
			for (std::map<IntersectionIDType, std::set<size_t> >::const_iterator mapIter=nodeIntersectionMap.begin();mapIter!=nodeIntersectionMap.end();++mapIter){
				assert(mapIter->second.size()==2 && "THERE SHOULD BE TWO NODES IN EACH INTERSECTION");
				const size_t i(*(mapIter->second.begin()));
				const size_t j(*(mapIter->second.rbegin()));
				
				
				const isNetworkNodeType Ni(DN.node(i));
				assert(Ni.first);
				assert(Ni.second->is_simple());
				const isNetworkNodeType Nj(DN.node(j));
				assert(Nj.first);
				assert(Nj.second->is_simple());
				
				const VectorDimD P1(Ni.second->get_P());
				const VectorDimD P2(Nj.second->get_P());
				const VectorDimD N1(Ni.second->openNeighborLink(0)->glidePlaneNormal); // CHANGE THIS
				const VectorDimD N2(Nj.second->openNeighborLink(0)->glidePlaneNormal); // CHANGE THIS
                const bool isIsessile(Ni.second->openNeighborLink(0)->sessilePlaneNormal.norm()>FLT_EPSILON);
                const bool isJsessile(Nj.second->openNeighborLink(0)->sessilePlaneNormal.norm()>FLT_EPSILON);
                
                
                
                if (!isIsessile && !isJsessile){
                    VectorDimD linePoint((P1+P2)*0.5);
                    const double denom(1.0-std::pow(N1.dot(N2),2));
                    if(std::fabs(denom)>FLT_EPSILON){ // planes are incident: make sure to place intersection on common line 
                        const double numer((P2-P1).dot(N2));
                        const double u(numer/denom);
                        linePoint = P1+(N2-N2.dot(N1)*N1)*u;
                    }
                    
                    
                    bool linePointInsideMesh(true);
                    if (DN.shared.boundary_type){
                        SearchData<dim> SD(linePoint);
                        DN.shared.domain.findIncludingTet(SD,Ni.second->meshID());
                        linePointInsideMesh*=(SD.nodeMeshLocation==insideMesh);
                    }
                    if(linePointInsideMesh){
                        std::cout<<"Glissile Junction"<<std::endl;
                        DN.contract(i,j,linePoint);
                    }

                }
                else if(isIsessile && !isJsessile){ // use P1 (which is on sessile segment) as the intersection point
                    //   const VectorDimD dir(Ni.second->openNeighborLink(0)->glidePlaneNormal.cross(Ni.second->openNeighborLink(0)->sessilePlaneNormal).normalized());
                    //   const double denom(dir.dot(N2));
                    //   if (std::fabs(denom)>FLT_EPSILON){
                    //       const double u((P2-P1).dot(N2)/denom);
                    // if (std::fabs(denom)<dx*0.25){
                    std::cout<<"First-Sessile Junction"<<std::endl;
                    //          DN.contract(i,j,P1+u*dir);
                    DN.contract(i,j,P1);
                    
                    //}
                    //}
                }
                else if(!isIsessile && isJsessile){ // use P2 (which is on sessile segment) as the intersection point
                    //  const VectorDimD dir(Nj.second->openNeighborLink(0)->glidePlaneNormal.cross(Nj.second->openNeighborLink(0)->sessilePlaneNormal).normalized());
                    //  const double denom(dir.dot(N1));
                    //  if (std::fabs(denom)>FLT_EPSILON){
                    //     const double u((P1-P2).dot(N1)/denom);
                    // if (std::fabs(denom)<dx*0.25){
                    std::cout<<"Second-Sessile Junction"<<std::endl;
                    //         DN.contract(i,j,P2+u*dir);
                    DN.contract(i,j,P2);
                    
                    // }
                    // }
                }
                else{
                    assert(0 && "CANNOT MAKE JUNCTION BETWEEN TWO SESSILE SEGMENTS.");
                    
                }
                
                
				
				
			}
            
            
            
            std::vector<std::pair<size_t,size_t> > nodeContractVector;
            
            for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
            {
                bool isSimple(nodeIter->second->is_simple());
                if(isSimple) // has two neighbors
                {
                    const NodeType* const pNi(nodeIter->second->openNeighborNode(0));
                    const NodeType* const pNj(nodeIter->second->openNeighborNode(1));
                    if(pNi->openOrder()>2 && pNj->openOrder()>2)
                    {
                        const isNetworkLinkType lIJ(DN.link(pNi->sID,pNj->sID));
                        const isNetworkLinkType lJI(DN.link(pNj->sID,pNi->sID));
                        if(lIJ.first || lJI.first)
                        {
                            nodeContractVector.push_back(std::make_pair(pNi->sID,nodeIter->second->sID));
                            
                        }
                    }
                }
                
            }
            
            
            for (unsigned int k=0;k<nodeContractVector.size();++k)
            {
                DN.contractSecond(nodeContractVector[k].first,nodeContractVector[k].second);
            }
            
            
            
		}
        
        

        
        
	};
    
    
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;

	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif






//        /* dontIntersect ******************************/
//        bool dontIntersect(const LinkType& L1, const LinkType& L2) const {
//            
//            bool areIncidentAtNodes(L1.source->sID==L2.source->sID || L1.source->sID==L2.sink->sID || L1.sink->sID==L2.source->sID || L1.sink->sID==L2.sink->sID); 
//            bool areOnDifferentPlanes((L1.glidePlaneNormal-L2.glidePlaneNormal).squaredNorm()>FLT_EPSILON);
//            bool L1isSessile(L1.sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
//            bool L2isSessile(L2.sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
//            bool eitherOneIsSessile(L1isSessile || L2isSessile);            
//            bool areIncidentAndOnDifferentPlanes(areIncidentAtNodes && areOnDifferentPlanes && eitherOneIsSessile);
//            
//            bool bothAreSessile(L1isSessile && L2isSessile);
//            
////            
////            
////            if (L1isSessile && !L2isSessile){
////                std::cout<<"L1 is "<<L1.source->sID<<"->"<<L1.sink->sID<<std::endl;
////                std::cout<<"L2 is "<<L2.source->sID<<"->"<<L2.sink->sID<<std::endl;
////                std::cout<<"L1.glidePlaneNormal="<<L1.glidePlaneNormal.transpose()<<std::endl;
////                std::cout<<"L1.sessilePlaneNormal="<<L1.sessilePlaneNormal.transpose()<<std::endl;
////                std::cout<<"L2.glidePlaneNormal="<<L2.glidePlaneNormal.transpose()<<std::endl;
////            }
////            
////            if (!L1isSessile && L2isSessile){
////                std::cout<<"L1 is "<<L1.source->sID<<"->"<<L1.sink->sID<<std::endl;
////                std::cout<<"L2 is "<<L2.source->sID<<"->"<<L2.sink->sID<<std::endl;
////                std::cout<<"L1.glidePlaneNormal="<<L1.glidePlaneNormal.transpose()<<std::endl;
////                std::cout<<"L2.glidePlaneNormal="<<L2.glidePlaneNormal.transpose()<<std::endl;
////                std::cout<<"L2.sessilePlaneNormal="<<L2.sessilePlaneNormal.transpose()<<std::endl;
////            }
////            
//            
//            bool L2NotOnL1Planes( L1isSessile && !L2isSessile && (L2.glidePlaneNormal-L1.glidePlaneNormal).norm()>FLT_EPSILON && (L2.glidePlaneNormal-L1.sessilePlaneNormal).norm()>FLT_EPSILON);
//            bool L1NotOnL2Planes(!L1isSessile &&  L2isSessile && (L1.glidePlaneNormal-L2.glidePlaneNormal).norm()>FLT_EPSILON && (L1.glidePlaneNormal-L2.sessilePlaneNormal).norm()>FLT_EPSILON);
//            
//            return areIncidentAndOnDifferentPlanes || bothAreSessile || L2NotOnL1Planes || L1NotOnL2Planes;
//            //            return areIncidentAndOnDifferentPlanes || bothAreSessile;
//            
//            //           return areIncident && areOnDifferentPlanes;
//        }
