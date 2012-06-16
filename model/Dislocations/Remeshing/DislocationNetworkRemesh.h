/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORKREMESH_H_
#define model_DISLOCATIONNETWORKREMESH_H_

#include <assert.h>

#include<Eigen/Dense>
#include <model/Network/Operations/EdgeFinder.h>

namespace model {
	
	template <typename DislocationNetworkType>
	class DislocationNetworkRemesh{
		
		typedef typename DislocationNetworkType::LinkType LinkType;
		
		enum{dim=3}; // CHANGE THIS
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		
		typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
		
		DislocationNetworkType& DN;
		
	public:
		
		/* Constructor ******************************/
		DislocationNetworkRemesh(DislocationNetworkType& DN_in) : DN(DN_in) {}
		
		/************************************************************/
		// remesh
		void remesh(const double& Lmax, const double& Lmin,const double& thetaDeg){
			remeshByContraction(Lmin);
			remeshByExpansion(Lmax,Lmin,thetaDeg); // CALL THIS AFTER CONTRACTION TO EXPAND 2-NODES SUBNETWORKS
		}
		
		/* remeshByContraction **************************************/
		void remeshByContraction(const double& Lmin){
			
			const double vTolcont=0.0;
			
			std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
			
			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){				
				
				VectorDimD chord(linkIter->second->chord()); // this is sink->get_P() - source->get_P()
				double chordLength(chord.norm());
				VectorDimD dv=linkIter->second->sink->get_V()-linkIter->second->source->get_V();				
				bool endsAreApproaching((chord.dot(dv))<-vTolcont*chordLength*dv.norm());				
				
				if (endsAreApproaching && chordLength<Lmin){// toBeContracted part
					toBeContracted.insert(std::make_pair(chordLength,linkIter->second->nodeIDPair));
				}
			}
			
			
			// Call Network::contract 
			for (std::set<std::pair<double,std::pair<size_t,size_t> > >::const_iterator smallIter=toBeContracted.begin(); smallIter!=toBeContracted.end(); ++smallIter) {
				const size_t i(smallIter->second.first);
				const size_t j(smallIter->second.second);
				typename EdgeFinder<LinkType>::isNetworkEdgeType Lij=DN.link(i,j);
				
				if (Lij.first ){
					bool   SinkRemovable(Lij.second->sink->is_removable());
					bool SourceRemovable(Lij.second->source->is_removable());
					
					
					if (SourceRemovable && SinkRemovable) {
						DN.contract(i,j,Lij.second->get_r(0.5));
					}
					else if (!SourceRemovable && SinkRemovable) {
						DN.contractSecond(i,j);						
					}
					else if (SourceRemovable && !SinkRemovable) {
						DN.contractSecond(j,i);
					}
					else { 
						
						
						
						// Store source and sink positions
						const VectorDimD P1(Lij.second->source->get_P());
						const VectorDimD P2(Lij.second-> sink->get_P());
						const VectorDimD C(P2-P1);
						const double cNorm(C.norm());
						if (cNorm<FLT_EPSILON){ // nodes are on to of each other
							DN.contract(i,j,0.5*(P1+P2)); 
						}
						else{ // only case when this can happen is when source moves on a line and sink moves on a line and the two lines intersect at one point
							
							const typename DislocationNetworkType::NodeType::VectorOfNormalsType CNsource=Lij.second->source->constraintNormals();
							const typename DislocationNetworkType::NodeType::VectorOfNormalsType CNsink  =Lij.second->  sink->constraintNormals();
							
							if (CNsource.size()==2 && CNsink.size()==2){
								// check if the lines X=P1+d1*u1 and X=P2+d2*u2 intersect at one point
								
								// Compute first direction
								VectorDimD d1(CNsource[0].cross(CNsource[1]));
								double d1norm(d1.norm());
								assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
								d1/=d1norm;
								
								// Compute second direction
								VectorDimD d2(CNsink[0].cross(CNsink[1]));
								double d2norm(d2.norm());
								assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
								d2/=d2norm;
								
								// Necessary condition is plannarity: C*(d1xd2)=0
								const VectorDimD d3(d1.cross(d2));
								const double d3Norm2(d3.squaredNorm());
								if (d3Norm2<FLT_EPSILON){ // colinear or parallel
									if (std::fabs((C/cNorm).dot(d1))<FLT_EPSILON){ // colinear
										DN.contract(i,j,0.5*(P1+P2)); 
									}
								}
								else{ // coplanar or no-intersection
									bool isPlanarConfiguration(std::fabs((C/cNorm).dot(d3))<FLT_EPSILON);
									if (isPlanarConfiguration){ // coplanar
										//std::cout<<"NEW CONSTRACT NONREMOVABLE"<<std::endl;
										const double u1=C.cross(d2).dot(d3)/d3Norm2;
										DN.contract(i,j,P1+d1*u1); 
									}
								}
								

								
							}
							
							
							
						}
						
						
						
						
						
						
						
						//						//const std::vector<VectorDimD> CNsource=Lij.second->source->constraintNormals();
						//						//const std::vector<VectorDimD> CNsink  =Lij.second->  sink->constraintNormals();
						//						
						//						const VectorDimD chord(Lij.second->chord());
						//						
						//						bool canContract(true);
						//						for (unsigned int kk=0;kk<CNsource.size();++kk){
						//							canContract*=(std::fabs(CNsource[kk].dot(chord))<FLT_EPSILON);
						//						}
						//						for (unsigned int kk=0;kk<CNsink.size();++kk){
						//							canContract*=(std::fabs(CNsink[kk].dot(chord))<FLT_EPSILON);
						//						}
						//						if(canContract){
						//							DN.contract(i,j,Lij.second->get_r(0.5)); 
						//							
						//						}
					}
					
				}
			}
			
			
		}
		
		
		
		/* remeshByExpansion ****************************************************************/
		void remeshByExpansion(const double& Lmax, const double& Lmin,const double& thetaDeg){
			
			
			double cos_theta_max_crit = std::cos(M_PI-thetaDeg*M_PI/180.0);  /*critical angle */
			std::set<std::pair<size_t,size_t> > toBeExpanded;
			
			
			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){
				
				VectorDimD chord(linkIter->second->chord()); // this is sink->get_P() - source->get_P()
				double chordLength(chord.norm());
				VectorDimD dv=linkIter->second->sink->get_V()-linkIter->second->source->get_V();				
				
				
				// Always expand single FR source segment 
				//				if (!linkIter->second->source->is_balanced() && !linkIter->second->sink->is_balanced()){
				//					assert(toBeExpanded.insert(linkIter->second->nodeIDPair).second && "COULD NOT INSERT LINK.");
				//				}
				
				if (linkIter->second->source->constraintNormals().size()>2 && linkIter->second->sink->constraintNormals().size()>2){
					assert(toBeExpanded.insert(linkIter->second->nodeIDPair).second && "COULD NOT INSERT LINK.");
				}
				
				
				
				if (chordLength>Lmax){
					toBeExpanded.insert(linkIter->second->nodeIDPair);
				}
				
				// Angle criterion
				double vTolexp=0.5;
				
				if (linkIter->second->source->is_simple()){ //check angle criterion at source
					VectorDimD c0(linkIter->second->source->openNeighborNode(0)->get_P()-linkIter->second->source->get_P());
					VectorDimD c1(linkIter->second->source->openNeighborNode(1)->get_P()-linkIter->second->source->get_P());
					double c0norm(c0.norm());
					double c1norm(c1.norm());
					if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm){
						if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/){
							toBeExpanded.insert(linkIter->second->source->openNeighborLink(0)->nodeIDPair);
						} 
						
						if (c1norm>3.0*Lmin /*&& c1.dot(v1)>vTolexp*c1norm*v1.norm()*/){
							toBeExpanded.insert(linkIter->second->source->openNeighborLink(1)->nodeIDPair);
						}
					}
				}
				if (linkIter->second->sink->is_simple()){ //check angle criterion at sink
					VectorDimD c0(linkIter->second->sink->openNeighborNode(0)->get_P()-linkIter->second->sink->get_P());
					VectorDimD c1(linkIter->second->sink->openNeighborNode(1)->get_P()-linkIter->second->sink->get_P());
					double c0norm(c0.norm());
					double c1norm(c1.norm());
					if(c0.dot(c1)>cos_theta_max_crit*c0norm*c1norm){
						if (c0norm>3.0*Lmin /*&& c0.dot(v0)>vTolexp*c0norm*v0.norm()*/){
							toBeExpanded.insert(linkIter->second->sink->openNeighborLink(0)->nodeIDPair);
						} 
						if (c1norm>3.0*Lmin/* && c1.dot(v1)>vTolexp*c1norm*v1.norm()*/){
							//														std::cout<<"Expanding 4"<<std::endl;
							toBeExpanded.insert(linkIter->second->sink->openNeighborLink(1)->nodeIDPair);
						}							
					}
				}
				
				if (!linkIter->second->source->is_simple() && !linkIter->second->sink->is_simple() 
					/*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/ && chordLength>3.0*Lmin){ // also expands a straight line to generate glissile segment 
					toBeExpanded.insert(linkIter->second->nodeIDPair);
				}
			}
			
			
			// Call Network::expand
			double expand_at(0.5);
			for (std::set<std::pair<size_t,size_t> >::const_iterator expIter=toBeExpanded.begin(); expIter!=toBeExpanded.end(); ++expIter) {
				const size_t i(expIter->first);
				const size_t j(expIter->second);				
				if(DN.link(i,j).first){
					DN.expand(i,j,expand_at);				
				}
			}
			
		}
		
		
		
		
		
		/* remeshByContraction **************************************/
		void contract0chordSegments(){
			std::set<std::pair<double,std::pair<size_t,size_t> > > toBeContracted; // order by increasing segment length
			
			for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){								
				VectorDimD chord(linkIter->second->chord()); // this is sink->get_P() - source->get_P()
				double chordLength(chord.norm());
				if (chordLength<=FLT_EPSILON){// toBeContracted part
					toBeContracted.insert(std::make_pair(chordLength,linkIter->second->nodeIDPair));
				}
			}
			
			// Call Network::contract 
			for (std::set<std::pair<double,std::pair<size_t,size_t> > >::const_iterator smallIter=toBeContracted.begin(); smallIter!=toBeContracted.end(); ++smallIter) {
				const size_t i(smallIter->second.first);
				const size_t j(smallIter->second.second);
				typename EdgeFinder<LinkType>::isNetworkEdgeType Lij=DN.link(i,j);				
				if (Lij.first ){
					DN.contractSecond(i,j);						
				}
			}
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif




//					if ( linkIter->second->source->is_removable() && linkIter->second->sink->is_removable() ){ // P0-P1->P2-P3
//						VectorDimD P0(linkIter->second->source->openNeighborNode(0)->get_P());
//						if(linkIter->second->source->openNeighborNode(0)->sID==linkIter->second->sink->sID){
//							P0=linkIter->second->source->openNeighborNode(1)->get_P();
//						}
//						VectorDimD P2(linkIter->second->sink->openNeighborNode(0)->get_P());
//						if(linkIter->second->sink->openNeighborNode(0)->sID==linkIter->second->source->sID){
//							P2=linkIter->second->sink->openNeighborNode(1)->get_P();
//						}
//						VectorDimD c02(P2-P0);
//						//VectorDimD c(linkIter->second->chord());
//						//double cNorm(c.norm());
//						if(std::fabs(chord.dot(c02))>cos_theta_min_crit*chord.norm()*chordLength){
//							toBeContracted.insert(std::make_pair(chordLength,linkIter->second->nodeIDPair));
//						}
//					} //both source and sink are removable




//					if (linkIter->second->source->is_simple()){ //check angle criterion at source
//						VectorDimD c0(linkIter->second->source->openNeighborNode(0)->get_P()-linkIter->second->source->get_P());
//						VectorDimD c1(linkIter->second->source->openNeighborNode(1)->get_P()-linkIter->second->source->get_P());
//						double c0norm(c0.norm());
//						double c1norm(c1.norm());
//						if(c0.dot(c1)<cos_theta_min_crit*c0norm*c1norm){
//							if (c0norm<=c1norm){
//	//							std::cout<<"Im here1"<<std::endl;
//								toBeContracted.insert(std::make_pair(c0norm,linkIter->second->source->openNeighborLink(0)->nodeIDPair));
//	//							std::cout<<"Im here2"<<std::endl;
//							}
//							else {
//	//							std::cout<<"Im here3"<<std::endl;
//								toBeContracted.insert(std::make_pair(c1norm,linkIter->second->source->openNeighborLink(1)->nodeIDPair));
//	//							std::cout<<"Im here4"<<std::endl;
//							}
//						}
//					}
//					if (linkIter->second->sink->is_simple()){ //check angle criterion at sink
//						VectorDimD c0(linkIter->second->sink->openNeighborNode(0)->get_P()-linkIter->second->sink->get_P());
//						VectorDimD c1(linkIter->second->sink->openNeighborNode(1)->get_P()-linkIter->second->sink->get_P());
//						double c0norm(c0.norm());
//						double c1norm(c1.norm());
//						if(c0.dot(c1)<cos_theta_min_crit*c0norm*c1norm){
//							if (c0norm<=c1norm){
//	//							std::cout<<"Im here5"<<std::endl;
//								toBeContracted.insert(std::make_pair(c0norm,linkIter->second->sink->openNeighborLink(0)->nodeIDPair));
//	//							std::cout<<"Im here6"<<std::endl;
//							}
//							else {
//	//							std::cout<<"Im here7"<<std::endl;
//								toBeContracted.insert(std::make_pair(c1norm,linkIter->second->sink->openNeighborLink(1)->nodeIDPair));
//	//							std::cout<<"Im here8"<<std::endl;
//							}
//						}
//					}
