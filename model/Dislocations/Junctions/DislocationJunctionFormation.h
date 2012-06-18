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

//#include <touple.h>


namespace model {
	
	template <typename DislocationNetworkType>
	class DislocationJunctionFormation{
		
		typedef typename DislocationNetworkType::LinkType LinkType;
		typedef typename DislocationNetworkType::isNetworkNodeType isNetworkNodeType;
		
		enum {dim=3};
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		
		typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
		
		typedef std::pair<const LinkType*, double> EdgeIntersectionType;
		typedef std::pair<EdgeIntersectionType,EdgeIntersectionType> EdgeIntersectionPairType;
		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
		
		
		
		// a reference to the dislocation network
		DislocationNetworkType& DN;
		
		//	std::touple<int,double> tp;
		
	public:
		
		
		/* Constructor ******************************/
		DislocationJunctionFormation(DislocationNetworkType& DN_in) : DN(DN_in) {}
		
		/* findIntersections **************************************************/
		EdgeIntersectionPairContainerType findIntersections(const double& avoidNodeIntersection, const double& tol=FLT_EPSILON) const { 
			
			EdgeIntersectionPairContainerType intersectionContainer;
			
			
			//	const double avoidNodeIntersection=0.01;
			
			//! 2- loop over all links and determine their intersections
			for (typename NetworkLinkContainerType::const_iterator linkIterA=DN.linkBegin();linkIterA!=DN.linkEnd();linkIterA++){
				for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.linkEnd();linkIterB++){
					if (linkIterA->second->sID!=linkIterB->second->sID){							
						//						std::cout<<"Intersecting link "<<linkIterA->second->nodeIDPair.first<<"->"<<linkIterA->second->nodeIDPair.second<<" and"<<
						//						linkIterB->second->nodeIDPair.first<<"->"<<linkIterB->second->nodeIDPair.second<<std::endl;
						std::set<std::pair<double,double> > temp ( linkIterA->second->intersectWith(linkIterB->second,tol));
						for (std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paramIter!=temp.end();++paramIter){					
							
							//					const double avoidNodeIntersection1(0.05);
							//					const double avoidNodeIntersection2(0.05);
							//					const double avoidNodeIntersection1(distanceFromNoodes/linkIterA->second->chordLength());
							//					const double avoidNodeIntersection2(distanceFromNoodes/linkIterB->second->chordLength());
							
							if (   paramIter->first >avoidNodeIntersection && paramIter-> first<1.0-avoidNodeIntersection
								&& paramIter->second>avoidNodeIntersection && paramIter->second<1.0-avoidNodeIntersection){
								
								EdgeIntersectionType intersectionOnA(std::make_pair(&(*linkIterA->second),paramIter->first ));
								EdgeIntersectionType intersectionOnB(std::make_pair(&(*linkIterB->second),paramIter->second));
								EdgeIntersectionPairType intersectionPair(std::make_pair(intersectionOnA,intersectionOnB));
								intersectionContainer.push_back(intersectionPair);
								
							}
//							else{
//								std::cout<<"avoidNodeIntersection between" <<linkIterA->second->nodeIDPair.first<<"->"<<linkIterA->second->nodeIDPair.second<<" and"<<
//														linkIterB->second->nodeIDPair.first<<"->"<<linkIterB->second->nodeIDPair.second<<std::endl;
//
//							}
						}
					}
				}
			}	
			return intersectionContainer;
		}
		
		
		
		
		
		
		
		/* formJunctions **************************************************/
		void formJunctions(const double& dl, const double& avoidNodeIntersection) {
			
			
			//! 1- Initialize intersectionContainer calling findIntersections
			EdgeIntersectionPairContainerType intersectionContainer(findIntersections(avoidNodeIntersection));
			std::cout<<intersectionContainer.size()<<" geometric) ";
			
			//! 2- Remove from intersectionContainer all intersections that don't satisfy Frank's rule
			std::vector<int> dirVector;			
			for (typename EdgeIntersectionPairContainerType::iterator iter=intersectionContainer.begin();iter!=intersectionContainer.end();){
				double u1(iter-> first.second);
				double u2(iter->second.second);
				VectorDimD b1(iter-> first.first->flow);
				VectorDimD b2(iter->second.first->flow);
				VectorDimD rl1(iter-> first.first->get_rl(u1));
				VectorDimD rl2(iter->second.first->get_rl(u2));
				const bool frankRule(b1.dot(b2)*rl1.dot(rl2)<=0.0);
//				const bool frankRule(true);
				if (!frankRule){
					iter=intersectionContainer.erase(iter);
				}
				else{
					int dir((rl1.dot(rl2) > 0.0) ? 1 : ((rl1.dot(rl2) < 0.0) ? -1 : 0)); // RETURNING 0 HERE MAY BE A PROBLEM
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
				
				EdgeIDType key1=intersectionContainer[interID]. first.first->nodeIDPair;
				EdgeIDType key2=intersectionContainer[interID].second.first->nodeIDPair;
				
				double du1(dl/DN.link(key1.first,key1.second).second->chordLength());
				double du2(dl/DN.link(key2.first,key2.second).second->chordLength());
				
				
				if (intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection &&
					intersectionContainer[interID].second.second-du2 > avoidNodeIntersection && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection){
					// Limit to 1 intersection per segment
					if(edgeIntersectionContainer.find(key1)==edgeIntersectionContainer.end() && edgeIntersectionContainer.find(key2)==edgeIntersectionContainer.end()){
//						std::cout<<"key1 is "<<key1.first<<" "<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
//						std::cout<<"key2 is "<<key2.first<<" "<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
						
						assert(edgeIntersectionContainer[key1].insert(std::make_pair(intersectionContainer[interID]. first.second-du1,2*interID)).second);
						assert(edgeIntersectionContainer[key2].insert(std::make_pair(intersectionContainer[interID].second.second-dirVector[interID]*du2,2*interID)).second);
						
						assert(edgeIntersectionContainer[key1].insert(std::make_pair(intersectionContainer[interID]. first.second+du1,2*interID+1)).second);
						assert(edgeIntersectionContainer[key2].insert(std::make_pair(intersectionContainer[interID].second.second+dirVector[interID]*du2,2*interID+1)).second);						
					}
				}
			}
			
			std::cout<<"Using "<< edgeIntersectionContainer.size()/2<<" intersections. ";
			
			
			std::map<IntersectionIDType, std::set<size_t> > nodeIntersectionMap;
			for (typename EdgeIntersectionContainerType::const_iterator edgeIter=edgeIntersectionContainer.begin();edgeIter!=edgeIntersectionContainer.end();++edgeIter){
				const size_t i(edgeIter->first. first);
				const size_t j(edgeIter->first.second);
				std::map<IntersectionIDType,size_t> multiExpandOut=DN.multiExpand(i,j,edgeIter->second);
				for (std::map<IntersectionIDType,size_t>::const_iterator iter=multiExpandOut.begin();iter!=multiExpandOut.end();++iter){
					nodeIntersectionMap[iter->first].insert(iter->second);
				}
			}
			
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

				VectorDimD linePoint((P1+P2)*0.5);
				
				const double denom= 1.0-std::pow(N1.dot(N2),2);
				if(std::fabs(denom)>FLT_EPSILON){ // planes are incident 
					const double numer= (P2-P1).dot(N2);
					double u=numer/denom;
					linePoint = P1+(N2-N2.dot(N1)*N1)*u;
				}
				
				DN.contract(i,j,linePoint);
				
			}			
		}
		
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif
