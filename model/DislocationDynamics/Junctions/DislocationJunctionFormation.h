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
#include <model/Geometry/PlanePlaneIntersection.h>

#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>
#include <model/DislocationDynamics/Remeshing/DislocationNetworkRemesh.h>
#include <model/BVP/SearchData.h>
#include <model/MPI/MPIcout.h>

//#include <model/Math/MatrixCompanion.h>


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
        //		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
		typedef std::deque<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
		
		
		
		//! A reference to the DislocationNetwork
		DislocationNetworkType& DN;
        
	public:
		
        //! The tolerance (in units of distance) used for collision detection
        static double collisionTol;
        static bool useVertexEdgeJunctions;
        
		
		/* Constructor ********************************************************/
		DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {}
		
		/* findIntersections **************************************************/
		EdgeIntersectionPairContainerType findEdgeEdgeIntersections(const double& avoidNodeIntersection) const
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
                        
			EdgeIntersectionPairContainerType intersectionContainer;
			
			
			//! 2- loop over all links and determine their intersections
			for (typename NetworkLinkContainerType::const_iterator linkIterA=DN.linkBegin();linkIterA!=DN.linkEnd();linkIterA++)
            {
                
                
                
                const DislocationSegmentIntersection<dim,pOrder> dsi(linkIterA->second-> hermiteCoefficients(),linkIterA->second->glidePlaneNormal);
                
				for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.linkEnd();linkIterB++)
                {
					if (linkIterA->second->sID!=linkIterB->second->sID) // don't intersect with itself
                    {
                        
                        //                model::cout<<"intersecting "<<linkIterA->second->source->sID<<"->"<<linkIterA->second->sink->sID<<" and "<<linkIterB->second->source->sID<<"->"<<linkIterB->second->sink->sID<<std::endl;
                        
                        //                        const bool areIncidentAtNodes(   linkIterA->second->source->sID==linkIterB->second->source->sID
                        //                                                      || linkIterA->second->source->sID==linkIterB->second->sink->sID
                        //                                                      || linkIterA->second->  sink->sID==linkIterB->second->source->sID
                        //                                                      || linkIterA->second->  sink->sID==linkIterB->second->sink->sID);
                        //
                        //                       const bool areOnDifferentPlanes((linkIterA->second->glidePlaneNormal-linkIterB->second->glidePlaneNormal).squaredNorm()>FLT_EPSILON);
                        const bool L1isSessile(linkIterA->second->sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
                        const bool L2isSessile(linkIterB->second->sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
                        
                        std::set<std::pair<double,double> > temp;
                        
                        
                        
                        if (!L1isSessile && !L2isSessile) // both are glissile
                        {
                            temp = dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                        }
                        
                        else if (!L1isSessile && L2isSessile) // L1 is glissile and L2 is sessile
                        {
                            const bool gnAgnB((linkIterA->second->glidePlaneNormal-linkIterB->second->glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                            const bool gnAsnB((linkIterA->second->glidePlaneNormal-linkIterB->second->sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                            
                            if (!gnAgnB && !gnAsnB){
                                // cannot intersect
                            }
                            else if(gnAgnB && !gnAsnB){ // use planeNormal of A and planeNormal of B
                                temp = dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else if (!gnAgnB && gnAsnB){ // use planeNormal of A and sessileNormal of B
                                temp = dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->sessilePlaneNormal,collisionTol);
                            }
                            else{
                                assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                            }
                        }
                        else if (L1isSessile && !L2isSessile) // L1 is sessile and L2 is glissile
                        {
                            const bool gnBgnA((linkIterB->second->glidePlaneNormal-linkIterA->second->glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                            const bool gnBsnA((linkIterB->second->glidePlaneNormal-linkIterA->second->sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                            
                            if (!gnBgnA && !gnBsnA){
                                // cannot intersect
                            }
                            else if(gnBgnA && !gnBsnA){ // use planeNormal of A and planeNormal of B
                                temp = dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else if (!gnBgnA && gnBsnA){ // use sessileNormal of A and use planeNormal of B
                                const DislocationSegmentIntersection<dim,pOrder> dsi2(linkIterA->second-> hermiteCoefficients(),linkIterA->second->sessilePlaneNormal);
                                temp = dsi2.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else{
                                assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                            }
                        }
                        else
                        { // both are sessile
                            // cannot intersect
                        }
                        
                        
                        
                        // std::set<std::pair<double,double> > temp ( dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol));
                        
                        
                        for (std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paramIter!=temp.end();++paramIter)
                        {
                            if (   paramIter->first >avoidNodeIntersection && paramIter-> first<1.0-avoidNodeIntersection
                                && paramIter->second>avoidNodeIntersection && paramIter->second<1.0-avoidNodeIntersection) // avoid node intersection
                            {
                                EdgeIntersectionType intersectionOnA(std::make_pair(&(*linkIterA->second),paramIter->first ));
                                EdgeIntersectionType intersectionOnB(std::make_pair(&(*linkIterB->second),paramIter->second));
                                EdgeIntersectionPairType intersectionPair(std::make_pair(intersectionOnA,intersectionOnB));
                                intersectionContainer.push_back(intersectionPair);
                            }
                        } // end for
                    } // end don't self-intersect
				} // end loop over second segment
			} // end loop over first segment
			return intersectionContainer;
		}
		
        /* formEdgeEdgeJunctions **********************************************/
        size_t formEdgeEdgeJunctions(const double& dx, const double& avoidNodeIntersection)
        {
            //! 1- Initialize intersectionContainer calling findIntersections
			EdgeIntersectionPairContainerType intersectionContainer(findEdgeEdgeIntersections(avoidNodeIntersection));
			model::cout<<intersectionContainer.size()<<" geometric) ";
			
			//! 2- Remove from intersectionContainer all intersections that don't satisfy Frank's rule
			std::vector<int> dirVector;
			for (typename EdgeIntersectionPairContainerType::iterator iter=intersectionContainer.begin();iter!=intersectionContainer.end();)
            {
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
                    if(commonLine.norm()>FLT_EPSILON)
                    { // planes are not parallel, intersection will be on common line
                        prjDir=commonLine.normalized();
                    }
                    else
                    {
                        const double rl1Norm(rl1.norm());
                        assert(rl1Norm>FLT_EPSILON && "TANGENT HAS ZERO NORM");
                        prjDir=rl1/rl1Norm;
                    }
                }
                else if(isIsessile && !isJsessile)
                { // use chord of I
                    prjDir=iter-> first.first->chord().normalized();
                }
                else if(!isIsessile && isJsessile)
                { // use chord of J
                    prjDir=iter->second.first->chord().normalized();
                }
                else
                {
                    assert(0 && "CANNOT DETERMINE COMMON LINE BETWEEN TWO SESSILE SEGMENTS.");
                }
                
                
                const double sgnrl1rl2(rl1.dot(prjDir)*rl2.dot(prjDir));
                
                
				//const bool frankRule(b1.dot(b2)*rl1.dot(rl2)<=0.0);
				const bool frankRule(b1.dot(b2)*sgnrl1rl2<=0.0);
                
                
                //				const bool frankRule(true);
				if (!frankRule)
                {
					iter=intersectionContainer.erase(iter);
				}
				else
                { // determine the relative orientation of the links
                    //					const int dir((rl1.dot(rl2) > 0.0) ? 1 : ((rl1.dot(rl2) < 0.0) ? -1 : 0));
					const int dir(( sgnrl1rl2 > 0.0) ? 1 : ((sgnrl1rl2 < 0.0) ? -1 : 0));
					dirVector.push_back(dir);
					++iter; // increment iterator only if didn't erase
				}
			}
			model::cout<<intersectionContainer.size()<<" physical intersections. ";
			assert(intersectionContainer.size()==dirVector.size());
            
			//! 3- Organize intersections by segments to use multiexpand
			typedef size_t IntersectionIDType;
			typedef std::map<double, IntersectionIDType> MultiExpandInputType;
			typedef std::pair<size_t,size_t> EdgeIDType;
			typedef std::map<EdgeIDType,MultiExpandInputType> EdgeIntersectionContainerType;
            
            const double dx2=std::pow(dx,2);
			
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
				
                
                
				if (
                    intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection
                    && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection
                    && intersectionContainer[interID].second.second-du2 > avoidNodeIntersection
                    && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection
                    
                    
                    //                       (L1.second->source->get_P()-L1.second->get_r(intersectionContainer[interID].first.second-du1) ).squaredNorm() > dx2
                    //                    && (L1.second->sink  ->get_P()-L1.second->get_r(intersectionContainer[interID].first.second+du1) ).squaredNorm() > dx2
                    //                    && (L2.second->source->get_P()-L2.second->get_r(intersectionContainer[interID].second.second-du2) ).squaredNorm() > dx2
                    //                    && (L2.second->sink  ->get_P()-L2.second->get_r(intersectionContainer[interID].second.second+du2) ).squaredNorm() > dx2
                    
                    && dirVector[interID]!=0)
                {
					// Limit to 1 intersection per segment
					if(   edgeIntersectionContainer.find(key1)==edgeIntersectionContainer.end()
                       && edgeIntersectionContainer.find(key2)==edgeIntersectionContainer.end())
                    {
                        //                        						model::cout<<"key1 is "<<key1.first<<" "<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                        //                        						model::cout<<"key2 is "<<key2.first<<" "<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
						
                        
                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        const double u2m(intersectionContainer[interID].second.second-dirVector[interID]*du2);
                        
                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.boundary_type){
//                            SearchData<dim> SD1(L1.second->get_r(u1m));
//                            DN.shared.domain.findIncludingTet(SD1,L1.second->source->meshID());
//                            firstIntersectionInsideMesh*=(SD1.nodeMeshLocation==insideMesh);
                            firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1m),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            if(firstIntersectionInsideMesh){
//                                SearchData<dim> SD2(L2.second->get_r(u2m));
//                                DN.shared.domain.findIncludingTet(SD2,L2.second->source->meshID());
//                                firstIntersectionInsideMesh*=(SD2.nodeMeshLocation==insideMesh);
                                firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2m),L2.second->source->includingSimplex(),FLT_EPSILON).first;

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
//                            SearchData<dim> SD1(L1.second->get_r(u1p));
//                            DN.shared.domain.findIncludingTet(SD1,L1.second->source->meshID());
//                            secondIntersectionInsideMesh*=(SD1.nodeMeshLocation==insideMesh);

                            secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1p),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            
                            if(secondIntersectionInsideMesh){
//                                SearchData<dim> SD2(L2.second->get_r(u2p));
//                                DN.shared.domain.findIncludingTet(SD2,L2.second->source->meshID());
//                                secondIntersectionInsideMesh*=(SD2.nodeMeshLocation==insideMesh);
                                
                                                            secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2p),L2.second->source->includingSimplex(),FLT_EPSILON).first;
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
			
			model::cout<<"Using "<< edgeIntersectionContainer.size()/2<<" intersections. ";
			
			// Call Multiexpand on segments in edgeIntersectionContainer. Store correspoinding new nodes in nodeIntersectionMap
			std::map<IntersectionIDType, std::set<size_t> > nodeIntersectionMap;
			for (typename EdgeIntersectionContainerType::const_iterator edgeIter=edgeIntersectionContainer.begin();
                 /*                                                  */ edgeIter!=edgeIntersectionContainer.end();
                 /*                                                */ ++edgeIter)
            {
				const size_t i(edgeIter->first. first);
				const size_t j(edgeIter->first.second);
				std::map<IntersectionIDType,size_t> multiExpandOut=DN.multiExpand(i,j,edgeIter->second);
				for (std::map<IntersectionIDType,size_t>::const_iterator iter=multiExpandOut.begin();iter!=multiExpandOut.end();++iter)
                {
					nodeIntersectionMap[iter->first].insert(iter->second);
				}
			}
			
            // Call Contract on nodes in nodeIntersectionMap
			for (std::map<IntersectionIDType, std::set<size_t> >::const_iterator mapIter=nodeIntersectionMap.begin();
                 /*                                                           */ mapIter!=nodeIntersectionMap.end();
                 /*                                                         */ ++mapIter)
            {
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
//                        SearchData<dim> SD(linePoint);
//                        DN.shared.domain.findIncludingTet(SD,Ni.second->meshID());
//                        linePointInsideMesh*=(SD.nodeMeshLocation==insideMesh);
                        
                        linePointInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(linePoint,Ni.second->includingSimplex(),FLT_EPSILON).first;

                        
                    }
                    if(linePointInsideMesh){
                        model::cout<<"Glissile Junction"<<std::endl;
                        DN.contract(i,j,linePoint);
                    }
                    
                }
                else if(isIsessile && !isJsessile){ // use P1 (which is on sessile segment) as the intersection point
                    model::cout<<"First-Sessile Junction"<<std::endl;
                    DN.contract(i,j,P1);
                }
                else if(!isIsessile && isJsessile){ // use P2 (which is on sessile segment) as the intersection point
                    model::cout<<"Second-Sessile Junction"<<std::endl;
                    DN.contract(i,j,P2);
                }
                else{
                    assert(0 && "CANNOT MAKE JUNCTION BETWEEN TWO SESSILE SEGMENTS.");
                    
                }
                
                
				
				
			}
            
            
            return nodeIntersectionMap.size();
        }
        
        
        /* formVertexEdgeJunctions ********************************************/
        size_t formVertexEdgeJunctions(const double& dx)
        {
            
            
            // Direct Node contraction
            std::vector<std::pair<size_t,std::pair<size_t,size_t> > > nodeContractVector;
            //            typename EdgeFinder<LinkType>::isNetworkEdgeType
            
            
            std::map<size_t,std::pair<std::pair<size_t,size_t>,VectorDimD> > vertexEdgeIntersectionContainer;
            
            for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
            {
                if(nodeIter->second->invertedMotion())
                {
                    if(nodeIter->second->outOrder()==1 && nodeIter->second->inOrder()==1 && nodeIter->second->is_balanced() ) // has one inflow and one outflow, therefore it makes sense to talk about B of that node
                    {
                        VectorDimD B(std::get<1>(nodeIter->second->outNeighborhood().begin()->second)->flow);
                        VectorDimD T(std::get<1>(nodeIter->second->outNeighborhood().begin()->second)->sourceT());
                        
                        const typename DislocationNetworkType::NodeType::VectorOfNormalsType vertexPN(GramSchmidt<dim>(nodeIter->second->constraintNormals()));
                        
                        if(vertexPN.size()==1) //node is not constrained, and a unique plane normal si associated with that node
                        {
                            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
                            {
                                
                                
                                
                                if (nodeIter->first != linkIter->first.first && nodeIter->first != linkIter->first.second)
                                {
                                    
                                    const double chord2(linkIter->second->chord().squaredNorm());
                                    if (  (linkIter->second->source->get_P()-nodeIter->second->get_P()).squaredNorm()<chord2    // vertex is not too far from edge source
                                        ||(linkIter->second->  sink->get_P()-nodeIter->second->get_P()).squaredNorm()<chord2 )  // vertex is not too far from edge sink
                                    {
                                        std::map<double,std::pair<double,VectorDimD> > rootMap( linkIter->second->closestPoint(nodeIter->second->get_P() ) );
                                        
                                        if (rootMap.size())
                                        {
                                            
                                            const double d(rootMap.begin()->first);
                                            const double u(rootMap.begin()->second.first);
                                            
                                            if (d<dx)
                                            {
                                                
                                                const bool frankRule(B.dot(linkIter->second->flow)*T.dot(linkIter->second->get_ru(u))<=0.0);
                                                
                                                if (frankRule)
                                                {
                                                    
                                                    //                                            vertexEdgeIntersectionContainer.insert(std::make_pair(nodeIter->second->sID,std::make_pair(,)));
                                                    
                                                    
                                                    const VectorDimD N1(vertexPN[0]);
                                                    const VectorDimD N2(linkIter->second->glidePlaneNormal);
                                                    const VectorDimD P1(nodeIter->second->get_P());
                                                    const VectorDimD P2(linkIter->second->source->get_P());
                                                    
                                                    const int ppt(PlanePlaneIntersection<dim>::planePlaneType(N1,N2,P1,P2));
                                                    
                                                    
                                                    switch (ppt)
                                                    {
                                                        case 0: // paralle planes
                                                            // do nothing, no intersection possile
                                                            break;
                                                            
                                                        case 1: // coincident planes
                                                        {
                                                            
                                                            break;
                                                        }
                                                            
                                                        default: // incident planes
                                                        {
                                                            const double denom(1.0-std::pow(N1.dot(N2),2));
                                                            if(std::fabs(denom)>FLT_EPSILON)
                                                            {
                                                                const double numer((P2-P1).dot(N2));
                                                                const double u(numer/denom);
                                                                const VectorDimD dP=(N2-N2.dot(N1)*N1)*u;
                                                                if (dP.norm()<dx)
                                                                {
                                                                    const VectorDimD linePoint = P1+dP;
                                                                    
                                                                    std::cout<<"formVertexEdgeJunctions: vertex "<<nodeIter->second->sID <<
                                                                    ", link "<<linkIter->second->source->sID<<"->"<<linkIter->second->sink->sID<<" @ "<< u <<std::endl;
                                                                }
                                                            }
                                                            break;
                                                        }
                                                            
                                                    }
                                                } // frank-rule appies
                                            } // root is closer than dx
                                        } // there is at least one root
                                    } // vertex is not too far from edge end points
                                } // vertex is not one of the edge end points
                            } // loop on edges
                        } // vertex has only one normal
                    } // vertex has one edge in, one edge out, and it is balanced
                } // vertex has inverted its motion
            } // loop on vertices
            
            return 1;
        }
        
		
		/* formJunctions ******************************************************/
		void formJunctions(const double& dx, const double& avoidNodeIntersection)
        {
			
			
            
            formEdgeEdgeJunctions(dx,avoidNodeIntersection);
            
            if (useVertexEdgeJunctions)
            {
                formVertexEdgeJunctions(dx);
            }
            
		} // end formJunctions
        
        
        
        
        
	};
    
    
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
    
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    bool DislocationJunctionFormation<DislocationNetworkType>::useVertexEdgeJunctions=true;
    
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif





//                    {
//                        for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
//                        {
//
//                            if (nodeIter->first != linkIter->first.first && nodeIter->first != linkIter->first.second) // avoid intersectio with ajecent edges
//                            {
//
//                                const VectorDimD N1(vertexPN[0]);
//                                const VectorDimD N2(linkIter->second->glidePlaneNormal);
//                                const VectorDimD P1(nodeIter->second->get_P());
//                                const VectorDimD P2(linkIter->second->source->get_P());
//
//                               const int ppt(PlanePlaneIntersection<dim>::planePlaneType(N1,N2,P1,P2));
//
//
//                                switch (ppt)
//                                {
//                                    case 0: // paralle planes
//                                        // do nothing, no intersection possile
//                                        break;
//
//                                    case 1: // coincident planes
//                                    {
//
//                                        break;
//                                    }
//
//                                    default: // =2 incident planes
//                                    {
//                                        const double denom(1.0-std::pow(N1.dot(N2),2));
//                                        if(std::fabs(denom)>FLT_EPSILON)
//                                        {
//                                            const double numer((P2-P1).dot(N2));
//                                            const double u(numer/denom);
//                                            const VectorDimD dP=(N2-N2.dot(N1)*N1)*u;
//                                            if (dP.norm()<dx)
//                                            {
//                                                const VectorDimD V = P1+dP; // the vertex point projected to the common line
//                                                const VectorDimD lineDir(N1.cross(N2).normalized());
//                                                double segNorm(linkIter->second->chord().norm());
////                                                VectorDimD V1(V+lineDir*2.0*segNorm);
////                                                VectorDimD V2(V-lineDir*2.0*segNorm);
//                                                VectorDimD V1(V+lineDir*2.0*dx);
//                                                VectorDimD V2(V-lineDir*2.0*dx);
//
//                                                enum {polyDegree=3};
//
//                                                const Eigen::Matrix<double,dim,dim> R(DislocationLocalReference<dim>::global2local(linkIter->second->chord(),N2));
//                                                PlanarSplineImplicitization<polyDegree> sli(linkIter->second->polynomialLocalCoeff());
////                                                PlanarSplineImplicitization<polyDegree>::physTol=physTol;
//
//
//                                                Eigen::Matrix<double,dim-1,2> H2L;
//                                                H2L.col(0)=(R*(V1-P2)).template segment<dim-1>(0); // coordinate of V1 in the local ref. system of the spline
//                                                H2L.col(1)=(R*(V2-P2)).template segment<dim-1>(0); // coordinate of V2 in the local ref. system of the spline
//
//
//                                                std::set<std::pair<double,double> > lineIntersectionParameters=sli.template intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L));
//
//                                                std::set<std::pair<double,double> > sortedDistance;
//
//                                                for (std::set<std::pair<double,double> >::const_iterator sIter=lineIntersectionParameters.begin(); sIter!=lineIntersectionParameters.end();++sIter)
//                                                {
//                                                    VectorDimD temp(V1*(1.0-sIter->second)+V2*sIter->second);
//
//                                                    const double tempNorm((V-temp).norm());
//                                                    if (tempNorm<dx)
//                                                    {
//                                                        sortedDistance.insert(std::make_pair(tempNorm,sIter->first));
//                                                    }
//
//                                                }
//
//                                                if (sortedDistance.size()>0)
//                                                {
//                                                    std::cout<<"Vertex-Edge intersection: node "<< nodeIter->second->sID
//                                                    <<", edge "<<linkIter->second->source->sID<<"->"<<linkIter->second->sink->sID<<" @ "<<sortedDistance.begin()->second<<std::endl;
//                                                }
//                                                else
//                                                {
//
//                                                }
//
//                                            }
//                                        }
//
//
//                                        break;
//                                    }
//                                }
//
//                            }
//
//
//                            //if (linkIterA->second->sID!=linkIterB->second->sID)
//                        }
//
//
//                    }




//            // Direct Node contraction
//            std::vector<std::pair<size_t,std::pair<size_t,size_t> > > nodeContractVector;
//            //            typename EdgeFinder<LinkType>::isNetworkEdgeType
//            for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
//            {
//                bool isSimple(nodeIter->second->is_simple());
//                if(isSimple) // has two neighbors
//                {
//                    const NodeType* const pNi(nodeIter->second->openNeighborNode(0));
//                    const NodeType* const pNj(nodeIter->second->openNeighborNode(1));
//                    if(pNi->openOrder()>2 && pNj->openOrder()>2)
//                    {
//                        nodeContractVector.push_back(std::make_pair(nodeIter->second->sID,std::make_pair(pNi->sID,pNj->sID)));
//                    }
//                }
//
//            }
//
//
//            for (unsigned int m=0;m<nodeContractVector.size();++m)
//            {
//                const size_t k(nodeContractVector[m].first);
//                const size_t i(nodeContractVector[m].second.first);
//                const size_t j(nodeContractVector[m].second.second);
//
//
//
//                const isNetworkLinkType Lki(DN.link(k,i));
//                if(Lki.first)
//                {
//                    DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(Lki);
//                    //                    nodeContractVector.push_back(Lki);
//                }
//                else
//                {
//                    DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(DN.link(i,k));
//                    //nodeContractVector.push_back(DN.link(i,k));
//                }
//
//                const isNetworkLinkType Lkj(DN.link(k,j));
//                if(Lkj.first)
//                {
//                    DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(Lkj);
//                    //nodeContractVector.push_back(Lkj);
//                }
//                else
//                {
//                    DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(DN.link(j,k));
//                    //                    nodeContractVector.push_back(DN.link(j,k));
//                }
//
//                const isNetworkLinkType Lij(DN.link(i,j));
//                if(Lij.first)
//                {
//                    DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(Lij);
//                    //                    nodeContractVector.push_back(Lij);
//                }
//                else
//                {
//                    DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(DN.link(j,i));
//                    //                    nodeContractVector.push_back(DN.link(j,i));
//                }
//
//
//                //                DislocationNetworkRemesh<DislocationNetworkType>(DN).singleEdgeContract(nodeContractVector[k]);
//                //DN.contractSecond(nodeContractVector[k].first,nodeContractVector[k].second);
//            }
