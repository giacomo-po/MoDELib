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
		
        //		typedef std::pair<const LinkType*, double> EdgeIntersectionType;
		typedef std::pair<std::pair<size_t,size_t>, double> EdgeIntersectionType;
        
		typedef std::pair<EdgeIntersectionType,EdgeIntersectionType> EdgeIntersectionPairType;
        //		typedef std::vector<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
		typedef std::deque<EdgeIntersectionPairType> EdgeIntersectionPairContainerType;
		
        
        /**********************************************************************/
        //		int junctionDir(const EdgeIntersectionPairType& intersectionPair) const
		int junctionDir(const LinkType& L1, const LinkType& L2,
                        const double& u1, const double& u2) const
        {
            //            const double u1(intersectionPair.first.second);
            //            const double u2(intersectionPair.second.second);
            const VectorDimD b1(L1.flow);
            const VectorDimD b2(L2.flow);
            const VectorDimD rl1(L1.get_rl(u1));
            const VectorDimD rl2(L2.get_rl(u2));
            
            
            const bool is1sessile(L1.sessilePlaneNormal.norm()>FLT_EPSILON);
            const bool is2sessile(L2.sessilePlaneNormal.norm()>FLT_EPSILON);
            
            VectorDimD prjDir(VectorDimD::Zero());
            if (!is1sessile && !is2sessile){
                const VectorDimD commonLine(L1.glidePlaneNormal.cross(L2.glidePlaneNormal));
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
            else if(is1sessile && !is2sessile)
            { // use chord of I
                prjDir=L1.chord().normalized();
            }
            else if(!is1sessile && is2sessile)
            { // use chord of J
                prjDir=L2.chord().normalized();
            }
            else
            {
                assert(0 && "CANNOT DETERMINE COMMON LINE BETWEEN TWO SESSILE SEGMENTS.");
            }
            
            
            const double sgnrl1rl2(rl1.dot(prjDir)*rl2.dot(prjDir));
            
            
            //const bool frankRule(b1.dot(b2)*rl1.dot(rl2)<=0.0);
            const bool frankRule(b1.dot(b2)*sgnrl1rl2<=0.0);
            const bool isValidJunction(frankRule || L1.is_boundarySegment() || L2.is_boundarySegment());
            
            int dir(0);
            if (isValidJunction)
            {
                dir=( sgnrl1rl2 > 0.0) ? 1 : ((sgnrl1rl2 < 0.0) ? -1 : 0);
            }
            return dir;
        }
        
        
        /**********************************************************************/
        size_t contractWithPlaneCheck(const size_t& i, const size_t& j)
        {
            
            size_t contracted(0);
            
            const auto N1=DN.node(i);
            const auto N2=DN.node(j);
            assert(N1.first && "NODE i DOES NOT EXIST.");
            assert(N2.first && "NODE j DOES NOT EXIST.");
            
            const VectorDimD P1=N1.second->get_P();
            const VectorDimD P2=N2.second->get_P();
            
            
            const typename DislocationNetworkType::NodeType::VectorOfNormalsType PN1(GramSchmidt<dim>(N1.second->constraintNormals()));
            const typename DislocationNetworkType::NodeType::VectorOfNormalsType PN2(GramSchmidt<dim>(N2.second->constraintNormals()));
            const size_t sizePN1(PN1.size());
            const size_t sizePN2(PN2.size());
            
            
            if(sizePN1==1 && sizePN2==1) // nodes constrained to move on planes
            {
//                const VectorDimD d(PN1[0].cross(PN2[0])); // direction of common line

                std::cout<<"contractWithPlaneCheck, case 1"<<std::endl;
                const double denom(1.0-std::pow(PN1[0].dot(PN2[0]),2));
                const double numer((P2-P1).dot(PN2[0]));

                if(denom<FLT_EPSILON)
                {
                    if(std::fabs(denom)<FLT_EPSILON) // planes are coincident
                    {
                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                    }
                    else // parallel planes
                    {
                        assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
                    }
                }
                else // incident planes
                {
                    contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P1+(PN2[0]-PN2[0].dot(PN1[0])*PN1[0])*numer/denom);
                }
            }
            else if(sizePN1==1 && sizePN2==2) // N1 moves on a plane, N2 moves on a line
            {
                std::cout<<"contractWithPlaneCheck, case 2"<<std::endl;
                const VectorDimD d2(PN2[0].cross(PN2[1]));
                const double den(d2.dot(PN1[0]));
                const double num((P1-P2).dot(PN1[0]));
                if(std::fabs(den)>FLT_EPSILON)
                {
                    contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P2+num/den*d2);
                }
                else
                {
                    if(std::fabs(num)<FLT_EPSILON && (P1-P2).norm()<FLT_EPSILON)
                    {
                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                    }
//                    else
//                    {
//                        assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
//                    }
                }
            }
            else if(sizePN1==2 && sizePN2==1) // N1 moves on a line, N2 moves on a plane
            {
                std::cout<<"contractWithPlaneCheck, case 3"<<std::endl;
                const VectorDimD d1(PN1[0].cross(PN1[1]));
                const double den(d1.dot(PN2[0]));
                const double num((P2-P1).dot(PN2[0]));
                if(std::fabs(den)>FLT_EPSILON)
                {
                    contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P1+num/den*d1);
                }
                else
                {
                    if(std::fabs(num)<FLT_EPSILON && (P1-P2).norm()<FLT_EPSILON)
                    {
                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                    }
//                    else
//                    {
//                        assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
//                    }
                }
            }
            else if(sizePN1==2 && sizePN2==2) // both N1 and N2 move on lines
            {
                std::cout<<"contractWithPlaneCheck, case 4"<<std::endl;
                
                const double P12norm((P1-P2).norm());
                if(P12norm<FLT_EPSILON)
                {
                    contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                }
                else
                {
                    const VectorDimD d1(PN1[0].cross(PN1[1]));
                    const VectorDimD d2(PN2[0].cross(PN2[1]));
                    if(std::fabs(d1.cross(d2).dot((P1-P2)/P12norm))<FLT_EPSILON) // planarity condition
                    {
                        const VectorDimD dOrth=d2-d2.dot(d1)*d1;
                        const double den=d2.dot(dOrth);
                        if(std::fabs(den)>FLT_EPSILON)
                        {
                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P2+(P1-P2).dot(dOrth)/den*d2);
                        }
                        else
                        {
                            assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
                        }
                    }
                    else
                    {
                        assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
                    }
                }
            }
            else if(sizePN1==1 && sizePN2==3)
            {
                if(std::fabs((P2-P1).dot(PN1[0]))<FLT_EPSILON) // P2 belongs to the plane of P1
                {// contract N1
                    contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                }
            }
            else if(sizePN1==3 && sizePN2==1)
            {
                if(std::fabs((P1-P2).dot(PN2[0]))<FLT_EPSILON) // P1 belongs to the plane of P2
                {// contract N2
                    contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                }
            }
            else if(sizePN1==2 && sizePN2==3)
            {
                
            }
            else if(sizePN1==3 && sizePN2==2)
            {
                
            }
            else if(sizePN1==3 && sizePN2==3)
            {
                
            }
            
            return contracted;
            

            
            
        }
		
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
		void findIntersections(EdgeIntersectionPairContainerType& intersectionContainer,
                                       std::deque<int>& dirVector,
                                       const double& avoidNodeIntersection) const
        
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
            
            //			EdgeIntersectionPairContainerType intersectionContainer;
			
			
			//! 2- loop over all links and determine their intersections
			for (typename NetworkLinkContainerType::const_iterator linkIterA=DN.linkBegin();linkIterA!=DN.linkEnd();linkIterA++)
            {
                const DislocationSegmentIntersection<dim,pOrder> dsi(linkIterA->second-> hermiteCoefficients(),linkIterA->second->glidePlaneNormal);
                
				for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.linkEnd();linkIterB++)
                {
					if (linkIterA->second->sID!=linkIterB->second->sID) // don't intersect with itself
                    {
                        const bool L1isSessile(linkIterA->second->sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
                        const bool L2isSessile(linkIterB->second->sessilePlaneNormal.squaredNorm()>FLT_EPSILON);
                        
                        std::set<std::pair<double,double> > temp; // the container of the roots
                        
                        if (!L1isSessile && !L2isSessile) // both are glissile
                        {
                            temp = dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                        }
                        
                        else if (!L1isSessile && L2isSessile) // L1 is glissile and L2 is sessile
                        {
                            const bool gnAgnB((linkIterA->second->glidePlaneNormal-linkIterB->second->glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                            const bool gnAsnB((linkIterA->second->glidePlaneNormal-linkIterB->second->sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                            
                            if (!gnAgnB && !gnAsnB)
                            {
                                // cannot intersect
                            }
                            else if(gnAgnB && !gnAsnB)
                            { // use planeNormal of A and planeNormal of B
                                temp = dsi.intersectWith(linkIterB->second-> hermiteCoefficients(),linkIterB->second->glidePlaneNormal,collisionTol);
                            }
                            else if (!gnAgnB && gnAsnB)
                            { // use planeNormal of A and sessileNormal of B
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
                            
                            if (!gnBgnA && !gnBsnA)
                            {
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
                            //                            if (   paramIter->first >avoidNodeIntersection && paramIter-> first<1.0-avoidNodeIntersection // THIS IS CHECKED LATER
                            //                                && paramIter->second>avoidNodeIntersection && paramIter->second<1.0-avoidNodeIntersection) // avoid node intersection
                            //                            {
                            
                            const bool intersectionIsSourceSource(   paramIter->first  <0.5
                                                                   && paramIter->second <0.5
                                                                   &&
                                                                  linkIterA->second->source->sID==linkIterB->second->source->sID);
                            
                            const bool   intersectionIsSourceSink(   paramIter->first  <0.5
                                                                   && paramIter->second >1.0-0.5
                                                                   &&
                                                                  linkIterA->second->source->sID==linkIterB->second->sink->sID);
                            
                            const bool   intersectionIsSinkSource(   paramIter->first  > 1.0-0.5
                                                                   && paramIter->second <0.5
                                                                   &&
                                                                  linkIterA->second->sink->sID==linkIterB->second->source->sID);
                            
                            const bool     intersectionIsSinkSink(   paramIter->first  > 1.0-0.5
                                                                   && paramIter->second > 1.0-0.5
                                                                   &&
                                                                  linkIterA->second->sink->sID==linkIterB->second->sink->sID);
                            
                            
                            if(!intersectionIsSourceSource && !intersectionIsSourceSink && !intersectionIsSinkSource && !intersectionIsSinkSink)
                            {
                                EdgeIntersectionType intersectionOnA(std::make_pair(linkIterA->second->nodeIDPair,paramIter->first ));
                                EdgeIntersectionType intersectionOnB(std::make_pair(linkIterB->second->nodeIDPair,paramIter->second));
                                const int dir(junctionDir(*linkIterA->second,*linkIterB->second,paramIter->first,paramIter->second));
                                if(dir!=0)
                                {
                                    intersectionContainer.emplace_back(intersectionOnA,intersectionOnB);
                                    dirVector.emplace_back(dir);
                                }
                            }
                        } // end for
                    } // end don't self-intersect
				} // end loop over second segment
			} // end loop over first segment
            //			return intersectionContainer;
		}
		
        /* formEdgeEdgeJunctions **********************************************/
        void formJunctions(const double& dx, const double& avoidNodeIntersection)
        {
            //! 1- Initialize intersectionContainer calling findIntersections deque<Pair<Pair<Link*double>,Pair<Link*,double>>>
            //			EdgeIntersectionPairContainerType intersectionContainer(findIntersections(avoidNodeIntersection));
			EdgeIntersectionPairContainerType intersectionContainer;
			std::deque<int> dirVector;
            findIntersections(intersectionContainer,dirVector,avoidNodeIntersection);
			assert(intersectionContainer.size()==dirVector.size());
            //			model::cout<<intersectionContainer.size()<<" geometric) ";
            model::cout<<intersectionContainer.size()<<" physical intersections. ";
            
            typedef std::pair<size_t,size_t> EdgeIDType;
            
			
            //			EdgeIntersectionContainerType edgeIntersectionContainer;
			for (size_t interID=0;interID!=intersectionContainer.size();++interID)
            {
                const EdgeIDType& key1(intersectionContainer[interID]. first.first);
				const EdgeIDType& key2(intersectionContainer[interID].second.first);
                
                const isNetworkLinkType L1(DN.link(key1.first,key1.second));
                const isNetworkLinkType L2(DN.link(key2.first,key2.second));
                
                if(L1.first && L2.first) // Links exist
                {
                    const NodeType& source1=*(L1.second->source);
                    const NodeType&   sink1=*(L1.second->sink);
                    const NodeType& source2=*(L2.second->source);
                    const NodeType&   sink2=*(L2.second->sink);
                    const double du1(dx/L1.second->chordLength());
                    const double du2(dx/L2.second->chordLength());
                    const VectorDimD N1(L1.second->glidePlaneNormal);
                    const VectorDimD N2(L2.second->glidePlaneNormal);
//                    const bool is1sessile(L1.second->sessilePlaneNormal.norm()>FLT_EPSILON);
//                    const bool is2sessile(L2.second->sessilePlaneNormal.norm()>FLT_EPSILON);
                    const Simplex<dim,dim>* S1(source1.includingSimplex());
                    const Simplex<dim,dim>* S2(source2.includingSimplex());
                    
                    
                    
                    if (   intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection     // intersection on first link is far enough from source
                        && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection   // intersection on first link is far enough from sink
                        && intersectionContainer[interID].second.second-du2 > avoidNodeIntersection     // intersection on second link is far enough from source
                        && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection   // intersection on second link is far enough from sink
                        //                        && dirVector[interID]!=0
                        )
                    {
                        
                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        const double u2m(intersectionContainer[interID].second.second-du2);
                        const VectorDimD P1m(L1.second->get_r(u1m));
                        const VectorDimD P2m(L2.second->get_r(u2m));
                        
                        const double u1p(intersectionContainer[interID]. first.second+du1);
                        const double u2p(intersectionContainer[interID].second.second+du2);
                        const VectorDimD P1p(L1.second->get_r(u1p));
                        const VectorDimD P2p(L2.second->get_r(u2p));
                        
                        
                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1m),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1m,S1).first;
                            if(firstIntersectionInsideMesh)
                            {
                                //                                firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2m),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                                firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2m,S2).first;
                            }
                        }
                        
                        bool secondIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1p),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1p,L1.second->source->includingSimplex()).first;
                            if(secondIntersectionInsideMesh)
                            {
                                //                                secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2p),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                                secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2p,L2.second->source->includingSimplex()).first;
                            }
                        }
                        
                        if(firstIntersectionInsideMesh && secondIntersectionInsideMesh)
                        {
                            
                            model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                            model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                            std::cout<<"Junction case 1"<<std::endl;

                            const std::pair<bool,size_t> success1m=DN.expand(source1.sID,sink1.sID,P1m); // now L1.second is invalid
                            assert(success1m.first && "COULD NOT EXPLAND LINK1 AT LOWER INTERSECTION");
                            const size_t& im=success1m.second; // id of the node obtained expanding L1
                            
                            const std::pair<bool,size_t> success2m=DN.expand(source2.sID,sink2.sID,P2m); // now L2.second is invalid
                            assert(success2m.first && "COULD NOT EXPLAND LINK2 AT LOWER INTERSECTION");
                            const size_t& jm=success2m.second; // id of the node obtained expanding L2
                            
                            const std::pair<bool,size_t> success1p=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
                            assert(success1p.first && "COULD NOT EXPLAND LINK1 AT UPPER INTERSECTION");
                            const size_t& ip=success1p.second; // id of the node obtained expanding L1
                            
                            const std::pair<bool,size_t> success2p=DN.expand(jm,sink2.sID,P2p); // now L2.second is invalid
                            assert(success2p.first && "COULD NOT EXPLAND LINK2 AT UPPER INTERSECTION");
                            const size_t& jp=success2p.second; // id of the node obtained expanding L2
                            
                            switch (dirVector[interID])
                            {
                                case +1:
                                {
                                    contractWithPlaneCheck(im,jm);
                                    contractWithPlaneCheck(ip,jp);
                                    break;
                                }
                                    
                                case -1:
                                {
                                    contractWithPlaneCheck(im,jp);
                                    contractWithPlaneCheck(ip,jm);
                                    break;
                                }
                                    
                                default:
                                    assert(0 && "DIR VECTOR CAN ONLY BE +1 or -1");
                                    break;
                            }
                            
                        }
                        
                    }
                    else if (   intersectionContainer[interID]. first.second-du1 < avoidNodeIntersection     // intersection on first link is near source
                             && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection        // intersection on first link is far enough from sink
                             && intersectionContainer[interID].second.second-du2 > avoidNodeIntersection          // intersection on second link is far enough from source
                             && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection        // intersection on second link is far enough from sink
                             )
                    {
                        //                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        const double u2m(intersectionContainer[interID].second.second-du2);
                        //                        const VectorDimD P1m(L1.second->get_r(u1m));
                        const VectorDimD P2m(L2.second->get_r(u2m));
                        
                        const double u1p(intersectionContainer[interID]. first.second+du1);
                        const double u2p(intersectionContainer[interID].second.second+du2);
                        const VectorDimD P1p(L1.second->get_r(u1p));
                        const VectorDimD P2p(L2.second->get_r(u2p));
                        
                        
                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2m,S2).first;
                        }
                        
                        bool secondIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1p,L1.second->source->includingSimplex()).first;
                            if(secondIntersectionInsideMesh)
                            {
                                secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2p,L2.second->source->includingSimplex()).first;
                            }
                        }
                        
                        if(firstIntersectionInsideMesh && secondIntersectionInsideMesh)
                        {
                            
                            
                            //                            const std::pair<bool,size_t> success1m=DN.expand(source1.sID,sink1.sID,P1m); // now L1.second is invalid
                            //                            assert(success1m.first && "COULD NOT EXPLAND LINK1 AT LOWER INTERSECTION");
                            const size_t& im=source1.sID; // id of the node obtained expanding L1
                            
                            model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                            model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                            std::cout<<"Junction case 2: node "<<im<<" survives"<<std::endl;
                            
                            const std::pair<bool,size_t> success2m=DN.expand(source2.sID,sink2.sID,P2m); // now L2.second is invalid
                            assert(success2m.first && "COULD NOT EXPLAND LINK2 AT LOWER INTERSECTION");
                            const size_t& jm=success2m.second; // id of the node obtained expanding L2
                            
                            const std::pair<bool,size_t> success1p=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
                            assert(success1p.first && "COULD NOT EXPLAND LINK1 AT UPPER INTERSECTION");
                            const size_t& ip=success1p.second; // id of the node obtained expanding L1
                            
                            const std::pair<bool,size_t> success2p=DN.expand(jm,sink2.sID,P2p); // now L2.second is invalid
                            assert(success2p.first && "COULD NOT EXPLAND LINK2 AT UPPER INTERSECTION");
                            const size_t& jp=success2p.second; // id of the node obtained expanding L2
                            
                            switch (dirVector[interID])
                            {
                                case +1:
                                {
                                    contractWithPlaneCheck(im,jm);
                                    contractWithPlaneCheck(ip,jp);
                                    break;
                                }
                                    
                                case -1:
                                {
                                    contractWithPlaneCheck(im,jp);
                                    contractWithPlaneCheck(ip,jm);
                                    break;
                                }
                                    
                                default:
                                    assert(0 && "DIR VECTOR CAN ONLY BE +1 or -1");
                                    break;
                            }
                            
                        }
                        
                    }
                    else if (   intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection     // intersection on first link is far enough from source
                             && intersectionContainer[interID]. first.second+du1 >1.0-avoidNodeIntersection   // intersection on first link is near sink
                             && intersectionContainer[interID].second.second-du2 > avoidNodeIntersection     // intersection on second link is far enough from source
                             && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection   // intersection on second link is far enough from sink
                             //                        && dirVector[interID]!=0
                             )
                    {
                        //                        model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                        //                        model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                        
                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        const double u2m(intersectionContainer[interID].second.second-du2);
                        const VectorDimD P1m(L1.second->get_r(u1m));
                        const VectorDimD P2m(L2.second->get_r(u2m));
                        
                        //                        const double u1p(intersectionContainer[interID]. first.second+du1);
                        const double u2p(intersectionContainer[interID].second.second+du2);
                        //                        const VectorDimD P1p(L1.second->get_r(u1p));
                        const VectorDimD P2p(L2.second->get_r(u2p));
                        
                        
                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1m),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1m,S1).first;
                            if(firstIntersectionInsideMesh)
                            {
                                //                                firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2m),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                                firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2m,S2).first;
                            }
                        }
                        
                        bool secondIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1p),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            //                            secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1p,L1.second->source->includingSimplex()).first;
                            //                            if(secondIntersectionInsideMesh)
                            //                            {
                            //                                secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2p),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                            secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2p,L2.second->source->includingSimplex()).first;
                            //                            }
                        }
                        
                        if(firstIntersectionInsideMesh && secondIntersectionInsideMesh)
                        {
                            
                            const std::pair<bool,size_t> success1m=DN.expand(source1.sID,sink1.sID,P1m); // now L1.second is invalid
                            assert(success1m.first && "COULD NOT EXPLAND LINK1 AT LOWER INTERSECTION");
                            const size_t& im=success1m.second; // id of the node obtained expanding L1
                            
                            const std::pair<bool,size_t> success2m=DN.expand(source2.sID,sink2.sID,P2m); // now L2.second is invalid
                            assert(success2m.first && "COULD NOT EXPLAND LINK2 AT LOWER INTERSECTION");
                            const size_t& jm=success2m.second; // id of the node obtained expanding L2
                            
                            //                            const std::pair<bool,size_t> success1p=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
                            //                            assert(success1p.first && "COULD NOT EXPLAND LINK1 AT UPPER INTERSECTION");
                            const size_t& ip=sink1.sID; // id of the node obtained expanding L1
                            
                            model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                            model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                            std::cout<<"Junction case 3: node "<<ip<<" survives"<<std::endl;
                            
                            const std::pair<bool,size_t> success2p=DN.expand(jm,sink2.sID,P2p); // now L2.second is invalid
                            assert(success2p.first && "COULD NOT EXPLAND LINK2 AT UPPER INTERSECTION");
                            const size_t& jp=success2p.second; // id of the node obtained expanding L2
                            
                            switch (dirVector[interID])
                            {
                                case +1:
                                {
                                    contractWithPlaneCheck(im,jm);
                                    contractWithPlaneCheck(ip,jp);
                                    break;
                                }
                                    
                                case -1:
                                {
                                    contractWithPlaneCheck(im,jp);
                                    contractWithPlaneCheck(ip,jm);
                                    break;
                                }
                                    
                                default:
                                    assert(0 && "DIR VECTOR CAN ONLY BE +1 or -1");
                                    break;
                            }
                            
                        }
                        
                    }
                    else if (   intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection     // intersection on first link is far enough from source
                             && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection   // intersection on first link is far enough from sink
                             && intersectionContainer[interID].second.second-du2 < avoidNodeIntersection     // intersection on second link near source
                             && intersectionContainer[interID].second.second+du2<1.0-avoidNodeIntersection   // intersection on second link is far enough from sink
                             //                        && dirVector[interID]!=0
                             )
                    {
                        //                        model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                        //                        model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                        
                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        //                        const double u2m(intersectionContainer[interID].second.second-du2);
                        const VectorDimD P1m(L1.second->get_r(u1m));
                        //                        const VectorDimD P2m(L2.second->get_r(u2m));
                        
                        const double u1p(intersectionContainer[interID]. first.second+du1);
                        const double u2p(intersectionContainer[interID].second.second+du2);
                        const VectorDimD P1p(L1.second->get_r(u1p));
                        const VectorDimD P2p(L2.second->get_r(u2p));
                        
                        
                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1m),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1m,S1).first;
                            //                            if(firstIntersectionInsideMesh)
                            //                            {
                            //                                //                                firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2m),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                            //                                firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2m,S2).first;
                            //                            }
                        }
                        
                        bool secondIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1p),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1p,L1.second->source->includingSimplex()).first;
                            if(secondIntersectionInsideMesh)
                            {
                                //                                secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2p),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                                secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2p,L2.second->source->includingSimplex()).first;
                            }
                        }
                        
                        if(firstIntersectionInsideMesh && secondIntersectionInsideMesh)
                        {
                            
                            const std::pair<bool,size_t> success1m=DN.expand(source1.sID,sink1.sID,P1m); // now L1.second is invalid
                            assert(success1m.first && "COULD NOT EXPLAND LINK1 AT LOWER INTERSECTION");
                            const size_t& im=success1m.second; // id of the node obtained expanding L1
                            
                            //                            const std::pair<bool,size_t> success2m=DN.expand(source2.sID,sink2.sID,P2m); // now L2.second is invalid
                            //                            assert(success2m.first && "COULD NOT EXPLAND LINK2 AT LOWER INTERSECTION");
                            const size_t& jm=source2.sID; // id of the node obtained expanding L2
                            
                            model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                            model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                            std::cout<<"Junction case 4: node "<<jm<<" survives"<<std::endl;
                            
                            
                            const std::pair<bool,size_t> success1p=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
                            assert(success1p.first && "COULD NOT EXPLAND LINK1 AT UPPER INTERSECTION");
                            const size_t& ip=success1p.second; // id of the node obtained expanding L1
                            
                            const std::pair<bool,size_t> success2p=DN.expand(jm,sink2.sID,P2p); // now L2.second is invalid
                            assert(success2p.first && "COULD NOT EXPLAND LINK2 AT UPPER INTERSECTION");
                            const size_t& jp=success2p.second; // id of the node obtained expanding L2
                            
                            switch (dirVector[interID])
                            {
                                case +1:
                                {
                                    contractWithPlaneCheck(im,jm);
                                    contractWithPlaneCheck(ip,jp);
                                    break;
                                }
                                    
                                case -1:
                                {
                                    contractWithPlaneCheck(im,jp);
                                    contractWithPlaneCheck(ip,jm);
                                    break;
                                }
                                    
                                default:
                                    assert(0 && "DIR VECTOR CAN ONLY BE +1 or -1");
                                    break;
                            }
                            
                        }
                        
                    }
                    else if (   intersectionContainer[interID]. first.second-du1 > avoidNodeIntersection     // intersection on first link is far enough from source
                             && intersectionContainer[interID]. first.second+du1<1.0-avoidNodeIntersection   // intersection on first link is far enough from sink
                             && intersectionContainer[interID].second.second-du2 > avoidNodeIntersection     // intersection on second link is far enough from source
                             && intersectionContainer[interID].second.second+du2>1.0-avoidNodeIntersection   // intersection on second link is near sink
                             //                        && dirVector[interID]!=0
                             )
                    {
                        //                        model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                        //                        model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                        
                        const double u1m(intersectionContainer[interID]. first.second-du1);
                        const double u2m(intersectionContainer[interID].second.second-du2);
                        const VectorDimD P1m(L1.second->get_r(u1m));
                        const VectorDimD P2m(L2.second->get_r(u2m));
                        
                        const double u1p(intersectionContainer[interID]. first.second+du1);
                        //                        const double u2p(intersectionContainer[interID].second.second+du2);
                        const VectorDimD P1p(L1.second->get_r(u1p));
                        //                        const VectorDimD P2p(L2.second->get_r(u2p));
                        //
                        
                        bool firstIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1m),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1m,S1).first;
                            if(firstIntersectionInsideMesh)
                            {
                                //                                firstIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2m),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                                firstIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2m,S2).first;
                            }
                        }
                        
                        bool secondIntersectionInsideMesh(true);
                        if (DN.shared.use_boundary)
                        {
                            //                            secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L1.second->get_r(u1p),L1.second->source->includingSimplex(),FLT_EPSILON).first;
                            secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P1p,L1.second->source->includingSimplex()).first;
                            //                            if(secondIntersectionInsideMesh)
                            //                            {
                            //                                //                                secondIntersectionInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(L2.second->get_r(u2p),L2.second->source->includingSimplex(),FLT_EPSILON).first;
                            //                                secondIntersectionInsideMesh*=DN.shared.mesh.searchWithGuess(P2p,L2.second->source->includingSimplex()).first;
                            //                            }
                        }
                        
                        if(firstIntersectionInsideMesh && secondIntersectionInsideMesh)
                        {
                            
                            const std::pair<bool,size_t> success1m=DN.expand(source1.sID,sink1.sID,P1m); // now L1.second is invalid
                            assert(success1m.first && "COULD NOT EXPLAND LINK1 AT LOWER INTERSECTION");
                            const size_t& im=success1m.second; // id of the node obtained expanding L1
                            
                            const std::pair<bool,size_t> success2m=DN.expand(source2.sID,sink2.sID,P2m); // now L2.second is invalid
                            assert(success2m.first && "COULD NOT EXPLAND LINK2 AT LOWER INTERSECTION");
                            const size_t& jm=success2m.second; // id of the node obtained expanding L2
                            
                            const std::pair<bool,size_t> success1p=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
                            assert(success1p.first && "COULD NOT EXPLAND LINK1 AT UPPER INTERSECTION");
                            const size_t& ip=success1p.second; // id of the node obtained expanding L1
                            
                            //                            const std::pair<bool,size_t> success2p=DN.expand(jm,sink2.sID,P2p); // now L2.second is invalid
                            //                            assert(success2p.first && "COULD NOT EXPLAND LINK2 AT UPPER INTERSECTION");
                            const size_t& jp=sink2.sID; // id of the node obtained expanding L2
                            
                            
                            model::cout<<"key1 is "<<key1.first<<"->"<<key1.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID]. first.second<<std::endl;
                            model::cout<<"key2 is "<<key2.first<<"->"<<key2.second<<" at "<<std::setprecision(15)<<intersectionContainer[interID].second.second<<std::endl;
                            std::cout<<"Junction case 5: node "<<jp<<" survives"<<std::endl;
                            
                            
                            
                            switch (dirVector[interID])
                            {
                                case +1:
                                {
                                    contractWithPlaneCheck(im,jm);
                                    contractWithPlaneCheck(ip,jp);
                                    break;
                                }
                                    
                                case -1:
                                {
                                    contractWithPlaneCheck(im,jp);
                                    contractWithPlaneCheck(ip,jm);
                                    break;
                                }
                                    
                                default:
                                    assert(0 && "DIR VECTOR CAN ONLY BE +1 or -1");
                                    break;
                            }
                            
                        }
                        
                    }
                    
                }
                
			}
			
        }
        
	};
    
    
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
    
    
    // Declare Static Data
    //    template <typename DislocationNetworkType>
    //    bool DislocationJunctionFormation<DislocationNetworkType>::useVertexEdgeJunctions=true;
    
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif



//            if (!is1sessile && !is2sessile)
//            {
//                VectorDimD linePoint((P1+P2)*0.5);
//                const double denom(1.0-std::pow(N1.dot(N2),2));
//                if(std::fabs(denom)>FLT_EPSILON)
//                { // planes are incident: make sure to place intersection on common line
//                    const double numer((P2-P1).dot(N2));
//                    const double u(numer/denom);
//                    linePoint = P1+(N2-N2.dot(N1)*N1)*u;
//                }
//
//                // If mesh is used, make cure that intersection point is inside mesh
//                bool linePointInsideMesh(true);
//                if (DN.shared.use_boundary)
//                {
//                    //                        linePointInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(linePoint,Ni.second->includingSimplex(),FLT_EPSILON).first;
//                    linePointInsideMesh*=DN.shared.mesh.searchWithGuess(linePoint,guess).first;
//                }
//                if(linePointInsideMesh)
//                {
//                    model::cout<<"Glissile Junction"<<std::endl;
//                    DN.contract(i,j,linePoint);
//                }
//
//            }
//            else if(is1sessile && !is2sessile){ // use P1 (which is on sessile segment) as the intersection point
//                model::cout<<"First-Sessile Junction"<<std::endl;
//                DN.contract(i,j,P1);
//            }
//            else if(!is1sessile && is2sessile){ // use P2 (which is on sessile segment) as the intersection point
//                model::cout<<"Second-Sessile Junction"<<std::endl;
//                DN.contract(i,j,P2);
//            }
//            else{
//                assert(0 && "CANNOT MAKE JUNCTION BETWEEN TWO SESSILE SEGMENTS.");
//
//            }


//        /* formVertexEdgeJunctions ********************************************/
//        size_t formVertexEdgeJunctions(const double& dx)
//        {
//
//
//            // Direct Node contraction
//            std::vector<std::pair<size_t,std::pair<size_t,size_t> > > nodeContractVector;
//            //            typename EdgeFinder<LinkType>::isNetworkEdgeType
//
//
//            std::map<size_t,std::pair<std::pair<size_t,size_t>,VectorDimD> > vertexEdgeIntersectionContainer;
//
//            for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
//            {
//                if(nodeIter->second->invertedMotion())
//                {
//                    if(nodeIter->second->outOrder()==1 && nodeIter->second->inOrder()==1 && nodeIter->second->is_balanced() ) // has one inflow and one outflow, therefore it makes sense to talk about B of that node
//                    {
//                        VectorDimD B(std::get<1>(nodeIter->second->outNeighborhood().begin()->second)->flow);
//                        VectorDimD T(std::get<1>(nodeIter->second->outNeighborhood().begin()->second)->sourceT());
//
//                        const typename DislocationNetworkType::NodeType::VectorOfNormalsType vertexPN(GramSchmidt<dim>(nodeIter->second->constraintNormals()));
//
//                        if(vertexPN.size()==1) //node is not constrained, and a unique plane normal si associated with that node
//                        {
//                            for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
//                            {
//
//
//
//                                if (nodeIter->first != linkIter->first.first && nodeIter->first != linkIter->first.second)
//                                {
//
//                                    const double chord2(linkIter->second->chord().squaredNorm());
//                                    if (  (linkIter->second->source->get_P()-nodeIter->second->get_P()).squaredNorm()<chord2    // vertex is not too far from edge source
//                                        ||(linkIter->second->  sink->get_P()-nodeIter->second->get_P()).squaredNorm()<chord2 )  // vertex is not too far from edge sink
//                                    {
//                                        std::map<double,std::pair<double,VectorDimD> > rootMap( linkIter->second->closestPoint(nodeIter->second->get_P() ) );
//
//                                        if (rootMap.size())
//                                        {
//
//                                            const double d(rootMap.begin()->first);
//                                            const double u(rootMap.begin()->second.first);
//
//                                            if (d<dx)
//                                            {
//
//                                                const bool frankRule(B.dot(linkIter->second->flow)*T.dot(linkIter->second->get_ru(u))<=0.0);
//
//                                                if (frankRule)
//                                                {
//
//                                                    //                                            vertexEdgeIntersectionContainer.insert(std::make_pair(nodeIter->second->sID,std::make_pair(,)));
//
//
//                                                    const VectorDimD N1(vertexPN[0]);
//                                                    const VectorDimD N2(linkIter->second->glidePlaneNormal);
//                                                    const VectorDimD P1(nodeIter->second->get_P());
//                                                    const VectorDimD P2(linkIter->second->source->get_P());
//
//                                                    const int ppt(PlanePlaneIntersection<dim>::planePlaneType(N1,N2,P1,P2));
//
//
//                                                    switch (ppt)
//                                                    {
//                                                        case 0: // paralle planes
//                                                            // do nothing, no intersection possile
//                                                            break;
//
//                                                        case 1: // coincident planes
//                                                        {
//
//                                                            break;
//                                                        }
//
//                                                        default: // incident planes
//                                                        {
//                                                            const double denom(1.0-std::pow(N1.dot(N2),2));
//                                                            if(std::fabs(denom)>FLT_EPSILON)
//                                                            {
//                                                                const double numer((P2-P1).dot(N2));
//                                                                const double u(numer/denom);
//                                                                const VectorDimD dP=(N2-N2.dot(N1)*N1)*u;
//                                                                if (dP.norm()<dx)
//                                                                {
//                                                                    const VectorDimD linePoint = P1+dP;
//
//                                                                    std::cout<<"formVertexEdgeJunctions: vertex "<<nodeIter->second->sID <<
//                                                                    ", link "<<linkIter->second->source->sID<<"->"<<linkIter->second->sink->sID<<" @ "<< u <<std::endl;
//                                                                }
//                                                            }
//                                                            break;
//                                                        }
//
//                                                    }
//                                                } // frank-rule appies
//                                            } // root is closer than dx
//                                        } // there is at least one root
//                                    } // vertex is not too far from edge end points
//                                } // vertex is not one of the edge end points
//                            } // loop on edges
//                        } // vertex has only one normal
//                    } // vertex has one edge in, one edge out, and it is balanced
//                } // vertex has inverted its motion
//            } // loop on vertices
//
//            return 1;
//        }


//		/* formJunctions ******************************************************/
//		void formJunctions(const double& dx, const double& avoidNodeIntersection)
//        {
//            formEdgeEdgeJunctions(dx,avoidNodeIntersection);
//
//		} // end formJunctions


//            //			model::cout<<"Using "<< edgeIntersectionContainer.size()/2<<" intersections. ";
//
//			// Call Multiexpand on segments in edgeIntersectionContainer. Store correspoinding new nodes in nodeIntersectionMap
//			std::map<IntersectionIDType, std::set<size_t> > nodeIntersectionMap;
//			for (typename EdgeIntersectionContainerType::const_iterator edgeIter=edgeIntersectionContainer.begin();
//                 /*                                                  */ edgeIter!=edgeIntersectionContainer.end();
//                 /*                                                */ ++edgeIter)
//            {
//				const size_t i(edgeIter->first. first);
//				const size_t j(edgeIter->first.second);
//				std::map<IntersectionIDType,size_t> multiExpandOut=DN.multiExpand(i,j,edgeIter->second);
//				for (std::map<IntersectionIDType,size_t>::const_iterator iter=multiExpandOut.begin();iter!=multiExpandOut.end();++iter)
//                {
//					nodeIntersectionMap[iter->first].insert(iter->second);
//				}
//			}
//
//            // Call Contract on nodes in nodeIntersectionMap
//			for (std::map<IntersectionIDType, std::set<size_t> >::const_iterator mapIter=nodeIntersectionMap.begin();
//                 /*                                                           */ mapIter!=nodeIntersectionMap.end();
//                 /*                                                         */ ++mapIter)
//            {
//				assert(mapIter->second.size()==2 && "THERE SHOULD BE TWO NODES IN EACH INTERSECTION");
//				const size_t i(*(mapIter->second.begin()));
//				const size_t j(*(mapIter->second.rbegin()));
//
//
//				const isNetworkNodeType Ni(DN.node(i));
//				assert(Ni.first);
//				assert(Ni.second->is_simple());
//				const isNetworkNodeType Nj(DN.node(j));
//				assert(Nj.first);
//				assert(Nj.second->is_simple());
//
//				const VectorDimD P1(Ni.second->get_P());
//				const VectorDimD P2(Nj.second->get_P());
//				const VectorDimD N1(Ni.second->openNeighborLink(0)->glidePlaneNormal); // CHANGE THIS
//				const VectorDimD N2(Nj.second->openNeighborLink(0)->glidePlaneNormal); // CHANGE THIS
//                const bool is1sessile(Ni.second->openNeighborLink(0)->sessilePlaneNormal.norm()>FLT_EPSILON);
//                const bool is2sessile(Nj.second->openNeighborLink(0)->sessilePlaneNormal.norm()>FLT_EPSILON);
//
//
//
//                if (!is1sessile && !is2sessile)
//                {
//                    VectorDimD linePoint((P1+P2)*0.5);
//                    const double denom(1.0-std::pow(N1.dot(N2),2));
//                    if(std::fabs(denom)>FLT_EPSILON)
//                    { // planes are incident: make sure to place intersection on common line
//                        const double numer((P2-P1).dot(N2));
//                        const double u(numer/denom);
//                        linePoint = P1+(N2-N2.dot(N1)*N1)*u;
//                    }
//
//                    // If mesh is used, make cure that intersection point is inside mesh
//                    bool linePointInsideMesh(true);
//                    if (DN.shared.use_boundary)
//                    {
//                        //                        linePointInsideMesh*=DN.shared.mesh.isStrictlyInsideMesh(linePoint,Ni.second->includingSimplex(),FLT_EPSILON).first;
//                        linePointInsideMesh*=DN.shared.mesh.searchWithGuess(linePoint,Ni.second->includingSimplex()).first;
//                    }
//                    if(linePointInsideMesh)
//                    {
//                        model::cout<<"Glissile Junction"<<std::endl;
//                        DN.contract(i,j,linePoint);
//                    }
//
//                }
//                else if(is1sessile && !is2sessile){ // use P1 (which is on sessile segment) as the intersection point
//                    model::cout<<"First-Sessile Junction"<<std::endl;
//                    DN.contract(i,j,P1);
//                }
//                else if(!is1sessile && is2sessile){ // use P2 (which is on sessile segment) as the intersection point
//                    model::cout<<"Second-Sessile Junction"<<std::endl;
//                    DN.contract(i,j,P2);
//                }
//                else{
//                    assert(0 && "CANNOT MAKE JUNCTION BETWEEN TWO SESSILE SEGMENTS.");
//
//                }
//
//
//
//
//			}



//			//! 2- Remove from intersectionContainer all intersections that don't satisfy Frank's rule
//			for (typename EdgeIntersectionPairContainerType::iterator iter=intersectionContainer.begin();iter!=intersectionContainer.end();)
//            {
//				const double u1(iter-> first.second);
//				const double u2(iter->second.second);
//				const VectorDimD b1(iter-> first.first->flow);
//				const VectorDimD b2(iter->second.first->flow);
//				const VectorDimD rl1(iter-> first.first->get_rl(u1));
//				const VectorDimD rl2(iter->second.first->get_rl(u2));
//
//
//                const bool is1sessile(iter-> first.first->sessilePlaneNormal.norm()>FLT_EPSILON);
//                const bool is2sessile(iter->second.first->sessilePlaneNormal.norm()>FLT_EPSILON);
//
//                VectorDimD prjDir(VectorDimD::Zero());
//                if (!is1sessile && !is2sessile){
//                    const VectorDimD commonLine(iter-> first.first->glidePlaneNormal.cross(iter->second.first->glidePlaneNormal));
//                    if(commonLine.norm()>FLT_EPSILON)
//                    { // planes are not parallel, intersection will be on common line
//                        prjDir=commonLine.normalized();
//                    }
//                    else
//                    {
//                        const double rl1Norm(rl1.norm());
//                        assert(rl1Norm>FLT_EPSILON && "TANGENT HAS ZERO NORM");
//                        prjDir=rl1/rl1Norm;
//                    }
//                }
//                else if(is1sessile && !is2sessile)
//                { // use chord of I
//                    prjDir=iter-> first.first->chord().normalized();
//                }
//                else if(!is1sessile && is2sessile)
//                { // use chord of J
//                    prjDir=iter->second.first->chord().normalized();
//                }
//                else
//                {
//                    assert(0 && "CANNOT DETERMINE COMMON LINE BETWEEN TWO SESSILE SEGMENTS.");
//                }
//
//
//                const double sgnrl1rl2(rl1.dot(prjDir)*rl2.dot(prjDir));
//
//
//				//const bool frankRule(b1.dot(b2)*rl1.dot(rl2)<=0.0);
//				const bool frankRule(b1.dot(b2)*sgnrl1rl2<=0.0);
//				const bool isValidJunction(frankRule || iter->first.first->is_boundarySegment() || iter->second.first->is_boundarySegment());
//
//
////				if (!frankRule)
//                    if (!isValidJunction)
//                {
//					iter=intersectionContainer.erase(iter);
//				}
//				else
//                { // determine the relative orientation of the links
//                    //					const int dir((rl1.dot(rl2) > 0.0) ? 1 : ((rl1.dot(rl2) < 0.0) ? -1 : 0));
//					const int dir(( sgnrl1rl2 > 0.0) ? 1 : ((sgnrl1rl2 < 0.0) ? -1 : 0));
////					dirVector.push_back(dir);
//					dirVector.emplace_back(dir);
//					++iter; // increment iterator only if didn't erase
//				}
//			}
//			model::cout<<intersectionContainer.size()<<" physical intersections. ";
//			assert(intersectionContainer.size()==dirVector.size());




//! 3- Organize intersections by segments to use multiexpand
//			typedef size_t IntersectionIDType;
//			typedef std::map<double, IntersectionIDType> MultiExpandInputType;
//			typedef std::pair<size_t,size_t> EdgeIDType;
//			typedef std::map<EdgeIDType,MultiExpandInputType> EdgeIntersectionContainerType;

//            const double dx2=std::pow(dx,2);


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
