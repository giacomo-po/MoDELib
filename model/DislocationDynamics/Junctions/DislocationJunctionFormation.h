/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONJUNCTIONFORMATION_H_
#define model_DISLOCATIONJUNCTIONFORMATION_H_

#include <iostream>
#include <utility> // for std::pair
#include <vector>
#include <Eigen/Dense>
#include <model/Network/Operations/EdgeFinder.h>
#include <model/DislocationDynamics/Junctions/DislocationSegmentIntersection.h>
#include <model/DislocationDynamics/Remeshing/DislocationNetworkRemesh.h>
#include <model/MPI/MPIcout.h>
#include <model/Threads/EqualIteratorRange.h>
#include <model/Threads/N2IteratorRange.h>
#define VerboseJunctions(N,x) if(verboseJunctions>=N){model::cout<<x;}

namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationJunctionFormation
    {
        
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
        
        typedef LatticeVector<dim> LatticeVectorType;
        
        //        /**********************************************************************/
        //        void bringBackToPlane(VectorDimD& P, const LinkType& L)
        //        {
        //            const VectorDimD& n(L.glidePlaneNormal);
        //            P -= (P-0.5*(L.source->get_P()+L.sink->get_P())).dot(n)*n;
        //        }
        
        /**********************************************************************/
        void bringBackToMesh(VectorDimD& P,
                             const std::pair<bool,const Simplex<dim,dim>*>& search,
                             const VectorDimD& c,
                             const VectorDimD& g) const
        {
            
            const VectorDimD dir(c.cross(g));
            assert(dir.squaredNorm()>FLT_EPSILON && "bringBackToMesh: direction has zero norm.");
            const VectorDimD d(dir.normalized());
            
            int faceID;
            search.second->pos2bary(P).minCoeff(&faceID); // find the ID of the face with minimum bary coordinate
            const bool isBoundaryFace(search.second->child(faceID).isBoundarySimplex());
            assert(isBoundaryFace && "bringBackToMesh: face is not a boundary face.");
            
            const VectorDimD V(search.second->child(faceID).vertices()[0]->P0);
            const VectorDimD N(search.second->nda.col(faceID).normalized());
            const double den(d.dot(N));
            if(std::fabs(den)<FLT_EPSILON)
            {
                model::cout<<"d="<<d.transpose()<<std::endl;
                model::cout<<"N="<<N.transpose()<<std::endl;
                assert(0 && "bringBackToMesh: direction is parallel to mesh face.");
            }
            
            const double u((V-P).dot(N)/den);
            
            P += u*d;
        }
        
        /**********************************************************************/
        //		int junctionDir(const EdgeIntersectionPairType& intersectionPair) const
        int junctionDir(const LinkType& L1, const LinkType& L2,
                        const double& u1, const double& u2) const
        {
            //            const double u1(intersectionPair.first.second);
            //            const double u2(intersectionPair.second.second);
            const VectorDimD b1(L1.Burgers);
            const VectorDimD b2(L2.Burgers);
            const VectorDimD rl1(L1.get_rl(u1));
            const VectorDimD rl2(L2.get_rl(u2));
            
            
            //            const bool L1.isSessile(L1.sessilePlaneNormal.norm()>FLT_EPSILON);
            //            const bool L2.isSessile(L2.sessilePlaneNormal.norm()>FLT_EPSILON);
            
            //            const bool L1.isSessile(L1.flow.dot()>FLT_EPSILON);
            //            const bool L2.isSessile(L2.sessilePlaneNormal.norm()>FLT_EPSILON);
            
            
            VectorDimD prjDir(VectorDimD::Zero());
            if (!L1.isSessile && !L2.isSessile)
            {
                //std::cout<<"junctionDir, case 1"<<std::endl;
                const VectorDimD commonLine(L1.glidePlaneNormal.cross(L2.glidePlaneNormal));
                if(commonLine.norm()>FLT_EPSILON)
                { // planes are not parallel, intersection will be on common line
                    //std::cout<<"junctionDir, case 1a"<<std::endl;
                    prjDir=commonLine.normalized();
                }
                else
                {
                    //std::cout<<"junctionDir, case 1b"<<std::endl;
                    const double rl1Norm(rl1.norm());
                    assert(rl1Norm>FLT_EPSILON && "TANGENT HAS ZERO NORM");
                    prjDir=rl1/rl1Norm;
                }
            }
            else if(L1.isSessile && !L2.isSessile)
            { // use chord of I
                //std::cout<<"junctionDir, case 2"<<std::endl;
                prjDir=L1.chord().normalized();
            }
            else if(!L1.isSessile && L2.isSessile)
            { // use chord of J
                //std::cout<<"junctionDir, case 2"<<std::endl;
                prjDir=L2.chord().normalized();
            }
            else
            {
                assert(0 && "CANNOT DETERMINE COMMON LINE BETWEEN TWO SESSILE SEGMENTS.");
            }
            
            
            const double sgnrl1rl2(rl1.dot(prjDir)*rl2.dot(prjDir));
            
            //std::cout<<"junctionDir, sgnrl1rl2="<<sgnrl1rl2<<std::endl;
            
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
        std::pair<size_t,size_t> junctionIDs(const LinkType& L,const double& u,const double& dx)
        {
            const int dx2=pow(dx,2);
            
            const NodeType& source=*(L.source);
            const NodeType&   sink=*(L.sink);
            
            
            size_t im = source.sID; // initialize first junciton point on segment1
            size_t ip = sink.sID;   // initialize second junciton point on segment1
            
            if(L.chordLength()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
            {// segment1 is not small, create new junction nodes
                
                const double du(dx/L.chordLength());
                
                double um(u-du);
                if(um<0.0)
                {
                    um=0.0;
                }
                if(um>1.0)
                {
                    um=1.0;
                }
                VectorDimD Pm(L.get_r(um));
                Pm=L.glidePlane.snapToLattice(Pm);
                
                double up(u+du);
                if(up<um)
                {
                    up=um;
                }
                if(up>1.0)
                {
                    up=1.0;
                }
                VectorDimD Pp(L.get_r(up));
                Pp=L.glidePlane.snapToLattice(Pp);
                
                bool insideMeshM=true;
                bool insideMeshP=true;
                if(DN.shared.use_boundary)
                {
                    const Simplex<dim,dim>* S(source.includingSimplex());
                    
                    insideMeshM=DN.shared.mesh.searchWithGuess(Pm,S).first;
                    if(!insideMeshM)
                    {
                        Pm=source.get_P()*(1.0-um)+sink.get_P()*um;
                        Pm=L.glidePlane.snapToLattice(Pm);
                        insideMeshM=DN.shared.mesh.searchWithGuess(Pm,S).first;
                    }
                    
                    insideMeshP=DN.shared.mesh.searchWithGuess(Pp,S).first;
                    if(!insideMeshP)
                    {
                        Pp=source.get_P()*(1.0-up)+sink.get_P()*up;
                        Pp=L.glidePlane.snapToLattice(Pp);
                        insideMeshP=DN.shared.mesh.searchWithGuess(Pp,S).first;
                    }
                }
                
                //                if(   (Pm-source.get_P()).squaredNorm()>dx2
                //                   && insideMeshM)
                //                {
                //                    std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(source.sID,sink.sID,Pm);
                //                    im=temp.first->first; // id of the node obtained expanding L1
                //                }
                //
                //                if(   (Pp-  sink.get_P()).squaredNorm()>dx2
                //                   && (Pm-Pp).squaredNorm()>dx2
                //                   && insideMeshP)
                //                {
                //                    std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(im,sink.sID,Pp); // now L1.second is invalid
                //                    ip=temp.first->first; // id of the node obtained expanding L1
                //                }
                
                if(   (Pm-source.get_P()).squaredNorm()>dx2
                   && (Pm-Pp).squaredNorm()>dx2
                   && (Pp-  sink.get_P()).squaredNorm()>dx2
                   && insideMeshM
                   && insideMeshP)
                {
                    //std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(source.sID,sink.sID,Pm);
                    im=DN.expand(source.sID,sink.sID,Pm).first->first; // id of the node obtained expanding L1
                    
                    //std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(im,sink.sID,Pp); // now L1.second is invalid
                    ip=DN.expand(im,sink.sID,Pp).first->first; // id of the node obtained expanding L1
                    
                }
                
            }
            
            
            return std::make_pair(im,ip);
        }
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        
    public:
        
        //! The tolerance (in units of distance) used for collision detection
        static double collisionTol;
        //        static bool useVertexEdgeJunctions;
        
        static int verboseJunctions;
        
        /* Constructor ********************************************************/
        DislocationJunctionFormation(DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in)
        {
            
        }
        
        /* findIntersections **************************************************/
        //		EdgeIntersectionPairContainerType findIntersections(const double& avoidNodeIntersection) const
        void findIntersections(std::deque<EdgeIntersectionPairContainerType>& intersectionContainer,
                               std::deque<std::deque<int>>& dirVector,
                               const size_t& nThreads) const
        
        {/*! @param[in]  avoidNodeIntersection
          *  Computes all the intersections between the edges of the DislocationNetwork
          */
            
            //                                        const double avoidNodeIntersection=0.2;
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Finding Junctions ("<<nThreads<<" threads)... "<<std::flush;
            
            // Create an EqualConstIteratorRange over links
            //            EqualConstIteratorRange<NetworkLinkContainerType> eir(DN.linkBegin(),DN.linkEnd(),nThreads);
            N2IteratorRange<typename NetworkLinkContainerType::const_iterator> eir(DN.linkBegin(),DN.linkEnd(),nThreads);
            
            
            //            std::cout<<"#links="<<DN.linkOrder()<<std::endl;
            //            for(auto pair : eir)
            //            {
            //                std::cout<<std::distance(pair.first,pair.second)<<std::endl;
            //            }
            //
            //            std::vector<int> threadVector(nThreads,0);
            
            
            const double avoidNodeIntersection=0.1;
            
            
            //! 2- loop over all links and determine their intersections
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t thread=0;thread<eir.size();thread++)
            {
                for (typename NetworkLinkContainerType::const_iterator linkIterA=eir[thread].first;linkIterA!=eir[thread].second;linkIterA++)
                {
                    const DislocationSegmentIntersection<LinkType> dsi(linkIterA->second,linkIterA->second.glidePlaneNormal);
                    
                    for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.linkEnd();linkIterB++)
                    {
                        if (linkIterA->second.sID!=linkIterB->second.sID) // don't intersect with itself
                        {
                            //                            threadVector[omp_get_thread_num()]++;
                            
                            VerboseJunctions(1,"Intersecting "<<linkIterA->second.nodeIDPair.first<<"->"<<linkIterA->second.nodeIDPair.second<<" " <<linkIterB->second.nodeIDPair.first<<"->"<<linkIterB->second.nodeIDPair.second<<std::flush)
                            //std::cout<< "Intersecting "<<linkIterA->second.nodeIDPair.first<<"->"<<linkIterA->second.nodeIDPair.second<<" " <<linkIterB->second.nodeIDPair.first<<"->"<<linkIterB->second.nodeIDPair.second<<std::flush;
                            
                            const bool& L1isSessile(linkIterA->second.isSessile);
                            const bool& L2isSessile(linkIterB->second.isSessile);
                            
                            std::set<std::pair<double,double> > temp; // the container of the roots
                            bool checkJunction(true);
                            if (linkIterA->second.is_boundarySegment() && linkIterB->second.is_boundarySegment())
                            {
								if ((linkIterA->second.glidePlaneNormal.cross(linkIterB->second.glidePlaneNormal)).norm()>FLT_EPSILON)
								{
									 checkJunction=false; // don't check junction between bnd segments on different planes
								}
								else
								{
									if((linkIterA->second.source->get_P()-linkIterB->second.source->get_P()).dot(linkIterA->second.glidePlaneNormal)>FLT_EPSILON)
									{
										checkJunction=false; // don't check junction between bnd segments on parallel offset planes
									}
								}
			              }
			              if (checkJunction)
			              {

                            if (!L1isSessile && !L2isSessile) // both are glissile
                            {
                                temp = dsi.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol,avoidNodeIntersection);
                            }
                            else if (!L1isSessile && L2isSessile) // L1 is glissile and L2 is sessile
                            {
                                const bool gnAgnB((linkIterA->second.glidePlaneNormal-linkIterB->second.glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                                const bool gnAsnB((linkIterA->second.glidePlaneNormal-linkIterB->second.sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                                
                                if (!gnAgnB && !gnAsnB)
                                {
                                    // cannot intersect
                                }
                                else if(gnAgnB && !gnAsnB)
                                { // use planeNormal of A and planeNormal of B
                                    temp = dsi.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol,avoidNodeIntersection);
                                }
                                else if (!gnAgnB && gnAsnB)
                                { // use planeNormal of A and sessileNormal of B
                                    temp = dsi.intersectWith(linkIterB->second,linkIterB->second.sessilePlaneNormal,collisionTol,avoidNodeIntersection);
                                }
                                else{
                                    assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                                }
                            }
                            else if (L1isSessile && !L2isSessile) // L1 is sessile and L2 is glissile
                            {
                                const bool gnBgnA((linkIterB->second.glidePlaneNormal-linkIterA->second.glidePlaneNormal  ).squaredNorm()<FLT_EPSILON);
                                const bool gnBsnA((linkIterB->second.glidePlaneNormal-linkIterA->second.sessilePlaneNormal).squaredNorm()<FLT_EPSILON);
                                
                                if (!gnBgnA && !gnBsnA)
                                {
                                    // cannot intersect
                                }
                                else if(gnBgnA && !gnBsnA)
                                { // use planeNormal of A and planeNormal of B
                                    temp = dsi.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol,avoidNodeIntersection);
                                }
                                else if (!gnBgnA && gnBsnA)
                                { // use sessileNormal of A and use planeNormal of B
                                    const DislocationSegmentIntersection<LinkType> dsi2(linkIterA->second,linkIterA->second.sessilePlaneNormal);
                                    temp = dsi2.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol,avoidNodeIntersection);
                                }
                                else{
                                    assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                                }
                            }
                            else
                            { // both are sessile
                                // cannot intersect
                            }
                            
                         }    
                            for (std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paramIter!=temp.end();++paramIter)
                            {
                                //                            if (   paramIter->first >avoidNodeIntersection && paramIter-> first<1.0-avoidNodeIntersection // THIS IS CHECKED LATER
                                //                                && paramIter->second>avoidNodeIntersection && paramIter->second<1.0-avoidNodeIntersection) // avoid node intersection
                                //                            {
                                
                                const bool intersectionIsSourceSource(   paramIter->first  <avoidNodeIntersection
                                                                      && paramIter->second <avoidNodeIntersection
                                                                      && linkIterA->second.source->sID==linkIterB->second.source->sID);
                                
                                const bool   intersectionIsSourceSink(   paramIter->first  <avoidNodeIntersection
                                                                      && paramIter->second >1.0-avoidNodeIntersection
                                                                      && linkIterA->second.source->sID==linkIterB->second.sink->sID);
                                
                                const bool   intersectionIsSinkSource(   paramIter->first  > 1.0-avoidNodeIntersection
                                                                      && paramIter->second <avoidNodeIntersection
                                                                      && linkIterA->second.sink->sID==linkIterB->second.source->sID);
                                
                                const bool     intersectionIsSinkSink(   paramIter->first  > 1.0-avoidNodeIntersection
                                                                      && paramIter->second > 1.0-avoidNodeIntersection
                                                                      && linkIterA->second.sink->sID==linkIterB->second.sink->sID);
                                
                                
                                if(!intersectionIsSourceSource && !intersectionIsSourceSink && !intersectionIsSinkSource && !intersectionIsSinkSink)
                                {
                                    EdgeIntersectionType intersectionOnA(std::make_pair(linkIterA->second.nodeIDPair,paramIter->first ));
                                    EdgeIntersectionType intersectionOnB(std::make_pair(linkIterB->second.nodeIDPair,paramIter->second));
                                    const int dir(junctionDir(linkIterA->second,linkIterB->second,paramIter->first,paramIter->second));
                                    
                                    //                               //std::cout<<paramIter->first<<" "<<paramIter->second<<" "<<dir<<std::endl;
                                    
                                    if(dir!=0)
                                    {
#ifdef _OPENMP
                                        intersectionContainer[omp_get_thread_num()].emplace_back(intersectionOnA,intersectionOnB);
                                        dirVector[omp_get_thread_num()].emplace_back(dir);
#else
                                        intersectionContainer[0].emplace_back(intersectionOnA,intersectionOnB);
                                        dirVector[0].emplace_back(dir);
#endif
                                    }
                                }
                            } // end for
                        } // end don't self-intersect
                    } // end loop over second segment
                } // end loop over first segment
            }// end loop ever threads
            
            
            //            for(int k=0;k<threadVector.size();++k)
            //            {
            //                std::cout<<"thread "<<k<<" computed "<<threadVector[k]<<" intersections"<<std::endl;
            //            }
            
            int nIntersections=0;
            for (size_t tt=0;tt<intersectionContainer.size();++tt)
            {
                assert(intersectionContainer[tt].size()==dirVector[tt].size());
                nIntersections+=intersectionContainer[tt].size();
            }
            model::cout<<nIntersections<<" physical intersections. ";
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        void formJunctions(const double& dx)
        {
            
            
            //! 1- Initialize intersectionContainer calling findIntersections deque<Pair<Pair<Link*double>,Pair<Link*,double>>>
            std::deque<EdgeIntersectionPairContainerType> intersectionContainer;
            std::deque<std::deque<int>> dirVector;
            
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            intersectionContainer.resize(nThreads);
            dirVector.resize(nThreads);
            findIntersections(intersectionContainer,dirVector,nThreads);
            
            
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Forming Junctions: "<<std::flush;
            
            typedef std::pair<size_t,size_t> EdgeIDType;
            
            
            for (size_t tt=0;tt<intersectionContainer.size();++tt)
            {
                for (size_t interID=0;interID!=intersectionContainer[tt].size();++interID)
                {
                    const EdgeIDType& key1(intersectionContainer[tt][interID]. first.first);
                    const EdgeIDType& key2(intersectionContainer[tt][interID].second.first);
                    
                    const isNetworkLinkType L1(DN.link(key1.first,key1.second));
                    const isNetworkLinkType L2(DN.link(key2.first,key2.second));
                    
                    //std::cout<<"forming Junction "<< key1.first<<"->"<<key1.second<<" and "<< key2.first<<"->"<<key2.second<<" @"<<intersectionContainer[tt][interID]. first.second<<","<<intersectionContainer[tt][interID]. second.second<<std::endl;
                    
                    if(L1.first && L2.first) // Links exist
                    {
                        //std::cout<<"I'm here 1"<<std::endl;
                        
                        const std::pair<size_t,size_t> I=junctionIDs(*L1.second,intersectionContainer[tt][interID]. first.second,dx);
                        const size_t im=I.first;
                        const size_t ip=I.second;
                        
                        
                        const std::pair<size_t,size_t> J=junctionIDs(*L2.second,intersectionContainer[tt][interID].second.second,dx);
                        const size_t jm=J.first;
                        const size_t jp=J.second;
                        
                        
                        switch (dirVector[tt][interID])
                        {
                            case +1:
                            {
                                //std::cout<<"+1: im="<<im<<", jm="<<jm<<std::endl;
                                //std::cout<<"+1: ip="<<ip<<", jp="<<jp<<std::endl;
                                if(im!=jm)
                                {
                                    const isNetworkNodeType N1=DN.node(im);
                                    const isNetworkNodeType N2=DN.node(jm);
                                    if(N1.first && N2.first)
                                    {
                                        //std::cout<<"first contract +1 "<<std::endl;
                                        DN.contractWithConstraintCheck(N1,N2);
                                    }
                                }
                                if(ip!=jp)
                                {
                                    const isNetworkNodeType N1=DN.node(ip);
                                    const isNetworkNodeType N2=DN.node(jp);
                                    if(N1.first && N2.first)
                                    {
                                        //std::cout<<"second contract +1 "<<std::endl;
                                        DN.contractWithConstraintCheck(N1,N2);
                                    }
                                }
                                break;
                            }
                                
                            case -1:
                            {
                                //std::cout<<"-1: im="<<im<<", jp="<<jp<<std::endl;
                                //std::cout<<"-1: ip="<<ip<<", jm="<<jm<<std::endl;
                                if(im!=jp)
                                {
                                    const isNetworkNodeType N1=DN.node(im);
                                    const isNetworkNodeType N2=DN.node(jp);
                                    if(N1.first && N2.first)
                                    {
                                        //std::cout<<"first contract -1 "<<std::endl;
                                        DN.contractWithConstraintCheck(N1,N2);
                                    }
                                }
                                if(ip!=jm)
                                {
                                    const isNetworkNodeType N1=DN.node(ip);
                                    const isNetworkNodeType N2=DN.node(jm);
                                    if(N1.first && N2.first)
                                    {
                                        //std::cout<<"second contract -1 "<<std::endl;
                                        DN.contractWithConstraintCheck(N1,N2);
                                    }
                                }
                                break;
                            }
                                
                            default:
                                assert(0 && "DIR VECTOR CAN ONLY BE +1 or -1");
                                break;
                        }
                        
                        
                    }
                    //std::cout<<"done forming Junction "<<std::endl;
                    
                }
            } // loop over threads
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            
        }
        
        /**********************************************************************/
        void breakZeroLengthJunctions()
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Breaking zero-length Junctions... "<<std::flush;
            
            std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecompFirst;
            std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecompSecond;
            
            
            // TO DO: PARALLELIZE THIS LOOP
            for (auto& nodePair : DN.nodes())
            {
                const auto temp=nodePair.second.edgeDecomposition();
                //std::deque<std::pair<size_t,size_t> > temp=nodePair.second.edgeDecomposition();
                
                if(temp.first.size() && temp.second.size())
                {
                    nodeDecompFirst.emplace_back(nodePair.second.sID,temp.first);
                    nodeDecompSecond.emplace_back(nodePair.second.sID,temp.second);
                }
            }
            
            int broken=0;
            
            for(size_t n=0;n<nodeDecompFirst.size();++n)
            {
                const size_t& i=nodeDecompFirst[n].first;
                auto Ni=DN.node(i);
                assert(Ni.first);
                //                size_t m=DN.insertVertex(Ni.second->get_P());
                
                // Check that Links still exist
                //std::deque<VectorDimD> linkDirs;
                
                bool linksFirstexist=true;
                VectorDimD avrFirst=VectorDimD::Zero();
                for(size_t d=0;d<nodeDecompFirst[n].second.size();++d)
                {
                    const size_t& j=nodeDecompFirst[n].second[d].first;
                    const size_t& k=nodeDecompFirst[n].second[d].second;
                    if(i==j)
                    {
                        auto Lik=DN.link(i,k);
                        if(Lik.first)
                        {
                            avrFirst+=Lik.second->chord().normalized();
                        }
                        linksFirstexist*=Lik.first;
                    }
                    else if(i==k)
                    {
                        auto Lji=DN.link(j,i);
                        if(Lji.first)
                        {
                            avrFirst-=Lji.second->chord().normalized();
                        }
                        linksFirstexist*=Lji.first;
                    }
                    else
                    {
                        assert(0 && "i must be equal to either j or k.");
                    }
                }
                const double avrFirstNorm=avrFirst.norm();
                if(avrFirstNorm>FLT_EPSILON)
                {
                    avrFirst/=avrFirstNorm;
                }
                
                bool linksSecondexist=true;
                VectorDimD avrSecond=VectorDimD::Zero();
                for(size_t d=0;d<nodeDecompSecond[n].second.size();++d)
                {
                    const size_t& j=nodeDecompSecond[n].second[d].first;
                    const size_t& k=nodeDecompSecond[n].second[d].second;
                    if(i==j)
                    {
                        auto Lik=DN.link(i,k);
                        if(Lik.first)
                        {
                            avrSecond+=Lik.second->chord().normalized();
                        }
                        linksSecondexist*=Lik.first;
                    }
                    else if(i==k)
                    {
                        auto Lji=DN.link(j,i);
                        if(Lji.first)
                        {
                            avrSecond-=Lji.second->chord().normalized();
                        }
                        linksSecondexist*=Lji.first;
                    }
                    else
                    {
                        assert(0 && "i must be equal to either j or k.");
                    }
                }
                const double avrSecondNorm=avrSecond.norm();
                if(avrSecondNorm>FLT_EPSILON)
                {
                    avrSecond/=avrSecondNorm;
                }
                
                
                if(linksFirstexist && linksSecondexist && avrSecond.dot(avrFirst)<-0.0)
                {
                    std::cout<<"NodeBreaking "<<Ni.second->sID<<" "<<avrSecond.dot(avrFirst)<<std::endl;
                    
                    size_t m=DN.insertVertex(Ni.second->get_P()).first->first;
                    
                    for(size_t d=0;d<nodeDecompFirst[n].second.size();++d)
                    {
                        const size_t& j=nodeDecompFirst[n].second[d].first;
                        const size_t& k=nodeDecompFirst[n].second[d].second;
                        if(i==j)
                        {
                            auto Lik=DN.link(i,k);
                            assert(Lik.first);
                            DN.connect(m,k,Lik.second->Burgers);
                            DN.template disconnect<0>(i,k);
                        }
                        else if(i==k)
                        {
                            auto Lji=DN.link(j,i);
                            assert(Lji.first);
                            DN.connect(j,m,Lji.second->Burgers);
                            DN.template disconnect<0>(j,i);
                        }
                        else
                        {
                            assert(0 && "i must be equal to either j or k.");
                        }
                    }
                    
                    
                    broken++;
                }
            }
            model::cout<<broken<<" broken."<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
    };
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;
    
    // Declare Static Data
    template <typename DislocationNetworkType>
    int DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions=0;
    
    
} // namespace model
#endif

