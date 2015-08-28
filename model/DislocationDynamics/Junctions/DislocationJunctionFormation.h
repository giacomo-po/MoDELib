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
#include <model/DislocationDynamics/Remeshing/DislocationNetworkRemesh.h>
#include <model/MPI/MPIcout.h>
#include <model/Threads/EqualIteratorRange.h>

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
                
//               //std::cout<<source.get_P().transpose()<<std::endl;
//               //std::cout<<Pm.transpose()<<std::endl;
//               //std::cout<<Pp.transpose()<<std::endl;
//               //std::cout<<sink.get_P().transpose()<<std::endl;
                
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
                   && insideMeshM)
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
            
            
            //			EdgeIntersectionPairContainerType intersectionContainer;
//#ifdef _OPENMP
//            const size_t nThreads = omp_get_max_threads();
//#else
//            const size_t nThreads = 1;
//#endif
            EqualConstIteratorRange<NetworkLinkContainerType> eir(DN.linkBegin(),DN.linkEnd(),nThreads);
            
            //! 2- loop over all links and determine their intersections
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t thread=0;thread<eir.size();thread++)
            {
                for (typename NetworkLinkContainerType::const_iterator linkIterA=eir[thread].first;linkIterA!=eir[thread].second;linkIterA++)
                    //            for (typename NetworkLinkContainerType::const_iterator linkIterA=DN.linkBegin();linkIterA!=DN.linkEnd();linkIterA++)
                {
                    //                typename NetworkLinkContainerType::const_iterator linkIterA=DN.linkBegin();
                    //                std::advance(linkIterA,ll);
                    const DislocationSegmentIntersection<LinkType> dsi(linkIterA->second,linkIterA->second.glidePlaneNormal);
                    
                    for (typename NetworkLinkContainerType::const_iterator linkIterB=linkIterA;linkIterB!=DN.linkEnd();linkIterB++)
                    {
                        if (linkIterA->second.sID!=linkIterB->second.sID) // don't intersect with itself
                        {
                            
                            //std::cout<< "Intersecting "<<linkIterA->second.nodeIDPair.first<<"->"<<linkIterA->second.nodeIDPair.second<<" " <<linkIterB->second.nodeIDPair.first<<"->"<<linkIterB->second.nodeIDPair.second<<std::flush;

                            const bool& L1isSessile(linkIterA->second.isSessile);
                            const bool& L2isSessile(linkIterB->second.isSessile);

//                            const bool L1isConfined(linkIterA->second.isSessile || linkIterA->second.isOnRegionBoundary());
//                            const bool L2isConfined(linkIterB->second.isSessile || linkIterB->second.isOnRegionBoundary());
                            
                            std::set<std::pair<double,double> > temp; // the container of the roots
                            
                            if (!L1isSessile && !L2isSessile) // both are glissile
                            {
                                temp = dsi.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol);
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
                                    temp = dsi.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol);
                                }
                                else if (!gnAgnB && gnAsnB)
                                { // use planeNormal of A and sessileNormal of B
                                    temp = dsi.intersectWith(linkIterB->second,linkIterB->second.sessilePlaneNormal,collisionTol);
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
                                    temp = dsi.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol);
                                }
                                else if (!gnBgnA && gnBsnA)
                                { // use sessileNormal of A and use planeNormal of B
                                    const DislocationSegmentIntersection<LinkType> dsi2(linkIterA->second,linkIterA->second.sessilePlaneNormal);
                                    temp = dsi2.intersectWith(linkIterB->second,linkIterB->second.glidePlaneNormal,collisionTol);
                                }
                                else{
                                    assert(0 && "GLISSILE AND SESSILE PLANE NORMALS OF B MUST BE DISTINCT.");
                                }
                            }
                            else
                            { // both are sessile
                                // cannot intersect
                            }
                            
                            
                            
                            //std::cout<<" ("<<temp.size()<<" intersections):"<<std::endl;

                            const double avoidNodeIntersection=0.1;
                            
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
            
            std::deque<std::pair<size_t,std::deque<std::pair<size_t,size_t> > > > nodeDecomp;
            
            //            for (auto& nodePair : DN.nodeContainer())
            for (auto& nodePair : DN.nodes())
            {
                std::deque<std::pair<size_t,size_t> > temp=nodePair.second.edgeDecomposition();
                if(temp.size())
                {
                    nodeDecomp.emplace_back(nodePair.second.sID,temp);
                }
            }
            
            int broken=0;
            
            for(size_t n=0;n<nodeDecomp.size();++n)
            {
                const size_t& i=nodeDecomp[n].first;
                auto Ni=DN.node(i);
                assert(Ni.first);
                //                size_t m=DN.insertVertex(Ni.second->get_P());
                
                // Check that Links still exist
                bool linksexist=true;
                for(size_t d=0;d<nodeDecomp[n].second.size();++d)
                {
                    const size_t& j=nodeDecomp[n].second[d].first;
                    const size_t& k=nodeDecomp[n].second[d].second;
                    if(i==j)
                    {
                        auto Lik=DN.link(i,k);
                        linksexist*=Lik.first;
                    }
                    else if(i==k)
                    {
                        auto Lji=DN.link(j,i);
                        linksexist*=Lji.first;
                    }
                    else
                    {
                        assert(0 && "i must be equal to either j or k.");
                    }
                }
                
                if(linksexist)
                {
                    size_t m=DN.insertVertex(Ni.second->get_P()).first->first;
                    
                    for(size_t d=0;d<nodeDecomp[n].second.size();++d)
                    {
                        const size_t& j=nodeDecomp[n].second[d].first;
                        const size_t& k=nodeDecomp[n].second[d].second;
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
    
} // namespace model
#endif


//                        const NodeType& source1=*(L1.second->source);
//                        const NodeType&   sink1=*(L1.second->sink);
//                        const NodeType& source2=*(L2.second->source);
//                        const NodeType&   sink2=*(L2.second->sink);
//                                                const Simplex<dim,dim>* S1(source1.includingSimplex());
//                                                const Simplex<dim,dim>* S2(source2.includingSimplex());

//                        const int dx2=pow(dx,2);
//
//
//                        size_t im = source1.sID; // initialize first junciton point on segment1
//                        size_t ip = sink1.sID;   // initialize second junciton point on segment1
//                        if(L1.second->chordLength()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
//                        {// segment1 is not small, create new junction nodes
//
//                            const double du1(dx/L1.second->chordLength());
//                            //std::cout<<"du1="<<du1<<std::endl;
//
//                            double u1m(intersectionContainer[tt][interID]. first.second-du1);
//                            if(u1m<0.0)
//                            {
//                                u1m=0.0;
//                            }
//                            if(u1m>1.0)
//                            {
//                                u1m=1.0;
//                            }
//                            VectorDimD P1m(L1.second->get_r(u1m));
//                            P1m=L1.second->glidePlane.snapToLattice(P1m);
//
//
//
//                            double u1p(intersectionContainer[tt][interID]. first.second+du1);
//                            if(u1p<u1m)
//                            {
//                                u1p=u1m;
//                            }
//                            if(u1p>1.0)
//                            {
//                                u1p=1.0;
//                            }
//                            VectorDimD P1p(L1.second->get_r(u1p));
//                            P1p=L1.second->glidePlane.snapToLattice(P1p);
//
////                            bool insideMesh1m=true;
////                            bool insideMesh1p=true;
////                            if(shared.use_boundary)
////                            {
////                                insideMesh1m=DN.shared.mesh.searchWithGuess(P1m,S1).first;
////                                if(!insideMesh1m)
////                                {
////                                    P1m=source1.get_P()*(1.0-u1m)+sink1.get_P()*u1m;
////                                    P1m=L1.second->glidePlane.snapToLattice(P1m);
////                                    insideMesh1m=DN.shared.mesh.searchWithGuess(P1m,S1).first;
////                                }
////
////                                insideMesh1p=DN.shared.mesh.searchWithGuess(P1p,S1).first;
////                                if(!insideMesh1p)
////                                {
////                                    P1p=source1.get_P()*(1.0-u1p)+sink1.get_P()*u1p;
////                                    P1p=L1.second->glidePlane.snapToLattice(P1p);
////                                    insideMesh1p=DN.shared.mesh.searchWithGuess(P1p,S1).first;
////                                }
////                            }
//
//                            if(  (P1m-source1.get_P()).squaredNorm()>dx2)
//                            {
//                                std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(source1.sID,sink1.sID,P1m);
//                                im=temp.first->first; // id of the node obtained expanding L1
//                            }
//
//                            if((P1p-  sink1.get_P()).squaredNorm()>dx2
//                               &&          (P1m-P1p).squaredNorm()>dx2)
//                            {
//                                std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
//                                ip=temp.first->first; // id of the node obtained expanding L1
//                            }
//                        }


//                        size_t jm = source2.sID;
//                        size_t jp = sink2.sID;
//                        if(L2.second->chordLength()>DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
//                        {
//                            const double du2(dx/L2.second->chordLength());
//                            //std::cout<<"du2="<<du2<<std::endl;
//
//                            double u2m(intersectionContainer[tt][interID].second.second-du2);
//                            if(u2m<0.0)
//                            {
//                                u2m=0.0;
//                            }
//                            if(u2m>1.0)
//                            {
//                                u2m=1.0;
//                            }
//                            VectorDimD P2m(L2.second->get_r(u2m));
//                            P2m=L2.second->glidePlane.snapToLattice(P2m);
//
//                            double u2p(intersectionContainer[tt][interID].second.second+du2);
//                            if(u2p<u2m)
//                            {
//                                u2p=u2m;
//                            }
//                            if(u2p>1.0)
//                            {
//                                u2p=1.0;
//                            }
//                            VectorDimD P2p(L2.second->get_r(u2p));
//                            P2p=L2.second->glidePlane.snapToLattice(P2p);
//
//                            if(  (P2m-source2.get_P()).squaredNorm()>dx2)
//                            {
//                                std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(source2.sID,sink2.sID,P2m);
//                                jm=temp.first->first; // id of the node obtained expanding L1
//                            }
//
//                            if((P2p-  sink2.get_P()).squaredNorm()>dx2
//                               &&          (P2m-P2p).squaredNorm()>dx2)
//                            {
//                                std::pair<typename NetworkNodeContainerType::iterator,bool> temp=DN.expand(jm,sink2.sID,P2p); // now L1.second is invalid
//                                jp=temp.first->first; // id of the node obtained expanding L1
//                            }
//
//                        }

//
//                        if(du1>0.5)
//                        {
//                            du1=0.5;
//                        }
//                        if(L1.second->chordLength()<DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
//                        {
//                            du1=1.0;
//                        }
//
//                        double du2(dx/L2.second->chordLength());
//                        if(du2>0.5)
//                        {
//                            du2=0.5;
//                        }
//                        if(L2.second->chordLength()<DislocationNetworkRemesh<DislocationNetworkType>::Lmin)
//                        {
//                            du2=1.0;
//                        }



//                        const VectorDimD C1(L1.second->chord());
//                        const VectorDimD C2(L2.second->chord());
//
//                        const VectorDimD N1(L1.second->glidePlaneNormal);
//                        const VectorDimD N2(L2.second->glidePlaneNormal);
//                        const Simplex<dim,dim>* S1(source1.includingSimplex());
//                        const Simplex<dim,dim>* S2(source2.includingSimplex());



//                        if (DN.shared.use_boundary)
//                        {
//                            std::pair<bool,const Simplex<dim,dim>*> search=DN.shared.mesh.searchWithGuess(P1m,S1);
//                            if(!search.first)
//                            {
//                                bringBackToMesh(P1m,search,C1,N1);
//                            }
//
//                            search=DN.shared.mesh.searchWithGuess(P1p,S1);
//                            if(!search.first)
//                            {
//                                bringBackToMesh(P1p,search,C1,N1);
//                            }
//
//                            search=DN.shared.mesh.searchWithGuess(P2m,S2);
//                            if(!search.first)
//                            {
//                                bringBackToMesh(P2m,search,C2,N2);
//                            }
//
//                            search=DN.shared.mesh.searchWithGuess(P2p,S2);
//                            if(!search.first)
//                            {
//                                bringBackToMesh(P2p,search,C2,N2);
//                            }
//                        }


// Prepare lower point on first link
//                        size_t im = source1.sID;
//                        if (u1m > avoidNodeIntersection)


// Prepare upper point on first link
//                        size_t ip = sink1.sID;
//                        if (u1p < 1.0-avoidNodeIntersection)
//                        {
//                            const std::pair<bool,size_t> success1p=DN.expand(im,sink1.sID,P1p); // now L1.second is invalid
//                            assert(success1p.first && "COULD NOT EXPLAND LINK1 AT UPPER INTERSECTION");
//                            ip=success1p.second; // id of the node obtained expanding L1
//                        }

//                        // Prepare lower point on second link
//                        if (u2m > avoidNodeIntersection)
//                        {
//                            const std::pair<bool,size_t> success2m=DN.expand(source2.sID,sink2.sID,P2m); // now L2.second is invalid
//                            assert(success2m.first && "COULD NOT EXPLAND LINK2 AT LOWER INTERSECTION");
//                            jm=success2m.second; // id of the node obtained expanding L2
//                        }

//                        // Prepare upper point on second link
//                        if (u2p < 1.0-avoidNodeIntersection)
//                        {
//                            const std::pair<bool,size_t> success2p=DN.expand(jm,sink2.sID,P2p); // now L2.second is invalid
//                            assert(success2p.first && "COULD NOT EXPLAND LINK2 AT UPPER INTERSECTION");
//                            jp=success2p.second; // id of the node obtained expanding L2
//                        }


//        /**********************************************************************/
//        //        size_t contractWithConstraintCheck(const size_t& i, const size_t& j)
//        size_t contractWithConstraintCheck(const isNetworkNodeType& N1, const isNetworkNodeType& N2)
//        {
//            size_t contracted(0);
//
//            const VectorDimD P1=N1.second->get_P();
//            const VectorDimD P2=N2.second->get_P();
//            const VectorDimD P12=P1-P2;
//            const double P12norm(P12.norm());
//
//            if(P12norm<FLT_EPSILON) // points are coincindent, just contract them
//            {
//                contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
//            }
//            else // points are distinct
//            {
//
////                const VectorDimD unitP12=P12/P12norm;
//
//
//                const typename DislocationNetworkType::NodeType::VectorOfNormalsType PN1(N1.second->constraintNormals());
//                const typename DislocationNetworkType::NodeType::VectorOfNormalsType PN2(N2.second->constraintNormals());
//                const size_t sizePN1(PN1.size());
//                const size_t sizePN2(PN2.size());
//
//                //std::cout<<"contractWithConstraintCheck "<<N1.second->sID<<" "<<N2.second->sID<<std::endl;
//                //std::cout<<"contractWithConstraintCheck: case "<<sizePN1<<" "<<sizePN2<<std::endl;
//
//                if(sizePN1==1 && sizePN2==1) // nodes constrained to move on planes
//                {
//                    //std::cout<<"contractWithConstraintCheck, case 1"<<std::endl;
////                    const double denom(1.0-std::pow(PN1[0].dot(PN2[0]),2));
////                    const double numer((P2-P1).dot(PN2[0]));
////
////                    if(denom<FLT_EPSILON) // parallel or coincident planes
////                    {
////                        if(std::fabs(denom)<FLT_EPSILON) // planes are coincident
////                        {
////                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
////                        }
////                        else // parallel planes
////                        {
////                            assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
////                        }
////                    }
////                    else // incident planes
////                    {
////                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P1+(PN2[0]-PN2[0].dot(PN1[0])*PN1[0])*numer/denom);
////                    }
//
////                    const VectorDimD d3(PN1[0].cross(PN2[0]));
////                    const double d3norm(d3.norm());
////                    if (d3norm<FLT_EPSILON) // parallel or coincident planes
//                    if ((PN1[0].cross(PN2[0])).norm()<FLT_EPSILON) // parallel or coincident planes
//                    {
//                        if(std::fabs((P1-P2).dot(PN1[0]))<FLT_EPSILON) // planes are coincident
//                        {
//                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
//                        }
//                        else // parallel planes
//                        {
//                            assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
//                        }
//                    }
//                    else // incident planes
//                    {
//                        // Contract at X, where X minimizes 0.5*(X-P1)^2+0.5*(X-P2)^2
//                        // under the constraints (X-P1)*N1=0 and (X-P2)*N2=0
//                        const VectorDimD P(P1+P2);
//                        const Eigen::Matrix<double,3,2> N((Eigen::Matrix<double,2,3>()<<PN1[0].transpose(),PN2[0].transpose()).finished().transpose());
//                        const Eigen::Matrix<double,2,1> A((Eigen::Matrix<double,2,1>()<<P1.dot(PN1[0]),P2.dot(PN2[0])).finished());
//                        const Eigen::Matrix<double,2,1> L=(N.transpose()*N).inverse()*(N.transpose()*P-2.0*A);
//                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P-N*L));
//                    }
//
//                }
//                else if(sizePN1==1 && sizePN2==2) // N1 moves on a plane, N2 moves on a line
//                {
//                    //std::cout<<"contractWithConstraintCheck, case 2"<<std::endl;
//                    const VectorDimD d2(PN2[0].cross(PN2[1]));
//                    const double den(d2.dot(PN1[0]));
//                    const double num(P12.dot(PN1[0]));
//                    if(std::fabs(den)>FLT_EPSILON) // line and plane are not parallel
//                    {
//                        //std::cout<<"contractWithConstraintCheck, case 2a"<<std::endl;
//                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P2+num/den*d2);
//                    }
//                    else
//                    {
//                        if(std::fabs(num)<FLT_EPSILON) // P2 belongs to the plane of P1
//                        {
//                            //std::cout<<"contractWithConstraintCheck, case 2b"<<std::endl;
//                            //                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2)); // HERE
//                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
//
//                        }
//                        //std::cout<<"contractWithConstraintCheck, case 2c"<<std::endl;
//                    }
//                }
//                else if(sizePN1==2 && sizePN2==1) // N1 moves on a line, N2 moves on a plane
//                {
//                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
//                    //                    //std::cout<<"contractWithConstraintCheck, case 3"<<std::endl;
//                    //                    const VectorDimD d1(PN1[0].cross(PN1[1]));
//                    //                    const double den(d1.dot(PN2[0]));
//                    //                    const double num((P2-P1).dot(PN2[0]));
//                    //                    if(std::fabs(den)>FLT_EPSILON) // line and plane are not parallel
//                    //                    {
//                    //                        //std::cout<<"contractWithConstraintCheck, case 3a"<<std::endl;
//                    //                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P1+num/den*d1);
//                    //                    }
//                    //                    else
//                    //                    {
//                    //                        if(std::fabs(num)<FLT_EPSILON) // P1 belongs to the plane of P2
//                    //                        {
//                    //                            //std::cout<<"contractWithConstraintCheck, case 3b"<<std::endl;
//                    //
//                    ////                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2)); // HERE
//                    //                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
//                    //                        }
//                    //                    }
//                }
//                else if(sizePN1==2 && sizePN2==2) // both N1 and N2 move on lines
//                {
//                    //std::cout<<"contractWithConstraintCheck, case 4"<<std::endl;
//
//                    VectorDimD d1(PN1[0].cross(PN1[1]));
//                    const double d1norm(d1.norm());
//                    assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
//                    d1/=d1norm;
//                    VectorDimD d2(PN2[0].cross(PN2[1]));
//                    const double d2norm(d2.norm());
//                    assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
//                    d2/=d2norm;
//
//                    const VectorDimD d3(d1.cross(d2));
//                    const double d3Norm(d3.norm());
//
//                    if(d3Norm<FLT_EPSILON) // d1 and d2 are aligned, this means colinear or no intersection
//                    {
//                        if(d1.cross(P12.normalized()).norm()<FLT_EPSILON) // colinear
//                        {
//                            //std::cout<<"contractWithConstraintCheck, case 4a"<<std::endl;
//
//                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
//                        }
//                        else // parallel (not colinear) lines, contraction not possible
//                        {
//                            //assert(0 && "COULD NOT CONTRACT JUNCTION POINTS. ALIGNMENT CONDITION FAILED.");
//
//                        }
//                    }
//                    else // d1 and d2 are not aligned
//                    {
//                        if(std::fabs((d3/d3Norm).dot(P12/P12norm))<FLT_EPSILON) // planarity condition
//                        {
//                            //std::cout<<"contractWithConstraintCheck, case 4b"<<std::endl;
//
//                            const VectorDimD dOrth=d2-d2.dot(d1)*d1; // component of d2 orthogonal to d1
//                            const double den=d2.dot(dOrth);
//                            assert(std::fabs(den)>FLT_EPSILON && "YOU SHOULD HAVE FOUND THIS ABOVE.");
//                            contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractWithCommonNeighborCheck(*N1.second,*N2.second,P2+P12.dot(dOrth)/den*d2);
//                        }
//                        else // planarity condition failed, contraction not possible
//                        {
//                            //assert(0 && "COULD NOT CONTRACT JUNCTION POINTS. PLANARITY CONDITION FAILED.");
//                        }
//                    }
//                }
//                else if(sizePN1==1 && sizePN2==3)
//                {
//                    //std::cout<<"contractWithConstraintCheck, case 5"<<std::endl;
//
//                    if(std::fabs(P12.normalized().dot(PN1[0]))<FLT_EPSILON) // P2 belongs to the plane of P1
//                    {// contract N1
//                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
//                    }
//                    else // contraction not possible
//                    {
//
//                    }
//                }
//                else if(sizePN1==3 && sizePN2==1)
//                {
//                    //std::cout<<"contractWithConstraintCheck, case 6"<<std::endl;
//
//                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
//
//                    //                    if(std::fabs(P12.normalized().dot(PN2[0]))<FLT_EPSILON) // P1 belongs to the plane of P2
//                    //                    {// contract N2
//                    //                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
//                    //                    }
//                    //                    else // contraction not possible
//                    //                    {
//                    //
//                    //                    }
//                }
//                else if(sizePN1==2 && sizePN2==3)
//                {
//                    //std::cout<<"contractWithConstraintCheck, case 7"<<std::endl;
//
//                    VectorDimD d1(PN1[0].cross(PN1[1]));
//                    const double d1norm(d1.norm());
//                    assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
//                    d1/=d1norm;
//                    if(P12.normalized().cross(d1).norm()<FLT_EPSILON) // P2 belongs to the line of P1
//                    {
//                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
//                    }
//                    else // contraction not possible
//                    {
//
//                    }
//                }
//                else if(sizePN1==3 && sizePN2==2)
//                {
//                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
//
//                    //                    VectorDimD d2(PN2[0].cross(PN2[1]));
//                    //                    const double d2norm(d2.norm());
//                    //                    assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
//                    //                    d2/=d2norm;
//                    //                    if(P12.normalized().cross(d2).norm()<FLT_EPSILON) // P2 belongs to the line of P1
//                    //                    {
//                    //                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
//                    //                    }
//                    //                    else // contraction not possible
//                    //                    {
//                    //
//                    //                    }
//                }
//                else if(sizePN1==3 && sizePN2==3)
//                {
//                    if(P12norm<FLT_EPSILON)
//                    {
//                        contracted+=DislocationNetworkRemesh<DislocationNetworkType>(DN).contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
//                    }
//                }
//            }
//
//            return contracted;
//        }


