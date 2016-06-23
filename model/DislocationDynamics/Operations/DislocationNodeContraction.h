/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNODECONTRACTION_H_
#define model_DISLOCATIONNODECONTRACTION_H_

#include <Eigen/Dense>
#include <model/Network/Operations/EdgeFinder.h>
#include <model/Network/Operations/VertexContraction.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/LatticeMath/LineMeshIntersection.h>
#include <model/DislocationDynamics/SimplexBndNormal.h>

namespace model
{
    
    /*! \brief Class template that handles the nodal remesh of the DislocationNetwork.
     */
    template <typename DislocationNetworkType>
    class DislocationNodeContraction
    // RENEABLE THIS LATER    : public VertexContraction<typename DislocationNetworkType::NodeType,typename DislocationNetworkType::LinkType>
    {
        
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef VertexContraction<NodeType,LinkType> VertexContractionType;
        
        typedef typename EdgeFinder<LinkType>::isNetworkEdgeType isNetworkLinkType;
        typedef typename DislocationNetworkType::isNetworkNodeType isNetworkNodeType;
        typedef typename DislocationNetworkType::VectorDimD VectorDimD;
        
        static constexpr int dim=DislocationNetworkType::dim;
        typedef LatticeVector<dim> LatticeVectorType;
        
        
        //        static double neighborRadius;
        
        DislocationNetworkType& DN;
        
        //        /**********************************************************************/
        //        std::pair<bool,const Simplex<dim,dim>*> pointIsInsideMesh(const VectorDimD& P0, const Simplex<dim,dim>* const guess)
        //        {/*!\param[in] P0 position vector
        //          * \param[in] guess pointer of the Simplex where the search starts
        //          * \returns true if P0 is inside the mesh
        //          */
        //            std::pair<bool,const Simplex<dim,dim>*> temp(true,NULL);
        //            if (DN.shared.use_boundary)
        //            {
        //                temp=DN.shared.mesh.searchWithGuess(P0,guess);
        //            }
        //            return temp;
        //        }
        
//        /**********************************************************************/
//        unsigned int contractWithCommonNeighborCheck(const NodeType& Ni,
//                                                     const NodeType& Nj,
//                                                     const LatticeVectorType& L0)
//        {/*!@param[in] Ni node i
//          * @param[in] Nj node j
//          * @param[in] P0 position of the new node onto which i and j will be
//          * contracted.
//          */
//            unsigned int temp(0);
//            const size_t i(Ni.sID); // StaticID of the source node
//            const size_t j(Nj.sID); // StaticID of the sink   node
//            
//            // collect all neighbors at P0 (but i and j)
//            std::set<size_t> neighbors;
//            Ni.neighborsAt(L0,neighbors);
//            Nj.neighborsAt(L0,neighbors);
//            //            neighbors.erase(i); // make sure
//            //            neighbors.erase(j); // make sure
//            
//            if(neighbors.size())
//            {
//                
//                // Remove all existing links among neighbors
//                for (std::set<size_t>::const_iterator nIter1=neighbors.begin();nIter1!=neighbors.end();++nIter1)
//                {
//                    for (std::set<size_t>::const_iterator nIter2=nIter1;nIter2!=neighbors.end();++nIter2)
//                    {
//                        if(nIter2!=nIter1)
//                        {
//                            if(DN.link(*nIter1,*nIter2).first || DN.link(*nIter2,*nIter1).first)
//                            {
//                                //std::cout<<"disconnecting "<<*nIter1<<" "<<*nIter2<<std::endl;
//                                DN.template disconnect<false>(*nIter1,*nIter2);
//                            }
//                        }
//                    }
//                }
//                
//                //std::cout<<"neighbors.size()="<<neighbors.size()<<std::endl;
//                
//                // Contract
//                for (std::set<size_t>::const_iterator nIter=neighbors.begin();nIter!=neighbors.end();++nIter)
//                {
//                    if(nIter!=neighbors.begin())
//                    {
//                        if(DN.node(*nIter).first)
//                        {
//                            //std::cout<<"neighbors contractingSecond "<<*neighbors.begin()<<" "<<*nIter<<std::endl;
//                            DN.contractSecond(*neighbors.begin(),*nIter);
//                            temp++;
//                        }
//                    }
//                }
//                
//                if(DN.node(i).first && DN.node(*neighbors.begin()).first && *neighbors.begin()!=i)
//                {
//                    //std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<i<<std::endl;
//                    
//                    DN.contractSecond(*neighbors.begin(),i);
//                    temp++;
//                }
//                
//                if(DN.node(j).first && DN.node(*neighbors.begin()).first && *neighbors.begin()!=j)
//                {
//                    //std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<j<<std::endl;
//                    
//                    DN.contractSecond(*neighbors.begin(),j);
//                    temp++;
//                }
//            }
//            else // neither i nor j has a neighbor at P0
//            {
//                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(L0.cartesian(),Ni.includingSimplex()));
//                if(guess.first) // check that P0 is inside mesh
//                {
//                    //std::cout<<"contracting "<<i<<" "<<j<<std::endl;
//                    
//                    //                    DN.contract(i,j,L0,guess.second); // 2 3 4 5
//                    DN.contract(i,j,L0); // 2 3 4 5
//                    temp++;
//                }
//            }
//            
//            return temp;
//        }
        
        
        /**********************************************************************/
        unsigned int contractSecondWithCommonNeighborCheck(const NodeType& Ni,
                                                           const NodeType& Nj)
        {/*!@param[in] i StaticID of the first node (vertex i remains)
          * @param[in] j StaticID of the second node (vertex j is contracted)
          *
          * Contracts  vertex j onto vertex i, making sure no other neighbors of j (but i)
          * occupies the position of i.
          */
            unsigned int temp(0);
            const int i(Ni.sID);
            const int j(Nj.sID);
            
            // collect all neighbors at Ni.get_P() (but j)
            std::set<size_t> neighbors;
            Nj.neighborsAt(Ni.get_L(),neighbors);
            Ni.neighborsAt(Ni.get_L(),neighbors);
            
            // Remove all existing links among neighbors
            for (std::set<size_t>::const_iterator nIter1=neighbors.begin();nIter1!=neighbors.end();++nIter1)
            {
                for (std::set<size_t>::const_iterator nIter2=nIter1;nIter2!=neighbors.end();++nIter2)
                {
                    if(nIter2!=nIter1)
                    {
                        if(DN.link(*nIter1,*nIter2).first || DN.link(*nIter2,*nIter1).first)
                        {
                            //std::cout<<"disconnecting "<<*nIter1<<" "<<*nIter2<<std::endl;
                            DN.template disconnect<false>(*nIter1,*nIter2);
                        }
                    }
                }
            }
            
            // Remove i from neighbors
            neighbors.erase(i);
            
            // Contract
            for (std::set<size_t>::const_iterator nIter=neighbors.begin();nIter!=neighbors.end();++nIter)
            {
                if(DN.node(*nIter).first)
                {
                    //std::cout<<"contractingSecond "<<i<<" "<<*nIter<<std::endl;
                    DN.contractSecond(i,*nIter);
                    temp++;
                }
            }
            
            if(DN.node(i).first && DN.node(j).first)
            {
                //std::cout<<"contractingSecond (j) "<<i<<" "<<j<<std::endl;
                DN.contractSecond(i,j);
                temp++;
            }
            
            return temp;
        }
        
        /**********************************************************************/
        unsigned int contractSecondWithCommonNeighborCheck(const int& i, const int& j)
        {/*! @param[in] i StaticID of the first node (vertex i remains)
          * @param[in] j StaticID of the second node (vertex j is contracted)
          *
          * Contracts  vertex j onto vertex i, making sure no other neighbors of j (but i)
          * occupies the position of i.
          */
            
            const typename DislocationNetworkType::isNetworkNodeType Ni(DN.node(i));
            assert(Ni.first && "NODE i DOES NOT EXIST");
            
            const typename DislocationNetworkType::isNetworkNodeType Nj(DN.node(j));
            assert(Nj.first && "NODE j DOES NOT EXIST");
            
            return contractSecondWithCommonNeighborCheck(*Ni.second,*Nj.second);
        }
        
        /**********************************************************************/
        std::pair<isNetworkNodeType,isNetworkNodeType>  selectSecond(const isNetworkNodeType& N1,
                                                                     const isNetworkNodeType& N2)
        {
            
            if(N1.second->isBoundaryNode() && !N2.second->isBoundaryNode())
            {
                return std::make_pair(N1,N2);
            }
            else if(!N1.second->isBoundaryNode() && N2.second->isBoundaryNode())
            {
                return std::make_pair(N2,N1);
            }
            else
            {
                if(N1.second->neighborhood().size()>N2.second->neighborhood().size())
                {
                    return std::make_pair(N1,N2);
                }
                else if(N1.second->neighborhood().size()<N2.second->neighborhood().size())
                {
                    return std::make_pair(N2,N1);
                }
                else
                {
                    if(N1.second->neighborLinksLength()<N2.second->neighborLinksLength())
                    {
                        return std::make_pair(N2,N1);
                    }
                    else if(N1.second->neighborLinksLength()>N2.second->neighborLinksLength())
                    {
                        return std::make_pair(N1,N2);
                    }
                    else
                    {
                        if(N1.second->sID<N2.second->sID)
                        {
                            return std::make_pair(N1,N2);
                        }
                        else
                        {
                            return std::make_pair(N2,N1);
                        }
                    }
                }
            }
        
        }
        
        
        
        
        
    public:
        
        /**********************************************************************/
        DislocationNodeContraction(DislocationNetworkType& DN_in) :
        // REENABLE THIS LATER        /* base init */ VertexContractionType(DN,DN),
        /* init list */ DN(DN_in)
        {}
        
        /**********************************************************************/
        //        size_t contractWithConstraintCheck(const size_t& i, const size_t& j)
        size_t contractWithConstraintCheck(const isNetworkNodeType& N1, const isNetworkNodeType& N2)
        {
            
            //std::cout<<"contractWithConstraintCheck "<<N1.second->sID<<" "<<N2.second->sID<<std::endl;
            
            size_t contracted(0);
            
            
            
            
            const LatticeVectorType& L1=N1.second->get_L();
            const LatticeVectorType& L2=N2.second->get_L();
            
            const bool isBoundaryNode1(N1.second->isBoundaryNode());
            const bool isBoundaryNode2(N2.second->isBoundaryNode());
            
            if((L1-L2).squaredNorm()==0) // points are coincindent, just contract them
            {
                //std::cout<<"contractWithConstraintCheck: chord= "<<(L1-L2).squaredNorm()<<std::endl;
                contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                
                //                contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,L1);
            }
            else // points are distinct
            {
                
                const typename DislocationNetworkType::NodeType::LatticePlaneContainerType& PN1=N1.second->confiningPlanes();
                const typename DislocationNetworkType::NodeType::LatticePlaneContainerType& PN2=N2.second->confiningPlanes();
                
                const size_t sizePN1(PN1.size());
                const size_t sizePN2(PN2.size());
                const VectorDimD P1=L1.cartesian();
                const VectorDimD P2=L2.cartesian();
                
                //std::cout<<"contractWithConstraintCheck: case "<<sizePN1<<" "<<sizePN2<<std::endl;
                
                if(sizePN1==1 && sizePN2==1)
                {// nodes constrained to move on planes, both inside mesh
                    PlanePlaneIntersection ppi(*PN1[0],*PN2[0]);
                    
                    if(!isBoundaryNode1 && !isBoundaryNode2)
                    {
                        const auto pr=selectSecond(N1,N2);
                        
                        //std::cout<<"contractWithConstraintCheck, case 1a"<<std::endl; // 2 3 5
                        switch (ppi.intersectionType)
                        {
                            case PlanePlaneIntersection::IntersectionType::coincident:
                            {
                                //std::cout<<"contractWithConstraintCheck, case 1a1"<<std::flush; // 2 3 5
                                LatticeVectorType L(PN1[0]->snapToLattice(0.5*(P1+P2)));

                                pr.first.second->set(L);
                                contracted+=contractSecondWithCommonNeighborCheck(*pr.first.second,*pr.second.second);
                                
//                                LatticeVectorType L(PN1[0]->snapToLattice(0.5*(P1+P2)));
//                                contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,L);
                                //std::cout<<" done"<<std::endl;
                                break;
                            }
                                
                            case PlanePlaneIntersection::IntersectionType::incident:
                            {
                                //std::cout<<"contractWithConstraintCheck, case 1a2"<<std::flush; // 2 3 5
                                // Find closest point on intersection line
                                const VectorDimD n1=PN1[0]->n.cartesian().normalized();
                                const VectorDimD n2=PN2[0]->n.cartesian().normalized();
                                const VectorDimD P(P1+P2);
                                const Eigen::Matrix<double,3,2> N((Eigen::Matrix<double,2,3>()<<n1.transpose(),n2.transpose()).finished().transpose());
                                const Eigen::Matrix<double,2,1> A((Eigen::Matrix<double,2,1>()<<P1.dot(n1),P2.dot(n2)).finished());
                                const Eigen::Matrix<double,2,1> lam=(N.transpose()*N).inverse()*(N.transpose()*P-2.0*A);
                                
                                const LatticeLine line(ppi.P,ppi.d);
                                LatticeVectorType L(line.snapToLattice(0.5*(P-N*lam)));
                                
                                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(L.cartesian(),N1.second->includingSimplex()));
                                if(guess.first)
                                {
                                    pr.first.second->set(L);
                                    contracted+=contractSecondWithCommonNeighborCheck(*pr.first.second,*pr.second.second);
                                }
                                
//                                contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,L);
                                //std::cout<<" done"<<std::endl;
                                
                                break;
                            }
                                
                            default:
                                break;
                        }
                    }
                    
                    else if(!isBoundaryNode1 && isBoundaryNode2)
                    {
                        contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                    }
                    else // case (isBoundaryNode1 && !isBoundaryNode2) and (isBoundaryNode1 && isBoundaryNode2)
                    { // N1 moves on its plane, constrained on the boundary, N2 moves on its plane, possibly on the boundary
                        //std::cout<<"contractWithConstraintCheck, case 1b"<<std::endl; // 2 3 5
                        switch (ppi.intersectionType)
                        {
                            case PlanePlaneIntersection::IntersectionType::coincident:
                            {
                                //std::cout<<"contractWithConstraintCheck, case 1b1"<<std::flush; // 2 3 5
                                contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                                //std::cout<<" done"<<std::endl;
                                break;
                            }
                                
                            case PlanePlaneIntersection::IntersectionType::incident:
                            {
                                //std::cout<<"contractWithConstraintCheck, case 1b2"<<std::flush; // 2 3 5
                                const LatticeLine line(ppi.P,ppi.d);
                                //                                const VectorDimD n1=PN1[0]->n.cartesian().normalized();
                                //                                const VectorDimD n2=PN2[0]->n.cartesian().normalized();
                                //                                const VectorDimD P(P1+P2);
                                //                                const Eigen::Matrix<double,3,2> N((Eigen::Matrix<double,2,3>()<<n1.transpose(),n2.transpose()).finished().transpose());
                                //                                const Eigen::Matrix<double,2,1> A((Eigen::Matrix<double,2,1>()<<P1.dot(n1),P2.dot(n2)).finished());
                                //                                const Eigen::Matrix<double,2,1> lam=(N.transpose()*N).inverse()*(N.transpose()*P-2.0*A);
                                
                                VectorDimD v0=line.snapToLattice(0.5*(P1+P2));
                                bool v0Inside=DN.shared.mesh.searchWithGuess(v0,N1.second->includingSimplex()).first;
                                if(!v0Inside)
                                {
                                    v0+=10.0*line.d.cartesian();
                                    v0Inside=DN.shared.mesh.searchWithGuess(v0,N1.second->includingSimplex()).first;
                                }
                                if(!v0Inside)
                                {
                                    v0-=20.0*line.d.cartesian();
                                    v0Inside=DN.shared.mesh.searchWithGuess(v0,N1.second->includingSimplex()).first;
                                }
                                if(v0Inside)
                                {
                                    const LatticeLine line2(LatticeVectorType(v0),line.d); // shift origin of line
                                    VectorDimD v1=line.snapToLattice(v0+10.0*line.d.cartesian());
                                    bool v1outside=!DN.shared.mesh.searchWithGuess(v1,N1.second->includingSimplex()).first;
                                    if(!v1outside)
                                    {
                                        v1=line.snapToLattice(v0-10.0*line.d.cartesian());
                                        v1outside=!DN.shared.mesh.searchWithGuess(v1,N1.second->includingSimplex()).first;
                                    }
                                    if(v1outside)
                                    {
                                        LineMeshIntersection lmi(line2,LatticeVectorType(v1),DN.shared.mesh,N1.second->includingSimplex());
                                        N1.second->set(lmi.L);
                                        contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                                        //contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,lmi.L);
                                        //std::cout<<" done"<<std::endl;
                                    }
                                }
                                //
                                //                                shared.mesh.searchWithGuess
                                ////                                LatticeVectorType L0(line.snapToLattice(N1.second->get_P()-5.0*line.d.cartesian().norm()*N1.second->bndNormal())); // here we could also say P2-f*P1.bndNormal, where f is a positive number
                                ////                                LatticeVectorType L1(line.snapToLattice(N1.second->get_P()+5.0*line.d.cartesian().norm()*N1.second->bndNormal())); // here we could also say P2-f*P1.bndNormal, where f is a positive number
                                //
                                //                                if((L0-L1).squaredNorm()!=0)
                                //                                {
                                ////                                    L=LatticeVectorType(line.snapToLattice(N2.second->get_P()-15.0*line.d.cartesian().norm()*N1.second->bndNormal()));
                                //                                    //                                contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,lineMeshIntersection(line2,N1.second->get_L(),N1.second->includingSimplex()));
                                //                                }
                                
                                break;
                            }
                                
                            default:
                                break;
                        }
                    }
                }
                else if(sizePN1==1 && sizePN2==2) // N1 moves on a plane, N2 moves on a line
                {
                    PlanePlaneIntersection ppi(*PN2[0],*PN2[1]);
                    assert(ppi.intersectionType==PlanePlaneIntersection::IntersectionType::incident);
                    const LatticeLine line(ppi.P,ppi.d);
                    PlaneLineIntersection pli(*PN1[0],line);

                    if(!isBoundaryNode1 && !isBoundaryNode2)
                    {
                        //std::cout<<"contractWithConstraintCheck, case 2a"<<std::endl; // 2 3 5
                        
                        // Find the LatticeLine over which node2 moves
                        
                        
                        // Find the intersection of the line and the plane of node1
                        switch (pli.intersectionType)
                        {
                            case PlaneLineIntersection::IntersectionType::intersecting:
                            {
                                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(pli.P.cartesian(),N1.second->includingSimplex()));
                                if(guess.first)
                                {
                                    N2.second->set(pli.P); // move node that is kept to intersection point
                                    contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                                }
 //                               contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,pli.P);
                                //std::cout<<"contractWithConstraintCheck, case 2a1"<<std::endl; // 2 3 5
                                
                                break;
                            }
                                
                            case PlaneLineIntersection::IntersectionType::coincident:
                                contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                                //std::cout<<"contractWithConstraintCheck, case 2a2"<<std::endl; // 2 3 5
                                
                                break;
                                
                            default:
                                //std::cout<<"contractWithConstraintCheck, case 2a3"<<std::endl; // 2 3 5
                                break;
                        }
                    }
//                    else if(isBoundaryNode1 && !isBoundaryNode2)
                    else
                    {// Intersection point must be on the boundary
                        //std::cout<<"contractWithConstraintCheck, case 2b"<<std::endl; // 2 3 5
//                        PlanePlaneIntersection ppi(*PN2[0],*PN2[1]);
//                        assert(ppi.intersectionType==PlanePlaneIntersection::IntersectionType::incident);
//                        const LatticeLine line(ppi.P,ppi.d);
//                        
//                        // Find the intersection of the line and the plane of node1
//                        PlaneLineIntersection pli(*PN1[0],line);
                        switch (pli.intersectionType)
                        {
                            case PlaneLineIntersection::IntersectionType::intersecting:
                            {// check if the intersection point is a boundary point
                                //std::cout<<"contractWithConstraintCheck, case 2b1"<<std::flush; // 2 3 5
                                
//                                std::pair<bool,const Simplex<dim,dim>*> temp=DN.shared.mesh.searchWithGuess(pli.P.cartesian(),N1.second->includingSimplex());
//                                std::pair<bool,const Simplex<dim,dim>*> tempp=DN.shared.mesh.searchWithGuess(pli.P.cartesian()+line.d.cartesian(),N1.second->includingSimplex());
//                                std::pair<bool,const Simplex<dim,dim>*> tempm=DN.shared.mesh.searchWithGuess(pli.P.cartesian()-line.d.cartesian(),N1.second->includingSimplex());
//                                if((temp.first && !tempp.first) || (temp.first && !tempm.first))
//                                {
//                                    contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,pli.P);
//                                    //std::cout<<" done"<<std::endl;
//                                    
//                                }
                                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(pli.P.cartesian(),N1.second->includingSimplex()));
                                if(guess.first && SimplexBndNormal::get_boundaryNormal(pli.P.cartesian(),*guess.second,NodeType::bndDistance).squaredNorm())
                                {
                                    N2.second->set(pli.P); // move node that is kept to intersection point
                                    contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                                }
                                
                                break;
                            }
                                
                            case PlaneLineIntersection::IntersectionType::coincident:
                            {// find the point where the line intersects the boundary
                                //std::cout<<"contractWithConstraintCheck, case 2b2"<<std::flush; // 2 3 5
                                
                                const LatticeLine line2(N2.second->get_L(),line.d);
                                const LatticeVectorType dL(line.d.snapToDirection(N1.second->bndNormal()*10.0));
                                if(dL.squaredNorm()>0)
                                {
                                    LineMeshIntersection lmi(line2,N2.second->get_L()+dL,DN.shared.mesh,N1.second->includingSimplex());
                                    //                                LatticeVectorType L=lineMeshIntersection(line2,N1.second->get_L(),N1.second->includingSimplex());
                                    if((lmi.L-N2.second->get_L()).squaredNorm()<2.0*(N1.second->get_L()-N2.second->get_L()).squaredNorm())
                                    {
                                        N2.second->set(lmi.L); // move node that is kept to intersection point
                                        contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
//                                        contracted+=contractWithCommonNeighborCheck(*N2.second,*N1.second,lmi.L);
                                    }
                                    //std::cout<<" done"<<std::endl;
                                    
                                }
                                break;
                            }
                                
                            default:
                                break;
                        }
                    }
//                    else if(!isBoundaryNode1 && isBoundaryNode2)
//                    {// N2 is stuck, N1 moves on a plane. Only option is that plane1 contains P2
//                        //std::cout<<"contractWithConstraintCheck, case 2c"<<std::flush; // 2 3 5
//                        if(PN1[0]->contains(L2))
//                        {
//                            contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
//                            //std::cout<<" done"<<std::endl;
//                            
//                        }
//                    }
//                    else
//                    {
//                        //std::cout<<"contractWithConstraintCheck, case 2d"<<std::endl; // 2 3 5
//                        
//                    }
                    
                }
                else if(sizePN1==2 && sizePN2==1) // N1 moves on a line, N2 moves on a plane
                {
                    //std::cout<<"contractWithConstraintCheck, case 3"<<std::endl; // 2 3 5
                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                }
                else if(sizePN1==2 && sizePN2==2) // both N1 and N2 move on lines
                {
                    
                    PlanePlaneIntersection ppi1(*PN1[0],*PN1[1]);
                    assert(ppi1.intersectionType==PlanePlaneIntersection::IntersectionType::incident);
                    const LatticeLine line1(ppi1.P,ppi1.d);
                    
                    
                    PlanePlaneIntersection ppi2(*PN2[0],*PN2[1]);
                    assert(ppi2.intersectionType==PlanePlaneIntersection::IntersectionType::incident);
                    const LatticeLine line2(ppi2.P,ppi2.d);
                    
                    LineLineIntersection lli(line1,line2);

                    
                    if(!isBoundaryNode1 && !isBoundaryNode2)
                    {
                        //std::cout<<"contractWithConstraintCheck, case 4a"<<std::endl; // 2 3 5

                        switch (lli.intersectionType)
                        {
                            case LineLineIntersection::IntersectionType::intersectingLines:
                            {
                                //std::cout<<"contractWithConstraintCheck, case 4a"<<std::flush; // 2 3 5
                                
                                const auto pr=selectSecond(N1,N2);
                                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(lli.P.cartesian(),pr.first.second->includingSimplex()));
                                if(guess.first)
                                {
                                    pr.first.second->set(lli.P); // move node that is kept to intersection point
                                    contracted+=contractSecondWithCommonNeighborCheck(*pr.first.second,*pr.second.second);//
                                }
                                
                                //                                contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,lli.P);
                                //std::cout<<" done"<<std::endl;
                                
                                break;
                            }
                                
                            case LineLineIntersection::IntersectionType::coincidentLines:
                            {
                                //std::cout<<"contractWithConstraintCheck, case 4b"<<std::flush; // 2 3 5
                                //std::cout<<"contractWithConstraintCheck, case 4a"<<std::endl; // 2 3 5
                                //std::cout<<" done"<<std::endl;
                                
                                //                                contracted+=contractSecondWithNodeSelection(N1,N2);
                                const auto pr=selectSecond(N1,N2);
                                contracted+=contractSecondWithCommonNeighborCheck(*pr.first.second,*pr.second.second);//
                                //                                if(N1.second->neighborhood().size()>N2.second->neighborhood().size())
                                //                                {
                                //                                    contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                                //
                                //                                }
                                //                                else if(N1.second->neighborhood().size()<N2.second->neighborhood().size())
                                //                                {
                                //                                    contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                                //                                }
                                //                                else
                                //                                {
                                //                                    contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,LatticeVectorType(line1.snapToLattice(0.5*(P1+P2))));
                                //                                }
                                
                                break;
                            }
                                
                                
                            default:
                                break;
                        }
                    }
                    else
                    { // Either nodes is on the boundary, so intersection must me on the boundary
                        switch (lli.intersectionType)
                        {
                            case LineLineIntersection::IntersectionType::intersectingLines:
                            {
                                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(lli.P.cartesian(),N1.second->includingSimplex()));
                                if(guess.first && SimplexBndNormal::get_boundaryNormal(lli.P.cartesian(),*guess.second,NodeType::bndDistance).squaredNorm())
                                {
                                    N1.second->set(lli.P); // move node that is kept to intersection point
                                    contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                                }
                                break;
                            }
                                
                            case LineLineIntersection::IntersectionType::coincidentLines:
                            {
                                const auto pr=selectSecond(N1,N2);
                                contracted+=contractSecondWithCommonNeighborCheck(*pr.first.second,*pr.second.second);
                                break;
                            }
                                
                                
                            default:
                                break;
                        }
                    }
//                    else if(isBoundaryNode1 && !isBoundaryNode2)
//                    {
//                        //std::cout<<"contractWithConstraintCheck, case 4c"<<std::flush; // 2 3 5
//                        
//                        // node 1 is stuck, while node2 moves on a line. Only option is thatl line2 contains P1
//                        PlanePlaneIntersection ppi2(*PN2[0],*PN2[1]);
//                        assert(ppi2.intersectionType==PlanePlaneIntersection::IntersectionType::incident);
//                        const LatticeLine line2(ppi2.P,ppi2.d);
//                        LineLineIntersection lli(line1,line2);

//
//                        
//                        
//                        //                        if(line2.contains(L1))
//                        //                        {
//                        //                            contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
//                        //                            //std::cout<<" done"<<std::endl;
//                        //
//                        //                        }
//                    }
//                    else if(!isBoundaryNode1 && isBoundaryNode2)
//                    {
//                        //                      //std::cout<<"contractWithConstraintCheck, case 4d"<<std::endl; // 2 3 5
//                        //std::cout<<"contractWithConstraintCheck, case 4 flip"<<std::flush; // 2 3 5
//                        
//                        contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
//                    }
//                    else // both are bonudary nodes
//                    {
//                        //std::cout<<"contractWithConstraintCheck, case 4d"<<std::endl; // 2 3 5
//                        
//                        // both nodes are stuck, so only option is then P1==P2, already considered above
//                    }
                }
                else if(sizePN1==1 && sizePN2==3)
                {// node2 is stuck,
                    if(!isBoundaryNode1)
                    {
                        //std::cout<<"contractWithConstraintCheck, case 5a"<<std::endl;
                        //std::cout<<isBoundaryNode1<<" "<<isBoundaryNode2<<std::endl;
                        if(PN1[0]->contains(L2))
                        {
                            //std::cout<<"contractWithConstraintCheck, case 5a1"<<std::flush; // 2 3 5
                            contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                            //std::cout<<" done"<<std::endl;
                            
                        }
                    }
                    else // node1 is a boundary node, so it can be contracted to N2 if N2 belongs to the plane of N1, and N2 is a boundary node
                    {
                        //std::cout<<"contractWithConstraintCheck, case 5b"<<std::endl;
                        if(PN1[0]->contains(L2) && isBoundaryNode2)
                        {
                            //std::cout<<"contractWithConstraintCheck, case 5b1"<<std::flush; // 2 3 5
                            contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                            //std::cout<<" done"<<std::endl;
                            
                        }
                    }
                }
                else if(sizePN1==3 && sizePN2==1)
                {
                    //std::cout<<"contractWithConstraintCheck, case 6 "<<std::flush;
                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                    //std::cout<<" done"<<std::endl;
                    
                }
                else if(sizePN1==2 && sizePN2==3)
                {// node2 is stuck,
                    //                    if(!isBoundaryNode1 || (isBoundaryNode1 && isBoundaryNode2))
                    if(isBoundaryNode1 && !isBoundaryNode2)
                    {
                        // do nothing, since N1 cannot be moved to the bulk
                    }
                    else
                    {// node1 is bulk and node2 is whatever, or both nodes are boundary
                        
                        //std::cout<<"contractWithConstraintCheck, case 7a"<<std::endl;
                        PlanePlaneIntersection ppi(*PN1[0],*PN1[1]);
                        assert(ppi.intersectionType==PlanePlaneIntersection::IntersectionType::incident);
                        const LatticeLine line(ppi.P,ppi.d);
                        if(line.contains(L2))
                        {
                            //std::cout<<"contractWithConstraintCheck, case 7aa"<<std::flush; // 2 3 5
                            contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                            //std::cout<<" done"<<std::endl;
                            
                        }
                    }
                    //                    else // node1 is a boundary node, so it is also stuck.
                    
                }
                else if(sizePN1==3 && sizePN2==2)
                {
                    //std::cout<<"contractWithConstraintCheck, case 8"<<std::flush; // 2 3 5
                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                    //std::cout<<" done"<<std::endl;
                    
                }
                //                else if(sizePN1==3 && sizePN2==3)
                //                {
                //                    if(P12norm<FLT_EPSILON)
                //                    {
                //                        contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                //                    }
                //                }
            }
            
            return contracted;
        }
        
    };
    
    
    //    template <typename DislocationNetworkType>
    //    double DislocationNodeContraction<DislocationNetworkType>::neighborRadius=0.001;
    
    
} // namespace model
#endif
