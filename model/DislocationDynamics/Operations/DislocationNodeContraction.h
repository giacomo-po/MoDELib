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
        
        static double neighborRadius;
        
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
        
        /**********************************************************************/
        unsigned int contractWithCommonNeighborCheck(const NodeType& Ni,
                                                     const NodeType& Nj,
                                                     const VectorDimD& P0)
        {/*!@param[in] Ni node i
          * @param[in] Nj node j
          * @param[in] P0 position of the new node onto which i and j will be
          * contracted.
          */
            unsigned int temp(0);
            const size_t i(Ni.sID); // StaticID of the source node
            const size_t j(Nj.sID); // StaticID of the sink   node
            
            // collect all neighbors at P0 (but i and j)
            std::set<size_t> neighbors;
            Ni.neighborsAt(P0,neighbors,neighborRadius);
            Nj.neighborsAt(P0,neighbors,neighborRadius);
            //            neighbors.erase(i); // make sure
            //            neighbors.erase(j); // make sure
            
            if(neighbors.size())
            {
                
                // Remove all existing links among neighbors
                for (std::set<size_t>::const_iterator nIter1=neighbors.begin();nIter1!=neighbors.end();++nIter1)
                {
                    for (std::set<size_t>::const_iterator nIter2=nIter1;nIter2!=neighbors.end();++nIter2)
                    {
                        if(nIter2!=nIter1)
                        {
                            if(DN.link(*nIter1,*nIter2).first || DN.link(*nIter2,*nIter1).first)
                            {
                                std::cout<<"disconnecting "<<*nIter1<<" "<<*nIter2<<std::endl;
                                DN.template disconnect<true>(*nIter1,*nIter2);
                            }
                        }
                    }
                }
                
                // Contract
                for (std::set<size_t>::const_iterator nIter=neighbors.begin();nIter!=neighbors.end();++nIter)
                {
                    if(nIter!=neighbors.begin())
                    {
                        if(DN.node(*nIter).first)
                        {
                            std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<*nIter<<std::endl;
                            DN.contractSecond(*neighbors.begin(),*nIter);
                            temp++;
                        }
                    }
                }
                
                if(DN.node(i).first && *neighbors.begin()!=i)
                {
                    std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<i<<std::endl;
                    
                    DN.contractSecond(*neighbors.begin(),i);
                    temp++;
                }
                
                if(DN.node(j).first && *neighbors.begin()!=j)
                {
                    std::cout<<"contractingSecond "<<*neighbors.begin()<<" "<<j<<std::endl;
                    
                    DN.contractSecond(*neighbors.begin(),j);
                    temp++;
                }
            }
            else // neither i nor j has a neighbor at P0
            {
                std::pair<bool,const Simplex<dim,dim>*> guess(DN.pointIsInsideMesh(P0,Ni.includingSimplex()));
                if(guess.first) // check that P0 is inside mesh
                {
                    std::cout<<"contracting "<<i<<" "<<j<<std::endl;
                    
                    DN.contract(i,j,P0,guess.second);
                    temp++;
                }
            }
            
            return temp;
        }
        
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
            Nj.neighborsAt(Ni.get_P(),neighbors,100.0*FLT_EPSILON);
            //            neighbors.erase(i);
            //            neighbors.erase(j);
            
            // Remove all existing links among neighbors
            for (std::set<size_t>::const_iterator nIter1=neighbors.begin();nIter1!=neighbors.end();++nIter1)
            {
                for (std::set<size_t>::const_iterator nIter2=nIter1;nIter2!=neighbors.end();++nIter2)
                {
                    if(nIter2!=nIter1)
                    {
                        if(DN.link(*nIter1,*nIter2).first || DN.link(*nIter2,*nIter1).first)
                        {
                            std::cout<<"disconnecting "<<*nIter1<<" "<<*nIter2<<std::endl;
                            DN.template disconnect<true>(*nIter1,*nIter2);
                        }
                    }
                }
            }
            
            // Contract
            neighbors.erase(i);
            
            //            if(neighbors.size())
            //            {
            for (std::set<size_t>::const_iterator nIter=neighbors.begin();nIter!=neighbors.end();++nIter)
            {
                if(DN.node(*nIter).first)
                {
                    std::cout<<"contractingSecond "<<i<<" "<<*nIter<<std::endl;
                    
                    DN.contractSecond(i,*nIter);
                    temp++;
                }
            }
            //            }
            
            if(DN.node(j).first)
            {
                std::cout<<"contractingSecond (j) "<<i<<" "<<j<<std::endl;
                
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
            
            std::cout<<"contractWithConstraintCheck "<<N1.second->sID<<" "<<N2.second->sID<<std::endl;

            
            assert(N1.first && "Non-existing node i. Change here");
            assert(N1.first && "Non-existing node j. Change here");
            
            size_t contracted(0);
            
            const VectorDimD P1=N1.second->get_P();
            const VectorDimD P2=N2.second->get_P();
            const VectorDimD P12=P1-P2;
            const double P12norm(P12.norm());
            
            if(P12norm<FLT_EPSILON) // points are coincindent, just contract them
            {
                std::cout<<"contractWithConstraintCheck: chord= "<<P12norm<<std::endl;

                contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
            }
            else // points are distinct
            {
                
                //                const VectorDimD unitP12=P12/P12norm;
                
                
                const typename DislocationNetworkType::NodeType::VectorOfNormalsType PN1(N1.second->constraintNormals());
                const typename DislocationNetworkType::NodeType::VectorOfNormalsType PN2(N2.second->constraintNormals());
                const size_t sizePN1(PN1.size());
                const size_t sizePN2(PN2.size());
                
                std::cout<<"contractWithConstraintCheck: case "<<sizePN1<<" "<<sizePN2<<std::endl;
                
                if(sizePN1==1 && sizePN2==1) // nodes constrained to move on planes
                {
                    std::cout<<"contractWithConstraintCheck, case 1"<<std::endl;
                    //                    const double denom(1.0-std::pow(PN1[0].dot(PN2[0]),2));
                    //                    const double numer((P2-P1).dot(PN2[0]));
                    //
                    //                    if(denom<FLT_EPSILON) // parallel or coincident planes
                    //                    {
                    //                        if(std::fabs(denom)<FLT_EPSILON) // planes are coincident
                    //                        {
                    //                            contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                    //                        }
                    //                        else // parallel planes
                    //                        {
                    //                            assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
                    //                        }
                    //                    }
                    //                    else // incident planes
                    //                    {
                    //                        contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,P1+(PN2[0]-PN2[0].dot(PN1[0])*PN1[0])*numer/denom);
                    //                    }
                    
                    //                    const VectorDimD d3(PN1[0].cross(PN2[0]));
                    //                    const double d3norm(d3.norm());
                    //                    if (d3norm<FLT_EPSILON) // parallel or coincident planes
                    if ((PN1[0].cross(PN2[0])).norm()<FLT_EPSILON) // parallel or coincident planes
                    {
                        std::cout<<"contractWithConstraintCheck, case 1a"<<std::endl;

                        if(std::fabs((P1-P2).dot(PN1[0]))<FLT_EPSILON) // planes are coincident
                        {
                            contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                        }
                        else // parallel planes
                        {
                            assert(0 && "COULD NOT CONTRACT JUNCTION POINTS.");
                        }
                    }
                    else // incident planes
                    {
                        std::cout<<"contractWithConstraintCheck, case 1b"<<std::endl;

                        // Contract at X, where X minimizes 0.5*(X-P1)^2+0.5*(X-P2)^2
                        // under the constraints (X-P1)*N1=0 and (X-P2)*N2=0
                        const VectorDimD P(P1+P2);
                        const Eigen::Matrix<double,3,2> N((Eigen::Matrix<double,2,3>()<<PN1[0].transpose(),PN2[0].transpose()).finished().transpose());
                        const Eigen::Matrix<double,2,1> A((Eigen::Matrix<double,2,1>()<<P1.dot(PN1[0]),P2.dot(PN2[0])).finished());
                        const Eigen::Matrix<double,2,1> L=(N.transpose()*N).inverse()*(N.transpose()*P-2.0*A);
                        contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P-N*L));
                    }
                    
                }
                else if(sizePN1==1 && sizePN2==2) // N1 moves on a plane, N2 moves on a line
                {
                    std::cout<<"contractWithConstraintCheck, case 2"<<std::endl;
                    const VectorDimD d2(PN2[0].cross(PN2[1]));
                    const double den(d2.dot(PN1[0]));
                    const double num(P12.dot(PN1[0]));
                    if(std::fabs(den)>FLT_EPSILON) // line and plane are not parallel
                    {
                        std::cout<<"contractWithConstraintCheck, case 2a"<<std::endl;
                        contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,P2+num/den*d2);
                    }
                    else
                    {
                        if(std::fabs(num)<FLT_EPSILON) // P2 belongs to the plane of P1
                        {
                            std::cout<<"contractWithConstraintCheck, case 2b"<<std::endl;
                            //                            contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2)); // HERE
                            contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                            
                        }
                        std::cout<<"contractWithConstraintCheck, case 2c"<<std::endl;
                    }
                }
                else if(sizePN1==2 && sizePN2==1) // N1 moves on a line, N2 moves on a plane
                {
                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                    //                    std::cout<<"contractWithConstraintCheck, case 3"<<std::endl;
                    //                    const VectorDimD d1(PN1[0].cross(PN1[1]));
                    //                    const double den(d1.dot(PN2[0]));
                    //                    const double num((P2-P1).dot(PN2[0]));
                    //                    if(std::fabs(den)>FLT_EPSILON) // line and plane are not parallel
                    //                    {
                    //                        std::cout<<"contractWithConstraintCheck, case 3a"<<std::endl;
                    //                        contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,P1+num/den*d1);
                    //                    }
                    //                    else
                    //                    {
                    //                        if(std::fabs(num)<FLT_EPSILON) // P1 belongs to the plane of P2
                    //                        {
                    //                            std::cout<<"contractWithConstraintCheck, case 3b"<<std::endl;
                    //
                    ////                            contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2)); // HERE
                    //                            contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                    //                        }
                    //                    }
                }
                else if(sizePN1==2 && sizePN2==2) // both N1 and N2 move on lines
                {
                    std::cout<<"contractWithConstraintCheck, case 4"<<std::endl;
                    
                    VectorDimD d1(PN1[0].cross(PN1[1]));
                    const double d1norm(d1.norm());
                    assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
                    d1/=d1norm;
                    VectorDimD d2(PN2[0].cross(PN2[1]));
                    const double d2norm(d2.norm());
                    assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
                    d2/=d2norm;
                    
                    const VectorDimD d3(d1.cross(d2));
                    const double d3Norm(d3.norm());
                    
                    if(d3Norm<FLT_EPSILON) // d1 and d2 are aligned, this means colinear or no intersection
                    {
                        if(d1.cross(P12.normalized()).norm()<FLT_EPSILON) // colinear
                        {
                            std::cout<<"contractWithConstraintCheck, case 4a"<<std::endl;
                            
                            contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,0.5*(P1+P2));
                        }
                        else // parallel (not colinear) lines, contraction not possible
                        {
                            //assert(0 && "COULD NOT CONTRACT JUNCTION POINTS. ALIGNMENT CONDITION FAILED.");
                            
                        }
                    }
                    else // d1 and d2 are not aligned
                    {
                        if(std::fabs((d3/d3Norm).dot(P12/P12norm))<FLT_EPSILON) // planarity condition
                        {
                            std::cout<<"contractWithConstraintCheck, case 4b"<<std::endl;
                            
                            const VectorDimD dOrth=d2-d2.dot(d1)*d1; // component of d2 orthogonal to d1
                            const double den=d2.dot(dOrth);
                            assert(std::fabs(den)>FLT_EPSILON && "YOU SHOULD HAVE FOUND THIS ABOVE.");
                            contracted+=contractWithCommonNeighborCheck(*N1.second,*N2.second,P2+P12.dot(dOrth)/den*d2);
                        }
                        else // planarity condition failed, contraction not possible
                        {
                            //assert(0 && "COULD NOT CONTRACT JUNCTION POINTS. PLANARITY CONDITION FAILED.");
                        }
                    }
                }
                else if(sizePN1==1 && sizePN2==3)
                {
                    std::cout<<"contractWithConstraintCheck, case 5"<<std::endl;
                    
                    if(std::fabs(P12.normalized().dot(PN1[0]))<FLT_EPSILON) // P2 belongs to the plane of P1
                    {// contract N1
                        contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                    }
                    else // contraction not possible
                    {
                        
                    }
                }
                else if(sizePN1==3 && sizePN2==1)
                {
                    std::cout<<"contractWithConstraintCheck, case 6"<<std::endl;
                    
                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                    
                    //                    if(std::fabs(P12.normalized().dot(PN2[0]))<FLT_EPSILON) // P1 belongs to the plane of P2
                    //                    {// contract N2
                    //                        contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                    //                    }
                    //                    else // contraction not possible
                    //                    {
                    //
                    //                    }
                }
                else if(sizePN1==2 && sizePN2==3)
                {
                    std::cout<<"contractWithConstraintCheck, case 7"<<std::endl;
                    
                    VectorDimD d1(PN1[0].cross(PN1[1]));
                    const double d1norm(d1.norm());
                    assert(d1norm>FLT_EPSILON && "DIRECTION d1 HAS ZERO NORM");
                    d1/=d1norm;
                    if(P12.normalized().cross(d1).norm()<FLT_EPSILON) // P2 belongs to the line of P1
                    {
                        contracted+=contractSecondWithCommonNeighborCheck(*N2.second,*N1.second);
                    }
                    else // contraction not possible
                    {
                        
                    }
                }
                else if(sizePN1==3 && sizePN2==2)
                {
                    contracted+=contractWithConstraintCheck(N2, N1); // call recursively switching N1 and N2
                    
                    //                    VectorDimD d2(PN2[0].cross(PN2[1]));
                    //                    const double d2norm(d2.norm());
                    //                    assert(d2norm>FLT_EPSILON && "DIRECTION d2 HAS ZERO NORM");
                    //                    d2/=d2norm;
                    //                    if(P12.normalized().cross(d2).norm()<FLT_EPSILON) // P2 belongs to the line of P1
                    //                    {
                    //                        contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                    //                    }
                    //                    else // contraction not possible
                    //                    {
                    //
                    //                    }
                }
                else if(sizePN1==3 && sizePN2==3)
                {
                    if(P12norm<FLT_EPSILON)
                    {
                        contracted+=contractSecondWithCommonNeighborCheck(*N1.second,*N2.second);
                    }
                }
            }
            
            return contracted;
        }
        
    };


    template <typename DislocationNetworkType>
    double DislocationNodeContraction<DislocationNetworkType>::neighborRadius=0.001;

    
} // namespace model
#endif
