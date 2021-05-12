/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationNetworkComponent_H_
#define model_DislocationNetworkComponent_H_

#include <vector>
#include <deque>
#include <float.h>
#include <Eigen/Dense>
#include <NetworkComponent.h>
#include <TextFileParser.h>

//#include <SchurComplementSolver.h>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/SparseCholesky> // simplicial DLDT

#ifdef _MODEL_PARDISO_SOLVER_
#include <Eigen/PardisoSupport>
#endif


namespace model
{
    
    template <typename NodeType,typename LinkType>
    class DislocationNetworkComponent
    {
        
        
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        typedef typename NetworkComponentType::NetworkComponentNodeContainerType NodeContainerType;
        typedef typename NetworkComponentType::NetworkComponentLinkContainerType LinkContainerType;
        typedef typename NodeType::VectorDim   VectorDim;
        typedef typename NodeType::MatrixDim   MatrixDim;
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        //        typedef std::vector<Eigen::Triplet<double> > TripletContainerType;
        typedef std::deque<Eigen::Triplet<double> > TripletContainerType;
        
        enum {dim=NodeType::dim};
        
        
        
        
        //! A reference to a NetworkComponentType
        NetworkComponentType& NC;
        
        /**********************************************************************/
        void storeNodeSolution(const Eigen::VectorXd& X)
        {
            size_t k=0;
            for (typename NodeContainerType::iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
            {
                nodeIter->second->set_V(X.segment(NdofXnode*k,NdofXnode));
                ++k;
            }
        }
        
        
        /************************************************************/
        //        template <bool symmetricConstraint>
        void assembleConstraints(TripletContainerType& vT, size_t& KPQ_row, const bool& symmetricConstraint) const
        {
            
            //            // Collect contraints
            //            std::vector<typename NodeType::VectorOfNormalsType> NNV; // DON'T NEED THIS ANYMORE WITH TRIPLETS
            //            for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
            //            {
            //                NNV.push_back(nodeIter->second->constraintNormals());
            //            }
            //
            //            for (size_t n=0 ;n<NNV.size();++n)
            //            {
            //                for (size_t c=0;c<NNV[n].size();++c)
            //                {
            //                    for(size_t d=0;d<dim;++d)
            //                    {
            //                        vT.push_back(Eigen::Triplet<double>(KPQ_row,n*dim+d,NNV[n][c](d)));
            //                        if (symmetricConstraint)
            //                        {
            //                            vT.push_back(Eigen::Triplet<double>(n*dim+d,KPQ_row,NNV[n][c](d)));
            //                        }
            //                    }
            //                    ++KPQ_row; // move to next line
            //                }
            //            }
            
            //            size_t n=0;
            //            for (const auto& node : NC.nodes())
            //            {
            //                const typename NodeType::VectorOfNormalsType vn=node.second->constraintNormals();
            //                for (size_t c=0;c<vn.size();++c)
            //                {
            //                    for(size_t d=0;d<dim;++d)
            //                    {
            ////                        vT.push_back(Eigen::Triplet<double>(KPQ_row,n*dim+d,vn[c](d)));
            //                        vT.emplace_back(KPQ_row,n*dim+d,vn[c](d));
            //                        if (symmetricConstraint)
            //                        {
            ////                            vT.push_back(Eigen::Triplet<double>(n*dim+d,KPQ_row,vn[c](d)));
            //                              vT.emplace_back(n*dim+d,KPQ_row,vn[c](d));
            //                        }
            //                    }
            //                    ++KPQ_row; // move to next line
            //                }
            //                n++;
            //            }
            
            //            //Constrain simple nodes to move normal to tangent
            //            for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
            //            {
            //                if(nodeIter->second->constraintNormals().size()==1)
            //                {
            //                    Eigen::Matrix<int,dim,1> node_dofID(nodeIter->second->node_dofID());
            //                    Eigen::Matrix<double,dim,1> T(nodeIter->second->get_T());
            //                    double normT(T.norm());
            //                    if (normT>FLT_EPSILON)
            //                    {
            //                        T/=normT;
            //                        for(size_t d=0;d<dim;++d)
            //                        {
            ////                            vT.push_back(Eigen::Triplet<double>(KPQ_row,node_dofID(d),T(d)));
            //                            vT.emplace_back(KPQ_row,node_dofID(d),T(d));
            //                            if (symmetricConstraint)
            //                            {
            ////                                vT.push_back(Eigen::Triplet<double>(node_dofID(d),KPQ_row,T(d)));
            //                                vT.emplace_back(node_dofID(d),KPQ_row,T(d));
            //                            }
            //                        }
            //                        ++KPQ_row;
            //                    }
            //                }
            //            }
            
            
            // loop over each segment and add segment contributions to kqqT and Fq
            //            for (typename LinkContainerType::const_iterator linkIter=NC.linkBegin();linkIter!=NC.linkEnd();++linkIter){
            //                const Eigen::VectorXd ortCS(linkIter->second->ortC*linkIter->second->Mseg);
            //                for (unsigned int i=0;i<linkIter->second->segmentDOFs.size();++i){
            //                    std::set<size_t>::const_iterator iterI(linkIter->second->segmentDOFs.begin());
            //                    std::advance(iterI,i);
            //                    //FQ(*iterI)+=tempFq(i);
            //                    kpqT.push_back(Eigen::Triplet<double>(KPQ_row,*iterI,ortCS(i)));
            //                }
            //                ++KPQ_row;
            //			}
            
            
            //size_t KPQ_row=0;
            
            
            
        }
   //Added by Yash
        //Improved Version 2
//        size_t assembleConstraintsforPeriodicSimulationsNULL(TripletContainerType& zT) const
//        {
//            // typedef std::map <std::pair< NodeType *, NodeType * >, std::tuple<NodeType *,double, NodeType *,double>> PeriodicConnectivityType;
//            size_t constrainedI = 0;
//            size_t unconstrainedNodes = 0;
//            std::map<size_t, size_t> correctedJPosition;
//            for (const auto &node : NC.nodes())
//            {
//
//                if (node.second->isBoundaryNode())
//                {
//                    const auto linkMapBNodes = node.second->getNodesMapforPeriodicAssembly();
//                    bool temp(std::get<0>(linkMapBNodes) && std::get<1>(linkMapBNodes) && std::get<2>(linkMapBNodes));
//                    if (temp)
//                    {
//                        constrainedI++;
//                    }
//                    else
//                    {
//                        const size_t ntempsnID(node.second->snID()); //Global position in the constraint matrix (j)
//                        correctedJPosition.emplace(ntempsnID, ntempsnID - constrainedI);
//                    }
//
//                }
//                else
//                {
//                    const size_t ntempsnID(node.second->snID()); //Global position in the constraint matrix (j)
//                    correctedJPosition.emplace(ntempsnID, ntempsnID - constrainedI);
//                }
//            }
//            for (const auto &node : NC.nodes())
//            {
//                const size_t ntempj(node.second->snID()); //Global position in the constraint matrix (j)
//                // std::cout<<"node.second->sID(is a Boundary node) "<<node.second->sID<<" ("<<node.second->isBoundaryNode()<<" )"<<std::endl;
//                if (node.second->isBoundaryNode())
//                {
//                    const auto linkMapBNodes = node.second->getNodesMapforPeriodicAssembly();
//                    bool temp(std::get<0>(linkMapBNodes) && std::get<1>(linkMapBNodes) && std::get<2>(linkMapBNodes) );
//                    if (temp)
//                    {
//                        assert(std::get<1>(linkMapBNodes)->sID == node.second->sID);
//                        const double lij(std::get<3>(linkMapBNodes));
//                        const double ljk(std::get<4>(linkMapBNodes));
//
//                        size_t ntempi(std::get<0>(linkMapBNodes)->snID()); //Global position in the constraint matrix corresponding to the original link (i)
//                        size_t ntempk(std::get<2>(linkMapBNodes)->snID()); //Global position in the constraint matrix corresponding to the neighbor link (k)
//                        const double lijk(lij + ljk);
//
//                        //First check if ntempi and ntempj are both non-boundary nodes, if not, the constraints are needed to be multiplied and ntempi and ntempk updated
//                        // std::cout << "Ntempi node " << std::get<0>(nodeinLinks.second)->sID << std::endl;
//                        // std::cout << "Ntempj node " << node.second->sID << std::endl;
//                        // std::cout << "Ntempk node " << std::get<2>(nodeinLinks.second)->sID << std::endl;
//
//                        assert(correctedJPosition.find(ntempi) != correctedJPosition.end());
//                        assert(correctedJPosition.find(ntempk) != correctedJPosition.end());
//
//                        const size_t correctedI(correctedJPosition.find(ntempi)->second);
//                        const size_t correctedK(correctedJPosition.find(ntempk)->second);
//
//                        for (size_t d = 0; d < dim; ++d)
//                        {
//                            // zT.emplace_back(dim*ntempj+d,dim*ntempi+d,ljk/lijk);
//
//                            zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
//                            zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
//                        }
//                    }
//                    else
//                    {
//                        assert(correctedJPosition.find(node.second->snID()) != correctedJPosition.end());
//
//                        const size_t correctedJ(correctedJPosition.find(node.second->snID())->second);
//                        for (size_t d = 0; d < dim; ++d)
//                        {
//                            zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
//                        }
//                        unconstrainedNodes++;
//                    }
//                }
//                else
//                {
//                    assert(correctedJPosition.find(node.second->snID()) != correctedJPosition.end());
//                    const size_t correctedJ(correctedJPosition.find(node.second->snID())->second);
//                    for (size_t d = 0; d < dim; ++d)
//                    {
//                        zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
//                    }
//                    unconstrainedNodes++;
//                }
//            }
//            return unconstrainedNodes;
//        }
        // //Added by Yash
        // //Improved Version
        // size_t assembleConstraintsforPeriodicSimulationsNULL(TripletContainerType& zT) const
        // {
        //     typedef std::map <std::pair< NodeType *, NodeType * >, std::tuple<NodeType *,double, NodeType *,double>> PeriodicConnectivityType;
        //     size_t constrainedI=0;
        //     size_t unconstrainedNodes=0;
        //     std::map<size_t,size_t> correctedJPosition;
        //     for (const auto &node : NC.nodes())
        //     {

        //         if (node.second->isBoundaryNode() && node.second->getNodesMapforPeriodicAssembly().size()>0)
        //         {
        //             constrainedI++;
        //         }
        //         else
        //         {
        //             const size_t ntempsnID(node.second->snID()); //Global position in the constraint matrix (j)
        //             correctedJPosition.emplace(ntempsnID,ntempsnID-constrainedI);
        //         }
                
        //     }
        //     for (const auto &node : NC.nodes())
        //     {
        //         const size_t ntempj(node.second->snID()); //Global position in the constraint matrix (j)
        //         // std::cout<<"node.second->sID(is a Boundary node) "<<node.second->sID<<" ("<<node.second->isBoundaryNode()<<" )"<<std::endl;
        //         if (node.second->isBoundaryNode())
        //         {
        //             const auto linkMapBNodes = node.second->getNodesMapforPeriodicAssembly();

        //             // if (linkMapBNodes.size()==1 || linkMapBNodes.size()==0)
        //             // {
        //             //     for (const auto &links : linkMapBNodes)
        //             //     {
        //             //         std::cout<<"link.first "<<links.first->source->sID<< "-----> "<<links.first->sink->sID<<std::endl;
        //             //         std::cout<<"link.second "<<links.second->source->sID<< "-----> "<<links.second->sink->sID<<std::endl; 

        //             //     }
        //             // }
        //             if (linkMapBNodes.size()>0)
        //             {
        //                 assert((linkMapBNodes.size() == 1 ) && "FINISH HERE");
        //                 for (const auto &nodeinLinks : linkMapBNodes)
        //                 {
        //                     const double ljk(std::get<3>(nodeinLinks.second));
        //                     const double lij(std::get<1>(nodeinLinks.second));

        //                     size_t ntempi(std::get<0>(nodeinLinks.second)->snID()); //Global position in the constraint matrix corresponding to the original link (i)
        //                     size_t ntempk(std::get<2>(nodeinLinks.second)->snID()); //Global position in the constraint matrix corresponding to the neighbor link (k)
        //                     const double lijk(lij + ljk);

        //                     //First check if ntempi and ntempj are both non-boundary nodes, if not, the constraints are needed to be multiplied and ntempi and ntempk updated
        //                     std::cout << "Ntempi node " << std::get<0>(nodeinLinks.second)->sID << std::endl;
        //                     std::cout << "Ntempj node " << node.second->sID << std::endl;
        //                     std::cout << "Ntempk node " << std::get<2>(nodeinLinks.second)->sID << std::endl;

        //                     assert(correctedJPosition.find(ntempi) != correctedJPosition.end());
        //                     assert(correctedJPosition.find(ntempk) != correctedJPosition.end());

        //                     const size_t correctedI(correctedJPosition.find(ntempi)->second);
        //                     const size_t correctedK(correctedJPosition.find(ntempk)->second);

        //                     for (size_t d = 0; d < dim; ++d)
        //                     {
        //                         // zT.emplace_back(dim*ntempj+d,dim*ntempi+d,ljk/lijk);

        //                         zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
        //                     }
        //                 }
        //             }
        //             else
        //             {
        //                 assert(correctedJPosition.find(node.second->snID()) != correctedJPosition.end());

        //                 const size_t correctedJ(correctedJPosition.find(node.second->snID())->second);
        //                 for (size_t d = 0; d < dim; ++d)
        //                 {
        //                     zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
        //                 }
        //                 unconstrainedNodes++;
        //             }
                    
                    
        //         }
        //         else
        //         {
        //             assert(correctedJPosition.find(node.second->snID()) != correctedJPosition.end());

        //             const size_t correctedJ(correctedJPosition.find(node.second->snID())->second);
        //             for (size_t d = 0; d < dim; ++d)
        //             {
        //                 zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
        //             }
        //             unconstrainedNodes++;
        //         }
        //     }
        //     return unconstrainedNodes;
        // }
        // size_t assembleConstraintsforPeriodicSimulationsNULL(TripletContainerType& zT) const
        // {
        //     size_t constrainedI=0;
        //     size_t unconstrainedNodes=0;
        //     std::map<size_t,size_t> correctedJPosition;
        //     for (const auto &node : NC.nodes())
        //     {

        //         if (node.second->isBoundaryNode())
        //         {
        //             constrainedI++;
        //         }
        //         else
        //         {
        //             const size_t ntempsnID(node.second->snID()); //Global position in the constraint matrix (j)
        //             correctedJPosition.emplace(ntempsnID,ntempsnID-constrainedI);
        //         }
                
        //     }
        //     for (const auto &node : NC.nodes())
        //     {
        //         const size_t ntempj(node.second->snID()); //Global position in the constraint matrix (j)
        //         // std::cout<<"node.second->sID(is a Boundary node) "<<node.second->sID<<" ("<<node.second->isBoundaryNode()<<" )"<<std::endl;
        //         if (node.second->isBoundaryNode())
        //         {
        //             const auto linkMapBNodes = node.second->getLinksMapsforBoundaryNodes();
        //             // if (linkMapBNodes.size()==1 || linkMapBNodes.size()==0)
        //             // {
        //             //     for (const auto &links : linkMapBNodes)
        //             //     {
        //             //         std::cout<<"link.first "<<links.first->source->sID<< "-----> "<<links.first->sink->sID<<std::endl;
        //             //         std::cout<<"link.second "<<links.second->source->sID<< "-----> "<<links.second->sink->sID<<std::endl; 

        //             //     }
        //             // }
        //             assert((linkMapBNodes.size()==1 || linkMapBNodes.size()==0) && "FINISH HERE");
        //             for (const auto& links : linkMapBNodes)
        //             {
        //                 const double ljk((links.first->source->get_P() - links.first->sink->get_P()).norm());
        //                 const double lij((links.second->source->get_P() - links.second->sink->get_P()).norm());
        //                 const double lijk(lij+ljk);

        //                 if (links.first->source->isBoundaryNode() && links.second->source->isBoundaryNode())
        //                 {
        //                     const size_t ntempi(links.second->sink->snID()); //Global position in the constraint matrix corresponding to the original link (i)
        //                     const size_t ntempk(links.first->sink->snID());  //Global position in the constraint matrix corresponding to the neighbor link (k)
        //                     std::cout<<"link.first source->sink "<<links.first->source->sID<<" ( "<<links.first->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.first->sink->sID<<" ( "<<links.first->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     std::cout<<"link.second source->sink "<<links.second->source->sID<<" ( "<<links.second->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.second->sink->sID<<" ( "<<links.second->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     assert(correctedJPosition.find(links.second->sink->snID())!=correctedJPosition.end());
        //                     assert(correctedJPosition.find(links.first->sink->snID())!=correctedJPosition.end());

        //                     const size_t correctedI(correctedJPosition.find(links.second->sink->snID())->second);
        //                     const size_t correctedK(correctedJPosition.find(links.first->sink->snID())->second);

        //                     for (size_t d = 0; d < dim; ++d)
        //                     {
        //                         // zT.emplace_back(dim*ntempj+d,dim*ntempi+d,ljk/lijk);

        //                         zT.emplace_back(dim*ntempj+d,dim*correctedI+d,ljk/lijk);
        //                         zT.emplace_back(dim*ntempj+d,dim*correctedK+d,lij/lijk);
        //                     }
        //                 }
        //                 else if (links.first->sink->isBoundaryNode() && links.second->source->isBoundaryNode())
        //                 {
        //                     const size_t ntempi(links.second->sink->snID());   //Global position in the constraint matrix corresponding to the original link (i)
        //                     const size_t ntempk(links.first->source->snID());    //Global position in the constraint matrix corresponding to the neighbor link (k)
        //                     std::cout<<"link.first source->sink "<<links.first->source->sID<<" ( "<<links.first->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.first->sink->sID<<" ( "<<links.first->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     std::cout<<"link.second source->sink "<<links.second->source->sID<<" ( "<<links.second->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.second->sink->sID<<" ( "<<links.second->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     assert(correctedJPosition.find(links.second->sink->snID()) != correctedJPosition.end());
        //                     assert(correctedJPosition.find(links.first->source->snID()) != correctedJPosition.end());

        //                     const size_t correctedI(correctedJPosition.find(links.second->sink->snID())->second);
        //                     const size_t correctedK(correctedJPosition.find(links.first->source->snID())->second);
        //                     for (size_t d = 0; d < dim; ++d)
        //                     {
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
        //                     }
        //                 }
        //                 else if (links.first->source->isBoundaryNode() && links.second->sink->isBoundaryNode())
        //                 {
        //                     // const size_t ntempj(links.second->sink->snID()); //Global position in the constraint matrix (j)
        //                     const size_t ntempi(links.second->source->snID());   //Global position in the constraint matrix corresponding to the original link (i)
        //                     const size_t ntempk(links.first->sink->snID());  //Global position in the constraint matrix corresponding to the neighbor link (k)
        //                     std::cout<<"link.first source->sink "<<links.first->source->sID<<" ( "<<links.first->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.first->sink->sID<<" ( "<<links.first->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     std::cout<<"link.second source->sink "<<links.second->source->sID<<" ( "<<links.second->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.second->sink->sID<<" ( "<<links.second->sink->isBoundaryNode()<<" )"<<std::endl;

        //                     assert(correctedJPosition.find(links.second->source->snID()) != correctedJPosition.end());
        //                     assert(correctedJPosition.find(links.first->sink->snID()) != correctedJPosition.end());

        //                     const size_t correctedI(correctedJPosition.find(links.second->source->snID())->second);
        //                     const size_t correctedK(correctedJPosition.find(links.first->sink->snID())->second);
                            
                            
        //                     for (size_t d = 0; d < dim; ++d)
        //                     {
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
        //                     }
        //                 }
        //                 else if (links.first->sink->isBoundaryNode() && links.second->sink->isBoundaryNode())
        //                 {
        //                     // const size_t ntempj(links.second->sink->snID()); //Global position in the constraint matrix (j)
        //                     const size_t ntempi(links.second->source->snID());   //Global position in the constraint matrix corresponding to the original link (i)
        //                     const size_t ntempk(links.first->source->snID());  //Global position in the constraint matrix corresponding to the neighbor link (k)
        //                     std::cout<<"link.first source->sink "<<links.first->source->sID<<" ( "<<links.first->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.first->sink->sID<<" ( "<<links.first->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     std::cout<<"link.second source->sink "<<links.second->source->sID<<" ( "<<links.second->source->isBoundaryNode()
        //                     <<" )"<<" ---> "<<links.second->sink->sID<<" ( "<<links.second->sink->isBoundaryNode()<<" )"<<std::endl;
        //                     assert(correctedJPosition.find(links.second->source->snID()) != correctedJPosition.end());
        //                     assert(correctedJPosition.find(links.first->source->snID()) != correctedJPosition.end());

        //                     const size_t correctedI(correctedJPosition.find(links.second->source->snID())->second);
        //                     const size_t correctedK(correctedJPosition.find(links.first->source->snID())->second);
                            
        //                     for (size_t d = 0; d < dim; ++d)
        //                     {
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
        //                         zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
        //                     }
        //                 }
        //             }
        //         }
        //         else
        //         {
        //             assert(correctedJPosition.find(node.second->snID()) != correctedJPosition.end());

        //             const size_t correctedJ(correctedJPosition.find(node.second->snID())->second);
        //             for (size_t d = 0; d < dim; ++d)
        //             {
        //                 zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
        //             }
        //             unconstrainedNodes++;
        //         }
        //     }
        //     return unconstrainedNodes;
        // }

        /************************************************************/
        //Added by Yash
//        void assembleConstraintsforPeriodicSimulations(TripletContainerType &cT, const double &pfactor) const
//        {
//            for (const auto &node : NC.nodes())
//            {
//                if (node.second->isBoundaryNode())
//                {
//                    const auto linkMapBNodes = node.second->getLinksMapsforBoundaryNodes();
//                    assert(linkMapBNodes.size() == 1 && "FINISH HERE");
//                    for (const auto &links : linkMapBNodes)
//                    {
//                        const double ljk((links.first->source->get_P() - links.first->sink->get_P()).norm());
//                        const double lij((links.second->source->get_P() - links.second->sink->get_P()).norm());
//                        const double lijk(lij + ljk);
//                        if (links.first->source->isBoundaryNode() && links.second->source->isBoundaryNode())
//                        {
//                            const size_t ntempj(links.second->source->snID()); //Global position in the constraint matrix (j)
//                            const size_t ntempi(links.second->sink->snID());   //Global position in the constraint matrix corresponding to the original link (i)
//                            const size_t ntempk(links.first->sink->snID());    //Global position in the constraint matrix corresponding to the neighbor link (k)
//                            //assemble from the sinks of both links
//                            for (size_t d = 0; d < dim; ++d)
//                            {
//                                cT.emplace_back(dim * ntempj + d, dim * ntempj + d, pfactor * 1.0);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempi + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempi + d, dim * ntempj + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempk + d, -pfactor * lij / lijk);
//                                cT.emplace_back(dim * ntempk + d, dim * ntempj + d, -pfactor * lij / lijk);
//                            }
//                        }
//                        else if (links.first->sink->isBoundaryNode() && links.second->source->isBoundaryNode())
//                        {
//                            const size_t ntempj(links.second->source->snID()); //Global position in the constraint matrix (j)
//                            const size_t ntempi(links.second->sink->snID());   //Global position in the constraint matrix corresponding to the original link (i)
//                            const size_t ntempk(links.first->source->snID());  //Global position in the constraint matrix corresponding to the neighbor link (k)
//                            //assemble from the sinks of both links
//                            for (size_t d = 0; d < dim; ++d)
//                            {
//                                cT.emplace_back(dim * ntempj + d, dim * ntempj + d, pfactor * 1.0);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempi + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempi + d, dim * ntempj + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempk + d, -pfactor * lij / lijk);
//                                cT.emplace_back(dim * ntempk + d, dim * ntempj + d, -pfactor * lij / lijk);
//                            }
//                        }
//                        else if (links.first->source->isBoundaryNode() && links.second->sink->isBoundaryNode())
//                        {
//                            const size_t ntempj(links.second->sink->snID());   //Global position in the constraint matrix (j)
//                            const size_t ntempi(links.second->source->snID()); //Global position in the constraint matrix corresponding to the original link (i)
//                            const size_t ntempk(links.first->sink->snID());    //Global position in the constraint matrix corresponding to the neighbor link (k)
//                            //assemble from the sinks of both links
//                            for (size_t d = 0; d < dim; ++d)
//                            {
//                                cT.emplace_back(dim * ntempj + d, dim * ntempj + d, pfactor * 1.0);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempi + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempi + d, dim * ntempj + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempk + d, -pfactor * lij / lijk);
//                                cT.emplace_back(dim * ntempk + d, dim * ntempj + d, -pfactor * lij / lijk);
//                            }
//                        }
//                        else if (links.first->sink->isBoundaryNode() && links.second->sink->isBoundaryNode())
//                        {
//                            const size_t ntempj(links.second->sink->snID());   //Global position in the constraint matrix (j)
//                            const size_t ntempi(links.second->source->snID()); //Global position in the constraint matrix corresponding to the original link (i)
//                            const size_t ntempk(links.first->source->snID());  //Global position in the constraint matrix corresponding to the neighbor link (k)
//                            //assemble from the sinks of both links
//                            for (size_t d = 0; d < dim; ++d)
//                            {
//                                cT.emplace_back(dim * ntempj + d, dim * ntempj + d, pfactor * 1.0);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempi + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempi + d, dim * ntempj + d, -pfactor * ljk / lijk);
//                                cT.emplace_back(dim * ntempj + d, dim * ntempk + d, -pfactor * lij / lijk);
//                                cT.emplace_back(dim * ntempk + d, dim * ntempj + d, -pfactor * lij / lijk);
//                            }
//                        }
//                    }
//                }
//            }
//        }
        /************************************************************/
        size_t assembleNCtriplets(TripletContainerType& kqqT, Eigen::VectorXd& Fq)
        {
//            const size_t nodeOrdr(NC.nodeOrder());
//            if (nodeOrdr==0 || nodeOrdr==1)
//            {
//                std::cout<<"DislocationNetworkComponent:" <<NC.sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<NC.linkOrder()<<std::endl;
//                typename NodeContainerType::const_iterator nodeIter1=NC.nodeBegin();
//                std::cout<<"Only DIslocatioNode is "<< nodeIter1->second->sID<<std::endl;
////                assert(0 && "DislocationSubNetework has less than 2 Nodes.");
//                std::cout<<"DislocationSubNetework has less than 2 Nodes."<<std::endl;
//            }
//            if (nodeOrdr==2)
//            {
//                typename NodeContainerType::const_iterator nodeIter1=NC.nodeBegin();
//                typename NodeContainerType::const_iterator nodeIter2(nodeIter1);
//                nodeIter2++;
//                if (nodeIter1->second->constraintNormals().size()>2 && nodeIter2->second->constraintNormals().size()>2)
//                {
//                    //                    assert(0 && "DislocationSubNetework has only 2 fixed Nodes.");
//                    std::cout<<"DislocationNetworkComponent:" <<NC.sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<NC.linkOrder()<<std::endl;
//                    std::cout<<"First DIslocatioNode is " <<nodeIter1->second->sID<<std::endl;
//                    std::cout<<"Second DIslocatioNode is "<<nodeIter2->second->sID<<std::endl;
//                    std::cout<<"WARNING: DislocationSubNetework has only 2 fixed Nodes"<<std::endl;
//                    
//                }
//            }
            
            
            // Assembly of Stiffness Matrix and force vector
            const size_t Ndof(NC.nodeOrder()*NdofXnode); // the total number of dof in the subnetwork
            //            size_t reserveSize(0);
            //            for (typename NodeContainerType::iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
            //            {
            //                reserveSize+=nodeIter->second->closedOrder();
            //            }
            //            //            TripletContainerType kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
            //            kqqT.reserve(reserveSize); // assume some fill percentage
            
            
            //            Eigen::VectorXd Fq(Eigen::VectorXd::Zero(Ndof));
            Fq.setZero(Ndof);
            
            for (typename LinkContainerType::const_iterator linkIter=NC.linkBegin();linkIter!=NC.linkEnd();++linkIter)
            {
                if (!linkIter->second->isBoundarySegment())
                {
                    linkIter->second->addToGlobalAssembly(kqqT,Fq); // loop over each segment and add segment contributions to kqqT and Fq
                }
            }
            
            return Ndof;
            
        }
        
        
    public:
        
        enum {NdofXnode=NodeType::NdofXnode};
        
        static bool outputKF;
        //        bool solvable;
        //       static bool useSchurComplementSolver;
        
        //        static bool use_directSolver;
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
//            EDR.readScalarInFile(fullName.str(),"outputDislocationStiffnessAndForce",DislocationNetworkComponentType::outputKF);
            outputKF=TextFileParser(fileName).readScalar<int>("outputDislocationStiffnessAndForce",true);            
        }
        
        /************************************************************/
        DislocationNetworkComponent(NetworkComponentType& NCin) :
        /* init list */ NC(NCin)
        //        /* init list */ solvable(true)
        {/*! Constructor initializes the reterence the NetworkComponent
          */
        }
        
        const NetworkComponentType& networkComponent() const
        {
            return NC;
        }
        
        /************************************************************/
        void directSolve()
        {
            TripletContainerType kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
            Eigen::VectorXd Fq; // the vector of nodal forces
            const size_t Ndof=assembleNCtriplets(kqqT,Fq);
            
            size_t newNdof=Ndof; // rows 0-(Ndof-1) are taken by Kqq. Start placing constraints at row Ndof
            assembleConstraints(kqqT,newNdof,true); // kqqT and KPQ_row are overwritten
            SparseMatrixType KQQ(newNdof,newNdof);
            
            KQQ.setFromTriplets(kqqT.begin(),kqqT.end()); // now KQQ is SPSD therefore conjugate-gradient cannot be used
            Eigen::VectorXd F(Eigen::VectorXd::Zero(newNdof));
            F.segment(0,Ndof)=Fq;
            
#ifdef _MODEL_PARDISO_SOLVER_
            Eigen::PardisoLDLT<SparseMatrixType>    solver(KQQ); // this fails
            //            Eigen::PardisoLU<SparseMatrixType>    solver(KQQ); // this fails
#else
            Eigen::SimplicialLDLT<SparseMatrixType> solver(KQQ);
#endif
            
            if(solver.info()==Eigen::Success)
            {
                Eigen::VectorXd x=solver.solve(F);
                storeNodeSolution(x.segment(0,Ndof));
                //                storeNodeSolution(F.segment(0,Ndof).cast<float>().cast<double>());
            }
            else
            {
                std::cout<<"DislocationNetworkComponent:" <<NC.sID<<" nodeOrder()="<<NC.nodeOrder()<<std::endl<<", linkOrder()="<<NC.linkOrder()<<std::endl;
                std::ofstream KQQfile("KQQfailed.txt");
//                KQQfile<<KQQ.toDense()<<std::endl;
                assert(0 && "RE-ENABLE OUTPUT");
                assert(0 && "LDLT DECOMPOSITION FAILED.");
            }
            
        }
        
        /**********************************************************************/
        void lumpedSolve(const size_t& runID)
        {
            if (NC.nodeOrder() > 0)
            {

                if (NC.nodeBegin()->second->network().simulationParameters.isPeriodicSimulation())
                {
//                    TripletContainerType kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
//                    Eigen::VectorXd Fq;        // the vector of nodal forces
//                    const size_t Ndof = assembleNCtriplets(kqqT, Fq);
//
//                    TripletContainerType zT;
//                    // std::cout<<"Node order size is "<<NC.nodeOrder()<<std::endl;
//                    size_t nUnconstrained = assembleConstraintsforPeriodicSimulationsNULL(zT);
//                    SparseMatrixType K(Ndof, Ndof);
//                    // std::cout<<"Ndof is "<<Ndof<<std::endl;
//                    K.setFromTriplets(kqqT.begin(), kqqT.end());
//                    SparseMatrixType Z(Ndof, dim * nUnconstrained);
//                    // std::cout<<"Setting up sparse matrix Z "<<std::endl;
//                    // std::cout<<"nUnconstrained is "<<nUnconstrained<<std::endl;
//                    // std::cout<<"Printing the triplets \n";
//                    // for (const auto& trip : zT)
//                    // {
//                    //     std::cout<<trip.row()<<", "<<trip.col()<<", "<<trip.value()<<std::endl;
//                    // }
//                    Z.setFromTriplets(zT.begin(), zT.end());
//                    // std::cout<<"Set up sparse matrix "<<std::endl;
//
//                    SparseMatrixType kqqZ(Z.transpose() * K * Z);
//                    // std::cout<<"KqqZ \n";
//                    // std::cout<<Eigen::MatrixXd(kqqZ)<<std::endl;
//
//                    // Zienkiewicz (See [1], section 16.2.4) discusses three methods for lumping the mass matrix
//                    TripletContainerType lumpedTriplets;
//                    for (int k = 0; k < kqqZ.outerSize(); ++k)
//                    {
//
//                        for (SparseMatrixType::InnerIterator it(kqqZ, k); it; ++it)
//                        {
//                            if (it.row() == it.col())
//                            {
//                                lumpedTriplets.emplace_back(it.row(), it.col(), it.value());
//                            }
//                            else
//                            {
//                                lumpedTriplets.emplace_back(it.row(), it.row(), 0.5 * it.value());
//                                lumpedTriplets.emplace_back(it.col(), it.col(), 0.5 * it.value());
//                            }
//                        }
//                    }
//                    SparseMatrixType kqq(nUnconstrained * dim, nUnconstrained * dim);
//                    kqq.setFromTriplets(lumpedTriplets.begin(), lumpedTriplets.end());
//                    Eigen::VectorXd Fqz(Z.transpose() * Fq);
//                    Eigen::VectorXd Kd(kqq.diagonal());
//                    // std::cout<<"Fqz "<<std::endl;
//                    // std::cout<<Fqz<<std::endl;
//                    // std::cout<<"Kd "<<std::endl;
//                    // std::cout<<Kd<<std::endl;
//
//                    Eigen::VectorXd x(Eigen::VectorXd::Zero(dim * nUnconstrained)); //Unconstrained Solution
//                    if (outputKF)
//                    {
//                        assert(0 && "RE-ENABLE OUTPUT");
//                        std::ofstream fileK("K_" + std::to_string(runID) + "_" + std::to_string(NC.sID) + ".txt");
//                        //                fileK<<Kd;
//                        std::ofstream fileF("F_" + std::to_string(runID) + "_" + std::to_string(NC.sID) + ".txt");
//                        //                fileF<<Fq;
//                    }
//                    // Check diagonal and force
//                    for (int k = 0; k < Kd.size(); ++k)
//                    {
//                        if (fabs(Kd(k)) > FLT_EPSILON)
//                        { // stiffness not zero
//                            x(k) = Fqz(k) / Kd(k);
//                        }
//                        else
//                        { // stiffness is zero
//                            if(fabs(Fqz(k)) > FLT_EPSILON)
//                            {
//                                std::cout<<"Kd(k)="<<Kd(k)<<std::endl;
//                                std::cout<<"Fqz(k)="<<Fqz(k)<<std::endl;
//                                assert(false && "if stiffness is zero also force must be zero.");
//                            }
//                        }
//                    }
//                    storeNodeSolution((Z * x).segment(0, Ndof));
                }
                else
                {
                    TripletContainerType kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
                    Eigen::VectorXd Fq;        // the vector of nodal forces
                    const size_t Ndof = assembleNCtriplets(kqqT, Fq);

                    // Zienkiewicz (See [1], section 16.2.4) discusses three methods for lumping the mass matrix
                    TripletContainerType lumpedTriplets;
                    for (const auto &t : kqqT)
                    {
                        if (t.col() == t.row())
                        {
                            lumpedTriplets.push_back(t);
                        }
                        else
                        {
                            lumpedTriplets.emplace_back(t.col(), t.col(), 0.5 * t.value());
                            lumpedTriplets.emplace_back(t.row(), t.row(), 0.5 * t.value());
                        }
                    }

                    SparseMatrixType kqq(Ndof, Ndof);
                    kqq.setFromTriplets(lumpedTriplets.begin(), lumpedTriplets.end());
                    Eigen::VectorXd Kd(kqq.diagonal());
                    // std::cout<<"Kd "<<std::endl;
                    // std::cout<<Kd<<std::endl;
                    // std::cout << "Fq " << std::endl;
                    // std::cout<<Fq<<std::endl;
                    Eigen::VectorXd x(Eigen::VectorXd::Zero(Ndof));

                    if (outputKF)
                    {
                        assert(0 && "RE-ENABLE OUTPUT");
                        std::ofstream fileK("K_" + std::to_string(runID) + "_" + std::to_string(NC.sID) + ".txt");
                        //                fileK<<Kd;
                        std::ofstream fileF("F_" + std::to_string(runID) + "_" + std::to_string(NC.sID) + ".txt");
                        //                fileF<<Fq;
                    }

                    // Check diagonal and force
                    for (int k = 0; k < Kd.size(); ++k)
                    {
                        if (fabs(Kd(k)) > FLT_EPSILON)
                        { // stiffness not zero
                            x(k) = Fq(k) / Kd(k);
                        }
                        else
                        { // stiffness is zero
                            if(fabs(Fq(k)) > FLT_EPSILON)
                            {
                                std::cout<<"k="<<k<<std::endl;
                                std::cout<<"k%NdofXnode="<<k%NdofXnode<<std::endl;
                                typename NodeContainerType::iterator nodeIter(NC.nodeBegin());
                                std::advance(nodeIter,k/NdofXnode);
                                std::cout<<"node "<<nodeIter->second->sID<<std::endl;
                                std::cout<<"Kd(k)="<<Kd(k)<<std::endl;
                                std::cout<<"Fq(k)="<<Fq(k)<<std::endl;


                                assert(false && "if stiffness is zero also force must be zero.");
                            }
//                            assert(fabs(Fq(k)) < FLT_EPSILON && "if stiffness is zero also force must be zero.");
                        }
                    }
                    storeNodeSolution(x.segment(0, Ndof));
                }
            }
        }

        /************************************************************/
        void iterativeSolve()
        {
            //            if(solvable)
            //            {
            if (NC.nodeBegin()->second->network().simulationParameters.isPeriodicSimulation())
            {
                //Currently the method implemented gives error in the velocity calculations using the penalty parameter
                assert(0 && "Periodic Simulations not supported with CG integration");
            }

            
            TripletContainerType kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
            Eigen::VectorXd Fq; // the vector of nodal forces
            const size_t Ndof=assembleNCtriplets(kqqT,Fq);
            
            // Find guess solution using current velocities
            Eigen::VectorXd x0(Ndof);
            size_t k=0;
            for (typename NodeContainerType::iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter){
                x0.segment(NdofXnode*k,NdofXnode)=nodeIter->second->get_V();
                ++k;
            }
            
            bool usePreSolver=true;
            if (NC.nodeBegin()->second->network().simulationParameters.isPeriodicSimulation())
            {
                assert(0 && "FINISH HERE");
                // if (usePreSolver)
                // {
                //     // kqq is SPD, find the unconstrained solution with ConjugateGradient (faster)
                //     size_t p_param=100000;
                //     assembleConstraintsforPeriodicSimulations(kqqT, p_param);
                //     SparseMatrixType kqq(Ndof, Ndof);
                //     kqq.setFromTriplets(kqqT.begin(), kqqT.end());
                //     Eigen::ConjugateGradient<SparseMatrixType> cg(kqq);
                //     if (cg.maxIterations() < 10000)
                //     {
                //         cg.setMaxIterations(10000);
                //     }
                //     x0 = cg.solveWithGuess(Fq, x0);
                // }
            }
            else
            {
                if (usePreSolver)
                {
                    // kqq is SPD, find the unconstrained solution with ConjugateGradient (faster)
                    SparseMatrixType kqq(Ndof, Ndof);
                    kqq.setFromTriplets(kqqT.begin(), kqqT.end());
                    Eigen::ConjugateGradient<SparseMatrixType> cg(kqq);

                    if (cg.maxIterations() < 10000)
                    {
                        cg.setMaxIterations(10000);
                    }
                    x0 = cg.solveWithGuess(Fq, x0);
                }
            }
            
            
            
            
            // Assemble constraints
            size_t KPQ_row=Ndof; // rows 0-(Ndof-1) are taken by Kqq. Start placing constraints at row Ndof
            assembleConstraints(kqqT,KPQ_row,true); // kqqT and KPQ_row are overwritten
            SparseMatrixType KQQ(KPQ_row,KPQ_row);
            KQQ.setFromTriplets(kqqT.begin(),kqqT.end()); // now KQQ is idefinite (NOT SPD) therefore conjugate-gradient cannot be used
            
            Eigen::VectorXd F(Eigen::VectorXd::Zero(KPQ_row));
            F.segment(0,Ndof)=Fq;
            Eigen::VectorXd x1(Eigen::VectorXd::Zero(KPQ_row));
            x1.segment(0,Ndof)=x0;
            
            bool useMINRES=true;
            if (useMINRES)
            {
                Eigen::MINRES<SparseMatrixType> solver(KQQ);
                solver.setTolerance(FLT_EPSILON);
                solver.setMaxIterations(10*F.rows());
                x1=solver.solveWithGuess(F,x1);
                if(solver.info()==Eigen::Success)
                {
                    storeNodeSolution(x1.segment(0,Ndof));
                }
                else
                {
                    std::cout<<"Solver did not converge\n";
                    std::cout<<"iterations="<<solver.iterations()<<"\n";
                    std::cout<<"maxIterations="<<solver.maxIterations()<<"\n";
                    std::cout<<"error="<<solver.error()<<"\n";
                    std::cout<<"tolerance="<<solver.tolerance()<<"\n";
                    assert(0 && "SOLVER DID NOT CONVERGE.");
                }
            }
            else
            {
                Eigen::BiCGSTAB<SparseMatrixType,Eigen::IdentityPreconditioner> bcs(KQQ);
                if(bcs.maxIterations()<10000)
                {
                    bcs.setMaxIterations(10000);
                }
                x1=bcs.solveWithGuess(F,x1);
                if(bcs.info()==Eigen::Success)
                {
                    storeNodeSolution(x1.segment(0,Ndof));
                    
                }
                else
                {
                    std::cout<<"Solver did not converge\n";
                    std::cout<<"iterations="<<bcs.iterations()<<"\n";
                    std::cout<<"maxIterations="<<bcs.maxIterations()<<"\n";
                    std::cout<<"error="<<bcs.error()<<"\n";
                    std::cout<<"tolerance="<<bcs.tolerance()<<"\n";
                    assert(0 && "SOLVER DID NOT CONVERGE.");
                }
            }
            //            }
            
            
        }
        
        //        /* isBoundarySubNetwork *********************************/
        //        bool isBoundarySubNetwork() const
        //        {
        //            bool temp(true);
        //            for (typename LinkContainerType::const_iterator linkIter=NC.linkBegin();linkIter!=NC.linkEnd();++linkIter)
        //            {
        //                temp*=linkIter->second->is_boundarySegment();
        //                if(!temp)
        //                {
        //                    break;
        //                }
        //            }
        //            return temp;
        //        }
        //
        //        /************************************************************/
        //        bool isLoop() const
        //        {
        //            bool temp=true;
        //            for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
        //            {
        //                temp*=nodeIter->second->is_simple();
        //            }
        //            return temp;
        //        }
        //
        //        /************************************************************/
        //        bool isPlanar() const
        //        {
        //            bool temp=true;
        //            VectorDim normal=NC.linkBegin()->second->glidePlaneNormal ;
        //            for (typename LinkContainerType::const_iterator linkIter=NC.linkBegin();linkIter!=NC.linkEnd();++linkIter)
        //            {
        //                temp*=(normal.cross(linkIter->second->glidePlaneNormal ).squaredNorm()<=FLT_EPSILON);
        //            }
        //            return temp;
        //        }
        
        //        /************************************************************/
        //        //		int findLinePoints(std::vector<VectorDim>& posVector) const {
        //        int findLinePoints(std::vector<typename NodeContainerType::const_iterator>& posVector) const
        //        {
        //            posVector.clear();
        //
        //            for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
        //            {
        //                if(nodeIter->second->constraintNormals().size()==2)
        //                { // this 2 should be dim-1 (moves on line)
        //                    //					posVector.push_back(nodeIter->second->get_P());
        //                    posVector.push_back(nodeIter);
        //
        //                }
        //            }
        //
        //            if (posVector.size()==0 && isPlanar() && isLoop())
        //            {
        //                typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();
        //                //				posVector.push_back(nodeIter->second->get_P());
        //                posVector.push_back(nodeIter);
        //                ++nodeIter;
        //                posVector.push_back(nodeIter);
        //                //				posVector.push_back(nodeIter->second->get_P());
        //            }
        //
        //            return posVector.size();
        //        }
        
        //        /************************************************************/
        //        bool loopInversion(const double& dt) const
        //        {
        //            bool temp(false);
        //            std::vector<typename NodeContainerType::const_iterator> posVector;
        //
        //            if(findLinePoints(posVector)==2){
        //                VectorDim L((posVector[1]->second->get_P()-posVector[0]->second->get_P()).normalized());
        //
        //                VectorDim A(VectorDim::Zero());
        //                VectorDim Anew(VectorDim::Zero());
        //                for (typename NodeContainerType::const_iterator nodeIter1=NC.nodeBegin();nodeIter1!=NC.nodeEnd();++nodeIter1){
        //                    if (nodeIter1->second->outOrder()==1){
        //                        NodeType* nextNode=std::get<0>(nodeIter1->second->outNeighborhood().begin()->second);
        //                        VectorDim Q1(nodeIter1->second->get_P());
        //                        VectorDim Q2(nextNode->get_P());
        //
        //                        VectorDim B((Q2-Q1).dot(L)*L);
        //                        VectorDim H((MatrixDim::Identity()-L*L.transpose())*(Q1+Q2));
        //                        A+=B.cross(H);
        //
        //                        // increment positions
        //                        Q1+=nodeIter1->second->get_V()*dt;
        //                        Q2+=nextNode->get_V()*dt;
        //                        B=(Q2-Q1).dot(L)*L;
        //                        H=(MatrixDim::Identity()-L*L.transpose())*(Q1+Q2);
        //                        Anew+=B.cross(H);
        //                    }
        //                }
        //
        //            }
        //            return temp;
        //        }
        
        //        /************************************************************/
        //        bool isSmall(const double& smallcritvalue,
        //                     const size_t& maxNodeSize) const
        //        {
        //            bool temp=false;
        //            const double critvelocity=0.02;
        //            if (NC.nodeOrder()<=maxNodeSize)
        //            {
        //                if( isLoop())
        //                {
        //                    VectorDim midpoint(VectorDim::Zero());
        //                    double velocityloop=0.0;
        //                    /////////////////to make sure that the loop is small is enough, distance to the middle point is small
        //                    for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
        //                    {
        //                        midpoint+=nodeIter->second->get_P();
        //                        double velocitynode=nodeIter->second->get_V().norm();
        //                        if (velocitynode>velocityloop)
        //                        {
        //                            velocityloop=velocitynode;
        //                        };
        //                    }
        //                    midpoint/=NC.nodeOrder();
        //                    double maxdis=0.0;
        //                    for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
        //                    {
        //                        double distance_to_mid=(nodeIter->second->get_P()-midpoint).norm();
        //                        if (distance_to_mid>maxdis)
        //                        {
        //                            maxdis=distance_to_mid;
        //                        }
        //                    }
        //
        //                    if(maxdis<smallcritvalue && velocityloop>critvelocity)
        //                    {
        //                        temp=true;
        //                    }
        //                }
        //
        //
        //            }
        //            return temp;
        //        }
        
    };
    
    template <typename NodeType,typename LinkType>
    bool DislocationNetworkComponent<NodeType,LinkType>::outputKF=false;
    
    //    //Static data
    //    template <typename NodeType,typename LinkType>
    //    bool DislocationNetworkComponent<NodeType,LinkType>::use_directSolver=true;
    
} // namespace model
#endif

