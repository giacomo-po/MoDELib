/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationGlideSolver_cpp_
#define model_DislocationGlideSolver_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationGlideSolver.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    DislocationGlideSolver<DislocationNetworkType>::DislocationGlideSolver(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    {
        
    }
    
    template <typename DislocationNetworkType>
    void DislocationGlideSolver<DislocationNetworkType>::solve(const size_t& runID)
    {
        lumpedSolve(runID);
    }
    
    
    
    template <typename DislocationNetworkType>
    size_t DislocationGlideSolver<DislocationNetworkType>::assembleNCtriplets(TripletContainerType& kqqT, Eigen::VectorXd& Fq)
    {
        const size_t Ndof(DN.networkNodes().size()*NdofXnode); // the total number of dof in the subnetwork
        kqqT.clear();
        Fq.setZero(Ndof);
        
        for (const auto& networkLink : DN.networkLinks())
        {
            if (!networkLink.second.lock()->isBoundarySegment())
            {
                networkLink.second.lock()->addToGlobalAssembly(kqqT,Fq); // loop over each segment and add segment contributions to kqqT and Fq
            }
        }
        
        return Ndof;
        
    }
    
    template <typename DislocationNetworkType>
    void DislocationGlideSolver<DislocationNetworkType>::storeNodeSolution(const Eigen::VectorXd& X)
    {
        size_t k=0;
        for (auto& networkNode : DN.networkNodes())
        {
        //    networkNode.second.lock()->set_V(X.segment(NdofXnode*k,NdofXnode).template cast<float>().template cast<double>()); // double cast to remove some numerical noise
            // std::cout<<"setting v for "<<networkNode.second.lock()->sID<<" "<<(X.segment(NdofXnode*k,NdofXnode)).transpose()<<std::endl;
            networkNode.second.lock()->set_V(X.segment(NdofXnode*k,NdofXnode)); // double cast to remove some numerical noise
            ++k;
        }
    }

//     template <typename DislocationNetworkType>
//     size_t DislocationGlideSolver<DislocationNetworkType>::assembleConstraintsforPeriodicSimulationsNULL(TripletContainerType &zT) const
//     {
//         // std::map<const NetworkNodeType *const, std::set<std::tuple<const NetworkNodeType *const,const double, const NetworkNodeType *const,const double>>> networkNodeContainer;
//         std::map<const NetworkNodeType *const, std::map<std::pair<const NetworkNodeType *const, const NetworkNodeType *const>, std::pair<const double, const double>>> networkNodeContainer;

//         size_t constrainedI = 0;
//         size_t unconstrainedNodes = 0;
//         std::map<size_t, size_t> correctedJPosition;

//         for (const auto &netNode : DN.networkNodes())
//         {
//             std::map<const LoopType *const, std::pair<std::pair<const LoopNodeType *const, const double>, std::pair<const LoopNodeType *const, const double>>> loopConnectivity;
//             if (netNode.second.lock()->isBoundaryNode())
//             {
//                 for (const auto &loopNode : netNode.second.lock()->loopNodes())
//                 {
//                     if (loopNode->periodicPrev() && loopNode->periodicNext())
//                     {
//                         const auto periodicPrev(loopNode->periodicPrev());
//                         const auto periodicNext(loopNode->periodicNext());
//                         const double lij = (loopNode->get_P() - periodicPrev->get_P()).norm();
//                         const double ljk = (loopNode->get_P() - periodicNext->get_P()).norm();

//                         loopConnectivity.emplace(loopNode->loop().get(), std::make_pair(std::make_pair(periodicPrev, lij),
//                                                                                         std::make_pair(periodicNext, ljk)));
//                     }
//                 }
                
//                 // std::cout<<"Current network node is "<<netNode.second.lock()->sID<<" [ "<<netNode.second.lock()->isBoundaryNode()<<" ]"<<std::endl;
//                 // std::cout<<"Loop Connectivity size "<<loopConnectivity.size()<<std::endl;
//                 assert((loopConnectivity.size() == 0 ||loopConnectivity.size() == 1|| loopConnectivity.size() == netNode.second.lock()->uniqueLoopNodes()) && "Junction Node at Boundary Not Defined Properly"); //Discuss this with Dr. Po
// // loopConnectivity.size() == 1 for the case of self annihilation
//                 if (loopConnectivity.size() == 0)
//                 {
//                     const size_t ntempsnID(netNode.second.lock()->networkID()); //Global position in the constraint matrix (j)
//                     correctedJPosition.emplace(ntempsnID, ntempsnID - constrainedI);
//                 }
//                 else
//                 {
//                     constrainedI++;
//                     std::map<std::pair<const NetworkNodeType *const, const NetworkNodeType *const>, std::pair<const double, const double>> innerMap;

//                     for (const auto &loopConn : loopConnectivity)
//                     {
//                         const auto nodeI(std::min(loopConn.second.first.first->networkNode->sID, loopConn.second.second.first->networkNode->sID) == loopConn.second.first.first->networkNode->sID ? std::make_pair(loopConn.second.first.first->networkNode.get(), loopConn.second.first.second) : std::make_pair(loopConn.second.second.first->networkNode.get(), loopConn.second.second.second));

//                         const auto nodeJ(std::max(loopConn.second.first.first->networkNode->sID, loopConn.second.second.first->networkNode->sID) == loopConn.second.first.first->networkNode->sID ? std::make_pair(loopConn.second.first.first->networkNode.get(), loopConn.second.first.second) : std::make_pair(loopConn.second.second.first->networkNode.get(), loopConn.second.second.second));

//                         assert(nodeI.first != nodeJ.first && "I and J network nodes cannot be same");
//                         innerMap.emplace(std::make_pair(nodeI.first, nodeJ.first), std::make_pair(nodeI.second, nodeJ.second));
//                     }
//                     networkNodeContainer.emplace(netNode.second.lock().get(), innerMap);
//                 }
//             }
//             else
//             {
//                 const size_t ntempsnID(netNode.second.lock()->networkID()); //Global position in the constraint matrix (j)
//                 correctedJPosition.emplace(ntempsnID, ntempsnID - constrainedI);
//             }
//         }

//         for (const auto &netNode : DN.networkNodes())
//         {
//             const size_t ntempj(netNode.second.lock()->networkID()); //Global position in the constraint matrix (j)
//             if (netNode.second.lock()->isBoundaryNode())
//             {
//                 const auto netNodeIter(networkNodeContainer.find(netNode.second.lock().get()));
//                 if (netNodeIter != networkNodeContainer.end())
//                 {
//                     //Constraints exist in the network node container
//                     for (const auto &innerMapIter : netNodeIter->second)
//                     {
//                         const double lij(innerMapIter.second.first );
//                         const double ljk(innerMapIter.second.second);
//                         const double lijk(lij + ljk);
//                         size_t ntempi(innerMapIter.first.first->networkID());  //Global position in the constraint matrix corresponding to the original link (i)
//                         size_t ntempk(innerMapIter.first.second->networkID()); //Global position in the constraint matrix corresponding to the neighbor link (k)

//                         assert(correctedJPosition.find(ntempi) != correctedJPosition.end());
//                         assert(correctedJPosition.find(ntempk) != correctedJPosition.end());

//                         const size_t correctedI(correctedJPosition.find(ntempi)->second);
//                         const size_t correctedK(correctedJPosition.find(ntempk)->second);

//                         for (size_t d = 0; d < dim; ++d)
//                         {
//                             // zT.emplace_back(dim*ntempj+d,dim*ntempi+d,ljk/lijk);

//                             zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
//                             zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
//                         }
//                     }
//                 }
//                 else
//                 {
//                     assert(correctedJPosition.find(netNode.second.lock()->networkID()) != correctedJPosition.end());

//                     const size_t correctedJ(correctedJPosition.find(netNode.second.lock()->networkID())->second);
//                     for (size_t d = 0; d < dim; ++d)
//                     {
//                         zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
//                     }
//                     unconstrainedNodes++;
//                 }
//             }
//             else
//             {
//                 //Constraints do not exist in the network node container
//                 assert(correctedJPosition.find(netNode.second.lock()->networkID()) != correctedJPosition.end());

//                 const size_t correctedJ(correctedJPosition.find(netNode.second.lock()->networkID())->second);
//                 for (size_t d = 0; d < dim; ++d)
//                 {
//                     zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1);
//                 }
//                 unconstrainedNodes++;
//             }
//         }
//         return unconstrainedNodes;
//     }
//Improved Implementation for the constrained nodes
    template <typename DislocationNetworkType>
    size_t DislocationGlideSolver<DislocationNetworkType>::assembleConstraintsforPeriodicSimulationsNULL(TripletContainerType &zT) const
    {
        // std::map<const NetworkNodeType *const, std::set<std::tuple<const NetworkNodeType *const,const double, const NetworkNodeType *const,const double>>> networkNodeContainer;
        std::map<const NetworkNodeType *const, std::map<std::pair<const NetworkNodeType *const, const NetworkNodeType *const>, std::pair<const double, const double>>> networkNodeContainer;

        size_t constrainedI = 0;
        size_t unconstrainedNodes = 0;
        std::map<size_t, size_t> correctedJPosition;

        for (const auto &netNode : DN.networkNodes())
        {
            std::set<std::pair<std::pair<const LoopNodeType *const, const double>, std::pair<const LoopNodeType *const, const double>>> loopConnectivity;

            if (netNode.second.lock()->isBoundaryNode() )
            {
                for (const auto &loopNode : netNode.second.lock()->loopNodes())
                {
                    
                    if (loopNode->periodicPrev() && loopNode->periodicNext())
                    {

                        const auto periodicPrev(loopNode->periodicPrev());
                        const auto periodicNext(loopNode->periodicNext());

                        const double lij = (loopNode->get_P() - periodicPrev->get_P()).norm();
                        const double ljk = (loopNode->get_P() - periodicNext->get_P()).norm();

                        loopConnectivity.emplace(std::make_pair(std::make_pair(periodicPrev, lij),
                                                               std::make_pair(periodicNext, ljk)));
                    }
                }
                
                // std::cout<<"Current network node is "<<netNode.second.lock()->sID<<" [ "<<netNode.second.lock()->isBoundaryNode()<<" ]"<<std::endl;
                // std::cout<<"Loop Connectivity size "<<loopConnectivity.size()<<std::endl;
                assert((loopConnectivity.size() == 0 ||loopConnectivity.size() == 1|| loopConnectivity.size() == netNode.second.lock()->loopNodes().size()) && "Junction Node at Boundary Not Defined Properly"); //Discuss this with Dr. Po
// loopConnectivity.size() == 1 for the case of self annihilation
                if (loopConnectivity.size() == 0)
                {
                    // std::cout<<"Coming in for "<<netNode.second.lock()->tag()<<std::endl;
                    const size_t ntempsnID(netNode.second.lock()->networkID()); //Global position in the constraint matrix (j)
                    correctedJPosition.emplace(ntempsnID, ntempsnID - constrainedI);
                }
                else
                {
                    constrainedI++;
                    std::map<std::pair<const NetworkNodeType *const, const NetworkNodeType *const>, std::pair<const double, const double>> innerMap;

                    for (const auto &loopConn : loopConnectivity)
                    {
                        // std::cout<<"Temp network IDs"<< loopConn.first.first->networkNode->sID <<"=>"<<loopConn.second.first->networkNode->sID<<std::endl;

                        const auto nodeI(std::min(loopConn.first.first->networkNode->sID, loopConn.second.first->networkNode->sID) == loopConn.first.first->networkNode->sID ? 
                        std::make_pair(loopConn.first.first->networkNode.get(), loopConn.first.second) : std::make_pair(loopConn.second.first->networkNode.get(), loopConn.second.second));

                        const auto nodeJ(std::max(loopConn.first.first->networkNode->sID, loopConn.second.first->networkNode->sID) == loopConn.first.first->networkNode->sID ? 
                        std::make_pair(loopConn.first.first->networkNode.get(), loopConn.first.second) : std::make_pair(loopConn.second.first->networkNode.get(), loopConn.second.second));

                        // if (nodeI.first == nodeJ.first)
                        // {
                        //     std::cout<<"I and J network Node" <<nodeI.first->sID<<"--->"<<nodeJ.first->sID<<std::endl;
                        // }
                        // assert(nodeI.first != nodeJ.first && "I and J network nodes cannot be same"); //They can be same (PeriodicPrev and PeriodicNext can point to same node)
                        
                        innerMap.emplace(std::make_pair(nodeI.first, nodeJ.first), std::make_pair(nodeI.second, nodeJ.second));
                    }
                    if (innerMap.size()!=1)
                    {
                        std::cout<<" For node "<<netNode.second.lock()->tag()<<" inner map size is "<<innerMap.size()<<std::endl;
                        for (const auto& inMap : innerMap)
                        {
                            std::cout<<inMap.first.first->tag()<<"=>"<<inMap.first.second->tag()<<std::endl;
                        }
                    }
                    assert(innerMap.size()==1 && "Unique constraint cannot be assigned");
                    networkNodeContainer.emplace(netNode.second.lock().get(), innerMap);
                }
            }
            else
            {
                const size_t ntempsnID(netNode.second.lock()->networkID()); //Global position in the constraint matrix (j)
                correctedJPosition.emplace(ntempsnID, ntempsnID - constrainedI);
            }
        }
        
        for (const auto &netNode : DN.networkNodes())
        {
            const size_t ntempj(netNode.second.lock()->networkID()); //Global position in the constraint matrix (j)
            if (netNode.second.lock()->isBoundaryNode()  )
            {
                const auto netNodeIter(networkNodeContainer.find(netNode.second.lock().get()));
                if (netNodeIter != networkNodeContainer.end())
                {
                    //Constraints exist in the network node container
                    for (const auto &innerMapIter : netNodeIter->second)
                    {
                        const double lij(innerMapIter.second.first );
                        const double ljk(innerMapIter.second.second);
                        const double lijk(lij + ljk);
                        size_t ntempi(innerMapIter.first.first->networkID());  //Global position in the constraint matrix corresponding to the original link (i)
                        size_t ntempk(innerMapIter.first.second->networkID()); //Global position in the constraint matrix corresponding to the neighbor link (k)

                        assert(correctedJPosition.find(ntempi) != correctedJPosition.end());
                        assert(correctedJPosition.find(ntempk) != correctedJPosition.end());

                        const size_t correctedI(correctedJPosition.find(ntempi)->second);
                        const size_t correctedK(correctedJPosition.find(ntempk)->second);

                        for (size_t d = 0; d < dim; ++d)
                        {
                            // zT.emplace_back(dim*ntempj+d,dim*ntempi+d,ljk/lijk);

                            zT.emplace_back(dim * ntempj + d, dim * correctedI + d, ljk / lijk);
                            zT.emplace_back(dim * ntempj + d, dim * correctedK + d, lij / lijk);
                        }
                    }
                }
                else
                {
                    assert(correctedJPosition.find(netNode.second.lock()->networkID()) != correctedJPosition.end());

                    const size_t correctedJ(correctedJPosition.find(netNode.second.lock()->networkID())->second);
                    for (size_t d = 0; d < dim; ++d)
                    {
                        zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1.0);
                    }
                    unconstrainedNodes++;
                }
            }
            else
            {
                //Constraints do not exist in the network node container
                assert(correctedJPosition.find(netNode.second.lock()->networkID()) != correctedJPosition.end());

                const size_t correctedJ(correctedJPosition.find(netNode.second.lock()->networkID())->second);
                for (size_t d = 0; d < dim; ++d)
                {
                    zT.emplace_back(dim * ntempj + d, dim * correctedJ + d, 1.0);
                }
                unconstrainedNodes++;
            }
        }
        
        if((constrainedI-(DN.networkNodes().size()-unconstrainedNodes))!=0)
        {
            throw std::runtime_error("Constraining nodes failed.");
        }
        
//        assert((constrainedI-(DN.networkNodes().size()-unconstrainedNodes))==0);
        return unconstrainedNodes;
    }

    template <typename DislocationNetworkType>
    void DislocationGlideSolver<DislocationNetworkType>::lumpedSolve(const size_t& )
    {
//        if (NC.nodeOrder() > 0)
//        {
//
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"GlideSolver (lumped) "<<std::flush;
        TripletContainerType kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
        Eigen::VectorXd Fq;        // the vector of nodal forces
        const size_t Ndof = assembleNCtriplets(kqqT, Fq);

        
            if (DN.simulationParameters.isPeriodicSimulation())
            {
                TripletContainerType zT;

               const size_t nUnconstrained = assembleConstraintsforPeriodicSimulationsNULL(zT);
                
                SparseMatrixType K(Ndof, Ndof);
                K.setFromTriplets(kqqT.begin(), kqqT.end());
                SparseMatrixType Z(Ndof, dim * nUnconstrained);
                Z.setFromTriplets(zT.begin(), zT.end());

                SparseMatrixType kqqZ(Z.transpose() * K * Z);

                // Zienkiewicz (See [1], section 16.2.4) discusses three methods for lumping the mass matrix
                TripletContainerType lumpedTriplets;
                for (int k = 0; k < kqqZ.outerSize(); ++k)
                {

                    for (SparseMatrixType::InnerIterator it(kqqZ, k); it; ++it)
                    {
                        if (it.row() == it.col())
                        {
                            lumpedTriplets.emplace_back(it.row(), it.col(), it.value());
                        }
                        else
                        {
                            lumpedTriplets.emplace_back(it.row(), it.row(), 0.5 * it.value());
                            lumpedTriplets.emplace_back(it.col(), it.col(), 0.5 * it.value());
                        }
                    }
                }
                                
                SparseMatrixType kqq(nUnconstrained * dim, nUnconstrained * dim);
                kqq.setFromTriplets(lumpedTriplets.begin(), lumpedTriplets.end());
                Eigen::VectorXd Fqz(Z.transpose() * Fq);
                Eigen::VectorXd Kd(kqq.diagonal());

                Eigen::VectorXd x(Eigen::VectorXd::Zero(dim * nUnconstrained)); //Unconstrained Solution
                // Check diagonal and force
                for (int k = 0; k < Kd.size(); ++k)
                {
                    if (fabs(Kd(k)) > FLT_EPSILON)
                    { // stiffness not zero
                        x(k) = Fqz(k) / Kd(k);
                    }
                    else
                    { // stiffness is zero
                        if(fabs(Fqz(k)) > FLT_EPSILON)
                        {
                            std::cout<<"Kd(k)="<<Kd(k)<<std::endl;
                            std::cout<<"Fqz(k)="<<Fqz(k)<<std::endl;
                            assert(false && "if stiffness is zero also force must be zero.");
                        }
                    }
                }

                storeNodeSolution((Z * x).segment(0, Ndof));
            }
            else
            {
//
//                // Zienkiewicz (See [1], section 16.2.4) discusses three methods for lumping the mass matrix
//                TripletContainerType lumpedTriplets;
//                for (const auto &t : kqqT)
//                {
//                    if (t.col() == t.row())
//                    {
//                        lumpedTriplets.push_back(t);
//                    }
//                    else
//                    {
//                        lumpedTriplets.emplace_back(t.col(), t.col(), 0.5 * t.value());
//                        lumpedTriplets.emplace_back(t.row(), t.row(), 0.5 * t.value());
//                    }
//                }
//
//                SparseMatrixType kqq(Ndof, Ndof);
//                kqq.setFromTriplets(lumpedTriplets.begin(), lumpedTriplets.end());
//                Eigen::VectorXd Kd(kqq.diagonal());
//                // std::cout<<"Kd "<<std::endl;
//                // std::cout<<Kd<<std::endl;
//                // std::cout << "Fq " << std::endl;
//                // std::cout<<Fq<<std::endl;
//                Eigen::VectorXd x(Eigen::VectorXd::Zero(Ndof));
//
//                if (outputKF)
//                {
//                    assert(0 && "RE-ENABLE OUTPUT");
//                    std::ofstream fileK("K_" + std::to_string(runID) + "_" + std::to_string(NC.sID) + ".txt");
//                    //                fileK<<Kd;
//                    std::ofstream fileF("F_" + std::to_string(runID) + "_" + std::to_string(NC.sID) + ".txt");
//                    //                fileF<<Fq;
//                }
//
//                // Check diagonal and force
//                for (int k = 0; k < Kd.size(); ++k)
//                {
//                    if (fabs(Kd(k)) > FLT_EPSILON)
//                    { // stiffness not zero
//                        x(k) = Fq(k) / Kd(k);
//                    }
//                    else
//                    { // stiffness is zero
//                        if(fabs(Fq(k)) > FLT_EPSILON)
//                        {
//                            std::cout<<"k="<<k<<std::endl;
//                            std::cout<<"k%NdofXnode="<<k%NdofXnode<<std::endl;
//                            typename NodeContainerType::iterator nodeIter(NC.nodeBegin());
//                            std::advance(nodeIter,k/NdofXnode);
//                            std::cout<<"node "<<nodeIter->second->sID<<std::endl;
//                            std::cout<<"Kd(k)="<<Kd(k)<<std::endl;
//                            std::cout<<"Fq(k)="<<Fq(k)<<std::endl;
//
//
//                            assert(false && "if stiffness is zero also force must be zero.");
//                        }
//                        //                            assert(fabs(Fq(k)) < FLT_EPSILON && "if stiffness is zero also force must be zero.");
//                    }
//                }
//                storeNodeSolution(x.segment(0, Ndof));
            }
//        }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

    }
    template class DislocationGlideSolver<DislocationNetwork<3,0>>;
    
}
#endif
