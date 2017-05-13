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
#include <model/Network/NetworkComponent.h>

//#include <model/Math/SchurComplementSolver.h>
#include <Eigen/Sparse>
//#include <model/Math/MINRES.h>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/SparseCholesky> // simplicial DLDT

//#include <model/Math/SparseNullSpace.h>

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

            size_t n=0;
            for (const auto& node : NC.nodes())
            {
                const typename NodeType::VectorOfNormalsType vn=node.second->constraintNormals();
                for (size_t c=0;c<vn.size();++c)
                {
                    for(size_t d=0;d<dim;++d)
                    {
//                        vT.push_back(Eigen::Triplet<double>(KPQ_row,n*dim+d,vn[c](d)));
                        vT.emplace_back(KPQ_row,n*dim+d,vn[c](d));
                        if (symmetricConstraint)
                        {
//                            vT.push_back(Eigen::Triplet<double>(n*dim+d,KPQ_row,vn[c](d)));
                              vT.emplace_back(n*dim+d,KPQ_row,vn[c](d));
                        }
                    }
                    ++KPQ_row; // move to next line
                }
                n++;
            }
            
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
        
        /************************************************************/
        size_t assembleNCtriplets(TripletContainerType& kqqT, Eigen::VectorXd& Fq)
        {
            const size_t nodeOrdr(NC.nodeOrder());
            if (nodeOrdr==0 || nodeOrdr==1)
            {
                std::cout<<"DislocationNetworkComponent:" <<NC.sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<NC.linkOrder()<<std::endl;
                typename NodeContainerType::const_iterator nodeIter1=NC.nodeBegin();
                std::cout<<"Only DIslocatioNode is "<< nodeIter1->second->sID<<std::endl;
                assert(0 && "DislocationSubNetework has less than 2 Nodes.");
            }
            if (nodeOrdr==2)
            {
                typename NodeContainerType::const_iterator nodeIter1=NC.nodeBegin();
                typename NodeContainerType::const_iterator nodeIter2(nodeIter1);
                nodeIter2++;
                if (nodeIter1->second->constraintNormals().size()>2 && nodeIter2->second->constraintNormals().size()>2)
                {
                    //					assert(0 && "DislocationSubNetework has only 2 fixed Nodes.");
                    std::cout<<"DislocationNetworkComponent:" <<NC.sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<NC.linkOrder()<<std::endl;
                    std::cout<<"First DIslocatioNode is " <<nodeIter1->second->sID<<std::endl;
                    std::cout<<"Second DIslocatioNode is "<<nodeIter2->second->sID<<std::endl;
                    std::cout<<"WARNING: DislocationSubNetework has only 2 fixed Nodes"<<std::endl;
                    
                }
            }
            
            
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
                linkIter->second->addToGlobalAssembly(kqqT,Fq); // loop over each segment and add segment contributions to kqqT and Fq
            }
            
            return Ndof;
            
        }
        
        
    public:
        //		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        enum {NdofXnode=NodeType::NdofXnode};
        //        bool solvable;
        //       static bool useSchurComplementSolver;
        
        static bool use_directSolver;
        
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
            Eigen::PardisoLDLT<SparseMatrixType>    solver(KQQ);
#else
            Eigen::SimplicialLDLT<SparseMatrixType> solver(KQQ);
#endif

            if(solver.info()==Eigen::Success)
            {
                Eigen::VectorXd x=solver.solve(F);
                storeNodeSolution(x.segment(0,Ndof));
            }
            else
            {
                std::cout<<"DislocationNetworkComponent:" <<NC.sID<<" nodeOrder()="<<NC.nodeOrder()<<std::endl<<", linkOrder()="<<NC.linkOrder()<<std::endl;
                std::ofstream KQQfile("KQQfailed.txt");
                KQQfile<<KQQ.toDense()<<std::endl;
                assert(0 && "LDLT DECOMPOSITION FAILED.");
            }

        }
        
        
        /************************************************************/
        void iterativeSolve()
        {
            //            if(solvable)
            //            {
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
            if (usePreSolver)
            {
                // kqq is SPD, find the unconstrained solution with ConjugateGradient (faster)
                SparseMatrixType kqq(Ndof,Ndof);
                kqq.setFromTriplets(kqqT.begin(),kqqT.end());
                Eigen::ConjugateGradient<SparseMatrixType> cg(kqq);
                
                if(cg.maxIterations()<10000)
                {
                    cg.setMaxIterations(10000);
                }
                x0=cg.solveWithGuess(Fq,x0);
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
        
        /* isBoundarySubNetwork *********************************/
        bool isBoundarySubNetwork() const
        {
            bool temp(true);
            for (typename LinkContainerType::const_iterator linkIter=NC.linkBegin();linkIter!=NC.linkEnd();++linkIter)
            {
                temp*=linkIter->second->is_boundarySegment();
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /************************************************************/
        bool isLoop() const
        {
            bool temp=true;
            for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
            {
                temp*=nodeIter->second->is_simple();
            }
            return temp;
        }
        
        /************************************************************/
        bool isPlanar() const
        {
            bool temp=true;
            VectorDim normal=NC.linkBegin()->second->glidePlaneNormal ;
            for (typename LinkContainerType::const_iterator linkIter=NC.linkBegin();linkIter!=NC.linkEnd();++linkIter)
            {
                temp*=(normal.cross(linkIter->second->glidePlaneNormal ).squaredNorm()<=FLT_EPSILON);
            }
            return temp;
        }
        
        /************************************************************/
        //		int findLinePoints(std::vector<VectorDim>& posVector) const {
        int findLinePoints(std::vector<typename NodeContainerType::const_iterator>& posVector) const
        {
            posVector.clear();
            
            for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
            {
                if(nodeIter->second->constraintNormals().size()==2)
                { // this 2 should be dim-1 (moves on line)
                    //					posVector.push_back(nodeIter->second->get_P());
                    posVector.push_back(nodeIter);
                    
                }
            }
            
            if (posVector.size()==0 && isPlanar() && isLoop())
            {
                typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();
                //				posVector.push_back(nodeIter->second->get_P());
                posVector.push_back(nodeIter);
                ++nodeIter;
                posVector.push_back(nodeIter);
                //				posVector.push_back(nodeIter->second->get_P());
            }
            
            return posVector.size();
        }
        
        /************************************************************/
        bool loopInversion(const double& dt) const
        {
            bool temp(false);
            std::vector<typename NodeContainerType::const_iterator> posVector;
            
            if(findLinePoints(posVector)==2){
                VectorDim L((posVector[1]->second->get_P()-posVector[0]->second->get_P()).normalized());
                
                VectorDim A(VectorDim::Zero());
                VectorDim Anew(VectorDim::Zero());
                for (typename NodeContainerType::const_iterator nodeIter1=NC.nodeBegin();nodeIter1!=NC.nodeEnd();++nodeIter1){
                    if (nodeIter1->second->outOrder()==1){
                        NodeType* nextNode=std::get<0>(nodeIter1->second->outNeighborhood().begin()->second);
                        VectorDim Q1(nodeIter1->second->get_P());
                        VectorDim Q2(nextNode->get_P());
                        
                        VectorDim B((Q2-Q1).dot(L)*L);
                        VectorDim H((MatrixDim::Identity()-L*L.transpose())*(Q1+Q2));
                        A+=B.cross(H);
                        
                        // increment positions
                        Q1+=nodeIter1->second->get_V()*dt;
                        Q2+=nextNode->get_V()*dt;
                        B=(Q2-Q1).dot(L)*L;
                        H=(MatrixDim::Identity()-L*L.transpose())*(Q1+Q2);
                        Anew+=B.cross(H);
                    }
                }
                
            }
            return temp;
        }
        
        /************************************************************/
        bool isSmall(const double& smallcritvalue,
                     const size_t& maxNodeSize) const
        {
            bool temp=false;
            const double critvelocity=0.02;
            if (NC.nodeOrder()<=maxNodeSize)
            {
                if( isLoop())
                {
                    VectorDim midpoint(VectorDim::Zero());
                    double velocityloop=0.0;
                    /////////////////to make sure that the loop is small is enough, distance to the middle point is small
                    for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
                    {
                        midpoint+=nodeIter->second->get_P();
                        double velocitynode=nodeIter->second->get_V().norm();
                        if (velocitynode>velocityloop)
                        {
                            velocityloop=velocitynode;
                        };
                    }
                    midpoint/=NC.nodeOrder();
                    double maxdis=0.0;
                    for (typename NodeContainerType::const_iterator nodeIter=NC.nodeBegin();nodeIter!=NC.nodeEnd();++nodeIter)
                    {
                        double distance_to_mid=(nodeIter->second->get_P()-midpoint).norm();
                        if (distance_to_mid>maxdis)
                        {
                            maxdis=distance_to_mid;
                        }
                    }                  
                    
                    if(maxdis<smallcritvalue && velocityloop>critvelocity)
                    {
                        temp=true;
                    }
                }
                
                
            }
            return temp;
        }
        
    };
    
    //Static data
    template <typename NodeType,typename LinkType>
    bool DislocationNetworkComponent<NodeType,LinkType>::use_directSolver=true;
    
} // namespace model
#endif

//        /************************************************************/
//        void nullSpaceSolve()
//        {
//            // Make some initial checks
//            TripletContainerType aT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
//            Eigen::VectorXd Fq; // the vector of nodal forces
//            const size_t Ndof=assembleNCtriplets(aT,Fq);
//            SparseMatrixType A(Ndof,Ndof);
//            A.setFromTriplets(aT.begin(),aT.end()); // now KQQ is idefinite (NOT SPD) therefore conjugate-gradient cannot be used
//
//            // testing null-space solver
//            TripletContainerType cT; // the vector of Eigen::Triplets of constraints
//            size_t rows=0; // start placing constraints at row=0
//            assembleConstraints<false>(cT,rows);
//            SparseMatrixType C(rows,Ndof);
//            std::cout<<"rows="<<rows<<std::endl;
//            std::cout<<"cols="<<Ndof<<std::endl;
//
//            int tempR=0;
//            int tempC=0;
//            for(const auto& t : cT)
//            {
//            if(t.row()>tempR)
//            {
//                tempR=t.row();
//            }
//                if(t.col()>tempC)
//                {
//                    tempC=t.col();
//                }
//            }
//            std::cout<<"cT_max_row="<<tempR<<std::endl;
//            std::cout<<"cT_max_col="<<tempC<<std::endl;
//
//
//            C.setFromTriplets(cT.begin(),cT.end());
//            SparseNullSpace<SparseMatrixType> ns(C);
//        }
