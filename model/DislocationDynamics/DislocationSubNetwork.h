/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DISLOCATIONSUBNETWORK_H_
#define model_DISLOCATIONSUBNETWORK_H_

#include <vector>
#include <float.h>
#include <Eigen/Dense>
#include <model/Network/SubNetwork.h>

#include <model/Math/SchurComplementSolver.h>
#include <model/Math/SchurComplementSparseSolver.h>
#include <model/Math/MINRES.h>

//#include <model/Geometry/Splines/SplineNetworkTraits.h>
//#include <model/Geometry/Splines/SplinesBase/SplineSubNetworkBase.h>

#include <model/Utilities/SequentialOutputFile.h>
#include <Eigen/Sparse>


namespace model {
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    //	class DislocationSubNetwork : public SplineSubNetworkBase<DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>,
    //	/*	                                                   */ dim,corder,alpha,InterpolationType> {		
	class DislocationSubNetwork : public SubNetwork<DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> >{		
		
		typedef DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> Derived;
#include <model/Network/NetworkTypedefs.h>
		
        //		typedef SplineSubNetworkBase<Derived,dim,corder,alpha,InterpolationType> SubNetworkBaseType;
		typedef SubNetwork<Derived> SubNetworkBaseType;
		
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		
        
        typedef Eigen::SparseMatrix<double> SparseMatrixType;

        
        /**********************************************************************/
        void storeNodeSolution(const Eigen::VectorXd& X){
            size_t k=0;
			for (typename SubNetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				nodeIter->second->set_V(X.segment(NdofXnode*k,NdofXnode));
				++k;
			}
        }
        
        
        /************************************************************/
        template <bool symmetricConstraint>
        void assembleConstraints(std::vector<Eigen::Triplet<double> >& vT, size_t& KPQ_row) const {
        
            // Collect contraints
            std::vector<typename NodeType::VectorOfNormalsType> NNV; // DON'T NEED THIS ANYMORE WITH TRIPLETS
			for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				NNV.push_back(nodeIter->second->constraintNormals());
			}
            
            // Constrain simple nodes to move normal to tangent
//            for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
//                if(nodeIter->second->constraintNormals().size()==1){
//                    Eigen::Matrix<int,dim,1> node_dofID(nodeIter->second->node_dofID());
//                    Eigen::Matrix<double,dim,1> T(nodeIter->second->get_T());
//
//                    for(size_t d=0;d<dim;++d){
//                        vT.push_back(Eigen::Triplet<double>(KPQ_row,node_dofID(d),T(d)));
//                    }
//                    ++KPQ_row;
//                }
//			}
            
            
            // loop over each segment and add segment contributions to kqqT and Fq
            //            for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
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
            
            for (size_t n=0 ;n<NNV.size();++n){
                for (size_t c=0;c<NNV[n].size();++c){
                    for(size_t d=0;d<dim;++d){
                        vT.push_back(Eigen::Triplet<double>(KPQ_row,n*dim+d,NNV[n][c](d)));
                        if (symmetricConstraint){
                            vT.push_back(Eigen::Triplet<double>(n*dim+d,KPQ_row,NNV[n][c](d)));
                        }
                    }
                    ++KPQ_row;
                }
            }
        
        }
        
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
		enum {NdofXnode=NodeType::NdofXnode};
        
        static bool useSchurComplementSolver;

		
		/************************************************************/
		/* Constructor with pointer to Node *************************/
		DislocationSubNetwork(NodeType* const pN) : SubNetworkBaseType::SubNetwork(pN){}
        //		DislocationSubNetwork(NodeType* const pN) : SubNetworkBaseType::SplineSubNetworkBase(pN){}
		
		/************************************************************/
		/* Constructor with pointer to Link *************************/
		DislocationSubNetwork(LinkType* const pL) : SubNetworkBaseType::SubNetwork(pL){}
        //		DislocationSubNetwork(LinkType* const pL) : SubNetworkBaseType::SplineSubNetworkBase(pL){}
		
		
        
        
        
        
        /************************************************************/
        void sparseSolve(){
            
            // Make some initial checks
            const size_t nodeOrdr(this->nodeOrder());
			if (nodeOrdr==0 || nodeOrdr==1){
                std::cout<<"DislocationSubNetwork:" <<this->sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<this->linkOrder()<<std::endl;
                assert(0 && "DislocationSubNetework has less than 2 Nodes.");
			}
			if (nodeOrdr==2){
				typename SubNetworkNodeContainerType::const_iterator nodeIter1=this->nodeBegin();
				typename SubNetworkNodeContainerType::const_iterator nodeIter2(nodeIter1);
				nodeIter2++;
				if (nodeIter1->second->constraintNormals().size()>2 && nodeIter2->second->constraintNormals().size()>2){
					assert(0 && "DislocationSubNetework has only 2 fixed Nodes.");
				}					
			}
            
            
            // Assembly of Stiffness Matrix and force vector
            const size_t Ndof(this->nodeOrder()*NdofXnode); // the total number of dof in the subnetwork
            size_t reserveSize(0); 
            for (typename SubNetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
                reserveSize+=nodeIter->second->closedOrder();
            }
            std::vector<Eigen::Triplet<double> > kqqT; // the vector of Eigen::Triplets corresponding to the matrix Kqq
            kqqT.reserve(reserveSize); // assume some fill percentage
            Eigen::VectorXd Fq(Eigen::VectorXd::Zero(Ndof));
            for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
                linkIter->second->addToGlobalAssembly(kqqT,Fq); // loop over each segment and add segment contributions to kqqT and Fq
			}
            
            if(useSchurComplementSolver){ // use SchurComplementSolver
                // Assemble Kqq and Kpq separately
                // Kqq
                SparseMatrixType KQQ(Ndof,Ndof);
                KQQ.setFromTriplets(kqqT.begin(),kqqT.end());
                // Kpq
                std::vector<Eigen::Triplet<double> > kpqT;
                size_t KPQ_row=0;
                assembleConstraints<false>(kpqT,KPQ_row);
                SparseMatrixType KPQ(KPQ_row,Ndof);
                KPQ.setFromTriplets(kpqT.begin(),kpqT.end());
                
                if (kqqT.size()<1){ // use direct solver (LLT decomposition) if the number of non-zeros in kqqT is less than a threshold
                    SchurComplementSparseSolver<SparseMatrixType,true> scs(KQQ);
                    scs.solve(KPQ,Fq);
                    storeNodeSolution(scs.X());
                }
                else{ // use iterative solver (Conjugate Gradient)
                    SchurComplementSparseSolver<SparseMatrixType,false> scs(KQQ);
                    scs.solve(KPQ,Fq,FLT_EPSILON*0.001); // FLT_EPSILON is the tolerance
                    storeNodeSolution(scs.X);
                }
                
            }
            else{ // use MINRES solver
                // Assemble Kqq and Kpq together
                size_t KPQ_row=Ndof; // start placing constraints at row Ndof
                assembleConstraints<true>(kqqT,KPQ_row); // kqqT and KPQ_row are overwritten
                SparseMatrixType KQQ(KPQ_row,KPQ_row);
                KQQ.setFromTriplets(kqqT.begin(),kqqT.end()); // now KQQ is idefinite (NOT SPD) therefore conjugate-gradient cannot be used
                Eigen::VectorXd F(Eigen::VectorXd::Zero(KPQ_row));
                F.segment(0,Ndof)=Fq;
                Eigen::VectorXd x0(Eigen::VectorXd::Zero(KPQ_row));
                MINRES<double> mRS(KQQ,F,x0,DBL_EPSILON*100.0);                
                storeNodeSolution(mRS.xMR.segment(0,Ndof));
            }
            
            

            
            
        }
        
        
        
        
		/************************************************************/
		void solve() {
            
            
            
            //#ifdef _OPENMP
            //            std::cout<<"Thread "<<omp_get_thread_num()<<" of "<<omp_get_num_threads()<<" solving SubNetwork "<<this->sID<<std::endl;
            //#endif	
            
			
			const size_t nodeOrdr(this->nodeOrder());
			
			if (nodeOrdr==0 || nodeOrdr==1){
                std::cout<<"DislocationSubNetwork:" <<this->sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<this->linkOrder()<<std::endl;
                assert(0 && "DislocationSubNetework has less than 2 Nodes.");
			}
			
			if (nodeOrdr==2){
				typename SubNetworkNodeContainerType::const_iterator nodeIter1=this->nodeBegin();
				typename SubNetworkNodeContainerType::const_iterator nodeIter2(nodeIter1);
				nodeIter2++;
				if (nodeIter1->second->constraintNormals().size()>2 && nodeIter2->second->constraintNormals().size()>2){
					assert(0 && "DislocationSubNetework has only 2 fixed Nodes.");
				}					
			}
            
            
			
            
			
			
			
			//! 1- Assemble KQQ and Fq
			const size_t NDOF = this->nodeOrder()*NdofXnode;
			
			
			
			Eigen::MatrixXd KQQ = Eigen::MatrixXd::Zero(NDOF,NDOF);
			Eigen::VectorXd Fq  = Eigen::VectorXd::Zero(NDOF);
			
			
			for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
                
                //		linkIter->second->assemble(); // moved to DislocationNetwork::solve
                
				Eigen::MatrixXd G2H(linkIter->second->get_G2H());
				KQQ+=G2H.transpose() * linkIter->second->get_Kqq() * G2H;   
				Fq +=G2H.transpose() * linkIter->second->get_Fq(); 
			}
			
			
			
			
			//! 1- Assemble KPQ and Fp
			std::vector<typename NodeType::VectorOfNormalsType> NNV;
			
			for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				NNV.push_back(nodeIter->second->constraintNormals());
			}
			
			size_t count=0;
			for (size_t k=0;k<NNV.size();++k){
				count+=NNV[k].size();
			}
			
			//assert(count<KQQ.rows() && "SUBNETWORK IS FULLY CONSTRAINED");
			
			Eigen::MatrixXd KPQ = Eigen::MatrixXd::Zero(count,KQQ.cols());
            //			Eigen::VectorXd Fp  = Eigen::VectorXd::Zero(count);
			
			size_t KPQ_row=0;
			
			for (size_t n=0 ;n<NNV.size();++n){	
				for (size_t c=0;c<NNV[n].size();++c){
					KPQ.block(KPQ_row,n*dim,1,dim)=NNV[n][c].transpose();
                    //					Fp(KPQ_row)=0.0;					
					++KPQ_row;
				}
			}
			
			
			
			//! Solve the system of equations
            //			const SchurComplementSolver LS(KQQ,KPQ,Fq,Fp);
			const SchurComplementSolver LS(KQQ,KPQ,Fq);
			
			//! Store velocities in DislocationNodes
            size_t k=0;
			for (typename SubNetworkNodeContainerType::iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				nodeIter->second->set_V(LS.X.segment(NdofXnode*k,NdofXnode));                
				++k;
			}
            
            
			
		}
		
		/* isBoundarySubNetwork *********************************/
		bool isBoundarySubNetwork() const {
			bool temp(true);
			for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp*=linkIter->second->is_boundarySegment();
				if(!temp){
					break;
				}
			}
			return temp;
		}
		
		/************************************************************/
		bool isLoop() const {
			bool temp=true;
			for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				temp*=nodeIter->second->is_simple();
			}
			return temp;
		}
        		
		/************************************************************/
		bool isPlanar() const {
			bool temp=true;
			VectorDim normal=this->linkBegin()->second->glidePlaneNormal ;
			for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
				temp*=(normal.cross(linkIter->second->glidePlaneNormal ).squaredNorm()<=FLT_EPSILON);
			}
			return temp;
		}
		
		/************************************************************/
        //		int findLinePoints(std::vector<VectorDim>& posVector) const {
		int findLinePoints(std::vector<typename SubNetworkNodeContainerType::const_iterator>& posVector) const {
			posVector.clear();
			
			for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
				if(nodeIter->second->constraintNormals().size()==2){ // this 2 should be dim-1 (moves on line)
                    //					posVector.push_back(nodeIter->second->get_P());
					posVector.push_back(nodeIter);
                    
				}
			}
			
			if (posVector.size()==0 && isPlanar() && isLoop()){
				typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();
                //				posVector.push_back(nodeIter->second->get_P());
				posVector.push_back(nodeIter);
				++nodeIter;
				posVector.push_back(nodeIter);
                //				posVector.push_back(nodeIter->second->get_P());
			}
			
			return posVector.size();
		}
		
		/************************************************************/
		bool loopInversion(const double& dt) const {
			bool temp(false);			
			std::vector<typename SubNetworkNodeContainerType::const_iterator> posVector;
			
			if(findLinePoints(posVector)==2){
				VectorDim L((posVector[1]->second->get_P()-posVector[0]->second->get_P()).normalized());
				
				VectorDim A(VectorDim::Zero());
				VectorDim Anew(VectorDim::Zero());
				for (typename SubNetworkNodeContainerType::const_iterator nodeIter1=this->nodeBegin();nodeIter1!=this->nodeEnd();++nodeIter1){
					if (nodeIter1->second->outOrder()==1){
						NodeType* nextNode=boost::tuples::get<0>(nodeIter1->second->outNeighborhood().begin()->second);						
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
		
	};
	
	// static data
    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	bool DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>::useSchurComplementSolver=false;

	
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif






//            // Collect contraints
//            std::vector<typename NodeType::VectorOfNormalsType> NNV; // DON'T NEED THIS ANYMORE WITH TRIPLETS
//			for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
//				NNV.push_back(nodeIter->second->constraintNormals());
//			}


//            // Constrain sinple nodes to move normal to tangent
//            for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
//                if(nodeIter->second->constraintNormals().size()==1){
//                    Eigen::Matrix<int,dim,1> node_dofID(nodeIter->second->node_dofID());
//                    Eigen::Matrix<double,dim,1> T(nodeIter->second->get_T());
//
//                    for(size_t d=0;d<dim;++d){
//                        kpqT.push_back(Eigen::Triplet<double>(KPQ_row,node_dofID(d),T(d)));
//                    }
//                    ++KPQ_row;
//                }
//			}


// loop over each segment and add segment contributions to kqqT and Fq
//            for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
//                const Eigen::VectorXd ortCS(linkIter->second->ortC*linkIter->second->Mseg);
//                for (unsigned int i=0;i<linkIter->second->segmentDOFs.size();++i){
//                    std::set<size_t>::const_iterator iterI(linkIter->second->segmentDOFs.begin());
//                    std::advance(iterI,i);
//                    //FQ(*iterI)+=tempFq(i);
//                    kpqT.push_back(Eigen::Triplet<double>(KPQ_row,*iterI,ortCS(i)));
//                }
//                ++KPQ_row;
//			}



// Create a sparse matrix from the triplets
//            SparseMatrixType KQQ(Ndof,Ndof);
//            KQQ.setFromTriplets(kqqT.begin(),kqqT.end());



//            std::vector<Eigen::Triplet<double> > kpqT;
//
//            size_t KPQ_row=0;
//
//            for (size_t n=0 ;n<NNV.size();++n){
//                for (size_t c=0;c<NNV[n].size();++c){
//                    for(size_t d=0;d<dim;++d){
//                        kpqT.push_back(Eigen::Triplet<double>(KPQ_row,n*dim+d,NNV[n][c](d)));
//                    }
//                    ++KPQ_row;
//                }
//            }
//
//            SparseMatrixType KPQ(KPQ_row,Ndof);
//            KPQ.setFromTriplets(kpqT.begin(),kpqT.end());
