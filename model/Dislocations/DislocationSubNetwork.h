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
//#include <model/Math/SchurComplementSparseSolver.h>

//#include <model/Geometry/Splines/SplineNetworkTraits.h>
//#include <model/Geometry/Splines/SplinesBase/SplineSubNetworkBase.h>

namespace model {
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
//	class DislocationSubNetwork : public SplineSubNetworkBase<DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>,
//	/*	                                                   */ dim,corder,alpha,InterpolationType> {		
	class DislocationSubNetwork : public SubNetwork<DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> >{		
		
		typedef DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType> Derived;
#include <model/Network/NetworkTypedefs.h>
		
//		typedef SplineSubNetworkBase<Derived,dim,corder,alpha,InterpolationType> SubNetworkBaseType;
		typedef SubNetwork<Derived> SubNetworkBaseType;
		
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		enum {NdofXnode=NodeType::NdofXnode};
		
		/************************************************************/
		/* Constructor with pointer to Node *************************/
		DislocationSubNetwork(NodeType* const pN) : SubNetworkBaseType::SubNetwork(pN){}
//		DislocationSubNetwork(NodeType* const pN) : SubNetworkBaseType::SplineSubNetworkBase(pN){}
		
		/************************************************************/
		/* Constructor with pointer to Link *************************/
		DislocationSubNetwork(LinkType* const pL) : SubNetworkBaseType::SubNetwork(pL){}
//		DislocationSubNetwork(LinkType* const pL) : SubNetworkBaseType::SplineSubNetworkBase(pL){}
		
		
		/************************************************************/
		/* Constructor with pointer to Link *************************/
		void solve() {
			
			size_t nodeOrdr(this->nodeOrder());
			
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
			
//			switch (nodeOrdr) {
//				case 0:
//					assert(0 && "DislocationSubNetework has 0 Nodes.");
//					break;
//				case 1:
//					assert(0 && "DislocationSubNetework has 1 Nodes.");
//					break;
//				case 2:
//					break;
//				default:
//					break;
//			}
			
//			if(nodeOrdr<=2){
//				std::cout<<"DislocationSubNetwork:" <<this->sID<<" nodeOrder()="<<nodeOrdr<<std::endl<<", linkOrder()="<<this->linkOrder()<<std::endl;
//				assert(0 && "DislocationSubNetework has less than 3 Nodes.");
//
//			}

			

			
			
			
			//! 1- Assemble KQQ and Fq
			size_t NDOF = this->nodeOrder()*NdofXnode;
			
			
			
			Eigen::MatrixXd KQQ = Eigen::MatrixXd::Zero(NDOF,NDOF);
			Eigen::VectorXd Fq  = Eigen::VectorXd::Zero(NDOF);
			
			
			for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){

		//		linkIter->second->assemble(); // moved to DislocationNetwork::solve

				Eigen::MatrixXd G2H(linkIter->second->get_G2H());
				KQQ+=G2H.transpose() * linkIter->second->get_Kqq() * G2H;   
				Fq+= G2H.transpose() * linkIter->second->get_Fq(); 
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
			Eigen::VectorXd Fp  = Eigen::VectorXd::Zero(count);
			
			size_t KPQ_row=0;
			
			for (size_t n=0 ;n<NNV.size();++n){	
				for (size_t c=0;c<NNV[n].size();++c){
					KPQ.block(KPQ_row,n*dim,1,dim)=NNV[n][c].transpose();
					Fp(KPQ_row)=0.0;					
					++KPQ_row;
				}
			}
			
			
			
			//! Solve the system of equations
			SchurComplementSolver LS(KQQ,KPQ,Fq,Fp);
			
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
	
	
	
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif


//		/********************************************************/
//		VectorDim displacement(const VectorDim & Rfield) const {
//			// compute the average source point
//			VectorDim rs(VectorDim::Zero());
////			double snLength(0.0);
////			for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
////				double linkLength(linkIter->second->template arcLength<qOrder,GaussLegendre>());
////				snLength+=linkLength;
////				rs+=linkIter->second->template rm<qOrder,GaussLegendre>()*linkLength;
////			}
////			rs/=snLength;
//
//			for (typename SubNetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
//				//double linkLength(linkIter->second->template arcLength<qOrder,GaussLegendre>());
////				snLength+=linkLength;
//				rs+=nodeIter->second->get_P();
//			}
//			rs/=this->nodeOrder();
//			
//			
//			// compute the vector S
//			VectorDim S((Rfield-rs).normalized());
//
//			// sum the displacements from DislocationSegments using  S			
//			VectorDim temp(VectorDim::Zero());
//			for (typename SubNetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
//				temp+=linkIter->second->displacement(Rfield,S);
//			}
//			
//			return temp;
//		}
