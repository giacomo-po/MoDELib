/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONENERGYRULES_H_
#define model_DISLOCATIONENERGYRULES_H_

#include <boost/tuple/tuple.hpp>
#include<Eigen/Dense>


namespace model {
	
	template <unsigned int N>
	struct EdgeConfigs{
		
		static Eigen::Matrix<int,Pow<2,N-1>::value,N> get_Ci(){
			Eigen::Matrix<int,Pow<2,N-1>::value,N> temp;
			temp.template block<Pow<2,N-2>::value,N-1>(0,1)=EdgeConfigs<N-1>::get_Ci();
			temp.template block<Pow<2,N-2>::value,N-1>(Pow<2,N-2>::value,1)=-EdgeConfigs<N-1>::get_Ci();
			temp.col(0).setOnes();
			return temp;
		}
		
		//	static const Eigen::Matrix<int,Pow<2,N-1>::value,N> Ci;
	};
	
	//! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
	//template<unsigned int N>
	//const Eigen::Matrix<int,Pow<2,N-1>::value,N> EdgeConfigs<N>::Ci=EdgeConfigs<N>::get_Ci();
	
	
	template <>
	struct EdgeConfigs<1>{
		enum {N=1};
		
		static Eigen::Matrix<int,Pow<2,N-1>::value,N> get_Ci(){
			//		Eigen::Matrix<int,Pow<2,N-1>::value,N> temp;
			//		temp<<1;
			//		return temp;
			return Eigen::Matrix<int,Pow<2,N-1>::value,N>::Ones();
		}
		
		//	static const Eigen::Matrix<int,Pow<2,N-1>::value,N> Ci;
	};
	
	
	
	template<short unsigned int dim>
	class DislocationEnergyRules {
		
		typedef Eigen::Matrix<double,dim,1> VectorDimD;
		
	public:
		
		
		static double interactionEnergy(const VectorDimD& b1, const VectorDimD& b2, const int& s1, const int& s2, const VectorDimD& t, const double& nu){
			const VectorDimD t1(static_cast<double>(s1)*t);
			const VectorDimD t2(static_cast<double>(s2)*t);
			//	   dot(b1,t1)*dot(b2,t2)+1/(1-nu)    *(dot(b1,b2)*dot(t1,t2)-dot(b1,t2)*dot(b2,t1)); MATLAB CODE
			return b1.dot(t1)*b2.dot(t2)+1.0/(1.0-nu)*(b1.dot(b2)*t1.dot(t2)-b1.dot(t2)*b2.dot(t1));
		}

		
		
		/* getCi *****************************************************/
		static Eigen::MatrixXi getCi(const unsigned int& n){
			Eigen::MatrixXi Ci;
			switch (n) {
				case 1:
					Ci=EdgeConfigs<1>::get_Ci();
					break;
				case 2:
					Ci=EdgeConfigs<2>::get_Ci();
					break;
				case 3:
					Ci=EdgeConfigs<3>::get_Ci();
					break;
				case 4:
					Ci=EdgeConfigs<4>::get_Ci();
					break;
				case 5:
					Ci=EdgeConfigs<5>::get_Ci();
					break;
				case 6:
					Ci=EdgeConfigs<6>::get_Ci();
					break;
				case 7:
					Ci=EdgeConfigs<7>::get_Ci();
					break;
				case 8:
					Ci=EdgeConfigs<8>::get_Ci();
					break;
				case 9:
					Ci=EdgeConfigs<9>::get_Ci();
					break;
				case 10:
					Ci=EdgeConfigs<10>::get_Ci();
					break;
				case 11:
					Ci=EdgeConfigs<11>::get_Ci();
					break;
				case 12:
					Ci=EdgeConfigs<12>::get_Ci();
					break;
					
					
				default:

					break;
			}
			return Ci;
		}
		
		
		/* findEdgeConfiguration *********************************************/
		template <typename DislocationNodeType>
		static void findEdgeConfiguration(DislocationNodeType& dN, const double& nu){
			
			if (dN.is_isolated()){
				dN.edgeConfiguration.setZero(0);
			}
			else{
				// 1- Collect all the Burgers from the neighbors considering them as out-neighbors.
				std::vector<typename DislocationNodeType::FlowType>  BsV; // the vector of Burgers as out-neighbors
				for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter){
					int dir=boost::tuples::get<2>(neighborIter->second);
					switch ( dir ) {
						case   1:
							BsV.push_back(+1.0*(boost::tuples::get<1>(neighborIter->second)->flow));
							break;
						case  -1:
							BsV.push_back(-1.0*(boost::tuples::get<1>(neighborIter->second)->flow));
							break;
						default:	// self
							break;
					}
				}
				
				// 2- Get all the possible edge configurations				
				Eigen::MatrixXi Ci(getCi(BsV.size()));
				if (Ci.cols()!=BsV.size()){
//					std::cout<<"Node "<<dN.sID<<": #edges="<<BsV.size()<<std::endl;
					for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter){
						int dir=boost::tuples::get<2>(neighborIter->second);
						switch ( dir ) {
							case   1:
	//							std::cout<<dN.sID<<"->"<<(boost::tuples::get<0>(neighborIter->second)->sID)<<"\n";
								//<<std::end;
								break;
							case  -1:
//								std::cout<<(boost::tuples::get<0>(neighborIter->second)->sID)<<"->"<<dN.sID<<"\n";			
								break;
							default:	// self
								break;
						}
					}
					assert(0 && "NEED TO IMPLEMENT MULTIPLE JUNCTIONS");				
				}
				
				// 3- Change the sign of Ci according to the direction of the first neighbor
				int sigCi=1;
				for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter){
					const int dir(boost::tuples::get<2>(neighborIter->second));
					if (dir!=0){
						sigCi=dir;
						break;
					}
				}
				Ci*=sigCi;
				
				// 4- Store total chord parametric length
				double gT(0.0);
				for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter){
					if (boost::tuples::get<2>(neighborIter->second)!=0){
						gT+=boost::tuples::get<1>(neighborIter->second)->chordParametricLength();
					}
				}
				
				int chi=1+!(dN.neighborhood().size()>2);
//				std::cout<<"chi="<<chi<<std::endl;
	//			std::cout<<"neigh.size="<<dN.neighborhood().size()<<std::endl;
				gT*=chi;			
				
				// 5- Compute and sort energy levels
				std::multimap<double,Eigen::VectorXi> ELev;	// MAP DOES NOT ACCEPT EQUAL VALUES SO multimap IS USED TO STORE  DEGENERATE STATES
				for (int k=0;k<Ci.rows();++k){
					
					// 5a- Compute the Catmull-Rom tangent using the current edge configuration
					VectorDimD tTemp(VectorDimD::Zero());
					int j(0);
					for (typename DislocationNodeType::constNeighborIteratorType neighborIter=dN.neighborhood().begin();neighborIter!=dN.neighborhood().end();++neighborIter){
						if (boost::tuples::get<2>(neighborIter->second)!=0){
							double gkj(boost::tuples::get<1>(neighborIter->second)->chordParametricLength());
//							std::cout<<Ci(k,j)*(boost::tuples::get<0>(neighborIter->second)->get_P()-dN.get_P()).transpose()/gkj*(gT-gkj)/gT<<std::endl;
//							std::cout<<boost::tuples::get<0>(neighborIter->second)->get_P().transpose()<<std::endl;
//							std::cout<<dN.get_P().transpose()<<std::endl;
							tTemp+=Ci(k,j)*(boost::tuples::get<0>(neighborIter->second)->get_P()-dN.get_P())/gkj*(gT-gkj)/gT;
							j++;
						}
					}

					const double tTempNorm(tTemp.norm());
//					assert(tTempNorm>FLT_EPSILON && "NORM OF NODE TANGENT MUST BE > 0");
//					tTemp/=tTempNorm;
					//tTemp=dN.prjM*tTemp; // this makes a dipolar loop crash
					
					// 5b- Compute the energy for the current configuration and current tangent
					double Ek(0.0);
					for(unsigned int a=0; a<BsV.size();++a){
						for(unsigned int b=0; b<BsV.size();++b){
							Ek+=interactionEnergy(BsV[a],BsV[b],Ci(k,a),Ci(k,b),tTemp,nu);
						}
					}
					ELev.insert(std::make_pair(Ek,Ci.row(k)));
				}
				
				
				// 5- reset edgeConfiguration
				dN.edgeConfiguration.setZero(Ci.cols());
				
				if (dN.is_balanced()){
					assert(ELev.size()>1 && "MORE THAN ONE ENERGY LEVELS MUST BE FOUND FOR A BALANCED NODE");	
					
					//std::map<double,Eigen::VectorXi>::const_iterator firstPositiveE=std::find_if(ELev.begin(),ELev.size.end(),isPositive)
					//				assert(ELev.begin()->first>=0.0 &&		  "FIRST ENERGY LEVEL MUST BE 0.");
					assert(std::fabs(ELev.begin()->first)<FLT_EPSILON && "FIRST ENERGY LEVEL MUST BE 0 FOR A BALANCED NODE.");
					
					std::multimap<double,Eigen::VectorXi>::const_iterator firstNonZero(ELev.lower_bound(FLT_EPSILON)); // the first element that compares >=FLT_EPSILON
					assert(firstNonZero!=ELev.end() && "AT LEAST ONE POSITIVE ENERGY LEVEL MUST EXIST");
					std::multimap<double,Eigen::VectorXi>::const_iterator nextNonZero(firstNonZero);
					++nextNonZero;
					if (nextNonZero==ELev.end()){	// MINIMUM ENERGY LEVEL IS LAST ONE SO IT IS UNIQUE
						dN.edgeConfiguration=firstNonZero->second;
					}
					else{
						if(std::fabs(nextNonZero->first-firstNonZero->first)>=FLT_EPSILON){ // MINIMUM ENERGY LEVEL IS UNIQUE
							dN.edgeConfiguration=firstNonZero->second;
						}
					}
					
					
					
					
				}
				
				
				
				
				// 7- store tangent coefficients
				unsigned int kk(0);
				typename DislocationNodeType::NeighborContainerType tempNeigh=dN.neighborhood(); // DESIGN FLAW: neighborhood should have a non-const version
				for (typename DislocationNodeType::NeighborIteratorType neighborIter=tempNeigh.begin();neighborIter!=tempNeigh.end();++neighborIter){
					int dir=boost::tuples::get<2>(neighborIter->second);
					switch ( dir ) {
						case   1:	// out
							boost::tuples::get<1>(neighborIter->second)->sourceTfactor=dN.edgeConfiguration(kk);
//							boost::tuples::get<1>(neighborIter->second)->sourceTfactor=sigCi*dN.edgeConfiguration(kk);
							
							kk++;
							break;
						case  -1:	// in
							boost::tuples::get<1>(neighborIter->second)->  sinkTfactor=dN.edgeConfiguration(kk);
//														boost::tuples::get<1>(neighborIter->second)->  sinkTfactor=sigCi*dN.edgeConfiguration(kk);
							kk++;
							break;
						default:	// self
							break;
					}
				}
				
				
				
				
				
				
			}
			
			
			
			
			
		}
		
		
		
		
	};
	
}
#endif






//				switch (BsV.size()) {
//					case 1:
//						Ci=EdgeConfigs<1>::get_Ci();
//						break;
//					case 2:
//						//					Ci.resize(Pow<2,N2-1>::value,N2);
//						Ci=EdgeConfigs<2>::get_Ci();
//						break;
//					case 3:
//						Ci=EdgeConfigs<3>::get_Ci();
//						break;
//					case 4:
//						Ci=EdgeConfigs<4>::get_Ci();
//						break;
//					case 5:
//						Ci=EdgeConfigs<5>::get_Ci();
//						break;
//					case 6:
//						Ci=EdgeConfigs<6>::get_Ci();
//						break;
//					case 7:
//						Ci=EdgeConfigs<7>::get_Ci();
//						break;
//					case 8:
//						Ci=EdgeConfigs<8>::get_Ci();
//						break;
//					case 9:
//						Ci=EdgeConfigs<9>::get_Ci();
//						break;
//					case 10:
//						Ci=EdgeConfigs<10>::get_Ci();
//						break;
//					case 11:
//						Ci=EdgeConfigs<11>::get_Ci();
//						break;
//					case 12:
//						Ci=EdgeConfigs<12>::get_Ci();
//						break;
//						
//						
//					default:
//						std::cout<<"Node "<<this->sID<<": #edges="<<BsV.size()<<std::endl;
//						for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
//							int dir=boost::tuples::get<2>(neighborIter->second);
//							switch ( dir ) {
//								case   1:
//									std::cout<<this->sID<<"->"<<(boost::tuples::get<0>(neighborIter->second)->sID)<<"\n";
//									//<<std::end;
//									break;
//								case  -1:
//									std::cout<<(boost::tuples::get<0>(neighborIter->second)->sID)<<"->"<<this->sID<<"\n";			
//									break;
//								default:	// self
//									break;
//							}
//						}
//						
//						assert(0 && "NEED TO IMPLEMENT MULTIPLE JUNCTIONS");
//						break;
//				}






//	bool isInsideMesh;		
//size_t constraintType;
//VectorDim constraintNormal;



//DislocationNode(const InitializationData<Derived>& idV) : NodeBaseType::SplineNodeBase(idV),
///*                                        */ constraintType(freeNode),
///*                                        */ constraintNormal(VectorDim::Zero()),


//velocity=this->get_nodeDof()*0.0;

//if (this->is_balanced()){

//}
//else{
//	constraintType=fixedNode;
//}

//			constraintNormal.setZero();
//			boundaryNormal.setZero();
//VectorDim projectedP;





//bvpfe::Tetrahedron::isTetrahedronType isT=shared.domain.findIncludingTet(this->get_P(),projectedP);
//			std::cout<<"OLD POSITION"<<this->P<<std::endl;
//		std::cout<<"NEW POSITION"<<newP<<std::endl;
//if (isT.first){

//}
//else{
//	currentMeshID=-1;
//	std::cout<< "now projectedP is "<<projectedP.transpose()<<std::endl;
//	this->P=projectedP; //! THIS SHOULD BE ENABLED BUT IT FAILS TO PUT THE NODE ON THE SLIP PLANE
//assert(0 && "NODE IS OUTSIDE DOMAIN");
//}


//		/////////////////////////////////////////////////////////////
//		// constraintNormals
//		std::vector<VectorDim> constraintNormals() const {
//			//	std::cout<<"I'm here 3"<<std::endl;
//			std::vector<VectorDim>  GS;
//			std::vector<VectorDim>  PN;
//			switch (constraintType) {
//				case fixedNode:
//					//			std::cout<<"I'm here 4"<<std::endl;
//					GS.push_back((VectorDim()<<1.0,0.0,0.0).finished());
//					GS.push_back((VectorDim()<<0.0,1.0,0.0).finished());
//					GS.push_back((VectorDim()<<0.0,0.0,1.0).finished());
//					break;
//				default:
//					//PN=planeNormals();
//					PN=planeNormals;
//					//PN.push_back(constraintNormal);
//					PN.push_back(boundaryNormal);
//					GramSchmidt<dim>( PN,  GS);
//					break;
//			}
//			assert(GS.size()>=1 && GS.size()<=3);
//			return GS;
//		}




//		/////////////////////////////////////////////////////////////
//		// planeNormals
//		std::vector<VectorDim> planeNormals() const{
//			std::vector<VectorDim> NV;
//			//LinkType* pL;
//		//	std::cout<<"Dislocation Node "<<this->sID<<std::endl;
//			
//			for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter){
//				
//		//		std::cout<<"Neighbors are:"<< (boost::tuples::get<0>(neighborIter->second))->sID<<std::endl;
//				
//				if (boost::tuples::get<2>(neighborIter->second)){
//					LinkType* pL=boost::tuples::get<1>(neighborIter->second); // use copy constructor here?
//					// insert only if glidePlaneNormal  and -glidePlaneNormal  are not found 
//					//					std::cout<<"Dislocation Node: I'm here 3"<<std::endl;
//					if (std::find(NV.begin(),NV.end(),pL->glidePlaneNormal )==NV.end() && std::find(NV.begin(),NV.end(),-pL->glidePlaneNormal )==NV.end()){
//						NV.push_back(pL->glidePlaneNormal );
//						//						std::cout<<"Dislocation Node: I'm here 4"<<std::endl;
//						
//					}
//				}				
//			}
//			
//			
//		//	std::cout<<"NV.size="<<NV.size()<<std::endl;
//		//	assert(NV.size()>0 && NV.size()<3);	// REWORK THIS FOR BCC, still the maximum on linear independent plane normals must be at max 2
//			return NV;
//		}



/////////////////////////////////////////////////////////////
// isBalanced // MARKED TO BE IN THE NETWORK LAYER
//		bool isBalanced() const {
//			return this->flowBalance().norm()==0.0;
//		}


//	template <short unsigned int dim>
//	class DislocationNodeData{
//		typedef Eigen::Matrix<double,dim,1> VectorDim;
//		
//	public:
//		const Eigen::VectorXd Q;
//		const size_t constraintType;
//		//const VectorDim constraintNormal;
//		//const VectorDim   boundaryNormal;
//		//const Eigen::VectorXd velocity;
//		
//		DislocationNodeData(const Eigen::VectorXd & Qin, const size_t & constraintTypeIn) : 
//		Q(Qin), constraintType(constraintTypeIn)
//		{}
//	
//	};





//		/////////////////////////////////////////////////////////////
//		// ????
//		DisconnectContainerType binaryBreakList() {
//			double tolNorm=1.0e-5;
//			
//			DisconnectContainerType DisconnectContainer;
//			
//			if (this->OutNeighborhood.size()>1 && this->InNeighborhood.size()>1){
//				
//				for (typename NeighborContainerType::iterator InIter=this->InNeighborhood.begin(); InIter!=this->InNeighborhood.end();++InIter){
//					for (typename NeighborContainerType::iterator OutIter=this->OutNeighborhood.begin(); OutIter!=this->OutNeighborhood.end();++OutIter){
//						if ( (boost::tuples::get<1>(OutIter->second)->flow()-boost::tuples::get<1>(InIter->second)->flow()).norm()<tolNorm /*&&
//																																			boost::tuples::get<1>(InIter ->second)->get_rll(1.0).dot(boost::tuples::get<1>(OutIter->second)->get_rll(0.0))>0.0*/ ){
//																																				DisconnectContainer.push_back(boost::tuples::make_tuple( boost::tuples::get<0>(InIter->second),
//																																																		this,
//																																																		boost::tuples::get<0>(OutIter->second) ));
//																																			}
//					}
//				}
//				
//				
//			}
//			
//			return DisconnectContainer;
//		}



/////////////////////////////////////////////////////////////
// GramSchmidt
//		void GramSchmidt(const std::vector<Eigen::Matrix<double,dim,1> > & NV, std::vector<Eigen::Matrix<double,dim,1> > & GS) const {
//			
//			GS=NV;
//			
//			for (size_t i=0;i<NV.size();++i){
//				
//				for (size_t j=0;j<i;++j){
//					GS[i]-= NV[i].dot(GS[j])*GS[j];
//				}
//				if (GS[i].squaredNorm()>0.0){
//					GS[i].normalize();
//				}
//				else{
//					GS[i].setZero();
//				}
//			}
//			
//			//			std::cout<<"Node "<<this->sID<< " size GS = "<<GS.size()<<std::endl;
//			
//			
//			for (typename std::vector<VectorDim>::iterator iter=GS.begin();iter!=GS.end();){
//				if (iter->squaredNorm()==0.0){
//					//					std::cout<<"I'm here 5"<<std::endl;
//					
//					iter=GS.erase(iter); // TEST THIS !!!
//					//					std::cout<<"I'm here 6"<<std::endl;
//					
//				}
//				else{
//					++iter;
//				}
//			}
//			
//			
//		}

