/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CRYSTAL_H_
#define model_CRYSTAL_H_

#include <set>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <Eigen/Core>
#include <model/Dislocations/Materials/SlipSystem.h>


namespace model {
	
	enum {FCC, FCCP};
	
	
	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// GENERAL CASE
	template<short unsigned int dim, short unsigned int Nslips>
	class CrystalBase : public std::set<model::SlipSystem<dim,Nslips> >{
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		typedef Eigen::Matrix<double,dim,Nslips> MatrixDimNslips;
		
		
	public:
		
		
		/////////////////////////////////////////////////////////////
		// find_slipSystem
		void find_slipSystem(const VectorDim  & chord, const VectorDim & Burgers,
							 std::set<SlipSystem<dim,Nslips> > & allowedSlipSystems, const double& tol = FLT_EPSILON){

			assert(chord.norm()>tol && "CHORD HAS ZERO NORM");
			assert(Burgers.norm()>tol && "BURGERS HAS ZERO NORM");
			
			
			VectorDim normalizedChord(chord.normalized());
	
			allowedSlipSystems.clear();
				
			for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=this->begin();iter!=this->end();++iter){
				if(	  std::fabs( iter->normal.dot(normalizedChord))<tol && std::fabs( iter->normal.dot(Burgers.normalized()))<tol){
					allowedSlipSystems.insert( *iter );
				}
			}
			
			unsigned int N(allowedSlipSystems.size());
			switch (N) {
				case 0: // CHECK FOR SESSILE
					for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=this->begin();iter!=this->end();++iter){
//						std::cout<<"c*n="<< std::fabs( iter->normal.dot(normalizedChord)) << " tol is "<<tol<<std::endl;
						if(	std::fabs( iter->normal.dot(normalizedChord))<tol ){
							allowedSlipSystems.insert( *iter );
						}
					}
					if (allowedSlipSystems.size()==0){
						std::cout<<" chord is"<<chord.transpose()<<std::endl;
						for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=this->begin();iter!=this->end();++iter){
							std::cout<<"n="<<iter->normal.transpose()<<" |c*n|="<< std::fabs( iter->normal.dot(normalizedChord)) << " tol is "<<tol<<std::endl;
						}
						assert(allowedSlipSystems.size() && "WRONG NUMBER OF SLIP SYSTEMS");					
					}
				//	std::cout<<"CRYSTAL.h: FOUND SESSILE SEGMENT"<<std::endl;
					break;
				case 1: // OK
					break;
				case 2: // CROSS-SLIP SEGMENT
					break;
				default:
					assert(0 && "NUMBER OF SLIP SYSTEMS>2");
					break;
			}

			
		}
		

//		/////////////////////////////////////////////////////////////
//		// find_slipSystem
//		void find_slipSystem(const VectorDim  & chord, const VectorDim & Burgers,
//							 /*const VectorDim  & T1, const VectorDim & T2,*/
//							 std::set<SlipSystem<dim,Nslips> > & allowedSlipSystems, const double& tol = FLT_EPSILON){
//			
//			
//			//;
//
//			
//			assert(chord.norm()>tol && "CHORD HAS ZERO NORM");
//			assert(Burgers.norm()>tol && "BURGERS HAS ZERO NORM");
//	//		assert(normalizedChord.cross(Burgers.normalized()).norm()>tol && "CHORD AND BURGERS ARE PARALLEL. SLIP PLANE CANNOT BE DETERMINED.");
//			
//			// TO DO: HERE IN ORDER TO FIND THE SLIP PLANE NORMAL USE THE CROSS
//			// OF CHORD AND BURGERS (unless is screw)!!
//			
//			//double tol=1.0e-10;
//			
//			
//			
//			
//			//VectorDim normalizedChord=chord/chord.norm();
//			//! Clear the content of planenormals
//			//			this->operator=(&FCCsss);
//			allowedSlipSystems.clear();
//			
//			//std::cout<<"normalized chord="<<normalizedChord<<std::endl;
//			
//			//! THIS DOES NOT ALLOW SESSILE SEGMENTS!!!!
//			for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=this->begin();iter!=this->end();++iter){
//				
////				std::cout<<"Checking vs n="<<iter->normal.transpose()<<std::endl;
////				std::cout<<"c*n="<< std::fabs( iter->normal.dot(normalizedChord)) << " tol is "<<tol<<std::endl;
////				std::cout<<"b*n="<< std::fabs( iter->normal.dot(Burgers.normalized())) << " tol is "<<tol<<std::endl;
//				
//				if(	  std::fabs( iter->normal.dot(normalizedChord))<tol 
//				   && std::fabs( iter->normal.dot(Burgers.normalized()))<tol
//				  /* && std::fabs( iter->normal.dot(T1))<tol
//				   && std::fabs( iter->normal.dot(T2))<tol*/
//				   ){
//					allowedSlipSystems.insert( *iter );
//				}
//			}
//			
//
//			
//			
//			if (allowedSlipSystems.size()<1 || allowedSlipSystems.size()>2){
//				/*std::cout<<"T1 is "<<T1.transpose()<<std::endl;
//				std::cout<<"T2 is "<<T2.transpose()<<std::endl;*/
//
//				std::cout<<"chord is "<<chord.transpose()<<std::endl;
//				std::cout<<"Burgers is "<< Burgers.transpose()<<std::endl; 
//				//	assert(allowedSlipSystems.size()>=1 && allowedSlipSystems.size()<=2);		// THE 2 IS ONLY FOR FCC, BOR BCC IS 3
//				assert(0 && "WRONG NUMBER OF SLIP SYSTEMS");
//			}
//			
//		}
		
		
		
	};
	
	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// GENERAL CASE
	template<short unsigned int dim, short unsigned int CrystalType>
	class Crystal {
		
	public:
		Crystal(){
			std::cout<<"Crystal Type not implemented"<<std::endl;
			assert(0);
		}
		
	};
	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// FCC CRYSTAL
	template<short unsigned int dim>
	class Crystal<dim,FCC> : public CrystalBase<dim,3> {
		
		
	public:
		
		
		enum{Nslips=3};
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		typedef Eigen::Matrix<double,dim,Nslips> MatrixDimNslips;
		
		
		Eigen::Matrix<double,dim,dim> C2G;
		
		//! The number of slip directions per slip plane
		
		
		//#include <model/Dislocations/Materials/Crystal_common.h>
		
		
		
		// THIS APPROACH IS WRONG BECAUSE IN FCC DIFFERENT PLANES CONTAIN
		// A DIFFERENT NUMBER OF SLIP DIRECTIONS SO NSLIPS IS NOT WELL DEFINED
		
		void rotate(const Eigen::Matrix<double,dim,dim>& C2G_in){
			
			
			
			
			// make sure that C2G is orthogonal
			assert((C2G_in*C2G_in.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && 
				   "CRYSTAL TO GLOBAL ROTATION MATRIX IS NOT ORTHOGONAL.");
			
			
			C2G=C2G_in;
			
			//! 0 clear existing slip planes
			this->clear();
			
			// 1- Define the 4 slip normals and directions
			VectorDim alpha (-1.0, 1.0,-1.0);
			MatrixDimNslips alpha_BC_CD_DB;
			/*                BC   CD   DB */
			alpha_BC_CD_DB<< 1.0,-1.0, 0.0,	
			/*            */ 0.0,-1.0, 1.0,	
			/*            */-1.0, 0.0, 1.0;
			
			
			VectorDim beta  (1.0,-1.0,-1.0);
			//			/**/VectorDim CA( 0.0,-1.0, 1.0);	// beta CA
			//			/**/VectorDim AD(-1.0, 0.0,-1.0);	// beta AD
			//			/**/VectorDim DC( 1.0, 1.0, 0.0);	// beta DC
			MatrixDimNslips beta_CA_AD_DC;
			beta_CA_AD_DC<< 0.0,-1.0, 1.0,	
			/*           */ -1.0, 0.0, 1.0,	
			/*           */ 1.0,-1.0, 0.0;
			
			
			VectorDim gamma (-1.0,-1.0, 1.0);
			//			/**/VectorDim AB(-1.0, 1.0, 0.0);	// gamma AB
			//			/**/VectorDim BD(0.0,-1.0,-1.0);	// gamma BD
			//			/**/VectorDim DA(1.0, 0.0, 1.0);	// gamma DA
			MatrixDimNslips gamma_AB_BD_DA;
			gamma_AB_BD_DA<< -1.0, 0.0, 1.0,	
			/*            */ 1.0,-1.0, 0.0,	
			/*            */ 0.0,-1.0, 1.0;
			
			
			VectorDim  delta(1.0, 1.0, 1.0);
			//			/**/VectorDim AC( 0.0, 1.0,-1.0);	// delta AC
			//			/**/VectorDim CB(-1.0, 0.0, 1.0);	// delta CB
			//			/**/VectorDim BA( 1.0,-1.0, 0.0);	// delta BA
			MatrixDimNslips delta_AC_CB_BA;
			/*                 AC  CB   BA */
			delta_AC_CB_BA<< 0.0,-1.0, 1.0,	
			/*            */ 1.0, 0.0,-1.0,	
			/*            */-1.0, 1.0, 0.0;	
			
			
//						Eigen::Matrix<double,dim,dim> RT;
//						
//			RT.col(2) = delta.normalized();
//			RT.col(0) = delta_AC_CB_BA.col(0).normalized();
//			RT.col(1) = RT.col(2).cross(RT.col(0));
//						
//			Eigen::Matrix<double,dim,dim> C2G=RT.transpose();
			//			Eigen::Matrix<double,dim,dim> R=Eigen::Matrix<double,dim,dim>::Identity();
			
			this->insert( SlipSystem<dim,Nslips>(C2G*alpha, C2G*alpha_BC_CD_DB));
			this->insert( SlipSystem<dim,Nslips>(C2G*beta,  C2G*beta_CA_AD_DC));
			this->insert( SlipSystem<dim,Nslips>(C2G*gamma, C2G*gamma_AB_BD_DA));
			this->insert( SlipSystem<dim,Nslips>(C2G*delta, C2G*delta_AC_CB_BA));
			
		}
		
		
		///////////////////////////////////////////////////////////////////////
		Crystal(){
			rotate(Eigen::Matrix<double,dim,dim>::Identity());
		}
		

		
		VectorDim conjugatePlaneNormal(const VectorDim& B, const VectorDim& N, const double& tol=FLT_EPSILON) const {
		
			int count(0);
			VectorDim temp;
			for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=this->begin();iter!=this->end();++iter){
				
				for (int k=0;k<Nslips;++k){
					if(	 B.normalized().cross(iter->slip.col(k)).norm()<tol 
					   && N.normalized().cross(iter->normal).norm()>tol){
						temp=iter->normal;
						++count;
					}
				}
			}
			assert(count==1 && "FOUND WRONG NUMBER OF CONJUGATE PLANES");
			return temp;
		}
		
		
//		VectorDim chooseCrossSlipNormal(const VectorDim& chord, const VectorDim& bOther, const VectorDim& nOther){
//			
//			if (normalizedChord.dot(nOther)<)
//				
//				shared.material.conjugatePlaneNormal(Burgers,this->glidePlaneNormal)
//				
//				}
		
		
	};
	
	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// FCCP CRYSTAL
	template<short unsigned int dim>
	class Crystal<dim,FCCP> : public CrystalBase<dim,2> {
		
		
		
	public:
		enum{Nslips=2};
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		typedef Eigen::Matrix<double,dim,Nslips> MatrixDimNslips;
		//#include <model/Dislocations/Materials/Crystal_common.h>
		
		
		
		
		
		Crystal(){
			//this->clear();
			
			
			// 1- Define the 2 slip normals and directions
			VectorDim  delta(0.0, 0.0, 1.0);
			MatrixDimNslips delta_Cd_dB;
			/*                               Cd                     dB */
			delta_Cd_dB<< std::pow(3.0,0.5)/2.0, std::pow(3.0,0.5)/2.0,
			/*                      */		0.5,				  -0.5,
			/*                      */		0.0,				   0.0;
			
			
			VectorDim alpha  (0,-2.0*std::pow(2.0,0.5)/3.0, -1.0/3.0);
			MatrixDimNslips alpha_Ca_aB;
			/*                               Ca                      aB */
			alpha_Ca_aB<< std::pow(3.0,0.5)/2.0,  std::pow(3.0,0.5)/2.0,
			/*                      */ -1.0/6.0,	            1.0/6.0,
			/*        */  std::pow(2.0,0.5)/3.0, -std::pow(2.0,0.5)/3.0;
			
			
			this->insert( SlipSystem<dim,Nslips>(delta,delta_Cd_dB));
			this->insert( SlipSystem<dim,Nslips>(alpha,alpha_Ca_aB));
			
		}
		
	};
	
	
	
	
	
	
	//////////////////////////////////////////////////////////////
} // namespace model 
#endif

