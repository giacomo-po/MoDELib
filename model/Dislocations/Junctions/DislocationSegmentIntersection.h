/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Tamer Crsoby <tamercrosby@gmail.com>, 
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DISLOCATIONSEGMENTINTERSECTION_H_
#define model_DISLOCATIONSEGMENTINTERSECTION_H_

#include <float.h>
#include <assert.h>
#include <set>
#include <map>
#include <utility>
#include <math.h>
#include <float.h>
#include <Eigen/Dense>
//#include <model/Geometry/LineLineIntersection.h>
#include <model/Geometry/Splines/Intersection/PlanarSplineImplicitization.h>
//#include <model/Geometry/Splines/SplineDegeneracy.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/Dislocations/DislocationLocalReference.h>

//#include <Eigen/Geometry>
//#include <Eigen/LU> 
//#include <algorithm>
//#include <iterator>
//#include <boost/math/special_functions/binomial.hpp>
//#include <model/Utilities/GeneralizedEigenSolver.h>
//#include <model/Geometry/Splines/Intersection/SplineImplicitization.h>
//#include <model/Geometry/Splines/Intersection/SplineLineIntersection.h>
//#include <model/Geometry/Splines/Intersection/SplineDislocationSegmentIntersection.h>
//#include "model/Patterns/Network/NetworkLink.h"
//#include "model/Geometry/ParametricCurve.h"
namespace model {
	
	
	//	Eigen::Matrix<double,3,3> rotMatrix(Eigen::Matrix<double,3,1> v1, Eigen::Matrix<double,3,1> n){
	//		v1.normalize();
	//		n.normalize();
	//		return (Eigen::Matrix<double,3,3>()<< v1, n.cross(v1), n).finished().transpose();
	//	}
	
	
	
	/******************************************************************************************************************/
	/* DislocationSegmentIntersection: general case  ******************************************************************************/
	/******************************************************************************************************************/
	template <short unsigned int dim, short unsigned int polyDegree>
	class DislocationSegmentIntersection {
		
	public:
		DislocationSegmentIntersection(){
			assert(0 && "DislocationSegmentIntersection: TEMPLATE SPECIALIZATION NOT IMPLEMENTED.");
			
		}
	};
	
	
	
	/******************************************************************************************************************/
	/* DislocationSegmentIntersection<3,3,1,1>: template specializatio (planar spline-spline intersection) ************************/
	/******************************************************************************************************************/
	template <>
	class DislocationSegmentIntersection<3,3> {
		
		enum {dim=3,polyDegree=3};
		enum {polyCoeff=polyDegree+1};
		typedef Eigen::Matrix<double,dim,polyCoeff> MatrixDimPolyCoeff;
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		
		
		
		/************************************************************************************/
		int planePlaneType(const VectorDim & normal1, const VectorDim & normal2, const VectorDim & pt1, const VectorDim & pt2,const double& tol) const {
			
			
			//			std::cout<<"normal1="<<normal1.transpose()<<std::endl;
			//			std::cout<<"normal2="<<normal2.transpose()<<std::endl;
			//			std::cout<<"normal cross normal"<<(normal1.cross(normal2).norm())<<std::endl;
			//			std::cout<<"tol"<<tol<<std::endl;
			//			std::cout<<"normal cross normal tol"<<(normal1.cross(normal2).norm()<tol)<<std::endl;
			//
			//			std::cout<<"points dot normal"<<(std::fabs((pt1-pt2).dot(normal1))>tol)<<std::endl;
			//			std::cout<<"points dist"<<(std::fabs((pt1-pt2).dot(normal1)))<<std::endl;
			bool areParallelNormals(normal1.cross(normal2).norm()<100.0*tol);
			bool areCoincidentPoints(std::fabs((pt1-pt2).dot(normal1))<tol);
			
			int i;
			if(areParallelNormals && !areCoincidentPoints){
				/*parallel planes: no intersection*/
				i= 0;	
			}else if(areParallelNormals && areCoincidentPoints){
				/*unique planes: planar intersection intersection*/
				i= 1;
			}else {
				/*angle planes: angular intersection*/
				i= 2;
			}
			
			return i;
		}
		
		
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		
		const MatrixDimPolyCoeff H1;
		const VectorDim n1;
		const VectorDim P0;	  // end point of the spline
		const VectorDim T0;	  // tangent of the spline source
		const VectorDim P1;   // end point of the spline
		const VectorDim T1;   // tangent of the spline sink
		const VectorDim chord1;	  // end point of the spline	
		//		const Eigen::Matrix<double,dim,dim> R1;
		
		
		
		/* Constructor *****************************************************/
		
		DislocationSegmentIntersection(const MatrixDimPolyCoeff& H1_in, 
		/*                          */ const VectorDim& n1_in) : H1(H1_in),
		/*                                                    */ n1(n1_in),
		/*                                                    */ P0 (H1.col(0)),   // end point of the spline
		/*                                                    */ T0 (H1.col(1)),   // tangent of the spline source
		/*                                                    */ P1 (H1.col(2)),   // end point of the spline
		/*                                                    */ T1 (H1.col(3)),   // tangent of the spline sink
		/*                                                    */ chord1 (P1-P0){   // end point of the spline	
			//		/*                                                    */ R1(DislocationLocalReference<dim>::global2local(chord1,n1)){
			// chech that H1 is planar
			//			std::cout<<"n1="<<n1.transpose()<<std::endl;
			//			std::cout<<"c1="<<chord1.transpose()<<std::endl;
			//			std::cout<<"|n1 x c1| ="<<n1.cross(chord1).norm()<<std::endl;
			assert(std::fabs(n1.dot(chord1))<FLT_EPSILON);
			assert(std::fabs(n1.dot(H1.col(1)))<FLT_EPSILON);
			assert(std::fabs(n1.dot(H1.col(3)))<FLT_EPSILON);
			
			
			
		}
		
		
		/* intersectWith ***************************************************/
		template <typename T>
		std::set<std::pair<double,double> > intersectWith(const MatrixDimPolyCoeff& H2, const VectorDim& n2, const double& tol,
														  const T* const pT){
			
			
			std::set<std::pair<double,double> > intersectionParameters;
			
			
			
			// chech that H1 is planar
			
			const double absN2dotC2(std::fabs(n2.dot(H2.col(2)-H2.col(0))));
			const double absN2dotT2source(std::fabs(n2.dot(H2.col(1))));
			const double   absN2dotT2sink(std::fabs(n2.dot(H2.col(3))));
			if (absN2dotC2>FLT_EPSILON || absN2dotT2source>FLT_EPSILON || absN2dotT2sink>FLT_EPSILON){
				std::cout<<"Dislocation Segment "<<pT->nodeIDPair.first<<"->"<<pT->nodeIDPair.second<<std::endl;
				std::cout<<"n2="<<n2.transpose()<<std::endl;
							std::cout<<"H2="<<H2<<std::endl;
				std::cout<<"absN2dotC2="<<absN2dotC2<<"\n";
				std::cout<<"absN2dotT2source="<<absN2dotT2source<<"\n";
				std::cout<<"absN2dotT2sink="<<absN2dotT2sink<<"\n";
//							std::cout<<"n2"std::fabs(n2.dot(H2.col(3)))<<std::endl;				
			}
			
			assert(      absN2dotC2<FLT_EPSILON);
			assert(absN2dotT2source<FLT_EPSILON);
			assert(  absN2dotT2sink<FLT_EPSILON);
			
			
			const VectorDim P2 = H2.col(0);	  // end point of the spline
			const VectorDim T2 = H2.col(1);	  // tangent of the spline source
			const VectorDim P3 = H2.col(2);   // end point of the spline
			const VectorDim T3 = H2.col(3);   // tangent of the spline sink
			const VectorDim chord2 = P3-P2;	  // end point of the spline	
			
			
			
			
			
			if((0.5*(P0+P1)-0.5*(P2+P3)).norm()<0.75*(chord1.norm()+chord2.norm())){
				
				double physTol=15.0;
				
				
				
				
				
				
				
				const int planesType=planePlaneType(n1,n2,P0,P2,tol);
				
				
				switch (planesType) {
					case 1:{ // coplanar planes
						
//											std::cout<<"Coplanar case"<<std::endl;
						const Eigen::Matrix<double,dim,dim> R1(DislocationLocalReference<dim>::global2local(chord1,n1));
						
						Eigen::Matrix<double,dim-1,polyCoeff> H1L;
						H1L.col(0)=(R1*(P0-P0)).segment<dim-1>(0);
						H1L.col(1)=(R1*T0).segment<dim-1>(0);
						H1L.col(2)=(R1*(P1-P0)).segment<dim-1>(0);
						H1L.col(3)=(R1*T1).segment<dim-1>(0);
						
						Eigen::Matrix<double,dim-1,polyCoeff> H2L;
						H2L.col(0)=(R1*(P2-P0)).segment<dim-1>(0);
						H2L.col(1)=(R1*T2).segment<dim-1>(0);
						H2L.col(2)=(R1*(P3-P0)).segment<dim-1>(0);
						H2L.col(3)=(R1*T3).segment<dim-1>(0);
						
						PlanarSplineImplicitization<polyDegree> sli(Coeff2Hermite<polyDegree>::h2c<dim-1>(H1L)); 
						intersectionParameters=sli.intersectWith<polyDegree>(Coeff2Hermite<polyDegree>::h2c<dim-1>(H2L),tol,physTol);
						
						break;
					}
						
					case 2:{ // incident planes
						
//						std::cout<<"Incident case"<<std::endl;

						
						std::set<std::pair<double,double> > lineIntersectionParameters1;
						std::set<std::pair<double,double> > lineIntersectionParameters2;
						
						
						//					std::cout<<"Incident case"<<std::endl;
						//					std::cout<<"n1 is"<<n1.transpose()<<std::endl;
						//					std::cout<<"n2 is"<<n2.transpose()<<std::endl;					
						// finding the common line
						const double denom(1.0-std::pow(n1.dot(n2),2));
						const double numer((P2-P0).dot(n2));
						
						//					std::cout<<"denom="<<denom<<std::endl;
						
						if(std::fabs(denom)>tol){ // planes are incident 
							double u=numer/denom;
							
							VectorDim linePoint = P0+(n2-n2.dot(n1)*n1)*u;
							VectorDim lineDir   = n1.cross(n2).normalized();
							
							Eigen::Matrix<double,4,1> projPoints;
							projPoints<<(P0-linePoint).dot(lineDir),
							/**********/(P1-linePoint).dot(lineDir),
							/**********/(P2-linePoint).dot(lineDir),
							/**********/(P3-linePoint).dot(lineDir);
							
							Eigen::Matrix<double,dim,1> Pmean=linePoint+0.5*(projPoints.minCoeff()+projPoints.maxCoeff())*lineDir;
							Eigen::Matrix<double,dim,1> P4=Pmean+(chord1.norm()+chord2.norm())*lineDir;
							Eigen::Matrix<double,dim,1> P5=Pmean-(chord1.norm()+chord2.norm())*lineDir;
							
							
							//						std::cout<<"P4 = "<<P4.transpose()<<std::endl;
							//						std::cout<<"P5 = "<<P5.transpose()<<std::endl;
							
							const Eigen::Matrix<double,dim,dim> R1(DislocationLocalReference<dim>::global2local(P5-P4,n1));
							const Eigen::Matrix<double,dim,dim> R2(DislocationLocalReference<dim>::global2local(P5-P4,n2));
							
							
							
							Eigen::Matrix<double,dim-1,polyCoeff> H1L;
							H1L.col(0)=(R1*(P0-P4)).segment<dim-1>(0);
							H1L.col(1)=(R1*T0).segment<dim-1>(0);
							H1L.col(2)=(R1*(P1-P4)).segment<dim-1>(0);
							H1L.col(3)=(R1*T1).segment<dim-1>(0);
							
							Eigen::Matrix<double,dim-1,2> H2L;
							H2L.col(0)=(R1*(P4-P4)).segment<dim-1>(0);
							H2L.col(1)=(R1*(P5-P4)).segment<dim-1>(0);
							
							//						std::cout<<"H1L="<<std::endl;
							//						std::cout<<H1L<<std::endl;
							//						std::cout<<"H2L="<<std::endl;
							//						std::cout<<H2L<<std::endl;
							
							PlanarSplineImplicitization<polyDegree> sli1(Coeff2Hermite<polyDegree>::h2c<dim-1>(H1L)); 
							lineIntersectionParameters1=sli1.intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L),tol,physTol);
							
//													std::cout<<"intersections of spline 1 and line are:"<<std::endl;
//													for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters1.begin();iter1!=lineIntersectionParameters1.end();++iter1){
//														std::cout<<iter1->first<<" "<<iter1->second<<std::endl;
//													}
							
							
							
							Eigen::Matrix<double,dim-1,polyCoeff> H3L;
							H3L.col(0)=(R2*(P2-P4)).segment<dim-1>(0);
							H3L.col(1)=(R2*T2).segment<dim-1>(0);
							H3L.col(2)=(R2*(P3-P4)).segment<dim-1>(0);
							H3L.col(3)=(R2*T3).segment<dim-1>(0);
							
							H2L.col(0)=(R2*(P4-P4)).segment<dim-1>(0);
							H2L.col(1)=(R2*(P5-P4)).segment<dim-1>(0);
							
							PlanarSplineImplicitization<polyDegree> sli2(Coeff2Hermite<polyDegree>::h2c<dim-1>(H3L)); 
							lineIntersectionParameters2=sli2.intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L),tol,physTol);
							
//												std::cout<<"intersections of spline 2 and line are:"<<std::endl;
//													for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters2.begin();iter1!=lineIntersectionParameters2.end();++iter1){
//														std::cout<<iter1->first<<" "<<iter1->second<<std::endl;
//													}
							
							
							// du = dl / j = dl/L for a line
							
							double tolU=physTol/(P5-P4).norm(); //! DislocationSegmentIntersection !!!! CHANGE THIS 5.0 LATER. NO ... FOR REAL						
							
							//					std::cout<<"tolU="<<tolU<<std::endl;
							
							
							for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters1.begin();iter1!=lineIntersectionParameters1.end();){
								
								// find the closest to iter1
								
								std::map<double,std::set<std::pair<double,double> >::const_iterator> compareTo1;
								
								for (std::set<std::pair<double,double> >::const_iterator iter2=lineIntersectionParameters2.begin();iter2!=lineIntersectionParameters2.end();++iter2){
									compareTo1.insert(std::make_pair(std::fabs(iter1->second-iter2->second),iter2));								
								}
								
								if(compareTo1.size()){
									if(compareTo1.begin()->first<tolU){
										double u1=iter1->first;
										double u2=compareTo1.begin()->second->first;
										intersectionParameters.insert(std::make_pair(u1,u2));
										std::set<std::pair<double,double> >::const_iterator toBeErased(iter1);
										++iter1;
										assert(lineIntersectionParameters1.erase(*toBeErased) && "NONE ERASED");
										assert(lineIntersectionParameters2.erase(*(compareTo1.begin()->second)) && "NONE DELETED");
									}
									else{
										++iter1;
									}
								}
								else{
									++iter1;
								}
								
								
							}
							
						}
						else{ // planes are parallel
							assert(0 && "SOMETHING WENT REALLY, REALLY WRONG HERE, YOU SHOULD HAVE FOUND THIS ABOVE");
						}
						
						
						break;
					}
						
					default:
						
						break;
				}
			}
			return intersectionParameters;
		} // close constructor
		
	};			
	
	
	
	
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif


//			std::cout<<"Doing planar spline-spline intersection"<<std::endl;

// Assert that spline1 is planar
//			VectorDim n1;
//			assert(( SplineDegeneracy<dim,polyDegree>::isPlanar(C1,n1,tol)) && "DislocationSegmentIntersection<3,3,1,1>: spline1 is not planar.");
//			assert((!SplineDegeneracy<dim,polyDegree>::isLine  (C1,   tol)) && "DislocationSegmentIntersection<3,3,1,1>: spline1 is a line. You should call DislocationSegmentIntersection<3,polyDegree,2,2>");

// Assert that spline1 is planar
//			VectorDim n2;
//			assert(( SplineDegeneracy<dim,polyDegree>::isPlanar(C2,n2,tol)) && "DislocationSegmentIntersection<3,3,1,1>: spline2 is not planar.");
//			assert((!SplineDegeneracy<dim,polyDegree>::isLine  (C2,   tol)) && "DislocationSegmentIntersection<3,3,1,1>: spline2 is a line. You should call DislocationSegmentIntersection<3,polyDegree,2,2>");





//assert(std::fabs(T2.normalized().cross(T3.normalized()).dot(chord2.normalized()))<tol && "DislocationSegmentIntersection<3,polyDegree,1,1>: spline2 is not planar.");
// Assert that spline1 is not degenerate
//			assert(T2.normalized().cross(chord2.normalized()).norm()>tol && T3.normalized().cross(chord2.normalized()).norm()>tol && 
//				   "DislocationSegmentIntersection<3,3,1,1>: spline2 is degenerate. You should call DislocationSegmentIntersection<3,polyDegree,2,2>");


//			std::cout<<"n1="<<n1.transpose()<<std::endl;
//			std::cout<<"n2="<<n2.transpose()<<std::endl;
//	const VectorDim n1 = get_normal(T0,T1,chord1,tol); 
//	const VectorDim n2 = get_normal(T2,T3,chord2,tol);


//						double physTol=10.0/(P5-P4).norm(); //! DislocationSegmentIntersection !!!! CHANGE THIS 5.0 LATER. NO ... FOR REAL

//						double physTol=0.01;
//						double physTol=20.0;
//						if(lineIntersectionParameters1.size() && lineIntersectionParameters2.size()){

//						for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters1.begin();iter1!=lineIntersectionParameters1.end();++iter1){


//						std::cout<<"lineIntersectionParameters1 are"<<std::endl;
//						for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters1.begin();iter1!=lineIntersectionParameters1.end();++iter1){
//							std::cout<<" "<<iter1->first<<" "<<iter1->second<<std::endl;
//						}
//						std::cout<<"lineIntersectionParameters2 are"<<std::endl;
//						for (std::set<std::pair<double,double> >::const_iterator iter2=lineIntersectionParameters2.begin();iter2!=lineIntersectionParameters2.end();++iter2){
//							std::cout<<" "<<iter2->first<<" "<<iter2->second<<std::endl;
//						}



//					}

//						std::cout<<"intersectionParameters are now"<<std::endl;
//						for (std::set<std::pair<double,double> >::const_iterator iter1=intersectionParameters.begin();iter1!=intersectionParameters.end();++iter1){
//							std::cout<<" "<<iter1->first<<" "<<iter1->second<<std::endl;
//						}





//							for (std::set<std::pair<double,double> >::const_iterator iter2=lineIntersectionParameters2.begin();iter2!=lineIntersectionParameters2.end();++iter2){
////								double u1(iter1->first);
////								VectorDim I1(C1.col(0)+C1.col(1)*u1+C1.col(2)*u1*u1+C1.col(3)*u1*u1*u1);
////								double u2(iter2->first);
////								VectorDim I2(C2.col(0)+C2.col(1)*u2+C2.col(2)*u2*u2+C2.col(3)*u2*u2*u2);
//
////								if((I1-I2).norm()<physTol){
//								if(std::fabs(iter1->second-iter2->second)<physTol){
////									std::cout<<I1.transpose()<<" "<<I2.transpose()<<std::endl;
//									intersectionParameters.insert(std::make_pair(iter1->first,iter2->first));
//								}
//								
//							}



//	/******************************************************************************************************************/
//	/* DislocationSegmentIntersection<3,polyDegree,2,2>: template specializatio (line-line intersection) ******************************/
//	/******************************************************************************************************************/
//	template <short unsigned int polyDegree>
//	class DislocationSegmentIntersection<3,polyDegree,2,2> : public LineLineIntersection<3>{
//		
//		enum{dim=3};
//		typedef Eigen::Matrix<double,dim,1> VectorDim;
//		
//		
//	public:
//		DislocationSegmentIntersection(const VectorDim& P0, const VectorDim& P1, 
//						   /*              */ const VectorDim& P2, const VectorDim& P3, 
//						   /*              */ const double& tol=FLT_EPSILON) : LineLineIntersection<3>(P0,P1,P2,P3,tol){}
//	};
//	
//	/******************************************************************************************************************/
//	/* DislocationSegmentIntersection<3,3,1,2>: template specializatio (planar spline-line intersection) **************************/
//	/******************************************************************************************************************/
//	template <short unsigned int polyDegree>
//	class DislocationSegmentIntersection<3,polyDegree,1,2> {
//		
//		enum {dim=3};
//		//		enum {polyDegree=3};
//		enum {polyCoeff=polyDegree+1};
//		typedef Eigen::Matrix<double,dim,polyCoeff> MatrixDimPolyCoeff;
//		typedef Eigen::Matrix<double,dim,1>   VectorDim;
//		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
//		
//		
//		const MatrixDimPolyCoeff H1;
//		
//	public:
//		
//		
//		std::set<std::pair<double,double> > intersectionParameters;
//		
//		DislocationSegmentIntersection(const MatrixDimPolyCoeff& C1, 
//						   /*              */ const VectorDim& P2, const VectorDim& P3, 
//						   /*              */ const double& tol=FLT_EPSILON) : H1(Coeff2Hermite<polyDegree>::template c2h<dim>(C1)){
//			
//			const VectorDim P0 = H1.col(0);	  // end point of the spline
//			const VectorDim T0 = H1.col(1);	  // tangent of the spline source
//			const VectorDim P1 = H1.col(2);   // end point of the spline
//			const VectorDim T1 = H1.col(3);;  // tangent of the spline sink
//			const VectorDim chord = P1-P0;	  // end point of the spline			
//			
//			assert(chord.norm()>tol && "DislocationSegmentIntersection<3,3,1,2>: spline chord too small.");
//			//			assert(T0.norm()>tol && "DislocationSegmentIntersection<3,3,1,2>: T0 too small.");  // NOT RIGHT FOR DIPOLAR LOOPS
//			//			assert(T1.norm()>tol && "DislocationSegmentIntersection<3,3,1,2>: T1 too small.");	// NOT RIGHT FOR DIPOLAR LOOPS
//			
//			// Assert that spline is planar
//			VectorDim N;
//			assert(( SplineDegeneracy<dim,polyDegree>::isPlanar(C1,N,tol)) && "DislocationSegmentIntersection<3,3,1,2>: is not planar.");
//			assert((!SplineDegeneracy<dim,polyDegree>::isLine  (C1,  tol)) && "DislocationSegmentIntersection<3,3,1,2>: spline is a line. You should call DislocationSegmentIntersection<3,polyDegree,2,2>");
//			
//			
//			double P23N=(P3-P2).dot(N);
//			double P20N=(P0-P2).dot(N);
//			
//			if (std::fabs(P23N)>tol){ // line intersects plane 
//				//				double u=P20N/P23N;
//				//				if (u>=0.0 && u<=1.0){ // line segment intersects plane at P2+(P3-P2)*u 
//				//					VectorDim X=P2+(P3-P2)*u;
//				//					
//				//					assert(0 && "FINISH HERE!!!! NEED TO ROTATE LOCALLY C AND X");
//				//					
//				//				}
//			}
//			else { // line segment is parallel to the plane
//				if (std::fabs(P20N)<tol){ // line is on the plane
//					MatrixDim R=rotMatrix(chord,N);
//					
//					
//					Eigen::Matrix<double,dim-1,polyCoeff> H1L;
//					H1L.col(0)=(R*(P0-P0)).template segment<dim-1>(0);
//					H1L.col(1)=(R*T0).template segment<dim-1>(0);
//					H1L.col(2)=(R*(P1-P0)).template segment<dim-1>(0);
//					H1L.col(3)=(R*T1).template segment<dim-1>(0);
//					
//					
//					Eigen::Matrix<double,dim-1,2> H2L;
//					H2L.col(0)=(R*(P2-P0)).template segment<dim-1>(0);
//					H2L.col(1)=(R*(P3-P0)).template segment<dim-1>(0);
//					
//					
//					PlanarSplineImplicitization<polyDegree> sli(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H1L)); 
//					intersectionParameters=sli.template intersectWith<1>(Coeff2Hermite<1>::template h2c<dim-1>(H2L),tol);
//					
//				} // close line is on the plane
//			} // close line segment is parallel to the plane
//			//			} // close spline is planar and not degenerate
//		} // close constructor
//		
//	};	
