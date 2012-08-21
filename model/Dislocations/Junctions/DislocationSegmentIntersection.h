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


namespace model {
	
    
	/******************************************************************************************************************/
	/* DislocationSegmentIntersection: general case  ******************************************************************************/
	/******************************************************************************************************************/
	template <short unsigned int dim, short unsigned int polyDegree>
	class DislocationSegmentIntersection {
		
	public:
        
        static double compTol;

        
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
		int planePlaneType(const VectorDim & normal1, const VectorDim & normal2, const VectorDim & pt1, const VectorDim & pt2) const {

			bool areParallelNormals(normal1.cross(normal2).norm()<FLT_EPSILON);
			bool areCoincidentPoints(std::fabs((pt1-pt2).dot(normal1))<FLT_EPSILON);
			
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
		
		static double compTol;
		
		/* Constructor *****************************************************/
		DislocationSegmentIntersection(const MatrixDimPolyCoeff& H1_in, 
        /*                          */ const VectorDim& n1_in) : H1(H1_in),
		/*                                                    */ n1(n1_in),
		/*                                                    */ P0 (H1.col(0)),   // end point of the spline
		/*                                                    */ T0 (H1.col(1)),   // tangent of the spline source
		/*                                                    */ P1 (H1.col(2)),   // end point of the spline
		/*                                                    */ T1 (H1.col(3)),   // tangent of the spline sink
		/*                                                    */ chord1 (P1-P0){   // end point of the spline	
			assert(std::fabs(n1.dot(chord1))<FLT_EPSILON);
			assert(std::fabs(n1.dot(H1.col(1)))<FLT_EPSILON);
			assert(std::fabs(n1.dot(H1.col(3)))<FLT_EPSILON);
		}
		
		
		/* intersectWith ***************************************************/
		//template <typename T>
		std::set<std::pair<double,double> > intersectWith(const MatrixDimPolyCoeff& H2, const VectorDim& n2, const double& physTol) const {
			
			
			std::set<std::pair<double,double> > intersectionParameters;
			
			
			
			// chech that H1 is planar
			
			const double absN2dotC2(std::fabs(n2.dot(H2.col(2)-H2.col(0))));
			const double absN2dotT2source(std::fabs(n2.dot(H2.col(1))));
			const double   absN2dotT2sink(std::fabs(n2.dot(H2.col(3))));
			if (absN2dotC2>FLT_EPSILON || absN2dotT2source>FLT_EPSILON || absN2dotT2sink>FLT_EPSILON){
				//std::cout<<"Dislocation Segment "<<pT->nodeIDPair.first<<"->"<<pT->nodeIDPair.second<<std::endl;
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
				
				
				
				
				
				
				
				
				
				const int planesType=planePlaneType(n1,n2,P0,P2);
				
				
				switch (planesType) {
					case 1:{ // coplanar planes
						
       //                 											std::cout<<"Coplanar case ";
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
                        PlanarSplineImplicitization<polyDegree>::physTol=physTol;
//						intersectionParameters=sli.intersectWith<polyDegree>(Coeff2Hermite<polyDegree>::h2c<dim-1>(H2L),physTol);
						intersectionParameters=sli.intersectWith<polyDegree>(Coeff2Hermite<polyDegree>::h2c<dim-1>(H2L));
						
						break;
					}
						
					case 2:{ // incident planes
						
      //                  						std::cout<<"Incident case ";
                        
						
						std::set<std::pair<double,double> > lineIntersectionParameters1;
						std::set<std::pair<double,double> > lineIntersectionParameters2;
						
						
						//					std::cout<<"Incident case"<<std::endl;
						//					std::cout<<"n1 is"<<n1.transpose()<<std::endl;
						//					std::cout<<"n2 is"<<n2.transpose()<<std::endl;					
						// finding the common line
						const double denom(1.0-std::pow(n1.dot(n2),2));
						const double numer((P2-P0).dot(n2));
						
						//					std::cout<<"denom="<<denom<<std::endl;
						
						if(std::fabs(denom)>compTol){ // planes are incident 
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
//							lineIntersectionParameters1=sli1.intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L),physTol);
                            PlanarSplineImplicitization<polyDegree>::physTol=physTol;
							lineIntersectionParameters1=sli1.intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L));
							
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
//							lineIntersectionParameters2=sli2.intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L),physTol);
                            PlanarSplineImplicitization<polyDegree>::physTol=physTol;
							lineIntersectionParameters2=sli2.intersectWith<1>(Coeff2Hermite<1>::h2c<dim-1>(H2L));
							
                            //												std::cout<<"intersections of spline 2 and line are:"<<std::endl;
                            //													for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters2.begin();iter1!=lineIntersectionParameters2.end();++iter1){
                            //														std::cout<<iter1->first<<" "<<iter1->second<<std::endl;
                            //													}
							
							
							// du = dl / j = dl/L for a line
							
							double tolU=physTol/(P5-P4).norm(); 						
//							
//                            std::cout<<"P4="<<P4.transpose()<<std::endl;
//                            std::cout<<"P5="<<P5.transpose()<<std::endl;
//                            for (std::set<std::pair<double,double> >::const_iterator iter1=lineIntersectionParameters1.begin();iter1!=lineIntersectionParameters1.end();++iter1){
//                                std::cout<<"spline1 and line: "<<iter1->first<<" "<<iter1->second<<std::endl;
//                            }
//                            for (std::set<std::pair<double,double> >::const_iterator iter2=lineIntersectionParameters2.begin();iter2!=lineIntersectionParameters2.end();++iter2){
//                                std::cout<<"spline2 and line: "<<iter2->first<<" "<<iter2->second<<std::endl;
//                            }
                            
							
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
	
    
    
    // Declare statica data members
    //template <short unsigned int dim, short unsigned int polyDegree>
	double DislocationSegmentIntersection<3,3>::compTol=FLT_EPSILON;
	
	
	
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif

