/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */



#ifndef model_LINELINEINTERSECTION_H_
#define model_LINELINEINTERSECTION_H_

#include <Eigen/Dense>
//#include <Eigen/Geometry>
#include <Eigen/LU> 
//#include <float.h>
#include <assert.h>
//#include <algorithm>
//#include <iterator>
//#include <map>
#include <set>
#include <utility> // defines std::pair
#include <math.h>  // defines std::fabs
#include <float.h> // defines FLT_EPSILON
//#include <boost/math/special_functions/binomial.hpp>
//#include <model/Utilities/GeneralizedEigenSolver.h>
#include <model/Geometry/ColinearLinesIntersection.h>
//#include "model/Patterns/Network/NetworkLink.h"
//#include "model/Geometry/ParametricCurve.h"

namespace model {
	
	/**************************************************************************************/
	/* LineLineIntersection<2>: general case **********************************************/
	/**************************************************************************************/
	/*! /brief Class template specialization for determining the intersection points between 
	 *  two staight line segments in two dimensions.
	 *
	 *  The first line segment is assumed to be prescribed in the form:
	 *	\f[
	 *		\mathbf{x}=\mathbf{P0}+\left(\mathbf{P1}-\mathbf{P0}\right)u
	 *	\f]
	 *  with \f$u\in [0,1]\f$. Analogously the second line segments is assumed to be in the form:
	 *	\f[
	 *		\mathbf{x}=\mathbf{P2}+\left(\mathbf{P3}-\mathbf{P2}\right)v
	 *	\f]
	 * with \f$v\in [0,1]\f$. Intersection points are the values
	 *	\f[
	 *		\left[\mathbf{P1}-\mathbf{P0} \mathbf{P2}-\mathbf{P3}\right]v
	 *	\f]
	 */
	template <short unsigned int dim>
	class LineLineIntersection {
		
	public:
		LineLineIntersection(){
			assert(0 && "LineLineIntersection: template specialization not implemented");
		}
		
	};
	
	
	/**************************************************************************************/
	/* LineLineIntersection<2>: template specializatio dim=2 ******************************/
	/**************************************************************************************/
	template <>
	class LineLineIntersection<2> {
		//		enum {dim=2, corder=0};
		enum {dim=2};
		
		
		
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		
		
		//	#include<model/Geometry/Splines/SplineEnums.h>
		//	typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		//	typedef Eigen::Matrix<double,dim,1>   VectorDim;
		
		
	public:
		
		//! The container of intersection points
		
		//		std::map<double,VectorDim> sortedIntersections1;
		//		std::map<double,VectorDim> sortedIntersections2;
		
		//	std::vector<std::pair<double,double> > intersectionParameters;
		std::set<std::pair<double,double> > intersectionParameters;
		//////////////////////////////////////////////////////////////
		// constructor
		LineLineIntersection(const VectorDim& P0,
							 /*                */ const VectorDim& P1,
							 /*                */ const VectorDim& P2,
							 /*                */ const VectorDim& P3,
							 /*                */ const double& tol=FLT_EPSILON) {
			//			intersectionParameters.clear();
			
			
			//			VectorDim P0=P0P1.col(0);
			//			VectorDim P1=P0P1.col(1);
			//			VectorDim P2=P2P3.col(0);
			//			VectorDim P3=P2P3.col(1);
			
			
			MatrixDim L;
			L<<P1-P0, P2-P3;
			VectorDim R;
			R<< P2-P0;
			
			if (std::fabs(L.determinant())>=tol){	// segments are non-parallel
				VectorDim sol=L.inverse()*R;
				//if (sol(0)>=0.0+dU && sol(0)<=1.0-dU && sol(1)>=0.0+dU && sol(1)<=1.0-dU){
				if (sol(0)>=0.0 && sol(0)<=1.0 && sol(1)>=0.0 && sol(1)<=1.0){ // take solutions only in [0,1] for both parameters
					//intersectionParameters.push_back(std::make_pair(sol(0),sol(1)));
					intersectionParameters.insert(std::make_pair(sol(0),sol(1)));
					//					sortedIntersections1.insert(std::make_pair(sol(0),P0+(P1-P0)*sol(0)));
					//					sortedIntersections2.insert(std::make_pair(sol(1),P0+(P1-P0)*sol(0)));
				}
				
			}
			else{ // segments are parallel
				if(R.squaredNorm()<tol){ // segments are colinear
					//assert(0 && "FINISH HERE, INTERSECTIONS SHOULD RETURN THE INNER POINTS");
					ColinearLinesIntersection<2> cli(P0,P1,P2,P3);
					intersectionParameters=cli.intersectionParameters;
					//					sortedIntersections1=cli.sortedIntersections1;
					//					sortedIntersections2=cli.sortedIntersections2;
					
					
					
				}
			}
		}			
		
	};
	
	
	
	
	/**************************************************************************************/
	/* LineLineIntersection<3>: template specializatio dim=3 ******************************/
	/**************************************************************************************/
	template <>
	class LineLineIntersection<3> {
		//		enum {dim=2, corder=0};
		enum {dim=3};
		
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		
		
		
		MatrixDim rotM(Eigen::Matrix<double,dim,1> v1, Eigen::Matrix<double,dim,1> n){
			v1.normalize();
			n.normalize();
			return (MatrixDim()<< v1, n.cross(v1), n).finished().transpose();
		}
		
	public:
		
		
		//		std::map<double,VectorDim> sortedIntersections1;
		//		std::map<double,VectorDim> sortedIntersections2;
		
		//	std::vector<std::pair<double,double> > intersectionParameters;
		std::set<std::pair<double,double> > intersectionParameters;
		
		LineLineIntersection(const VectorDim& P0, const VectorDim& P1, 
							 const VectorDim& P2, const VectorDim& P3, 
							 const double& tol=FLT_EPSILON){
			
			VectorDim L01=P1-P0;
			VectorDim L23=P3-P2;			
			assert(L01.norm()>=tol && "POINTS P0 and P1 are the same for prescribed tolerance");
			assert(L23.norm()>=tol && "POINTS P3 and P2 are the same for prescribed tolerance");
			VectorDim N=L01.cross(L23);
			
			if (std::abs(N.dot(P2-P0))<=tol){ // segments are co-planar because either |N|=0 (parrallel) or (P2-P0) has no component on N 
				if (N.norm()<tol ){ // co-planar and parallel 
					if (L01.cross(P2-P0).norm()<tol){ // co-planar, parallel and colinear
						ColinearLinesIntersection<3> cli(P0,P1,P2,P3);
						intersectionParameters=cli.intersectionParameters;
						//						sortedIntersections1=cli.sortedIntersections1;
						//						sortedIntersections2=cli.sortedIntersections2;
					}
				}
				else{  // co-planar and non parallel 
					MatrixDim R = rotM(L01,N);
					LineLineIntersection<2> lli((R*(P0-P0)).segment<dim-1>(0),
					/*                       */ (R*(P1-P0)).segment<dim-1>(0),
					/*                       */ (R*(P2-P0)).segment<dim-1>(0),
					/*                       */ (R*(P3-P0)).segment<dim-1>(0));
					intersectionParameters=lli.intersectionParameters;
					//					sortedIntersections1=lli.sortedIntersections1;
					//					sortedIntersections2=lli.sortedIntersections2;
				}
			}
			else { // segments are skew (non co-planar)
//				std::cout<<"SplineIntersection<3,porder,2,2>: no intersection"<<std::endl;
			}
		}
		
		
		
	};
	
	
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif




//					double u2=(P1-P0).dot(P2-P0)/(P1-P0).squaredNorm();
//					if (u2>=0.0 && u2<=1.0){
//						intersectionParameters.insert(std::make_pair(u2,0.0));
//					}
//					
//					double u3=(P1-P0).dot(P3-P0)/(P1-P0).squaredNorm();
//					if (u3>=0.0 && u3<=1.0){
//						intersectionParameters.insert(std::make_pair(u3,1.0));
//					}
//					
//					double v0=(P3-P2).dot(P0-P2)/(P3-P2).squaredNorm();
//					if (v0>=0.0 && v0<=1.0){
//						intersectionParameters.insert(std::make_pair(0.0,v0));
//					}
//					
//					double v1=(P3-P2).dot(P1-P2)/(P3-P2).squaredNorm();
//					if (v1>=0.0 && v1<=1.0){
//						intersectionParameters.insert(std::make_pair(1.0,v1));
//					}
