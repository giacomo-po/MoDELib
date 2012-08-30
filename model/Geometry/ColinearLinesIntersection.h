/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_COLINEARLINESINTERSECTION_H_
#define model_COLINEARLINESINTERSECTION_H_

#include <assert.h>
//#include <map>
#include <set>
#include <utility> // defines std::pair
#include <math.h>  // defines std::fabs
#include <float.h> // defines FLT_EPSILON
#include <Eigen/Dense>

namespace model {
	
	/*****************************************************************************/
	/* ColinearLinesIntersection<dim> ********************************************/
	/*****************************************************************************/
	/*! /brief Class template for determining the intersection points between two 
	 *  clinear line segments in two dimensions.
	 *
	 *  The first line segment is assumed to be prescribed in the form:
	 *	\f[
	 *		\mathbf{x}=\mathbf{P0}+\left(\mathbf{P1}-\mathbf{P0}\right)u
	 *	\f]
	 *  with \f$u\in [0,1]\f$. Analogously the second line segments is assumed to 
	 *  be in the form:
	 *	\f[
	 *		\mathbf{x}=\mathbf{P2}+\left(\mathbf{P3}-\mathbf{P2}\right)v
	 *	\f]
	 * with \f$v\in [0,1]\f$. Intersection points are the values
	 *	\f[
	 *		\left[\mathbf{P1}-\mathbf{P0} \mathbf{P2}-\mathbf{P3}\right]v
	 *	\f]
	 */
	template <short unsigned int dim>
	class ColinearLinesIntersection {
		
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		

		
	public:
		
		//! The container of intersection points
		
//		std::map<double,VectorDim> sortedIntersections1;
//		std::map<double,VectorDim> sortedIntersections2;
		std::set<std::pair<double,double> > intersectionParameters; // use set to avoid inserting repeated intersection points
		
		/*****************************************************************************/
		ColinearLinesIntersection(const VectorDim& P0,
		/*                     */ const VectorDim& P1,
		/*                     */ const VectorDim& P2,
		/*                     */ const VectorDim& P3,
		/*                     */ const double& tol=FLT_EPSILON){
			
			const double snP0P1=(P1-P0).squaredNorm();
			const double snP2P3=(P3-P2).squaredNorm();
			
			// Check that line segments are not too small (need to divide by squard norm)
			assert(snP0P1>tol && "Vector P1-P0 is too small.");
			assert(snP2P3>tol && "Vector P3-P2 is too small.");
			// Check that the two line segments are actually colinear
			assert(std::fabs(std::fabs((P1-P0).dot(P2-P0)) - (P1-P0).norm()*(P2-P0).norm())<tol && "(P1-P0) and (P2-P0) are not colinear.");
			assert(std::fabs(std::fabs((P1-P0).dot(P3-P0)) - (P1-P0).norm()*(P3-P0).norm())<tol && "(P1-P0) and (P3-P0) are not colinear.");
			
			/*            u2                      */
			/* P0 +>------+--------->+ P1         */
			/*         P2 +>---------------->+ P3 */
			const double u2=(P1-P0).dot(P2-P0)/snP0P1;
			if (u2>=0.0 && u2<=1.0){
				intersectionParameters.insert(std::make_pair(u2,0.0));
//				sortedIntersections1.insert(std::make_pair(u2,P0+(P1-P0)*u2));
//				sortedIntersections2.insert(std::make_pair(0.0,P2));
			}
			
			/*                       u3           */
			/*         P0 +>---------+------>+ P1 */
			/* P2 +>---------------->+ P3         */
			const double u3=(P1-P0).dot(P3-P0)/snP0P1;
			if (u3>=0.0 && u3<=1.0){
				intersectionParameters.insert(std::make_pair(u3,1.0));
//				sortedIntersections1.insert(std::make_pair(u3,P0+(P1-P0)*u3));
//				sortedIntersections2.insert(std::make_pair(1.0,P3));
			}
			
			const double v0=(P3-P2).dot(P0-P2)/snP2P3;
			if (v0>=0.0 && v0<=1.0){
				intersectionParameters.insert(std::make_pair(0.0,v0));
//				sortedIntersections1.insert(std::make_pair(0.0,P0));
//				sortedIntersections2.insert(std::make_pair(v0,P2+(P3-P2)*v0));
			}
			
			const double v1=(P3-P2).dot(P1-P2)/snP2P3;
			if (v1>=0.0 && v1<=1.0){
				intersectionParameters.insert(std::make_pair(1.0,v1));
//				sortedIntersections1.insert(std::make_pair(1.0,P1));
//				sortedIntersections2.insert(std::make_pair(v1,P2+(P3-P2)*v1));
			}
			
			assert(intersectionParameters.size()<=2 && "SOMETHING WENT REALLY WRONG HERE.");
			
		}			
		
	};
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif