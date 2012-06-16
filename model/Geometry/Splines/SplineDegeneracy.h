/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 Benjamin Ramirez<ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SPLINEDEGENERACY_H_
#define model_SPLINEDEGENERACY_H_

#include <Eigen/Dense>
#include <model/Geometry/Splines/Coeff2Hermite.h>

//#include <Eigen/Geometry>


namespace model {
	
	/***************************************************************************/
	/* SplineDegeneracy<dim,polyDegree> general case ***************************/
	/***************************************************************************/
	template <short unsigned int dim, short unsigned int polyDegree>
	struct SplineDegeneracy {
		
	};
	
	
	/***************************************************************************/
	/* SplineDegeneracy<2,3> template specialization ***************************/
	/***************************************************************************/
	template <>
	struct SplineDegeneracy<2,3> {
		
		enum {dim=2};
		enum {polyDegree=3};
		enum {polyCoeff=polyDegree+1};
		
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		
		static bool isLine(const Eigen::Matrix<double,dim,polyCoeff> & Coef, const double& tol){
			
			//! 
			//			assert(polyDegree==3 && "THIS IS TEMPORARY!!! FIND CRITERIA FOR DEGENERACY FOR ARBITRARY POLYDEGREE");
			
			Eigen::Matrix<double,dim,polyCoeff> Hermite = Coeff2Hermite<polyDegree>::c2h<dim>(Coef);
			
			
			VectorDim T1    = Hermite.col(1).normalized();
			VectorDim T2    = Hermite.col(3).normalized();
			VectorDim chord = Hermite.col(2)-Hermite.col(0);
			
			//			assert(T1.squaredNorm()>tol && "T1 is too small.");
			//			assert(T2.squaredNorm()>tol && "T2 is too small.");
			assert(chord.squaredNorm()>tol && "chord is too small.");
			//T1.cross(T2);
			
			MatrixDim T1T2, T1chord;
			T1T2<<T1,T2;
			T1chord<< T1,chord;
			
			return std::fabs(T1T2.determinant())<tol &&  std::fabs(T1chord.determinant())<tol;
		}
		
	};
	
	
	/***************************************************************************/
	/* SplineDegeneracy<2,3> template specialization ***************************/
	/***************************************************************************/
	template <>
	struct SplineDegeneracy<3,3> {
		
		enum {dim=3};
		enum {polyDegree=3};
		enum {polyCoeff=polyDegree+1};
		
		typedef Eigen::Matrix<double,dim,1>   VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		
		
		/******************/
		static bool isPlanar(const Eigen::Matrix<double,dim,polyCoeff> & Coef, VectorDim& N, const double& tol){
			
			/*! IS PLANAR SHOULD RETURN 0 IF IS A LINE!!!!! CHANGE ALSO IN SPLINESEGMENTBASE
			 *
			 */
			
			bool temp(true);
			N.setZero();

			if (!isLine(Coef, tol)){
				Eigen::Matrix<double,dim,polyCoeff> Hermite (Coeff2Hermite<polyDegree>::c2h<dim>(Coef));
				const VectorDim T0(Hermite.col(1));
				const VectorDim T1(Hermite.col(3));
				const VectorDim chord (Hermite.col(2)-Hermite.col(0));
				assert(chord.norm()>tol && "chord is too small.");

				VectorDim N0(T0.cross(chord.normalized())); // N0 can be 0 either because T0=0 or because T0 and chord are parallel
				VectorDim N1(T1.cross(chord.normalized())); // N1 can be 0 either because T1=0 or because T1 and chord are parallel
				
				if (N0.norm()>=tol && N1.norm()>=tol){
//					N0.normalize();
//					N1.normalize();
					//temp=(N0.normalized().cross(N1.normalized()).norm()<tol);
					temp=(N0.cross(N1).norm()<tol*N0.norm()*N1.norm());

					if(temp){
						N=N0.normalized();
					}

				}
				else if (N0.norm()>=tol && N1.norm()<tol) { // if one of the two tangets is 0 then the spline is automatically planar
					temp=true;
					N=N0.normalized();
				}
				else if (N1.norm()>=tol || N0.norm()<tol) { // if one of the two tangets is 0 then the spline is automatically planar
					temp=true;
					N=N1.normalized();
				}
				else{
					assert(0 && "THE SPLINE SHOULD BE A LINE.");
				}
			}

			return temp;
		}
		
		
		
		/******************/
		static bool isLine(const Eigen::Matrix<double,dim,polyCoeff> & Coef, const double& tol){
			
			Eigen::Matrix<double,dim,polyCoeff> Hermite = Coeff2Hermite<polyDegree>::c2h<dim>(Coef);
			assert(polyDegree==3 && "THIS IS TEMPORARY!!! FIND CRITERIA FOR DEGENERACY FOR ARBITRARY POLYDEGREE");
			VectorDim T0    = Hermite.col(1);
			VectorDim T1    = Hermite.col(3);
			VectorDim chord = Hermite.col(2)-Hermite.col(0);
			
			
			assert(chord.norm()>tol && "chord is too small.");
			//	 assert(T0.squaredNorm()>tol && "T1 is too small.");
			//	 assert(T1.squaredNorm()>tol && "T2 is too small.");
			
			bool temp(false);
			if (T0.norm()>=tol && T1.norm()>=tol){
				temp=(T0.normalized().cross(chord.normalized()).norm()<tol && T1.normalized().cross(chord.normalized()).norm()<tol);
			}
			else if (T0.norm()>=tol && T1.norm()<tol){
				temp=T0.normalized().cross(chord.normalized()).norm()<tol;
			}
			else if (T1.norm()>=tol && T0.norm()<tol){
				temp=T1.normalized().cross(chord.normalized()).norm()<tol;
			}
			else{
				temp=true;
			}
			
			return temp;
		}
	
};


//////////////////////////////////////////////////////////////s
} // namespace model
#endif
