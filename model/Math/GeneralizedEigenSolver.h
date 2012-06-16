/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
 
#ifndef model_GENERALIZEDEIGENSOLVER_H_
#define model_GENERALIZEDEIGENSOLVER_H_


#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <Eigen/Dense>


namespace model {

template <typename T, int n>
class GeneralizedEigenSolver{
	
	
public:
	
	Eigen::Matrix<std::complex<double>,n,1>  D;
	Eigen::Matrix<std::complex<double>,n,n>  V;
	
	
	///////////////////////////////////////////////////////////////////
	// Constructor
	GeneralizedEigenSolver(){
		D.setZero();
		V.setZero();
		
	}
	
	
	///////////////////////////////////////////////////////////////////
	// solve beta*A*x = alpha*B*x
	void solve(const Eigen::Matrix<T,n,n> & MA, const  Eigen::Matrix<T,n,n> & MB){
		
		D.setZero();
		V.setZero();
		
		
//		std::cout<<MA<<std::endl<<std::endl;
//		std::cout<<MB<<std::endl;
		
		Eigen::Matrix<T,n,n> MAT=MA.transpose();
		Eigen::Matrix<T,n,n> MBT=MB.transpose();
		
		// 1 - allocates a workspace for computing eigenvalues of n-by-n real generalized nonsymmetric eigensystems
		gsl_eigen_genv_workspace* w = gsl_eigen_genv_alloc (n);
		
		gsl_matrix_view A = gsl_matrix_view_array (MAT.data(), n, n);
		gsl_matrix_view B = gsl_matrix_view_array (MBT.data(), n, n);
		gsl_vector_complex * alpha = gsl_vector_complex_alloc (n);
		gsl_vector * beta = gsl_vector_alloc (n);
//		gsl_vector_complex * beta = gsl_vector_complex_alloc (n);
		gsl_matrix_complex * evec = gsl_matrix_complex_alloc(n,n);
		
		
		gsl_eigen_genv (&A.matrix, &B.matrix, alpha, beta, evec, w);
		
		
		for (int col = 0; col < n; ++col) {
			gsl_complex alpha_i = gsl_vector_complex_get (alpha, col);
			double beta_i       = gsl_vector_get(beta,col);
			D(col) = std::complex<double>(GSL_REAL(alpha_i)/beta_i,GSL_IMAG(alpha_i)/beta_i);
			
			for (int row = 0; row < n; ++row) {
				gsl_complex Vij = gsl_matrix_complex_get(evec,row,col);
				V(row,col) = std::complex<double>(GSL_REAL(Vij),GSL_IMAG(Vij));
			}
		}
		
		// 2 - frees the memory associated with the workspace w
		gsl_eigen_genv_free (w);
		gsl_vector_complex_free(alpha);
		gsl_vector_free(beta);
		gsl_matrix_complex_free(evec);
	}
	
};

} // namespace model
#endif
