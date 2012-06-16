// This file is part of mmdl, a C++ template library
// for network dislocation dynamics.
//
// Copyright (C) 2009 Giacomo Po <giacomopo@gmail.com>,
//
//
// mmdl is distributed without any warranty under the 
// GNU Lesser General Public License <http://www.gnu.org/licenses/>.


#ifndef mmdl_LAGRANGESPARSESOLVER_H_
#define mmdl_LAGRANGESPARSESOLVER_H_

#include <Eigen/Dense>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <Eigen/Sparse>
#include <Eigen/SparseExtra>
//#include <Eigen/CholmodSupport>
//#include <Eigen/UmfPackSupport>

namespace mmdl {
	
	class LagrangeSparseSolver  {
		
		/*! \brief Solves the block-matrix system of equations:
		 *	\f[
		 *		\left[
		 *		\begin{array}{ll}
		 *		K_{qq} & K_{pq}^T\\
		 *		K_{pq} & 0\\
		 *		\end{array}
		 *		\right]
		 *		\left[
		 *		\begin{array}{l}
		 *		X_q\\
		 *		L_p
		 *		\end{array}
		 *		\right]
		 *		=
		 *		\left[
		 *		\begin{array}{l}
		 *		F_q\\
		 *		F_p
		 *		\end{array}
		 *		\right]
		 *	\f]
		 * where the matrix Kqq is symmetric and positive definite. 
		 *
		 * The full system is symmetric but not necessarely poisitive-definite, therefore Cholesky decomposition (LLT) cannot be used. 
		 * However, since \f$K_{qq}\f$ is symmetric and positive-definite, LLT can be used for the decomposition of \f$K_{qq}\f$. 
		 * In order to take advantage of LLT and reduce the size of the overall matrix decomposition, the solution is found in the form:
		 *	\f[ \left\{
		 *		\begin{array}{l}
		 *		L_p=\left(K_{pq}K_{qq}^{-1}K_{pq}^T\right)^{-1}\left(K_{pq}K_{qq}^{-1}F_q-F_p\right)\\
		 *		X_q=K_{qq}^{-1}\left(F_q-K_{pq}^TL_p\right)
		 *		\end{array}
		 *		\right.
		 *	\f]
		 */
		
	public:
		//! The LLT decomposition of KQQ
		//const Eigen::LLT<Eigen::MatrixXd> llt;
		//const Eigen::SparseLLT<SparseMatrix<double>, Eigen::Cholmod>  llt;
		//Eigen::SparseLLT<Eigen::SparseMatrix<double> > llt;
		Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > llt;
		//, SparseLDLT<SparseMatrixType>
		//const Eigen::SparseLU<SparseMatrix<double>,Eigen::UmfPack> llt;
		
		//! the partial-pivot LU decomposition of KPQ*inv(KQQ)*KPQ^T
		const Eigen::SparseLU<Eigen::SparseMatrix<double> > pplu;
		//const Eigen::SparseLU<SparseMatrix<double>,Eigen::UmfPack> pplu;
		//const Eigen::PartialPivLU<Eigen::MatrixXd> pplu;
		
		//! The vector of Lagrange multipliers
		const Eigen::VectorXd L;
		
		//! The solution vector
		 Eigen::VectorXd X;
		Eigen::MatrixXd Y;
		
		 
		
		//Eigen::SparseMatrix<double> XXXX;
		
		Eigen::VectorXd XXXX;
		
		LagrangeSparseSolver(const Eigen::SparseMatrix<double>& KQQ, const Eigen::SparseMatrix<double>& KPQ, 
		/*                */ const Eigen::VectorXd& Fq,  const Eigen::VectorXd& Fp) : llt(KQQ)
		//pplu()
		{
			
			Eigen::SparseMatrix<double> KQP (50,50);
			Y.resize(50,50);
			
			
			llt.solve(KQP);
			
//			int  = ;
//			int m = ;
//			
//			// Finally, the fastest way to fill a SparseMatrix object is to insert the elements in purely increasing order 
//			// (increasing inner index per outer index, and increasing outer index) using the insertBack() function:
//			Eigen::SparseMatrix<double> A(N+M,N+M);
//			//mat.reserve(estimated_number_of_non_zero);  // optional
//			for(int j=0; j<N+M; ++j)
//			{
//				A.startVec(j);                          // optional for a DynamicSparseMatrix
//				for each i interacting with j             // with increasing i
//					A.insertBack(i,j) = foo(i,j);
//			}
//			A.finalize();                             // optional for a DynamicSparseMatrix
			
			
			
			//KPQ*llt.solve(KPQ.transpose());
			//: llt(KQQ), 
		///*                                                                   */ pplu(KPQ*llt.solve(KPQ.transpose())), 
		///*                                                                   */ L(pplu.solve(KPQ*llt.solve(Fq)-Fp)),
		///*                                                                   */ X(llt.solve(Fq-KPQ.transpose()*L)){
			/*! Numerically this is implemented as:
			 * \code
			 *	Eigen::LLT<Eigen::MatrixXd> llt(KQQ);
			 *	Eigen::PartialPivLU<Eigen::MatrixXd> pplu(KPQ*llt.solve(KPQ.transpose()));
			 *	Eigen::VectorXd Lp = pplu.solve(KPQ*llt.solve(Fq)-Fp);
			 *	Xq  = llt.solve(Fq-KPQ.transpose()*Lp);
			 * \endcode
			 */
		}
		
	};
	
} // namespace mmdl
#endif
