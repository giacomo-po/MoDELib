/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SchurComplementSolver_H_
#define model_SchurComplementSolver_H_

#include <assert.h>
#include <Eigen/Dense>

namespace model {
	
	class SchurComplementSolver  {
	
		/*! \brief Solves the system of equations in block-form:
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
		 * where the matrix \f$K_{qq}\f$ is a symmetric  positive definite matrix and \f$K_{pq}\f$ is a full row-rank matrix.
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
		const Eigen::LLT<Eigen::MatrixXd> llt;
		
		//! the partial-pivot LU decomposition of KPQ*inv(KQQ)*KPQ^T
		const Eigen::PartialPivLU<Eigen::MatrixXd> pplu;
		
		//! The vector of Lagrange multipliers
		const Eigen::VectorXd L;
		
		//! The solution vector
		const Eigen::VectorXd X;
		
		
		SchurComplementSolver(const Eigen::MatrixXd& KQQ, const Eigen::MatrixXd& KPQ, 
		/*          */ const Eigen::VectorXd& Fq,  const Eigen::VectorXd& Fp) : llt(KQQ), 
		/*                                                                   */ pplu(KPQ*llt.solve(KPQ.transpose())), 
		/*                                                                   */ L(pplu.solve(KPQ*llt.solve(Fq)-Fp)),
		/*                                                                   */ X(llt.solve(Fq-KPQ.transpose()*L)){
			/*! Numerically this is implemented as:
			 * \code
			 *	Eigen::LLT<Eigen::MatrixXd> llt(KQQ);
			 *	Eigen::PartialPivLU<Eigen::MatrixXd> pplu(KPQ*llt.solve(KPQ.transpose()));
			 *	Eigen::VectorXd Lp = pplu.solve(KPQ*llt.solve(Fq)-Fp);
			 *	Xq  = llt.solve(Fq-KPQ.transpose()*Lp);
			 * \endcode
			 */		
//			std::cout<<"KQQ.determinant()="<<KQQ.determinant()<<std::endl;
//			std::cout<<"llt.info()="<<llt.info()<<std::endl;
//			std::cout<<KQQ<<std::endl;
			
//			assert(llt.info() && "LLT FAILED.");			
		}
		
	};
	
} // namespace model
#endif
