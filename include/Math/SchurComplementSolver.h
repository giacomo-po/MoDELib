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
     *
     * Numerically this is implemented as:
     * \code
     *	Eigen::LLT<Eigen::MatrixXd> llt(KQQ);
     *	Eigen::PartialPivLU<Eigen::MatrixXd> pplu(KPQ*llt.solve(KPQ.transpose()));
     *	Eigen::VectorXd Lp = pplu.solve(KPQ*llt.solve(Fq)-Fp);
     *	Xq  = llt.solve(Fq-KPQ.transpose()*Lp);
     * \endcode
     */	
	class SchurComplementSolver  {
        
        
        /* constructor ********************************************************/
        void checkDecomposition() const {
            switch (llt.info()) {
                case Eigen::Success:
                    break;
                    
                case Eigen::NumericalIssue:
                    assert(0 && "SchurComplementSolver found Eigen::NumericalIssue.");
                    break;
                    
                case Eigen::NoConvergence:
                    assert(0 && "SchurComplementSolver found Eigen::NoConvergence.");
                    break;
                    
                case Eigen::InvalidInput:
                    assert(0 && "SchurComplementSolver found Eigen::InvalidInput.");
                    break;
                    
                default:
                    assert(0 && "SchurComplementSolver failed.");
                    break;
            }        
        }
	
	public:
		//! The LLT decomposition of KQQ
		const Eigen::LLT<Eigen::MatrixXd> llt;
		
		//! the partial-pivot LU decomposition of KPQ*inv(KQQ)*KPQ^T
		const Eigen::PartialPivLU<Eigen::MatrixXd> pplu;
		
		//! The vector of Lagrange multipliers
        Eigen::VectorXd L;
		
		//! The solution vector
        Eigen::VectorXd X;
		
        /* constructor with Fp ************************************************/
		SchurComplementSolver(const Eigen::MatrixXd& KQQ, 
        /*                 */ const Eigen::MatrixXd& KPQ, 
		/*                 */ const Eigen::VectorXd& Fq,  
        /*                 */ const Eigen::VectorXd& Fp) : 
        /* init list                                    */ llt(KQQ), 
		/* init list                                    */ pplu(KPQ*llt.solve(KPQ.transpose())){
            checkDecomposition();
            //solve(KPQ,Fq,Fp);
            L=pplu.solve(KPQ*llt.solve(Fq)-Fp);
            X=llt.solve(Fq-KPQ.transpose()*L);
		}
        
        /* constructor without Fp *********************************************/
		SchurComplementSolver(const Eigen::MatrixXd& KQQ, 
        /*                 */ const Eigen::MatrixXd& KPQ, 
        /*                 */ const Eigen::VectorXd&  Fq) : 
        /* init list                                     */ llt(KQQ), 
		/* init list                                     */ pplu(KPQ*llt.solve(KPQ.transpose())){
            checkDecomposition();
            //solve(KPQ,Fq);
            L=pplu.solve(KPQ*llt.solve(Fq));
            X=llt.solve(Fq-KPQ.transpose()*L);
		}
        
//        /* constructor without Fp *********************************************/
//		SchurComplementSolver(const Eigen::MatrixXd& KQQ, 
//        /*                 */ const Eigen::MatrixXd& KPQ) : 
//        /* init list                                     */ llt(KQQ), 
//		/* init list                                     */ pplu(KPQ*llt.solve(KPQ.transpose())){
//            checkDecomposition();
//		}

        
//        /* constructor without Fp *********************************************/
//		void solve(const Eigen::VectorXd& Fq) {
//            L=pplu.solve(KPQ*llt.solve(Fq));
//            X=llt.solve(Fq-KPQ.transpose()*L);
//		}
//            
//        /* constructor without Fp *********************************************/
//        void solve(const Eigen::VectorXd& KPQ, const Eigen::VectorXd& Fq, const Eigen::VectorXd& Fp) {
//                L=pplu.solve(KPQ*llt.solve(Fq)-Fp);
//                X=llt.solve(Fq-KPQ.transpose()*L);
//        }

        

		
	};
	
} // namespace model
#endif
