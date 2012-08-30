/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SchurComplementSparseSolver_H_
#define model_SchurComplementSparseSolver_H_
#include <assert.h>
#include <float.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/CholmodSupport>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#include <model/Utilities/AreSameType.h>

namespace model {

    
    /* constructor ********************************************************/
    template <typename EigenSolverType>
    void checkSolution(const EigenSolverType& solver) {
        switch (solver.info()) {
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
    
    /*! \brief Solves the sparse system of equations in block-form:
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
    
//    template <typename SparseMatrixType, template<typename T, int UL> class SPDsolverType>
//    template <typename SPDsolverType,typename SPsDsolverType>
    
    
    /* SchurComplementSparseSolver: general case ******************************/
    template<typename SparseMatrixType, bool useDirect>
    struct SchurComplementSparseSolver{};

    /**************************************************************************/
    /* SchurComplementSparseSolver: DIRECT FACTORIZATION **********************/
    template<typename SparseMatrixType>
    struct SchurComplementSparseSolver<SparseMatrixType,true>{
    
        typedef  Eigen::SimplicialLLT<SparseMatrixType> SPDsolverType;
        typedef Eigen::SimplicialLDLT<SparseMatrixType> SPsDsolverType;
        
        const SPDsolverType llt;
        
        //! The vector of Lagrange Multipliers
        Eigen::VectorXd L;
		
		//! The solution vector
        Eigen::VectorXd X;
        
        /* constructor ********************************************************/
        SchurComplementSparseSolver(const SparseMatrixType& KQQ) : llt(KQQ){
            checkSolution(llt);
        }
        
        /* solve ********************************************************/ // DO NOT ERASE, WAIT FOR NEXT EIGEN
        void solve(const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq){            
            
            const SparseMatrixType KPQT(KPQ.transpose());  // store the tranpose
            
            const SparseMatrixType temp(llt.solve(KPQT));  // find temp=(Kqq*KpqT)^-1
            checkSolution(llt);
            const Eigen::VectorXd temp1(llt.solve(Fq));
            checkSolution(llt);
            
            const SPsDsolverType ldlt((KPQ*temp).pruned(FLT_EPSILON)); // decompose temp
            checkSolution(ldlt);
            const Eigen::VectorXd temp2(KPQ*temp1);
            L= ldlt.solve(temp2);
            checkSolution(ldlt);
            X=temp1-temp*L;
//            X=llt.solve(Fq-KPQT*L);
//            checkSolution<(llt);
        }
        
    
    };
    
    
    /**************************************************************************/
    /* SchurComplementSparseSolver: ITERATIVE SOLVER **************************/
    template<typename SparseMatrixType>
    struct SchurComplementSparseSolver<SparseMatrixType,false>{
        
        typedef Eigen::ConjugateGradient<SparseMatrixType>  SPDsolverType;
        typedef SPDsolverType          SPsDsolverType;

        
        SPDsolverType llt;
        
        //! The vector of Lagrange Multipliers
        Eigen::VectorXd L;
		
		//! The solution vector
        Eigen::VectorXd X;
        
        /* constructor ********************************************************/
        SchurComplementSparseSolver(const SparseMatrixType& KQQ) : llt(KQQ){
            checkSolution(llt);
        }
        
        /* solve ********************************************************/
        void solve(const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq, const double& tol=FLT_EPSILON){
 
            const SparseMatrixType KPQT(KPQ.transpose());  // store the tranpose
            
            llt.setTolerance(tol);
            const SparseMatrixType temp(llt.solve(KPQT));  // find temp=(Kqq*KpqT)^-1
            checkSolution(llt);
            const Eigen::VectorXd temp1(llt.solve(Fq));
            checkSolution(llt);
            
            
            const SparseMatrixType temp3(KPQ*temp);   // without this temporary, next "compute" goes crazy            
            SPsDsolverType ldlt;//((KPQ*temp).pruned(FLT_EPSILON)); // decompose temp
            ldlt.setTolerance(tol);

            ldlt.compute(temp3); // decompose temp
            checkSolution(ldlt);
            //ldlt.setTolerance(tol);

            const Eigen::VectorXd temp2(KPQ*temp1);
            
            L= ldlt.solve(temp2);
            checkSolution(ldlt);
            X=temp1-temp*L;
            //            X=llt.solve(Fq-KPQT*L);
            //            checkSolution<(llt);
        }
        
    };
    
  
	
} // namespace model
#endif










/**/
/**/
/**/
/**/
/**/

//    template <typename SPDsolverType>
////	struct SchurComplementSparseSolver : public SPDsolverType, public SPsDsolverType {
//        struct SchurComplementSparseSolver : public SPDsolverType {
//
//            typedef SPDsolverType SPDsolverType;
//
//        static_assert(AreSameType<typename SPDsolverType::MatrixType,typename SPsDsolverType::MatrixType>::value, "INCONSISTENT SPD and SPsD SOLVERS ARE USED.");
//
//        typedef typename SPDsolverType::MatrixType SparseMatrixType;
//
//
//       // const SPDsolverType spdSolver;
//
//
//   //     const SparseMatrixType& KQQ; // a const reference to the KQQ matrix
//
//		//! The vector of Lagrange multipliers
//        Eigen::VectorXd L;
//
//		//! The solution vector
//        Eigen::VectorXd X;
//
//
//
//
//        /* constructor ********************************************************/
//        SchurComplementSparseSolver(const SparseMatrixType& KQQ) {
//            SPDsolverType::compute(KQQ);
//            checkSolution<SPDsolverType>(*this);
//        }
//
//
//
//
//
//        /* solve ********************************************************/
//        void solve(const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq){
//            std::cout<<"I'm here 1"<<std::endl;
//
//            const SparseMatrixType KPQT(KPQ.transpose());
//            std::cout<<"I'm here 2"<<std::endl;
//
//            const SparseMatrixType temp(SPDsolverType::solve(KPQT));
//            checkSolution<SPDsolverType>(*this);
//            std::cout<<"I'm here 3"<<std::endl;
//
//            SPsDsolverType::compute((KPQ*temp).pruned(FLT_EPSILON));
//            checkSolution<SPsDsolverType>(*this);
//
//            std::cout<<"I'm here 4"<<std::endl;
//
////            const SPsDsolverType sqSolver((KPQ*temp).pruned(FLT_EPSILON));
//            //checkSolution(sqSolver);
//            const Eigen::VectorXd temp1(SPDsolverType::solve(Fq));
//            checkSolution<SPDsolverType>(*this);
//
//            std::cout<<"I'm here 5"<<std::endl;
//
////            L= sqSolver.solve(KPQ*temp1);
//            const Eigen::VectorXd temp2(KPQ*temp1);
//            std::cout<<"temp2="<<std::endl<<temp2<<std::endl;
//            L= SPsDsolverType::solve(temp2);
//            checkSolution<SPsDsolverType>(*this);
//
//            //L= SPsDsolverType::solve(KPQ*temp1);
//            std::cout<<"I'm here 6"<<std::endl;
//
//            X=SPDsolverType::solve(Fq-KPQT*L);
//            std::cout<<"I'm here 7"<<std::endl;
//            checkSolution<SPDsolverType>(*this);
//
//
//
//        }





//	};


//        /* constructor ********************************************************/
//        SchurComplementSparseSolver(const SparseMatrixType& KQQ_in) : spdSolver(KQQ) {
//            //checkSolution(llt);
//
//        }

//        void compute(){
//            SPDsolverType::compute(KQQ);
//            //            //checkSolution(llt);
//
//        }


//        /* directSolve ********************************************************/
//        void iterativeSolve(const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq){
//            const SparseMatrixType KPQT(KPQ.transpose());
//            const SparseMatrixType temp(spdSolver.solve(KPQT));
//            const Eigen::BiCGSTAB<SparseMatrixType> sqSolver((KPQ*temp).pruned(FLT_EPSILON));
//            //checkSolution(sqSolver);
//            const Eigen::VectorXd temp1(spdSolver.solve(Fq));
//            L= sqSolver.solve(KPQ*temp1);
//            X=spdSolver.solve(Fq-KPQT*L);
//        }


//        SchurComplementSparseSolver(const SparseMatrixType& KQQ, const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq) {
//
//            //checkSolution(llt);
//            const SparseMatrixType KPQT(KPQ.transpose());
//            const SparseMatrixType temp(spdSolver.solve(KPQT));
//            const SQsolverType<SparseMatrixType> sqSolver((KPQ*temp).pruned(FLT_EPSILON));
//            //checkSolution(ldlt);
//            const Eigen::VectorXd temp1(spdSolver.solve(Fq));
//            L= sqSolver.solve(KPQ*temp1);
//            X=spdSolver.solve(Fq-KPQT*L);
//		}




//        SchurComplementSparseSolver(const SparseMatrixType& KQQ, const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq) {
//
//            const Eigen::SimplicialLLT <SparseMatrixType> llt(KQQ);
//            checkSolution(llt);
//            const SparseMatrixType KPQT(KPQ.transpose());
//            const SparseMatrixType temp(llt.solve(KPQT));
//            const Eigen::SimplicialLDLT<SparseMatrixType> ldlt((KPQ*temp).pruned(FLT_EPSILON));
//            checkSolution(ldlt);
//            const Eigen::VectorXd temp1(llt.solve(Fq));
//            L=ldlt.solve(KPQ*temp1);
//            X=llt.solve(Fq-KPQT*L);
//		}
//
//        SchurComplementSparseSolver(const SparseMatrixType& KQQ, const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq) {
//
////            const Eigen::SimplicialLLT <SparseMatrixType> llt;
//            const Eigen::ConjugateGradient<SparseMatrixType> cg(KQQ);
//
//           // checkSolution(llt);
//            const SparseMatrixType KPQT(KPQ.transpose());
//            const SparseMatrixType temp(cg.solve(KPQT));
//            const Eigen::BiCGSTAB<SparseMatrixType> bic((KPQ*temp).pruned(FLT_EPSILON));
//           // checkSolution(ldlt);
//            const Eigen::VectorXd temp1(cg.solve(Fq));
//            L=bic.solve(KPQ*temp1);
//            X=cg.solve(Fq-KPQT*L);
//		}