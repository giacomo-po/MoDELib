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
                assert(0 && "Solver found Eigen::NumericalIssue.");
                break;
                
            case Eigen::NoConvergence:
                assert(0 && "Solver found Eigen::NoConvergence.");
                break;
                
            case Eigen::InvalidInput:
                assert(0 && "Solver found Eigen::InvalidInput.");
                break;
                
            default:
                assert(0 && "Solver failed.");
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
     */
    
    
    /* SchurComplementSparseSolver: general case ******************************/
    template<typename SparseMatrixType, bool useDirect>
    class SchurComplementSparseSolver{};

    /**************************************************************************/
    /* SchurComplementSparseSolver: DIRECT FACTORIZATION **********************/
    template<typename SparseMatrixType>
    class SchurComplementSparseSolver<SparseMatrixType,true>{
    
        
        //! The vector of Lagrange Multipliers
        Eigen::VectorXd _L;
		
		//! The solution vector
        Eigen::VectorXd _X;
        
        bool solutionIsUpToDate;
        
    public:
        
        typedef  Eigen::SimplicialLLT<SparseMatrixType> SPDsolverType;
        const SPDsolverType solver1;
        
        /* constructor ********************************************************/
        SchurComplementSparseSolver(const SparseMatrixType& KQQ) : solutionIsUpToDate(false), solver1(KQQ){
            //! Computes and stores the LLT decomposition of \f$K_{qq}\f$
            checkSolution(solver1);
        }
        
        /* X ******************************************************************/
        const Eigen::VectorXd& X() const {
            //! a const-reference to the vector of unknowns
            assert(solutionIsUpToDate && "METHOD solve(...) HAS NOT BEEN CALLED YET.");
            return _X;
        };
        
        /* L ******************************************************************/
        const Eigen::VectorXd& L() const {
            //! a const-reference to the vector of Lagrange multipliers
            assert(solutionIsUpToDate && "METHOD solve(...) HAS NOT BEEN CALLED YET.");
            return _L;
        };
        
        /* solve ********************************************************/
        void solve(const Eigen::MatrixXd& KPQ, const Eigen::VectorXd& Fq){            
            //! Using the stored LLT decomposition of \f$K_{qq}\f$:
            
            //!     1- compute and store the solution \f$temp=K_{qq}^{-1}*K_{pq}^T\f$
            const Eigen::MatrixXd temp(solver1.solve(KPQ.transpose()));  // find temp=(Kqq*KpqT)^-1
            checkSolution(solver1);
            //!     2- compute and store the solution \f$temp1=K_{qq}^{-1}*F_{q}\f$
            const Eigen::VectorXd temp1(solver1.solve(Fq));
            checkSolution(solver1);

            // Convert to sparse
//          const SparseMatrixType Kpp((KPQ*temp).sparseView());
//          SPDsolverType ldlt(Kpp);

            //!     3- compute and store the LLT decomposition of \f$K_{pp}=K_{pq}*K_{qq}^{-1}*K_{pq}^T=K_{pq}*temp\f$
            const Eigen::LLT<Eigen::MatrixXd> ldlt(KPQ*temp);
            checkSolution(ldlt);
            //!     4- compute and store the solution \f$L_p=K_{pp}^{-1}*F_{q}\f$
            _L= ldlt.solve(KPQ*temp1);
            checkSolution(ldlt);
            //!     5- compute and store the solution \f$X_q=K_{qq}^{-1}*(F_{q}-K_{pq}^T*L_p)=temp1-temp*L_p\f$
            _X=temp1-temp*_L;
            //      6- update bool "solutionIsUpToDate"
            solutionIsUpToDate=true;
        }
    
    };
    
    
    /**************************************************************************/
    /* SchurComplementSparseSolver: ITERATIVE SOLVER **************************/
    template<typename SparseMatrixType>
    class SchurComplementSparseSolver<SparseMatrixType,false>{
        
    public:
        
        typedef Eigen::ConjugateGradient<SparseMatrixType>  SPDsolverType;
        
        
        const SparseMatrixType& KQQ;
        
        SPDsolverType solver1;
        
        //! The vector of Lagrange Multipliers
        Eigen::VectorXd L;
		
		//! The solution vector
        Eigen::VectorXd X;
        
        /* constructor ********************************************************/
        SchurComplementSparseSolver(const SparseMatrixType& _KQQ) : KQQ(_KQQ), solver1(KQQ){
            checkSolution(solver1);
        }
        
        /* solve ********************************************************/
        void solve(const SparseMatrixType& KPQ, const Eigen::VectorXd& Fq, const double& tol=FLT_EPSILON){
            //!     1- compute and store the solution \f$temp=K_{qq}^{-1}*K_{pq}^T\f$ by Conjugate Gradient
            
            
            
            solver1.setTolerance(tol);
            const Eigen::VectorXd temp1(solver1.solve(Fq));

            const Eigen::VectorXd rhs(KPQ*temp1);
            const double threshold ( tol*tol*rhs.squaredNorm() );
            //const double threshold ( FLT_EPSILON*FLT_EPSILON*rhs.squaredNorm() );

            
            L=Eigen::VectorXd::Zero(KPQ.rows());
            Eigen::VectorXd residual(rhs-KPQ*KQQ*(L.transpose()*KPQ).transpose());
            
            Eigen::VectorXd p(residual);
            
            Eigen::VectorXd gamma(Eigen::VectorXd::Zero(KQQ.rows()));
            
            
            int kttMaxiterations=2*KPQ.rows();
            
            int i = 0;
            while(i<kttMaxiterations){
                const Eigen::Matrix<double,1,Eigen::Dynamic> ptK(p.transpose()*KPQ); // a row-vector
                gamma=solver1.solveWithGuess(ptK.transpose(),gamma);
//                gamma=solver1.solve(ptK.transpose());
                checkSolution(solver1);
                const double residualNorm2(residual.squaredNorm());
                const double den(ptK.transpose().dot(gamma));
                const double alpha( residualNorm2 / den );
                L+= alpha * p;
                residual -= alpha * KPQ*gamma;
                const double newResidualNorm2(residual.squaredNorm());
                if (newResidualNorm2<threshold){
                    break;
                } 
                const double beta(newResidualNorm2 / residualNorm2);
                p = residual + beta * p;
                i++;
                if (i==kttMaxiterations){ // no convergence
                    std::cout<<"No conversion after "<<i<<" iterations"<<std::endl;
                    std::cout<<"newResidualNorm2= "<<newResidualNorm2<<">"<<threshold<<std::endl;
                    assert(0 && "KKT solver did not converge");
                }
            }
            

            


            X=solver1.solve(Fq-(L.transpose()*KPQ).transpose());
            //X=solver1.solveWithGuess(Fq-(L.transpose()*KPQ).transpose(),temp1);
            checkSolution(solver1);

        }
        
    };
    
} // namespace model
#endif


//            assert(i<);

//solver1.setTolerance(tol);
//const Eigen::MatrixXd temp(solver1.solve(KPQ.transpose()));
//checkSolution(solver1);
//!     2- compute and store the solution \f$temp1=K_{qq}^{-1}*F_{q}\f$ by Conjugate Gradient
//const Eigen::VectorXd temp1(solver1.solve(Fq));
//checkSolution(solver1);
//!     3- store the (sparse) matrix \f$K_{pp}=K_{pq}*K_{qq}^{-1}*K_{pq}^T=K_{pq}*temp\f$ and dense vector \f$temp2=K_{pq}*K_{qq}^{-1}*F_q=K_{pq}*temp1\f$
//            const SparseMatrixType Kpp((KPQ*temp).sparseView());
//            const Eigen::VectorXd temp2(KPQ*temp1);
//            //!     4- compute and store the solution \f$L_p=K_{pp}^{-1}*F_{q}\f$ by Conjugate Gradient
//            SPDsolverType ldlt;
//            ldlt.setTolerance(tol);
//            L= ldlt.compute(Kpp).solve(temp2);
//            checkSolution(ldlt);
//            //!     5- compute and store the solution \f$X_q=K_{qq}^{-1}*(F_{q}-K_{pq}^T*L_p)=temp1-temp*L_p\f$
//            X=temp1-temp*L; // ok


//            solver1.setTolerance(tol);
