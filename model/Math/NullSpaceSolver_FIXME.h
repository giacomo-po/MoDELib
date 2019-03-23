
#ifndef _model_NullSpaceSolver_h_
#define _model_NullSpaceSolver_h_

#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <SparseNullSpace.h>

//#ifdef _MODEL_PARDISO_SOLVER_
#include <Eigen/PardisoSupport>
//#endif

namespace model
{

    template<typename SparseMatrixType>
    class NullSpaceSolver
    {
        
        typedef typename SparseMatrixType::Scalar Scalar;
        
        SparseMatrixType CY;
        SparseMatrixType ZTAZ;

        Eigen::PardisoLLT<SparseMatrixType> solver1;
        Eigen::PardisoLLT<SparseMatrixType> solver2;
        
    public:
        
        /********************************************************/
        NullSpaceSolver(const SparseMatrixType& A,
                        const SparseMatrixType& C)
        {
            SparseNullSpace<SparseMatrixType> sn(C);

            CY=C*sn.matrixY();
            solver1.factorize(CY);
            
            ZTAZ=sn.matrixY().transpose()*A*sn.matrixZ();
            solver2.factorize(ZTAZ);
        }
        
        /********************************************************/
        Eigen::VectorXd solve(const Eigen::VectorXd& b, const Eigen::VectorXd& d) const
        {/*! Solves for the current constraint matrix
          */
            const auto t1= std::chrono::system_clock::now();
//#ifdef _MODEL_PARDISO_SOLVER_
//            model::cout<<"Solving (PardisoLLT)..."<<std::flush;
            const Eigen::VectorXd YxY=sn.matrixY()*solver1.solve(d);
            return YxY+ns.matrixZ()*solver2.solve(ns.matrixZ().transpose()*(b-A*YxY));

//#else
//            model::cout<<"Solving (ConjugateGradient)..."<<std::flush;
//            Eigen::ConjugateGradient<SparseMatrixType> solver(CY);
//            //                    Eigen::ConjugateGradient<SparseMatrixType> solver(T.transpose()*A*T); // this gives a segmentation fault
//            solver.setTolerance(tol);
//            //                    x=T*solver.solve(T.transpose()*(b-A*g))+g;
//            x=T*solver.solveWithGuess(b1,guess)+g;
//            //                    x=T*solver.solve(b1)+g;
//            model::cout<<" (relative error ="<<solver.error()<<", tolerance="<<solver.tolerance();
//#endif
        
            
        }

    };
    
}

#endif
