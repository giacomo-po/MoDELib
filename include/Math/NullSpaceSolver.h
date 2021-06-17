
#ifndef _model_NullSpaceSolver_h_
#define _model_NullSpaceSolver_h_

#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <SparseNullSpace.h>

#ifdef _MODEL_PARDISO_SOLVER_
#include <Eigen/PardisoSupport>
#endif

namespace model
{
    
    template<typename SparseMatrixType>
    class NullSpaceSolver
    {
        
        typedef typename SparseMatrixType::Scalar Scalar;
        
        
        const SparseNullSpace<SparseMatrixType> sNs;
        const SparseMatrixType KY;
        //        SparseMatrixType CY;
        //        SparseMatrixType ZTKZ;
        
        Eigen::VectorXd _y;
        Eigen::VectorXd _z;
        

        
#ifdef _MODEL_PARDISO_SOLVER_
        Eigen::PardisoLLT<SparseMatrixType> solver1;
        Eigen::PardisoLLT<SparseMatrixType> solver2;
#else
        Eigen::SimplicialLLT<SparseMatrixType> solver1;
        Eigen::SimplicialLLT<SparseMatrixType> solver2;
#endif
        
    public:
        
        /********************************************************/
        NullSpaceSolver(const SparseMatrixType& K,
                        const SparseMatrixType& C) :
        /* init */ sNs(C)
        /* init */,KY(K*sNs.matrixY())
        {
            
            
            
            //            CY=C*sNs.matrixY();
            //            ZTKZ=sNs.matrixZ().transpose()*K*sNs.matrixZ();
            //
            //            solver1.factorize(CY);
            //            solver2.factorize(ZTKZ);
            
            solver1.factorize(C*sNs.matrixY());
            solver2.factorize(sNs.matrixZ().transpose()*K*sNs.matrixZ());
            
        }
        
        /********************************************************/
        Eigen::VectorXd solve(const Eigen::VectorXd& f, const Eigen::VectorXd& g)
        {/*! Solves for the current constraint matrix
          */            _y=solver1.solve(g);
            _z=solver2.solve((sNs.matrixZ().transpose()*(f-KY*_y)).eval());
            return sNs.matrixY()*_y+sNs.matrixZ()*_z;
        }
        
        /********************************************************/
        Eigen::VectorXd solve(const Eigen::VectorXd& f)
        {/*! Solves for the current constraint matrix
          */
            _y.setZero(sNs.matrixY().cols());
            _z=solver2.solve((sNs.matrixZ().transpose()*f).eval());
            return sNs.matrixZ()*_z;
        }
        
        /********************************************************/
        const Eigen::VectorXd& y() const
        {/*!\returns the vector y which is used to compute to solution as
          * x=Y*y+Z*z
          */
            return _y;
        }
        
        /********************************************************/
        const Eigen::VectorXd& z() const
        {/*!\returns the vector z which is used to compute to solution as
          * x=Y*y+Z*z
          */
            return _z;
        }
        
        const SparseMatrixType& matrixZ() const
        {/*!\returns the matrix Z which is used to compute to solution as
          * x=Y*y+Z*z
          * The matrix Z is the null-space of the constaint matrix C
          */
            return sNs.matrixZ();
        }
        
        const SparseMatrixType& matrixY() const
        {/*!\returns the matrix Y which is used to compute to solution as
          * x=Y*y+Z*z
          */
            return sNs.matrixY();
        }
    };
    
}

#endif
