/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_MINRES_H_
#define model_MINRES_H_
#include <assert.h>
#include <float.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace model {
    
    /**************************************************************************/
    template<typename _ScalarType>
    class MINRES{
        /*! \brief The MINimum RESidual norm solver for sparse symmetric (indefinite) problems.
         *
         *  Adapted from B. Fischer, "Polynomial Based Iteration Methods for Symmetric Linear Systems", SIAM 2011
         */
        
    public:
        
       // typedef Eigen::ConjugateGradient<SparseMatrixType>  SPDsolverType;
        typedef _ScalarType ScalarType;
        
        typedef Eigen::SparseMatrix<ScalarType> SparseMatrixType;
        typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> DenseVectorType;
        
		
		//! The solution vector
        DenseVectorType xMR;
        ScalarType norm_rMR;
        
        DenseVectorType xOR;        
        ScalarType norm_rOR;
        
        /* constructor ********************************************************/
        MINRES(const SparseMatrixType& A, const DenseVectorType& f, const DenseVectorType& x0, const ScalarType& tol){
            
            // initialize
            int n=0;
            const int N(A.rows());
            const int n_max(10*N);
            assert(A.cols()==N && "A is not square");
            assert(f.rows()==N && "A and F have different # of rows");
            assert(x0.rows()==N && "A and x0 have different # of rows");
            DenseVectorType v(DenseVectorType::Zero(N));
            DenseVectorType v_hat(f-A*x0);
            ScalarType beta(v_hat.norm());
            ScalarType c(1.0); // the cosine of the Givens rotation
            ScalarType c_old(1.0); 
            ScalarType s(0.0); // the sine of the Givens rotation
            ScalarType s_old(0.0); // the sine of the Givens rotation
            DenseVectorType w(DenseVectorType::Zero(N));
            DenseVectorType w_old(w);
            ScalarType eta(beta);
            xMR=x0;
            norm_rMR=beta;
            const ScalarType norm_r0(beta);

//            while ( n<n_max && norm_rMR/norm_r0 > tol){
                while ( n<n_max ){
                
                
                // Lanczos
                DenseVectorType v_old(v);
                v=v_hat/beta;
                DenseVectorType Av(A*v);
                ScalarType alpha(v.transpose()*Av);
                v_hat=Av-alpha*v-beta*v_old;
                ScalarType beta_old(beta);
                beta=v_hat.norm();
                
                // QR
                ScalarType c_oold(c_old);
                c_old=c;
                ScalarType s_oold(s_old);
                s_old=s;
                ScalarType r1_hat=c_old *alpha-c_oold*s_old *beta_old;
                    ScalarType r1 =std::pow(std::pow(r1_hat,2)+std::pow(beta,2),0.5);
                ScalarType r2 =s_old *alpha+c_oold*c_old*beta_old;
                ScalarType r3 =s_oold*beta_old;
                
                // Givens rotation
                c=r1_hat/r1;
                s=beta/r1;
                
                // update
                DenseVectorType w_oold(w_old);
                w_old=w;
                w=(v-r3*w_oold-r2*w_old) /r1;
                xMR=xMR+c*eta*w;
                norm_rMR=norm_rMR*std::fabs(s);
                eta=-s*eta;
                    
                    if(norm_rMR/norm_r0 < tol){
                        break;
                    }
                
                n++;
                    
                    if(n==n_max){
                        std::cout<<"MINRES iteration "<< n << "of "<<n_max<<": realtive error = "<< norm_rMR/norm_r0 << ", tol="<<tol<<std::endl;
                        assert(0 && "MINRES DID NOT CONVERGE.");
                    }
            
            }
            if (std::fabs(c)>DBL_EPSILON){
                xOR =xMR-s*eta*w/c;
                norm_rOR=norm_rMR/abs(c);
            }
            
        }
        

        
    };
    
} // namespace model
#endif

