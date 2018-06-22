/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DSCL_h_
#define model_DSCL_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
#include <model/Math/SmithDecomposition.h>
#include <model/LatticeMath/Lattice.h>
#include <model/LatticeMath/RationalMatrix.h>
#include <model/LatticeMath/LLL.h>
#include <model/LatticeMath/RLLL.h>


namespace model
{
    /*!Class template that computes the Displacement Shift Complete Lattice (DSCL)
     * of two parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class DSCL : public Lattice<dim>
    {
        static_assert(dim>0,"dim must be > 0.");
        //        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Lattice<dim> LatticeType;
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,dim,dim> MatrixInt;

        
        IntValueType _sigma;
        
        /**********************************************************************/
        static IntValueType gcd(const IntValueType& a,const IntValueType& b)
        {
            return b>0? gcd(b, a % b) : a;
        }

    public:
        
        const LatticeType& A;
        const LatticeType& B;
        
        /**********************************************************************/
        DSCL(const LatticeType& A_in,
            const LatticeType& B_in,
             const bool& useRLLL=true) :
        /* init */ _sigma(1),
        /* init */ A(A_in),
        /* init */ B(B_in)
        {
            update(useRLLL);
        }

        /**********************************************************************/
        void update(const bool& useRLLL)
        {/*!Suppose that lattice B is obtained from A through a rotaion R,
          * that is B=R*A. Then R=B*inv(A)=B.covBasis()*A.contraBasis().transpose()
          */
            const MatrixDimD R(B.covBasis()*A.contraBasis().transpose());
            
            // Check that R is a proper rotation
            const MatrixDimD RRT=R*R.transpose();
            const double RRTmInorm=(RRT-Eigen::Matrix<double,dim,dim>::Identity()).norm()/Eigen::Matrix<double,dim,dim>::Identity().norm();
            if(RRTmInorm>FLT_EPSILON)
            {
                std::cout<<"R="<<std::endl<<R<<std::endl;
                std::cout<<"R*R^T="<<std::endl<<RRT<<std::endl;
                std::cout<<"norm(R*R^T-I)/norm(I)="<<RRTmInorm<<", tol="<<FLT_EPSILON<<std::endl;
                assert(0 && "R IS NOT ORTHOGONAL.");
            }
            // make sure that C2G is proper
            assert(std::fabs(R.determinant()-1.0) < FLT_EPSILON && "R IS NOT PROPER.");

            
            // Compute the transition matrix T=inv(A)*B
            const MatrixDimD T=A.contraBasis().transpose()*B.covBasis();
            
            // For the two lattices to have coincident sites, R must be a rational matrix
            // Compute the integer matrix P and the integer sigma such that T=P/sigma
            RationalMatrix<dim> rm(T);
            const MatrixInt& P=rm.integerMatrix;
            _sigma=rm.den;
            
            // The integer matrix P can be decomposed as P=X*D*Y using the Smith decomposition.
            // X and Y are unimodular, and D is diagonal with D(k,k) dividing D(k+1,k+1)
            // The decopmosition also computes the matices U and V such that D=U*P*V
            SmithDecomposition<dim> sd(P);
            
            // From T=inv(A)*B=P/sigma=X*D*Y/sigma=X*D*inv(V)/sigma, we have
            // B1*(sigma*I)=A1*D
            // where
            // B1=B*V
            // A1=A*X
            // Since V and X are unimodular matrices, B1 and A1 are new bases
            // of the lattices B and A, respectively. Moreover, since
            // (sigma*I) and D are diagonal, the columns of B1 and A1 are
            // proportional, with rational proportionality factors different for each column.
            // For example, columns "i" read
            // b1_i*sigma=a1_i*D(i,i)
            // Therefore, the i-th primitive vectors of the DSCL is
            // c_i=b1_i*sigma/gcd(sigma,D(i,i))=a1_i*D(i,i)/gcd(sigma,D(i,i))
            // or, in matrix form
            // C=B1*N=A1*M, that is
            // C=B*V*N=A*X*M
            // where M=diag(D(i,i)/gcd(sigma,D(i,i))) and
            //       N=diag(sigma/gcd(sigma,D(i,i))) and
            
            MatrixInt M(MatrixInt::Identity());
            MatrixInt N(MatrixInt::Identity());
            for(int i=0;i<dim;++i)
            {
                const IntValueType& dii=sd.matrixD()(i,i);
                M(i,i)=dii/gcd(_sigma,dii);
                N(i,i)=_sigma/gcd(_sigma,dii);
            }
            
            const MatrixDimD D1=A.covBasis()*sd.matrixX().template cast<double>()*N.template cast<double>().inverse();
            const MatrixDimD D2=B.covBasis()*sd.matrixV().template cast<double>()*M.template cast<double>().inverse();
            
            assert((D1-D2).norm()<FLT_EPSILON && "DSCL calculation failed.");
            
            if(useRLLL)
            {
                this->setLatticeBasis(RLLL(0.5*(D1+D2),0.75).reducedBasis());
            }
            else
            {
                this->setLatticeBasis(0.5*(D1+D2));
            }
        }
        
        /**********************************************************************/
        const IntValueType& sigma() const
        {
            return _sigma;
        }
    };
    
    
}
#endif


