/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_BiCrystal_cpp_
#define model_BiCrystal_cpp_

#include <LatticeModule.h>

namespace model
{
        
        template <int dim>
        typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getM(const RationalMatrix<dim>& rm,
                            const SmithDecomposition<dim>& sd)
        {
            typename BiCrystal<dim>::MatrixDimI M=BiCrystal<dim>::MatrixDimI::Identity();
            for(int i=0;i<dim;++i)
            {
                const auto& dii(sd.matrixD()(i,i));
                M(i,i)=dii/IntegerMath<IntScalarType>::gcd(rm.mu,dii);
            }
            return M;
        }
        
        template <int dim>
        typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getN(const RationalMatrix<dim>& rm,
                            const SmithDecomposition<dim>& sd)
        {
            typename BiCrystal<dim>::MatrixDimI N=BiCrystal<dim>::MatrixDimI::Identity();
            for(int i=0;i<dim;++i)
            {
                const auto& dii=sd.matrixD()(i,i);
                N(i,i)=rm.mu/IntegerMath<IntScalarType>::gcd(rm.mu,dii);
            }
            return N;
        }
        
        template <int dim>
        typename BiCrystal<dim>::MatrixDimD BiCrystal<dim>::getCSLBasis(const Lattice<dim>& A,
                                                                                 const Lattice<dim>& B,
                                                                                 const SmithDecomposition<dim>& sd,
                                                                                 const typename BiCrystal<dim>::MatrixDimI& M,
                                                                                 const typename BiCrystal<dim>::MatrixDimI& N,
                                                                                 const bool& useRLLL)
        {
            // The transition matrix is T=P/sigma, where P=rm.integerMatrix is
            // an integer matrix and sigma=rm.sigma is an integer
            // The integer matrix P can be decomposed as P=X*D*Y using the Smith decomposition.
            // X and Y are unimodular, and D is diagonal with D(k,k) dividing D(k+1,k+1)
            // The decomposition also computes the matices U and V such that D=U*P*V
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
            // Therefore, the i-th primitive vectors of the CSL is
            // c_i=b1_i*sigma/gcd(sigma,D(i,i))=a1_i*D(i,i)/gcd(sigma,D(i,i))
            // or, in matrix form
            // C=B1*N=A1*M, that is
            // C=B*V*N=A*X*M
            // where M=diag(D(i,i)/gcd(sigma,D(i,i))) and
            //       N=diag(sigma/gcd(sigma,D(i,i))) and

            const auto C1(A.latticeBasis*(sd.matrixX()*M).template cast<double>());
            const auto C2(B.latticeBasis*(sd.matrixV()*N).template cast<double>());
            if ((C1-C2).norm()>FLT_EPSILON)
            {
                throw std::runtime_error("CSL calculation failed.\n");
            }

            if(useRLLL)
            {
                return RLLL(0.5*(C1+C2),0.75).reducedBasis();
            }
            else
            {
                return 0.5*(C1+C2);
            }
        }
        
        template <int dim>
        typename BiCrystal<dim>::MatrixDimD BiCrystal<dim>::getDSCLBasis(const Lattice<dim>& A,
                                       const Lattice<dim>& B,
                                       const SmithDecomposition<dim>& sd,
                                       const typename BiCrystal<dim>::MatrixDimI& M,
                                       const typename BiCrystal<dim>::MatrixDimI& N,
                                       const bool& useRLLL)
        {
                        
            const auto D1(A.latticeBasis*sd.matrixX().template cast<double>()*N.template cast<double>().inverse());
            const auto D2(B.latticeBasis*sd.matrixV().template cast<double>()*M.template cast<double>().inverse());
            if ((D1-D2).norm()>FLT_EPSILON)
            {
                throw std::runtime_error("DSCL calculation failed.\n");
            }
            if(useRLLL)
            {
                return RLLL(0.5*(D1+D2),0.75).reducedBasis();
            }
            else
            {
                return 0.5*(D1+D2);
            }
        }
        
        
        template <int dim>
        BiCrystal<dim>::BiCrystal(const Lattice<dim>& A_in,
                  const Lattice<dim>& B_in,
                  const bool& useRLLL) :
        /* init */ RationalMatrix<dim>(A_in.reciprocalBasis.transpose()*B_in.latticeBasis)
        /* init */,SmithDecomposition<dim>(this->integerMatrix)
        /* init */,A(A_in)
        /* init */,B(B_in)
        /* init */,M(getM(*this,*this))
        /* init */,N(getN(*this,*this))
        /* init */,sigmaA(M.determinant())
        /* init */,sigmaB(N.determinant())
        /* init */,sigma(std::abs(sigmaA)==std::abs(sigmaB)? std::abs(sigmaA) : 0)
        /* init */, csl(getCSLBasis (A,B,*this,M,N,useRLLL),MatrixDimD::Identity())
        /* init */,dscl(getDSCLBasis(A,B,*this,M,N,useRLLL),MatrixDimD::Identity())
        /* init */,Ap(A.latticeBasis*this->matrixX().template cast<double>())
        /* init */,Bp(B.latticeBasis*this->matrixV().template cast<double>())
        {
            
            if(true)
            {//verify that CSL can be obtained as multiple of A and B
                
                const MatrixDimD tempA(A.reciprocalBasis.transpose()*csl.latticeBasis);
                if ((tempA-tempA.array().round().matrix()).norm()>FLT_EPSILON)
                {
                    throw std::runtime_error("CSL is not a multiple of lattice A\n");
                }
                
                const MatrixDimD tempB(B.reciprocalBasis.transpose()*csl.latticeBasis);
                if ((tempB-tempB.array().round().matrix()).norm()>FLT_EPSILON)
                {
                    throw std::runtime_error("CSL is not a multiple of lattice B\n");
                }
            }
            
            if(true)
            {//verify that A and B are multiples of DSCL
                
                const MatrixDimD tempA(dscl.reciprocalBasis.transpose()*A.latticeBasis);
                if ((tempA-tempA.array().round().matrix()).norm()>FLT_EPSILON)
                {
                    throw std::runtime_error("Lattice A is not a multiple of the DSCL\n");
                }
                
                const MatrixDimD tempB(dscl.reciprocalBasis.transpose()*B.latticeBasis);
                if ((tempB-tempB.array().round().matrix()).norm()>FLT_EPSILON)
                {
                    throw std::runtime_error("Lattice B is not a multiple of the DSCL\n");
                }
            }

            //            update(useRLLL);
        }
        
//        template<int dim>
//        LatticeDirection<dim> BiCrystal<dim>::AtoCSLvector(const LatticeVector<dim>& n) const
//        {
//            return csl.latticeDirection(M.adjoint()*this->matrixU()*n); // NOOO, adjoint in Eigen means conjugate transpose
//        }
    
    
//    template class BiCrystal<1>;
    template class BiCrystal<2>;
    template class BiCrystal<3>;

} // end namespace
#endif

