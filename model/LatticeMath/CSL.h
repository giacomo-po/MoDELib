/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CSL_h_
#define model_CSL_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
#include <model/Math/SmithDecomposition.h>
#include <model/LatticeMath/Lattice.h>
#include <model/LatticeMath/RationalRotation.h>


namespace model
{
    /*!Class template that computes the Coincident Site Lattice (CSL) of two 
     * parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class CSL : public Lattice<dim>
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
        CSL(const LatticeType& A_in,
            const LatticeType& B_in) :
        /* init */ _sigma(1),
        /* init */ A(A_in),
        /* init */ B(B_in)
        {
            update();
        }
        
        /**********************************************************************/
        void update()
        {
            // Suppose that lattice B is obtained from A through a rotaion R, that is B=R*A
            // The rotation matrix is R=B*inv(A)
            const MatrixDimD R(B.covBasis()*A.contraBasis().transpose());

            // For the two lattices to have coincident sites, R must be a rational matrix
            // Compute the integer matrix P and the integer sigma such that R=P/sigma
            RationalRotation<dim> rr(R);
            const MatrixInt& P=rr.integerMatrix();
            _sigma=rr.sigma();
            
            // The integer matrix P can be decomposed as P=X*D*Y using the Smith decomposition.
            // X and Y are unimodular, and D is diagonal with D(k,k) dividing D(k+1,k+1)
            // The decopmosition also computes the matices U and V such that D=U*P*V
            SmithDecomposition<dim> sd(P);

            // From R=B*inv(A)=P/sigma=X*D*Y/sigma=inv(U)*D*inv(V)/sigma, we have
            // U*B*inv(A)*V=U*B*inv(A)*inv(Y)=(U*B)*inv(Y*A)=D/sigma
            // Since U and Y are unimodular matrices, they define a change in basis,
            // and therefore the bases B1=U*B and A1=Y*A satisfy
            // B1*inv(A1)=D/sigma
            // Since D is diagonal, B1 and A1 are simply scaled bases, with different
            // scaling in different directions. Therefore their primitive vectors
            // (columns of B1 and A1) satisfy
            // sigma*b1_i=D(i,i)*a1_i
            // Therefore, the primitive vectors of the CSL are
            // c_i=sigma/gcd(sigma,D(i,i))*b1_i=D(i,i)/gcd(sigma,D(i,i))*a1_i
            MatrixInt M(MatrixInt::Identity());
            MatrixInt N(MatrixInt::Identity());
            for(int i=0;i<dim;++i)
            {
                const IntValueType& dii=sd.matrixD()(i,i);
                M(i,i)=dii/gcd(_sigma,dii);
                N(i,i)=_sigma/gcd(_sigma,dii);
            }
            
            const MatrixDimD C1=(M*sd.matrixY()).template cast<double>()*A.covBasis();
            const MatrixDimD C2=(N*sd.matrixU()).template cast<double>()*B.covBasis();
            assert((C1-C2).norm()<FLT_EPSILON && "CSL calculation failed.");
            this->setLatticeBasis(0.5*(C1+C2));
        }

        /**********************************************************************/
        const IntValueType& sigma() const
        {
            return _sigma;
        }
    };
    
    
} // end namespace
#endif


