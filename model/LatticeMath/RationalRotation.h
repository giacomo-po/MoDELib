/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RationalRotation_h_
#define model_RationalRotation_h_

#include <cfloat> // FLT_EPSILON
#include <assert.h> // FLT_EPSILON
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
#include <model/Math/BestRationalApproximation.h>

//#include <model/Math/SmithDecomposition.h>
//#include <model/Math/BestRationalApproximation.h>
//#include <model/LatticeMath/Lattice.h>


namespace model
{
    
    template <int dim>
    class RationalRotation
    {
        static_assert(dim>0,"dim must be > 0.");
//        //        static constexpr double roundTol=FLT_EPSILON;
//        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
//        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef long long int IntValueType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Eigen::Matrix<int,dim,dim> MatrixDimI;
        typedef Eigen::Matrix<IntValueType,dim,dim> MatrixInt;

        static constexpr int maxDen=1000;
        
        /**********************************************************************/
        static IntValueType gcd(const IntValueType& a,const IntValueType& b)
        {
            return b>0? gcd(b, a % b) : a;
        }
        
        /**********************************************************************/
        static IntValueType lcm(const IntValueType& a,const IntValueType& b)
        {
            return a*b/gcd(a,b);
        }
        
        IntValueType _sigma;
        MatrixInt im;
        
    public:
        
//        const LatticeType& parentLattice1;
//        const LatticeType& parentLattice2;
        
        /**********************************************************************/
        RationalRotation(const MatrixDimD& R) :
        /* init */ _sigma(1),
        /* init */ im(MatrixInt::Identity())
        {
            assert((R*R.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && "R IS NOT ORTHOGONAL.");
            // make sure that C2G is proper
            assert(std::fabs(R.determinant()-1.0) < FLT_EPSILON && "R IS NOT PROPER.");
         
            
            // Find the BestRationalApproximation of each entry
            MatrixDimI nums(MatrixDimI::Zero());
            MatrixDimI dens(MatrixDimI::Ones());
            
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    BestRationalApproximation bra(R(i,j),maxDen);
                    nums(i,j)=bra.num;
                    dens(i,j)=bra.den;
                    _sigma=lcm(_sigma,bra.den);
                }
            }
            
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    im(i,j)=nums(i,j)*_sigma/dens(i,j);
                }
            }
         
            assert((im.template cast<double>()/_sigma-R).norm()<2.0*DBL_EPSILON*dim*dim && "Rational Matrix failed, check maxDen.");
            
        }
        
        /**********************************************************************/
        const IntValueType& sigma() const
        {
            return _sigma;
        }
        
        /**********************************************************************/
        const MatrixInt& integerMatrix() const
        {
            return im;
        }
        
    };
    
    
} // end namespace
#endif


