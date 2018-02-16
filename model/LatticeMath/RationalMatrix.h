/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RationalMatrix_h_
#define model_RationalMatrix_h_

#include <iomanip>
#include <cfloat> // FLT_EPSILON
#include <assert.h> // FLT_EPSILON
#include <utility>
#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
#include <model/Math/BestRationalApproximation.h>

namespace model
{
    
    template <int dim>
    class RationalMatrix
    {
        static_assert(dim>0,"dim must be > 0.");
        typedef long long int IntValueType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Eigen::Matrix<int,dim,dim> MatrixDimI;
        typedef Eigen::Matrix<IntValueType,dim,dim> MatrixInt;

        static constexpr int64_t maxDen=100000;
        
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
        
        /**********************************************************************/
        static std::pair<MatrixInt,IntValueType> compute(const MatrixDimD& R)
        {
                        
            // Find the BestRationalApproximation of each entry
            MatrixDimI nums(MatrixDimI::Zero());
            MatrixDimI dens(MatrixDimI::Ones());
            
            IntValueType sigma=1;
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    BestRationalApproximation bra(R(i,j),maxDen);
                    nums(i,j)=bra.num;
                    dens(i,j)=bra.den;
                    sigma=lcm(sigma,bra.den);
                }
            }
            
            MatrixInt im(MatrixInt::Zero());
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    im(i,j)=nums(i,j)*sigma/dens(i,j);
                }
            }
            
            const double error=(im.template cast<double>()/sigma-R).norm()/(dim*dim);
            if(error>100.0*DBL_EPSILON)
            {
                std::cout<<"error="<<error<<std::endl;
                std::cout<<"maxDen="<<maxDen<<std::endl;
                std::cout<<"im=\n"<<std::setprecision(15)<<std::scientific<<im.template cast<double>()/sigma<<std::endl;
                std::cout<<"= 1/"<<sigma<<"*\n"<<std::setprecision(15)<<std::scientific<<im<<std::endl;
                std::cout<<"R=\n"<<std::setprecision(15)<<std::scientific<<R<<std::endl;
                assert(false && "Rational Matrix failed, check maxDen.");
            }
            

            
            return std::make_pair(im,sigma);
        }
        
//        IntValueType _sigma;
//        MatrixInt im;
        
        const std::pair<MatrixInt,IntValueType> returnPair;
        
    public:
        
        const MatrixInt& integerMatrix;
        const IntValueType& den;
        
        /**********************************************************************/
        RationalMatrix(const MatrixDimD& R) :
        /* init */ returnPair(compute(R)),
        /* init */ integerMatrix(returnPair.first),
        /* init */ den(returnPair.second)
        {
            
        }
        
//        /**********************************************************************/
//        const IntValueType& sigma() const
//        {
//            return _sigma;
//        }
//        
//        /**********************************************************************/
//        const MatrixInt& integerMatrix() const
//        {
//            return im;
//        }
//        
    };
    
    
} // end namespace
#endif


