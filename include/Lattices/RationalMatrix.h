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
//#include <RoundEigen.h>
#include <BestRationalApproximation.h>
#include <iostream>

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
        static IntValueType gcd(const IntValueType& a,const IntValueType& b);
        
        /**********************************************************************/
        static IntValueType lcm(const IntValueType& a,const IntValueType& b);
        
        /**********************************************************************/
        static std::pair<MatrixInt,IntValueType> compute(const MatrixDimD& R);
        
        const std::pair<MatrixInt,IntValueType> returnPair;
        
    public:
        
        const MatrixInt& integerMatrix;
        const IntValueType& sigma;
        
        /**********************************************************************/
        RationalMatrix(const MatrixDimD& R) ;
    };
    
    
} // end namespace
#endif


