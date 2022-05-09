/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_LatticeCore_h_
#define model_LatticeCore_h_

#include <iostream>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <IntegerMath.h>


namespace model
{
    template <int dim>
    struct LatticeCore
    {
        static_assert(dim>0,"dim must be > 0.");
        static constexpr double roundTol=FLT_EPSILON;
        
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef long long int IntScalarType;
        typedef Eigen::Matrix<IntScalarType,dim,1> VectorDimI;
        typedef Eigen::Matrix<IntScalarType,dim,dim> MatrixDimI;
        
        static VectorDimI rationalApproximation(VectorDimD v);
        static VectorDimI integerCoordinates(const VectorDimD& d,const MatrixDimD& invA);
    };
}
#endif
