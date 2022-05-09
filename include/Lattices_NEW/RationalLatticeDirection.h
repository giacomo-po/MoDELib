/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_RationalLatticeDirection_h_
#define model_RationalLatticeDirection_h_

#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <Rational.h>
#include <BestRationalApproximation.h>

namespace model
{
    template <int dim>
    struct RationalLatticeDirection
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        const Rational<IntScalarType> rat;
        const LatticeDirection<dim> dir;
        
    public:
        
        
        RationalLatticeDirection(const Rational<IntScalarType>& _rat, const LatticeDirection<dim>& _dir) ;
        RationalLatticeDirection(const Rational<IntScalarType>& _rat, const LatticeVector<dim>& v) ;
        RationalLatticeDirection(const LatticeVector<dim>& v) ;
        RationalLatticeDirection(const RationalLatticeDirection<dim>& other) = default;
        RationalLatticeDirection(RationalLatticeDirection<dim>&& other) =default;
        VectorDimD cartesian() const;
        Rational<IntScalarType> dot(const ReciprocalLatticeVector<dim>& other) const;
        RationalLatticeDirection<dim> operator*(const IntScalarType& scalar) const;
        RationalLatticeDirection<dim> operator/(const IntScalarType& scalar) const;
        RationalLatticeDirection<dim> operator+(const RationalLatticeDirection<dim>& other) const;
        RationalLatticeDirection<dim> operator-(const RationalLatticeDirection<dim>& other) const;
        RationalLatticeDirection<dim> operator+(const LatticeVector<dim>& other) const;
        RationalLatticeDirection<dim> operator-(const LatticeVector<dim>& other) const;
        double squaredNorm() const;
        
    };
    
    template<int dim>
    RationalLatticeDirection<dim> operator*(const typename RationalLatticeDirection<dim>::IntScalarType& scalar, const RationalLatticeDirection<dim>& L);
    
} // end namespace
#endif
