/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_RationalLatticeDirection_cpp_
#define model_RationalLatticeDirection_cpp_

#include <LatticeModule.h>
#include <RationalLatticeDirection.h>

namespace model
{
    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim>::RationalLatticeDirection(const Rational<IntScalarType> &_rat, const LatticeDirection<dim> &_dir) :
/* init */ rat(_rat)
/* init */,dir(_dir)
    {
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim>::RationalLatticeDirection(const Rational<IntScalarType> &_rat, const LatticeVector<dim> &v) :
    /* init */ RationalLatticeDirection(_rat, LatticeDirection<dim>(v))
    {
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim>::RationalLatticeDirection(const LatticeVector<dim> &v) :
    /* init */ rat(Rational<IntScalarType>(IntegerMath<IntScalarType>::gcd(v), 1)),
    /* init */ dir(v)
    {
    }

    /**********************************************************************/

    template <int dim>
    typename RationalLatticeDirection<dim>::VectorDimD RationalLatticeDirection<dim>::cartesian() const
    {
        return dir.cartesian() * rat.asDouble();
    }

    /**********************************************************************/
    template <int dim>
    Rational<typename RationalLatticeDirection<dim>::IntScalarType> RationalLatticeDirection<dim>::dot(const ReciprocalLatticeVector<dim> &other) const
    {
        return rat * dir.dot(other);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator*(const IntScalarType &scalar) const
    {
        return RationalLatticeDirection<dim>(rat * scalar, dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator/(const IntScalarType &scalar) const
    {
        return RationalLatticeDirection<dim>(rat / scalar, dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator+(const RationalLatticeDirection<dim> &other) const
    {
        assert(&dir.lattice == &other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const VectorDimI temp(rat.n * other.rat.d * dir + other.rat.n * rat.d * other.dir);
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
        const LatticeVector<dim> v(temp / gcd, dir.lattice);
        return RationalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), LatticeDirection<dim>(v));
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator-(const RationalLatticeDirection<dim> &other) const
    {
        assert(&dir.lattice == &other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const VectorDimI temp(rat.n * other.rat.d * dir - other.rat.n * rat.d * other.dir);
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
        const LatticeVector<dim> v(temp / gcd, dir.lattice);
        return RationalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), LatticeDirection<dim>(v));
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator+(const LatticeVector<dim> &other) const
    {
        assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
        return this->operator+(RationalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), LatticeDirection<dim>(other)));
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator-(const LatticeVector<dim> &other) const
    {
        assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
        return this->operator-(RationalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), LatticeDirection<dim>(other)));
    }

    /**********************************************************************/
    template <int dim>
    double RationalLatticeDirection<dim>::squaredNorm() const
    {
        return dir.squaredNorm() * std::pow(rat.asDouble(), 2);
    }
    
    template<int dim>
    RationalLatticeDirection<dim> operator*(const typename RationalLatticeDirection<dim>::IntScalarType& scalar, const RationalLatticeDirection<dim>& L)
    {
        return L*scalar;
    }
    
    template struct RationalLatticeDirection<1>;
    template RationalLatticeDirection<1> operator*(const typename RationalLatticeDirection<1>::IntScalarType& scalar, const RationalLatticeDirection<1>& L);
    template struct RationalLatticeDirection<2>;
    template RationalLatticeDirection<2> operator*(const typename RationalLatticeDirection<2>::IntScalarType& scalar, const RationalLatticeDirection<2>& L);
    template struct RationalLatticeDirection<3>;
    template RationalLatticeDirection<3> operator*(const typename RationalLatticeDirection<3>::IntScalarType& scalar, const RationalLatticeDirection<3>& L);

} // end namespace
#endif
