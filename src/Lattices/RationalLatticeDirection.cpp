/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RationalLatticeDirection_cpp_
#define model_RationalLatticeDirection_cpp_

#include <LatticeModule.h>
#include <RationalLatticeDirection.h>

namespace model
{
    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim>::RationalLatticeDirection(const Rational &_rat, const LatticeDirection<dim> &_dir) : /* init */ rat(_rat)
                                                                                        /* init */,
                                                                                        dir(_dir)
    {
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim>::RationalLatticeDirection(const Rational &_rat, const LatticeVector<dim> &v) : /* init */ RationalLatticeDirection(_rat, LatticeDirection<dim>(v))
    {
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim>::RationalLatticeDirection(const LatticeVector<dim> &v) : /* init */ rat(Rational(LatticeGCD<dim>::gcd(v), 1))
                                                            /* init */,
                                                            dir(v)
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
    Rational RationalLatticeDirection<dim>::dot(const ReciprocalLatticeVector<dim> &other) const
    {
        return rat * dir.dot(other);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator*(const long int &scalar) const
    {
        return RationalLatticeDirection<dim>(rat * scalar, dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator/(const long int &scalar) const
    {
        return RationalLatticeDirection<dim>(rat / scalar, dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator+(const RationalLatticeDirection<dim> &other) const
    {
        assert(&dir.lattice == &other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const VectorDimI temp(rat.n * other.rat.d * dir + other.rat.n * rat.d * other.dir);
        const long int gcd(LatticeGCD<dim>::gcd(temp));
        const LatticeVector<dim> v(temp / gcd, dir.lattice);
        return RationalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), LatticeDirection<dim>(v));
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator-(const RationalLatticeDirection<dim> &other) const
    {
        assert(&dir.lattice == &other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const VectorDimI temp(rat.n * other.rat.d * dir - other.rat.n * rat.d * other.dir);
        const long int gcd(LatticeGCD<dim>::gcd(temp));
        const LatticeVector<dim> v(temp / gcd, dir.lattice);
        return RationalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), LatticeDirection<dim>(v));
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator+(const LatticeVector<dim> &other) const
    {
        assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const long int gcd(LatticeGCD<dim>::gcd(other));
        return this->operator+(RationalLatticeDirection<dim>(Rational(gcd, 1), LatticeDirection<dim>(other)));
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator-(const LatticeVector<dim> &other) const
    {
        assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const long int gcd(LatticeGCD<dim>::gcd(other));
        return this->operator-(RationalLatticeDirection<dim>(Rational(gcd, 1), LatticeDirection<dim>(other)));
    }

    /**********************************************************************/
    template <int dim>
    double RationalLatticeDirection<dim>::squaredNorm() const
    {
        return dir.squaredNorm() * std::pow(rat.asDouble(), 2);
    }
    
    template<int dim>
    RationalLatticeDirection<dim> operator*(const long int& scalar, const RationalLatticeDirection<dim>& L)
    {
        return L*scalar;
    }
    
    template struct RationalLatticeDirection<3>;
    template RationalLatticeDirection<3> operator*(const long int& scalar, const RationalLatticeDirection<3>& L);

} // end namespace
#endif
