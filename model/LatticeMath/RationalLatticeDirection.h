/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RationalLatticeDirection_h_
#define model_RationalLatticeDirection_h_

//#include <LatticeBase.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeGCD.h>
#include <Rational.h>
#include <BestRationalApproximation.h>
//#include <ReciprocalLatticeVector.h>

namespace model
{
    template <int dim>
    struct RationalLatticeDirection
    {
//        typedef LatticeGCD<dim> LatticeGCDType;
//        typedef LatticeBase<dim> LatticeBaseType;
//        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        
        const Rational rat;
        const LatticeDirection<dim> dir;
        
    public:
        
        /**********************************************************************/
        RationalLatticeDirection(const Rational& _rat, const LatticeDirection<dim>& _dir) :
        /* init */ rat(_rat)
        /* init */,dir(_dir)
        {
        }
        
        /**********************************************************************/
        RationalLatticeDirection(const Rational& _rat, const LatticeVector<dim>& v) :
        /* init */ RationalLatticeDirection(_rat,LatticeDirection<dim>(v))
        {
        }
        
        /**********************************************************************/
        RationalLatticeDirection(const LatticeVector<dim>& v) :
        /* init */ rat(Rational(LatticeGCD<dim>::gcd(v),1))
        /* init */,dir(v)
        {
        }
        
        /**********************************************************************/
        RationalLatticeDirection(const RationalLatticeDirection& other) = default;
        RationalLatticeDirection(RationalLatticeDirection&& other) =default;
        
        VectorDimD cartesian() const
        {
            return dir.cartesian()*rat.asDouble();
        }
        
        /**********************************************************************/
        Rational dot(const ReciprocalLatticeVector<dim>& other) const
        {
            return rat*dir.dot(other);
        }
        
        /**********************************************************************/
        RationalLatticeDirection operator*(const long int& scalar) const
        {
            return RationalLatticeDirection(rat*scalar,dir);
        }
        
        /**********************************************************************/
        RationalLatticeDirection operator/(const long int& scalar) const
        {
            return RationalLatticeDirection(rat/scalar,dir);
        }
        
        /**********************************************************************/
        RationalLatticeDirection operator+(const RationalLatticeDirection& other) const
        {
            assert(&dir.lattice==&other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
            const VectorDimI temp(rat.n*other.rat.d*dir+other.rat.n*rat.d*other.dir);
            const long int gcd(LatticeGCD<dim>::gcd(temp));
            const LatticeVector<dim> v(temp/gcd,dir.lattice);
            return RationalLatticeDirection(Rational(gcd,rat.d*other.rat.d),LatticeDirection<dim>(v));
        }

        /**********************************************************************/
        RationalLatticeDirection operator-(const RationalLatticeDirection& other) const
        {
            assert(&dir.lattice==&other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
            const VectorDimI temp(rat.n*other.rat.d*dir-other.rat.n*rat.d*other.dir);
            const long int gcd(LatticeGCD<dim>::gcd(temp));
            const LatticeVector<dim> v(temp/gcd,dir.lattice);
            return RationalLatticeDirection(Rational(gcd,rat.d*other.rat.d),LatticeDirection<dim>(v));
        }
        
        /**********************************************************************/
        RationalLatticeDirection operator+(const LatticeVector<dim>& other) const
        {
            assert(&dir.lattice==&other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
            const long int gcd(LatticeGCD<dim>::gcd(other));
            return this->operator+(RationalLatticeDirection(Rational(gcd,1),LatticeDirection<dim>(other)));
        }
        
        /**********************************************************************/
        RationalLatticeDirection operator-(const LatticeVector<dim>& other) const
        {
            assert(&dir.lattice==&other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
            const long int gcd(LatticeGCD<dim>::gcd(other));
            return this->operator-(RationalLatticeDirection(Rational(gcd,1),LatticeDirection<dim>(other)));
        }
        
        /**********************************************************************/
        double squaredNorm() const
        {
            return dir.squaredNorm()*std::pow(rat.asDouble(),2);
        }
        
    };
    
    template<int dim>
    RationalLatticeDirection<dim> operator*(const long int& scalar, const RationalLatticeDirection<dim>& L)
    {
        return L*scalar;
    }
    
} // end namespace
#endif
