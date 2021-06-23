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
        RationalLatticeDirection(const Rational& _rat, const LatticeDirection<dim>& _dir) ;
        
        /**********************************************************************/
        RationalLatticeDirection(const Rational& _rat, const LatticeVector<dim>& v) ;
        
        /**********************************************************************/
        RationalLatticeDirection(const LatticeVector<dim>& v) ;
        
        /**********************************************************************/
        RationalLatticeDirection(const RationalLatticeDirection<dim>& other) = default;
        RationalLatticeDirection(RationalLatticeDirection<dim>&& other) =default;
        
        VectorDimD cartesian() const;
        
        /**********************************************************************/
        Rational dot(const ReciprocalLatticeVector<dim>& other) const;
        
        /**********************************************************************/
        RationalLatticeDirection<dim> operator*(const long int& scalar) const;
        
        /**********************************************************************/
        RationalLatticeDirection<dim> operator/(const long int& scalar) const;
        
        /**********************************************************************/
        RationalLatticeDirection<dim> operator+(const RationalLatticeDirection<dim>& other) const;

        /**********************************************************************/
        RationalLatticeDirection<dim> operator-(const RationalLatticeDirection<dim>& other) const;
        
        /**********************************************************************/
        RationalLatticeDirection<dim> operator+(const LatticeVector<dim>& other) const;
        
        /**********************************************************************/
        RationalLatticeDirection<dim> operator-(const LatticeVector<dim>& other) const;
        
        /**********************************************************************/
        double squaredNorm() const;
        
    };
    
    template<int dim>
    RationalLatticeDirection<dim> operator*(const long int& scalar, const RationalLatticeDirection<dim>& L);
    
} // end namespace
#endif
