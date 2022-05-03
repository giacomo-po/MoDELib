/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeDirection_cpp_
#define model_LatticeDirection_cpp_

#include <LatticeModule.h>
namespace model
{

    
    /**********************************************************************/
    template <int dim>
    LatticeDirection<dim>::LatticeDirection(const LatticeVectorType& v) :
    /* base init */ LatticeVectorType(((v.squaredNorm()==0)? v : (v/LatticeGCD<dim>::gcd(v)).eval()),v.lattice)
    {
    }
        
    /**********************************************************************/
    template <int dim>
    template<typename Derived>
    LatticeDirection<dim>::LatticeDirection(const Eigen::MatrixBase<Derived>& v,
                     const Lattice<dim>& lat) :
    //        /* base init */ LatticeGCDType(v),
    //        /* base init */ ReciprocalLatticeVectorType(((v.squaredNorm()==0)? v : (v/this->gCD).eval()),v.covBasis,v.contraBasis)
    //        /* base init */ ReciprocalLatticeVectorType(((v.squaredNorm()==0)? v : (v/this->gCD).eval()),v.lattice)
    /* base init */ LatticeVectorType(((v.squaredNorm()==0)? v : (v/LatticeGCD<dim>::gcd(v)).eval()),lat)
    {
    //            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
    }
        
    /**********************************************************************/
    template <int dim>
    LatticeDirection<dim>::LatticeDirection(const ReciprocalLatticeVectorType& r1,const ReciprocalLatticeVectorType& r2) :
    /* delegating */ LatticeDirection(LatticeVectorType(r1.cross(r2)))
    {
    }
        

    /**********************************************************************/
    template <int dim>
    typename LatticeDirection<dim>::LatticeVectorType LatticeDirection<dim>::snapToLattice(const VectorDimD &dP) const
    { /*!@param[in] dP a cartesian input vector
          *\returns the closest lattice vector corresponding to the projection
          * of dP on *this LatticeDirection.
          */
        const VectorDimD dc(this->cartesian());
        return LatticeVectorType((round(dP.dot(dc) / dc.squaredNorm()) * dc).eval(),
                                 this->lattice);
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeDirection<dim>::VectorDimD LatticeDirection<dim>::snapToDirection(const VectorDimD &dP) const
    {
        VectorDimD dc(this->cartesian());
        const double dNorm(dc.norm());
        assert(dNorm > FLT_EPSILON);
        dc /= dNorm;
        return dP.dot(dc) * dc;
    }
    template struct LatticeDirection<3>;

} // end namespace
#endif
