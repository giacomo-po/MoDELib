/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeDirection_h_
#define model_LatticeDirection_h_

//#include <LatticeBase.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeGCD.h>

namespace model
{
    template <int dim>
    struct LatticeDirection :
//    /* inherits */ public LatticeGCD<dim>,
    /* inherits */ public LatticeVector<dim>
    {
        typedef LatticeGCD<dim> LatticeGCDType;
//        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
//    public:
        
        /**********************************************************************/
        LatticeDirection(const LatticeVectorType& v) :
        /* base init */ LatticeVectorType(((v.squaredNorm()==0)? v : (v/LatticeGCD<dim>::gcd(v)).eval()),v.lattice)
        {
        }
        
        /**********************************************************************/
        template<typename Derived>
        LatticeDirection(const Eigen::MatrixBase<Derived>& v,
                         const Lattice<dim>& lat) :
        //        /* base init */ LatticeGCDType(v),
        //        /* base init */ ReciprocalLatticeVectorType(((v.squaredNorm()==0)? v : (v/this->gCD).eval()),v.covBasis,v.contraBasis)
        //        /* base init */ ReciprocalLatticeVectorType(((v.squaredNorm()==0)? v : (v/this->gCD).eval()),v.lattice)
        /* base init */ LatticeVectorType(((v.squaredNorm()==0)? v : (v/LatticeGCD<dim>::gcd(v)).eval()),lat)
        {
            //            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
        }
        
        /**********************************************************************/
        LatticeDirection(const ReciprocalLatticeVectorType& r1,const ReciprocalLatticeVectorType& r2) :
        /* delegating */ LatticeDirection(LatticeVectorType(r1.cross(r2)))
        {
        }
        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const VectorDimD& dP) const
        {/*!@param[in] dP a cartesian input vector
          *\returns the closest lattice vector corresponding to the projection
          * of dP on *this LatticeDirection.
          */
            const VectorDimD dc(this->cartesian());
            return LatticeVectorType((round(dP.dot(dc)/dc.squaredNorm())*dc).eval(),
                                     this->lattice);
        }
        
        /**********************************************************************/
        VectorDimD snapToDirection(const VectorDimD& dP) const
        {
            VectorDimD dc(this->cartesian());
            const double dNorm(dc.norm());
            assert(dNorm>FLT_EPSILON);
            dc/=dNorm;
            return dP.dot(dc)*dc;
        }
        
    };
    
} // end namespace
#endif

//        LatticeDirection(const VectorDimD& d) :
//        /* delegating */ LatticeDirection(LatticeVectorType(LatticeBaseType::latticeDirection(d)))
//        {
//            //            assert(this->squaredNorm() && "LatticeDirection has Zero Norm");
//        }

//        LatticeDirection(const VectorDimD& d):
//        /* delegating */ LatticeDirection(LatticeVectorType(d))
//        {}

//        /**********************************************************************/
//        VectorDimD snapToDirection(const VectorDimD& dP) const
//        {
//            const VectorDimD dc(this->cartesian());
//            return round(dP.dot(dc)/dc.squaredNorm())*dc;
//        }
