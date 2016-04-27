/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeDirection_h_
#define model_LatticeDirection_h_

#include <model/LatticeMath/LatticeBase.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>

namespace model
{
    template <int dim>
    class LatticeDirection :
    /* inherits */ public LatticeGCD<dim>,
    /* inherits */ public LatticeVector<dim>
    {
        typedef LatticeGCD<dim> LatticeGCDType;
        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
    public:
        
        LatticeDirection(const LatticeVectorType& v) :
        /* base init */ LatticeGCDType(v),
        /* base init */ LatticeVectorType((v.squaredNorm()==0)? v : (v/this->gCD).eval())
        {
//            assert(this->squaredNorm() && "LatticeDirection has Zero Norm");
        }
        
        LatticeDirection(const ReciprocalLatticeVectorType& r1,const ReciprocalLatticeVectorType& r2) :
        /* delegating */ LatticeDirection(LatticeVectorType(r1.cross(r2)))
        {
//            assert(this->squaredNorm() && "LatticeDirection has Zero Norm");
        }
        
        LatticeDirection(const VectorDimD& d) :
        /* delegating */ LatticeDirection(LatticeVectorType(LatticeBaseType::latticeDirection(d)))
        {
            //            assert(this->squaredNorm() && "LatticeDirection has Zero Norm");
        }
        
//        LatticeDirection(const VectorDimD& d):
//        /* delegating */ LatticeDirection(LatticeVectorType(d))
//        {}
        
        /**********************************************************************/
        VectorDimD snapToDirection(const VectorDimD& dP) const
        {
            const VectorDimD dc(this->cartesian());
            return round(dP.dot(dc)/dc.squaredNorm())*dc;
        }
        
    };
    
} // end namespace
#endif
