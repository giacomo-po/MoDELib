/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeDirection_h_
#define model_ReciprocalLatticeDirection_h_

#include <model/LatticeMath/LatticeBase.h>
#include <model/LatticeMath/LatticeGCD.h>
#include <model/LatticeMath/LatticeVector.h>
//#include <model/LatticeMath/LatticeDirection.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>

namespace model
{
    template <int dim>
    class ReciprocalLatticeDirection :
    /* inherits */ public LatticeGCD<dim>,
    /* inherits */ public ReciprocalLatticeVector<dim>
    {
        typedef LatticeGCD<dim> LatticeGCDType;
        typedef LatticeBase<dim> LatticeBaseType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;

        
    public:
        
        ReciprocalLatticeDirection(const ReciprocalLatticeVectorType& v) :
        /* base init */ LatticeGCDType(v),
        /* base init */ ReciprocalLatticeVectorType( (v.squaredNorm()==0)? v : (v/this->gCD).eval() )
        {
//            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
        }
        
        ReciprocalLatticeDirection(const LatticeVectorType& v1,const LatticeVectorType& v2) :
        /* delegating */ ReciprocalLatticeDirection(ReciprocalLatticeVectorType(v1.cross(v2)))
        {
//            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
        }
        
        ReciprocalLatticeDirection(const VectorDimD& d) :
        /* delegating */ ReciprocalLatticeDirection(ReciprocalLatticeVectorType(LatticeBaseType::reciprocalLatticeDirection(d)))
        {
            //            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
        }
        
//        Eigen::Matrix<double,dim,1> cartesian() const
//        {
//            return (this->squaredNorm()>0)? ReciprocalLatticeVectorType::cartesian().normalized() : Eigen::Matrix<double,dim,1>::Zero();
//        }
        
    };
    
} // end namespace
#endif
