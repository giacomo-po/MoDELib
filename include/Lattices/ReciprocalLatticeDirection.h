/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeDirection_h_
#define model_ReciprocalLatticeDirection_h_

//#include <LatticeBase.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeGCD.h>

namespace model
{
    template <int dim>
    struct ReciprocalLatticeDirection :
    /* inherits */ public ReciprocalLatticeVector<dim>
    {
        typedef LatticeGCD<dim> LatticeGCDType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        /**********************************************************************/
        ReciprocalLatticeDirection(const ReciprocalLatticeVectorType& v) ;
        
        /**********************************************************************/
        template<typename Derived>
        ReciprocalLatticeDirection(const Eigen::MatrixBase<Derived>& v,
                                   const Lattice<dim>& lat) :
        /* base init */ ReciprocalLatticeVectorType(((v.squaredNorm()==0)? v : (v/LatticeGCD<dim>::gcd(v)).eval()),lat)
        {
            //            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
        }
        
        /**********************************************************************/
        ReciprocalLatticeDirection(const LatticeVectorType& v1,const LatticeVectorType& v2);
        
    };
    
} // end namespace
#endif
