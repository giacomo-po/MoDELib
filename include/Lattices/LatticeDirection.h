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
        LatticeDirection(const LatticeVectorType& v) ;
        /**********************************************************************/
        template<typename Derived>
        LatticeDirection(const Eigen::MatrixBase<Derived>& v,
                         const Lattice<dim>& lat) ;
        
        /**********************************************************************/
        LatticeDirection(const ReciprocalLatticeVectorType& r1,const ReciprocalLatticeVectorType& r2);
        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const VectorDimD& dP) const;
        
        /**********************************************************************/
        VectorDimD snapToDirection(const VectorDimD& dP) const;
        
    };
    
} // end namespace
#endif
