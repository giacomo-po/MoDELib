/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlaneBase.h>

namespace model
{
    struct LatticePlane
    {
        
        typedef LatticeVector<3>    LatticeVectorType;
        
        const LatticeVectorType P;
        const LatticePlaneBase& n;
        
        LatticePlane(const LatticeVectorType& P_in,const LatticePlaneBase& n_in) :
        /* init */ P(P_in),
        /* init */ n(n_in)
        {}
        
        Eigen::Matrix<double,3,1> snapToLattice(const Eigen::Matrix<double,3,1>& P0) const
        {
            return P.cartesian()+n.snapToLattice(P0-P.cartesian());
        }
        
        bool contains(const LatticeVectorType& L) const
        {
            return (L-P).dot(n)==0;
        }
        
    };
    
} // end namespace
#endif
