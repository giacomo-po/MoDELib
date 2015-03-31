/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeDirection.h>

namespace model
{
    struct LatticePlane
    {
        
        typedef LatticeVector<3>    LatticeVectorType;
        typedef ReciprocalLatticeVector<3> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<3> ReciprocalLatticeDirectionType;
        
        LatticeVectorType P;
        ReciprocalLatticeDirectionType n;
        
        LatticePlane(const LatticeVectorType& P_in,const ReciprocalLatticeVectorType n_in) :
        /* init */ P(P_in),
        /* init */ n(n_in)
        {}
        
        
        LatticePlane(const LatticeVectorType& P_in,const ReciprocalLatticeDirectionType n_in) :
        /* init */ P(P_in),
        /* init */ n(n_in)
        {}
        
//        /**********************************************************************/
//        PlaneLineIntersection intersectWith(const LatticeLine& line) const
//        {
//            return PlaneLineIntersection(*this,line);
//        }
        
    };
    
} // end namespace
#endif
