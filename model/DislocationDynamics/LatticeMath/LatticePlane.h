/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <model/DislocationDynamics/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/LatticeMath/LatticeDirection.h>

namespace model
{
    struct LatticePlane
    {
        
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3> LatticeDirectionType;
        
        LatticeVectorType P;
        LatticeDirectionType n;
        
        LatticePlane(const LatticeVectorType& P_in,const LatticeVectorType n_in) :
        /* init */ P(P_in),
        /* init */ n(n_in)
        {}
        
        
        LatticePlane(const LatticeVectorType& P_in,const LatticeDirectionType n_in) :
        /* init */ P(P_in),
        /* init */ n(n_in)
        {}
        
    };
    
} // end namespace
#endif
