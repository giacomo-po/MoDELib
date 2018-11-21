/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_H_
#define model_SLIPSYSTEM_H_

#include <assert.h>
#include <LatticePlaneBase.h>
#include <LatticeVector.h>

namespace model
{
    
    struct SlipSystem
    {
        

        const LatticePlaneBase n;
        const LatticeVector<3>  s;
        const Eigen::Matrix<double,3,1>  unitNormal;
        
        
        SlipSystem(const LatticeVector<3>& a1,
                   const LatticeVector<3>& a2,
                   const LatticeVector<3>& slip_in):
        /* init list */ n(a1,a2),
        /* init list */ s(slip_in),
        /* init list */ unitNormal(n.cartesian().normalized())
        {
            assert(std::fabs(n.dot(s))==0 && "PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL.");
        }
        
    };

}
#endif
