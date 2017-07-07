/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_H_
#define model_DislocationLoop_H_

#include <model/LoopNetwork/Loop.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlaneBase.h>
#include <model/LatticeMath/LatticePlane.h>

namespace model
{
    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    class DislocationLoop : public Loop<DislocationLoop<dim,corder,InterpolationType,QuadratureRule>>
    {
    
        typedef DislocationLoop<dim,corder,InterpolationType,QuadratureRule> DislocationLoopType;
        typedef Loop<DislocationLoopType> LoopType;
        typedef typename TypeTraits<DislocationLoopType>::FlowType FlowType;
        
        
    public:
        
        const LatticePlane glidePlane;
        
        /**********************************************************************/
        DislocationLoop(const FlowType& flow,
                        const LatticeVector<dim>& P,
                        const LatticePlaneBase& N) :
        /* base init */ LoopType(flow),
        /* base init */ glidePlane(P,N)
        {
        
        }
    
    };
    
}

#endif
