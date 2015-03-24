/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePoint_h_
#define model_LatticePoint_h_

#include <model/DislocationDynamics/LatticeMath/LatticeVector.h>

namespace model
{
    template <int dim>
    struct LatticePoint : public LatticeVector<dim>
    {
        
    };
    
} // end namespace
#endif
