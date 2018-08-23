/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicElement_H_
#define model_PeriodicElement_H_

//#include <list>
#include <deque>
#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
#include <model/DislocationDynamics/Materials/FCClattice.h>
#include <model/DislocationDynamics/Materials/BCClattice.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int Z, typename SymmetryType>
    struct PeriodicElement
    {
        // static_assert(0,"PERIODIC ELEMENT NOT IMPLEMENTED");
    };
    
}

#include <model/DislocationDynamics/Materials/Elements/Al.h>
#include <model/DislocationDynamics/Materials/Elements/Ni.h>
#include <model/DislocationDynamics/Materials/Elements/Cu.h>
#include <model/DislocationDynamics/Materials/Elements/W.h>

#endif
