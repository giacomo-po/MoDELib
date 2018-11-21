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
#include <MaterialSymmetry.h>
#include <FCClattice.h>
#include <BCClattice.h>
#include <DislocationMobility.h>
#include <GrainBoundaryType.h>

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

#include <Al.h>
#include <Ni.h>
#include <Cu.h>
#include <W.h>

#endif
