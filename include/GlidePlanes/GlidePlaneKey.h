/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpF@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneKey_H
#define model_GlidePlaneKey_H

#include <array>
#include <Eigen/Dense>
#include <LatticeModule.h>


namespace model
{
    
    template <int dim>
    using GlidePlaneKey = LatticePlaneKey<dim>;
    
}
#endif

