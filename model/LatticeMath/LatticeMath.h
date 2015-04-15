/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeMath_
#define model_LatticeMath_

namespace model
{
    template <int dim>
    class LatticeVector;

    template <int dim>
    class ReciprocalLatticeVector;

    class LatticePlane;
    class LatticeLine;
    class PlaneLineIntersection;

} // end namespace

#include <model/LatticeMath/LatticeGCD.h>
#include <model/LatticeMath/LatticeBase.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticeDirection.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeDirection.h>

#include <model/LatticeMath/LatticePlaneBase.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/LatticeMath/LatticeLine.h>
#include <model/LatticeMath/PlaneLineIntersection.h>
#include <model/LatticeMath/PlanePlaneIntersection.h>
#include <model/LatticeMath/LineLineIntersection.h>

#endif
