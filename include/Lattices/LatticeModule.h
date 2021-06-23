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
    class Lattice;
    
    template <int dim>
    class LatticeVector;

    template <int dim>
    class ReciprocalLatticeVector;

    template <int dim>
    struct LatticeDirection;
    
    template <int dim>
    struct ReciprocalLatticeDirection;
    
    template <int dim>
    class LatticeVector ;
    
    
    struct LatticePlane;
    struct LatticeLine;
//    class PlaneLineIntersection;

} // end namespace

//#include <LatticeGCD.h>
//#include <LatticeBase.h>
#include <Lattice.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>

#include <LatticePlaneBase.h>
#include <LatticePlane.h>
#include <LatticeLine.h>
#include <RationalLatticeDirection.h>

//#include <PlaneLineIntersection.h>
//#include <PlanePlaneIntersection.h>
//#include <LineLineIntersection.h>

#endif
