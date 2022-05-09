/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_LatticeMath_
#define model_LatticeMath_

namespace model
{

template <int dim>
    class Lattice;

    template <int dim>
    class LatticeVectorBase;

    template <int dim>
    class LatticeVector;

    template <int dim>
    class ReciprocalLatticeVectorBase;

    template <int dim>
    class ReciprocalLatticeVector;

    template <int dim>
    struct LatticeDirection;
    
    template <int dim>
    struct ReciprocalLatticeDirection;

    template <int dim>
    struct RationalLatticeDirection;

    template <int dim>
    class BiCrystal ;
        
//    struct LatticePlane;
//    struct LatticeLine;

} // end namespace

//#include <LatticeGCD.h>
//#include <LatticeBase.h>
#include <LatticeCore.h>
#include <Lattice.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>
#include <RationalLatticeDirection.h>
#include <RationalMatrix.h>
#include <SmithDecomposition.h>
#include <BiCrystal.h>
#include <LatticePlaneBase.h>
#include <LatticePlane.h>
#include <LLL.h>
#include <RLLL.h>


#endif
