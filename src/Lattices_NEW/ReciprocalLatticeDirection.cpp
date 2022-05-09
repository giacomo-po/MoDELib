/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_ReciprocalLatticeDirection_cpp_
#define model_ReciprocalLatticeDirection_cpp_

#include <LatticeModule.h>

namespace model
{
    template <int dim>
    ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v) :
    /* init */ ReciprocalLatticeVector<dim>(((v.squaredNorm()==0)? v : (v/IntegerMath<IntScalarType>::gcd(v)).eval()),v.lattice)
    {
    }


    template struct ReciprocalLatticeDirection<1>;
    template struct ReciprocalLatticeDirection<2>;
    template struct ReciprocalLatticeDirection<3>;

} // end namespace
#endif
