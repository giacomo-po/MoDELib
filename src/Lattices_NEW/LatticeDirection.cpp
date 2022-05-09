/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_LatticeDirection_cpp_
#define model_LatticeDirection_cpp_

#include <LatticeModule.h>
namespace model
{

    
    template <int dim>
    LatticeDirection<dim>::LatticeDirection(const LatticeVector<dim>& v) :
    /* base init */ LatticeVector<dim>(((v.squaredNorm()==0)? v : (v/IntegerMath<IntScalarType>::gcd(v)).eval()),v.lattice)
    {
    }
        

    
    template class LatticeDirection<1>;
    template class LatticeDirection<2>;
    template class LatticeDirection<3>;

} // end namespace
#endif
