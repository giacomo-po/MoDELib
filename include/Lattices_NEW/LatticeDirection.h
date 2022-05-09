/* This file is part of model.
 *
 * model is distributed without any warranty under the MIT License.
 */


#ifndef model_LatticeDirection_h_
#define model_LatticeDirection_h_

#include <LatticeModule.h>

namespace model
{
    template <int dim>
    struct LatticeDirection : public LatticeVector<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        
        LatticeDirection(const LatticeVector<dim>& v) ;
        
    };
    
} // end namespace
#endif
