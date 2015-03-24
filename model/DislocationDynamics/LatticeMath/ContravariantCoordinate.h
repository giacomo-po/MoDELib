/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ContravariantCoordinate_h_
#define model_ContravariantCoordinate_h_

#include <Eigen/Dense>


namespace model
{
    template <typename T,int dim>
    struct ContravariantCoordinate : public Eigen::Matrix<T,dim,1>
    {
        
    };
    
} // end namespace
#endif
