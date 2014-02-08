/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_H_
#define model_SLIPSYSTEM_H_

#include <cfloat> // FLT_EPSILON
#include <assert.h>
#include <vector>
#include <Eigen/Dense>

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct SlipSystem {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;

        const VectorDimD normal;
        const VectorDimD slip;
        
        SlipSystem(const VectorDimD& normal_in,const VectorDimD& slip_in):
        /* init list */ normal(normal_in),
        /* init list */ slip(slip_in)
        {
            assert(std::fabs(normal.dot(slip))<FLT_EPSILON && "PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL.");
        }
        
    };
    /**************************************************************************/
} // namespace model
#endif
