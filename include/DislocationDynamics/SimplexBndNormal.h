/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexBndNormal_H_
#define model_SimplexBndNormal_H_

#include <Eigen/Dense>
#include <Simplex.h>

namespace model
{
    
    class SimplexBndNormal
    {
        
    public:
        
        // make dim available outside class
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        /**********************************************************************/
        static VectorDim get_boundaryNormal(const VectorDim& P,
                                     const Simplex<dim,dim>& simplex,
                                     const double& dmax);
        
        
    };
    
    
    
    
} // close namespace
#endif

