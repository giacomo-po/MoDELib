/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Fix_H_
#define model_Fix_H_

#include <Eigen/Dense>

namespace model
{

    struct Fix
    {        
        /**************************************/
        template <typename NodeType>
        static double at(const NodeType&)
        {
            return 0.0;
        }
        
    };
    
    
}	// close namespace
#endif

