/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RuntimeError_H_
#define model_RuntimeError_H_

#include <iostream>
#include <string>

namespace model
{
    struct RuntimeError
    {
        
        RuntimeError(const bool& dontCrash, const std::string& msg)
        {
            if(!dontCrash) // similar to assert
            {
                std::cerr<<msg<<std::endl;
            }
        }
        
    };
    
    
} // namespace model
#endif
