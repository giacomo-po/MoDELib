/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DirichletCondition_H_
#define model_DirichletCondition_H_

#include <Eigen/Dense>

namespace model
{

    template <typename _TrialFunctionType>
    struct DirichletCondition
    {
        
        typename _TrialFunctionType TrialFunctionType;

        
        DirichletCondition(const TrialFunctionType& trial_in) :trial(trial_in)
        /* init list */ value(vale_in)
        {
        
        }
        
        
        /**************************************/
        template <typename NodeType>
        double at(const NodeType& node) const
        {
            return value;
        }
        
    };
    
    
}	// close namespace
#endif

