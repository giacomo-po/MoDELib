/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AtPoint_H_
#define model_AtPoint_H_

#include <Eigen/Dense>

namespace model
{
    class AtPoint
    {
        

        
    public:
        
        const Eigen::VectorXd P0;
        
        /**************************************/
        template <typename FiniteElementType>
        AtPoint(const FiniteElementType&,
                const Eigen::VectorXd& P0_in) :
        P0(P0_in)
        {

        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return (node.P0-P0).norm()<FLT_EPSILON;
            
        }
        
    };
    
    
}	// close namespace
#endif

