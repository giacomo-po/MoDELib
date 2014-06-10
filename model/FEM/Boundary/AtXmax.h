/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AtXmax_H_
#define model_AtXmax_H_

#include <Eigen/Dense>

namespace model
{
    template <int i>
    class   AtXmax
    {
        
        const double max;
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        AtXmax(const FiniteElementType& fe) :
        max(fe.xMax()(i))
        {
            
        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return std::abs(node.p0(i)-max)<FLT_EPSILON;
            
        }
        
    };
    
    
}	// close namespace
#endif

