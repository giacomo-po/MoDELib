/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AtxMin_H_
#define model_AtxMin_H_

#include <Eigen/Dense>

namespace model
{
    template <int i>
    class AtXmin
    {
        
        const double min;
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        AtXmin(const FiniteElementType& fe) :
        min(fe.xMin()(i))
        {

        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return std::abs(node.p0(i)-min)<FLT_EPSILON;
            
        }
        
    };
    
    
}	// close namespace
#endif

