/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FixMin_H_
#define model_FixMin_H_

#include <Eigen/Dense>

namespace model
{
    template <int i>
    class FixMin
    {
        
        const double val;
        double min;
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        FixMin(const FiniteElementType& fe) :
        val(0.0),
        min(DBL_MAX)
        {
            for (int n=0;n<fe.nodeSize();++n)
            {
                if (fe.node(n).p0(i)<min)
                {
                    min=fe.node(n).p0(i);
                }
            }
        }
        
        /**************************************/
        template <typename NodeType>
        std::pair<bool,double> operator()(const NodeType& node) const
        {
            return std::pair<bool,double>(node.p0(i)==min,val);
            
        }
        
    };
    
    
}	// close namespace
#endif

