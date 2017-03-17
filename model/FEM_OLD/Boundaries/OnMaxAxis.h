/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_OnMaxAxis_H_
#define model_OnMaxAxis_H_

#include <Eigen/Dense>

namespace model
{
    template <int i>
    class OnMaxAxis
    {
        
        const Eigen::MatrixXd axisPoint;
        
        
        Eigen::MatrixXd get_axisPoint(const Eigen::MatrixXd& xMax) const
        {
            Eigen::MatrixXd temp(xMax*0.0);
            temp(i)=xMax(i);
            return temp;
        }
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        OnMaxAxis(const FiniteElementType& fe) :
        axisPoint(get_axisPoint(fe.xMin()))
        {

        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return (node.P0-axisPoint).norm()<FLT_EPSILON;
            
        }
        
    };
    
    
}	// close namespace
#endif

