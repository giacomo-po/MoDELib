/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AxisRotation_H_
#define model_AxisRotation_H_

#include <Eigen/Dense>


namespace model
{


    struct AxisRotation
    {        

        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        
        const VectorDim P;
//        const VectorDim a;
//        const double theta;
        const Eigen::AngleAxisd R;
        
        /**************************************/
        AxisRotation(const VectorDim& P_in,
                     const VectorDim& axis,
                     const double& theta) :
        P(P_in),
        R(theta,axis)
        {
        
        }
        
        /**************************************/
        template <typename NodeType,int dofPerNode>
        Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType& node,
                                                       Eigen::Matrix<double,dofPerNode,1>& val) const
        {
            val=P-node.P0+R*(node.P0-P);
            return val;
        }
        
    };
    
    
}	// close namespace
#endif
