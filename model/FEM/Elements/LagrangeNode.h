/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LagrangeNode_H_
#define model_LagrangeNode_H_

#include <Eigen/Dense>

namespace model
{
    
	
    /**************************************************************************/
	/**************************************************************************/
	template<int _dim>
	struct LagrangeNode
    {
        static constexpr int dim=_dim;
        typedef Eigen::Matrix<double,dim,1> PositionType;
        
        const PositionType p0;
        const size_t gID; // global ID
        
        LagrangeNode(const PositionType& p, const size_t& gid) :
        /* init list */ p0(p),
        /* init list */ gID(gid)
        {}
        
    };
    
    
}	// close namespace
#endif