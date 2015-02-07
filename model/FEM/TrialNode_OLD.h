/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_TrialNode_H_
#define  model_TrialNode_H_


namespace model
{
	
    
    /**************************************************************************/
	/**************************************************************************/
	template<typename _TrialFunctionType>
	struct TrialNode
    {
        typedef _TrialFunctionType TrialFunctionType;
        typedef typename TrialFunctionType::NodeType NodeType;
        constexpr static int dofPerNode=TrialFunctionType::dofPerNode;

        const NodeType& node;
        
        TrialNode(const NodeType& node_in) : node(node_in)
        {
        
        }
        
	};
    
	
}
#endif
