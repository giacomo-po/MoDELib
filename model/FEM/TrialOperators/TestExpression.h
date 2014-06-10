/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TestExpression_H_
#define model_TestExpression_H_


namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
	template<typename TrialExpressionType>
	struct TestExpression : public TrialExpressionType
    {
        
        typedef TestExpression<TrialExpressionType> TestExpressionType ;

        /**********************************************************************/
        TestExpression(const TrialExpressionType& trialExp) :
        /* init list */ TrialExpressionType(trialExp)
        {
            
        }
        
    };
    
}	// close namespace
#endif

