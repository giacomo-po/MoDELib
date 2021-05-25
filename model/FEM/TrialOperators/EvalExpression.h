/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EvalExpression_H_
#define model_EvalExpression_H_

#include <utility> // for std::move
#include <EvalFunction.h>
#include <ExpressionRef.h>


namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename TrialExpressionType>
    struct EvalExpression : public EvalFunction<EvalExpression<TrialExpressionType>>
    {
        typedef typename TrialExpressionType::TrialFunctionType TrialFunctionType;
        constexpr static int rows=TrialExpressionType::rows;
        constexpr static int cols=TrialExpressionType::cols;

        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;
        
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        
        typedef TrialBase<TrialFunctionType> TrialBaseType;
                
        ExpressionRef<TrialExpressionType> trialExp;
        
        /**********************************************************************/
        EvalExpression(const TrialExpressionType& te) :
        /* init */ trialExp(te)
        {
        }
        
        /**********************************************************************/
        EvalExpression(TrialExpressionType&& te) :
        /* init */ trialExp(std::move(te))
        {
        }
        
        /**********************************************************************/
        EvalMatrixType operator()(const ElementType& ele,
                                  const BaryType& bary) const
        {/*!@param[in] ele the element
          * @param[in] bary the vector of barycentric coordinates
          * \returns the value of trialExp in the element ele at baricentric coordinate bary.
          */
            return trialExp().sfm(ele,bary)*TrialBaseType::dofs(ele);
        }
        
        /**********************************************************************/
        EvalMatrixType operator()(const Eigen::Matrix<double,dim,1>& P,
                                  const Simplex<dim,dim>* guess) const
        {/*!@param[in] P the position vector
          * @param[in] guess the Simplex where the search starts
          * \returns the value of trialExp expression at position P. If P is
          * ouside the mesh, a zero matrix is returned.
          */
            const std::pair<bool,const ElementType*> temp=TrialBaseType::fe().searchWithGuess(P,guess);
            EvalMatrixType val(EvalMatrixType::Zero());
            if(temp.first)
            {
                val=trialExp().sfm(*(temp.second),temp.second->simplex.pos2bary(P))*TrialBaseType::dofs(*(temp.second));
            }
            else
            {
                std::cout<<"WARNING: EVALUATING EXPRESSION OUTSIDE DOMAIN"<<std::endl;
            }
            return val;
        }
        
        /**********************************************************************/
        EvalMatrixType operator()(const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P the position vector
          * \returns the value of trialExp expression at position P. The search 
          * P starts by default in the first simplex of the mesh. If P is
          * ouside the mesh, a zero matrix is returned.
          */
            return this->operator()(P,&(TrialBaseType::mesh().simplices().begin()->second));
        }
        
    };
    
}
#endif
