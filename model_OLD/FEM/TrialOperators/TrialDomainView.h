/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialDomainView_H_
#define model_TrialDomainView_H_

//#include <TypeTraits.h>
#include <EvalExpression.h>
#include <ExpressionRef.h>


namespace model
{
    
//    template<typename TrialExpressionType>
//    struct TrialExpressionBase;
    
    /**************************************************************************/
	/**************************************************************************/
	template<typename EvalType,int DimMinusDomainDim>
	struct TrialDomainView
    {
        
//        const EvalExpression<EvalType>& evalExp;
        
        ExpressionRef<EvalExpression<EvalType>> evalExp;
        
        TrialDomainView(const EvalExpression<EvalType>& ee) :
        evalExp(ee)
        {
        
        }
        
        TrialDomainView(EvalExpression<EvalType>&& ee) :
        evalExp(std::move(ee))
        {
            
        }
        

    };
    
    

    
    /**************************************************************************/
    template <typename T, typename TrialExpressionType>
    T& operator << (T& os, const TrialDomainView<TrialExpressionType,0>& dv)
    {/*!@param[in] os the stream object
      * @param[in] expr the expression
      *
      * Outputs the value of TrialExpressionType on the mesh
      */
        constexpr int dim=TrialExpressionType::TrialFunctionType::dim;
        const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
        
        
        
        for(const auto& ele : TrialBase<typename TrialExpressionType::TrialFunctionType>::fe().elements())
        {
            for (int v=0;v<dim+1;++v)
            {
                os<<ele.second.position(vertexBary.col(v)).transpose()<<" "
                <<dv.evalExp()(ele.second,vertexBary.col(v)).transpose()<<"\n";
            }
        }
        
        return os;
    }
    
    /**************************************************************************/
    template <typename T, typename TrialExpressionType>
    T& operator << (T& os, const TrialDomainView<TrialExpressionType,1>& dv)
    {/*!@param[in] os the stream object
      * @param[in] expr the expression
      *
      * Outputs the value of Derived on the faces of the mesh
      */
        constexpr int dim=TrialExpressionType::TrialFunctionType::dim;
        const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
        for(const auto& ele : TrialBase<typename TrialExpressionType::TrialFunctionType>::fe().elements())
        {
            if(ele.second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=ele.second.boundaryFaces();
                for (size_t f=0;f<boundaryFaces.size();++f)
                {
                    for (int v=0;v<dim+1;++v)
                    {
                        if (v!=boundaryFaces[f])
                        {
                            os<<std::setprecision(15)<<std::scientific<<ele.second.position(vertexBary.col(v)).transpose()<<" "
                            <<dv.evalExp()(ele.second,vertexBary.col(v)).transpose()<<"\n";
                        }
                    }
                }
            }
        }
        return os;
    }

    
    
}	// close namespace
#endif

