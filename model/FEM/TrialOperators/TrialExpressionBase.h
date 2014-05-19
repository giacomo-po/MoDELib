/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialExpressionBase_H_
#define model_TrialExpressionBase_H_

#include <model/FEM/TrialOperators/TestExpression.h>

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    /*!\brief A class template that provides the base for all the expressions
     * involving a TrialFunction. The expression-template mechanism is based
     * on the CRTP pattern.
     */
	template<typename Derived>
	struct TrialExpressionBase
    {
        /**********************************************************************/
        Derived& derived()
        {/*!\returns A reference to the Derived object
          */
            return *static_cast<Derived*>(this);
        }
        
        /**********************************************************************/
        Derived* p_derived()
        {/*!\returns A pointer to the Derived object
          */
            return  static_cast<Derived*>(this);
        }
        
        /**********************************************************************/
        const Derived& derived() const
        {/*! A const reference to the Derived object
          */
            return *static_cast<const Derived*>(this);
        }
        
        /**********************************************************************/
        const Derived* p_derived() const
        {/*!\returns A pointer to the const Derived object
          */
            return  static_cast<const Derived*>(this);
        }
        
        /**********************************************************************/
        operator const Derived&() const
        {/*!\returns the const Derived object (cast operator)
          */
            return derived();
        }
        
        /**********************************************************************/
        TestExpression<Derived> test() const
        {/*!\returns A TestExpression of the Derived expression.
          */
            return TestExpression<Derived>(derived());
        }
        
        /**********************************************************************/
        template <typename ElementType, typename BaryType>
        Eigen::Matrix<double,Eigen::Dynamic,1> operator()(const ElementType& ele, const BaryType& bary) const
        {/*!@param[in] ele the element
          * @param[in] bary the vector of barycentric coordinates
          * \returns the value of the Derived expression at bary.
          */
            return derived().sfm(ele,bary)*derived().trial().dofs(ele);
        }
        
        
        /**********************************************************************/
//        void onExternalNodes() const
		/**********************************************************************/
		template <class T>
		friend T& operator << (T& os, const TrialExpressionBase<Derived>& expr)
        {
//            const int dim=expr.derived().trial().dim;
            const int dim=3;
            const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
//            const Eigen::MatrixXd vertexBary(Eigen::MatrixXd::Identity(dim,dim));

            
            for (int n=0;n<expr.derived().trial().fe.elementSize();++n)
            {
                if(expr.derived().trial().fe.element(n).isBoundaryElement())
                {
                    const std::vector<int> boundaryFaces=expr.derived().trial().fe.element(n).boundaryFaces();
                    for (int f=0;f<boundaryFaces.size();++f)
                    {
                        for (int v=0;v<dim+1;++v)
                        {
                        if (v!=boundaryFaces[f])
                        {
                        os<<expr.derived().trial().fe.element(n).position(vertexBary.col(v)).transpose()<<" "
                          <<expr(expr.derived().trial().fe.element(n),vertexBary.col(v)).transpose()<<"\n";
                        }
                        }
                    }
                }
            }
            return os;
        }
        
        
    };
    
    
    
}	// close namespace
#endif