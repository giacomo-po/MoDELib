/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialExpressionBase_H_
#define model_TrialExpressionBase_H_

#include <model/Utilities/TypeTraits.h>
#include <model/FEM/TrialOperators/TestExpression.h>
#include <model/Mesh/Simplex.h>

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
        size_t gSize() const
        {
            return derived().trial().fe.nodeSize()*derived().trial().dofPerNode;
        }

    };
    
    
    /**************************************************************************/
    template <typename T, typename Derived>
    T& operator << (T& os, const TrialExpressionBase<Derived>& expr)
    {/*!@param[in] os the stream object
      * @param[in] expr the expression
      *
      * Outputs the value of Derived on the faces of the mesh
      */
        constexpr int dim=Derived::TrialFunctionType::dim;
        const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
        
        
        for (typename Derived::FiniteElementType::ElementContainerType::const_iterator eIter =expr.derived().trial().fe.elementBegin();
             /*                                                                     */ eIter!=expr.derived().trial().fe.elementEnd();
             /*                                                                   */ ++eIter)
        {
            if(eIter->second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
                for (size_t f=0;f<boundaryFaces.size();++f)
                {
                    for (int v=0;v<dim+1;++v)
                    {
                        if (v!=boundaryFaces[f])
                        {
                            os<<eIter->second.position(vertexBary.col(v)).transpose()<<" "
                            <<expr.derived()(eIter->second,vertexBary.col(v)).transpose()<<"\n";
                        }
                    }
                }
            }
        }
        return os;
    }
    
    
}	// close namespace
#endif
