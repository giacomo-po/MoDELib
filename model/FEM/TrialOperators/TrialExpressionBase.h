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
        
//        /**********************************************************************/
//        template <typename ElementType, typename BaryType>
//        Eigen::Matrix<double,Eigen::Dynamic,1> operator()(const ElementType& ele, const BaryType& bary) const
//        {/*!@param[in] ele the element
//          * @param[in] bary the vector of barycentric coordinates
//          * \returns the value of the Derived expression at bary.
//          *
//          * \todo: in order to be optimized, this function should be Derived-specific
//          */
//            return derived().sfm(ele,bary)*derived().trial().dofs(ele);
//        }
        
        //        /**********************************************************************/
        //        template <int dim>
        //        Eigen::Matrix<double,Eigen::Dynamic,1> operator()(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* guess) const
        //        {/*!@param[in] ele the element
        //          * @param[in] bary the vector of barycentric coordinates
        //          * \returns the value of the Derived expression at bary.
        //          */
        //            const std::pair<bool,const ElementType*> temp=derived().fe.searchWithGuess(P,guess);
        //
        //
        //            Eigen::Matrix<double,Eigen::Dynamic,1> temp1(derived().sfm(*temp.second,bary)*derived().trial().dofs(*temp.second));
        //
        //            return temp.first?  : ;
        //        }
        
        
    };
    
    
    /**************************************************************************/
    template <typename T, typename Derived>
    T& operator << (T& os, const TrialExpressionBase<Derived>& expr)
    {/*!@param[in] os the stream object
      * @param[in] expr the expression
      *
      * Outputs the value of Derived on the faces of the mes
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
                for (int f=0;f<boundaryFaces.size();++f)
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






/**********************************************************************/
//        void onExternalNodes() const
/**********************************************************************/
//		template <class T>
//		friend T& operator << (T& os, const TrialExpressionBase<Derived>& expr)
//        {
////            const int dim=expr.derived().trial().dim;
//            const int dim=3;
//            const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
////            const Eigen::MatrixXd vertexBary(Eigen::MatrixXd::Identity(dim,dim));
//
//
//            for (int n=0;n<expr.derived().trial().fe.elementSize();++n)
//            {
//                if(expr.derived().trial().fe.element(n).isBoundaryElement())
//                {
//                    const std::vector<int> boundaryFaces=expr.derived().trial().fe.element(n).boundaryFaces();
//                    for (int f=0;f<boundaryFaces.size();++f)
//                    {
//                        for (int v=0;v<dim+1;++v)
//                        {
//                        if (v!=boundaryFaces[f])
//                        {
//                        os<<expr.derived().trial().fe.element(n).position(vertexBary.col(v)).transpose()<<" "
//                          <<expr(expr.derived().trial().fe.element(n),vertexBary.col(v)).transpose()<<"\n";
//                        }
//                        }
//                    }
//                }
//            }
//            return os;
//        }