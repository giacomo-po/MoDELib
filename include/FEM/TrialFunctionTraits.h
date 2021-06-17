/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialFunctionTraits_H_
#define model_TrialFunctionTraits_H_

#include <Eigen/Dense>
#include <TypeTraits.h>


namespace model
{
    
    template<char name,int nComponents, typename _FiniteElementType>
	class TrialFunction;
    
    template<typename TrialExpressionType>
	struct TestExpression;
    
    template<typename Derived>
	struct TrialExpressionBase;
    
    /*!
     * NOTE: if dofPerElement is too large, you may get the Eigen error OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG
     * #define EIGEN_STACK_ALLOCATION_LIMIT 1000000
     */
    template<char _name, int _nComponents, typename _FiniteElementType>
	struct TypeTraits<TrialFunction<_name,_nComponents,_FiniteElementType> >
    {
        typedef _FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename FiniteElementType::NodeType NodeType;
        typedef typename FiniteElementType::NodeContainerType NodeContainerType;
        
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int nComponents=_nComponents;
        constexpr static int dim=FiniteElementType::dim;
        constexpr static int dofPerNode=ElementType::dofPerNode(nComponents);
        constexpr static int dofPerElement=dofPerNode*nodesPerElement;
        
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        typedef Eigen::Matrix<double,nComponents,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,nComponents*dim,dofPerElement> ShapeFunctionGradMatrixType;
        typedef Eigen::Matrix<double,(dim*(dim+1))/2,dofPerElement> ShapeFunctionDefMatrixType;
    };
    
	   
}	// close namespace
#endif
