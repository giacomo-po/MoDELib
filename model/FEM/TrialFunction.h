/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialFunction_H_
#define model_TrialFunction_H_

#include <utility>      // std::pair, std::make_pair
#include <deque>

#include <model/Utilities/TerminalColors.h>
#include <model/FEM/TrialFunctionTraits.h>
#include <model/FEM/Constant.h>
#include <model/FEM/TrialOperators/TrialSum.h>
#include <model/FEM/TrialOperators/TrialProd.h>
#include <model/FEM/TrialOperators/TrialGrad.h>
#include <model/FEM/TrialOperators/TrialDef.h>
#include <model/FEM/DirichletCondition.h>



namespace model
{
    
	template<int _nComponents, typename _FiniteElementType>
	struct TrialFunction : public TrialExpressionBase<TrialFunction<_nComponents,_FiniteElementType> >
    {
        
        typedef _FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename FiniteElementType::NodeType NodeType;

        typedef TrialFunction<_nComponents,FiniteElementType> TrialFunctionType;
        constexpr static int rows=_nComponents;

        typedef typename TypeTraits<TrialFunctionType>::NodeContainerType NodeContainerType;
        typedef std::pair<size_t,double> DirichletConditionType;
        typedef std::deque<DirichletConditionType> DirichletConditionContainerType;

        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionMatrixType ShapeFunctionMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionGradMatrixType ShapeFunctionGradMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionDefMatrixType ShapeFunctionDefMatrixType;
        
        const FiniteElementType& fe;
        Eigen::Matrix<double,Eigen::Dynamic,1> dofContainer;
        
        DirichletConditionContainerType dcContainer;
        
        
        /**********************************************************************/
        TrialFunction(const FiniteElementType& fe_in) :
        /* init list */ fe(fe_in)
        {
            std::cout<<greenColor<<"Creating TrialFunction "<<defaultColor<<std::endl;
            
            dofContainer.setZero(fe.nodeSize()*dofPerNode);
            
            std::cout<<"    #dof "<<dofContainer.size()<<std::endl;
            
        }
        
        /**********************************************************************/
        const TrialFunctionType& trial() const
        {
            return *this;
        }
        
        /**********************************************************************/
        size_t elementSize() const
        {/*!\returns the number of elements in the FiniteElement
          */
            return fe.elementSize();
        }
        
        /**********************************************************************/
        const ElementType& element(const size_t& n) const
        {/*!@param[in] n the node ID
          * \returns a const reference to the n-th node in the FiniteElement
          */
            return fe.element(n);
        }
        
//        /**********************************************************************/
//        ElementType& element(const size_t& n)
//        {/*!@param[in] n the node ID
//          * \returns a reference to the n-th node in the FiniteElement
//          */
//            return fe.element(n);
//        }
        
        //        /**********************************************************************/
        //        typename NodeContainerType::const_iterator nodeBegin() const
        //        {
        //            return fe.nodeBegin();
        //        }
        //
        //        /**********************************************************************/
        //        typename NodeContainerType::const_iterator nodeEnd() const
        //        {
        //            return fe.nodeEnd();
        //        }
        
        /**********************************************************************/
        size_t nodeSize() const
        {
            return fe.nodeSize();
        }
        
        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {
            return fe.node(n);
        }

        /**********************************************************************/
        ShapeFunctionMatrixType sfm(const ElementType& ele, const BaryType& bary) const
        {/*!\param[in] ele a finite element
          * \param[in] bary the vector of barycentric coordinates
          * \returns The matrix of shape functions for the element, evaluated at
          * bary.
          */
            return ele.template sfm<TrialFunctionType>(bary);
        }
        
        /**********************************************************************/
        ShapeFunctionGradMatrixType sfmGrad(const ElementType& ele, const BaryType& bary) const
        {/*!\param[in] ele a finite element
          * \param[in] bary the vector of barycentric coordinates
          * \returns The matrix of shape functions derivatives for the element,
          * evaluated at bary.
          */
            return ele.template sfmGrad<TrialFunctionType>(bary);
        }
        
        /**********************************************************************/
        ShapeFunctionDefMatrixType sfmDef(const ElementType& ele, const BaryType& bary) const
        {/*!\param[in] ele a finite element
          * \param[in] bary the vector of barycentric coordinates
          * \returns The matrix of symmetric shape functions derivatives for the element,
          * evaluated at bary.
          */
            return ele.template sfmDef<TrialFunctionType>(bary);
        }
        
        /**********************************************************************/
        template<typename Condition>
        void addDirechletCondition(const Condition& dc, const int& dof)
        {/*!@param[in] dc the Dirichlet condition
          * @param[in] the nodal dof to be constrained
          */
            
            assert(dof>=0 && "dof MUST BE >=0");
            assert(dof<dofPerNode && "dof MUST BE < dofPerNode");
            
            for(typename NodeContainerType::const_iterator nIter=fe.nodeBegin();nIter!=fe.nodeEnd();++nIter)
            {
                std::pair<bool,double> temp(dc(*nIter));
                if(temp.first) // condition applies o current node
                {
                    dcContainer.emplace_back(dofPerNode*nIter->gID+dof,temp.second);
                }
            }
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele) const
        {/*!\returns
          */
            Eigen::Matrix<double,dofPerElement,1> temp;
            for(int n=0;n<nodesPerElement;++n)
            {
                temp.template segment<dofPerNode>(n*dofPerNode)=dofContainer.template segment<dofPerNode>(ele.node(n).gID*dofPerNode);
            }
            return temp;
        }
        

        
    };
    
}	// close namespace
#endif
