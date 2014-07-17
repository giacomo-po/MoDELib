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
#include <map>
#include <list>

#include <model/Utilities/TerminalColors.h>
#include <model/FEM/TrialFunctionTraits.h>
#include <model/FEM/Constant.h>
#include <model/FEM/TrialOperators/TrialSum.h>
#include <model/FEM/TrialOperators/TrialProd.h>
#include <model/FEM/TrialOperators/TrialGrad.h>
#include <model/FEM/TrialOperators/TrialDef.h>
//#include <model/FEM/DirichletCondition.h>
#include <model/FEM/Boundaries/NodeList.h>
#include <model/MPI/MPIcout.h>



namespace model
{
    
	template<int _nComponents, typename _FiniteElementType>
	class TrialFunction : public TrialExpressionBase<TrialFunction<_nComponents,_FiniteElementType> >,
    /*                 */ private Eigen::Matrix<double,Eigen::Dynamic,1>
    {
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        typedef _FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename FiniteElementType::NodeType NodeType;
        
        typedef TrialFunction<_nComponents,FiniteElementType> TrialFunctionType;
        constexpr static int rows=_nComponents;
        
        typedef typename TypeTraits<TrialFunctionType>::NodeContainerType NodeContainerType;
        typedef std::map<size_t,double> DirichletConditionContainerType;
        
        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionMatrixType ShapeFunctionMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionGradMatrixType ShapeFunctionGradMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionDefMatrixType ShapeFunctionDefMatrixType;
        
        const FiniteElementType& fe;
        
//        private:
//        Eigen::Matrix<double,Eigen::Dynamic,1> dofContainer; // \todo make this private inheritance
        
        
//        public:
        DirichletConditionContainerType dcContainer;

        
        /**********************************************************************/
        TrialFunction(const FiniteElementType& fe_in) :
        /* init list */ fe(fe_in)
        {
//            dofContainer.setZero(fe.nodeSize()*dofPerNode);
//             model::cout<<greenColor<<"Creating TrialFunction: "<<dofContainer.size()<<" #dof."<<defaultColor<<std::endl;
            this->setZero(fe.nodeSize()*dofPerNode);
             model::cout<<greenColor<<"Creating TrialFunction: "<<this->size()<<" #dof."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        const TrialFunction& operator=(const Eigen::VectorXd& temp)
        {
            assert(temp.size()==this->size() && "DOF SIZE MISMATCH");
//            dofContainer=temp;
            static_cast<Eigen::Matrix<double,Eigen::Dynamic,1>*>(this)->operator=(temp);
//            dofContainer=temp;
            return *this;
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
        void addDirichletCondition(const Condition& cond, const NodeList<FiniteElementType>& nodeList, const int& dof)
        {/*!@param[in] dc the Dirichlet condition
          * @param[in] the nodal dof to be constrained
          */
            
            assert(dof>=0 && "dof MUST BE >=0");
            assert(dof<dofPerNode && "dof MUST BE < dofPerNode");
            
            //Make sure that (&this->fe)==&(nodeList.fe)
            assert((&this->fe)==(&nodeList.fe) && "USING NODE LIST CREATED FROM A DIFFERENT FINITE ELEMENT.");
            
            for(typename NodeList<FiniteElementType>::const_iterator nIter=nodeList.begin();nIter!=nodeList.end();++nIter)
            {
                const auto temp=dcContainer.emplace(dofPerNode*(*nIter)->gID+dof,cond.at(**nIter));
                assert(temp.second && "UNABLE TO INSERT DIRICHLET CONDITION");
            }
        }
        
        /**********************************************************************/
        void addDirichletCondition(const double& val, const NodeList<FiniteElementType>& nodeList, const int& dof)
        {/*!@param[in] dc the Dirichlet condition
          * @param[in] the nodal dof to be constrained
          */
            
            assert(dof>=0 && "dof MUST BE >=0");
            assert(dof<dofPerNode && "dof MUST BE < dofPerNode");
            
            //Make sure that (&this->fe)==&(nodeList.fe)
            assert((&this->fe)==(&nodeList.fe) && "USING NODE LIST CREATED FROM A DIFFERENT FINITE ELEMENT.");
            
            for(typename NodeList<FiniteElementType>::const_iterator nIter=nodeList.begin();nIter!=nodeList.end();++nIter)
            {
                const auto temp=dcContainer.emplace(dofPerNode*(*nIter)->gID+dof,val);
                assert(temp.second && "UNABLE TO INSERT DIRICHLET CONDITION");
            }
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele) const
        {/*!@param[in] ele the element
          * \returns A vector of the degrees of freedom of *this on ele
          */
            Eigen::Matrix<double,dofPerElement,1> temp;
            for(int n=0;n<nodesPerElement;++n)
            {
                temp.template segment<dofPerNode>(n*dofPerNode)=this->template segment<dofPerNode>(ele.node(n).gID*dofPerNode);
            }
            return temp;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const ElementType& ele, const BaryType& bary) const
        {/*!@param[in] ele the element
          * @param[in] bary the vector of barycentric coordinates
          * \returns the value of the Derived expression at bary.
          *
          * \todo: in order to be optimized, this function should be Derived-specific
          */
            return sfm(ele,bary)*dofs(ele);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* guess) const
        {/*!@param[in] P the position vector
          * @param[in] guess the Simplex where the search starts
          * \returns the value of the Derived expression at P.
          */
            const std::pair<bool,const ElementType*> temp=fe.searchWithGuess(P,guess);
            //            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
            if(temp.first)
            {
                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*dofs(*(temp.second));
            }
            return val;
        }

    };
    
}	// close namespace
#endif


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
