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
#include <assert.h>

#include <TerminalColors.h>
#include <NonCopyable.h>
#include <InstanceCount.h>

#include <TrialFunctionTraits.h>
#include <Constant.h>
//#include <ExpressionEvaluator.h>
//#include <TrialSum.h>
#include <TrialProd.h>
#include <TrialGrad.h>
#include <TrialDef.h>
//#include <NodeList.h>


// READ THIS
// http://stackoverflow.com/questions/39267388/best-practices-for-dependency-injection-via-constructor

// std::enable_shared_from_this
// http://shaharmike.com/cpp/shared-ptr/

namespace model
{
    
    template<char _name, int _nComponents, typename _FiniteElementType>
    class TrialFunction : public TrialExpressionBase<TrialFunction<_name,_nComponents,_FiniteElementType> >,
    /*                 */ //public ExpressionEvaluator<TrialFunction<_name,_nComponents,_FiniteElementType>,TrialFunction<_name,_nComponents,_FiniteElementType>>,
//    /*                 */ private std::map<size_t,double>, // DirichletConditionContainer
//    /*                 */ private std::map<size_t,std::array<bool,TypeTraits<TrialFunction<_name,_nComponents,_FiniteElementType>>::dofPerNode>>, // DirichletNodeMap
    /*                 */ public InstanceCount<TrialFunction<_name,_nComponents,_FiniteElementType>>
    {
     
        
    public:
        
        
        typedef double Scalar;
        typedef _FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename FiniteElementType::NodeType NodeType;
        
        typedef TrialFunction<_name,_nComponents,FiniteElementType> TrialFunctionType;
        
        constexpr static int rows=_nComponents;
        constexpr static int cols=1;

        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
//        typedef typename TypeTraits<TrialFunctionType>::NodeContainerType NodeContainerType;
        
        typedef std::map<size_t,double> DirichletConditionContainerType;
        typedef std::map<size_t,std::array<bool,dofPerNode>>  DirichletNodeMapType;
        
        
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionMatrixType     ShapeFunctionMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionGradMatrixType ShapeFunctionGradMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionDefMatrixType  ShapeFunctionDefMatrixType;
        
//        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;
        
        //! A const reference of the FiniteElement over which *this is constructed.
//        const FiniteElementType& fe;
//        FiniteElementType& fe; // made non-const only to allow fe.createNodeList

    private:

//        Eigen::Matrix<double,Eigen::Dynamic,1> dofvector;

    public:
        
        /**********************************************************************/
        TrialFunction(FiniteElementType& fe)
//        /* init list */ fe(fe_in)
        {
            std::cout<<greenBoldColor<<"Creating TrialFunction "<< _name <<" "<<std::flush;
            assert(this->counter()==1 && "More of one TrialFunction of the same type exist. Change char template parameter.");
            TrialBase<TrialFunctionType>::init(fe);
//            dofvector.setZero(this->gSize());
//            std::cout<<"("<<this->gSize()<<" dof)."<<defaultColor<<std::endl;
        }
        
        TrialFunction(const TrialFunctionType& other) = delete; // do not copy trial functions
        TrialFunction(TrialFunctionType&& other) = default; // do not copy trial functions

        /**********************************************************************/
        FiniteElementType& fe() const
        {
            return TrialBase<TrialFunctionType>::fe();
        }
        
        /**********************************************************************/
        TrialFunctionType& operator=(const double& c)
        {
            TrialBase<TrialFunctionType>::dofVector().setConstant(c);
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& setNodeDof(const Eigen::Matrix<double,_nComponents,1>& c)
        {
            for(int n=0; n< nodeSize();++n)
            {
                TrialBase<TrialFunctionType>::dofVector().template segment<_nComponents>(n*_nComponents)=c;
            }
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& operator=(const Eigen::VectorXd& temp)
        {
//            assert(size_t(temp.size())==this->gSize() && "DOF SIZE MISMATCH");
//            dofvector=temp;
            assert(size_t(temp.size())==TrialBase<TrialFunctionType>::gSize() && "DOF SIZE MISMATCH");
            TrialBase<TrialFunctionType>::dofVector()=temp;
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& operator+=(const Eigen::VectorXd& temp)
        {
//            assert(size_t(temp.size())==this->gSize() && "DOF SIZE MISMATCH");
//            dofvector+=temp;
            assert(size_t(temp.size())==TrialBase<TrialFunctionType>::gSize() && "DOF SIZE MISMATCH");
            TrialBase<TrialFunctionType>::dofVector()+=temp;
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& operator-=(const Eigen::VectorXd& temp)
        {
//            assert(temp.size()==this->gSize() && "DOF SIZE MISMATCH");
//            dofvector-=temp;
            assert(size_t(temp.size())==TrialBase<TrialFunctionType>::gSize() && "DOF SIZE MISMATCH");
            TrialBase<TrialFunctionType>::dofVector()-=temp;
            return *this;
        }
        
//        /**********************************************************************/
//        const TrialFunctionType& trial() const
//        {
//            return *this;
//        }

        /**********************************************************************/
        size_t gSize() const
        {/*!\returns the total number of dof of this TrialFunction
          */
            return TrialBase<TrialFunctionType>::gSize();
        }
        
        /**********************************************************************/
        size_t elementSize() const
        {/*!\returns the number of elements in the FiniteElement
          */
            return TrialBase<TrialFunctionType>::elementSize();
        }
        
        /**********************************************************************/
        size_t nodeSize() const
        {
            return TrialBase<TrialFunctionType>::nodeSize();
        }
        
        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {
            return TrialBase<TrialFunctionType>::node(n);
        }
        
        /**********************************************************************/
        static ShapeFunctionMatrixType sfm(const ElementType& ele, const BaryType& bary)
        {/*!\param[in] ele a finite element
          * \param[in] bary the vector of barycentric coordinates
          * \returns The matrix of shape functions for the element, evaluated at
          * bary.
          */
            return ele.template sfm<TrialFunctionType>(bary);
        }
        
        /**********************************************************************/
        static ShapeFunctionGradMatrixType sfmGrad(const ElementType& ele, const BaryType& bary)
        {/*!\param[in] ele a finite element
          * \param[in] bary the vector of barycentric coordinates
          * \returns The matrix of shape functions derivatives for the element,
          * evaluated at bary.
          */
            return ele.template sfmGrad<TrialFunctionType>(bary);
        }
        
        /**********************************************************************/
        static ShapeFunctionDefMatrixType sfmDef(const ElementType& ele, const BaryType& bary)
        {/*!\param[in] ele a finite element
          * \param[in] bary the vector of barycentric coordinates
          * \returns The matrix of symmetric shape functions derivatives for the element,
          * evaluated at bary.
          */
            return ele.template sfmDef<TrialFunctionType>(bary);
        }

        /**********************************************************************/
        template<typename Condition>
        void addDirichletCondition(const size_t& nodeListID,
                                   const Condition& cond,
                                   const std::array<bool,dofPerNode>& constrainDof)
        {/*!@param[in] nodeListID ID of the FiniteElement nodeList to which this DirichletCondition applies
          * @param[in] cond the condition object
          * @param[in] constrainDof array of booleans indicating which dofs are to be constrained
          */
            
            TrialBase<TrialFunctionType>::addDirichletCondition(nodeListID,cond,constrainDof);
        }
        
//        /**********************************************************************/
//        template<typename Condition>
//        void addDirichletCondition(const size_t& nodeListID,
//                                   const Condition& cond,
//                                   const std::array<bool,dofPerNode>& constrainDof)
//        {/*!@param[in] nodeListID ID of the FiniteElement nodeList to which this DirichletCondition applies
//          * @param[in] cond the condition object
//          * @param[in] constrainDof array of booleans indicating which dofs are to be constrained
//          */
//            
//            const auto t0= std::chrono::system_clock::now();
//            std::cout<<"Adding Dirichlet condition..."<<std::flush;
//            // Check that at least one constrainDof is true
//            bool isConstrained(false);
//            for(int dof=0;dof<dofPerNode;++dof)
//            {
//                isConstrained+=constrainDof[dof];
//            }
//            
//            if(isConstrained)
//            {
//                // Add to dirichletConditions
//                for(const auto& pNode :  TrialBase<TrialFunctionType>::fe().nodeList(nodeListID))
//                {
//                    const Eigen::Matrix<Scalar,dofPerNode,1> value(cond(*pNode,*this));
//                    // DirichletBoundaryCondition<TrialFunctionType> value(*pNode); // TO DO , implement this
//                    // cond(value); // TO DO , implement this
//                    for(int dof=0;dof<dofPerNode;++dof)
//                    {
//                        dirichletNodeMap()[pNode->gID][dof]+=constrainDof[dof];
//                        
//                        if(constrainDof[dof])
//                        {
//                            const auto temp=dirichletConditions().emplace(dofPerNode*(pNode->gID)+dof,value(dof));
//                            // assert that the condition has been inserted, or that an equivalent condition already existed
//                            if(!(temp.second || std::abs(temp.first->second-value(dof))<FLT_EPSILON ))
//                            {
//                                std::cout<<"FEM node "<<pNode->gID<<", dof "<<dof<<std::endl;
//                                std::cout<<"exising dirichlet-condition= "<<temp.first->second<<std::endl;
//                                std::cout<<"    new dirichlet-condition= "<<value(dof)<<std::endl;
//                                assert(0 && "CONFLICTING DIRICHLET CONDITIONS");
//                            }
//                        }
//                    }
//                }
//            }
//
//            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//        }
        
        /**********************************************************************/
        void clearDirichletConditions()
        {
            TrialBase<TrialFunctionType>::clearDirichletConditions();
//            TrialBase<TrialFunctionType>::dirichletConditions().clear();
        }
        
//        /**********************************************************************/
//        DirichletConditionContainerType& dirichletConditions()
//        {
//            return *this;
//        }
//        
        /**********************************************************************/
        DirichletConditionContainerType& dirichletConditions() const
        {
            return TrialBase<TrialFunctionType>::dirichletConditions();
        }
        
//        /**********************************************************************/
//        DirichletNodeMapType& dirichletNodeMap()
//        {/*!\returns A map of the nodes having a DirichletCondition. The key is
//          * the global ID of the node (gID), and the value is an array of bool
//          * indicating which component is constrained.
//          */
//            return *this;
//        }
//        
        /**********************************************************************/
        DirichletNodeMapType& dirichletNodeMap() const
        {/*!\returns A map of the nodes having a DirichletCondition. The key is
          * the global ID of the node (gID), and the value is an array of bool
          * indicating which component is constrained.
          */
            return TrialBase<TrialFunctionType>::dirichletNodeMap();
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,Eigen::Dynamic,1>& dofVector() const
        {/*!\returns A column vector of the DOFs of *this.
          */
            return TrialBase<TrialFunctionType>::dofVector();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerNode,1> dofs(const NodeType& node) const
        {/*!@param[in] node the FEM node
          * \returns A vector of the degrees of freedom of the node
          */
//            return dofvector.template segment<dofPerNode>(node.gID*dofPerNode);
            return TrialBase<TrialFunctionType>::dofs(node);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerNode,1> dofs(const size_t& gID) const
        {/*!@param[in] node the FEM node
          * \returns A vector of the degrees of freedom of the node
          */
//            return dofvector.template segment<dofPerNode>(gID*dofPerNode);
            return TrialBase<TrialFunctionType>::dofs(gID);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele) const
        {/*!@param[in] ele the element
          * \returns A vector of the degrees of freedom of *this on ele
          */
//            Eigen::Matrix<double,dofPerElement,1> temp;
//            for(int n=0;n<nodesPerElement;++n)
//            {
//                temp.template segment<dofPerNode>(n*dofPerNode)=dofs(ele.node(n));
//            }
//            return temp;
            return TrialBase<TrialFunctionType>::dofs(ele);

        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const ElementType& ele, const BaryType& bary) const
//        {/*!@param[in] ele the element
//          * @param[in] bary the vector of barycentric coordinates
//          * \returns the value of the Derived expression at bary.
//          */
//            return sfm(ele,bary)*dofs(ele);
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            const std::pair<bool,const ElementType*> temp=TrialBase<TrialFunctionType>::fe().searchWithGuess(P,guess);
//            //            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
//            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
//            if(temp.first)
//            {
//                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*dofs(*(temp.second));
//            }
//            return val;
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
////            return operator()(P,&(fe.mesh.begin()->second));
//            return operator()(P,&(TrialBase<TrialFunctionType>::mesh().begin()->second));
//
//        }
        
        
    };
    
}	// close namespace
#endif
