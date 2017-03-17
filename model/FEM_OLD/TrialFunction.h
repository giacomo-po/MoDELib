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

#include <model/Utilities/TerminalColors.h>
#include <model/Utilities/NonCopyable.h>

#include <model/FEM/TrialFunctionTraits.h>
#include <model/FEM/Constant.h>
#include <model/FEM/TrialOperators/TrialSum.h>
#include <model/FEM/TrialOperators/TrialProd.h>
#include <model/FEM/TrialOperators/TrialGrad.h>
#include <model/FEM/TrialOperators/TrialDef.h>
//#include <model/FEM/Boundaries/NodeList.h>
#include <model/MPI/MPIcout.h>




namespace model
{
    
    template<int _nComponents, typename _FiniteElementType>
    class TrialFunction : public TrialExpressionBase<TrialFunction<_nComponents,_FiniteElementType> >,
    /*                 */ private std::map<size_t,double>, // DirichletConditionContainer
    /*                 */ private std::map<size_t,std::array<bool,TypeTraits<TrialFunction<_nComponents,_FiniteElementType>>::dofPerNode>> // DirichletNodeMap
    {
     
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        typedef double Scalar;
        typedef _FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename FiniteElementType::NodeType NodeType;
        
        typedef TrialFunction<_nComponents,FiniteElementType> TrialFunctionType;
        
        constexpr static int rows=_nComponents;
        constexpr static int cols=1;

        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
        typedef typename TypeTraits<TrialFunctionType>::NodeContainerType NodeContainerType;
        typedef std::map<size_t,double> DirichletConditionContainerType;
        typedef std::map<size_t,std::array<bool,dofPerNode>>  DirichletNodeMapType;
        
        
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionMatrixType     ShapeFunctionMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionGradMatrixType ShapeFunctionGradMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionDefMatrixType  ShapeFunctionDefMatrixType;
        
        //! A const reference of the FiniteElement over which *this is constructed.
//        const FiniteElementType& fe;
        FiniteElementType& fe; // made non-const only to allow fe.createNodeList

    private:

        Eigen::Matrix<double,Eigen::Dynamic,1> dofvector;

    public:
        
        /**********************************************************************/
        TrialFunction(FiniteElementType& fe_in) : // made non-const only to allow fe.createNodeList
        /* init list */ fe(fe_in)
        {
            model::cout<<greenColor<<"Creating TrialFunction..."<<std::flush;
            dofvector.setZero(this->gSize());
            model::cout<<"("<<this->gSize()<<" dof)."<<defaultColor<<std::endl;
        }
        
        TrialFunction(const TrialFunction& other) = delete; // do not copy trial functions
        TrialFunction(TrialFunction&& other) = default; // do not copy trial functions

        /**********************************************************************/
        TrialFunctionType& operator=(const double& c)
        {
            dofvector.setConstant(c);
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& operator=(const Eigen::VectorXd& temp)
        {
            assert(size_t(temp.size())==this->gSize() && "DOF SIZE MISMATCH");
            dofvector=temp;
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& operator+=(const Eigen::VectorXd& temp)
        {
            assert(size_t(temp.size())==this->gSize() && "DOF SIZE MISMATCH");
            dofvector+=temp;
            return *this;
        }
        
        /**********************************************************************/
        TrialFunctionType& operator-=(const Eigen::VectorXd& temp)
        {
            assert(temp.size()==this->gSize() && "DOF SIZE MISMATCH");
            dofvector-=temp;
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
        void addDirichletCondition(const size_t& nodeListID,
                                   const Condition& cond,
                                   const std::array<bool,dofPerNode>& constrainDof)
        {/*!@param[in] nodeListID ID of the FiniteElement nodeList to which this DirichletCondition applies
          * @param[in] cond the condition object
          * @param[in] constrainDof array of booleans indicating which dofs are to be constrained
          */
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"Adding Dirichlet condition..."<<std::flush;
            // Check that at least one constrainDof is true
            bool isConstrained(false);
            for(int dof=0;dof<dofPerNode;++dof)
            {
                isConstrained+=constrainDof[dof];
            }
            
            if(isConstrained)
            {
                // Add to dirichletConditions
                for(const auto& pNode :  fe.nodeList(nodeListID))
                {
                    const Eigen::Matrix<Scalar,dofPerNode,1> value(cond(*pNode,*this));
                    // DirichletBoundaryCondition<TrialFunctionType> value(*pNode); // TO DO , implement this
                    // cond(value); // TO DO , implement this
                    for(int dof=0;dof<dofPerNode;++dof)
                    {
                        dirichletNodeMap()[pNode->gID][dof]+=constrainDof[dof];
                        
                        if(constrainDof[dof])
                        {
                            const auto temp=dirichletConditions().emplace(dofPerNode*(pNode->gID)+dof,value(dof));
                            // assert that the condition has been inserted, or that an equivalent condition already existed
                            if(!(temp.second || std::abs(temp.first->second-value(dof))<FLT_EPSILON ))
                            {
                                model::cout<<"FEM node "<<pNode->gID<<", dof "<<dof<<std::endl;
                                model::cout<<"exising dirichlet-condition= "<<temp.first->second<<std::endl;
                                model::cout<<"    new dirichlet-condition= "<<value(dof)<<std::endl;
                                assert(0 && "CONFLICTING DIRICHLET CONDITIONS");
                            }
                        }
                    }
                }
            }

            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void clearDirichletConditions()
        {
            dirichletNodeMap().clear();
            dirichletConditions().clear();
        }
        
        /**********************************************************************/
        DirichletConditionContainerType& dirichletConditions()
        {
            return *this;
        }
        
        /**********************************************************************/
        const DirichletConditionContainerType& dirichletConditions() const
        {
            return *this;
        }
        
        /**********************************************************************/
        DirichletNodeMapType& dirichletNodeMap()
        {/*!\returns A map of the nodes having a DirichletCondition. The key is
          * the global ID of the node (gID), and the value is an array of bool
          * indicating which component is constrained.
          */
            return *this;
        }
        
        /**********************************************************************/
        const DirichletNodeMapType& dirichletNodeMap() const
        {/*!\returns A map of the nodes having a DirichletCondition. The key is
          * the global ID of the node (gID), and the value is an array of bool
          * indicating which component is constrained.
          */
            return *this;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,Eigen::Dynamic,1>& dofVector() const
        {/*!\returns A column vector of the DOFs of *this.
          */
            return dofvector;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerNode,1> dofs(const NodeType& node) const
        {/*!@param[in] node the FEM node
          * \returns A vector of the degrees of freedom of the node
          */
            return dofvector.template segment<dofPerNode>(node.gID*dofPerNode);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerNode,1> dofs(const size_t& gID) const
        {/*!@param[in] node the FEM node
          * \returns A vector of the degrees of freedom of the node
          */
            return dofvector.template segment<dofPerNode>(gID*dofPerNode);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele) const
        {/*!@param[in] ele the element
          * \returns A vector of the degrees of freedom of *this on ele
          */
            Eigen::Matrix<double,dofPerElement,1> temp;
            for(int n=0;n<nodesPerElement;++n)
            {
                temp.template segment<dofPerNode>(n*dofPerNode)=dofs(ele.node(n));
            }
            return temp;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const ElementType& ele, const BaryType& bary) const
        {/*!@param[in] ele the element
          * @param[in] bary the vector of barycentric coordinates
          * \returns the value of the Derived expression at bary.
          */
            return sfm(ele,bary)*dofs(ele);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess) const
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
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P the position vector
          * @param[in] guess the Simplex where the search starts
          * \returns the value of the Derived expression at P.
          */
            return operator()(P,&(fe.mesh.begin()->second));
        }
        
        
    };
    
}	// close namespace
#endif
