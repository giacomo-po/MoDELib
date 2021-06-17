/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialBase_H_
#define model_TrialBase_H_

//#include <TestExpression.h>

#include <type_traits> // std::is_same
#include <map>
#include <Eigen/Dense>
#include <TypeTraits.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    /*!\brief A class template that provides the base for all the expressions
     * involving a TrialFunction. The expression-template mechanism is based
     * on the CRTP pattern.
     */
    template<typename _TrialFunctionType>
    class TrialBase
    {
        
    public:
        typedef double Scalar;
        typedef typename TypeTraits<_TrialFunctionType>::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::MeshType MeshType;
        typedef typename TypeTraits<_TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<_TrialFunctionType>::BaryType BaryType;
        
        typedef typename TypeTraits<_TrialFunctionType>::NodeType NodeType;
        typedef _TrialFunctionType TrialFunctionType;
        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
        typedef std::map<size_t,std::array<bool,TypeTraits<TrialFunctionType>::dofPerNode>> DirichletNodeMapType;
        typedef std::map<size_t,double> DirichletConditionContainerType;
        
    private:
        static Eigen::Matrix<double,Eigen::Dynamic,1> dofvector;
        static FiniteElementType* p_fe;
        static DirichletNodeMapType _dirichletNodeMap; // DirichletNodeMap
        static DirichletConditionContainerType _dirichletConditionContainer;
    public:
        
        /**********************************************************************/
        static void init(FiniteElementType& fe_in)
        {
            p_fe=&fe_in;
            dofvector.setZero(gSize());
            std::cout<<"("<<dofvector.size()<<" dof)."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        static FiniteElementType& fe()
        {
            return *p_fe;
        }
        
        /**********************************************************************/
        static size_t gSize()
        {
            return nodeSize()*dofPerNode;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,Eigen::Dynamic,1>& dofVector()
        {/*!\returns A column vector of the DOFs of *this.
          */
            return dofvector;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dofPerNode,1> dofs(const NodeType& node)
        {/*!@param[in] node the FEM node
          * \returns A vector of the degrees of freedom of the node
          */
            return dofvector.template segment<dofPerNode>(node.gID*dofPerNode);
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dofPerNode,1> dofs(const size_t& gID)
        {/*!@param[in] node the FEM node
          * \returns A vector of the degrees of freedom of the node
          */
            return dofvector.template segment<dofPerNode>(gID*dofPerNode);
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele)
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
        static const ElementType& element(const size_t& n)
        {/*!@param[in] n the node ID
          * \returns a const reference to the n-th node in the FiniteElement
          */
            return fe().element(n);
        }
        
        /**********************************************************************/
        static const MeshType& mesh()
        {
            return fe().mesh;
        }
        
        /**********************************************************************/
        static size_t elementSize()
        {/*!\returns the number of elements in the FiniteElement
          */
            return fe().elementSize();
        }
        
        /**********************************************************************/
        static size_t nodeSize()
        {
            return fe().nodeSize();
        }
        
        /**********************************************************************/
        static const NodeType& node(const size_t& n)
        {
            return fe().node(n);
        }
        
        /**********************************************************************/
        static void clearDirichletConditions()
        {
            _dirichletNodeMap.clear();
            _dirichletConditionContainer.clear();
        }
        
        //        /**********************************************************************/
        //        static DirichletNodeMapType& dirichletNodeMap()
        //        {/*!\returns A map of the nodes having a DirichletCondition. The key is
        //          * the global ID of the node (gID), and the value is an array of bool
        //          * indicating which component is constrained.
        //          */
        //            return _dirichletNodeMap;
        //        }
        //
        /**********************************************************************/
        static DirichletNodeMapType& dirichletNodeMap()
        {/*!\returns A map of the nodes having a DirichletCondition. The key is
          * the global ID of the node (gID), and the value is an array of bool
          * indicating which component is constrained.
          */
            return _dirichletNodeMap;
        }
        
        //        /**********************************************************************/
        //        static DirichletConditionContainerType& dirichletConditions()
        //        {
        //            return _dirichletConditionContainer;
        //        }
        
        /**********************************************************************/
        static DirichletConditionContainerType& dirichletConditions()
        {
            return _dirichletConditionContainer;
        }
        
        /**********************************************************************/
        template<typename Condition>
        static void addDirichletCondition(const size_t& nodeListID,
                                          const Condition& cond,
                                          const std::array<bool,dofPerNode>& constrainDof)
        {/*!@param[in] nodeListID ID of the FiniteElement nodeList to which this DirichletCondition applies
          * @param[in] cond the condition object
          * @param[in] constrainDof array of booleans indicating which dofs are to be constrained
          */
            
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Adding Dirichlet condition..."<<std::flush;
            // Check that at least one constrainDof is true
            bool isConstrained(false);
            for(int dof=0;dof<dofPerNode;++dof)
            {
                isConstrained+=constrainDof[dof];
            }
            
            if(isConstrained)
            {
                // Add to dirichletConditions
                for(const auto& pNode :  TrialBase<TrialFunctionType>::fe().nodeList(nodeListID))
                {
                    //const Eigen::Matrix<double,dofPerNode,1> value(cond(*pNode,*this));
                    Eigen::Matrix<double,dofPerNode,1> value(Eigen::Matrix<double,dofPerNode,1>::Zero());
                    cond(*pNode,value);
                    //assert(0 && "RE-ENABLE ABOVE");
                    // DirichletBoundaryCondition<TrialFunctionType> value(*pNode); // TO DO , implement this
                    // cond(value); // TO DO , implement this
                    for(int dof=0;dof<dofPerNode;++dof)
                    {
                        _dirichletNodeMap[pNode->gID][dof]+=constrainDof[dof];
                        
                        if(constrainDof[dof])
                        {
                            const auto temp=_dirichletConditionContainer.emplace(dofPerNode*(pNode->gID)+dof,value(dof));
                            // assert that the condition has been inserted, or that an equivalent condition already existed
                            if(!(temp.second || std::abs(temp.first->second-value(dof))<FLT_EPSILON ))
                            {
                                std::cout<<"FEM node "<<pNode->gID<<", dof "<<dof<<std::endl;
                                std::cout<<"exising dirichlet-condition= "<<temp.first->second<<std::endl;
                                std::cout<<"    new dirichlet-condition= "<<value(dof)<<std::endl;
                                assert(0 && "CONFLICTING DIRICHLET CONDITIONS");
                            }
                        }
                    }
                }
            }
            
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
    };
    
    template<typename TrialFunctionType>
    Eigen::Matrix<double,Eigen::Dynamic,1> TrialBase<TrialFunctionType>::dofvector;
    
    template<typename TrialFunctionType>
    typename TypeTraits<TrialFunctionType>::FiniteElementType* TrialBase<TrialFunctionType>::p_fe;
    
    template<typename TrialFunctionType>
    typename std::map<size_t,std::array<bool,TypeTraits<TrialFunctionType>::dofPerNode>> TrialBase<TrialFunctionType>::_dirichletNodeMap;
    
    template<typename TrialFunctionType>
    std::map<size_t,double> TrialBase<TrialFunctionType>::_dirichletConditionContainer;
    
    
    
    
    
}	// close namespace
#endif


