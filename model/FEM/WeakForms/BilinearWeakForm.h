/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakForm_H_
#define model_BilinearWeakForm_H_

#include <chrono>
#include <vector>

//#ifdef _OPENMP
//#include <omp.h>
//#include <Model/Threads/EqualIteratorRange.h>
//#endif

#include <Eigen/Sparse>

#include <model/FEM/WeakForms/LinearWeakForm.h>
#include <model/FEM/WeakForms/WeakProblem.h>
#include <model/Utilities/TerminalColors.h>
#include <model/FEM/WeakForms/BilinearForm.h>
#include <model/MPI/MPIcout.h>



namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename _BilinearFormType,typename _IntegrationDomainType>
	struct BilinearWeakForm
    {
        
        typedef _BilinearFormType BilinearFormType;
        typedef _IntegrationDomainType IntegrationDomainType;
        typedef BilinearWeakForm<BilinearFormType,IntegrationDomainType> BilinearWeakFormType;
        typedef typename BilinearFormType::TrialFunctionType TrialFunctionType;
        typedef typename BilinearFormType::TestExpressionType TestExpressionType;
        typedef typename BilinearFormType::TrialExpressionType TrialExpressionType;
        typedef typename TrialFunctionType::ElementType ElementType;
        typedef typename IntegrationDomainType::QuadratureType QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;
        
        
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int nodesPerElement=TypeTraits<TrialFunctionType>::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
        
        typedef Eigen::Matrix<double,dofPerElement,dofPerElement> ElementMatrixType;
        
        const BilinearFormType bilinearForm;
        const IntegrationDomainType& domain;
        const TestExpressionType testExpr;
        const TrialExpressionType trialExpr;
        const size_t gSize;
        
        
        
        /**********************************************************************/
        BilinearWeakForm(const BilinearFormType& bf, const IntegrationDomainType& dom) :
        /* init list */ bilinearForm(bf), // cast testE to its base T2 type
        /* init list */ domain(dom), // cast trialE to its derived T1 type
        /* init list */ testExpr(bf.testExpr),
        /* init list */ trialExpr(bf.trialExpr),
        /* init list */ gSize(bilinearForm.gSize)
        {
             model::cout<<greenColor<<"Creating BilinearWeakForm: gSize="<<gSize<<defaultColor<<std::endl;
        }
        
//        const TrialFunctionType& trialExpr() const
//        {
//            return trialExprr();
//        }
//
//        const TrialFunctionType& testExpr() const
//        {
//            return trialExprr();
//        }
        
        /**********************************************************************/
        //template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        std::vector<Eigen::Triplet<double> >  assembleOnDomain(double& maxAbsValue) const
        {
            
             model::cout<<"Assembling BilinearWeakForm on domain..."<<std::flush;
            
            std::vector<Eigen::Triplet<double> > globalTriplets;
            //            globalTriplets.clear();
            globalTriplets.reserve(dofPerElement*dofPerElement*trialExpr.elementSize());
            maxAbsValue=0.0;
            
            
            const auto t0= std::chrono::system_clock::now();
            for (int k=0;k<domain.size();++k)
            {
                const ElementType& ele(domain.element(k));
                ElementMatrixType ke(ElementMatrixType::Zero());
                QuadratureType::integrate(this,ke,&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,ele);
                //Assemble into global array of triplets
                for (int i=0;i<dofPerElement;++i)
                {
                    for (int j=0;j<dofPerElement;++j)
                    {
                        if(ke(i,j)!=0.0)
                        {
                            const size_t  nodeID_I(i/dofPerNode);
                            const size_t nodeDof_I(i%dofPerNode);
                            const size_t gI= ele.node(nodeID_I).gID*dofPerNode+nodeDof_I;
                            
                            const size_t  nodeID_J(j/dofPerNode);
                            const size_t nodeDof_J(j%dofPerNode);
                            const size_t gJ=ele.node(nodeID_J).gID*dofPerNode+nodeDof_J;
                            
                            globalTriplets.emplace_back(gI,gJ,ke(i,j));
                            
                            //                            A.coeffRef(gI,gJ) += ke(i,j);
                            
                            if (std::fabs(ke(i,j))>maxAbsValue)
                            {
                                maxAbsValue=std::fabs(ke(i,j));
                            }
                        }
                    }
                }
            }
            
             model::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            return globalTriplets;
        }
        
        /**********************************************************************/
		template <typename AbscissaType>
        ElementMatrixType elementAssemblyKernel(const AbscissaType& a, const ElementType& ele) const
        {
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
            return testExpr.sfm(ele,bary).transpose()*trialExpr.sfm(ele,bary)*ele.absJ(bary);
		}
        
        /**********************************************************************/
        template< typename T3>
        WeakProblem<BilinearWeakFormType,T3> operator=(const LinearWeakExpression<T3>& lwf) const
        {
            return WeakProblem<BilinearWeakFormType,T3>(*this,lwf);
        }
        
    };
    
    
    /**************************************************************************/
    /**************************************************************************/
    // Operator *
    template <typename T1, typename T2, typename FiniteElementType, int qOrder, int dimMinusDomainDim,
    /*     */ template <short unsigned int, short unsigned int> class QuadratureRule>
    BilinearWeakForm<BilinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> > operator*(const BilinearForm<T1,T2>& bilinearForm,
                                                                                                                                  const IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule>& domain)
    {
        return BilinearWeakForm<BilinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> >(bilinearForm,domain);
    }
    
}	// close namespace
#endif









//        /**********************************************************************/
//        WeakProblem<BilinearWeakFormType> operator=(const int& a) const
//        {
//            return WeakProblem<BilinearWeakFormType>(*this);
//        }

//        /**********************************************************************/
//        template< typename T3,typename T4>
//        WeakProblem<BilinearWeakFormType,LinearWeakForm<T3,T4> > operator=(const LinearWeakForm<T3,T4>& lwf) const
//        {
//            return WeakProblem<BilinearWeakFormType,LinearWeakForm<T3,T4> >(*this,lwf);
//        }


//        /**********************************************************************/
//        BilinearWeakForm(BilinearWeakForm&&) = default; // explicit move constructor (required since copy constructor is private)
//
//        /**********************************************************************/
//        BilinearWeakForm& operator=(BilinearWeakForm&&) = default; // explicit a move assignment (required since assignment is private)


//        /**********************************************************************/
//        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//        const BilinearWeakFormType& operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
//        {
//            assembleOnDomain<qOrder,QuadratureRule>();
//            return *this;
//        }


//        /**********************************************************************/
//        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//        BilinearWeakFormType&& operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
//        {/*
//          * See section "Returning an explicit rvalue-reference from a function" in
//          * http://www.cprogramming.com/c++11/rvalue-references-and-move-semantics-in-c++11.html
//          */
//            assembleOnDomain<qOrder,QuadratureRule>();
//            return std::move(*this);
//        }

//        /**********************************************************************/
//        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//        BilinearWeakFormType operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
//        {
//            assembleOnDomain<qOrder,QuadratureRule>();
//            return *this;
//        }






//            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =trialExp.elementBegin();
//                 /*                                                            */ eIter!=trialExp.elementEnd();
//                 /*                                                            */ ++eIter)
//            {
//                // Compute element stiffness Matrix
//                ElementMatrixType ke(ElementMatrixType::Zero());
//                QuadratureType::integrate(this,ke,&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,eIter->second);
//
//                //Assemble into global array of triplets
//                //                for (int i=0;i<dofPerElement;++i)
//                //                {
//                //                    for (int j=0;j<dofPerElement;++j)
//                //                    {
//                //                        if(ke(i,j)!=0.0)
//                //                        {
//                //                            const size_t  nodeID_I(i/dofPerNode);
//                //                            const size_t nodeDof_I(i%dofPerNode);
//                //                            const size_t gI= eIter->second.node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                //
//                //                            const size_t  nodeID_J(j/dofPerNode);
//                //                            const size_t nodeDof_J(j%dofPerNode);
//                //                            const size_t gJ=eIter->second.node(nodeID_J).gID*dofPerNode+nodeDof_J;
//                //
//                //                            globalTriplets.emplace_back(gI,gJ,ke(i,j));
//                //
//                //                            //                            A.coeffRef(gI,gJ) += ke(i,j);
//                //
//                //                            if (std::fabs(ke(i,j))>maxAbsValue)
//                //                            {
//                //                                maxAbsValue=std::fabs(ke(i,j));
//                //                            }
//                //                        }
//                //                    }
//                //                }
//
//                //                    std::map<int,const int> iMap;
//                //                    for (int i=0;i<dofPerElement;++i)
//                //                    {
//                //                        const size_t  nodeID_I(i/dofPerNode);
//                //                        const size_t nodeDof_I(i%dofPerNode);
//                //                        const size_t gI= eIter->second.node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                //                        iMap.emplace(gI,i);
//                //                    }
//                //
//                //                    std::map<int,const int> jMap;
//                //                    for (int j=0;j<dofPerElement;++j)
//                //                    {
//                //                        const size_t  nodeID_J(j/dofPerNode);
//                //                        const size_t nodeDof_J(j%dofPerNode);
//                //                        const size_t gJ=eIter->second.node(nodeID_J).gID*dofPerNode+nodeDof_J;
//                //                        jMap.emplace(gJ,j);
//                //                    }
//                //
//                //                    for (std::map<int,const int>::const_iterator jIter=jMap.begin();jIter!=jMap.end();++jIter)
//                //                        for (std::map<int,const int>::const_iterator iIter=iMap.begin();iIter!=iMap.end();++iIter)
//                //                    {
//                //                        {
//                //                            A.insert(iIter->first,jIter->first) += ke(iIter->second,jIter->second);
//                //                        }
//                //                    }
//
//
//            }

//            A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());

//            A.makeCompressed();


