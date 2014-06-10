/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakForm_H_
#define model_BilinearWeakForm_H_

//#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <chrono>
#include <vector>
#include <iterator>     // std::advance
#include <utility>      // std::move

#ifdef _OPENMP
#include <omp.h>
#include <Model/Threads/EqualIteratorRange.h>
#endif

#include <Eigen/Sparse>

#include <model/FEM/TrialOperators/TrialExpressionBase.h>
#include <model/FEM/WeakForms/LinearWeakForm.h>
#include <model/FEM/WeakForms/WeakProblem.h>
#include <model/Quadrature/Quadrature.h>
#include <model/Utilities/AreSameType.h>
#include <model/Utilities/TerminalColors.h>
#include <model/FEM/WeakForms/ElementaryDomain.h>



namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T1,typename T2>
	class BilinearWeakForm
    {
        
        static_assert(AreSameType<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE CREATING A BilinearWeakForm OF DIFFERENT TRIALFUNCTIONS.");
        static_assert((T1::rows-T2::rows)==0,"YOU ARE CREATING A BilinearWeakForm BETWEEN A TRIALEXPRESSION AND A TESTEXPRESSION WITH DIFFERENT NUMBER OF ROWS");
        
        typedef typename T2::TrialFunctionType TF;
        typedef BilinearWeakForm<T1,T2> BilinearWeakFormType;
        typedef typename TypeTraits<TF>::ElementType ElementType;
        typedef typename TypeTraits<TF>::FiniteElementType FiniteElementType;
        
        /**********************************************************************/
        BilinearWeakForm(const BilinearWeakFormType&) = default; // prevent copy (too expensive)
        
        /**********************************************************************/
        BilinearWeakForm& operator=(const BilinearWeakFormType&) = default; // prevent assignment (too expensive)
        
    public:
        
        enum {dim=TypeTraits<TF>::dim};
        enum {nodesPerElement=TypeTraits<TF>::nodesPerElement};
        enum {dofPerNode=TypeTraits<TF>::dofPerNode};
        enum {dofPerElement=TypeTraits<TF>::dofPerElement};
        
        
        typedef Eigen::Matrix<double,dofPerElement,dofPerElement> ElementMatrixType;
        
        const T1  testExp;
        const T2 trialExp;
        
        const size_t gSize;
        
        std::vector<Eigen::Triplet<double> > globalTriplets; // this should be private
        
        
        double maxAbsValue; // this should be private
        
        
        
        /**********************************************************************/
        BilinearWeakForm(const TestExpression<T1>& testE, const TrialExpressionBase<T2>& trialE) :
        /* init list */ testExp(testE), // cast testE to its base T2 type
        /* init list */ trialExp(trialE.derived()), // cast trialE to its derived T1 type
        /* init list */ gSize(trialExp.nodeSize()*dofPerNode),
        /* init list */ maxAbsValue(0.0)
        {
            
            std::cout<<greenColor<<"Creating BilinearWeakForm: gSize="<<gSize<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        BilinearWeakForm(BilinearWeakForm&&) = default; // explicit move constructor (required since copy constructor is private)
        
        /**********************************************************************/
        BilinearWeakForm& operator=(BilinearWeakForm&&) = default; // explicit a move assignment (required since assignment is private)
        
        
        //        /**********************************************************************/
        //        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        //        const BilinearWeakFormType& operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
        //        {
        //            assembleOnDomain<qOrder,QuadratureRule>();
        //            return *this;
        //        }
        
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        BilinearWeakFormType&& operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
        {/*
          * See section "Returning an explicit rvalue-reference from a function" in
          * http://www.cprogramming.com/c++11/rvalue-references-and-move-semantics-in-c++11.html
          */
            assembleOnDomain<qOrder,QuadratureRule>();
            return std::move(*this);
        }
        
        //        /**********************************************************************/
        //        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        //        BilinearWeakFormType operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
        //        {
        //            assembleOnDomain<qOrder,QuadratureRule>();
        //            return *this;
        //        }
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        void assembleOnDomain()
        {
            
            
//            CHANGE THIS CLASS AS THE LINEARWEAKEXPRESSION, AND THEN CHANGE THIS FUNCTION to
//            RETURN THE GLOBAL TRIPLETS.
            
            typedef Quadrature<TF::dim,qOrder,QuadratureRule> QuadratureType;
            typedef typename QuadratureType::VectorDim AbscissaType;
            
            //            double t0(clock());
            auto t0= std::chrono::system_clock::now();
            
            
            std::cout<<"Assembling on domain..."<<std::flush;

            globalTriplets.clear();
            globalTriplets.reserve(dofPerElement*dofPerElement*trialExp.elementSize());
            maxAbsValue=0.0;

            
            
//#ifdef _OPENMP
//            // Equally divide iterator range
//            EqualIteratorRange<FiniteElementType::ElementContainerType> eir(trialExp.elementBegin(),trialExp.elementEnd(),omp_get_max_threads());
//            
//#pragma omp parallel for
//            for (unsigned int k = 0; k < eir.size(); ++k)
//            {
//                for (typename FiniteElementType::ElementContainerType::const_iterator eIter =eir[k].first;
//                     /*                                                            */ eIter!=eir[k].second;
//                     /*                                                            */ ++eIter)
//                {
//                    ElementMatrixType ke(ElementMatrixType::Zero());
//                    QuadratureType::integrate(this,ke,&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,eIter->second);
//
//                    std::map<int,const int> iMap;
//                    for (int i=0;i<dofPerElement;++i)
//                    {
//                        const size_t  nodeID_I(i/dofPerNode);
//                        const size_t nodeDof_I(i%dofPerNode);
//                        const size_t gI= eIter->second.node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                        iMap.insert(gI,i);
//                    }
//                    
//                    std::map<int,const int> jMap;
//                    for (int j=0;j<dofPerElement;++j)
//                    {
//                        const size_t  nodeID_J(j/dofPerNode);
//                        const size_t nodeDof_J(j%dofPerNode);
//                        const size_t gJ=eIter->second.node(nodeID_J).gID*dofPerNode+nodeDof_J;
//                    }
//                    
//                    for (std::map<int,const int>::const_iterator jIter=jMap.begin();jIter!=jMap.end();++jIter)
//                        for (std::map<int,const int>::const_iterator iIter=iMap.begin();iIter!=iMap.end();++iIter)
//                    {
//                        {
////                            A.coeffRef(iIter->first,jIter->first) += ke(iIter->second,jIter->second);
//                        }
//                    }
//                
//                }
//                
//                FINISH HERE
//                
//                // For each Element:
//                // 1 Compute dense element matrix
//                // 2 Find globad row and column IDs
//                // 3 Insert each element in a global sparse matrix (for smart insertion see http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html)
//                // keep summing the globl matrices
//                
//                // OR
//                
//                // For each Element:
//                // 1 Compute dense element matrix
//                // 2 Find globad row and column IDs. Create two map of type (gI,i) and (gJ,j), so that gI and gJ are sorted, see below
//                // 3 Insert each element in a global sparse matrix using insert(gI,gJ). For a column-major matrix (Eigen default)
//                // this requires looping over gJ first (inner vector) and then loop for increasing gI (inner index).
//                // (see http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html)
//                // keep summing the globl matrices
//                
//            }
//            
//#else
//            Eigen::SparseMatrix<double> A(gSize,gSize);

            
            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =trialExp.elementBegin();
                 /*                                                            */ eIter!=trialExp.elementEnd();
                 /*                                                            */ ++eIter)
            {
                // Compute element stiffness Matrix
                ElementMatrixType ke(ElementMatrixType::Zero());
                QuadratureType::integrate(this,ke,&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,eIter->second);
                
                //Assemble into global array of triplets
                for (int i=0;i<dofPerElement;++i)
                {
                    for (int j=0;j<dofPerElement;++j)
                    {
                        if(ke(i,j)!=0.0)
                        {
                            const size_t  nodeID_I(i/dofPerNode);
                            const size_t nodeDof_I(i%dofPerNode);
                            const size_t gI= eIter->second.node(nodeID_I).gID*dofPerNode+nodeDof_I;
                            
                            const size_t  nodeID_J(j/dofPerNode);
                            const size_t nodeDof_J(j%dofPerNode);
                            const size_t gJ=eIter->second.node(nodeID_J).gID*dofPerNode+nodeDof_J;
                            
                            globalTriplets.emplace_back(gI,gJ,ke(i,j));
                            
                            //                            A.coeffRef(gI,gJ) += ke(i,j);
                            
                            if (std::fabs(ke(i,j))>maxAbsValue)
                            {
                                maxAbsValue=std::fabs(ke(i,j));
                            }
                        }
                    }
                }
                
//                    std::map<int,const int> iMap;
//                    for (int i=0;i<dofPerElement;++i)
//                    {
//                        const size_t  nodeID_I(i/dofPerNode);
//                        const size_t nodeDof_I(i%dofPerNode);
//                        const size_t gI= eIter->second.node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                        iMap.emplace(gI,i);
//                    }
//
//                    std::map<int,const int> jMap;
//                    for (int j=0;j<dofPerElement;++j)
//                    {
//                        const size_t  nodeID_J(j/dofPerNode);
//                        const size_t nodeDof_J(j%dofPerNode);
//                        const size_t gJ=eIter->second.node(nodeID_J).gID*dofPerNode+nodeDof_J;
//                        jMap.emplace(gJ,j);
//                    }
//
//                    for (std::map<int,const int>::const_iterator jIter=jMap.begin();jIter!=jMap.end();++jIter)
//                        for (std::map<int,const int>::const_iterator iIter=iMap.begin();iIter!=iMap.end();++iIter)
//                    {
//                        {
//                            A.insert(iIter->first,jIter->first) += ke(iIter->second,jIter->second);
//                        }
//                    }

                
            }
            
//            A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());
            
//            A.makeCompressed();
            
//#endif
            std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            //            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
        }
        
        /**********************************************************************/
		template <typename AbscissaType>
        ElementMatrixType elementAssemblyKernel(const AbscissaType& a, const ElementType& ele) const
        {
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
            return testExp.sfm(ele,bary).transpose()*trialExp.sfm(ele,bary)*ele.absJ(bary);
		}
        
        
        
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
        
        /**********************************************************************/
        template< typename T3>
        WeakProblem<BilinearWeakFormType,T3> operator=(const LinearWeakExpression<T3>& lwf) const
        {
            return WeakProblem<BilinearWeakFormType,T3>(*this,lwf);
        }
        
    };
    
    
    template <typename T1,typename T2>
    BilinearWeakForm<T1,T2> operator, (const TestExpression<T1>& testE, const TrialExpressionBase<T2>& trialE)
    {
        return BilinearWeakForm<T1,T2>(testE,trialE);
    }
    
}	// close namespace
#endif

