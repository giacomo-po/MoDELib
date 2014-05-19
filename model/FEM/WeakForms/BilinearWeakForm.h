/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakForm_H_
#define model_BilinearWeakForm_H_

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <vector>
#include <Eigen/SparseCore>

#include <model/FEM/TrialOperators/TrialExpressionBase.h>
#include <model/FEM/WeakForms/LinearWeakForm.h>
#include <model/FEM/WeakForms/WeakProblem.h>
#include <model/Quadrature/Quadrature.h>
#include <model/Utilities/AreSameType.h>
#include <model/Utilities/TerminalColors.h>
#include <model/FEM/WeakForms/ElementaryDomain.h>


#ifdef _OPENMP
#include <omp.h>
//#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T1,typename T2>
	class BilinearWeakForm
    {
        
        static_assert(AreSameType<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE CREATING A BilinearWeakForm OF DIFFERENT TRIALFUNCTIONS.");
        static_assert((T1::rows-T2::rows)==0,"YOU ARE CREATING A BilinearWeakForm BETWEEN A TRIALEXPRESSION AND A TESTEXPRESSION WITH DIFFERENT NUMBER OF ROWS");
        
        typedef typename T1::TrialFunctionType TF;
        typedef BilinearWeakForm<T1,T2> BilinearWeakFormType;
        typedef typename TypeTraits<TF>::ElementType ElementType;
        
        
        
//        /**********************************************************************/
//        void assembleTriplets()
//        {
//        
//        }
        
        
    public:
        
        enum {dim=TypeTraits<TF>::dim};
        enum {nodesPerElement=TypeTraits<TF>::nodesPerElement};
        enum {dofPerNode=TypeTraits<TF>::dofPerNode};
        enum {dofPerElement=TypeTraits<TF>::dofPerElement};
        
        
        typedef Eigen::Matrix<double,dofPerElement,dofPerElement> ElementMatrixType;
        
        const T2  testExp;
        const T1 trialExp;

        const size_t gSize;
        
        std::vector<Eigen::Triplet<double> > globalTriplets; // this should be private
        
        
        double maxAbsValue; // this should be private
        

        
        /**********************************************************************/
        BilinearWeakForm(const TestExpression<T2>& testE, const TrialExpressionBase<T1>& trialE) :
        /* init list */ testExp(testE), // cast testE to its base T2 type
        /* init list */ trialExp(trialE.derived()), // cast trialE to its derived T1 type
        /* init list */ gSize(trialExp.nodeSize()*dofPerNode),
        /* init list */ maxAbsValue(0.0)
        {
            
            std::cout<<greenColor<<"Creating BilinearWeakForm: gSize="<<gSize<<defaultColor<<std::endl;

        }
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        const BilinearWeakFormType& operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
        {
            assembleOnDomain<qOrder,QuadratureRule>();
            return *this;
        }
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        void assembleOnDomain()
        {
            typedef Quadrature<TF::dim,qOrder,QuadratureRule> QuadratureType;
            typedef typename QuadratureType::VectorDim AbscissaType;
            
            double t0(clock());
            std::cout<<"Assembling on domain..."<<std::flush;
            
            globalTriplets.clear();
            maxAbsValue=0.0;
            globalTriplets.reserve(dofPerElement*dofPerElement*trialExp.elementSize());
            
#ifdef _OPENMP
            std::cout<<"\n    computing element matrices..."<<std::flush;
            
            std::vector<ElementMatrixType> kv;
            kv.resize(trialExp.elementSize());
            
#pragma omp parallel for
            for (int n=0;n<trialExp.elementSize();++n)
            {
                // Compute element stiffness Matrix
                ElementMatrixType ke(ElementMatrixType::Zero());
                kv[n].setZero();
                QuadratureType::integrate(this,kv[n],&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,trialExp.element(n));
            }
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            
            // Assemble in globalTriplets
            t0=clock();
            std::cout<<"    assembling in global vector of triplets..."<<std::flush;
            
            for (int n=0;n<trialExp.elementSize();++n)
            {
                //Assemble into global array of triplets
                for (int i=0;i<dofPerElement;++i)
                {
                    for (int j=0;j<dofPerElement;++j)
                    {
                        if(kv[n](i,j)!=0.0)
                        {
                            const size_t  nodeID_I(i/dofPerNode);
                            const size_t nodeDof_I(i%dofPerNode);
                            const size_t gI= testExp.element(n).node(nodeID_I).gID*dofPerNode+nodeDof_I;
                            
                            const size_t  nodeID_J(j/dofPerNode);
                            const size_t nodeDof_J(j%dofPerNode);
                            const size_t gJ=trialExp.element(n).node(nodeID_J).gID*dofPerNode+nodeDof_J;
                            
                            globalTriplets.emplace_back(gI,gJ,kv[n](i,j));
                            
                            if (std::fabs(kv[n](i,j))>maxAbsValue)
                            {
                                maxAbsValue=std::fabs(kv[n](i,j));
                            }
                        }
                    }
                }
            }
            
#else
//            A.resize(gSize,gSize);
//            A.reserve(Eigen::VectorXi::Constant(gSize,1000));

            
            for (int n=0;n<trialExp.elementSize();++n)
            {
                // Compute element stiffness Matrix
                ElementMatrixType ke(ElementMatrixType::Zero());
                QuadratureType::integrate(this,ke,&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,trialExp.element(n));
                
                //Assemble into global array of triplets
                for (int i=0;i<dofPerElement;++i)
                {
                    for (int j=0;j<dofPerElement;++j)
                    {
                        if(ke(i,j)!=0.0)
                        {
                            const size_t  nodeID_I(i/dofPerNode);
                            const size_t nodeDof_I(i%dofPerNode);
                            const size_t gI= testExp.element(n).node(nodeID_I).gID*dofPerNode+nodeDof_I;
                            
                            const size_t  nodeID_J(j/dofPerNode);
                            const size_t nodeDof_J(j%dofPerNode);
                            const size_t gJ=trialExp.element(n).node(nodeID_J).gID*dofPerNode+nodeDof_J;
                            
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
            
//            A.makeCompressed();
            
#endif
            
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
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
        
        /**********************************************************************/
        template< typename T3,typename T4>
        WeakProblem<BilinearWeakFormType,LinearWeakForm<T3,T4> > operator=(const LinearWeakForm<T3,T4>& lwf) const
        {
            return WeakProblem<BilinearWeakFormType,LinearWeakForm<T3,T4> >(*this,lwf);
        }
        
    };
    
    
    template <typename T1,typename T2>
    BilinearWeakForm<T1,T2> operator, (const TestExpression<T2>& testE, const TrialExpressionBase<T1>& trialE)
    {
        return BilinearWeakForm<T1,T2>(testE,trialE);
    }
    
}	// close namespace
#endif

