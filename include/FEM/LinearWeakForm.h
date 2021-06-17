/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakForm_H_
#define model_LinearWeakForm_H_

#include <chrono>
#include <vector>
#include <Eigen/SparseCore>

//#include <EvalExpression.h>
//#include <TestExpression.h>
//#include <Quadrature.h>
//#include <AreSameType.h>
#include <TerminalColors.h>
#include <IntegrationDomain.h>
#include <LinearWeakExpression.h>
#include <LinearForm.h>
#include <LinearWeakSum.h>
#include <LinearWeakDiff.h>
#include <LinearWeakAssembler.h>
//#include <JGNselector.h>

#include <ExpressionRef.h>

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename _LinearFormType, typename _IntegrationDomainType>
	struct LinearWeakForm : public LinearWeakExpression<LinearWeakForm<_LinearFormType,_IntegrationDomainType> >
    {
        
        typedef _LinearFormType LinearFormType;
        typedef _IntegrationDomainType IntegrationDomainType;
        typedef LinearWeakForm<LinearFormType,IntegrationDomainType> LinearWeakFormType;

        
        typedef typename LinearFormType::TrialFunctionType TrialFunctionType;
        constexpr static int dim=IntegrationDomainType::dim;
        static_assert(dim==TypeTraits<TrialFunctionType>::dim, "dim MUST BE THE SAME");
        constexpr static int domainDim=IntegrationDomainType::domainDim;
        static_assert(domainDim==dim || domainDim==dim-1, "DOMAIN DIMENSIONALITY MUST BE EITHER dim (volume integration) OR dim-1 (boundary integration)");

        ExpressionRef<LinearFormType> linearForm;
//        const LinearFormType& linearForm;
        const IntegrationDomainType& domain;
        
        /**********************************************************************/
        LinearWeakForm(const LinearFormType& lF,
                       const IntegrationDomainType& dom) :
        /*init list */ linearForm(lF),
        /*init list */ domain(dom)
        {
             std::cout<<greenColor<<"Creating LinearWeakForm 1"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        LinearWeakForm(LinearFormType&& lF,
                       const IntegrationDomainType& dom) :
        /*init list */ linearForm(std::move(lF)),
        /*init list */ domain(dom)
        {
            std::cout<<greenColor<<"Creating LinearWeakForm 2"<<defaultColor<<std::endl;
        }
        
//        /**********************************************************************/
//        size_t gSize() const
//        {
//            return linearForm().testExp().gSize();
//        }

        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            return LinearWeakAssembler<LinearWeakFormType,dim-domainDim>(*this).assemble();
        }

    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Operator *
    template <typename T1, typename T2, typename FiniteElementType, int qOrder, int dimMinusDomainDim,
    /*     */ template <short unsigned int, size_t> class QuadratureRule>
    LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> > operator*(const LinearForm<T1,T2>& linearForm,
                             const IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule>& domain)
    {
        return LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> >(linearForm,domain);
    }
    
    template <typename T1, typename T2, typename FiniteElementType, int qOrder, int dimMinusDomainDim,
    /*     */ template <short unsigned int, size_t> class QuadratureRule>
    LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> > operator*(LinearForm<T1,T2>&& linearForm,
                                                                                                                              const IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule>& domain)
    {
        return LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> >(std::move(linearForm),domain);
    }
    
}	// close namespace
#endif





//        /**********************************************************************/
//        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
//        {
//
//            Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(linearForm.testExp.gSize()));
//            switch (domainDim)
//            {
//                case dim:
//                    assembleOnDomain(_globalVector);
//                    break;
//                case dim-1:
//                    assembleOnBoundary(_globalVector);
//                    break;
//
//                default:
//                    assert(0 && "CASE NOT IMPLEMENTED");
//                    break;
//            }
//
//            return _globalVector;
//        }
//
//        /**********************************************************************/
////		template <typename AbscissaType>
//        ElementVectorType domainAssemblyKernel(const AbscissaType& a, const ElementType& ele) const
//        {
//            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
////            return linearForm.testExp.sfm(ele,bary).transpose()*linearForm.evalExp(ele,bary)*ele.absJ(bary);
//            return linearForm(ele,bary)*ele.absJ(bary);
//		}
//
//        /**********************************************************************/
////		template <typename AbscissaType>
//        ElementVectorType boundaryAssemblyKernel(const AbscissaType& a1, const ElementType& ele, const int& boundaryFace) const
//        {
//            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
//            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::face2domainBary(b1,boundaryFace));
////            return linearForm.testExp.sfm(ele,bary).transpose()*linearForm.evalExp(ele,bary)*JGNselector<evalCols>::jGN(ele.jGN(bary,boundaryFace));
//            return linearForm(ele,bary)*JGNselector<evalCols>::jGN(ele.jGN(bary,boundaryFace));
//
//		}

