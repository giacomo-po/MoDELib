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

#include <model/FEM/TrialOperators/EvalExpression.h>
#include <model/FEM/TrialOperators/TestExpression.h>
#include <model/Quadrature/Quadrature.h>
#include <model/Utilities/AreSameType.h>
#include <model/Utilities/TerminalColors.h>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/FEM/WeakForms/LinearWeakExpression.h>
#include <model/FEM/WeakForms/LinearForm.h>
#include <model/FEM/WeakForms/LinearWeakSum.h>
#include <model/FEM/WeakForms/LinearWeakDiff.h>
#include <model/FEM/WeakForms/JGNselector.h>
#include <model/MPI/MPIcout.h>


namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename _LinearFormType, typename _IntegrationDomainType>
	class LinearWeakForm : public LinearWeakExpression<LinearWeakForm<_LinearFormType,_IntegrationDomainType> >
    {
        
    public:
        
        typedef _LinearFormType LinearFormType;
        typedef typename LinearFormType::TrialFunctionType TrialFunctionType;
        typedef _IntegrationDomainType IntegrationDomainType;

        constexpr static int dim=IntegrationDomainType::dim;
        static_assert(dim==TypeTraits<TrialFunctionType>::dim, "dim MUST BE THE SAME");
        constexpr static int domainDim=IntegrationDomainType::domainDim;
        static_assert(domainDim==dim || domainDim==dim-1, "DOMAIN DIMENSIONALITY MUST BE EITHER dim (volume integration) OR dim-1 (boundary integration)");

    private:
        typedef LinearWeakForm<LinearFormType,IntegrationDomainType> LinearWeakFormType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        typedef Eigen::Matrix<double,dofPerElement,1> ElementVectorType;
        constexpr static int evalCols=LinearFormType::evalCols;
        typedef typename IntegrationDomainType::QuadratureType QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;
        
        /**********************************************************************/
        void assembleOnDomain(Eigen::Matrix<double,Eigen::Dynamic,1>& _globalVector) const
        {

            assert(0 && "FINISH HERE, AbscissaType is the wrong type in this function in case of boundary integration (and viceversa).");
//             model::cout<<"Assembling LinearWeakForm on domain..."<<std::flush;
//            const auto t0= std::chrono::system_clock::now();
//            for (int k=0;k<domain.size();++k)
//            {
//                // Integrate element vector
//                ElementVectorType ve(ElementVectorType::Zero());
//                QuadratureType::integrate(this,ve,&LinearWeakFormType::template domainAssemblyKernel<AbscissaType>,domain.element(k));
//                
//                //Assemble element vector into global vector
//                for (int i=0;i<dofPerElement;++i)
//                {
//                    const size_t  nodeID_I(i/dofPerNode);
//                    const size_t nodeDof_I(i%dofPerNode);
//                    const size_t gI= domain.element(k).node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                    _globalVector(gI) += ve(i);
//                }
//            }
//             model::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
        }
        
        /**********************************************************************/
        void assembleOnBoundary(Eigen::Matrix<double,Eigen::Dynamic,1>& _globalVector) const
        {
             model::cout<<"Assembling LinearWeakForm on faces ("<<domain.size()<<" faces) ..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            for (int k=0;k<domain.size();++k)
            {
                const ElementType& ele(domain.element(k));  // element ID
                const int f(domain[k].second); //    face ID
                ElementVectorType ve(ElementVectorType::Zero());
                QuadratureType::integrate(this,ve,&LinearWeakFormType::boundaryAssemblyKernel<AbscissaType>,ele,f);
                for (int i=0;i<dofPerElement;++i)
                {
                    const size_t  nodeID_I(i/dofPerNode);
                    const size_t nodeDof_I(i%dofPerNode);
                    const size_t gI= ele.node(nodeID_I).gID*dofPerNode+nodeDof_I;
                    _globalVector(gI) += ve(i);
                }
            }
             model::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
        }
        
    public:
        
        const LinearFormType linearForm;
        const IntegrationDomainType& domain;
        
        /**********************************************************************/
        LinearWeakForm(const LinearFormType& lF,
                       const IntegrationDomainType& dom) :
        /*init list */ linearForm(lF),
        /*init list */ domain(dom)
        {
             model::cout<<greenColor<<"Creating LinearWeakForm "<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        size_t gSize() const
        {
            return linearForm.testExp.gSize();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            
            Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(linearForm.testExp.gSize()));
            switch (domainDim)
            {
                case dim:
                    assembleOnDomain(_globalVector);
                    break;
                case dim-1:
                    assembleOnBoundary(_globalVector);
                    break;
                    
                default:
                    assert(0 && "CASE NOT IMPLEMENTED");
                    break;
            }
            
            return _globalVector;
        }
        
        /**********************************************************************/
		template <typename AbscissaType>
        ElementVectorType domainAssemblyKernel(const AbscissaType& a, const ElementType& ele) const
        {
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
//            return linearForm.testExp.sfm(ele,bary).transpose()*linearForm.evalExp(ele,bary)*ele.absJ(bary);
            return linearForm(ele,bary)*ele.absJ(bary);
		}
        
        /**********************************************************************/
		template <typename AbscissaType>
        ElementVectorType boundaryAssemblyKernel(const AbscissaType& a1, const ElementType& ele, const int& boundaryFace) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::face2domainBary(b1,boundaryFace));
//            return linearForm.testExp.sfm(ele,bary).transpose()*linearForm.evalExp(ele,bary)*JGNselector<evalCols>::jGN(ele.jGN(bary,boundaryFace));
            return linearForm(ele,bary)*JGNselector<evalCols>::jGN(ele.jGN(bary,boundaryFace));

		}
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Operator *
    template <typename T1, typename T2, typename FiniteElementType, int qOrder, int dimMinusDomainDim,
    /*     */ template <short unsigned int, short unsigned int> class QuadratureRule>
    LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> > operator*(const LinearForm<T1,T2>& linearForm,
                             const IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule>& domain)
    {
        return LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FiniteElementType,dimMinusDomainDim,qOrder,QuadratureRule> >(linearForm,domain);
    }
    
}	// close namespace
#endif
