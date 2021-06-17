/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakAssembler_H_
#define model_LinearWeakAssembler_H_

#include <chrono>
#include <vector>
//#include <Eigen/SparseCore>
#include <Eigen/Dense>

#include <Quadrature.h>
//#include <AreSameType.h>
#include <TerminalColors.h>
#include <TrialBase.h>
//#include <IntegrationDomain.h>
//#include <LinearWeakExpression.h>
//#include <LinearForm.h>
//#include <LinearWeakSum.h>
//#include <LinearWeakDiff.h>
#include <JGNselector.h>



namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename LinearWeakFormType, short unsigned int DimMinusDomainDim>
    struct LinearWeakAssembler
    {
    
        const LinearWeakFormType& lwf;
        
        /**********************************************************************/
        LinearWeakAssembler(const LinearWeakFormType& lwf_in) :
        /* init list */ lwf(lwf_in)
        {
            assert(0 && "LinearWeakAssembler: DimMinusDomainDim must be 0 or 1.");
        }
        
    };
    

    /**************************************************************************/
    /**************************************************************************/
    template<typename LinearWeakFormType>
    struct LinearWeakAssembler<LinearWeakFormType,0>
    {
        typedef LinearWeakAssembler<LinearWeakFormType,0> LinearWeakAssemblerType;
        typedef typename LinearWeakFormType::LinearFormType LinearFormType;
        typedef typename LinearWeakFormType::IntegrationDomainType IntegrationDomainType;
        constexpr static int dim=IntegrationDomainType::dim;
        
        typedef typename LinearFormType::TrialFunctionType TrialFunctionType;
        
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        typedef Eigen::Matrix<double,dofPerElement,1> ElementVectorType;
        constexpr static int evalCols=LinearFormType::evalCols;
        typedef typename IntegrationDomainType::QuadratureType QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;
        
        
        const LinearWeakFormType& lwf;
        
        /**********************************************************************/
        LinearWeakAssembler(const LinearWeakFormType& lwf_in) :
        /* init list */ lwf(lwf_in)
        {
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> assemble() const
        {
            std::cout<<"Assembling LinearWeakForm on domain ("<<lwf.domain.size()<<" elements) ..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
//            std::cout<<"lwf.gSize()="<<lwf.gSize()<<std::endl;

//            std::cout<<"lwf.domain.size()="<<lwf.domain.size()<<std::endl;
////            std::cout<<"lwf.linearForm.testExp.gSize()="<<lwf.linearForm.testExp.gSize()<<std::endl;
//
//            std::cout<<"lwf.linearForm.evalExp.rows="<<lwf.linearForm.evalExp.rows<<std::endl;
//
//            std::cout<<"lwf.linearForm.testExp.trial.gSize()"<<lwf.linearForm.testExp.trial().gSize()<<std::endl;
            
            Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(TrialBase<TrialFunctionType>::gSize()));
            
            for (size_t k=0;k<lwf.domain.size();++k)
            {
                // Integrate element vector
                ElementVectorType ve(ElementVectorType::Zero());
                QuadratureType::integrate(this,ve,&LinearWeakAssemblerType::assemblyKernel,lwf.domain.element(k));
                
                //Assemble element vector into global vector
                for (int i=0;i<dofPerElement;++i)
                {
                    const size_t  nodeID_I(i/dofPerNode);
                    const size_t nodeDof_I(i%dofPerNode);
                    const size_t gI= lwf.domain.element(k).node(nodeID_I).gID*dofPerNode+nodeDof_I;
                    _globalVector(gI) += ve(i);
                }
            }
            std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;

            return _globalVector;
        }
        
        /**********************************************************************/
        ElementVectorType assemblyKernel(const AbscissaType& a, const ElementType& ele) const
        {
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
            return lwf.linearForm()(ele,bary)*ele.absJ(bary);
        }

        
    };
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename LinearWeakFormType>
    struct LinearWeakAssembler<LinearWeakFormType,1>
    {
        typedef LinearWeakAssembler<LinearWeakFormType,1> LinearWeakAssemblerType;
        typedef typename LinearWeakFormType::LinearFormType LinearFormType;
        typedef typename LinearWeakFormType::IntegrationDomainType IntegrationDomainType;
        constexpr static int dim=IntegrationDomainType::dim;

        typedef typename LinearFormType::TrialFunctionType TrialFunctionType;
        
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        typedef Eigen::Matrix<double,dofPerElement,1> ElementVectorType;
        constexpr static int evalCols=LinearFormType::evalCols;
        typedef typename IntegrationDomainType::QuadratureType QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;

        
        const LinearWeakFormType& lwf;
        
        /**********************************************************************/
        LinearWeakAssembler(const LinearWeakFormType& lwf_in) :
        /* init list */ lwf(lwf_in)
        {
        
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> assemble() const
        {
            std::cout<<"Assembling LinearWeakForm on faces ("<<lwf.domain.size()<<" faces) ..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(TrialBase<TrialFunctionType>::gSize()));

            for (size_t k=0;k<lwf.domain.size();++k)
            {
                const ElementType& ele(lwf.domain.element(k));  // element ID
                const int f(lwf.domain[k].second); //    face ID
                ElementVectorType ve(ElementVectorType::Zero());
                QuadratureType::integrate(this,ve,&LinearWeakAssemblerType::assemblyKernel,ele,f);
                for (int i=0;i<dofPerElement;++i)
                {
                    const size_t  nodeID_I(i/dofPerNode);
                    const size_t nodeDof_I(i%dofPerNode);
                    const size_t gI= ele.node(nodeID_I).gID*dofPerNode+nodeDof_I;
                    _globalVector(gI) += ve(i);
                }
            }
            std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            return _globalVector;
        }
        
        /**********************************************************************/
        ElementVectorType assemblyKernel(const AbscissaType& a1,
                            const ElementType& ele,
                            const int& boundaryFace) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::face2domainBary(b1,boundaryFace));
            return lwf.linearForm()(ele,bary)*JGNselector<evalCols>::jGN(ele.jGN(bary,boundaryFace));
        }
        
    };
    
    
}	// close namespace
#endif
