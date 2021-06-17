/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IntegrationDomain_H_
#define model_IntegrationDomain_H_

#include <deque>
#include <utility>      // std::pair, std::make_pair
#include <chrono>
#include <Quadrature.h>

#include <BarycentricTraits.h>
#include <IntegrationList.h>
#include <JGNselector.h>


namespace model
{
    
    
    /**************************************************************************/
	/**************************************************************************/
	template <typename FiniteElementType, int dimMinusDomainDim, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
	struct IntegrationDomain
    {
        
        /**********************************************************************/
        IntegrationDomain()
        {
            assert(0 && "IntegrationDomain not implemented");
        }
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename FiniteElementType, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
	struct IntegrationDomain<FiniteElementType,0,qOrder,QuadratureRule> : public std::deque<const typename FiniteElementType::ElementType*>
    {// Volume integration
        constexpr static int dim=FiniteElementType::dim;
        constexpr static int domainDim=dim;
        typedef Quadrature<domainDim,qOrder,QuadratureRule> QuadratureType;
        
        /**********************************************************************/
        const typename FiniteElementType::ElementType& element(const size_t& k) const
        {
            return *(this->operator[](k));
        }
    };

    /**************************************************************************/
    /**************************************************************************/
    template <typename FiniteElementType, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
	struct IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> : public std::deque<std::pair<const typename FiniteElementType::ElementType* const,int> >
    {// Boundary integration
        constexpr static int dim=FiniteElementType::dim;
        constexpr static int dimMinusDomainDim=1;
        constexpr static int domainDim=dim-dimMinusDomainDim;
        typedef Quadrature<domainDim,qOrder,QuadratureRule> QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> IntegrationDomainType;
        
        /**********************************************************************/
        const typename FiniteElementType::ElementType& element(const size_t& k) const
        {
            return *(this->operator[](k).first);
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        void integrate(const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const AbscissaType&,const ElementType&, const int&, const Args&...) const, const Args&...args) const
        {
//            std::cout<<"Integrating on boundary ("<<this->size()<<" faces) ..."<<std::flush;
//            const auto t0= std::chrono::system_clock::now();
            for (size_t k=0;k<this->size();++k)
            {
                QuadratureType::integrate(C,intgrl,mfp,*(this->operator[](k)).first,this->operator[](k).second,args...);
            }
//            std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
        }
        
//        /**********************************************************************/
//        double volumeKernel(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace) const
//        {
//            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
//            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
//            return JGNselector<1>::jGN(ele.jGN(bary,boundaryFace));
//        }
        
        /**********************************************************************/
        double volume() const
        {
            double temp=0.0;
            for(const auto& pair : *this)
            {
                temp+=pair.first->simplex.child(pair.second).vol0;
            }
//            integrate(this,temp,&IntegrationDomainType::volumeKernel);
            return temp;
        }
//        /**********************************************************************/
//        template <typename PointType>
//        void integrationList(std::deque<PointType>& deq) const
//        {
//            std::cout<<"populating integration List (size "<<deq.size()<<"->"<<std::flush;
//            const auto t0= std::chrono::system_clock::now();
//            for (int k=0;k<this->size();++k)
//            {
//                const ElementType& ele(*(this->operator[](k).first));  // element ID
//                const int f(this->operator[](k).second); //    face ID
//                //                ElementVectorType ve(ElementVectorType::Zero());
//                //                QuadratureType::integrate(this,ve,&LinearWeakFormType::boundaryAssemblyKernel<AbscissaType>,ele,f);
//                for(unsigned int q=0;q<qOrder;++q)
//                {
//                    const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(QuadratureType::abscissa(q)));
//                    const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::face2domainBary(b1,f));
//                    deq.emplace_back(ele.position(bary));
//                }
//                
//            }
//            std::cout<<deq.size()<<") ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
//            
//        }
        
        /**********************************************************************/
        template <typename PointType>
        IntegrationList<dimMinusDomainDim,PointType> integrationList() const
        {
            return IntegrationList<dimMinusDomainDim,PointType>(*this);
        }
        
    };
    
    
    
    
}	// close namespace
#endif
