/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IntegrationList_H_
#define model_IntegrationList_H_

#include <deque>
#include <utility>      // std::pair, std::make_pair
#include <chrono>
#include <model/Quadrature/Quadrature.h>
#include <model/MPI/MPIcout.h>
#include <model/FEM/BarycentricTraits.h>
//#include <model/FEM/TrialOperators/TestExpression.h>
//#include <model/FEM/Domains/LinearWeakList.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dimMinusDomainDim, typename PointType>
	struct IntegrationList;
    
    template <typename PointType>
	struct IntegrationList<0,PointType> : public std::deque<PointType,Eigen::aligned_allocator<PointType>>
    {
//        static_assert(0,"NOT IMPLEMENTED YET");
    };
    
    
    template <typename PointType>
	struct IntegrationList<1,PointType> : public std::deque<PointType,Eigen::aligned_allocator<PointType>>
    {
        
        /**********************************************************************/
        template<typename IntegrationDomainType>
        IntegrationList(const IntegrationDomainType& domain)
        {
            
            typedef typename IntegrationDomainType::ElementType ElementType;
            typedef typename IntegrationDomainType::QuadratureType QuadratureType;
            
            model::cout<<"Creating IntegrationList (size "<<this->size()<<"->"<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            for (size_t k=0;k<domain.size();++k)
            {
                const ElementType& ele(*domain[k].first);  // element ID
                const int f(domain[k].second); //    face ID
                for(unsigned int q=0;q<QuadratureType::quadratureOrder;++q)
                {
                    const Eigen::Matrix<double,ElementType::dim,1>   faceBary(BarycentricTraits<ElementType::dim-1>::x2l(QuadratureType::abscissa(q)));
                    const Eigen::Matrix<double,ElementType::dim+1,1> domainBary(BarycentricTraits<ElementType::dim>::face2domainBary(faceBary,f));
                    this->emplace_back(ele,domainBary,f,QuadratureType::weight(q));
                }
                
            }
            model::cout<<this->size()<<") ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
        }
        
    };
    
    template <typename PointType>
    struct IntegrationList<2,PointType> : public std::deque<PointType,Eigen::aligned_allocator<PointType>>
    {
        
//        /**********************************************************************/
//        template<typename IntegrationDomainType>
//        IntegrationList(const IntegrationDomainType& domain)
//        {
//            
//            typedef typename IntegrationDomainType::ElementType ElementType;
//            typedef typename IntegrationDomainType::QuadratureType QuadratureType;
//            
//            model::cout<<"Creating IntegrationList (size "<<this->size()<<"->"<<std::flush;
//            const auto t0= std::chrono::system_clock::now();
//            for (size_t k=0;k<domain.size();++k)
//            {
//                const ElementType& ele(*domain[k].first);  // element ID
//                const int f(domain[k].second); //    face ID
//                for(unsigned int q=0;q<QuadratureType::quadratureOrder;++q)
//                {
//                    const Eigen::Matrix<double,ElementType::dim,1>   faceBary(BarycentricTraits<ElementType::dim-1>::x2l(QuadratureType::abscissa(q)));
//                    const Eigen::Matrix<double,ElementType::dim+1,1> domainBary(BarycentricTraits<ElementType::dim>::face2domainBary(faceBary,f));
//                    this->emplace_back(ele,domainBary,f,QuadratureType::weight(q));
//                }
//                
//            }
//            model::cout<<this->size()<<") ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
//        }
        
    };
    
    

    
    
}	// close namespace
#endif
