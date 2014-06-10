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
#include <model/Quadrature/Quadrature.h>


namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
	template <typename FiniteElementType, int dimMinusDomainDim, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct IntegrationDomain 
    {
        
        /**********************************************************************/
        IntegrationDomain()
        {
            assert(0 && "IntegrationDomain not implemented");
        }
    };

    
    /**************************************************************************/
    template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct IntegrationDomain<FiniteElementType,0,qOrder,QuadratureRule> : public std::deque<const typename FiniteElementType::ElementType*>
    {// Volume ntegration
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
    template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> : public std::deque<std::pair<const typename FiniteElementType::ElementType* const,int> >
    {// Boundary ntegration
        constexpr static int dim=FiniteElementType::dim;
        constexpr static int domainDim=dim-1;
        typedef Quadrature<domainDim,qOrder,QuadratureRule> QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;
        
        /**********************************************************************/
        const typename FiniteElementType::ElementType& element(const size_t& k) const
        {
            return *(this->operator[](k).first);
        }
    };


    
    
}	// close namespace
#endif