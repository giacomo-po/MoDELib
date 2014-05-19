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
//#include <Eigen/Dense>
//#include <model/Mesh/Simplex.h>


namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
	template <int dim, int minusDim, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct IntegrationDomain 
    {
        typedef QuadratureRule<dim-minusDim,qOrder> QuadratureType;
        
        IntegrationDomain()
        {
            assert(0 && "IntegrationDomain not implemented");
        }
 //       static_assert(0,"IntegrationDomain not implemented");
        
    };

    
    template <int dim, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct IntegrationDomain<dim,0,qOrder,QuadratureRule> : public std::deque<int>
    {// Volume ntegration
     
        typedef QuadratureRule<dim,qOrder> QuadratureType;

        
    };
    
    template <int dim, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct IntegrationDomain<dim,1,qOrder,QuadratureRule> : public std::deque<std::pair<int,int> >
    {// Boundary ntegration
        
        typedef QuadratureRule<dim-1,qOrder> QuadratureType;
        
        
    };


    
    
}	// close namespace
#endif