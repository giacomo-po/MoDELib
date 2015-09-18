/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryNearDislocations_H_
#define model_BoundaryNearDislocations_H_

#include <deque>
#include <utility>      // std::pair, std::make_pair
#include <Eigen/Dense>
#include <model/Mesh/Simplex.h>
#include <model/FEM/Domains/IntegrationDomain.h>

namespace model
{
    
 
    /**************************************************************************/
	/**************************************************************************/
	struct BoundaryNearDislocations //: public IntegrationDomain<dim,1,qOrder,QuadratureRule>
    {
        
        template <typename FiniteElementType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        static IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary(const FiniteElementType& fe)
        {
            IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> temp;

            for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
                 /*                                                            */ eIter!=fe.elementEnd();
                 /*                                                            */ eIter++)
            {
                if(eIter->second.isBoundaryElement())
                {
                    
                    HERE CREATE A PARTICLE AND CHECK WHAT IS THE NUMBER OF DISLOCATION PARTICLES NEAR IT!!!!
                    
                    const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
                    for (int f=0;f<boundaryFaces.size();++f)
                    {
                        temp.emplace_back(&eIter->second,boundaryFaces[f]);
                    }
                }
            }
            
            return temp;
        }
        
    };

}	// close namespace
#endif