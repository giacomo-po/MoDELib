/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_OnExternalBoundary_H_
#define model_OnExternalBoundary_H_

#include <vector>
#include <array>
#include <cfloat>
#include <Eigen/Dense>
#include <Simplex.h>
#include <IntegrationDomain.h>

namespace model
{
    class   OnExternalBoundary
    {
        
    public:
        
        /**************************************/
        template <typename FiniteElementType>
        OnExternalBoundary(const FiniteElementType&)
        {
            
        }
        
        /**************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return node.outNormal().squaredNorm()>FLT_EPSILON;
        }
        
        /**********************************************************************/
        template <typename FiniteElementType, size_t qOrder, template <short unsigned int, size_t> class QuadratureRule>
        static IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary(const FiniteElementType& fe)
        {
            IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> temp;
            for (const auto& eIter : fe.elements())
            {
                if(eIter.second.isBoundaryElement())
                {
                    const std::vector<int> boundaryFaces=eIter.second.boundaryFaces();
                    for (size_t f=0;f<boundaryFaces.size();++f)
                    {
                        temp.emplace_back(&(eIter.second),boundaryFaces[f]);
                    }
                }
            }
            
            return temp;
        }
        
        
    };
    
    
}	// close namespace
#endif

