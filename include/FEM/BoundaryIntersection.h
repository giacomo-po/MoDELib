/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryIntersection_H_
#define model_BoundaryIntersection_H_

#include <vector>
#include <array>
#include <cfloat>
#include <Eigen/Dense>
#include <Simplex.h>
#include <IntegrationDomain.h>

namespace model
{
    template <typename A,typename B>
    class   BoundaryIntersection
    {
        const A a;
        const B b;
        
    public:
        
        /**********************************************************************/
        template <typename FiniteElementType>
        BoundaryIntersection(const FiniteElementType& fe) :
        /* init list */ a(fe),
        /* init list */ b(fe)
        {
            
        }
        
        /**********************************************************************/
        template <typename NodeType>
        bool operator()(const NodeType& node) const
        {
            return a(node) && b(node);
        }
        
    };
    
    
}	// close namespace
#endif

