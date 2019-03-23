/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryIntegrationPoint_H_
#define model_BoundaryIntegrationPoint_H_

#include <Eigen/Dense>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename ElementType>
	struct BoundaryIntegrationPoint
    {
        constexpr static int dim=ElementType::dim;
        const ElementType& ele;
        const Eigen::Matrix<double,dim+1,1> domainBary;
        const int boundaryFace;
        const double& weight;
        
        /**********************************************************************/
        BoundaryIntegrationPoint(const ElementType& _ele,
                                 const Eigen::Matrix<double,dim+1,1>& _domainBary,
                                 const int& _boundaryFace,
                                 const double& _weight) :
        /*                    */ ele(_ele),
        /*                    */ domainBary(_domainBary),
        /*                    */ boundaryFace(_boundaryFace),
        /*                    */ weight(_weight)
        {

        }
        
    };
    
}	// close namespace
#endif