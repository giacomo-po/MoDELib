/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryQuadraturePoint_H_
#define model_BoundaryQuadraturePoint_H_

#include <Eigen/Dense>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename ElementType,int rows,int cols>
    struct BoundaryQuadraturePoint : public Eigen::Matrix<double,rows,cols>
    {
        constexpr static int dim=ElementType::dim;
        const ElementType& ele;
        const Eigen::Matrix<double,dim+1,1> domainBary;
        const int boundaryFace;
        const double& weight;
        const Eigen::Matrix<double,dim,1> P;
        
        /**********************************************************************/
        BoundaryQuadraturePoint(const ElementType& _ele,
                                 const Eigen::Matrix<double,dim+1,1>& _domainBary,
                                 const int& _boundaryFace,
                                 const double& _weight) :
        /* init */ Eigen::Matrix<double,rows,cols>(Eigen::Matrix<double,rows,cols>::Zero()),
        /* init */ ele(_ele),
        /* init */ domainBary(_domainBary),
        /* init */ boundaryFace(_boundaryFace),
        /* init */ weight(_weight),
        /* init */ P(ele.position(domainBary))
        {

        }
        
    };
    
}	// close namespace
#endif
