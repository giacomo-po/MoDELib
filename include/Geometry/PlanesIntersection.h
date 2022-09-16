/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanesIntersection_H_
#define model_PlanesIntersection_H_

#include <cfloat>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace model
{
    
    template <int dim>
    struct PlanesIntersection
    {
        
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixDimDynamic;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::JacobiSVD<MatrixDimDynamic> SVDsolverType;

        const std::pair<MatrixDimDynamic,MatrixDimDynamic> NP;
        const SVDsolverType svd;
        
        
        PlanesIntersection(const MatrixDimDynamic& N, const MatrixDimDynamic& P, const double& tol=FLT_EPSILON);
        
        std::pair<bool,VectorDim> snap(const VectorDim& x) const;


    };
    


}
#endif
