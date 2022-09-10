/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanesIntersection_H_
#define model_PlanesIntersection_H_

//#include <Eigen/LU>
#include <Eigen/SVD>
#include <Plane.h>

namespace model
{
    
    template <int dim>
    struct PlanesIntersection
    {
        
        typedef Plane<dim>::VectorDim VectorDim;
        typedef Eigen::JacobiSVD<MatrixXd, ComputeThinU | ComputeThinV> SVDsolverType;
//        typedef Eigen::BDCSVD   <MatrixXd, ComputeThinU | ComputeThinV> SVDsolverType;

        
        
        static Eigen::MatrixXd getN(const std:vector<const Plane<dim>* const>& planes)
        {
            Eigen::MatrixXd N(Eigen::MatrixXd::Zero(dim,planes.size()));
            for(int c=0;c<planes.size();++c)
            {
                N.col(c)=planes[c].unitNormal;
            }
            return N;
        }
        
        const Eigen::MatrixXd N;
        const SVDsolverType svd;

//        const size_t rank; // note: rank(N)=rank(N^T*N)
        
        PlanesIntersection(const std:vector<const Plane<dim>* const>& planes, const double& tol=FLT_EPSILON) :
        /* init */ N(getN(planes))
        /* init */ svd(N)
//        /* init */ rank(N.fullPivLu().setThreshold(tol).rank())
        {
            
        }
        
        VectorDim snap(const VectorDim& x) const
        {
            const size_t N(planes.size());
            Eigen::MatrixXd A(Eigen::MatrixXd::Identity(N,N));
            for(size_t i=0;i<N;++i)
            {
                for(size_t j=i+1;j<N;++j)
                {
                    A(i,j)=planes[i].unitNormaldot(planes[j]);
                    A(j,i)=A(i,j);
                }
            }
        }


    };
    
}
#endif
