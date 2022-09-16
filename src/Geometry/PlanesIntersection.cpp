/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanesIntersection_cpp_
#define model_PlanesIntersection_cpp_

#include <iostream>
#include <PlanesIntersection.h>

//#include <Plane.h>

namespace model
{
    
    
        
        template<int dim>
        PlanesIntersection<dim>::PlanesIntersection(const MatrixDimDynamic& N, const MatrixDimDynamic& P, const double& tol) :
        /* init */ NP(std::make_pair(N,P))
        /* init */,svd(SVDsolverType().setThreshold(tol).compute(NP.first,Eigen::ComputeThinU | Eigen::ComputeThinV))
        {// Ax=b
        }
        
        template<int dim>
        std::pair<bool,typename PlanesIntersection<dim>::VectorDim> PlanesIntersection<dim>::snap(const VectorDim& x) const
        {/*!\param [in] x point to be snapped to the intersection between the planes
          * \returns a pair where pair.first is true if snapping was successful, and pair.second is the position of the snapped point
          *
          * The snapped point y minimizes the distance from x constrained to y belonging to the planes.
          * Using lgrange multipliers, the solution is y=x-N*lam, where lam is a vector of Lagrance multipliers
          * solving N^T*N *lam = b, and b_k = n_k*(x-P_k), where n_k and P_k are the normal and origin of the k0th plane.
          * We solve this using svd, hence we get the "best solution", which may not actually satisfy the contraints if the planes are not intersecting.
          * Hence we check that the solution actually belongs to the intersectoin between the planes before returning it.
          */

            Eigen::Matrix<double,Eigen::Dynamic,1> b(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(NP.first.cols(),1));
            Eigen::Matrix<double,Eigen::Dynamic,1> b1(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(NP.first.cols(),1));
            for(int k=0;k<NP.first.cols();++k)
            {
                b(k)=NP.first.col(k).dot(x-NP.second.col(k));
                b1(k)=NP.first.col(k).dot(NP.second.col(k));
            }
            
            const int r(svd.rank());
            const Eigen::MatrixXd Vr(svd.matrixV().block(0,0,NP.first.cols(),r));
            const Eigen::VectorXd Sinv2(svd.singularValues().segment(0,r).array().pow(-2).matrix());
            const Eigen::VectorXd lam(Vr*Sinv2.asDiagonal()*Vr.transpose()*b);
            const VectorDim y(x-NP.first*lam);
                        
            const double b1norm2(b1.squaredNorm());
            if(b1norm2>svd.threshold())
            {
                if((NP.first.transpose()*y-b1).squaredNorm()/b1norm2<svd.threshold())
                {// y belongs to planes
                    return std::make_pair(true,y);
                }
                else
                {
                    return std::make_pair(false,y);
                }
            }
            else
            {
                if((NP.first.transpose()*y-b1).squaredNorm()<svd.threshold())
                {// y belongs to planes
                    return std::make_pair(true,y);
                }
                else
                {
                    return std::make_pair(false,y);
                }
            }
        }
    
template struct PlanesIntersection<1>;
template struct PlanesIntersection<2>;
template struct PlanesIntersection<3>;

}
#endif
