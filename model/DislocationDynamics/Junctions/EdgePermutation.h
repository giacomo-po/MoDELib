/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EdgePermutation_h_
#define model_EdgePermutation_h_


#include <iostream>
#include <deque>
#include <float.h>
#include <Eigen/Dense>

#include <model/Math/CompileTimeMath/Pow.h>

namespace model
{
    
    /**************************************************************************/
    template <int n>
    struct EdgePermutation
    {
        static_assert(n>0,"n MUST BE >0");
        static const Eigen::Matrix<double,Pow<3,n>::value,n> edgePerm;
    };
    
    // Static data
    template <int n>
    const Eigen::Matrix<double,Pow<3,n>::value,n> EdgePermutation<n>::edgePerm=(Eigen::Matrix<double,n,Pow<3,n>::value>()<<
                                                                                +1.0*Eigen::Matrix<double,1,Pow<3,n-1>::value>::Ones(),
                                                                                -1.0*Eigen::Matrix<double,1,Pow<3,n-1>::value>::Ones(),
                                                                                +0.0*Eigen::Matrix<double,1,Pow<3,n-1>::value>::Ones(),
                                                                                  EdgePermutation<n-1>::edgePerm.transpose(),
                                                                                  EdgePermutation<n-1>::edgePerm.transpose(),
                                                                                  EdgePermutation<n-1>::edgePerm.transpose()
                                                                                ).finished().transpose();
    
    /**************************************************************************/
    template <>
    struct EdgePermutation<1>
    {
        static const Eigen::Matrix<double,3,1> edgePerm;
    };
    
    // Static data
    const Eigen::Matrix<double,3,1> EdgePermutation<1>::edgePerm=(Eigen::Matrix<double,3,1>()<<1.0,-1.0,0.0).finished();
    
    
    /**************************************************************************/
    struct EdgePermutations
    {
        
        /**********************************************************************/
        static Eigen::MatrixXd edgePerm(const int& n)
        {
            switch (n)
            {
                case 1:
                    return EdgePermutation< 1>::edgePerm;
                    break;
                case 2:
                    return EdgePermutation< 2>::edgePerm;
                    break;
                case 3:
                    return EdgePermutation< 3>::edgePerm;
                    break;
                case 4:
                    return EdgePermutation< 4>::edgePerm;
                    break;
                case 5:
                    return EdgePermutation< 5>::edgePerm;
                    break;
                case 6:
                    return EdgePermutation< 6>::edgePerm;
                    break;
                case 7:
                    return EdgePermutation< 7>::edgePerm;
                    break;
//                case 8:
//                    return EdgePermutation< 8>::edgePerm;
//                    break;
//                case 9:
//                    return EdgePermutation< 9>::edgePerm;
//                    break;
//                case 10:
//                    return EdgePermutation<10>::edgePerm;
//                    break;
                    
                default:
//                    assert(0 && "CASE NOT IMPLEMENTED");
                    model::cout<<"WARNING: EdgePermutations, case not implemented"<<std::endl;
                    return Eigen::MatrixXd::Zero(1,n);
                    
            }
            
        }
        
        
        /**********************************************************************/
        static std::deque<Eigen::Matrix<double,1,Eigen::Dynamic>> edgeStats(const Eigen::MatrixXd& m)
        {
            std::deque<Eigen::Matrix<double,1,Eigen::Dynamic>> vec;
            const Eigen::MatrixXd perm(edgePerm(m.rows()));
            const Eigen::MatrixXd temp=(perm*m).rowwise().squaredNorm();
            
            for (int k=0;k<temp.rows();++k)
            {
                if(temp(k)<FLT_EPSILON && !perm.row(k).all() && perm.row(k).squaredNorm()>FLT_EPSILON)
                {
                    bool isNewVector=true;
                    for (auto v : vec)
                    {
                        isNewVector *=  (perm.row(k)-v).squaredNorm()>FLT_EPSILON &&
                                        (perm.row(k)+v).squaredNorm()>FLT_EPSILON;
                    }
                    if(isNewVector)
                    {
                        vec.emplace_back(perm.row(k));

                    }
                }
            }
            
            return vec;
        }
        
    };
    
}
#endif
