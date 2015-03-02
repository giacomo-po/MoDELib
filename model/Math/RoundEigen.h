/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_roundEigen_H_
#define model_roundEigen_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <Eigen/Dense>

namespace model
{
    
    /**************************************************************************/
    //TO DO: chenge floorEigen into a class like RoundEigen
    template <short unsigned int dim>
    Eigen::Matrix<int,dim,1> floorEigen(const Eigen::Matrix<double,dim,1>& P)
    {
        return (Eigen::Matrix<int,dim,1>()<< (int)std::floor(P(0)), floorEigen<dim-1>(P.template segment<dim-1>(1))).finished();
    }
    
    template <>
    Eigen::Matrix<int,1,1> floorEigen<1>(const Eigen::Matrix<double,1,1>& P)
    {
        return (Eigen::Matrix<int,1,1>()<< (int)std::floor(P(0)) ).finished();
    }
    
    /**************************************************************************/
    template<typename T,short unsigned int dim>
    struct RoundEigen
    {
        static Eigen::Matrix<T,dim,1> round(const Eigen::Matrix<double,dim,1>& P)
        {
            return (Eigen::Matrix<T,dim,1>()<< (T)std::round(P(0)), RoundEigen<T,dim-1>::round(P.template segment<dim-1>(1))).finished();
        }
    };
    
    /**************************************************************************/
    template<typename T>
    struct RoundEigen<T,1>
    {
        static Eigen::Matrix<T,1,1> round(const Eigen::Matrix<double,1,1>& P)
        {
            return (Eigen::Matrix<T,1,1>()<< (T)std::round(P(0)) ).finished();
        }
    };
    
} // close namespace model
#endif