/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CatmullRomSplineSegment_H_
#define model_CatmullRomSplineSegment_H_

#include <math.h>
#include <Eigen/Dense>




#include <ParametricCurve.h>

namespace model
{
    template <int dim>
    class CatmullRomSplineSegment : public ParametricCurve <CatmullRomSplineSegment<dim>, dim>
    {
        typedef Eigen::Matrix<double, dim, 1> VectorDim;
        typedef Eigen::Matrix<double, dim, 4> PointMatrixType;

        const PointMatrixType P;
        const Eigen::Matrix<double, 4, 4> M;
        const PointMatrixType PM;
        

        static Eigen::Matrix<double, 4, 4> get_M(const VectorDim& P0,const VectorDim& P1,const VectorDim& P2,const VectorDim& P3)
        {
            Eigen::Matrix<double, 4, 4> temp;
            const double g0(std::pow((P1-P0).norm(),0.5));
            const double g1(std::pow((P2-P1).norm(),0.5));
            const double g2(std::pow((P3-P2).norm(),0.5));
            
            temp<<0, 1.0, 0, 0,
                -std::pow(g1,2)/g0/(g0+g1), g1/g0-1.0, g0/(g1+g0), 0,
                2.0*std::pow(g1,2)/g0/(g0+g1), -2.0*g1/g0-g1/(g1+g2), 2.0*g1/(g1+g0)+g1/g2, -std::pow(g1,2)/g2/(g1+g2),
                -std::pow(g1,2)/g0/(g0+g1), g1/g0+g1/(g1+g2), -g1/(g0+g1)-g1/g2, std::pow(g1,2)/g2/(g1+g2);
                
            return temp.transpose();
        }

        public:
        CatmullRomSplineSegment(const VectorDim& P0,const VectorDim& P1,const VectorDim& P2,const VectorDim& P3) : 
        /*init*/  P  ((PointMatrixType()<<P0,P1,P2,P3).finished())
        /*init*/, M  (get_M(P0,P1,P2,P3))
        /*init*/, PM (P*M)
        {
            // std::cout<<P<<std::endl;
        }


        VectorDim get_r(const double& u) const
        {
            return PM*(Eigen::Matrix<double, 4,1>()<<1.0,u,std::pow(u,2),std::pow(u,3)).finished();
        }

        VectorDim get_ru(const double& u) const
        {
            return PM*(Eigen::Matrix<double, 4,1>()<<0.0,1.0,2.0*u,3.0*std::pow(u,2)).finished();
        }

        VectorDim get_ruu(const double& u) const
        {
            return PM*(Eigen::Matrix<double, 4,1>()<<0.0,0.0,2.0,6.0*u).finished();
        }
    };
    
}
#endif

