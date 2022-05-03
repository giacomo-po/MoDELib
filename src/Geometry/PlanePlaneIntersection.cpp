/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanePlaneIntersection_CPP_
#define model_PlanePlaneIntersection_CPP_

#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Plane.h>

#include <iostream>
#include<iomanip>
#include <PlanePlaneIntersection.h>

namespace model
{
    
    /**********************************************************************/
    template <int dim>
    bool PlanePlaneIntersection<dim>::checkIntersection(const VectorDimD& C,
                                  const VectorDimD& p1,
                                  const VectorDimD& n1,
                                  const VectorDimD& p2,
                                  const VectorDimD& n2,
                                  const double& tol)
    {
    
        // Check that segments C-p1 and C-p2 are orthogonal to n1 and n2, respectively
        const double Cp1((C-p1).norm());
        const double Cp2((C-p2).norm());
        const bool Cp1OK((Cp1<tol)? true : fabs((C-p1).dot(n1))<tol*Cp1);
        const bool Cp2OK((Cp2<tol)? true : fabs((C-p2).dot(n2))<tol*Cp2);
        const bool success(Cp1OK && Cp2OK);
        if(!success)
        {
            std::cout<<"PlanePlaneIntersection FAILED"<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"C="<<C.transpose()<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"p1="<<p1.transpose()<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"n1="<<n1.transpose()<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"p2="<<p2.transpose()<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"n2="<<n2.transpose()<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"Cp1="<<Cp1<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"Cp2="<<Cp2<<std::endl;
            std::cout<<std::setprecision(15)<<std::scientific<<"tol="<<tol<<std::endl;
        }
        
        return success;
    }
    
    /**********************************************************************/
    template <int dim>
    typename PlanePlaneIntersection<dim>::SolutionType PlanePlaneIntersection<dim>::findIntersection(const VectorDimD& p1,
                                         const VectorDimD& N1,
                                         const VectorDimD& p2,
                                         const VectorDimD& N2,
                                         const double& tol)
    {
        
        const double normN1(N1.norm());
        const double normN2(N2.norm());
        assert(normN1>tol);
        assert(normN2>tol);
        
        const VectorDimD n1(N1/normN1);
        const VectorDimD n2(N2/normN2);
        const VectorDimD D(n1.cross(n2));
        const double normD(D.norm());
          
        if(normD>tol)
        {
            // Find the Cartesian point P which minimizes (C-p1)^2+(C-p2)^2 under
            // the constraints (C-p1)*n1=0 and (C-p2)*n2=0
            Eigen::Matrix<double,dim+2,dim+2> M(Eigen::Matrix<double,dim+2,dim+2>::Zero());
            M.template block<dim,dim>(0,0)=2.0*Eigen::Matrix<double,dim,dim>::Identity();
            M.template block<dim,1>(0,dim+0)=n1;
            M.template block<dim,1>(0,dim+1)=n2;
            M.template block<1,dim>(dim+0,0)=n1.transpose();
            M.template block<1,dim>(dim+1,0)=n2.transpose();
            
            Eigen::Matrix<double,dim+2,1> b(Eigen::Matrix<double,dim+2,1>::Zero());
            b.template segment<dim>(0)=p1+p2;
            b(dim+0)=p1.dot(n1);
            b(dim+1)=p2.dot(n2);
            
            const VectorDimD C=M.ldlt().solve(b).template segment<dim>(0);
            assert(checkIntersection(C,p1,n1,p2,n2,tol) && "PlanePlaneIntersection FAILED.");

            return std::make_tuple(INCIDENT,C,D/normD);
        }
        else
        {
            if(fabs((p1-p2).dot(n1))<tol)
            {// coincident planes
                return std::make_tuple(COINCIDENT,0.5*(p1+p2),VectorDimD::Zero());
            }
            else
            {// paralle planes
                return std::make_tuple(PARALLEL,0.5*(p1+p2),VectorDimD::Zero());
            }
        }
    }

    /**********************************************************************/
    template <int dim>
    PlanePlaneIntersection<dim>::PlanePlaneIntersection(const VectorDimD& p1,
                          const VectorDimD& n1,
                          const VectorDimD& p2,
                          const VectorDimD& n2,
                          const double& tolerance) :
    /* init */ sol(findIntersection(p1,n1,p2,n2,tolerance)),
    /* init */ type(std::get<0>(sol)),
    /* init */ P(std::get<1>(sol)),
    /* init */ d(std::get<2>(sol))
    {
    }

    /**********************************************************************/
    template <int dim>
    PlanePlaneIntersection<dim>::PlanePlaneIntersection(const Plane<dim>& plane0,
                          const Plane<dim>& plane1,
                          const double& tolerance) :
    /* init */ sol(findIntersection(plane0.P,plane0.unitNormal,plane1.P,plane1.unitNormal,tolerance)),
    /* init */ type(std::get<0>(sol)),
    /* init */ P(std::get<1>(sol)),
    /* init */ d(std::get<2>(sol))
    {
    }
    
    template struct PlanePlaneIntersection<3>;
    
}
#endif
