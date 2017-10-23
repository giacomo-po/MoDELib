/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanePlaneIntersection_H_
#define model_PlanePlaneIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

namespace model
{
    
    template <int dim>
    struct PlanePlaneIntersection
    {
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT};
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        const SolutionType sol;
        
        /**********************************************************************/
        static SolutionType findIntersection(const VectorDimD& p1,
                                             const VectorDimD& N1,
                                             const VectorDimD& p2,
                                             const VectorDimD& N2)
        {
            
            const double normN1(N1.norm());
            const double normN2(N2.norm());
            assert(normN1>FLT_EPSILON);
            assert(normN2>FLT_EPSILON);
            
            const VectorDimD n1(N1/normN1);
            const VectorDimD n2(N2/normN2);
            const VectorDimD D(n1.cross(n2));
            const double normD(D.norm());
            
            if(normD>FLT_EPSILON)
            {
                // Find the Cartesian point P which minimizes (P-p1)^2+(P-p2)^2 under
                // the constraints (P-P1)*n1=0 and (P-P2)*n2=0
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
                
                assert(fabs((C-p1).dot(n1))<FLT_EPSILON);
                assert(fabs((C-p2).dot(n2))<FLT_EPSILON);
                
                return std::make_tuple(INCIDENT,C,D/normD);
            }
            else
            {
                if(fabs((p1-p2).dot(n1))<FLT_EPSILON)
                {// coincident planes
                    return std::make_tuple(COINCIDENT,0.5*(p1+p2),VectorDimD::Zero());
                }
                else
                {// paralle planes
                    return std::make_tuple(PARALLEL,0.5*(p1+p2),VectorDimD::Zero());
                }
            }
        }
        
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        const IntersectionType& type;
        const VectorDimD& P;
        const VectorDimD& d;
        
        /**********************************************************************/
        PlanePlaneIntersection(const VectorDimD& p1,
                               const VectorDimD& n1,
                               const VectorDimD& p2,
                               const VectorDimD& n2) :
        /* init */ sol(findIntersection(p1,n1,p2,n2)),
        /* init */ type(std::get<0>(sol)),
        /* init */ P(std::get<1>(sol)),
        /* init */ d(std::get<2>(sol))
        {
        }
        
    };
    
}
#endif
