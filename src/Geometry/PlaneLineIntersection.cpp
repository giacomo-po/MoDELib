/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneLineIntersection_CPP_
#define model_PlaneLineIntersection_CPP_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <PlaneLineIntersection.h>

namespace model
{
    /**********************************************************************/
    template <int dim>
    typename PlaneLineIntersection<dim>::SolutionType PlaneLineIntersection<dim>::findIntersection(const VectorDimD& pP,
                                         const VectorDimD& N,
                                         const VectorDimD& pL,
                                         const VectorDimD& D)
    {
        
        const double normN(N.norm());
        assert(normN>FLT_EPSILON && "PLANE MUST HAVE NON-ZERO NORMAL");
        const VectorDimD n(N/normN); // unit plane normal
        
        const double normD(D.norm());
        if(normD<FLT_EPSILON)
        {// line is degenerate (a point)
            if(fabs((pL-pP).dot(n))<FLT_EPSILON)
            {// degenerate incident
                return std::make_tuple(INCIDENT,pL,VectorDimD::Zero());
            }
            else
            {// degenerate parallel (non-intersection)
                return std::make_tuple(PARALLEL,pL-(pL-pP).dot(n)*n,VectorDimD::Zero());
            }
        }
        else
        {
            const VectorDimD d(D/normD); // unit line direction
            const double den(d.dot(n));
            const double num((pP-pL).dot(n));
            
            if(fabs(den)>FLT_EPSILON)
            {// incident
                const double u=num/den;
                return std::make_tuple(INCIDENT,pL+u*d,VectorDimD::Zero());
                
            }
            else
            {// parallel or conincident
                if(fabs(num)<FLT_EPSILON)
                {// coincident
                    return std::make_tuple(COINCIDENT,pL,d);
                    
                }
                else
                {// parallel
                    return std::make_tuple(PARALLEL,pL,d);
                }
            }
            
        }
        
    }
    
    
    /**********************************************************************/
    template <int dim>
    PlaneLineIntersection<dim>::PlaneLineIntersection(const Eigen::Matrix<double,dim,1>& pP,
                          const Eigen::Matrix<double,dim,1>& N,
                          const Eigen::Matrix<double,dim,1>& pL,
                          const Eigen::Matrix<double,dim,1>& D) :
    /* init */ sol(findIntersection(pP,N,pL,D)),
    /* init */ type(std::get<0>(sol)),
    /* init */ P(std::get<1>(sol)),
    /* init */ d(std::get<2>(sol))
    {
        
    }

    template struct PlaneLineIntersection<3>;
        
} /* namespace model */
#endif
