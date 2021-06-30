/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineLineIntersection_CPP_
#define model_LineLineIntersection_CPP_

#include <cfloat>
#include <tuple>
#include <map>
#include <stdexcept>
#include <Eigen/Dense>

#include <LineLineIntersection.h>

namespace model
{
    /**********************************************************************/
    template <int dim>
    typename LineLineIntersection<dim>::SolutionType LineLineIntersection<dim>::findIntersections(const VectorDimD& A0,
                                          const VectorDimD& D0,
                                          const VectorDimD& A1,
                                          const VectorDimD& D1)
    {/*!\param[in] A0 start point of first line
      *\param[in] D0  line direction of first line
      *\param[in] A1 start point of second line
      *\param[in] D1   line direction of second line
      *\returns a tuple, where the first element is an elmenet of IntersectionType (the type of intersection).
      * Second and third elements depend of the type of intersection:
      * - For COINCIDENT lines
      * - the second element is A0
      * - the third elements is A1
      * - For SKEW lines
      * - the second element is the closest point on line 1
      * - the third elements is the closest point on line 2
      * - For INCIDENT lines
      * - the second element is the point of intersection
      * - the third elements is the point of intersection
      * - For PARALLEL lines
      * - the second element is A0
      * - the third elements is A1
      */
        
        const double normD0(D0.norm());
        const double normD1(D1.norm());
        if(normD0<FLT_EPSILON)
        {
            throw std::runtime_error("Norm of line direction D0 is close to singular. norm(D0)="+std::to_string(normD0)+"\n");
        }
        if(normD1<FLT_EPSILON)
        {
            throw std::runtime_error("Norm of line direction D1 is close to singular. norm(D1)="+std::to_string(normD1)+"\n");
        }
        
        const VectorDimD d0(D0/normD0);
        const VectorDimD d1(D1/normD1);
        const VectorDimD A0A1(A1-A0);
        
        if(normD0<FLT_EPSILON && normD1<FLT_EPSILON)
        {// both lines are degenerate
            if(A0A1.squaredNorm()<FLT_EPSILON)
            {
                return std::make_tuple(INCIDENT,0.5*(A0+A1),VectorDimD::Zero());
            }
            else
            {
                return std::make_tuple(SKEW,A0,A1);
            }
        }
        else if(normD0<FLT_EPSILON && normD1>=FLT_EPSILON)
        {// line0 is degenerate
            const double u1 = -A0A1.dot(d1);
            const VectorDimD x(A1+u1*d1);
            if((x-A0).squaredNorm()<FLT_EPSILON)
            {
                return std::make_tuple(INCIDENT,0.5*(A0+x),VectorDimD::Zero());
            }
            else
            {
                return std::make_tuple(SKEW,A0,x);
            }
        }
        else if(normD0>=FLT_EPSILON && normD1<FLT_EPSILON)
        {// line1 is degenerate
            const double u0 = A0A1.dot(d0);
            const VectorDimD x(A0+u0*d0);
            if((x-A1).squaredNorm()<FLT_EPSILON)
            {
                return std::make_tuple(INCIDENT,0.5*(A1+x),VectorDimD::Zero());
            }
            else
            {
                return std::make_tuple(SKEW,x,A1);
            }
        }
        else
        {/* Both lines are non-degenerate. We search for the points x0 and x1
          * that minimize the distance between the two lines:
          * x0=A0+u0*d0
          * x1=A1+u1*d1
          * The scalar parameters u0 and u1 are found as
          * [u0,u1]=argmin{1/2*(A0+u0*d0-A1-u1*d1)^2}
          */
            const double d0d1(d0.dot(d1));
            const double det(1.0-d0d1*d0d1);
            const double num0=d0.dot(A0A1)-d0d1*d1.dot(A0A1);
            const double num1=d0d1*d0.dot(A0A1)-d1.dot(A0A1);
            if(fabs(det)>FLT_EPSILON)
            {
                const double u0=num0/det;
                const double u1=num1/det;
                const VectorDimD x0(A0+u0*d0);
                const VectorDimD x1(A1+u1*d1);
                if((x0-x1).squaredNorm()<FLT_EPSILON)
                {
                    return std::make_tuple(INCIDENT,x0,x1);
                }
                else
                {
                    return std::make_tuple(SKEW,x0,x1);
                }
            }
            else
            {// A solution could not be found. This means parallel or coincident lines
                const double A0A1norm(A0A1.norm());
                if(A0A1norm<FLT_EPSILON)
                {
                    return std::make_tuple(COINCIDENT,A0,A1);
                }
                else if(fabs(fabs(d0.dot(A0A1))-A0A1norm)<FLT_EPSILON)
                {// coi
                    return std::make_tuple(COINCIDENT,A0,A1);
                }
                else
                {
                    return std::make_tuple(PARALLEL,A0,A1);
                }
//                        if(fabs(num0)<FLT_EPSILON && fabs(num1)<FLT_EPSILON)
//                        {
//                            return std::make_tuple(COINCIDENT,0.5*(A0+A1),D0);
//                        }
//                        else
//                        {// parallel segments
//
//                            return std::make_tuple(PARALLEL,A0,A1);
//                        }
            }
            
        }
    }
    
    /**********************************************************************/
    template<int dim>
    LineLineIntersection<dim>::LineLineIntersection(const VectorDimD& A0,
                         const VectorDimD& D0,
                         const VectorDimD& A1,
                         const VectorDimD& D1) :
    /* init */ sol(findIntersections(A0,D0,A1,D1)),
    /* init */ type(std::get<0>(sol)),
    /* init */ x0(std::get<1>(sol)),
    /* init */ x1(std::get<2>(sol))
    {
        
    }

    template class LineLineIntersection<3>;

} /* namespace model */
#endif
