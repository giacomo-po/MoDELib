/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineLineIntersection_H_
#define model_LineLineIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>

namespace model
{
    
    template <int dim>
    class LineLineIntersection
    {
        
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT,SKEW};
        
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        
        const SolutionType sol;
        
        /**********************************************************************/
        static SolutionType findntersections(const VectorDimD& A0,
                                             const VectorDimD& D0,
                                             const VectorDimD& A1,
                                             const VectorDimD& D1)
        {/*!\param[in] A0 start point of first segment
          *\param[in] D0   line direction
          *\param[in] A1 start point of second segment
          *\param[in] D1   end point of second segment
          *\returns a tuple, where the first element is the type of intersection.
          * Second and third elements depend of the type of intersection:
          * - For COINCIDENT lines
          * - the second element is an intersection point
          * - the third elements is the direction of the line of intersection
          * - For SKEW lines
          * - the second element is the closest point on line 1
          * - the third elements is the closest point on line 2
          * - For INCIDENT lines
          * - the second element is the point of intersection
          * - the third elements is set to zero (not used)
          * - For PARALLEL lines
          * - the second element is the origin of line 1
          * - the third elements is the origin of line 2
          */
            
            const double normD0(D0.norm());
            const VectorDimD d0(D0/normD0);
            const double normD1(D1.norm());
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
            {// both line and segment are non-degenerate
                
                const double d0d1(d0.dot(d1));
                const double det(1.0-d0d1*d0d1);
                const double num0=d0.dot(A0A1)-d0d1*d1.dot(A0A1);
                const double num1=d0d1*d0.dot(A0A1)-d1.dot(A0A1);
                if(fabs(det)>FLT_EPSILON)
                {
                    const double u0=num0/det;
                    const double u1=num1/det;
                    const VectorDimD x0=A0+u0*d0;
                    const VectorDimD x1=A1+u1*d1;
                    if((x0-x1).squaredNorm()<FLT_EPSILON)
                    {
                        return std::make_tuple(INCIDENT,0.5*(x0+x1),VectorDimD::Zero());
                    }
                    else
                    {
                        return std::make_tuple(SKEW,x0,x1);
                    }
                }
                else
                {
                    if(fabs(num0)<FLT_EPSILON && fabs(num0)<FLT_EPSILON)
                    {
                        return std::make_tuple(COINCIDENT,0.5*(A0+A1),D0);
                    }
                    else
                    {// parallel segments
                        
                        return std::make_tuple(PARALLEL,A0,A1);
                    }
                }
                
            }
        }
        
    public:
        
        const IntersectionType& type;
        const VectorDimD& x0;
        const VectorDimD& x1;
        
        /**********************************************************************/
        LineLineIntersection(const VectorDimD& A0,
                             const VectorDimD& D0,
                             const VectorDimD& A1,
                             const VectorDimD& D1) :
        /* init */ sol(findIntersections(A0,D0,A1,D1)),
        /* init */ type(std::get<0>(sol)),
        /* init */ x0(std::get<1>(sol)),
        /* init */ x1(std::get<2>(sol))
        {
            
        }
        
    };
    
    /******************************************************************/
    /******************************************************************/
} /* namespace model */
#endif
