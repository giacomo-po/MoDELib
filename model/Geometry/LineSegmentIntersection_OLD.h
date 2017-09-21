/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineSegmentIntersection_H_
#define model_LineSegmentIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>

namespace model
{
    
    template <int dim>
    class LineSegmentIntersection
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        std::tuple<size_t,VectorDimD,VectorDimD> _tup;
        
        /**********************************************************************/
        static std::tuple<size_t,VectorDimD,VectorDimD> oneSidedIntersections(const VectorDimD& A0,
                                                                          const VectorDimD& D,
                                                                          const VectorDimD& A1,
                                                                          const VectorDimD& B1)
        {/*!\param[in] A0 start point of first segment
          *\param[in] d   line direction
          *\param[in] A1 start point of second segment
          *\param[in] B1   end point of second segment
          *\returns a tuple where the first element is the number of intersections found (0,1, or two)
          * the second and third elements are the intersection points. Second and
          * third element contain the same values in case of one intersection.
          */

                        const double normD(D.norm());
            const VectorDimD d(D/normD);
            const VectorDimD L1(B1-A1);
            const VectorDimD A0A1(A1-A0);
            

            const double norm2L1(L1.squaredNorm());
            
            if(normD<FLT_EPSILON && norm2L1<FLT_EPSILON)
            {// both segments are degenerate
                if(A0A1.squaredNorm()<FLT_EPSILON)
                {
                    return std::make_tuple(1,A0,A0);
                }
                else
                {
                    return std::make_tuple(0,A0,A1);
                }
            }
            else if(normD<FLT_EPSILON && norm2L1>=FLT_EPSILON)
            {// line is degenerate
                const double u1 = -A0A1.dot(L1)/norm2L1;
                if(u1>-FLT_EPSILON && u1<1.0+FLT_EPSILON)
                {
                    return std::make_tuple(1,A0,A0);
                }
                else
                {
                    return std::make_tuple(0,A0,A1);
                }
            }
            else if(normD>=FLT_EPSILON && norm2L1<FLT_EPSILON)
            {// segment is degenerate
                const double u0 = A0A1.dot(d);
                if(u0>-FLT_EPSILON && u0<1.0+FLT_EPSILON)
                {
                    return std::make_tuple(1,A1,A1);
                }
                else
                {
                    return std::make_tuple(0,A0,A1);
                }
            }
            else
            {// both line and segment are non-degenerate
                
                const double dL1(d.dot(L1));
                const double det(norm2L1-dL1*dL1);
                const double num0=norm2L1*d.dot(A0A1)-dL1*L1.dot(A0A1);
                const double num1=dL1*d.dot(A0A1)-L1.dot(A0A1);
                if(fabs(det)>FLT_EPSILON)
                {
                    const double u0=num0/det;
                    const double u1=num1/det;
                    if(   u0>-FLT_EPSILON
                       && u1>-FLT_EPSILON && u1<1.0+FLT_EPSILON)
                    {
                        const VectorDimD x0=A0+u0*d;
                        const VectorDimD x1=A1+u1*L1;
                        if((x0-x1).squaredNorm()<FLT_EPSILON)
                        {
                            const VectorDimD x=0.5*(x0+x1);
                            return std::make_tuple(1,x,x);
                        }
                        else
                        {
                            return std::make_tuple(0,x0,x1);
                        }
                    }
                    else
                    {
                        return std::make_tuple(0,A0,A1);
                    }
                }
                else
                {
                    if(fabs(num0)<FLT_EPSILON && fabs(num0)<FLT_EPSILON)
                    {
                        std::multimap<double,VectorDimD> ms;
                        
                        ms.emplace(0.0,A0);
                        const double u1=(A1-A0).dot(d);
                        ms.emplace(u1,A0+u1*d);
                        const double u2=(B1-A0).dot(d);
                        ms.emplace(u2,A0+u2*d);
                        
                        auto iter1=ms.begin();
                        std::advance(iter1,1);
                        auto iter2=ms.begin();
                        std::advance(iter2,2);
                        
                        return std::make_tuple(2,iter1->second,iter2->second);
                    }
                    else
                    {// parallel segments
                        
                        return std::make_tuple(0,A0,A1);
                    }
                }
                
            }
        }
        
    public:
        
        const size_t& size;
        const VectorDimD& x0;
        const VectorDimD& x1;
        
        /**********************************************************************/
        LineSegmentIntersection(const VectorDimD& A0,
                                   const VectorDimD& D,
                                   const VectorDimD& A1,
                                   const VectorDimD& B1) :
        /* init */ _tup(findIntersections(A0,D,A1,B1)),
        /* init */ size(std::get<0>(_tup)),
        /* init */ x0(std::get<1>(_tup)),
        /* init */ x1(std::get<2>(_tup))
        {
            
        }
        
    };
    
    /******************************************************************/
    /******************************************************************/
} /* namespace model */
#endif
