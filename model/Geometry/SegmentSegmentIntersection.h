/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SegmentSegmentIntersection_H_
#define model_SegmentSegmentIntersection_H_

#include <tuple>
#include <map>
#include <Eigen/Dense>

namespace model
{
    
    template <int dim>
    class SegmentSegmentIntersection
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        std::tuple<size_t,VectorDimD,VectorDimD> _tup;
        
        /**********************************************************************/
        static std::tuple<size_t,VectorDimD,VectorDimD> findIntersections(const VectorDimD& A0,
                                                                          const VectorDimD& B0,
                                                                          const VectorDimD& A1,
                                                                          const VectorDimD& B1)
        {/*!\param[in] A0 start point of first segment
          *\param[in] B0   end point of first segment
          *\param[in] A1 start point of second segment
          *\param[in] B1   end point of second segment
          *\returns a tuple where the first element is the number of intersections found (0,1, or two)
          * the second and third elements are the intersection points. Second and
          * third element contain the same values in case of one intersection.
          */
            
            const VectorDimD L0(B0-A0);
            const VectorDimD L1(B1-A1);
            const VectorDimD A0A1(A1-A0);
            
            const double norm2L0(L0.squaredNorm());
            const double norm2L1(L1.squaredNorm());
            
            if(norm2L0<FLT_EPSILON && norm2L1<FLT_EPSILON)
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
            else if(norm2L0<FLT_EPSILON && norm2L1>=FLT_EPSILON)
            {// first segment is degenerate
                const double u1 = -A0A1.dot(L1)/norm2L1;
                const VectorDimD x(A1+u1*L1);
                if((x-A0).squaredNorm()<FLT_EPSILON)
                {// point of min distance is an intersection point
                    if(u1<=-FLT_EPSILON)
                    {// intersection point beyond A1
                        return std::make_tuple(0,A0,A1);
                    }
                    else if(u1>=1.0+FLT_EPSILON)
                    {// intersection point beyond B1
                        return std::make_tuple(0,A0,B1);
                    }
                    else
                    {// intersection point is within [A1,B1]
                        return std::make_tuple(1,0.5*(A0+x),0.5*(A0+x));
                    }
                }
                else
                {// no intersection
                    return std::make_tuple(0,A0,A1);
                }
            }
            else if(norm2L0>=FLT_EPSILON && norm2L1<FLT_EPSILON)
            {// second segment is degenerate
                const double u0 = A0A1.dot(L0)/norm2L0;
                const VectorDimD x(A0+u0*L0);
                if( (x-A1).squaredNorm()<FLT_EPSILON)
                {// point of min distance is an intersection point
                    if(u0<=-FLT_EPSILON)
                    {// intersection point beyond A0
                        return std::make_tuple(0,A0,A1);
                    }
                    else if(u0>=1.0+FLT_EPSILON)
                    {// intersection point beyond B0
                        return std::make_tuple(0,B0,A1);
                    }
                    else
                    {// intersection point is within [A0,B0]
                        return std::make_tuple(1,0.5*(A1+x),0.5*(A1+x));
                    }
                }
                else
                {// no intersection
                    return std::make_tuple(0,A0,A1);
                }
            }
            else
            {// both segments are non-degenerate
                
                const double L0L1(L0.dot(L1));
                const double det(norm2L0*norm2L1-L0L1*L0L1);
                const double num0=norm2L1*L0.dot(A0A1)-L0L1*L1.dot(A0A1);
                const double num1=L0L1*L0.dot(A0A1)-norm2L0*L1.dot(A0A1);
                if(fabs(det)>FLT_EPSILON)
                {
                    const double u0=num0/det;
                    const double u1=num1/det;
                    const VectorDimD x0=A0+u0*L0;
                    const VectorDimD x1=A1+u1*L1;
                    if((x0-x1).squaredNorm()<FLT_EPSILON)
                    {// points of minimum distance coincide, so intersection exists
                        if(   u0>-FLT_EPSILON && u0<1.0+FLT_EPSILON
                           && u1>-FLT_EPSILON && u1<1.0+FLT_EPSILON)
                        {
                                const VectorDimD x=0.5*(x0+x1);
                                return std::make_tuple(1,x,x);
                        }
                        else
                        {// intersection exists but outside segments
                            return std::make_tuple(0,x0,x1);
                        }
                    }
                    else
                    {// no intersections
                        return std::make_tuple(0,A0,A1);
                    }
                    
                }
                else
                {// segments are coincident or parallel
                    if(fabs(num0)<FLT_EPSILON && fabs(num0)<FLT_EPSILON)
                    {// segments are coincident, keep innermost points
                        std::multimap<double,VectorDimD> ms;
                        
                        ms.emplace(0.0,A0);
                        ms.emplace(1.0,B0);
                        const double u1=(A1-A0).dot(L0)/norm2L0;
                        const double u2=(B1-A0).dot(L0)/norm2L0;

                        
//                        THIS IS WRONG THERE COULD BE NON-INTERSECTIONG COINCIDENT LINES! FINISH THE 9 CASES BELOW

                        
                        if(u1<0.0)
                        {
                            if(u2<0.0)
                            {/*         A0----B0
                              * A1---B1
                              */
                                return std::make_tuple(0,A0,A1);
                            }
                            else if(u2>=0.0 && u2<=1.0)
                            {/*         A0----B0
                              *       A1---B1
                              */
                                return std::make_tuple(2,A0,B1);
                            }
                            else // u2>1.0
                            {/*         A0----B0
                              *      A1----------B1
                              */
                                return std::make_tuple(2,A0,B0);
                            }
                        }
                        else if(u1>=0.0 && u1<=1.0)
                        {
                            if(u2<0.0)
                            {/*         A0----B0
                              *      B1----A1
                              */
                                return std::make_tuple(2,A0,A1);
                            }
                            else if(u2>=0.0 && u2<=1.0)
                            {/*      A0----------B0
                              *         A1----B1
                              */
                                return std::make_tuple(2,A1,B1);
                            }
                            else // u2>1.0
                            {/*      A0------B0
                              *           A1-----B1
                              */
                                return std::make_tuple(2,A1,B0);
                            }
                        }
                        else // u1>1.0
                        {
                            if(u2<0.0)
                            {/*      A0------B0
                              *   B1-------------A1
                              */
                                return std::make_tuple(2,A0,B0);
                            }
                            else if(u2>=0.0 && u2<=1.0)
                            {/*      A0------B0
                              *          B1-------A1
                              */
                                return std::make_tuple(2,B1,B0);
                            }
                            else // u2>1.0
                            {/*      A0----B0
                              *               B1----A1
                              */
                                return std::make_tuple(0,B0,B1);
                            }
                        }
                        
//                        ms.emplace(u1,A0+u1*(B0-A0));
//
//                        ms.emplace(u2,A0+u2*(B0-A0));
//                        
//                        auto iter1=ms.begin();
//                        std::advance(iter1,1);
//                        auto iter2=ms.begin();
//                        std::advance(iter2,2);
//                        return std::make_tuple(2,iter1->second,iter2->second);
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
        SegmentSegmentIntersection(const VectorDimD& A0,
                                   const VectorDimD& B0,
                                   const VectorDimD& A1,
                                   const VectorDimD& B1) :
        /* init */ _tup(findIntersections(A0,B0,A1,B1)),
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
