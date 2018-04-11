/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SegmentSegmentDistance_H_
#define model_SegmentSegmentDistance_H_

#include <cfloat>
#include <deque>
#include <map>
#include <tuple>
#include <Eigen/Dense>

namespace model
{
    /*! Class template which computes the distance between two lines segments
     * in dimension dim. Algorithm adapted from
     * V. LUMELSKY. ON FAST COMPUTATION OF DISTANCE BETWEEN LINE SEGMENTS.
     * Information Processing Letters 21 (1985) 55-61.
     */
    template <int dim>
    struct SegmentSegmentDistance
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        typedef std::tuple<VectorDim,double,double> IntersectionPointType;
        typedef std::deque<IntersectionPointType,Eigen::aligned_allocator<IntersectionPointType>> IntersectionContainerType;
        
        static constexpr double tol=FLT_EPSILON;
        
        /**********************************************************************/
        static double limitRange(const double& v)
        {
            if(v<0.0)
            {
                return 0.0;
            }
            else if(v>1.0)
            {
                return 1.0;
            }
            else
            {
                return v;
            }
        }
        
        /**********************************************************************/
        std::pair<double,double> step4U(const double& T) const
        {/*! Compute t from eq. (10) and limit by eq. (12)
          */
            return std::make_pair(T,limitRange((T*R-S2)/D2));
        }
        
        /**********************************************************************/
        std::pair<double,double> step4(const double& U) const
        {/*! Compute t from eq. (10) and limit by eq. (12)
          */
            return std::make_pair(limitRange((U*R+S1)/D1),U);
        }
        
        /**********************************************************************/
        std::pair<double,double> step3(const double& T) const
        {
            const double U=(T*R-S2)/D2;       // compute U from eq. (10)
            if(U<0.0 || U>1.0)
            {// u is not in [0,1]
                return step4(limitRange(U));
            }
            else
            {// u is in [0,1]
                return std::make_pair(T,U);
            }
        }
        
        /**********************************************************************/
        std::pair<double,double> step2() const
        {
            const double T=limitRange((S1*D2-S2*R)/den); // compute t from eq. (11) and limit by eq. (12)
            return step3(T);
        }
        
        /**********************************************************************/
        std::pair<double,double> getTU() const
        {
            if(D1<tol && D2<tol)
            {// Step 1b: both segments are degenerate
                return std::make_pair(degeneracyValue,degeneracyValue);
            }
            else if(D1<tol && D2>=tol)
            {// Step 1a: first segment is degenerate
                return step4U(degeneracyValue);
            }
            else if(D1>=tol && D2<tol)
            {// Step 1a: second segment is degenerate
                return step4(degeneracyValue);
            }
            else
            {// both segments are not degenerate
                if(fabs(den)>tol) // MAY NEED fabs(den)>tol*D1*D2
                {// Step 1d: skew segments
                    return step2();
                }
                else
                {//  Step 1c: parallel or coincident segments
                    return step3(degeneracyValue);
                }
            }
        }
        
        
        const VectorDim A;
        const VectorDim B;
        const VectorDim C;
        const VectorDim D;
        const double degeneracyValue;
        
        const VectorDim d1;
        const VectorDim d2;
        const VectorDim d12;
        
        const double R;
        const double S1;
        const double S2;
        const double D1;
        const double D2;
        const double den;
        
        const std::pair<double,double> tu;
        const double t;
        const double u;
        const double dMin;
        const VectorDim x0;
        const VectorDim x1;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /**********************************************************************/
        SegmentSegmentDistance(const VectorDim& A0,
                               const VectorDim& B0,
                               const VectorDim& C0,
                               const VectorDim& D0,
                               const double& degeneracyValue_in=0.0) :
        /* init */ A(A0),
        /* init */ B(B0),
        /* init */ C(C0),
        /* init */ D(D0),
        /* init */ degeneracyValue(degeneracyValue_in),
        /* init */ d1(B-A),
        /* init */ d2(D-C),
        /* init */ d12(C-A),
        /* init */ R(d1.dot(d2)),
        /* init */ S1(d1.dot(d12)),
        /* init */ S2(d2.dot(d12)),
        /* init */ D1(d1.squaredNorm()), // length of AB squared
        /* init */ D2(d2.squaredNorm()), // length of CD squared
        /* init */ den(D1*D2-R*R),
        /* init */ tu(getTU()),
        /* init */ t(tu.first),
        /* init */ u(tu.second),
        /* init */ dMin((d1*t-d2*u-d12).norm()),
        /* init */ x0(A+t*d1),
        /* init */ x1(C+u*d2)
        {/*!\param[in] A start point of first segment
          * \param[in] B   end point of first segment
          * \param[in] C start point of second segment
          * \param[in] D   end point of second segment
          */
            
        }
        
        /**********************************************************************/
        IntersectionContainerType intersectionSegment() const
        {/*\returns the point(s) of intersection between the segments, if any.
          * The return type is deque<tuple<x,t,u>>, where x is the intersection
          * point, t the parameter on the first segment, u the parameter on the
          * second segment. The deque is empty if no intersection was found.
          * If one intersection is found, the deque contains that intersection.
          * In case of coincident segments, the deque contains two points. These
          * points define the line element over which the original segments overlap.
          */
            
            IntersectionContainerType temp;
            
            if(dMin<FLT_EPSILON)
            {// an intersection exists
                
                if(D1<tol || D2<tol)
                {// At least one segment is degenerate
                    temp.emplace_back(0.5*(x0+x1),t,u);
                }
                else
                {// both segments are not degenerate
                    
//                    std::cout<<"den="<<den<<std::endl;
//                    std::cout<<"tol="<<tol<<std::endl;
                    
                    if(fabs(den)>tol*D1*D2)
                    {// Step 1d: skew segments
                        temp.emplace_back(0.5*(x0+x1),t,u);
                    }
                    else
                    {//  coincident segments. Since there must be an intersection, the two innermost points are the overlapping segment

                        const double tC=(C-A).dot(d1)/D1; // value of t for point C
                        const double tD=(D-A).dot(d1)/D1; // value of t for point D

                        std::multimap<double,VectorDim> ms;
                        ms.emplace(0.0,A);
                        ms.emplace(1.0,B);
                        ms.emplace(tC,C);
                        ms.emplace(tD,D);
                        
                        auto iter1=ms.begin();
                        std::advance(iter1,1);
                        auto iter2=ms.begin();
                        std::advance(iter2,2);

                        temp.emplace_back(iter1->second,iter1->first,(iter1->second-C).dot(d2)/D2);
                        temp.emplace_back(iter2->second,iter2->first,(iter2->second-C).dot(d2)/D2);
                        
                    }
                }
                
            }
            
            return temp;
        }
        
        
    };
    
}
#endif
