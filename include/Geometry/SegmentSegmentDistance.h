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

        typedef Eigen::Matrix<long double,dim,1> VectorDimLD;
        
        typedef std::tuple<VectorDim,double,double> IntersectionPointType;
        typedef std::deque<IntersectionPointType> IntersectionContainerType;
        
        static constexpr double tol=100*DBL_EPSILON;
        
        /**********************************************************************/
        static long double limitRange(const long double& v);
        
        
        /**********************************************************************/
        std::pair<long double,long double> step4U(const long double& T) const;
        
        
        /**********************************************************************/
        std::pair<long double,long double> step4(const long double& U) const;
       
        
        /**********************************************************************/
        std::pair<long double,long double> step3(const long double& T) const;
        
        
        /**********************************************************************/
        std::pair<long double,long double> step2() const;
        
        
        /**********************************************************************/
        std::pair<long double,long double> getTU() const;
       
        
        
        const VectorDimLD A;
        const VectorDimLD B;
        const VectorDimLD C;
        const VectorDimLD D;
        const double degeneracyValue;
        
        const VectorDimLD d1;
        const VectorDimLD d2;
        const VectorDimLD d12;
        
        const long double R;
        const long double S1;
        const long double S2;
        const long double D1;
        const long double D2;
        const long double den;
        
        const std::pair<long double,long double> tu;
        const long double t;
        const long double u;
        const long double dMin;
        const VectorDim x0;
        const VectorDim x1;
        
        
        /**********************************************************************/
        SegmentSegmentDistance(const VectorDim& A0,
                               const VectorDim& B0,
                               const VectorDim& C0,
                               const VectorDim& D0,
                               const double& degeneracyValue_in=0.0) ;
        
        /**********************************************************************/
        IntersectionContainerType intersectionSegment() const;

        bool isIntersecting(const double tolerance = FLT_EPSILON) const;
    };
}
#endif