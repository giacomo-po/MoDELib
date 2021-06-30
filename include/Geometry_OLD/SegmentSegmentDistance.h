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
        typedef std::deque<IntersectionPointType> IntersectionContainerType;
        
        static constexpr double tol=FLT_EPSILON;
        
        /**********************************************************************/
        static double limitRange(const double& v);
        
        /**********************************************************************/
        std::pair<double,double> step4U(const double& T) const;
        
        /**********************************************************************/
        std::pair<double,double> step4(const double& U) const;
        
        /**********************************************************************/
        std::pair<double,double> step3(const double& T) const;
        
        /**********************************************************************/
        std::pair<double,double> step2() const;
        
        /**********************************************************************/
        std::pair<double,double> getTU() const;
        
        
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
        
        
        /**********************************************************************/
        SegmentSegmentDistance(const VectorDim& A0,
                               const VectorDim& B0,
                               const VectorDim& C0,
                               const VectorDim& D0,
                               const double& degeneracyValue_in=0.0);
        
        /**********************************************************************/
        IntersectionContainerType intersectionSegment() const;
        
        
    };
    
}
#endif
