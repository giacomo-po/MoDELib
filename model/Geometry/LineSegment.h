/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineSegment_H_
#define model_LineSegment_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>

#include <model/Geometry/SegmentSegmentDistance.h>

namespace model
{
    
    template <int dim>
    struct LineSegment
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        VectorDim P0;
        VectorDim P1;
        
        /**********************************************************************/
        LineSegment(const VectorDim& p0,
                    const VectorDim& p1) :
        /* init */ P0(p0),
        /* init */ P1(p1)
        {
            
        }
        
        /**********************************************************************/
        VectorDim snap(const VectorDim& x) const
        {/*!\param[in] x point to be snapped to the LineSegment
          * \returns the closest point to x on the LineSegment
          *
          * miminime distance to the line:
          * u=argMin{ (P0+u*(P1-P0)-x)^2 }
          * (P0+u*(P1-P0)-x)*(P1-P0)=0
          * u=(x-P0)*(P1-P0)/(P1-P0)^2
          */
            
            const VectorDim chord(P1-P0);
            const double chordNorm2(chord.squaredNorm());
            if(chordNorm2>FLT_EPSILON)
            {
                const double u=(x-P0).dot(chord)/chordNorm2;
                if(u<0.0)
                {
                    return P0;
                }
                else if(u>1.0)
                {
                    return P1;
                }
                else
                {
                    return P0+u*chord;
                }
            }
            else
            {
                return 0.5*(P0+P1);
            }
            
        }
        
        
        /**********************************************************************/
        bool contains(const VectorDim& x) const
        {
            return (x-snap(x)).squaredNorm()<FLT_EPSILON;
        }
        
        /**********************************************************************/
        SegmentSegmentDistance<dim> distanceTo(const LineSegment<dim>& other,
                                               const double& degeneracyValue=0.0) const
        {
            return SegmentSegmentDistance<dim>(P0,P1,other.P0,other.P1,degeneracyValue);
        }
        
    };
    
}
#endif
