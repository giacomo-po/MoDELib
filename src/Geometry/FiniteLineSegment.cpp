/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FiniteLineSegment_CPP_
#define model_FiniteLineSegment_CPP_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>

#include <SegmentSegmentDistance.h>
#include <FiniteLineSegment.h>

namespace model
{
    
    /**********************************************************************/
    template <int dim>
    FiniteLineSegment<dim>::FiniteLineSegment(const VectorDim& p0,
                const VectorDim& p1) :
    /* init */ P0(p0),
    /* init */ P1(p1)
    {
        
    }
    
    /**********************************************************************/
    template <int dim>
    typename FiniteLineSegment<dim>::VectorDim FiniteLineSegment<dim>::snap(const VectorDim& x) const
    {/*!\param[in] x point to be snapped to the FiniteLineSegment
      * \returns the closest point to x on the FiniteLineSegment
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
    template <int dim>
    typename FiniteLineSegment<dim>::VectorDim FiniteLineSegment<dim>::snapToInfiniteLine(const VectorDim& x) const
    {/*!\param[in] x point to be snapped to the FiniteLineSegment
      * \returns the closest point to x on the FiniteLineSegment
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
            return P0 + u * chord;
        }
        else
        {
            assert(false && "Not possible to snap for an infinite line");
            return 0.5*(P0+P1);
        }
        
    }

    /**********************************************************************/
    template <int dim>
    SegmentSegmentDistance<dim> FiniteLineSegment<dim>::distanceTo(const FiniteLineSegment<dim>& other,
                                           const double& degeneracyValue) const
    {
        return SegmentSegmentDistance<dim>(P0,P1,other.P0,other.P1,degeneracyValue);
    }
    
    
    /**********************************************************************/
    template <int dim>
    double FiniteLineSegment<dim>::distanceTo(const VectorDim& x) const
    {
        return (x-snap(x)).norm();
    }
    
    /**********************************************************************/
    template <int dim>
    bool FiniteLineSegment<dim>::contains(const VectorDim& x) const
    {
        return distanceTo(x)<FLT_EPSILON;
    }

    /**********************************************************************/
    template <int dim>
    double FiniteLineSegment<dim>::length() const
    {
        return (P1-P0).norm();
    }
    
    /**********************************************************************/
    template <int dim>
    bool FiniteLineSegment<dim>::hasZeroLength() const
    {
        return length()<FLT_EPSILON;
    }

    /**********************************************************************/
    template <int dim>
    typename FiniteLineSegment<dim>::VectorDim FiniteLineSegment<dim>::center() const
    {
        return 0.5*(P0+P1);
    }

//    template class SegmentSegmentDistance<1>;
//    template class SegmentSegmentDistance<2>;
//    template class SegmentSegmentDistance<3>;
    template struct FiniteLineSegment<1>;
    template struct FiniteLineSegment<2>;
    template struct FiniteLineSegment<3>;
}
#endif
