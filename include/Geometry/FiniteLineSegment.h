/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FiniteLineSegment_H_
#define model_FiniteLineSegment_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>

#include <SegmentSegmentDistance.h>

namespace model
{
    
    template <int dim>
    struct FiniteLineSegment
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
                
        VectorDim P0;
        VectorDim P1;
        
        /**********************************************************************/
        FiniteLineSegment(const VectorDim& p0,
                          const VectorDim& p1);
        
        /**********************************************************************/
        VectorDim snap(const VectorDim& x) const;

        /**********************************************************************/
        VectorDim snapToInfiniteLine(const VectorDim& x) const;

        

        
        /**********************************************************************/
        SegmentSegmentDistance<dim> distanceTo(const FiniteLineSegment<dim>& other,
                                               const double& degeneracyValue=0.0) const;
        
        /**********************************************************************/
        double distanceTo(const VectorDim& x) const;
        
        /**********************************************************************/
        bool contains(const VectorDim& x) const;
        
        /**********************************************************************/
        double length() const;
        
        /**********************************************************************/
        bool hasZeroLength() const;
        
        /**********************************************************************/
        VectorDim center() const;
        
    };
    
}
#endif
