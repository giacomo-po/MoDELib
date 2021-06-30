/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneSegmentIntersection_H_
#define model_PlaneSegmentIntersection_H_

#include <tuple>
#include <Eigen/Dense>
#include <Plane.h>
#include <FiniteLineSegment.h>

namespace model
{
    
    template <int dim>
    struct PlaneSegmentIntersection
    {
        
        enum IntersectionType {PARALLEL,COINCIDENT,INCIDENT,OFFSET};
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        typedef std::tuple<IntersectionType,VectorDimD,VectorDimD> SolutionType;
        const SolutionType sol;
        
        /**********************************************************************/
        static SolutionType findIntersection(const VectorDimD& P0,
                                             const VectorDimD& Normal,
                                             const VectorDimD& v0,
                                             const VectorDimD& v1);
        
        
    public:
        
        const IntersectionType& type;
        const VectorDimD x0;
        const VectorDimD x1;
        
        /**********************************************************************/
        PlaneSegmentIntersection(const VectorDimD& P0,
                                 const VectorDimD& N,
                                 const VectorDimD& v0,
                                 const VectorDimD& v1) :
        /* init */ sol(findIntersection(P0,N,v0,v1)),
        /* init */ type(std::get<0>(sol)),
        /* init */ x0(std::get<1>(sol)),
        /* init */ x1(std::get<2>(sol))
        {
            
        }
        
        /**********************************************************************/
        PlaneSegmentIntersection(const Plane<dim>& plane,
                                 const FiniteLineSegment<dim>& seg) :
        /* init */ sol(findIntersection(plane.P,plane.unitNormal,seg.P0,seg.P1)),
        /* init */ type(std::get<0>(sol)),
        /* init */ x0(std::get<1>(sol)),
        /* init */ x1(std::get<2>(sol))
        {
            
        }
        
    };

}
#endif
