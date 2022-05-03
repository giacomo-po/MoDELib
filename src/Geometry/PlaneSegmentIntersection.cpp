/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneSegmentIntersection_CPP_
#define model_PlaneSegmentIntersection_CPP_

#include <tuple>
#include <Eigen/Dense>
#include <Plane.h>
#include <FiniteLineSegment.h>

#include <PlaneSegmentIntersection.h>

namespace model
{

    /**********************************************************************/
    template <int dim>
    typename PlaneSegmentIntersection<dim>::SolutionType PlaneSegmentIntersection<dim>::findIntersection(const VectorDimD& P0,
                                         const VectorDimD& Normal,
                                         const VectorDimD& v0,
                                         const VectorDimD& v1)
    {
        const double nNorm(Normal.norm());
        assert(nNorm>FLT_EPSILON && "PLANE MUST HAVE NON-ZERO NORMAL");
        const VectorDimD n(Normal/nNorm);

        // check intersection of v0->v1 with plane
        // x=v0+u(v1-v0)
        // (x-P0).n=0
        // (v0+u(v1-v0)-P0).n=0
        // u=(P0-v0).n/(v1-v0).n;
        const double edgeNorm=(v1-v0).norm();
        if(edgeNorm<FLT_EPSILON)
        {
            if(fabs((P0-0.5*(v0+v1)).dot(n))<FLT_EPSILON)
            {
                return std::make_tuple(INCIDENT,v0,v1);
            }
            else
            {
                return std::make_tuple(OFFSET,VectorDimD::Zero(),VectorDimD::Zero());
            }
        }
        else
        {
            const double den=(v1-v0).dot(n);
            const double num=(P0-v0).dot(n);

            if (fabs(den/edgeNorm)>FLT_EPSILON)
            {
                // edge intersects plane
                const double u=num/den;
                
                if(fabs(u)<FLT_EPSILON)
                {
                    return std::make_tuple(INCIDENT,v0,v0);
                }
                else if (u>=FLT_EPSILON && u<=1.0-FLT_EPSILON)
                {
                    const VectorDimD x((1.0-u)*v0 + u*v1);
                    return std::make_tuple(INCIDENT,x,x);
                }
                else if (fabs(1.0-u)<FLT_EPSILON)
                {
                    return std::make_tuple(INCIDENT,v1,v1);
                }
                else
                {// no roots
                    return std::make_tuple(OFFSET,VectorDimD::Zero(),VectorDimD::Zero());
                }
                
            }
            else
            {
                const double P0v0norm=(P0-v0).norm();
                const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
                if (fabs(numCheck)>FLT_EPSILON)
                {// edge is parallel to plane, no intersection
                    return std::make_tuple(PARALLEL,VectorDimD::Zero(),VectorDimD::Zero());

                }
                else
                {// edge is coplanar
                    return std::make_tuple(COINCIDENT,v0,v1);
                }
            }
        }
    }

    /**********************************************************************/
    template <int dim>
    PlaneSegmentIntersection<dim>::PlaneSegmentIntersection(const VectorDimD& P0,
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
    template <int dim>
    PlaneSegmentIntersection<dim>::PlaneSegmentIntersection(const Plane<dim>& plane,
                             const FiniteLineSegment<dim>& seg) :
    /* init */ sol(findIntersection(plane.P,plane.unitNormal,seg.P0,seg.P1)),
    /* init */ type(std::get<0>(sol)),
    /* init */ x0(std::get<1>(sol)),
    /* init */ x1(std::get<2>(sol))
    {
        
    }

    template struct PlaneSegmentIntersection<3>;
}
#endif
