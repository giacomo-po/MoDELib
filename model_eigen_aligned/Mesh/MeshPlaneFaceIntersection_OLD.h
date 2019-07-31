/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshPlaneFaceIntersection_H_
#define model_MeshPlaneFaceIntersection_H_

#include <cfloat>
#include <tuple>
#include <Eigen/Dense>
#include <SimplicialMeshFace.h>
#include <LineSegment.h>

namespace model
{
    
    template <int dim>
    struct MeshPlaneIntersection : public LineSegment<dim>
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef LineSegment<dim> LineSegmentType;
        typedef SimplicialMeshFace<dim> SimplicialMeshFaceType;

        
        const SimplicialMeshFaceType* const face;
        
        MeshPlaneFaceIntersection(const VectorDim& p0,
                            const VectorDim& p1,
                            const SimplicialMeshFaceType* const face_in) :
        /* init */ LineSegmentType(p0,p1)
        /* init */ face(face_in)
        {
            
            
        }
        
        
    };
    

    
}
#endif
