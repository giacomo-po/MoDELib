/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshBoundarySegment_H_
#define model_MeshBoundarySegment_H_

#include <iomanip>
#include <cfloat>
#include <set>
#include <Eigen/Dense>
#include <FiniteLineSegment.h>
#include <PlanarMeshFace.h>
#include <LineLineIntersection.h>
#include <PlanePlaneIntersection.h>

namespace model
{
    /*!\brief Class template representing a straight line at the intersection
     * of PlanarMeshFaces
     */
    template <int dim>
    class MeshBoundarySegment : public FiniteLineSegment<dim>
    {
        
    public:
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef FiniteLineSegment<dim> FiniteLineSegmentType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFace<dim>*> FaceContainerType;
        
    private:
        
        static VectorDim getBoundaryNormal(const FaceContainerType& fcs);
        
    public:
        
        const FaceContainerType faces;
        const VectorDim outNormal;
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const PlanarMeshFaceType* const face);
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const FaceContainerType& faces_in);
        
        VectorDim periodicShift() const;

    };

    template <int dim,class T>
    T& operator << (T& os, const MeshBoundarySegment<dim>& seg)
    {
        os<<std::setprecision(15)<<std::scientific<<seg.P0.transpose()<<","<<seg.P1.transpose()<<"\n";
        return os;
    }

    
    template <int dim>
    struct BoundingMeshSegments : public std::vector<std::shared_ptr<MeshBoundarySegment<dim>>>
    {
        
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::set<Eigen::Matrix<double,dim,1>,CompareVectorsByComponent<double,dim,float> > UniquePointContainer;
        typedef MeshBoundarySegment<dim> MeshBoundarySegmentType;
        typedef typename MeshBoundarySegmentType::FaceContainerType FaceContainerType;
        typedef std::vector<std::shared_ptr<MeshBoundarySegment<dim>>> MeshBoundarySegmentContainerType;
        
        void emplace_unique(const VectorDim& P0,const VectorDim& P1,const PlanarMeshFace<dim>* const face);
        void computeFaceIntersections(const Plane<dim>& plane,
                                      const std::shared_ptr<PlanarMeshFace<dim>>& face);
        BoundingMeshSegments();
        BoundingMeshSegments(const SimplicialMesh<dim>& mesh,
                             const Plane<dim>& plane);
        BoundingMeshSegments(const SimplicialMesh<dim>& mesh,
                             const size_t& rID,
                             const Plane<dim>& plane);
        std::set<const MeshBoundarySegment<dim>*> containingSegments(const VectorDim& P) const;
        bool contains(const VectorDim& P) const;
        const BoundingMeshSegments<dim>& boundingBoxSegments() const;
        BoundingMeshSegments<dim>& boundingBoxSegments();
        
    };

    template <int dim,class T>
    T& operator << (T& os, const BoundingMeshSegments<dim>& bls)
    {
        for(const auto& seg : bls)
        {
            os<<*seg;
        }
        return os;
    }

    
}
#endif
