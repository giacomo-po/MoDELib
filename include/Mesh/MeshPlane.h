/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshPlane_H_
#define model_MeshPlane_H_


#include <cfloat>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <SimplicialMesh.h>
#include <PlanarMeshFace.h>
#include <Plane.h>
#include <StaticID.h>
#include <PlaneSegmentIntersection.h>
#include <PlanePlaneIntersection.h>
#include <MeshBoundarySegment.h>
#include <LineLineIntersection.h>
//#include <EmbeddedPolygonTriangulation.h>

namespace model
{
    
    
    template <int dim>
    struct MeshPlane : public StaticID<MeshPlane<dim>>
    /*              */,public Plane<dim>
    {
        
        typedef std::array<long int,dim+3> MeshPlaneKeyType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::pair<VectorDim,const Simplex<dim,1>* const> RootType;
        typedef std::deque<RootType> RootContainerType;
        typedef std::set<Eigen::Matrix<double,dim,1>,CompareVectorsByComponent<double,dim,float> > UniquePointContainer;
        
        static BoundingMeshSegments<dim> getFaceBoundary(const PlanarMeshFace<dim>& face);
        static void checkPlaneIntersections(const BoundingMeshSegments<dim>& temp);
        BoundingMeshSegments<dim> sortMeshIntersections(const BoundingMeshSegments<dim>& mshInt) const;
        
        const std::pair<int,int> regionIDs;
        const BoundingMeshSegments<dim> meshIntersections;

        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const VectorDim& p,
                  const VectorDim& n);
        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const size_t& rID,
                  const VectorDim& p,
                  const VectorDim& n);
        MeshPlane(const PlanarMeshFace<dim>& face,
                  const size_t& rID1,
                  const size_t& rID2);
    };
    
    template <int dim,class T>
    T& operator << (T& os, const MeshPlane<dim>& gp)
    {
        for (const auto& x : gp.meshIntersections)
        {
            os<<gp.sID<<" "<<*x;
        }
        return os;
    }
    
}
#endif
