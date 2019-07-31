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

namespace model
{
    
    template <int dim>
    struct MeshPlane : public StaticID<MeshPlane<dim>>
    /*              */,public Plane<dim>
    {
        
        typedef std::array<long int,dim+3> MeshPlaneKeyType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::pair<VectorDim,const Simplex<dim,1>* const> RootType;
        typedef std::deque<RootType> RootContainerType;
        typedef std::set<Eigen::Matrix<double,dim,1>,CompareVectorsByComponent<double,dim,float> > UniquePointContainer;
        
        /**********************************************************************/
        static BoundingMeshSegments<dim> getFaceBoundary(const PlanarMeshFace<dim>& face)
        {
            BoundingMeshSegments<dim> temp;
            for(size_t k=0;k<face.convexHull().size();++k)
            {
                const size_t k1(k==face.convexHull().size()-1? 0 : k+1);
                temp.emplace_back(new MeshBoundarySegment<dim>(face.convexHull()[k]->P0,face.convexHull()[k1]->P0,&face));
            }
            return temp;
        }
        
        /**********************************************************************/
        static void checkPlaneIntersections(const BoundingMeshSegments<dim>& temp)
        {
            
            if(temp.size()<3)
            {
                std::cout<<"meshIntersections.size()="<<temp.size()<<std::endl;
                assert(false && "meshIntersections FAILED");
            }
            
            for(const auto& seg : temp)
            {
                assert(!seg->hasZeroLength() && "Plane-Face intersection has zero length");
            }
        }
        
        const std::pair<int,int> regionIDs;
        const BoundingMeshSegments<dim> meshIntersections;
        
        /**********************************************************************/
        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const int& rID,
                  const VectorDim& p,
                  const VectorDim& n) :
        /* init */ Plane<dim>(p,n)
        /* init */,regionIDs(rID,rID)
        /* init */,meshIntersections(mesh,rID,*this)
        {/*!\param[in] mesh
          * \param[in] rID the region ID where the plane is defined
          * \param[in] p position of the plane
          * \param[in] n normal to the plane
          * Constructor for plane internal to a mesh region
          */
            checkPlaneIntersections(meshIntersections);
        }
        
        /**********************************************************************/
        MeshPlane(const PlanarMeshFace<dim>& face,
                  const int& rID1,
                  const int& rID2) :
        /* init */ Plane<dim>(face.asPlane()),
        /* init */ regionIDs(rID1,rID2),
        /* init */ meshIntersections(getFaceBoundary(face)) // WARNING: CALLING meshIntersections with rID1
        {
            checkPlaneIntersections(meshIntersections);
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const MeshPlane<dim>& gp)
        {
            for (const auto& x : gp.meshIntersections)
            {
                os<<gp.sID<<" "<<*x;
            }
            return os;
        }
        
    };
    
}
#endif
