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
//#include <map>
#include <Eigen/Dense>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Geometry/Plane.h>
#include <model/Utilities/StaticID.h>
#include <model/Mesh/PlaneMeshIntersection.h>

namespace model
{
    
    template <int dim>
    struct MeshPlane : public StaticID<MeshPlane<dim>>,
    /*              */ public Plane<dim>
    {
        
        
        
        typedef std::array<long int,dim+3> MeshPlaneKeyType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;

        
        /**********************************************************************/
        static Plane<dim> getPlane(const SimplicialMesh<dim>& mesh,
                                        const int& rID1,
                                        const int& rID2)
        {
            
            const auto& regionBnd(mesh.regionBoundary(rID1,rID2));
            const Simplex<dim,dim-1>& triangle(**regionBnd.begin());
            const Simplex<dim,dim>& tet(**triangle.parents().begin());
            const size_t faceID=tet.childOrder(triangle.xID);
            const VectorDim N=tet.nda.col(faceID);
            const VectorDim P=triangle.vertexPositionMatrix().col(0);
            
            Plane<dim> plane(P,N);
            
            //            for(const auto& tet : triangle.parents())
            //            {
            //                const size_t faceID=tet->childOrder(triangle.xID);
            //                const VectorDim outNormal=tet->nda.col(faceID);
            //                const double outNorm(outNormal.norm());
            //                assert(outNorm>FLT_EPSILON && "Simplex normal has zero norm.");
            //                storeLatticePlane(dn,mesh,grain(tet->region->regionID),outNormal/outNorm);
            //            }
            //
            //            assert(glidePlanes().size()==2);
            //            assert((glidePlanes().begin()->second->unitNormal+glidePlanes().rbegin()->second->unitNormal).norm()<FLT_EPSILON && "plane normals are not opposite to each other.");
            
            // Check that all tringle vertices are contained by both GB planes
            for(const auto& triangle : regionBnd)
            {
                const auto Ps=triangle->vertexPositionMatrix();
                for(int j=0;j<Ps.cols();++j)
                {
                    assert(plane.contains(Ps.col(j)) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
                    //                    for(const auto& grain : grains())
                    //                    {
                    //                        assert(glidePlane(grain.second->grainID).contains(Ps.col(j)) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
                    //                    }
                }
            }
            
            return plane;
        }
        
//        const MeshPlaneKeyType planeKey;
//        const int regionID;
        const std::pair<int,int> regionIDs;
        const PlaneMeshIntersectionContainerType meshIntersections;
        
        /**********************************************************************/
        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const int& rID,
                  const VectorDim& p,
                  const VectorDim& n) :
        /* init */ Plane<dim>(p,n),
        /* init */ regionIDs(rID,rID),
//        /* init */ regionID(rID),
        /* init */ meshIntersections(PlaneMeshIntersection<dim>(mesh,this->P,this->unitNormal,rID))

        {/*!\param[in] mesh
          * \param[in] rID the region ID where the plane is defined
          * \param[in] p position of the plane
          * \param[in] n normal to the plane
          * Constructor for plane internal to a mesh region
          */
        }
        
        /**********************************************************************/
        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const int& rID1,
                  const int& rID2) :
        /* init */ Plane<dim>(getPlane(mesh,rID1,rID2)),
//        /* init */ regionID(rID),
        /* init */ regionIDs(rID1,rID2),
        /* init */ meshIntersections(PlaneMeshIntersection<dim>(mesh,this->P,this->unitNormal,rID1)) // WARNING: CALLING meshIntersections with rID1
        
        {/*!\param[in] mesh
          * \param[in] rID the region ID where the plane is defined
          * \param[in] p position of the plane
          * \param[in] n normal to the plane
          * Constructor for plane between two mesh regions
          */
        }

    };
    
}
#endif



//        /**********************************************************************/
//        static int sign(const long int& i)
//        {
//            if(i>0)
//            {
//                return 1;
//            }
//            else if(i<0)
//            {
//                return -1;
//            }
//            else
//            {
//                return 0;
//            }
//
//        }

//        /**********************************************************************/
//        static MeshPlaneKeyType getMeshPlaneKey(const Lattice<dim>& lattice,
//                                                  const int& grainID1,
//                                                  //                                                  const int& grainID2,
//                                                  const VectorDim& P,
//                                                  const VectorDim& N)
//        /**********************************************************************/
//        static MeshPlaneKeyType getMeshPlaneKey(const int& grainID,
//                                                const LatticePlane& lp
////                                                const VectorDim& N
//        )
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
////            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
////            //            return (MeshPlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
////            //            return (MeshPlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
////            const long int h(LatticePlane::height(r,P));
//            MeshPlaneKeyType temp;
//            temp[0]=grainID;
//            temp[1]=grainID;
////            const int signh(sign(h));
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=lp.n(d);
//            }
//            temp[2+dim]=lp.h;
//            return temp;
//            //            return (MeshPlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }
//
//        /**********************************************************************/
//        static MeshPlaneKeyType getMeshPlaneKey(const int& grainID1,
//                                                  const int& grainID2)
//        {/*!\param[in] grain the grain on which the GlidePlane is defined
//          * \param[in] P a point on the plane
//          * \param[in] N the normal to the plane
//          * \returns the key which uniquely identifies the plane.
//          * The type of the key is a tuple with entries (grainID,r,h), where r
//          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
//          * is an integer indicating the "heigth" of the plane from the origin,
//          * in integer multiples of the interplanar distance d=1/|r|.
//          */
//            //            const ReciprocalLatticeDirection<dim> r(lattice.reciprocalLatticeDirection(N));
//            //            return (MeshPlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            return (MeshPlaneKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
//            //            const long int h(LatticePlane::height(r,P));
//            MeshPlaneKeyType temp;
//            temp[0]=grainID1;
//            temp[1]=grainID2;
//            //            const int signh(sign(h));
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=0;
//            }
//            temp[2+dim]=0;
//            return temp;
//            //            return (MeshPlaneKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }
