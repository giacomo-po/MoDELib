/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshPlaneIntersection_H_
#define model_MeshPlaneIntersection_H_

#include <cfloat>
#include <tuple>
//#include <map>
#include <Eigen/Dense>
#include <SimplicialMesh.h>
#include <Plane.h>
#include <StaticID.h>
//#include <PlaneMeshIntersection.h>
//#include <MeshPlaneIntersectionFaceIntersection.h>
namespace model
{
    
    
    
    template <int dim>
    class MeshPlaneIntersection
    {
        
        
    public:
        MeshPlaneIntersection(const SimplicialMesh<dim>& mesh,
                              const int& rID1,
                              const int& rID2)
        {
            
            assert(rID1!=rID2);
            
            mesh.regionBonudary(rID1,rID2)
            
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
//        static MeshPlaneIntersectionKeyType getMeshPlaneIntersectionKey(const Lattice<dim>& lattice,
//                                                  const int& grainID1,
//                                                  //                                                  const int& grainID2,
//                                                  const VectorDim& P,
//                                                  const VectorDim& N)
//        /**********************************************************************/
//        static MeshPlaneIntersectionKeyType getMeshPlaneIntersectionKey(const int& grainID,
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
////            //            return (MeshPlaneIntersectionKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
////            //            return (MeshPlaneIntersectionKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
////            const long int h(LatticePlane::height(r,P));
//            MeshPlaneIntersectionKeyType temp;
//            temp[0]=grainID;
//            temp[1]=grainID;
////            const int signh(sign(h));
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=lp.n(d);
//            }
//            temp[2+dim]=lp.h;
//            return temp;
//            //            return (MeshPlaneIntersectionKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }
//
//        /**********************************************************************/
//        static MeshPlaneIntersectionKeyType getMeshPlaneIntersectionKey(const int& grainID1,
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
//            //            return (MeshPlaneIntersectionKeyType()<<grainID1,grainID2,r,LatticePlane::height(LatticePlane::computeHeight(r,P))).finished();
//            //            return (MeshPlaneIntersectionKeyType()<<grainID1,grainID2,r,LatticePlane::height(r,P)).finished();
//            //            const long int h(LatticePlane::height(r,P));
//            MeshPlaneIntersectionKeyType temp;
//            temp[0]=grainID1;
//            temp[1]=grainID2;
//            //            const int signh(sign(h));
//            for(int d=0;d<dim;++d)
//            {
//                temp[2+d]=0;
//            }
//            temp[2+dim]=0;
//            return temp;
//            //            return (MeshPlaneIntersectionKeyType()<<grainID1,grainID2,r*sign(h),h*sign(h)).finished(); // make sure that key heigh is always positive for uniqueness
//
//        }
