/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshBoundarySegment_H_
#define model_MeshBoundarySegment_H_

#include <cfloat>
#include <set>
#include <Eigen/Dense>
#include <LineSegment.h>
#include <PlanarMeshFace.h>

namespace model
{
    /*!\brief Class template representing a straight line at the intersection
     * of PlanarMeshFaces
     */
    template <int dim>
    struct MeshBoundarySegment : public LineSegment<dim>
//    /*                        */,public std::set<const PlanarMeshFace<dim>*>
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef LineSegment<dim> LineSegmentType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFace<dim>*> FaceContainerType;
        
        
        static VectorDim getBoundaryNormal(const FaceContainerType& fcs)
        {
            assert(fcs.size() && "EMPY FACE CONTAINER");
            VectorDim temp(VectorDim::Zero());
            for(const auto& face : fcs)
            {
                temp+=face->outNormal();
            }
            const double tempNorm(temp.norm());
            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
        }
        
        const FaceContainerType faces;
        const VectorDim boundaryNormal;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        //        const PlanarMeshFaceType* const face;
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const PlanarMeshFaceType* const face) :
        /* init */ LineSegmentType(p0,p1)
        /* init */,faces(FaceContainerType({face}))
        /* init */,boundaryNormal(getBoundaryNormal(faces))
        //        /* init */,face(face_in)
        {
//            assert(face.size() && "EMPY FACE CONTAINER");
        }
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const FaceContainerType& faces_in) :
        /* init */ LineSegmentType(p0,p1)
        /* init */,faces(faces_in)
        /* init */,boundaryNormal(getBoundaryNormal(faces))
        {
//            assert(faces().size() && "EMPY FACE CONTAINER");
        }
        
//        FaceContainerType& faces()
//        {
//            return *this;
//
//        }
//
//        const FaceContainerType& faces() const
//        {
//            return *this;
//
//        }
        
//        /**********************************************************************/
//        VectorDim boundaryNormal() const
//        {
//            VectorDim temp(VectorDim::Zero());
//            for(const auto& face : faces())
//            {
//                temp+=face->outNormal();
//            }
//            const double tempNorm(temp.norm());
//            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
//        }
        
        
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const MeshBoundarySegment<dim>& seg)
        {
            os<<seg.faces.size()<<", "<<std::setprecision(15)<<std::scientific<<seg.P0.transpose()<<","<<seg.P1.transpose();
            return os;
        }
        
        
    };
    
    template <int dim>
    struct BoundingMeshSegments : public std::vector<MeshBoundarySegment<dim>, Eigen::aligned_allocator<MeshBoundarySegment<dim>>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        
//        BoundingMeshSegments(){}
        
        //        BoundingMeshSegments(const BoundingMeshSegments<dim>& bb1,const BoundingMeshSegments<dim>& bb2)
        //        {
        //            assert(0 && "FINISH HERE");
        //        }
        
        /**********************************************************************/
        std::set<const MeshBoundarySegment<dim>*> containingSegments(const VectorDim& P) const
        {
            std::set<const MeshBoundarySegment<dim>*> temp;
            
            for(const auto& seg : *this)
            {
                //                std::cout<<pair.first<<", d="<<pair.second.distanceTo(P)<<std::endl;
                if(seg.contains(P))
                {
                    temp.insert(&seg);
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool contains(const VectorDim& P) const
        {
            return containingSegments(P).size();
        }
        
        /**********************************************************************/
        const BoundingMeshSegments<dim>& boundingBoxSegments() const
        {
            return *this;
        }
        
        BoundingMeshSegments<dim>& boundingBoxSegments()
        {
            return *this;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const BoundingMeshSegments<dim>& bls)
        {
            for(const auto& seg : bls)
            {
                os<<seg<<std::endl;
            }
            return os;
        }
        
    };
    
    
//    template<int dim>
//    struct MeshBoundarySegmentIntersection : public std::unique_ptr<MeshBoundarySegment<dim>>
//    {
//
//
//        MeshBoundarySegmentIntersection(const MeshBoundarySegment<dim>& s1,const MeshBoundarySegment<dim>& s2)
//        {// Constructs the intersection of two MeshBoundarySegment(s)
//            const SegmentSegmentDistance<dim> ssd(s1.P0,s1.P1,s2.P0,s2.P1);
//            const auto iSeg(ssd.intersectionSegment());
//
//            switch (iSeg.size())
//            {
//                    std::set<const PlanarMeshFace<dim>*> allFaces;
//                    for(const auto& tup : iSeg)
//                    {// Make sure that each face cointais end points, and add to allFaces
//                        const auto& x(std::get<0>(tup));
//                        for(const auto& face : s1.faces)
//                        {
//                            assert(face->asPlane().contains(x) && "FACE DOES NOT CONTAIN INTERSECTION POINT");
//                            allFaces.insert(face);
//                        }
//                        for(const auto& face : s2.faces)
//                        {
//                            assert(face->asPlane().contains(x) && "FACE DOES NOT CONTAIN INTERSECTION POINT");
//                            allFaces.insert(face);
//                        }
//                    }
//
//                case 1:
//                {// Single intersection point. This point belongs must be common to all faces of the origina segments
//                    const auto& x(std::get<0>(iSeg[0]));
//                    this->reset(new MeshBoundarySegment<dim>(x,x,allFaces));
//                    break;
//                }
//
//                case 2:
//                {// extended intersection segment
//                    const auto& x0(std::get<0>(iSeg[0]));
//                    const auto& x1(std::get<0>(iSeg[1]));
//                    this->reset(new MeshBoundarySegment<dim>(x0,x1,allFaces));
//                    break;
//                }
//
//                default:
//                    break;
//            }
//
//
//
//        }
//
//        MeshBoundarySegmentIntersection(const MeshBoundarySegment<dim>& seg,const PlanarMeshFace<dim>& face)
//        {// Constructs the intersection of a MeshBoundarySegment and a PlanarMeshFace
//
//            const PlaneSegmentIntersection<dim> psi(face.asPlane(),seg);
//            if(psi.type==PlaneSegmentIntersection<dim>::COINCIDENT || psi.type==PlaneSegmentIntersection<dim>::INCIDENT)
//            {
//                std::set<const PlanarMeshFace<dim>*> allFaces(seg.faces);
//                allFaces.insert(face);
//                this->reset(new MeshBoundarySegment<dim>(psi.x0,psi.x1,allFaces));
//            }
////
////            switch (psi.type)
////            {
////                case PlaneSegmentIntersection::COINCIDENT:
////                {
////                    std::set<const PlanarMeshFace<dim>*> allFaces(seg.faces);
////                    allFaces.insert(face);
////                    this->reset(new MeshBoundarySegment<dim>(psi.x0,psi.x1,allFaces));
////                    break;
////                }
////
////                case PlaneSegmentIntersection::INCIDENT:
////                {
////                    this->reset(new MeshBoundarySegment<dim>(psi.x0,psi.x1,allFaces));
////                    break;
////                }
////
////                default:
////                    break;
////            }
//        }
//
//    };
    
    
}
#endif

