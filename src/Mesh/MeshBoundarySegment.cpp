/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshBoundarySegment_cpp_
#define model_MeshBoundarySegment_cpp_

#include <cfloat>
#include <set>
#include <Eigen/Dense>
#include <FiniteLineSegment.h>

#include <MeshModule.h>

namespace model
{

template <int dim>
typename MeshBoundarySegment<dim>::VectorDim MeshBoundarySegment<dim>::periodicShift() const
{
    VectorDim temp(VectorDim::Zero());
    for(const auto& face : faces)
    {
        if(face->periodicFacePair.second)
        {
            temp+=face->periodicFacePair.first;
        }
    }
    return temp;
}


    template <int dim>
    typename MeshBoundarySegment<dim>::VectorDim MeshBoundarySegment<dim>::getBoundaryNormal(const typename MeshBoundarySegment<dim>::FaceContainerType& fcs)
    {
        assert(fcs.size() && "EMPY FACE CONTAINER");
        typename MeshBoundarySegment<dim>::VectorDim temp(MeshBoundarySegment<dim>::VectorDim::Zero());
        for(const auto& face : fcs)
        {
            temp+=face->outNormal();
        }
        const double tempNorm(temp.norm());
        return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
    }
    
    template <int dim>
    MeshBoundarySegment<dim>::MeshBoundarySegment(const typename MeshBoundarySegment<dim>::VectorDim& p0,
                                                  const typename MeshBoundarySegment<dim>::VectorDim& p1,
                                                  const typename MeshBoundarySegment<dim>::PlanarMeshFaceType* const face) :
    /* init */ FiniteLineSegmentType(p0,p1)
    /* init */,faces(FaceContainerType({face}))
    /* init */,outNormal(getBoundaryNormal(faces))
    {
    }
    
    template <int dim>
    MeshBoundarySegment<dim>::MeshBoundarySegment(const typename MeshBoundarySegment<dim>::VectorDim& p0,
                                                  const typename MeshBoundarySegment<dim>::VectorDim& p1,
                                                  const typename MeshBoundarySegment<dim>::FaceContainerType& faces_in) :
    /* init */ FiniteLineSegmentType(p0,p1)
    /* init */,faces(faces_in)
    /* init */,outNormal(getBoundaryNormal(faces))
    {
    }
        
    /**********************************************************************/
    template <int dim>
    void BoundingMeshSegments<dim>::emplace_unique(const VectorDim& P0,const VectorDim& P1,const PlanarMeshFace<dim>* const face)
    {
        if((P0-P1).norm()>FLT_EPSILON)
        {// ignore degenerate segments
            int merged(0);
            for(auto& seg : *this)
            {
                if(   ((seg->P0-P0).norm()<FLT_EPSILON && (seg->P1-P1).norm()<FLT_EPSILON)
                   || ((seg->P0-P1).norm()<FLT_EPSILON && (seg->P1-P0).norm()<FLT_EPSILON)
                   )
                {// coincident segments on different faces, need to merge faces
                    FaceContainerType mergedFaces(seg->faces);
                    mergedFaces.insert(face);
                    seg.reset(new MeshBoundarySegment<dim>(P0,P1,mergedFaces));
                    merged++;
                }
            }
            if (merged)
            {// P0->P1 was found to be not unique. Swap vectors
                assert(merged==1 && "MORE THAN ONE SEGMENTS FOUND ON LINE P0 P1");
            }
            else
            {// P0->P1 was found to be unique, add it
                this->emplace_back(new MeshBoundarySegment<dim>(P0,P1,face));
            }
        }
    }
    
    template <int dim>
    void BoundingMeshSegments<dim>::computeFaceIntersections(const Plane<dim>& plane,
                                                             const std::shared_ptr<PlanarMeshFace<dim>>& face)
    {
        PlanePlaneIntersection<dim> ppi(plane,face->asPlane());
                
        if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
        {// plane and mesh-face are incident
            
            UniquePointContainer roots;
            for(size_t k=0;k<face->convexHull().size();++k)
            {
                
                // Pick face boundary segment
                const size_t k1(k==face->convexHull().size()-1? 0 : k+1);
                const VectorDim& P0(face->convexHull()[k]->P0);
                const VectorDim& P1(face->convexHull()[k1]->P0);
                const double segLength((P1-P0).norm());
                const VectorDim D0((P1-P0)/segLength);
                // Compute intersection between boundary segment and line of intersection of the two planes
                LineLineIntersection<dim> lli(P0,D0,ppi.P,ppi.d);
                if(lli.type==LineLineIntersection<dim>::INCIDENT)
                {
                    const double u0((lli.x0-P0).dot(D0));
                    if(u0>0.0-FLT_EPSILON && u0<segLength+FLT_EPSILON)
                    {// intersection within segment
                        roots.insert(lli.x0);
                    }
                }
                else if(lli.type==LineLineIntersection<dim>::COINCIDENT)
                {// a coincident line was found, which means that the glide planes intersec on a boundary face
                    roots.insert(P0);
                    roots.insert(P1);
                }
            }
            
            switch (roots.size())
            {
                    case 0:
                {// no intersaction between plane and face
                    break;
                }
                    
                    case 1:
                {// single-point intersection, not a valid boundary segment, so don't consider it
                    break;
                }
                    
                    case 2:
                {// an intersection segment
                    const VectorDim& P0(*roots.begin());
                    const VectorDim& P1(*roots.rbegin());
                    emplace_unique(P0,P1,face.get());
                    break;
                }
                    
                default:
                {
                    std::cout<<"FAILED TO FIND A BOUNDARY PERIMETER FOR PLANE"<<std::endl;
                    std::cout<<"plane.P="<<std::setprecision(15)<<std::scientific<<plane.P.transpose()<<std::endl;
                    std::cout<<"plane.unitNormal="<<std::setprecision(15)<<std::scientific<<plane.unitNormal.transpose()<<std::endl;
                    std::cout<<"IN INTERSECTING PLANE AND FACE"<<std::endl;
                    std::cout<<face->sID<<", n="<<std::setprecision(15)<<std::scientific<<face->outNormal()<<std::endl;
                    std::cout<<"ROOTS ARE"<<std::endl;
                    for(const auto& root : roots)
                    {
                        std::cout<<std::setprecision(15)<<std::scientific<<root.transpose()<<std::endl;
                    }
                    assert(false && "FAILED TO FIND A BOUNDARY PERIMETER FOR PLANE");
                    break;
                }
            }
            
        }
        else if (ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
        {
            for(size_t k=0;k<face->convexHull().size();++k)
            {
                //std::cout<<"I'm here A8"<<std::endl;
                
                const size_t k1(k==face->convexHull().size()-1? 0 : k+1);
                const VectorDim& P0(face->convexHull()[k]->P0);
                const VectorDim& P1(face->convexHull()[k1]->P0);
                emplace_unique(P0,P1,face.get());
            }
        }
    }
    
    /**********************************************************************/
    template <int dim>
    BoundingMeshSegments<dim>::BoundingMeshSegments()
    {
        
    }
    
    
    
    /**********************************************************************/
    template <int dim>
    BoundingMeshSegments<dim>::BoundingMeshSegments(const SimplicialMesh<dim>& mesh,
                                                    const Plane<dim>& plane)
    {
        for(const auto& region : mesh.regions())
        {
            for(const auto& face : region.second->faces())
            {
                if(face.second->isExternal())
                {
                    computeFaceIntersections(plane,face.second);
                }
            }
        }
    }
    
    /**********************************************************************/
    template <int dim>
    BoundingMeshSegments<dim>::BoundingMeshSegments(const SimplicialMesh<dim>& mesh,
                                                    const size_t& rID,
                                                    const Plane<dim>& plane)
    {
        //            const MatrixDim R(plane.localRotationMatrix());
        
        for(const auto& face : mesh.region(rID)->faces())
        {
            computeFaceIntersections(plane,face.second);
        }
        

        
        // Now sort segments
        ConvexHull<2,std::shared_ptr<MeshBoundarySegment<dim>>> finalHull;
        //std::cout<<"unsorted hull"<<std::endl;
        for(const auto& pt : *this)
        {
            const auto x(plane.localPosition(0.5*(pt->P0+pt->P1)));
            //                VectorDim x(R*(0.5*(pt->P0+pt->P1)-plane.P));
            finalHull.emplace(std::array<double,2>{x[0],x[1]},&pt);
            // THE PROBLEM HERE IS THAT IF COINCIDENT POINTS FROM DIFFERENCE FACES EXIST, THEN ONLY ONE OF THEM IS KEPT. E.G. A PLANE CUTTING AT THE INTERSECTION OF TWO FACES. IF WE HAD UNIQUE FACE EDGES WITH POINTERS TO THE ADJECENT FACES WE COULD SOLVE THIS
        }
        

        
        const auto hullPts=finalHull.getPoints();
        BoundingMeshSegments<dim> sortedTemp;
        for(size_t k=0;k<hullPts.size();++k)
        {
            const auto& seg(*hullPts[k].t);
            if(k==0)
            {
                const auto& seg1(*hullPts[1].t);
                
                if((seg->P1-seg1->P0).norm()<FLT_EPSILON || (seg->P1-seg1->P1).norm()<FLT_EPSILON)
                {// first segment is oriented in the right direction
                    sortedTemp.push_back(seg);
                }
                else
                {
                    sortedTemp.emplace_back(new MeshBoundarySegment<dim>(seg->P1,seg->P0,seg->faces));
                }
            }
            else
            {
                if((sortedTemp.back()->P1-seg->P0).norm()<FLT_EPSILON)
                {
                    sortedTemp.push_back(seg);
                }
                else if((sortedTemp.back()->P1-seg->P1).norm()<FLT_EPSILON)
                {
                    sortedTemp.emplace_back(new MeshBoundarySegment<dim>(seg->P1,seg->P0,seg->faces));
                }
                else
                {
                    std::cout<<(sortedTemp.back()->P1-seg->P0).norm()<<std::endl;
                    std::cout<<(sortedTemp.back()->P1-seg->P1).norm()<<std::endl;
                    for(size_t k=0;k<hullPts.size();++k)
                    {
                        const auto& segTemp(*hullPts[k].t);
                        
                        std::cout<<segTemp->P0.transpose()<<" "<<segTemp->P1.transpose()<<std::endl;
                    }
                    assert(false && "DISCONNECTED FACE BOUNDARY");
                }
            }
        }

        
        if(sortedTemp.size())
        {
            
            assert((sortedTemp.back()->P1-sortedTemp.front()->P0).norm()<FLT_EPSILON && "OPEN FACE BOUNDARY");

            
            //                VectorDim nA(VectorDim::Zero());
            //                const VectorDim P0(sortedTemp.front().P0);
            //                for(const auto& seg : sortedTemp)
            //                {
            //                    nA+= 0.5*(seg.P0-P0).cross(seg.P1-seg.P0);
            //                }
            //                if(!isRightHandedBoundary(sortedTemp,plane))
            //                {// boundary makes a left-handed loop with respect to plane normal. We need a right-handed loop
            //                    BoundingMeshSegments<dim> flippedTemp;
            //                    for (auto it = sortedTemp.rbegin(); it != sortedTemp.rend(); ++it)
            //                    {
            //                        flippedTemp.emplace_back(it->P1,it->P0,it->faces);
            //                    }
            //                    sortedTemp.swap(flippedTemp);
            //                }
            //                assert(isRightHandedBoundary(sortedTemp,plane));
        }
        
        
        assert(sortedTemp.size()==this->size());

        this->swap(sortedTemp);

    }
    
    /**********************************************************************/
    template <int dim>
    std::set<const MeshBoundarySegment<dim>*> BoundingMeshSegments<dim>::containingSegments(const VectorDim& P) const
    {
        std::set<const MeshBoundarySegment<dim>*> temp;
        
        for(const auto& seg : *this)
        {
            if(seg->contains(P))
            {
                temp.insert(seg.get());
            }
        }
        return temp;
    }
    
    /**********************************************************************/
    template <int dim>
    bool BoundingMeshSegments<dim>::contains(const VectorDim& P) const
    {
        return containingSegments(P).size();
    }
    
    /**********************************************************************/
    template <int dim>
    const BoundingMeshSegments<dim>& BoundingMeshSegments<dim>::boundingBoxSegments() const
    {
        return *this;
    }
    
    template <int dim>
    BoundingMeshSegments<dim>& BoundingMeshSegments<dim>::boundingBoxSegments()
    {
        return *this;
    }
    
    
    //
    //
    //
    //
    //        /**********************************************************************/
    //        template <class T>
    //        friend T& operator << (T& os, const BoundingMeshSegments<dim>& bls)
    //        {
    //            for(const auto& seg : bls)
    //            {
    //                os<<*seg;
    //            }
    //            return os;
    //        }
    
    template class MeshBoundarySegment<3>;
    template struct BoundingMeshSegments<3>;

    
}
#endif
