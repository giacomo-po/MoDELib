/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundingLineSegments_H_
#define model_BoundingLineSegments_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/StdDeque>
#include <deque>
#include <map>
#include <utility>
//#include <GlidePlane.h>
//#include <Plane.h>
#include <SegmentSegmentDistance.h>
#include <PlaneSegmentIntersection.h>
#include <LineSegment.h>
#include <MeshPlane.h>

//#include <SegmentSegmentIntersection.h>

namespace model
{
    
    template <int dim>
    class BoundingLineSegments : public std::map<int,MeshBoundarySegment<dim>>
    {
        
    public:
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        //        typedef std::pair<VectorDim,VectorDim> MeshBoundarySegment<dim>;
        //typedef std::deque<MeshBoundarySegment<dim>/*,Eigen::aligned_allocator<MeshBoundarySegment<dim>>*/> LineSegmentContainerType;
        typedef std::map<int,MeshBoundarySegment<dim>> LineSegmentContainerType;
//        typedef std::tuple<VectorDim,size_t,double> SnapReturnType;
        
        //        /**********************************************************************/
        //        const LineSegmentContainerType& segments() const
        //        {
        //            return *this;
        //        }
        
    private:
        
        
        /**********************************************************************/
        void emplaceFromIntersection(const MeshBoundarySegment<dim>& s1,
                                     const MeshBoundarySegment<dim>& s2)
        {
            if(s1.face->sID==s2.face->sID)
            {// Lines on the same face
                SegmentSegmentDistance<dim> ssd(s1.P0,s1.P1,
                                                s2.P0,s2.P1);
                const auto iSeg=ssd.intersectionSegment();
                
                
                if(this->find(s1.face->sID)==this->end())
                {// no MeshBoundarySegment exists on current face
                    const auto iSeg=ssd.intersectionSegment();
                    if(iSeg.size()==1)
                    {
                        this->emplace(s1.face->sID,MeshBoundarySegment<dim>(std::get<0>(iSeg[0]),std::get<0>(iSeg[0]),s1.face));
                    }
                    else if(iSeg.size()==2)
                    {
                        this->emplace(s1.face->sID,MeshBoundarySegment<dim>(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]),s1.face));
                    }
                    else
                    {// do nothing
                        
                    }
                }
                else
                {// another MeshBoundarySegment exists on current face
                    assert(false && "FINISH HERE");
                }
                
            }
            
        }
        
    public:
        
        /**********************************************************************/
        BoundingLineSegments()
        {/* Empty constructor
          */
            
        }
        
        /**********************************************************************/
        BoundingLineSegments(const MeshPlane<dim>& gp)
        {/* Empty constructor
          */
            for(const auto& seg : gp.meshIntersections)
            {
                this->emplace(seg.face->sID,seg);
            }
        }

        /**********************************************************************/
        BoundingLineSegments(const BoundingLineSegments<dim>& bls1,
                             const BoundingLineSegments<dim>& bls2)
        {/*! Constructs the BoundingLineSegments from the intersection of two BoundingLineSegments
          */
            for (const auto& s1pair : bls1)
            {
                for (const auto& s2pair : bls2)
                {
                    
                    const MeshBoundarySegment<dim>& s1(s1pair.second);
                    const MeshBoundarySegment<dim>& s2(s2pair.second);
                    
                    emplaceFromIntersection(s1,s2);

                }
            }
            
            assert(this->size()<=std::max(bls1.size(),bls2.size()) && "intersections cannot be more than original segments");
        }
        
        /**********************************************************************/
        void updateWithMeshPlane(const MeshPlane<dim>& gp)
        {
            if(this->size())
            {
                BoundingLineSegments<dim> temp(*this,BoundingLineSegments(gp));
                this->swap(temp);
            }
            else
            {
                BoundingLineSegments<dim> temp(gp);
                this->swap(temp);
            }
        }

        /**********************************************************************/
        std::set<const MeshBoundarySegment<dim>*> containingSegments(const VectorDim& P) const
        {
            std::set<const MeshBoundarySegment<dim>*> temp;
            
            for(const auto& pair : *this)
            {
                if(pair.second.contains(P))
                {
                    temp.insert(&pair.second);
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
        VectorDim boundaryNormal(const VectorDim& P) const
        {
            VectorDim temp(VectorDim::Zero());
            for(const auto& seg : containingSegments(P))
            {
                temp+=seg->face->outNormal();
            }
            const double tempNorm(temp.norm());
            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const BoundingLineSegments<dim>& bls)
        {
            for(const auto& pair : bls)
            {
                os<<std::setprecision(15)<<std::scientific<<pair.second.P0.transpose()<<" "<<pair.second.P1.transpose()<<std::endl;
            }
            return os;
        }
        
    };
    
}
#endif


//        /**********************************************************************/
//        std::pair<bool,size_t> contains(const VectorDim& P) const
//        {
//
//            //            std::cout<<"Point="<<P.transpose()<<std::endl;
//            std::pair<bool,size_t> temp=std::make_pair(false,0);
//            for(size_t segID=0;segID<this->size();++segID)
//                //            for(const auto& vertexPair : *this)
//            {
//                const auto& vertexPair(this->operator[](segID));
//                const VectorDim segm(vertexPair.P1-vertexPair.P0);
//                const double segmNorm2(segm.squaredNorm());
//                if(segmNorm2>FLT_EPSILON)
//                {
//                    double u((P-vertexPair.P0).dot(segm)/segmNorm2);
//                    if(u<0.0)
//                    {
//                        u=0.0;
//                    }
//                    if(u>1.0)
//                    {
//                        u=1.0;
//                    }
//                    const VectorDim x(vertexPair.P0+u*segm);
//                    //                    std::cout<<vertexPair.P0.transpose()<<" "<<vertexPair.P1.transpose()<<"      case A: P-x="<<(P-x).squaredNorm()<<std::endl;
//                    if((P-x).squaredNorm()<FLT_EPSILON)
//                    {
//                        temp=std::make_pair(true,segID);
//                        break;
//                    }
//                }
//                else
//                {
//                    const VectorDim x(0.5*(vertexPair.P0+vertexPair.P1));
//                    //                    std::cout<<vertexPair.P0.transpose()<<" "<<vertexPair.P1.transpose()<<"      case B: P-x="<<(P-x).squaredNorm()<<std::endl;
//                    if((P-x).squaredNorm()<FLT_EPSILON)
//                    {
//                        temp=std::make_pair(true,segID);
//                        break;
//                    }
//                }
//            }
//            return temp;
//        }

//        /**********************************************************************/
//        bool contains(const VectorDim& P,const size_t& k) const
//        {
//            bool temp(false);
//            if(k<this->size())
//            {
//                const auto& vertexPair(this->operator[](k));
//                const VectorDim segm(vertexPair.P1-vertexPair.P0);
//                const double segmNorm2(segm.squaredNorm());
//                if(segmNorm2>FLT_EPSILON)
//                {
//                    double u((P-vertexPair.P0).dot(segm)/segmNorm2);
//                    if(u<0.0)
//                    {
//                        u=0.0;
//                    }
//                    if(u>1.0)
//                    {
//                        u=1.0;
//                    }
//                    const VectorDim x(vertexPair.P0+u*segm);
//                    if((P-x).squaredNorm()<FLT_EPSILON)
//                    {
//                        temp=true;
//                    }
//                }
//                else
//                {
//                    const VectorDim x(0.5*(vertexPair.P1+vertexPair.P0));
//                    if((P-x).squaredNorm()<FLT_EPSILON)
//                    {
//                        temp=true;
//                    }
//                }
//            }
//            return temp;
//        }

//        /**********************************************************************/
//        SnapReturnType snap(const VectorDim& P) const __attribute__ ((deprecated))
//        {
//
//            assert(this->size() && "CANNOT SNAP TO EMPTY BoundingLineSegments");
//
//            std::map<double,std::pair<VectorDim,size_t>,std::less<double>/*,Eigen::aligned_allocator<std::pair<double,std::pair<VectorDim,size_t>>>*/> snapMap;
//
//            //            for(const auto& vertexPair : *this)
//            for(size_t k=0;k<this->size();++k)
//            {
//                const MeshBoundarySegment<dim>& vertexPair(this->operator[](k));
//                const VectorDim segm(vertexPair.P1-vertexPair.P0);
//                const double segmNorm2(segm.squaredNorm());
//                if(segmNorm2>FLT_EPSILON)
//                {
//                    double u((P-vertexPair.P0).dot(segm)/segmNorm2);
//                    if(u<0.0)
//                    {
//                        u=0.0;
//                    }
//                    if(u>1.0)
//                    {
//                        u=1.0;
//                    }
//                    const VectorDim x(vertexPair.P0+u*segm);
//                    //snapMap.emplace((P-x).squaredNorm(),x);
//                    snapMap.emplace(std::piecewise_construct,
//                                    std::make_tuple((P-x).squaredNorm()),
//                                    std::make_tuple(x,k)
//                                    );
//
//                }
//                else
//                {
//                    const VectorDim x(0.5*(vertexPair.P1+vertexPair.P0));
//                    //                    snapMap.emplace((P-x).squaredNorm(),x);
//                    snapMap.emplace(std::piecewise_construct,
//                                    std::make_tuple((P-x).squaredNorm()),
//                                    std::make_tuple(x,k)
//                                    );
//
//                }
//            }
//
//            assert(contains(snapMap.begin()->second.first).first && "BoundingLineSegments does not contains snapped point.");
//
//            //            return snapMap.begin()->second;
//            return std::make_tuple(snapMap.begin()->second.first,snapMap.begin()->second.second,snapMap.begin()->first);
//        }

//        /**********************************************************************/
//        std::pair<double,VectorDim> snapToVertex(const VectorDim& P) const __attribute__ ((deprecated))
//        {
//
//            assert(this->size() && "CANNOT SNAP TO EMPTY BoundingLineSegments");
//
//            std::map<double,VectorDim,std::less<double>,Eigen::aligned_allocator<std::pair<double,VectorDim>>> snapMap;
//
//            for(const auto& vertexPair : *this)
//            {
//                snapMap.emplace((P-vertexPair.P0).norm(),vertexPair.P0);
//            }
//
//            return std::make_pair(snapMap.begin()->first,snapMap.begin()->second);
//
//        }

//        /**********************************************************************/
//        static void emplaceUnique(LineSegmentContainerType& temp,const VectorDim& P1,const VectorDim& P2)
//        {
//            bool isUnique=true;
//
//            std::cout<<"emplaceUnique "<<P1.transpose()<<" ### "<<P2.transpose()<<std::endl;
//
//            if((P1-P2).squaredNorm()>FLT_EPSILON)
//            {// P1-P2 is not a degenerate line
//                for(int p=0;p<temp.size();++p)
//                {
//                    VectorDim& E1(temp[p].first);
//                    VectorDim& E2(temp[p].second);
//                    if((E1-E2).squaredNorm()>FLT_EPSILON)
//                    {// existing pair is not a degenerate line
////                        isUnique*=((E1-P1).squaredNorm()>FLT_EPSILON && (E2-P2).squaredNorm()>FLT_EPSILON);
//                        isUnique*=((E1-P1).squaredNorm()>FLT_EPSILON || (E2-P2).squaredNorm()>FLT_EPSILON);
//                        isUnique*=((E1-P2).squaredNorm()>FLT_EPSILON || (E2-P1).squaredNorm()>FLT_EPSILON);
//                        std::cout<<"    A existing "<<E1.transpose()<<" ### "<<E2.transpose()<<" ### "<<isUnique<<std::endl;
//                        if(!isUnique)
//                        {
//                            break;
//                        }
//                    }
//                    else
//                    {// existing pair is a degenerate line
//                        const VectorDim x=0.5*(E1+E2);
//                        isUnique*=((P1-x).squaredNorm()>FLT_EPSILON && (P2-x).squaredNorm()>FLT_EPSILON);
//                        std::cout<<"    B existing "<<E1.transpose()<<" ### "<<E2.transpose()<<" ### "<<isUnique<<std::endl;
//                        if(!isUnique)
//                        {// keep non-defenerate line
//                            temp[p].first=P1;
//                            temp[p].second=P2;
//                            break;
//                        }
//                    }
//
//                }
//            }
//            else
//            {// P1-P2 is a degenerate line (a point)
//                const VectorDim x=0.5*(P1+P2);
//                for(const auto& pair : temp)
//                {
//                    isUnique*=((pair.first-x).squaredNorm()>FLT_EPSILON && (pair.second-x).squaredNorm()>FLT_EPSILON);
//                    std::cout<<"    C existing "<<pair.first.transpose()<<" ### "<<pair.second.transpose()<<" ### "<<isUnique<<std::endl;
//                    if(!isUnique)
//                    {
//                        break;
//                    }
//                }
//
//            }
//
//
//            if(isUnique)
//            {
//                temp.emplace_back(P1,P2);
//            }
//            else
//            {
////                std::cout<<"BoundingLineSegments: not unique"<<std::endl;
//            }
//        }

//        /**********************************************************************/
//        static void emplaceUnique(LineSegmentContainerType& temp,const VectorDim& P1,const VectorDim& P2)
//        {
//
//            //            std::cout<<"emplaceUnique "<<P1.transpose()<<" ### "<<P2.transpose()<<std::endl;
//
//            bool isUnique=true;
//
//
//            if((P1-P2).squaredNorm()>FLT_EPSILON)
//            {// P1-P2 is not a degenerate line
//                for(typename LineSegmentContainerType::iterator iter=temp.begin(); iter!=temp.end();)
//                {
//
//                    VectorDim& E1(iter->P0);
//                    VectorDim& E2(iter->P1);
//                    if((E1-E2).squaredNorm()>FLT_EPSILON)
//                    {// existing pair is not a degenerate line
//                        //                        isUnique*=((E1-P1).squaredNorm()>FLT_EPSILON && (E2-P2).squaredNorm()>FLT_EPSILON);
//                        isUnique*=((E1-P1).squaredNorm()>FLT_EPSILON || (E2-P2).squaredNorm()>FLT_EPSILON);
//                        isUnique*=((E1-P2).squaredNorm()>FLT_EPSILON || (E2-P1).squaredNorm()>FLT_EPSILON);
//                        //                        std::cout<<"    A existing "<<E1.transpose()<<" ### "<<E2.transpose()<<" ### "<<isUnique<<std::endl;
//                        if(!isUnique)
//                        {
//                            break;
//                        }
//                        else
//                        {
//                            iter++;
//                        }
//                    }
//                    else
//                    {// existing pair is a degenerate line
//                        const VectorDim x=0.5*(E1+E2);
//                        if((P1-x).squaredNorm()>FLT_EPSILON && (P2-x).squaredNorm()>FLT_EPSILON)
//                        {// distinct point
//                            //                            isUnique*=true;
//                            //                            std::cout<<"    B existing "<<E1.transpose()<<" ### "<<E2.transpose()<<" ### "<<isUnique<<std::endl;
//                            iter++;
//                        }
//                        else
//                        {// existing degenerate point is erased since P1-P2 will be appended later
//                            //                            isUnique*=true;
//                            //                            std::cout<<"    C existing "<<E1.transpose()<<" ### "<<E2.transpose()<<" ### "<<isUnique<<std::endl;
//                            iter=temp.erase(iter);
//                        }
//                    }
//                }
//            }
//            else
//            {// P1-P2 is a degenerate line (a point)
//                const VectorDim x=0.5*(P1+P2);
//                for(const auto& pair : temp)
//                {
//                    isUnique*=((pair.P0-x).squaredNorm()>FLT_EPSILON && (pair.P1-x).squaredNorm()>FLT_EPSILON);
//                    //                    std::cout<<"    D existing "<<pair.P0.transpose()<<" ### "<<pair.P1.transpose()<<" ### "<<isUnique<<std::endl;
//                    if(!isUnique)
//                    {
//                        break;
//                    }
//                }
//
//            }
//
//            if(isUnique)
//            {
//                temp.emplace_back(P1,P2);
//            }
//            else
//            {
//                //                std::cout<<"BoundingLineSegments: not unique"<<std::endl;
//            }
//
//        }


//        /**********************************************************************/
//        BoundingLineSegments(const BoundingLineSegments<dim>& bls1,
//                             const BoundingLineSegments<dim>& bls2)
//        {/* Constructor from intersectin of BoundingLineSegments
//          */
//            for (const auto& s1 : bls1)
//            {
//                for (const auto& s2 : bls2)
//                {
//                    //                    std::cout<<"Comparison"<<std::endl;
//                    //
//                    //                    SegmentSegmentIntersection<dim> ssi(s1.P0,s1.P1,
//                    //                                                        s2.P0,s2.P1);
//                    //
//                    //                    if(ssi.size)
//                    //                    {
//                    //                        std::cout<<ssi.x0.transpose()<<" -A- "<<ssi.x1.transpose()<<std::endl;
//                    //                                                this->emplace_back(ssi.x0,ssi.x1);
//                    //                    }
//
//                    SegmentSegmentDistance<dim> ssd(s1.P0,s1.P1,
//                                                    s2.P0,s2.P1);
//                    const auto iSeg=ssd.intersectionSegment();
//                    if(iSeg.size()==1)
//                    {
//                        //                        std::cout<<std::get<0>(iSeg[0]).transpose()<<" -B- "<<std::get<0>(iSeg[0]).transpose()<<std::endl;
//                        //                        this->emplace_back(std::get<0>(iSeg[0]),std::get<0>(iSeg[0]));
//                        emplaceUnique(*this,std::get<0>(iSeg[0]),std::get<0>(iSeg[0]));
//
//                    }
//                    else if(iSeg.size()==2)
//                    {
//                        //                        std::cout<<std::get<0>(iSeg[0]).transpose()<<" -C- "<<std::get<0>(iSeg[1]).transpose()<<std::endl;
//                        //                        this->emplace_back(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]));
//                        emplaceUnique(*this,std::get<0>(iSeg[0]),std::get<0>(iSeg[1]));
//                    }
//                    else
//                    {// do nothing
//
//                    }
//                }
//            }
//
//            assert(this->size()<=std::max(bls1.size(),bls2.size()) && "intersections cannot be more than original segments");
//            //            assert(this->size()<=bls2.size() && "intersections cannot be more than original segments");
//            //            std::cout<<"BoundingLineSegments "<<this<<", size="<<this->size()<<std::endl;
//
//        }

//        /**********************************************************************/
//        //        template <typename LoopType>
//        void updateWithMeshPlane(const MeshPlane<dim>& gp)
//        {
//
//            if(this->size())
//            {
//                *this=BoundingLineSegments(*this,BoundingLineSegments(gp));
//            }
//            else
//            {
//                *this=BoundingLineSegments(gp);
//            }
//
////            if(this->size())
////            {
////                LineSegmentContainerType temp;
////
////                for(const auto& oldPair : *this)
////                {// loop over current segments
////
////                    PlaneSegmentIntersection<dim> psi(gp.P,
////                                                      gp.unitNormal,
////                                                      oldPair.P0,
////                                                      oldPair.P1);
////
////                    //                    if(psi.size())
////                    if(psi.type==PlaneSegmentIntersection<dim>::INCIDENT || psi.type==PlaneSegmentIntersection<dim>::COINCIDENT)
////                    {// plane and current segment intersect
////                        for(const auto& seg : gp.meshIntersections)
////                        {
////
////                            SegmentSegmentDistance<dim> ssd(seg.P0,
////                                                            seg.P1,
////                                                            oldPair.P0,
////                                                            oldPair.P1);
////
////                            const auto iSeg=ssd.intersectionSegment();
////                            if(iSeg.size()==1)
////                            {
////                                emplaceUnique(temp,std::get<0>(iSeg[0]),std::get<0>(iSeg[0]));
////                            }
////                            else if(iSeg.size()==2)
////                            {
////                                emplaceUnique(temp,std::get<0>(iSeg[0]),std::get<0>(iSeg[1]));
////
////                            }
////                            else
////                            {// do nothing
////                                //                                std::cout<<"BoundingLineSegments: ssd no intersection"<<std::endl;
////
////                            }
////
////                        }
////                    }
////                    else
////                    {// plane does not contain segment, do nothing
////
////                    }
////                }
////
////                assert(temp.size()<=std::max(this->size(),gp.meshIntersections.size()) && "intersections cannot be more than original segments");
////                //assert(temp.size()<=gp.meshIntersections.size() && "intersections cannot be more than original segments");
////
////
////                this->swap(temp);
////            }
////            else
////            {// this is empty, so just copy segments from MashPlane
////                for(const auto& seg : gp.meshIntersections)
////                {
////                    this->push_back(seg);
////                }
////            }
//
//        }
