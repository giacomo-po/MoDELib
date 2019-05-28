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
#include <Eigen/StdVector>
#include <SimplicialMesh.h>
#include <PlanarMeshFace.h>
#include <Plane.h>
#include <StaticID.h>
//#include <PlaneMeshIntersection.h>
#include <PlaneSegmentIntersection.h>
#include <PlanePlaneIntersection.h>

//#include <MeshPlaneIntersection.h>

//#include <MeshPlaneFaceIntersection.h>
namespace model
{
    /*!\brief Class template representing a straight line at the intersection
     * of PlanarMeshFaces
     */
    template <int dim>
    struct MeshBoundarySegment : public LineSegment<dim>
    /*                        */,public std::set<const PlanarMeshFace<dim>*>
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef LineSegment<dim> LineSegmentType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFace<dim>*> FaceContainerType;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
//        const PlanarMeshFaceType* const face;
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const PlanarMeshFaceType* const face_in) :
        /* init */ LineSegmentType(p0,p1)
//        /* init */,face(face_in)
        {
            faces().insert(face_in);
            
        }
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const FaceContainerType& faces_in) :
        /* init */ LineSegmentType(p0,p1)
        /* init */,FaceContainerType(faces_in)
        {
        }
        
        FaceContainerType& faces()
        {
            return *this;
            
        }

        const FaceContainerType& faces() const
        {
            return *this;
            
        }
        
        /**********************************************************************/
        VectorDim boundaryNormal() const
        {
            VectorDim temp(VectorDim::Zero());
            for(const auto& face : faces())
            {
                temp+=face->outNormal();
            }
            const double tempNorm(temp.norm());
            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
        }

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const MeshBoundarySegment<dim>& seg)
        {
                os<<seg.faces().size()<<", "<<std::setprecision(15)<<std::scientific<<seg.P0.transpose()<<","<<seg.P1.transpose();
            return os;
        }
        
        
    };
    
    template <int dim>
    struct MeshPlane : public StaticID<MeshPlane<dim>>
    /*              */,public Plane<dim>
    {
        
        typedef std::array<long int,dim+3> MeshPlaneKeyType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        //        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        typedef std::pair<VectorDim,const Simplex<dim,1>* const> RootType;
        typedef std::deque<RootType> RootContainerType;
        typedef std::vector<MeshBoundarySegment<dim>, Eigen::aligned_allocator<MeshBoundarySegment<dim>>> MeshBoundaryContainerType;
        
        /**********************************************************************/
        static MeshBoundaryContainerType getFaceBoundary(const PlanarMeshFace<dim>& face)
        {
            MeshBoundaryContainerType temp;
            for(size_t k=0;k<face.convexHull().size();++k)
            {
                const size_t k1(k==face.convexHull().size()-1? 0 : k+1);
                temp.emplace_back(face.convexHull()[k]->P0,face.convexHull()[k1]->P0,&face);
            }
            return temp;
        }
        
//        static std::set<const PlanarMeshFace<dim>*> getFaces(const MeshBoundaryContainerType& meshBoundaryContainer)
//        {
//            std::set<const PlanarMeshFace<dim>*> temp;
//            for(const auto& seg : meshBoundaryContainer)
//            {
//                temp.insert(seg.face);
//            }
//            return temp;
//        }

        
        /**********************************************************************/
        static MeshBoundaryContainerType getMeshBoundarySegments(const SimplicialMesh<dim>& mesh,
                                                                 const int& rID,
                                                                 const Plane<dim>& plane)
        {/*!\todo MeshFace and this function use ConvexHull and therefore will yield a intersection
          * for non-convex faces. In order to support non-convex faces
          * we need to have a structure with face edges.
          */
            
            //
            //            std::cout<<"getMeshBoundarySegments start"<<std::endl;
            // Collect all intersection points
            //            Plane<dim> plane(p,n);
            MeshBoundaryContainerType temp;
            const MatrixDim R(plane.localRotationMatrix());
            
            for(const auto& face : mesh.region(rID)->faces())
            {
                //std::cout<<face.second->outNormal()<<std::endl;
                PlanePlaneIntersection<dim> ppi(plane,face.second->asPlane());
                
                if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                {// Collect all intersection points with lines defined by convexHull of face
                                        //std::cout<<"incident"<<std::endl;
                    std::vector<VectorDim> pts;
                    //std::cout<<"face.second->convexHull().size()="<<face.second->convexHull().size()<<std::endl;

                    for(size_t k=0;k<face.second->convexHull().size();++k)
                    {
                        const size_t k1(k==face.second->convexHull().size()-1? 0 : k+1);
                        PlaneSegmentIntersection<dim> psi(plane.P,plane.unitNormal,face.second->convexHull()[k]->P0,face.second->convexHull()[k1]->P0);
                        
                        if(psi.type==PlaneSegmentIntersection<dim>::INCIDENT)
                        {
                            pts.push_back(psi.x0);
                        }
                        else if(psi.type==PlaneSegmentIntersection<dim>::COINCIDENT)
                        {
                            pts.push_back(psi.x0);
                            pts.push_back(psi.x1);
                            break;
                        }
                    }

                    // Compute ConvexHull of intersection points
                    ConvexHull<2,VectorDim> hull;
                    for(const auto& pt : pts)
                    {
                        VectorDim x(R*(pt-plane.P));
//                        HullPoint<2,VectorDim> hp({x[0],x[1]},&pt);
//                        hull.push_back(hp);
                        hull.emplace(std::array<double,2>{x[0],x[1]},&pt);
                    }

                    const auto hullPts=hull.getPoints();
                    if(hullPts.size()==2)
                    {
                        //std::cout<<"a"<<std::endl;
                        temp.emplace_back(*hullPts[0].t,*hullPts[1].t,face.second.get());
                    }
                    else if(hullPts.size()>2)
                    {
                                                //std::cout<<"b"<<std::endl;
                        for(size_t k=0;k<hullPts.size();++k)
                        {
                            const size_t k1(k==hullPts.size()-1? 0 : k+1);
                            temp.emplace_back(*hullPts[k].t,*hullPts[k1].t,face.second.get());
                        }
                    }
                    else
                    {// don't push points
                                                //std::cout<<"c"<<std::endl;
                    }
                    
                }
                else if (ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                {
                                                            //std::cout<<"coincident"<<std::endl;
                    temp.clear();
                    for(size_t k=0;k<face.second->convexHull().size();++k)
                    {
                        const size_t k1(k==face.second->convexHull().size()-1? 0 : k+1);
                        temp.emplace_back(face.second->convexHull()[k]->P0,face.second->convexHull()[k1]->P0,face.second.get());
                    }
                    break;
                }
                
                
            }
            

            
            ConvexHull<2,MeshBoundarySegment<dim>> finalHull;
            //std::cout<<"unsorted hull"<<std::endl;
            for(const auto& pt : temp)
            {
                
                VectorDim x(R*(0.5*(pt.P0+pt.P1)-plane.P));
                finalHull.emplace(std::array<double,2>{x[0],x[1]},&pt);
                //std::cout<<pt.P0.transpose()<<","<<pt.P1.transpose()<<std::endl;
            }
            const auto hullPts=finalHull.getPoints();

            assert(hullPts.size()==temp.size());

            MeshBoundaryContainerType sortedTemp;
//            for(const auto& seg : hullPts)
//            {
//                sortedTemp.push_back(*seg.t);
//            }
                        //std::cout<<"sorted hull"<<std::endl;
            for(size_t k=0;k<hullPts.size();++k)
            {
                const auto& seg(*hullPts[k].t);
                if(k==0)
                {
                    const auto& seg1(*hullPts[1].t);

                    if((seg.P1-seg1.P0).norm()<FLT_EPSILON || (seg.P1-seg1.P1).norm()<FLT_EPSILON)
                    {// first segment is oriented in the right direction
                        sortedTemp.push_back(seg);
                    }
                    else
                    {
                        sortedTemp.emplace_back(seg.P1,seg.P0,seg.faces());
                    }
                }
                else
                {
                    if((sortedTemp.back().P1-seg.P0).norm()<FLT_EPSILON)
                    {
                        sortedTemp.push_back(seg);
                    }
                    else if((sortedTemp.back().P1-seg.P1).norm()<FLT_EPSILON)
                    {
                        sortedTemp.emplace_back(seg.P1,seg.P0,seg.faces());
                    }
                    else
                    {
                        //std::cout<<"k="<<k<<" of "<<hullPts.size()<<std::endl;
                        assert(false && "DISCONNECTED FACE BOUNDARY");
                    }
                }
                //std::cout<<sortedTemp.back().P0.transpose()<<","<<sortedTemp.back().P1.transpose()<<std::endl;

            }
            assert((sortedTemp.back().P1-sortedTemp.front().P0).norm()<FLT_EPSILON && "OPEN FACE BOUNDARY");

            assert(sortedTemp.size()==temp.size());

            
            return sortedTemp;
            
        }
        
        const std::pair<int,int> regionIDs;
        //        const PlaneMeshIntersectionContainerType meshIntersections;
        const MeshBoundaryContainerType meshIntersections;
//        const std::set<const PlanarMeshFace<dim>*> meshFaces;
        
        /**********************************************************************/
        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const int& rID,
                  const VectorDim& p,
                  const VectorDim& n) :
        /* init */ Plane<dim>(p,n)
        /* init */,regionIDs(rID,rID)
        /* init */,meshIntersections(getMeshBoundarySegments(mesh,rID,*this))
//        /* init */,meshFaces(getFaces(meshIntersections))
        {/*!\param[in] mesh
          * \param[in] rID the region ID where the plane is defined
          * \param[in] p position of the plane
          * \param[in] n normal to the plane
          * Constructor for plane internal to a mesh region
          */
            //            std::cout<<"Constructing MeshPlane "<<std::endl;
            if(meshIntersections.size()<3)
            {
                std::cout<<"p="<<p.transpose()<<std::endl;
                std::cout<<"n="<<n.transpose()<<std::endl;
                std::cout<<"meshIntersections.size()="<<meshIntersections.size()<<std::endl;
                assert(false && "meshIntersections FAILED");
            }
        }
        
        /**********************************************************************/
        MeshPlane(const PlanarMeshFace<dim>& face,
                  const int& rID1,
                  const int& rID2) :
        /* init */ Plane<dim>(face.asPlane()),
        //        /* init */ regionID(rID),
        /* init */ regionIDs(rID1,rID2),
        /* init */ meshIntersections(getFaceBoundary(face)) // WARNING: CALLING meshIntersections with rID1
        //        ,meshIntersections2(getMeshBoundarySegments(mesh,rID),*this)
        {
            
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const MeshPlane<dim>& gp)
        {
            for (const auto& x : gp.meshIntersections)
            {
                os<<gp.sID<< " "<<x;
            }
            return os;
        }
        
    };
    
}
#endif


//        /**********************************************************************/
//        static PlaneMeshIntersectionContainerType getBoundaryBetweenRegions(const SimplicialMesh<dim>& mesh,
//                                                                            const int& rID1,
//                                                                            const int& rID2) __attribute__ ((deprecated))
//        {
//            assert(dim==3 && "ALGORITHM ONLY VALID IN dim=3");
//            // Normal to plane
//
//            // Get unsorted boundary edges
//            std::set<const Simplex<dim,dim-2>*> ub=mesh.regionBoundary(rID1,rID2).unsortedBoundary();
//
//            // Compute center of plane
//            VectorDim c(VectorDim::Zero());
//            for(const auto& edge : ub)
//            {
//                c+=edge->child(0).P0;
//                c+=edge->child(1).P0;
//            }
//            c/=(2*ub.size()); //center of the plane
//
//            // Rotation matrix to local reference system
//            const VectorDim x3((*mesh.regionBoundary(rID1,rID2).simplices().begin())->outNormal(rID1).normalized());
//            const VectorDim x1(((*ub.begin())->child(0).P0-c).normalized());
//            Eigen::Matrix<double,dim,dim> R;
//            R.col(0)=x1;
//            R.col(1)=x3.cross(x1);
//            R.col(2)=x3;
//
//            // Sort the edges
//            std::map<double,std::pair<const Simplex<dim,dim-2>* const,VectorDim>> sortedEdges;
//            for(const auto& edge : ub)
//            {
//                const VectorDim x0=R.transpose()*(edge->child(0).P0-c);
//                const double angle0=atan2(x0(1),x0(0));
//                sortedEdges.emplace(angle0,std::make_pair(edge,edge->child(0).P0));
//
//                const VectorDim x1=R.transpose()*(edge->child(1).P0-c);
//                const double angle1=atan2(x1(1),x1(0));
//                sortedEdges.emplace(angle1,std::make_pair(edge,edge->child(1).P0));
//            }
//
//            assert(sortedEdges.size()==ub.size());
//
//            PlaneMeshIntersectionContainerType temp;
//            for(const auto& pair : sortedEdges)
//            {
//                temp.push_back(pair.second);
//            }
//
//            return PlaneMeshIntersection<dim>::reducedPlaneMeshIntersection(temp);
//        }


//        /**********************************************************************/
//        MeshPlane(const SimplicialMesh<dim>& mesh,
//                  const int& rID1,
//                  const int& rID2) :
//        /* init */ Plane<dim>(getPlaneBetweenRegions(mesh,rID1,rID2)),
//        //        /* init */ regionID(rID),
//        /* init */ regionIDs(rID1,rID2),
//        /* init */ meshIntersections(getBoundaryBetweenRegions(mesh,rID1,rID2)) // WARNING: CALLING meshIntersections with rID1
////        ,meshIntersections2(getMeshBoundarySegments(mesh,rID),*this)
//        {/*!\param[in] mesh
//          * \param[in] rID the region ID where the plane is defined
//          * \param[in] p position of the plane
//          * \param[in] n normal to the plane
//          * Constructor for plane between two mesh regions
//          */
//
//            if(meshIntersections.size()<3)
//            {
//                std::cout<<"rID1="<<rID1<<std::endl;
//                std::cout<<"rID2="<<rID2<<std::endl;
//                std::cout<<"meshIntersections.size()="<<meshIntersections.size()<<std::endl;
//                assert(false && "meshIntersections FAILED");
//            }
//
//        }

//         typedef std::deque<std::pair<const Simplex<dim,dim-2>* const,VectorDim>> PlaneMeshIntersectionContainerType;

//        PlaneMeshIntersectionContainerType getRegionBndPerimeter


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


//        /**********************************************************************/
//        static Plane<dim> getPlaneBetweenRegions(const SimplicialMesh<dim>& mesh,
//                                                 const int& rID1,
//                                                 const int& rID2) __attribute__ ((deprecated))
//        {
//
//            const auto& regionBnd(mesh.regionBoundary(rID1,rID2));
//
//            VectorDim N(VectorDim::Zero());
//            VectorDim P(VectorDim::Zero());
//            size_t k=0;
//            for(const auto& triangle : regionBnd.simplices())
//            {
//
//                const Simplex<dim,dim>& tet(*triangle->parents().begin()->second);
//                const size_t faceID=tet.childOrder(triangle->xID);
//                const VectorDim n=tet.nda.col(faceID);
//
//                if(n.dot(N)>=0.0)
//                {
//                    N+=n.normalized();
//                }
//                else
//                {
//                    N-=n.normalized();
//                }
//
//                const auto vertices= triangle->vertices();
//                for(const auto& vertex :  vertices)
//                {
//                    P+=vertex->P0;
//                    k++;
//                }
//            }
//            Plane<dim> plane(P/k,N);
//
//            // Check that all tringle vertices are contained by both GB planes
//            for(const auto& triangle : regionBnd.simplices())
//            {
//                const auto vertices= triangle->vertices();
//                for(const auto& vertex :  vertices)
//                {
//                    assert(plane.contains(vertex->P0) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
//                }
//            }
//
//            return plane;
//        }

//        /**********************************************************************/
//        static PlaneMeshIntersectionContainerType getPlaneIntersection(const SimplicialMesh<dim>& mesh,
//                                                                       const int& rID,
//                                                                       const VectorDim& P0,
//                                                                       const VectorDim& N) __attribute__ ((deprecated))
//        {
//
//            const double nNorm(N.norm());
//            assert(nNorm>FLT_EPSILON);
//            const VectorDim n(N/nNorm);
//
//            RootContainerType rootDeq;
//            for(const auto& edge : mesh.template observer<1>())
//            {
//                if(   edge.second->isBoundarySimplex()
//                   || edge.second->isRegionBoundarySimplex())
//                {
//                    if(edge.second->isInRegion(rID))
//                    {
//                        const VectorDim& v0(edge.second->child(0).P0);
//                        const VectorDim& v1(edge.second->child(1).P0);
//
//                        // check intersection of v0->v1 with plane
//                        // x=v0+u(v1-v0)
//                        // (x-P0).n=0
//                        // (v0+u(v1-v0)-P0).n=0
//                        // u=(P0-v0).n/(v1-v0).n;
//                        const double edgeNorm=(v1-v0).norm();
//                        assert(edgeNorm>FLT_EPSILON && "mesh edge has zero norm.");
//                        const double den=(v1-v0).dot(n);
//                        const double num=(P0-v0).dot(n);
//                        const double P0v0norm=(P0-v0).norm();
//
//                        const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
//
//                        if (fabs(den/edgeNorm)>FLT_EPSILON)
//                        {
//                            // edge intersects plane
//                            const double u=num/den;
//
//                            if(fabs(u)<FLT_EPSILON)
//                            {
//                                rootDeq.emplace_back(v0,edge.second);
//                            }
//                            else if (u>=FLT_EPSILON && u<=1.0-FLT_EPSILON)
//                            {
//                                rootDeq.emplace_back((1.0-u)*v0 + u*v1,edge.second);
//                            }
//                            else if (fabs(1.0-u)<FLT_EPSILON)
//                            {
//                                rootDeq.emplace_back(v1,edge.second);
//                            }
//                            else
//                            {// no roots
//
//                            }
//
//                        }
//                        else
//                        {
//                            if (fabs(numCheck)>FLT_EPSILON)
//                            {// edge is parallel to plane, no intersection
//
//                            }
//                            else
//                            {// edge is coplanar
//                                rootDeq.emplace_back(v0,edge.second);
//                                rootDeq.emplace_back(v1,edge.second);
//                            }
//                        }
//                    }
//                }
//            }
////            std::cout<<"MeshPlane rootDeq: "<<rootDeq.size()<<std::endl;
//
//            PlaneMeshIntersectionContainerType temp;
//
//            if(rootDeq.size())
//            {
//
//                // compute center
//                VectorDim c(VectorDim::Zero());
//                for(const auto& pair : rootDeq)
//                {
//                    c+=pair.first;
//                }
//                c/=rootDeq.size(); //center of the plane
//
//                // Check that points belong to a plane
//                for(const auto& pair : rootDeq)
//                {
//                    assert(fabs(n.dot(c-pair.first))<FLT_EPSILON);
//                }
//
//                const VectorDim refDir(rootDeq[0].first-c);
//                const double refDirNorm(refDir.norm());
//                assert(refDirNorm>FLT_EPSILON);
//
//                // Local rotation matrix
//                Eigen::Matrix<double,dim,dim> R;
//                R.col(0)=refDir/refDirNorm;
//                R.col(2)=n;
//                R.col(1)=R.col(2).cross(R.col(0));
//
//
//                assert((R*R.transpose()-MatrixDim::Identity()).norm()<FLT_EPSILON);
//                assert(fabs(R.determinant()-1.0)<FLT_EPSILON);
//
//                std::map<double,std::pair<const Simplex<dim,dim-2>* const,VectorDim>> sortedEdges;
//
//                for(const auto& pair : rootDeq)
//                {
//                    const VectorDim x0=R.transpose()*(pair.first-c);
//                    const double angle0=atan2(x0(1),x0(0));
//                    sortedEdges.emplace(angle0,std::make_pair(pair.second,pair.first));
//                }
//
//                for(const auto& pair : sortedEdges)
//                {
//                    temp.push_back(pair.second);
//                }
//
//            }
//
//
//            return PlaneMeshIntersection<dim>::reducedPlaneMeshIntersection(temp);
//        }

//        static MeshBoundaryContainerType sortMeshBoundarySegments(const MeshBoundaryContainerType& vin)
//        {/*!\todo Remove this function when getMeshBoundarySegments uses FaceEdge to sort intersections
//          */
//
//            MeshBoundaryContainerType vout;
//            if(vin.size())
//            {
//                for(const auto& seg : vin)
//                {
//
//                }
//            }
//            return vout;
//        }
