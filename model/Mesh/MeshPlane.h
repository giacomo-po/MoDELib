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
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        typedef std::pair<VectorDim,const Simplex<dim,1>* const> RootType;
        typedef std::deque<RootType> RootContainerType;
        
        
        /**********************************************************************/
        static Plane<dim> getPlaneBetweenRegions(const SimplicialMesh<dim>& mesh,
                                                 const int& rID1,
                                                 const int& rID2)
        {
            
            const auto& regionBnd(mesh.regionBoundary(rID1,rID2));
            
            VectorDim N(VectorDim::Zero());
            VectorDim P(VectorDim::Zero());
            size_t k=0;
            for(const auto& triangle : regionBnd)
            {
                
                const Simplex<dim,dim>& tet(**triangle->parents().begin());
                const size_t faceID=tet.childOrder(triangle->xID);
                const VectorDim n=tet.nda.col(faceID);

                if(n.dot(N)>=0.0)
                {
                    N+=n.normalized();
                }
                else
                {
                    N-=n.normalized();
                }
                
                const auto vertices= triangle->vertices();
                for(const auto& vertex :  vertices)
                {
                    P+=vertex->P0;
                    k++;
                }
            }
            Plane<dim> plane(P/k,N);
            
            // Check that all tringle vertices are contained by both GB planes
            for(const auto& triangle : regionBnd)
            {
                const auto vertices= triangle->vertices();
                for(const auto& vertex :  vertices)
                {
                    assert(plane.contains(vertex->P0) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
                }
            }
            
            return plane;
        }
        
        /**********************************************************************/
        static PlaneMeshIntersectionContainerType getBoundaryBetweenRegions(const SimplicialMesh<dim>& mesh,
                                                                            const int& rID1,
                                                                            const int& rID2)
        {
            assert(dim==3 && "ALGORITHM ONLY VALID IN dim=3");
            // Normal to plane
            
            // Get unsorted boundary edges
            std::set<const Simplex<dim,dim-2>*> ub=mesh.regionBoundary(rID1,rID2).unsortedBoundary();
            
            // Compute center of plane
            VectorDim c(VectorDim::Zero());
            for(const auto& edge : ub)
            {
                c+=edge->child(0).P0;
                c+=edge->child(1).P0;
            }
            c/=(2*ub.size()); //center of the plane
            
            // Rotation matrix to local reference system
            const VectorDim x3((*mesh.regionBoundary(rID1,rID2).begin())->outNormal(rID1).normalized());
            const VectorDim x1(((*ub.begin())->child(0).P0-c).normalized());
            Eigen::Matrix<double,dim,dim> R;
            R.col(0)=x1;
            R.col(1)=x3.cross(x1);
            R.col(2)=x3;
            
            // Sort the edges
            std::map<double,std::pair<const Simplex<dim,dim-2>* const,VectorDim>> sortedEdges;
            for(const auto& edge : ub)
            {
                const VectorDim x0=R.transpose()*(edge->child(0).P0-c);
                const double angle0=atan2(x0(1),x0(0));
                sortedEdges.emplace(angle0,std::make_pair(edge,edge->child(0).P0));
                
                const VectorDim x1=R.transpose()*(edge->child(1).P0-c);
                const double angle1=atan2(x1(1),x1(0));
                sortedEdges.emplace(angle1,std::make_pair(edge,edge->child(1).P0));
            }
            
            assert(sortedEdges.size()==ub.size());
            
            PlaneMeshIntersectionContainerType temp;
            for(const auto& pair : sortedEdges)
            {
                temp.push_back(pair.second);
            }
            
            return PlaneMeshIntersection<dim>::reducedPlaneMeshIntersection(temp);
        }
        
        /**********************************************************************/
        static PlaneMeshIntersectionContainerType getPlaneIntersection(const SimplicialMesh<dim>& mesh,
                                                                       const int& rID,
                                                                       const VectorDim& P0,
                                                                       const VectorDim& N)
        {
            
            const double nNorm(N.norm());
            assert(nNorm>FLT_EPSILON);
            const VectorDim n(N/nNorm);
            
            RootContainerType rootDeq;
            for(const auto& edge : mesh.template observer<1>())
            {
                if(   edge.second->isBoundarySimplex()
                   || edge.second->isRegionBoundarySimplex())
                {
                    if(edge.second->isInRegion(rID))
                    {
                        const VectorDim& v0(edge.second->child(0).P0);
                        const VectorDim& v1(edge.second->child(1).P0);
                        
                        // check intersection of v0->v1 with plane
                        // x=v0+u(v1-v0)
                        // (x-P0).n=0
                        // (v0+u(v1-v0)-P0).n=0
                        // u=(P0-v0).n/(v1-v0).n;
                        const double edgeNorm=(v1-v0).norm();
                        assert(edgeNorm>FLT_EPSILON && "mesh edge has zero norm.");
                        const double den=(v1-v0).dot(n);
                        const double num=(P0-v0).dot(n);
                        const double P0v0norm=(P0-v0).norm();
                        
                        const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
                        
                        if (fabs(den/edgeNorm)>FLT_EPSILON)
                        {
                            // edge intersects plane
                            const double u=num/den;
                            
                            if(fabs(u)<FLT_EPSILON)
                            {
                                rootDeq.emplace_back(v0,edge.second);
                            }
                            else if (u>=FLT_EPSILON && u<=1.0-FLT_EPSILON)
                            {
                                rootDeq.emplace_back((1.0-u)*v0 + u*v1,edge.second);
                            }
                            else if (fabs(1.0-u)<FLT_EPSILON)
                            {
                                rootDeq.emplace_back(v1,edge.second);
                            }
                            else
                            {// no roots
                                
                            }
                            
                        }
                        else
                        {
                            if (fabs(numCheck)>FLT_EPSILON)
                            {// edge is parallel to plane, no intersection
                                
                            }
                            else
                            {// edge is coplanar
                                rootDeq.emplace_back(v0,edge.second);
                                rootDeq.emplace_back(v1,edge.second);
                            }
                        }
                    }
                }
            }
//            std::cout<<"MeshPlane rootDeq: "<<rootDeq.size()<<std::endl;
            
            PlaneMeshIntersectionContainerType temp;

            if(rootDeq.size())
            {
            
                // compute center
                VectorDim c(VectorDim::Zero());
                for(const auto& pair : rootDeq)
                {
                    c+=pair.first;
                }
                c/=rootDeq.size(); //center of the plane
                
                // Check that points belong to a plane
                for(const auto& pair : rootDeq)
                {
                    assert(fabs(n.dot(c-pair.first))<FLT_EPSILON);
                }
                
                const VectorDim refDir(rootDeq[0].first-c);
                const double refDirNorm(refDir.norm());
                assert(refDirNorm>FLT_EPSILON);

                // Local rotation matrix
                Eigen::Matrix<double,dim,dim> R;
                R.col(0)=refDir/refDirNorm;
                R.col(2)=n;
                R.col(1)=R.col(2).cross(R.col(0));
                
                
                assert((R*R.transpose()-MatrixDim::Identity()).norm()<FLT_EPSILON);
                assert(fabs(R.determinant()-1.0)<FLT_EPSILON);
                
                std::map<double,std::pair<const Simplex<dim,dim-2>* const,VectorDim>> sortedEdges;
                
                for(const auto& pair : rootDeq)
                {
                    const VectorDim x0=R.transpose()*(pair.first-c);
                    const double angle0=atan2(x0(1),x0(0));
                    sortedEdges.emplace(angle0,std::make_pair(pair.second,pair.first));
                }
                
                for(const auto& pair : sortedEdges)
                {
                    temp.push_back(pair.second);
                }
                
            }
            
            
            return PlaneMeshIntersection<dim>::reducedPlaneMeshIntersection(temp);
        }
        
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
        meshIntersections(getPlaneIntersection(mesh,rID,p,this->unitNormal))
        //        /* init */ meshIntersections(PlaneMeshIntersection<dim>(mesh,this->P,this->unitNormal,rID))
        
        {/*!\param[in] mesh
          * \param[in] rID the region ID where the plane is defined
          * \param[in] p position of the plane
          * \param[in] n normal to the plane
          * Constructor for plane internal to a mesh region
          */
            if(meshIntersections.size()<3)
            {
                std::cout<<"p="<<p.transpose()<<std::endl;
                std::cout<<"n="<<n.transpose()<<std::endl;
                std::cout<<"meshIntersections.size()="<<meshIntersections.size()<<std::endl;
                assert(false && "meshIntersections FAILED");
            }

        }
        
        /**********************************************************************/
        MeshPlane(const SimplicialMesh<dim>& mesh,
                  const int& rID1,
                  const int& rID2) :
        /* init */ Plane<dim>(getPlaneBetweenRegions(mesh,rID1,rID2)),
        //        /* init */ regionID(rID),
        /* init */ regionIDs(rID1,rID2),
        /* init */ meshIntersections(getBoundaryBetweenRegions(mesh,rID1,rID2)) // WARNING: CALLING meshIntersections with rID1
        
        {/*!\param[in] mesh
          * \param[in] rID the region ID where the plane is defined
          * \param[in] p position of the plane
          * \param[in] n normal to the plane
          * Constructor for plane between two mesh regions
          */
            
            if(meshIntersections.size()<3)
            {
                std::cout<<"rID1="<<rID1<<std::endl;
                std::cout<<"rID2="<<rID2<<std::endl;
                std::cout<<"meshIntersections.size()="<<meshIntersections.size()<<std::endl;
                assert(false && "meshIntersections FAILED");
            }
            
        }
        
        //         typedef std::deque<std::pair<const Simplex<dim,dim-2>* const,VectorDim>> PlaneMeshIntersectionContainerType;
        
        //        PlaneMeshIntersectionContainerType getRegionBndPerimeter
        
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
