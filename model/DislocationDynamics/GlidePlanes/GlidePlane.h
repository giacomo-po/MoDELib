/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GLIDEPLANE_H
#define model_GLIDEPLANE_H

//#include <math.h>
#include <deque>
#include <chrono>
#include <map>
#include <assert.h>

#include <Eigen/Core>
//#include <Eigen/StdVector>

#include <model/Utilities/StaticID.h>
//#include <model/DislocationDynamics/DislocationNetworkTraits.h>
//#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
//#include <model/DislocationDynamics/DislocationSharedObjects.h>
//#include <model/DislocationDynamics/BVP/VirtualBoundarySlipSurface.h>
//#include <model/MPI/MPIcout.h>
//#include <model/Mesh/SimplexObserver.h>
//#include <model/LatticeMath/LatticeMath.h>

#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/MPI/MPIcout.h>

namespace model
{
    
    /*************************************************************/
    /*************************************************************/
    template <typename LoopType>
    class GlidePlane : public StaticID<GlidePlane<LoopType> >,
    /* base class   */ public LatticePlane,
    /* base class   */ private std::map<size_t,const LoopType* const>
    //	/* base class   */ private std::set<const VirtualBoundarySlipSurface<SegmentType>*>,
    //	/* base class   */ private GlidePlaneObserver<SegmentType>
    {
        
        
    public:
        
        constexpr static int dim=LoopType::dim;
        typedef typename TypeTraits<LoopType>::LinkType LinkType;
        typedef GlidePlane<LoopType> GlidePlaneType;
        typedef GlidePlaneObserver<LoopType> GlidePlaneObserverType;
//        typedef LatticeVector<dim> LatticeVectorType;
        //        typedef ReciprocalLatticeDi<dim> LatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename GlidePlaneObserverType::GlidePlaneKeyType GlidePlaneKeyType;
        //        typedef std::pair<VectorDim,VectorDim> segmentMeshCollisionPair;
        //        typedef std::deque<segmentMeshCollisionPair> SegmentMeshCollisionPairContainerType;
        typedef std::deque<std::pair<const Simplex<dim,dim-2>* const,VectorDim>> PlaneMeshIntersectionContainer;
        
        //        const LatticePlane glidePlane;
        
        static double meshIntersectionTol;
        const Grain<dim>& grain;
        const GlidePlaneKeyType glidePlaneKey;
        //! A container of the intersection lines between *this and the mesh boundary
        //const SegmentMeshCollisionPairContainerType segmentMeshCollisionPairContainer;
        const PlaneMeshIntersectionContainer meshIntersections;
        //                //! A container of the intersection lines between *this and the internal region mesh boundaries
        //                const SegmentMeshCollisionPairContainerType segmentRegionCollisionPairContainer;
        
        
        /**********************************************************************/
        GlidePlane(const Grain<dim>& grain_in,
                   const VectorDim& P,
                   const VectorDim& N) :
        /* init */ LatticePlane(grain_in.latticeVector(P),grain_in.reciprocalLatticeDirection(N)), // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
        /* init */ grain(grain_in),
        /* init */ glidePlaneKey(GlidePlaneObserverType::getGlidePlaneKey(grain,P,N)),
        /* init */ meshIntersections(planeMeshIntersection())
        //        /* init list */ segmentMeshCollisionPairContainer(getPlaneMeshIntersection(planeNormal*height,planeNormal)),
        //        /* init list */ segmentRegionCollisionPairContainer(getPlaneRegionIntersection(planeNormal*height,planeNormal))
        {
            model::cout<<"Creating GlidePlane "<<glidePlaneKey.transpose()<<std::endl;
            GlidePlaneObserverType::addGlidePlane(this);
        }
        
        /**********************************************************************/
        ~GlidePlane()
        {
            GlidePlaneObserverType::removeGlidePlane(this);
        }
        
        /**********************************************************************/
        const std::map<size_t,const LoopType* const>& loops() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addLoop(const LoopType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Adds pS to *this GLidePlane
          */
            const bool success=this->emplace(pL->sID,pL).second;
            assert( success && "COULD NOT INSERT LOOP POINTER IN GLIDE PLANE.");
        }
        
        /**********************************************************************/
        void removeLoop(const LoopType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Removes pS from *this GLidePlane
          */
            const int success=this->erase(pL->sID);
            assert(success==1 && "COULD NOT ERASE LOOP POINTER FROM GLIDE PLANE.");
        }
        
        /**********************************************************************/
        PlaneMeshIntersectionContainer planeMeshIntersection() const
        {
            std::cout<<"Computing mesh intersection"<<std::flush;
            const auto t0=std::chrono::system_clock::now();
            
            
            PlaneMeshIntersectionContainer temp;
            
            if (DislocationSharedObjects<dim>::use_boundary)
            {
                const VectorDim P0=this->P.cartesian();
                const VectorDim N=this->n.cartesian().normalized();
                
//                std::cout<<"P0="<<P0.transpose()<<std::endl;
//                std::cout<<"N="<<N.transpose()<<std::endl;
               
                // Find initial intersection
                std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> pEI=getFirstPlaneEdgeIntersection(P0,N);
//                std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> pEI;

                
                std::set<const Simplex<dim,dim-2>*> tested;
                
                tested.insert(pEI.first);
                for(const auto& x: pEI.second)
                {
                    temp.emplace_back(pEI.first,x);
                }
                
//                std::cout<<""
                
//                bool keepTesting=true;
                while(true)
                {
//                    keepTesting=false;
                    
                    std::map<double,std::pair<const Simplex<dim,dim-2>*,VectorDim>> intersectionMap;

                    
                    
                    for(const auto& sibling : temp.rbegin()->first->siblings())
                    {
//                        std::cout<<this->sID<<" "
//                        <<(*tested.rbegin())->child(0).xID<<" "<<(*tested.rbegin())->child(1).xID<<" "
//                        <<sibling->child(0).xID<<" "<<sibling->child(1).xID<<" "
//                        << sibling->isBoundarySimplex()<<" "
//                        << (tested.find(sibling)==tested.end())<<std::endl;
                        
                        
                        if(sibling->isBoundarySimplex())
                        {
                            if(tested.find(sibling)==tested.end())
                            {
                                tested.insert(sibling);
                                
                                std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> nextIntersections=planeEdgeIntersection(P0,N,*sibling);
                                
//                                std::cout<<"# intersections "<<nextIntersections.second.size()<<std::endl;
                                
                                
                                for(const auto& x: nextIntersections.second)
                                {
//                                    if((x-temp.rbegin()->second).norm()>FLT_EPSILON)
//                                    {
                                    
                                    //double w=1.0;
                                    
                                    const VectorDim dX=x-temp[temp.size()-1].second;
                                    VectorDim dXold(dX);
                                
                                    if(temp.size()>1)
                                    {
                                        dXold=temp[temp.size()-1].second-temp[temp.size()-2].second;
                                    }

                                    
                                    intersectionMap.emplace(dX.dot(dXold),std::make_pair(sibling,x));

                                    
                                        //temp.emplace_back(sibling,x);
//                                    }
                                    //                                        else
                                    //                                        {
                                    //                                            temp.rbegin()->first=sibling;
                                    //                                        }
                                }
                                
                                
//                                if(nextIntersections.second.size())
//                                {
//                                    for(const auto& x: nextIntersections.second)
//                                    {
//                                        if((x-temp.rbegin()->second).norm()>FLT_EPSILON)
//                                        {
//                                            temp.emplace_back(sibling,x);
//                                        }
////                                        else
////                                        {
////                                            temp.rbegin()->first=sibling;
////                                        }
//                                    }
//                                    
//                                    keepTesting=true;
//                                    break;
//                                }
                            }
                        }
                        
                        
                    }
                    
                    
                    if(intersectionMap.size())
                    {
                        if(intersectionMap.rbegin()->first<FLT_EPSILON)
                        {
                            temp.pop_back();
                            
                        }
                        //                    else
                        //                    {
                        //                        //temp.rbegin()->first=intersectionMap.rbegin()->second.first;
                        //                    }
                        
                        temp.emplace_back(intersectionMap.rbegin()->second.first,intersectionMap.rbegin()->second.second);
                    }
                    else
                    {
                        break;
                    }
                    
                    

                    
                }
                
                
                
                //                std::cout
                
            }
            
            std::cout<<" boundary perimeter has "<<temp.size()<<" points"<<std::flush;

            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            
            return temp;
        }
        
        /**********************************************************************/
        static std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> getFirstPlaneEdgeIntersection(const VectorDim& P0,
                                                                                                        const VectorDim& n)
        {
            std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> temp;
            for(const auto& edge : SimplexObserver<dim,dim-2>::simplices())
            {// loop over edges
                if(edge.second->isBoundarySimplex())
                {
                    temp=planeEdgeIntersection(P0,n,*edge.second);
                    if(temp.second.size())
                    { // first intersection found
                        break;
                    }
                }
            }
            
            assert(temp.second.size() && "plane does not intersect mesh.");
            
            return temp;
        }
        
        /**********************************************************************/
        static std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> planeEdgeIntersection(const VectorDim& P0,
                                                                                                const VectorDim& n,
                                                                                                const Simplex<dim,dim-2>& edge)
        {
            std::deque<VectorDim> temp;
            
            const VectorDim& v0(edge.child(0).P0);
            const VectorDim& v1(edge.child(1).P0);
            
            // check intersection of v0->v1 with plane
            // x=v0+u(v1-v0)
            // (x-P0).n=0
            // (v0+u(v1-v0)-P0).n=0
            // u=(P0-v0).n/(v1-v0).n;
            const double edgeNorm=(v1-v0).norm();
            assert(edgeNorm>FLT_EPSILON && "mesh edge has zero norm.");
            const double den=(v1-v0).dot(n);
            const double num=(P0-v0).dot(n);
            
            if (fabs(den/edgeNorm)>meshIntersectionTol)
            {
                // edge intersects plane
                const double u=num/den;
                if(u>=0.0 && u<=1.0)
                {
                    temp.emplace_back(v0 + u*(v1-v0));
                }
            }
            else
            {
                if (fabs(num)>meshIntersectionTol)
                {// edge is parallel to plane, no intersection
                    
                }
                else
                {// edge is coplanar
                    temp.emplace_back(v0);
                    temp.emplace_back(v1);
                }
            }
            
            return std::make_pair(&edge,temp);
        }
        
        
        /* friend T& operator << **********************************************/
        template <class T>
        friend T& operator << (T& os, const GlidePlaneType& gp)
        {
            size_t kk=0;
            for (const auto& x : gp.meshIntersections)
            {
                os<<gp.sID<< " "<<kk<<" "<< x.first->child(0).xID<< " "<< x.first->child(1).xID <<" "<<x.second.transpose()<<"\n";
                kk++;
            }

            return os;
        }
        
    };
    
    template <typename LoopType>
    double GlidePlane<LoopType>::meshIntersectionTol=100.0*DBL_EPSILON;
    
    /*************************************************************/
    /*************************************************************/
} // namespace model
#endif

