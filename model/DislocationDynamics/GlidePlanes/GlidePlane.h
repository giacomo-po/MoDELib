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
#include <model/Mesh/PlaneMeshIntersection.h>

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
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;

        //        typedef std::pair<VectorDim,VectorDim> segmentMeshCollisionPair;
        //        typedef std::deque<segmentMeshCollisionPair> SegmentMeshCollisionPairContainerType;
        const Grain<dim>& grain;
        const GlidePlaneKeyType glidePlaneKey;
        //! A container of the intersection lines between *this and the mesh boundary
        //const SegmentMeshCollisionPairContainerType segmentMeshCollisionPairContainer;
        const PlaneMeshIntersectionContainerType meshIntersections;
        //                //! A container of the intersection lines between *this and the internal region mesh boundaries
        //                const SegmentMeshCollisionPairContainerType segmentRegionCollisionPairContainer;
        
        
        /**********************************************************************/
        GlidePlane(const Grain<dim>& grain_in,
                   const VectorDim& P,
                   const VectorDim& N) :
        /* init */ LatticePlane(grain_in.latticeVector(P),grain_in.reciprocalLatticeDirection(N)), // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
        /* init */ grain(grain_in),
        /* init */ glidePlaneKey(GlidePlaneObserverType::getGlidePlaneKey(grain,P,N)),
        /* init */ meshIntersections(PlaneMeshIntersection<dim>::planeMeshIntersection(this->P.cartesian(),this->n.cartesian().normalized(),grain.grainID))
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
    
//    template <typename LoopType>
//    double GlidePlane<LoopType>::meshIntersectionTol=FLT_EPSILON;
    
} // namespace model
#endif



//                for(const auto& x: pEI.second)
//                {
//                    temp.emplace_back(pEI.first,x);
//                }

//                std::cout<<""

//                //                bool keepTesting=true;
//                while(true)
//                {
//                    //                    keepTesting=false;
//
////                    std::map<double,std::pair<const Simplex<dim,dim-2>*,VectorDim>> intersectionMap;
//                    std::map<std::tuple<int,double,double>,std::pair<const Simplex<dim,dim-2>*,VectorDim>> intersectionMap;
//
//
//                    //                    std::cout<<"gID="<<grain.grainID<<std::endl;
//                    //                    std::cout<<"regionSiblings.size="<<temp.rbegin()->first->regionSiblings(grain.grainID)<<std::endl;
//
//                    //                 for(const auto& sibling : temp.rbegin()->first->siblings())
//                    //                        for(const auto& sibling : temp.rbegin()->first->siblingsInRegion(grain.grainID))
//
//
//
//                    //std::set<const SimplexType*,SimplexCompare<dim,order> >
//
//
//                    SiblingsContainerType siblings;
//                    for(const auto& parent : temp.rbegin()->first->parents())
//                    {
//                        if(  (parent->isBoundarySimplex() || parent->isRegionBoundarySimplex())
//                           && parent->isInRegion(grain.grainID))
//                        {
//                            for(int c=0; c<ParentSimplexType::nFaces;++c)
//                            {
//                                siblings.insert(&parent->child(c));
//                            }
//                        }
//                    }
//
//
////                    for(const auto& sibling : temp.rbegin()->first->boundarySiblingsInRegion(grain.grainID))
//                        for(const auto& sibling : siblings)
//
//                    {
//                        std::cout<<this->sID<<" "
//                        <<temp.rbegin()->first->child(0).xID<<" "<<temp.rbegin()->first->child(1).xID<<" "
//                        <<sibling->child(0).xID<<" "<<sibling->child(1).xID<<" "
//                        << sibling->isBoundarySimplex()<<" "
//                        << sibling->isRegionBoundarySimplex()<<" "
//                        << (tested.find(sibling)==tested.end())<<std::endl;
//
//
////                        if(//sibling->isBoundarySimplex() ||
////                           sibling->isRegionBoundarySimplex())
////                        {
//                            if(tested.find(sibling)==tested.end())
//                            {
//
//                                std::pair<const Simplex<dim,dim-2>*,std::deque<VectorDim>> nextIntersections=planeEdgeIntersection(P0,N,*sibling);
//
//                                std::cout<<"# intersections "<<nextIntersections.second.size()<<std::endl;
//
//
//                                if(nextIntersections.second.size()==0)
//                                {
//                                    tested.insert(sibling);
//                                }
//
//                                for(const auto& x: nextIntersections.second)
//                                {
//                                    //                                    if((x-temp.rbegin()->second).norm()>FLT_EPSILON)
//                                    //                                    {
//
//                                    //double w=1.0;
//
////                                    const VectorDim dX=x-temp[temp.size()-1].second;
////                                    VectorDim dXold(dX);
//
//                                    const double b=(x-temp[temp.size()-1].second).norm();
//                                    double a=b;
//
//
//                                    if(temp.size()>1)
//                                    {
////                                        dXold=temp[temp.size()-1].second-temp[temp.size()-2].second;
//                                        a=(x-temp[temp.size()-2].second).norm();
//                                       // if (a<FLT_EPSILON)
//                                       // {
//
//                                      //  }
//                                    }
//
//
//                                    //intersectionMap.emplace(dX.dot(dXold),std::make_pair(sibling,x));
////                                    intersectionMap.emplace(std::make_pair(a+b,a),std::make_pair(sibling,x));
//                                    intersectionMap.emplace(std::make_tuple(nextIntersections.second.size(),a+b,a),std::make_pair(sibling,x));
//
//
//                                    //                                    intersectionMap.emplace((x-temp[temp.size()-1]).squaredNorm(),std::make_pair(sibling,x));
//
//
//                                    //temp.emplace_back(sibling,x);
//                                    //                                    }
//                                    //                                        else
//                                    //                                        {
//                                    //                                            temp.rbegin()->first=sibling;
//                                    //                                        }
//                                }
//
//
//                                //                                if(nextIntersections.second.size())
//                                //                                {
//                                //                                    for(const auto& x: nextIntersections.second)
//                                //                                    {
//                                //                                        if((x-temp.rbegin()->second).norm()>FLT_EPSILON)
//                                //                                        {
//                                //                                            temp.emplace_back(sibling,x);
//                                //                                        }
//                                ////                                        else
//                                ////                                        {
//                                ////                                            temp.rbegin()->first=sibling;
//                                ////                                        }
//                                //                                    }
//                                //
//                                //                                    keepTesting=true;
//                                //                                    break;
//                                //                                }
//                            }
//                        //}
//
//
//                    }
//
//
//                    if(intersectionMap.size())
//                    {
//
//                        tested.insert(intersectionMap.rbegin()->second.first);
////                                                if((intersectionMap.rbegin()->first.first-intersectionMap.rbegin()->first.second)<FLT_EPSILON)
////                                                {
////                                                    temp.pop_back();
////
////                                                }
//                        if((std::get<1>(intersectionMap.rbegin()->first)-std::get<2>(intersectionMap.rbegin()->first))<FLT_EPSILON)
//                        {
//                            temp.pop_back();
//
//                        }
//                        temp.emplace_back(intersectionMap.rbegin()->second.first,intersectionMap.rbegin()->second.second);
//
//
//                        //                        tested.insert(intersectionMap.begin()->second.first);
//                        //
//                        //                        if(intersectionMap.begin()->first<FLT_EPSILON)
//                        //                        {
//                        //                            temp.pop_back();
//                        //
//                        //                        }
//                        //
//                        //
//                        //                        temp.emplace_back(intersectionMap.begin()->second.first,intersectionMap.begin()->second.second);
//
//                        //                    else
//                        //                    {
//                        //                        //temp.rbegin()->first=intersectionMap.rbegin()->second.first;
//                        //                    }
//                    }
//                    else
//                    {
//                        break;
//                    }
//
//
//
//
//                }



//                std::cout

