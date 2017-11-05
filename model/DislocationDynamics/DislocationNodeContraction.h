/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_H_
#define model_DislocationNodeContraction_H_

#include <memory>
#include <tuple>
#include <Eigen/Dense>
#include <model/Utilities/TypeTraits.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/Geometry/LineSegment.h>

#ifndef NDEBUG
#define VerboseNodeContraction(N,x) if(verboseNodeContraction>=N){model::cout<<x;}
#else
#define VerboseNodeContraction(N,x) 
#endif

namespace model
{
    template <typename DislocationNetworkType>
    struct DislocationNodeContraction
    {
        
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        static int verboseNodeContraction;

        
        DislocationNetworkType& DN;
        
        /**********************************************************************/
        DislocationNodeContraction(DislocationNetworkType& DN_in) : DN(DN_in)
        {
            
        }
        
        /**********************************************************************/
        bool contractYoungest(std::shared_ptr<NodeType> nA,
                              std::shared_ptr<NodeType> nB)
        {
            return nA->sID<nB->sID? DN.contractSecond(nA,nB) : DN.contractSecond(nB,nA);
            
        }
        
        
        //        /**********************************************************************/
        //        static bool isMovable(const std::shared_ptr<NodeType>& nA,
        //                              const std::shared_ptr<NodeType>& nB,
        //                              const VectorDim& X) // THIS SHOULD BE A DislocationNode::isMovableTo(X)
        //        {
        //            bool nAisMovable=true;
        //            for(const auto& pair : nA->neighbors())
        //            {
        //                //                            std::cout<<"neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first<<std::endl;
        //                if(pair.first!=nB->sID && std::get<1>(pair.second)->isBoundarySegment())
        //                {// boundary segments other than A->B must remain boundary
        //                    nAisMovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X)).first;
        //                }
        //
        //                if(pair.first!=nB->sID && std::get<1>(pair.second)->isGrainBoundarySegment())
        //                {// boundary segments other than A->B must remain boundary
        //                    for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
        //                    {
        //                        //                                        std::cout<<"grainBoundary neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nB->get_P())<<std::endl;
        //                        nAisMovable*=gb->glidePlanes().begin()->second.contains(X);
        //                    }
        //                }
        //            }
        //            return nAisMovable;
        //        }
        
        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            
            VerboseNodeContraction(1,"DislocationNodeContraction::contract "<<nA->sID<<" "<<nB->sID<<std::endl;);
            
            if(nA->isGlissile() && nB->isGlissile())
            {// both nodes are glissile
                
                if(nA->isOnBoundingBox() || nB->isOnBoundingBox())
                {// either one of the nodes is a boundary node. Therefore the contraction point must be a boundary node
                
                    BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
                    VerboseNodeContraction(1,"temp.size="<<temp.size()<<std::endl;);
                    if(temp.size())
                    {// a necessay condition is that the bounding boxes intersect
                        if(nA->isOnBoundingBox() && nB->isOnBoundingBox())
                        {// both nodes on bounding boxes. Intersect bounding boxes'
                            const bool nAisMovable=nA->isMovableTo(nB->get_P());
                            const bool nBisMovable=nB->isMovableTo(nA->get_P());
                            
                            if(nAisMovable && nBisMovable)
                            {
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1a"<<std::endl;);
                                return contractYoungest(nA,nB);
                            }
                            else if(nAisMovable && !nBisMovable)
                            {
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1b"<<std::endl;);
                                return DN.contractSecond(nB,nA);
                            }
                            else if(!nAisMovable && nBisMovable)
                            {
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1c"<<std::endl;);
                                return DN.contractSecond(nA,nB);
                            }
                            else
                            {
                                VerboseNodeContraction(1,"DislocationNodeContraction case 1d"<<std::endl;);
                                return false;
                            }
                        }
                        else if(nA->isOnBoundingBox() && !nB->isOnBoundingBox())
                        {// a on box, b is not
                            const VectorDim X(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P())))); // point at which we should contract
                            if(nA->isMovableTo(X))
                            {// nA can be moved to X
                                VerboseNodeContraction(1,"DislocationNodeContraction case 2a"<<std::endl;);
                                nA->set_P(X);
                                return DN.contractSecond(nA,nB);
                            }
                            else
                            {// nA cannot be moved to X
                                VerboseNodeContraction(1,"nA cannot be moved to X"<<std::endl;);
                                return false;
                            }
                        }
                        else
                        {// b is on box, a is not
                            const VectorDim X(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P())))); // point at which we should contract
                            if(nB->isMovableTo(X))
                            {// nA can be moved to X
                                VerboseNodeContraction(1,"DislocationNodeContraction case 3a"<<std::endl;);
                                nB->set_P(X);
                                return DN.contractSecond(nB,nA);
                            }
                            else
                            {// nB cannot be moved to X
                                VerboseNodeContraction(1,"nB cannot be moved to X"<<std::endl;);
                                return false;
                            }
                        }
                    }
                    else
                    {
                        VerboseNodeContraction(1,"bounding boxes don't intersect."<<std::endl;);
                        return false;
                    }
                    
                }
                else
                {// neither a nor b are on bounding box
                    //std::cout<<"CHECK IF NODES ARE ON GB!"<<std::endl;
                    if(nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
                    {// both nodes confined by more then one plane
                        
                        SegmentSegmentDistance<dim> ssd(nA->glidePlaneIntersections()[0].first,nA->glidePlaneIntersections()[0].second,
                                                        nB->glidePlaneIntersections()[0].first,nB->glidePlaneIntersections()[0].second);
                        
                        const auto iSeg=ssd.intersectionSegment();
                        if(iSeg.size()==1)
                        {// incident intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5a"<<std::endl;);
                            nA->set_P(std::get<0>(iSeg[0]));
                            return DN.contractSecond(nA,nB);
                        }
                        else if(iSeg.size()==2)
                        {// coincident intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5b"<<std::endl;);
                            LineSegment<dim> ls(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]));
                            nA->set_P(ls.snap(0.5*(nA->get_P()+nB->get_P())));
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {// parallel or skew intersection
                            VerboseNodeContraction(1,"DislocationNodeContraction case 5c"<<std::endl;);
                            return false;
                        }
                        
//                        BoundingLineSegments<dim> temp(nA->glidePlaneIntersections(),nB->glidePlaneIntersections());
//                        if(temp.size())
//                        {
//                            std::cout<<"contractPoint case 8"<<std::endl;
//                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
//                            return DN.contractSecond(nA,nB);
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 9"<<std::endl;
//                            return false;
//                        }
                    }
                    else if(nA->glidePlaneIntersections().size() && !nB->glidePlaneIntersections().size())
                    {// nA confined by more then one plane, nB confined by only one plane
                        PlaneLineIntersection<dim> pli(nB->glidePlane(0).P.cartesian(),
                                                       nB->glidePlane(0).unitNormal,
                                                       nA->glidePlaneIntersections()[0].first, // origin of line
                                                       nA->glidePlaneIntersections()[0].second-nA->glidePlaneIntersections()[0].first // line direction
                                                       );
                        
                        if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                        {// nothing to do, _glidePlaneIntersections remains unchanged
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6a"<<std::endl;);
                            nA->set_P(nA->snapToGlidePlaneIntersection(0.5*(nA->get_P()+nB->get_P())));
                            return DN.contractSecond(nA,nB);
                        }
                        else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                        {// _glidePlaneIntersections becomes a singular point
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6b"<<std::endl;);
                            nA->set_P(pli.P);
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {// parallel planes, cannot contract
                            VerboseNodeContraction(1,"DislocationNodeContraction case 6c"<<std::endl;);
                            return false;
                        }
                        
                        //                        BoundingLineSegments<dim> temp=nA->glidePlaneIntersections();
                        //                        temp.updateWithGlidePlane(nB->glidePlane(0));
                        //
                        //                        if(temp.size())
                        //                        {
                        //                            std::cout<<"contractPoint case 10"<<std::endl;
                        //                            nA->set_P(temp.snap(0.5*(nA->get_P()+nB->get_P()))); THE PROBLE HERE IS THAT WE NEED TO SNAP TO GLIDEPLANE INTERSECTION!!!
                        //                            return DN.contractSecond(nA,nB);
                        //                        }
                        //                        else
                        //                        {
                        //                            std::cout<<"contractPoint case 11"<<std::endl;
                        //                            return false;
                        //                        }
                    }
                    else if(!nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
                    {
                        VerboseNodeContraction(1,"DislocationNodeContraction case 7a"<<std::endl;);
                        return contract(nB,nA);
                    }
                    else
                    {// both nodes confined by only one plane
                        
                        assert(nA->glidePlanes().size()==1);
                        assert(nB->glidePlanes().size()==1);
                        
                        GlidePlaneObserver<DislocationNetworkType>* const gpo(nA->glidePlane(0).glidePlaneObserver);
                        const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(&nA->glidePlane(0),&nB->glidePlane(0)));
                        
                        //                        std::cout<<"gpA: "<<nA->glidePlane(0).P.cartesian().transpose()<<std::endl;
                        //                        std::cout<<"gpA: "<<nA->glidePlane(0).n.cartesian().normalized().transpose()<<std::endl;
                        //
                        //                        std::cout<<"gpB: "<<nB->glidePlane(0).P.cartesian().transpose()<<std::endl;
                        //                        std::cout<<"gpB: "<<nB->glidePlane(0).n.cartesian().normalized().transpose()<<std::endl;
                        
                        
                        if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 8a"<<std::endl;);
                            nA->set_P(nA->snapToGlidePlaneIntersection(0.5*(nA->get_P()+nB->get_P()))); // this snaps to the glide planes to kill numerical errors
                            return DN.contractSecond(nA,nB);
                        }
                        else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                        {
                            //                            std::cout<<nA->sID<<" "<<nB->sID<<std::endl;
                            VerboseNodeContraction(1,"DislocationNodeContraction case 8b"<<std::endl;);
                            //                            std::cout<<"norm d="<<ppi.d.norm()<<std::endl;
                            const double u=(0.5*(nA->get_P()+nB->get_P())-ppi.P).dot(ppi.d);
                            nA->set_P(ppi.P+u*ppi.d);
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 8c"<<std::endl;);
                            return false;
                        }
                    }
                }
            }
            else if(nA->isGlissile() && !nB->isGlissile())
            {// a is glissile, b is sessile
                if(nA->isOnBoundingBox())
                {// a is a boundary node. Check is the bonuding box contains b
                    if(nA->boundingBoxSegments().contains(nB->get_P()).first)
                    {
                        VerboseNodeContraction(1,"DislocationNodeContraction case 9a"<<std::endl;);
                        return DN.contractSecond(nB,nA);
                    }
                    else
                    {
                        VerboseNodeContraction(1,"DislocationNodeContraction case 9b"<<std::endl;);
                        return false;
                    }
                }
                else
                {// a is an inner node. Check if _glidePlaneIntersections contains b
                    if(nA->glidePlaneIntersections().size())
                    {// nA confined by multiple planes
                        if(nA->glidePlaneIntersections().contains(nB->get_P()).first)
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 9c"<<std::endl;);
                            return DN.contractSecond(nB,nA);
                        }
                        else
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 9d"<<std::endl;);
                            return false;
                        }
                    }
                    else
                    {// nA confined by only one plane
                        assert(nA->glidePlanes().size()==1);
                        
                        bool found=false;
                        VectorDim Y(VectorDim::Zero());

                        for(const auto& pair : nB->neighbors())
                        {
                            PlaneSegmentIntersection<dim> psi(nA->glidePlane(0).P.cartesian(),
                                                              nA->glidePlane(0).unitNormal,
                                                              std::get<1>(pair.second)->source->get_P(),
                                                              std::get<1>(pair.second)->sink->get_P()
                                                              );
//                            std::cout<<"psi.type="<<psi.type<<std::endl;
                            Y=0.5*(psi.x0+psi.x1);
//                            std::cout<<"movable="<<nB->isMovableTo(Y)<<std::endl;
//                            std::cout<<"contains="<<nA->glidePlane(0).contains(Y)<<std::endl;

                            if(nB->isMovableTo(Y) && nA->glidePlane(0).contains(Y))
                            {
                                found=true;
                                break;
                            }
                        }
                        
                        if(found)
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 9e"<<std::endl;);
                            nA->set_P(Y);
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {
                            VerboseNodeContraction(1,"DislocationNodeContraction case 9f"<<std::endl;);
                            return false;
                        }
                        
//                        if(nA->glidePlane(0).contains(nB->get_P()))
//                        {
//                            VerboseNodeContraction(1,"DislocationNodeContraction case 9e"<<std::endl;);
//                            return DN.contractSecond(nB,nA);
//                        }
//                        else
//                        {
//                            VerboseNodeContraction(1,"DislocationNodeContraction case 9f"<<std::endl;);
//                            return false;
//                        }
                    }
                }
            }
            else if(!nA->isGlissile() && nB->isGlissile())
            {
                VerboseNodeContraction(1,"DislocationNodeContraction case 10a"<<std::endl;);
                return contract(nB,nA);
            }
            else
            {// both nodes are sessile
                if((nA->get_P()-nB->get_P()).squaredNorm()<FLT_EPSILON)
                {// contract only if coincident
                    VerboseNodeContraction(1,"DislocationNodeContraction case 11a"<<std::endl;);
                    nA->set_P(0.5*(nA->get_P()+nB->get_P()));
                    return DN.contractSecond(nA,nB);
                }
                else
                {
                    VerboseNodeContraction(1,"DislocationNodeContraction case 11b"<<std::endl;);
                    return false;
                }
            }
            
        }
    };
    
    // Static data
    template <typename DislocationNetworkType>
    int DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction=0;

    
}
#endif





//                if(nA->isOnBoundingBox() && nB->isOnBoundingBox())
//                {// both nodes on bounding boxes. Intersect bounding boxes'
//                    std::cout<<"contractPoint case 1aa"<<std::endl;
//                    const bool nAisMovable=nA->isMovableTo(nB->get_P());
//                    const bool nBisMovable=nB->isMovableTo(nA->get_P());
//
//                    if(nAisMovable && nBisMovable)
//                    {
//                        VerboseNodeContraction(1,"DislocationNodeContraction case 1a"<<std::endl;);
//                        return contractYoungest(nA,nB);
//                    }
//                    else if(nAisMovable && !nBisMovable)
//                    {
//                        VerboseNodeContraction(1,"DislocationNodeContraction case 1b"<<std::endl;);
//                        return DN.contractSecond(nB,nA);
//                    }
//                    else if(!nAisMovable && nBisMovable)
//                    {
//                        VerboseNodeContraction(1,"DislocationNodeContraction case 1c"<<std::endl;);
//                        return DN.contractSecond(nA,nB);
//                    }
//                    else
//                    {
//                        VerboseNodeContraction(1,"DislocationNodeContraction case 1d"<<std::endl;);
//                        return false;
//                    }
//                }
//                else if(nA->isOnBoundingBox() && !nB->isOnBoundingBox())
//                {// a on box, b is not
//                    BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
//                    VerboseNodeContraction(1,"temp.size="<<temp.size()<<std::endl;);
//                    if(temp.size())
//                    {
//                        const VectorDim X(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P())))); // point at which we should contract
//                        if(nA->isMovableTo(X))
//                        {// nA can be moved to X
//                            VerboseNodeContraction(1,"DislocationNodeContraction case 2a"<<std::endl;);
//                            nA->set_P(X);
//                            return DN.contractSecond(nA,nB);
//                        }
//                        else
//                        {// nA cannot be moved to X
//                            VerboseNodeContraction(1,"DislocationNodeContraction case 2b"<<std::endl;);
//                            return false;
//
//                        }
//                    }
//                    else
//                    {
//                        VerboseNodeContraction(1,"DislocationNodeContraction case 2c"<<std::endl;);
//                        return false;
//                    }
//
//
//
////                    if(nB->glidePlaneIntersections().size())
////                    {// other is confined internally by two or more planes
////                        BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
////                        if(temp.size())
////                        {// An intersection of the bounding boxes exists
////                            VerboseNodeContraction(1,"DislocationNodeContraction case 2a"<<std::endl;);
////                            const VectorDim X(std::get<0>(temp.snap(nA->get_P()))); // point at which we should contract
////
////                            if(nA->isMovableTo(X))
////                            {// nA can be moved to X
////                                VerboseNodeContraction(1,"DislocationNodeContraction case 2b"<<std::endl;);
////                                nA->set_P(X);
////                                return DN.contractSecond(nA,nB);
////                            }
////                            else
////                            {// nA cannot be moved to X
////                                VerboseNodeContraction(1,"DislocationNodeContraction case 2c"<<std::endl;);
////                                return false;
////
////                            }
////                        }
////                        else
////                        {
////                            VerboseNodeContraction(1,"DislocationNodeContraction case 2d"<<std::endl;);
////                            return false;
////                        }
////                    }
////                    else
////                    {// other is confined by only one glide plane
////                        BoundingLineSegments<dim> temp=nA->boundingBoxSegments();
////                        temp.updateWithGlidePlane(nB->glidePlane(0));
////                        VerboseNodeContraction(1,"temp.size="<<temp.size()<<std::endl;);
////                        if(temp.size())
////                        {
////                            const VectorDim X(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P())))); // point at which we should contract
////                            if(nA->isMovableTo(X))
////                            {// nA can be moved to X
////                                VerboseNodeContraction(1,"DislocationNodeContraction case 3a"<<std::endl;);
////                                nA->set_P(X);
////                                return DN.contractSecond(nA,nB);
////                            }
////                            else
////                            {// nA cannot be moved to X
////                                VerboseNodeContraction(1,"DislocationNodeContraction case 3a1"<<std::endl;);
////                                return false;
////
////                            }
////                        }
////                        else
////                        {
////                            VerboseNodeContraction(1,"DislocationNodeContraction case 3b"<<std::endl;);
////                            return false;
////                        }
////                    }
//                }
//                else if(!nA->isOnBoundingBox() && nB->isOnBoundingBox())
//                {// b on box, a is not
//                    VerboseNodeContraction(1,"DislocationNodeContraction case 4a"<<std::endl;);
//                    return contract(nB,nA);
//                }
