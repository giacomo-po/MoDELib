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


namespace model
{
    template <typename DislocationNetworkType>
    struct DislocationNodeContraction
    {
        
        static constexpr int dim=DislocationNetworkType::dim;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
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
            
            bool success=false;
            
            
            std::cout<<"DislocationNetwork::contract "<<nA->sID<<" "<<nB->sID<<std::endl;
            
            
            
            if(nA->isGlissile() && nB->isGlissile())
            {// both nodes are glissile
                if(nA->isOnBoundingBox() && nB->isOnBoundingBox())
                {// both nodes on bounding boxes. Intersect bounding boxes'
                    
                    std::cout<<"contractPoint case 1aa"<<std::endl;

//                    const bool nAisRemovable=isMovable(nA,nB,nB->get_P());
//                    const bool nBisRemovable=isMovable(nB,nA,nA->get_P());

                    const bool nAisRemovable=nA->isMovableTo(nB->get_P());
                    const bool nBisRemovable=nB->isMovableTo(nA->get_P());

                    
                    
//                    std::cout<<"nA neighbors:"<<std::endl;
//                    bool nAisRemovable=true;
//                    for(const auto& pair : nA->neighbors())
//                    {
//                        //                            std::cout<<"neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first<<std::endl;
//                        if(pair.first!=nB->sID && std::get<1>(pair.second)->isBoundarySegment())
//                        {// boundary segments other than A->B must remain boundary
//                            nAisRemovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first;
//                        }
//                        
//                        if(pair.first!=nB->sID && std::get<1>(pair.second)->isGrainBoundarySegment())
//                        {// boundary segments other than A->B must remain boundary
//                            for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
//                            {
//                                std::cout<<"grainBoundary neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nB->get_P())<<std::endl;
//                                nAisRemovable*=gb->glidePlanes().begin()->second.contains(nB->get_P());
//                            }
//                        }
//                    }
//                    
//                    std::cout<<"nB neighbors:"<<std::endl;
//                    bool nBisRemovable=true;
//                    for(const auto& pair : nB->neighbors())
//                    {
//                        //                            std::cout<<"neighbor "<<nB->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first<<std::endl;
//                        if(pair.first!=nA->sID && std::get<1>(pair.second)->isBoundarySegment())
//                        {// boundary segments other than A->B must remain boundary
//                            //                                std::cout<<"boundary neighbor "<<nB->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first<<std::endl;
//                            //                                nBisRemovable*=temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first;
//                            nBisRemovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+nA->get_P())).first;
//                        }
//                        
//                        if(pair.first!=nA->sID && std::get<1>(pair.second)->isGrainBoundarySegment())
//                        {// boundary segments other than A->B must remain boundary
//                            for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
//                            {
//                                std::cout<<"grainBoundary neighbor "<<nB->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nA->get_P())<<std::endl;
//                                nBisRemovable*=gb->glidePlanes().begin()->second.contains(nA->get_P());
//                            }
//                        }
//                        
//                    }
                    
                    if(nAisRemovable && nBisRemovable)
                    {
                        std::cout<<"contractPoint case 1a"<<std::endl;
                        //                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
                        //                            return DN.contractSecond(nA,nB);
                        return contractYoungest(nA,nB);
                    }
                    else if(nAisRemovable && !nBisRemovable)
                    {
                        std::cout<<"contractPoint case 1b"<<std::endl;
                        return DN.contractSecond(nB,nA);
                    }
                    else if(!nAisRemovable && nBisRemovable)
                    {
                        std::cout<<"contractPoint case 1c"<<std::endl;
                        return DN.contractSecond(nA,nB);
                    }
                    else
                    {
                        std::cout<<"contractPoint case 1d"<<std::endl;
                        return false;
                    }
                    
                    
                    
                }
                else if(nA->isOnBoundingBox() && !nB->isOnBoundingBox())
                {// a on box, b is not
                    if(nB->glidePlaneIntersections().size())
                    {// other is confined internally by two or more planes
                        std::cout<<"contractPoint case 3aa"<<std::endl;
                        //                        BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->glidePlaneIntersections()); // NO THIS IS WRONG. USE nB->boundingBoxSegments()
                        BoundingLineSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
                        if(temp.size())
                        {// An intersection of the bonuding boxes exists
                            std::cout<<"contractPoint case 3"<<std::endl;
                            //                            const VectorDim X(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P())))); // point at which we will contract
                            const VectorDim X(std::get<0>(temp.snap(nA->get_P()))); // point at which we should contract
                            
//                            if(isMovable(nA,nB,X))
                                if(nA->isMovableTo(X))
                            {// nA can be moved to X
                                std::cout<<"contractPoint case 3a1"<<std::endl;
                                nA->set_P(X);
                                return DN.contractSecond(nA,nB);
                            }
                            else
                            {// nA cannot be moved to X
                                std::cout<<"contractPoint case 3a2"<<std::endl;
                                return false;

                            }
                            
//                            const size_t A_ID(std::get<1>(nA->boundingBoxSegments().snap(nA->get_P())));    // bonuding box line ID containing A
//                            if(nA->boundingBoxSegments().contains(X,A_ID))
//                            {
//                                std::cout<<"contractPoint case 3a1"<<std::endl;
//                                nA->set_P(X);
//                                return DN.contractSecond(nA,nB);
//                            }
//                            else
//                            {
//                                                                std::cout<<"contractPoint case 3a2"<<std::endl;
//                                return false;
//                            }
                        }
                        else
                        {
                            std::cout<<"contractPoint case 4"<<std::endl;
                            return false;
                        }
                    }
                    else
                    {// other is confined by only one glide plane
                        std::cout<<"contractPoint case 5aa"<<std::endl;
                        BoundingLineSegments<dim> temp=nA->boundingBoxSegments();
                        temp.updateWithGlidePlane(nB->glidePlane(0));
                        std::cout<<"temp.size="<<temp.size()<<std::endl;
                        if(temp.size())
                        {
                            std::cout<<"contractPoint case 5"<<std::endl;
                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {
                            std::cout<<"contractPoint case 6"<<std::endl;
                            return false;
                        }
                    }
                }
                else if(!nA->isOnBoundingBox() && nB->isOnBoundingBox())
                {// b on box, a is not
                    std::cout<<"contractPoint case 7"<<std::endl;
                    return contract(nB,nA);
                }
                else
                {// neither a nor b on bounding box
                    
                    std::cout<<"CHECK IF NODES ARE ON GB!"<<std::endl;
                    
                    if(nA->glidePlaneIntersections().size() && nB->glidePlaneIntersections().size())
                    {// both nodes confined by more then one plane
                        std::cout<<"contractPoint case 8aa"<<std::endl;
                        BoundingLineSegments<dim> temp(nA->glidePlaneIntersections(),nB->glidePlaneIntersections());
                        if(temp.size())
                        {
                            std::cout<<"contractPoint case 8"<<std::endl;
                            nA->set_P(std::get<0>(temp.snap(0.5*(nA->get_P()+nB->get_P()))));
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {
                            std::cout<<"contractPoint case 9"<<std::endl;
                            return false;
                        }
                    }
                    else if(nA->glidePlaneIntersections().size() && !nB->glidePlaneIntersections().size())
                    {// nA confined by more then one plane, nB confined by only one plane
                        std::cout<<"contractPoint case 10aa"<<std::endl;
                        
                        PlaneLineIntersection<dim> pli(nB->glidePlane(0).P.cartesian(),
                                                       nB->glidePlane(0).unitNormal,
                                                       nA->glidePlaneIntersections()[0].first, // origin of line
                                                       nA->glidePlaneIntersections()[0].second-nA->glidePlaneIntersections()[0].first // line direction
                                                       );
                        
                        if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                        {// nothing to do, _glidePlaneIntersections remains unchanged
                            std::cout<<"contractPoint case 10aa1"<<std::endl;
                            nA->set_P(nA->snapToGlidePlaneIntersection(0.5*(nA->get_P()+nB->get_P())));
                            return DN.contractSecond(nA,nB);
                        }
                        else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                        {// _glidePlaneIntersections becomes a singular point
                            std::cout<<"contractPoint case 10aa2"<<std::endl;
                            nA->set_P(pli.P);
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {// parallel planes, cannot contract
                            std::cout<<"contractPoint case 10aa3"<<std::endl;
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
                        std::cout<<"contractPoint case 12"<<std::endl;
                        return contract(nB,nA);
                    }
                    else
                    {// both nodes confined by only one plane
                        std::cout<<"contractPoint case 13aa"<<std::endl;
                        
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
                            std::cout<<"contractPoint case 13"<<std::endl;
                            nA->set_P(nA->snapToGlidePlaneIntersection(0.5*(nA->get_P()+nB->get_P()))); // this snaps to the glide planes to kill numerical errors
                            return DN.contractSecond(nA,nB);
                        }
                        else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                        {
                            //                            std::cout<<nA->sID<<" "<<nB->sID<<std::endl;
                            std::cout<<"contractPoint case 13a"<<std::endl;
                            //                            std::cout<<"norm d="<<ppi.d.norm()<<std::endl;
                            const double u=(0.5*(nA->get_P()+nB->get_P())-ppi.P).dot(ppi.d);
                            nA->set_P(ppi.P+u*ppi.d);
                            return DN.contractSecond(nA,nB);
                        }
                        else
                        {
                            std::cout<<"contractPoint case 14"<<std::endl;
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
                        std::cout<<"contractPoint case 15"<<std::endl;
                        return DN.contractSecond(nB,nA);
                    }
                    else
                    {
                        std::cout<<"contractPoint case 16"<<std::endl;
                        return false;
                    }
                }
                else
                {// a is an inner node. Check if _glidePlaneIntersections contains b
                    if(nA->glidePlaneIntersections().size())
                    {
                        if(nA->glidePlaneIntersections().contains(nB->get_P()).first)
                        {
                            std::cout<<"contractPoint case 17"<<std::endl;
                            return DN.contractSecond(nB,nA);
                        }
                        else
                        {
                            std::cout<<"contractPoint case 17a"<<std::endl;
                            return false;
                        }
                    }
                    else
                    {
                        assert(nA->glidePlanes().size()==1);
                        if(nA->glidePlane(0).contains(nB->get_P()))
                        {
                            std::cout<<"contractPoint case 17b"<<std::endl;
                            return DN.contractSecond(nB,nA);
                        }
                        else
                        {
                            std::cout<<"contractPoint case 17c"<<std::endl;
                            return false;
                        }
                    }
                }
            }
            else if(!nA->isGlissile() && nB->isGlissile())
            {
                std::cout<<"contractPoint case 18"<<std::endl;
                return contract(nB,nA);
            }
            else
            {// both nodes are sessile
                if((nA->get_P()-nB->get_P()).squaredNorm()<FLT_EPSILON)
                {// contract only if coincident
                    std::cout<<"contractPoint case 19"<<std::endl;
                    nA->set_P(0.5*(nA->get_P()+nB->get_P()));
                    return DN.contractSecond(nA,nB);
                }
                else
                {
                    std::cout<<"contractPoint case 20"<<std::endl;
                    return false;
                }
            }
            
        }
    };
    
}
#endif
