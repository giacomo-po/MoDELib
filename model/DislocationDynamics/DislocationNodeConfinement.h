/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeConfinement_H_
#define model_DislocationNodeConfinement_H_

#include <memory>
#include <set>

#include <model/Utilities/TypeTraits.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>


namespace model
{
    template <typename NodeType>
    class DislocationNodeConfinement :
    /*         */ private std::set<const GlidePlane<typename TypeTraits<NodeType>::LoopType>*>,
    /*         */ private BoundingLineSegments<TypeTraits<NodeType>::dim>,
    /*         */ private std::set<const Grain<TypeTraits<NodeType>::dim>*>
    {
        
        static constexpr int dim=TypeTraits<NodeType>::dim;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef GlidePlane<LoopType> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef std::set<const Grain<dim>*> GrainContainerType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        /**********************************************************************/
        void updateGlidePlaneIntersections(const GlidePlaneType& lastGlidePlane)
        {
            BoundingLineSegments<dim> temp;
            
            switch (glidePlanes().size())
            {
                case 0:
                {// there must be at least one glide plane
                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                    break;
                }
                    
                case 1:
                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
                    _glidePlaneIntersections.clear();
                    break;
                }
                    
                case 2:
                {// a second plane is being added, so we must have no _glidePlaneIntersections
                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
                    
                    // Grab the infinite line of intersection between the two planes
                    GlidePlaneObserver<LoopType>* const gpo(glidePlane(0).glidePlaneObserver);
                    const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(&glidePlane(0),&glidePlane(1)));
                    
                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                    {/* Two distinct glide planes can be coincident only if they belong to different grains
                      * In that case, the intersection of their bounding boxes should be one line segment
                      */
                        assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
                        _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
                    }
                    else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                    {/* If the two planes are incident then the intersection of
                      * their bounding boxes is either a pair of singluar segments (2 points)
                      * or a line segment on the boundary
                      */
                        switch (boundingBoxSegments().size())
                        {
                            case 1:
                            {// the bounding boxes of the two planes intersect on a boundary segment. Add end points to _glidePlaneIntersections
                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
                                break;
                            }
                                
                            case 2:
                            {// The two intersections must be degenerate (2 boundary points)
                                assert((boundingBoxSegments()[0].first-boundingBoxSegments()[0].second).squaredNorm()<FLT_EPSILON);
                                assert((boundingBoxSegments()[1].first-boundingBoxSegments()[1].second).squaredNorm()<FLT_EPSILON);
                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[1].first);
                                break;
                            }
                                
                            default:
                            {
                                model::cout<<"DislocationNode "<<node->sID<<" boundingBoxSegments() are:"<<std::endl;
                                for(const auto& pair : boundingBoxSegments())
                                {
                                    model::cout<<"("<<pair.first.transpose()<<","<<pair.second.transpose()<<")"<<std::endl;
                                }
                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
                            }
                        }
                    }
                    else
                    {
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                    // Now we must have exactly one _glidePlaneIntersections
                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must have size 1");
                    
                    break;
                }
                    
                default:
                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must exist");
                    
                    // intersect the _glidePlaneIntersections with the new plane
                    PlaneLineIntersection<dim> pli(lastGlidePlane.P.cartesian(),
                                                   lastGlidePlane.n.cartesian(),
                                                   _glidePlaneIntersections[0].first, // origin of line
                                                   _glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first // line direction
                                                   );
                    
                    if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                    {// nothing to do, _glidePlaneIntersections remains unchanged
                        
                    }
                    else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                    {// _glidePlaneIntersections becomes a singular point
                        _glidePlaneIntersections[0].first =pli.P;
                        _glidePlaneIntersections[0].second=pli.P;
                    }
                    else
                    {
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                }
                    
            }
            
        }
        
        /**********************************************************************/
        BoundingLineSegments<dim> _glidePlaneIntersections; //
    
        const NodeType* const node;
        
    public:
        
        /**********************************************************************/
        DislocationNodeConfinement(const NodeType* const pN) :
        /* init */node(pN)
        {
            
        }
        
        /**********************************************************************/
        void clear()
        {
            glidePlanes().clear();
            boundingBoxSegments().clear();
            _glidePlaneIntersections.clear();
            grains().clear();
        }
        
        /**********************************************************************/
        bool addGlidePlane(const GlidePlaneType& gp)
        {
            const bool success=glidePlanes().insert(&gp).second;
            if(success)
            {
                assert(gp.contains(node->get_P()) && "Glide Plane does not contain DislocationNode");
                //                _isGlissile*=pL->loop()->isGlissile;
                boundingBoxSegments().updateWithGlidePlane(gp); // Update _boundingBoxSegments. This must be called before updateGlidePlaneIntersections
                updateGlidePlaneIntersections(gp);
                grains().insert(&(gp.grain)); // Insert new grain in grainSet
                if(grains().size()>1)
                {
                    std::cout<<"WARNING: CHECK THAT NODE IS ON REGION BND"<<std::endl;
                }
                //
                //                const VectorDim bbP(_boundingBoxSegments.snap(this->get_P()));
                //                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                //                {
                //                    set_P(bbP);
                //                    _isOnBoundingBox=true;
                //                }
            }
            return success;
        }
        
        /**********************************************************************/
        const GlidePlaneType& glidePlane(const size_t& n) const
        {
            assert(n<glidePlanes().size());
            auto iter=glidePlanes().begin();
            std::advance(iter,n);
            return **iter;
        }
        
        /**********************************************************************/
        const GlidePlaneContainerType& glidePlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GlidePlaneContainerType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainContainerType& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GrainContainerType& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& boundingBoxSegments() const
        {
            return *this;
        }
        
        /**********************************************************************/
        BoundingLineSegments<dim>& boundingBoxSegments()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& glidePlaneIntersections() const
        {
            
            return _glidePlaneIntersections;
        }
        
        /**********************************************************************/
        VectorDim snapToGlidePlaneIntersection(const VectorDim& P)
        {
            
            switch (_glidePlaneIntersections.size())
            {
                case 0:
                {
                    //                    assert(nodeConfinement().glidePlanes().size()>0);
                    assert(glidePlanes().size()==1);
                    return glidePlane(0).snapToPlane(P);
                    break;
                }
                    
                case 1:
                {
                    const VectorDim D=_glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first;
                    const double normD2(D.squaredNorm());
                    return normD2>FLT_EPSILON? _glidePlaneIntersections[0].first+(P-_glidePlaneIntersections[0].first).dot(D)*D/normD2 : _glidePlaneIntersections[0].first;
                    break;
                }
                    
                default:
                {
                    assert(0 && "THERE CAN BE AT MOST ONE LINE OF INTERSECTION");
                    return VectorDim::Zero();
                    break;
                }
            }
            
        }
        

        
//        /**********************************************************************/
//        std::tuple<bool,VectorDim,size_t,size_t> contractPoint(const NodeType* const other) const
//        {
//            
//            MOVE THIS FUNCTION TO DISLOCATION NETWORK
//            
//            //std::pair<bool,VectorDim> c=std::make_pair(false,VectorDim::Zero());
//            
//            std::cout<<node->sID<<" "<<other->sID<<std::endl;
//            
//            if(node->isGlissile() && other->isGlissile())
//            {// both nodes are glissile
//                if(node->isOnBoundingBox() && other->isOnBoundingBox())
//                {// both nodes on bounding boxes. Intersect bounding boxes'
//                    std::cout<<"contractPoint case 1aa"<<std::endl;
//                    BoundingLineSegments<dim> temp(node->boundingBoxSegments(),other->boundingBoxSegments());
//                    if(temp.size())
//                    {// a common portion of the boundary exists
//                        std::cout<<"contractPoint case 1"<<std::endl;
//                        return std::make_tuple(true,temp.snap(0.5*(node->get_P()+other->get_P())),node->sID,other->sID);
//                    }
//                    else
//                    {
//                                                std::cout<<"contractPoint case 2"<<std::endl;
//                        return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                    }
//                }
//                else if(node->isOnBoundingBox() && !other->isOnBoundingBox())
//                {// a on box, b is not
//                    if(other->nodeConfinement().glidePlaneIntersections().size())
//                    {// other is confined internally by two or more planes
//                        std::cout<<"contractPoint case 3aa"<<std::endl;
//                        BoundingLineSegments<dim> temp(node->boundingBoxSegments(),other->nodeConfinement().glidePlaneIntersections());
//                        if(temp.size())
//                        {
//                                                    std::cout<<"contractPoint case 3"<<std::endl;
//                            return std::make_tuple(true,temp.snap(0.5*(node->get_P()+other->get_P())),node->sID,other->sID);
//                        }
//                        else
//                        {
//                                                    std::cout<<"contractPoint case 4"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                    else
//                    {// other is confined by only one glide plane
//                        std::cout<<"contractPoint case 5aa"<<std::endl;
//                        BoundingLineSegments<dim> temp=node->boundingBoxSegments();
//                        temp.updateWithGlidePlane(other->glidePlane(0));
//                        if(temp.size())
//                        {
//                                                    std::cout<<"contractPoint case 5"<<std::endl;
//                            return std::make_tuple(true,temp.snap(0.5*(node->get_P()+other->get_P())),node->sID,other->sID);
//                        }
//                        else
//                        {
//                                                    std::cout<<"contractPoint case 6"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                }
//                else if(!node->isOnBoundingBox() && other->isOnBoundingBox())
//                {// b on box, a is not
//                                            std::cout<<"contractPoint case 7"<<std::endl;
//                    return other->contractPoint(node);
//                }
//                else
//                {// neither a nor b on bounding box
//                    if(glidePlaneIntersections().size() && other->nodeConfinement().glidePlaneIntersections().size())
//                    {
//                        std::cout<<"contractPoint case 8aa"<<std::endl;
//                        BoundingLineSegments<dim> temp(node->nodeConfinement().glidePlaneIntersections(),other->nodeConfinement().glidePlaneIntersections());
//                        if(temp.size())
//                        {
//                                                    std::cout<<"contractPoint case 8"<<std::endl;
//                            return std::make_tuple(true,temp.snap(0.5*(node->get_P()+other->get_P())),node->sID,other->sID);
//                        }
//                        else
//                        {
//                                                    std::cout<<"contractPoint case 9"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                    else if(glidePlaneIntersections().size() && !other->nodeConfinement().glidePlaneIntersections().size())
//                    {
//                        std::cout<<"contractPoint case 10aa"<<std::endl;
//                        BoundingLineSegments<dim> temp=node->boundingBoxSegments();
//                        temp.updateWithGlidePlane(other->glidePlane(0));
//                        if(temp.size())
//                        {
//                                                    std::cout<<"contractPoint case 10"<<std::endl;
//                            return std::make_tuple(true,temp.snap(0.5*(node->get_P()+other->get_P())),node->sID,other->sID);
//                        }
//                        else
//                        {
//                                                    std::cout<<"contractPoint case 11"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                    else if(!glidePlaneIntersections().size() && other->nodeConfinement().glidePlaneIntersections().size())
//                    {
//                                                std::cout<<"contractPoint case 12"<<std::endl;
//                        return other->contractPoint(node);
//                    }
//                    else
//                    {// both nodes confined by only one plane
////                        std::cout<<"contractPoint case 13aa"<<std::endl;
//
//                        assert(glidePlanes().size()==1);
//                        assert(other->nodeConfinement().glidePlanes().size()==1);
//                        
//                        GlidePlaneObserver<LoopType>* const gpo(glidePlane(0).glidePlaneObserver);
//                        const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(&glidePlane(0),&other->nodeConfinement().glidePlane(0)));
//                        
//                        if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
//                        {
//                                                    std::cout<<"contractPoint case 13"<<std::endl;
//                            return std::make_tuple(true,0.5*(node->get_P()+other->get_P()),node->sID,other->sID);
//                        }
//                        else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
//                        {
//                            std::cout<<node->sID<<" "<<other->sID<<std::endl;
//                            std::cout<<"contractPoint case 13a"<<std::endl;
//                            
//                            const double u=(0.5*(node->get_P()+other->get_P())+ppi.P).dot(ppi.d);
//                            const VectorDim Pmin=ppi.P+u*ppi.d;
//                            node->nodeConfinement().updateGlidePlaneIntersections(other->nodeConfinement().glidePlane(0));
////                            std::cout<<"u="<<u<<std::endl;
////                            std::cout<<"NEED TO SNAP P TO THE GLIDEPLANE INTERSECTION"
//                            return std::make_tuple(true,snapToGlidePlaneIntersection(Pmin),node->sID,other->sID);
////
////                            
////                            assert(0 && "FINISH HERE, RETURN POINT ON INTTERSECTIOn LINE WHICH MINIMIZES DISTANCES");
//                        }
//                        else
//                        {
//                                                    std::cout<<"contractPoint case 14"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                }
//            }
//            else if(node->isGlissile() && !other->isGlissile())
//            {// a is glissile, b is sessile
//                if(node->isOnBoundingBox())
//                {// a is a boundary node. Check is the bonuding box contains b
//                    if(node->boundingBoxSegments().contains(other->get_P()))
//                    {
//                                                std::cout<<"contractPoint case 15"<<std::endl;
//                        return std::make_tuple(true,other->get_P(),other->sID,node->sID);
//                    }
//                    else
//                    {
//                                                std::cout<<"contractPoint case 16"<<std::endl;
//                        return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                    }
//                }
//                else
//                {// a is an inner node. Check if _glidePlaneIntersections contains b
//                    if(glidePlaneIntersections().size())
//                    {
//                        if(glidePlaneIntersections().contains(other->get_P()))
//                        {
//                            std::cout<<"contractPoint case 17"<<std::endl;
//                            return std::make_tuple(true,other->get_P(),other->sID,node->sID);
//
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 17a"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                    else
//                    {
//                        assert(glidePlanes().size()==1);
//                        if(glidePlane(0).contains(other->get_P()))
//                        {
//                            std::cout<<"contractPoint case 17b"<<std::endl;
//                            return std::make_tuple(true,other->get_P(),other->sID,node->sID);
//
//                        }
//                        else
//                        {
//                            std::cout<<"contractPoint case 17c"<<std::endl;
//                            return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                        }
//                    }
//                    //assert(0 && "FINISH HERE");
//                }
//                
//            }
//            else if(!node->isGlissile() && other->isGlissile())
//            {
//                                        std::cout<<"contractPoint case 18"<<std::endl;
//                return other->contractPoint(node);
//            }
//            else
//            {// both nodes are sessile
//                if((node->get_P()-other->get_P()).squaredNorm()<FLT_EPSILON)
//                {// contract only if coincident
//                                            std::cout<<"contractPoint case 19"<<std::endl;
//                    return std::make_tuple(true,0.5*(node->get_P()+other->get_P()),node->sID,other->sID);
//                }
//                else
//                {
//                                            std::cout<<"contractPoint case 20"<<std::endl;
//                    return std::make_tuple(false,VectorDim::Zero(),node->sID,other->sID);
//                }
//            }
//            
//        }
        
        
    };
}

#endif
