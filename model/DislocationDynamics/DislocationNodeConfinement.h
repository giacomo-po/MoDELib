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
    /*         */ private std::set<const GlidePlane<typename TypeTraits<NodeType>::LoopNetworkType>*>,
    /*         */ private BoundingLineSegments<TypeTraits<NodeType>::dim>,
    /*         */ private std::set<const Grain<typename TypeTraits<NodeType>::LoopNetworkType>*>
    {
        
        static constexpr int dim=TypeTraits<NodeType>::dim;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef typename TypeTraits<NodeType>::LoopNetworkType NetworkType;
        typedef GlidePlane<NetworkType> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef std::set<const Grain<NetworkType>*> GrainContainerType;
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
                    GlidePlaneObserver<NetworkType>* const gpo(glidePlane(0).glidePlaneObserver);
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
            assert(_glidePlaneIntersections.size()<=1 && "_glidePlaneIntersections can have at the most size 1");
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
        std::deque<size_t> bndIDs() const
        {/*!\returns The set of IDs into BoundingLineSegments on which node is allowed to be
          */
            std::deque<size_t> temp;
            if(node->isBoundaryNode())
            {
            
            }
            else
            {
                for(size_t k=0;k<boundingBoxSegments().size();++k)
                {
                    temp.push_back(k);
                }
            }
            return temp;
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
                    if(normD2>FLT_EPSILON)
                    {
                        const double u=(P-_glidePlaneIntersections[0].first).dot(D)/normD2;
                        if(u<0.0)
                        {
                            return _glidePlaneIntersections[0].first;
                        }
                        else if(u>1.0)
                        {
                            return _glidePlaneIntersections[0].second;
                        }
                        else
                        {
                            return _glidePlaneIntersections[0].first+u*D;
                        }
                    }
                    else
                    {
                        return _glidePlaneIntersections[0].first;
                    }

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
        
    };

}
#endif
