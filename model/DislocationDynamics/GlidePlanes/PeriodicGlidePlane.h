/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_H_
#define model_PeriodicGlidePlane_H_


#include <memory>
#include <string>


#include <GlidePlane.h>
#include <GlidePlaneFactory.h>
#include <TerminalColors.h>




namespace model
{
    
    template<int dim>
    class PeriodicGlidePlane;
    
    template<int dim>
    class PeriodicPlanePatch;
    
    template<int dim>
    class PeriodicPlaneNode;
    
    template<int dim>
    class PeriodicPlaneEdge;
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlaneEdge
    {
        
        const PeriodicPlanePatch<dim>* const patch;
        const std::shared_ptr<PeriodicPlaneNode<dim>> source;
        const std::shared_ptr<PeriodicPlaneNode<dim>>   sink;
        const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection;
        PeriodicPlaneEdge<dim>* next;
        PeriodicPlaneEdge<dim>* prev;
        PeriodicPlaneEdge<dim>* twin;
        
        /**********************************************************************/
        PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                          const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in)
        ;
        ~PeriodicPlaneEdge();
        
        std::string tag() const;
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct NodalConnectivity
    {
        typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
        
        PeriodicPlaneEdgeType* inEdge;
        PeriodicPlaneEdgeType* outEdge;
        
        NodalConnectivity() :
        /* init */ inEdge(nullptr)
        /* init */,outEdge(nullptr)
        {
            
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class PeriodicPlaneNode : public StaticID<PeriodicPlaneNode<dim>>
    /*                   */,public Eigen::Matrix<double,dim-1,1>
    {
        
    public:
        
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef PeriodicPlaneNode<dim> PeriodicPlaneNodeType;
        typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
        typedef NodalConnectivity<dim> NodalConnectivityType;
        typedef std::map<size_t,NodalConnectivity<dim>> NodalConnectivityContainerType;
        typedef std::set<PeriodicPlaneEdgeType*> InOutEdgeContainerType;
        
    private:
        
        NodalConnectivityContainerType _patchConnectivities;
        NodalConnectivityContainerType _neighborConnectivities;
        
        
    public:
        
        /**********************************************************************/
        PeriodicPlaneNode(const VectorLowerDim& pos) :
        /* init*/ VectorLowerDim(pos)
        {

        }
        

        
        /**********************************************************************/
        void addLink(PeriodicPlaneEdgeType* const link)
        {
            
            if(link->sink->sID==this->sID)
            {// edge ends at this node, so link is an inLink
                
                // Update _patchConnectivities
                NodalConnectivityType& patchConnectivity(patchConnectivities()[link->patch->sID]);
                assert(patchConnectivity.inEdge==nullptr);
                patchConnectivity.inEdge=link;
                if(patchConnectivity.outEdge)
                {// and outEdge of the same patch is connected to this node
                    assert(patchConnectivity.outEdge->prev==nullptr);
                    assert(link->next==nullptr);
                    patchConnectivity.outEdge->prev=link;
                    link->next=patchConnectivity.outEdge;
                }
                
                // Update _neighborConnectivities
                NodalConnectivityType& neighborConnectivity(neighborConnectivities()[link->source->sID]);
                assert(neighborConnectivity.inEdge==nullptr || neighborConnectivity.inEdge==link);
                neighborConnectivity.inEdge=link;
                if(neighborConnectivity.outEdge)
                {// and outEdge of the same patch is connected to this node
                    assert(neighborConnectivity.outEdge->twin==nullptr || neighborConnectivity.outEdge->twin==link);
                    assert(link->twin==nullptr || link->twin==neighborConnectivity.outEdge);
                    neighborConnectivity.outEdge->twin=link;
                    link->twin=neighborConnectivity.outEdge;
                }
                
            }
            else if(link->source->sID==this->sID)
            {// edge starts at this node, so link is an outLink
                
                // Update _patchConnectivities
                NodalConnectivityType& patchConnectivity(patchConnectivities()[link->patch->sID]);
                assert(patchConnectivity.outEdge==nullptr);
                patchConnectivity.outEdge=link;
                if(patchConnectivity.inEdge)
                {
                    assert(patchConnectivity.inEdge->next==nullptr);
                    assert(link->prev==nullptr);
                    patchConnectivity.inEdge->next=link;
                    link->prev=patchConnectivity.inEdge;
                }
                
                // Update _neighborConnectivities
                NodalConnectivityType& neighborConnectivity(neighborConnectivities()[link->sink->sID]);
                assert(neighborConnectivity.outEdge==nullptr || neighborConnectivity.outEdge==link);
                neighborConnectivity.outEdge=link;
                if(neighborConnectivity.inEdge)
                {// and outEdge of the same patch is connected to this node
                    assert(neighborConnectivity.inEdge->twin==nullptr || neighborConnectivity.inEdge->twin==link);
                    assert(link->twin==nullptr || link->twin==neighborConnectivity.inEdge);
                    neighborConnectivity.inEdge->twin=link;
                    link->twin=neighborConnectivity.inEdge;
                }
            }
            else
            {
                assert(false && "CONNECTING LINK TO NON-INCIDENT NODE");
            }
            
        }
        
        /**********************************************************************/
        void removeLink(PeriodicPlaneEdgeType* const link)
        {
            
            if(link->sink->sID==this->sID)
            {// edge ends at this node, so link is an inLink
                
                // Update _patchConnectivities
                auto patchIter(patchConnectivities().find(link->patch->sID));
                assert(patchIter!=patchConnectivities().end() && "PATCH NOT FOUND IN patchConnectivities");
                NodalConnectivityType& patchConnectivity(patchIter->second);
                assert(patchConnectivity.inEdge==link);
                patchConnectivity.inEdge=nullptr;
                if(patchConnectivity.inEdge==nullptr && patchConnectivity.outEdge==nullptr)
                {
                    patchConnectivities().erase(patchIter);
                }
                
                // Update _neighborConnectivities
                auto neighborIter(neighborConnectivities().find(link->source->sID));
                assert(neighborIter!=neighborConnectivities().end() && "NEIGHBOR NOT FOUND IN neighborConnectivity");
                NodalConnectivityType& neighborConnectivity(neighborIter->second);
                assert(neighborConnectivity.inEdge==link);
                neighborConnectivity.inEdge=nullptr;
                if(neighborConnectivity.inEdge==nullptr && neighborConnectivity.outEdge==nullptr)
                {
                    neighborConnectivities().erase(neighborIter);
                }
                
            }
            else if(link->source->sID==this->sID)
            {// edge starts at this node, so link is an outLink
                
                // Update _patchConnectivities
                auto patchIter(patchConnectivities().find(link->patch->sID));
                assert(patchIter!=patchConnectivities().end() && "PATCH NOT FOUND IN patchConnectivities");
                NodalConnectivityType& patchConnectivity(patchIter->second);
                assert(patchConnectivity.outEdge==link);
                patchConnectivity.outEdge=nullptr;
                if(patchConnectivity.inEdge==nullptr && patchConnectivity.outEdge==nullptr)
                {
                    patchConnectivities().erase(patchIter);
                }
                
                // Update _neighborConnectivities
                auto neighborIter(neighborConnectivities().find(link->sink->sID));
                assert(neighborIter!=neighborConnectivities().end() && "NEIGHBOR NOT FOUND IN neighborConnectivity");
                NodalConnectivityType& neighborConnectivity(neighborIter->second);
                assert(neighborConnectivity.outEdge==link);
                neighborConnectivity.outEdge=nullptr;
                if(neighborConnectivity.inEdge==nullptr && neighborConnectivity.outEdge==nullptr)
                {
                    neighborConnectivities().erase(neighborIter);
                }
            }
            else
            {
                assert(false && "DISCONNECTING LINK FROM NON-INCIDENT NODE");
            }
            
            
        }
        
        /**********************************************************************/
        const NodalConnectivityContainerType& patchConnectivities() const
        {
            return _patchConnectivities;
        }
        
        /**********************************************************************/
        NodalConnectivityContainerType& patchConnectivities()
        {
            return _patchConnectivities;
        }
        
        /**********************************************************************/
        const NodalConnectivityContainerType& neighborConnectivities() const
        {
            return _neighborConnectivities;
        }
        
        /**********************************************************************/
        NodalConnectivityContainerType& neighborConnectivities()
        {
            return _neighborConnectivities;
        }
        
        
        
        /**********************************************************************/
        PeriodicPlaneEdgeType* inEdge() const
        {
            InOutEdgeContainerType temp;
            for(const auto& connectivity : neighborConnectivities())
            {
                if(connectivity.second.inEdge && !connectivity.second.outEdge)
                {
                    assert(connectivity.second.inEdge->twin==nullptr);
                    temp.insert(connectivity.second.inEdge);
                }
            }
            return temp.size()==1? *temp.begin() : nullptr;
        }
        
        /**********************************************************************/
        PeriodicPlaneEdgeType* outEdge() const
        {
            InOutEdgeContainerType temp;
            for(const auto& connectivity : neighborConnectivities())
            {
                if(connectivity.second.outEdge && !connectivity.second.inEdge)
                {
                    assert(connectivity.second.outEdge->twin==nullptr);
                    temp.insert(connectivity.second.outEdge);
                }
            }
            return temp.size()==1? *temp.begin() : nullptr;
        }
        
        
        
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlanePatch : public StaticID<PeriodicPlanePatch<dim>>
    /*                      */, private std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        PeriodicGlidePlane<dim>* const periodicPlane;
        const VectorDim shift;
        const std::shared_ptr<GlidePlane<dim>> glidePlane;
        typedef std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>> PeriodicPlaneEdgeContainerType;
        
        /**********************************************************************/
        PeriodicPlanePatch(PeriodicGlidePlane<dim>& periodicPlane_in,
                           const VectorDim& shift_in) :
        /* init */ periodicPlane(&periodicPlane_in)
        /* init */,shift(shift_in)
        /* init */,glidePlane(periodicPlane->getGlidePlane(shift))
        {
//            std::cout<<cyanColor<<"Creating patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
            
            if(isRightHandedBoundary(glidePlane->meshIntersections,*periodicPlane->referencePlane))
            {
                addMeshIntersections(glidePlane->meshIntersections);
            }
            else
            {
                BoundingMeshSegments<dim> flippedTemp;
                for (auto it = glidePlane->meshIntersections.rbegin(); it != glidePlane->meshIntersections.rend(); ++it)
                {
                    flippedTemp.emplace_back(new MeshBoundarySegment<dim>((*it)->P1,(*it)->P0,(*it)->faces));
                }
                addMeshIntersections(flippedTemp);
            }
            
            periodicPlane->updateBoundaries();
            
        }
        
        /**********************************************************************/
        ~PeriodicPlanePatch()
        {
//            std::cout<<cyanColor<<"Destroying patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
            this->clear(); // call before updateBoundaries
            periodicPlane->updateBoundaries();
        }

        /**********************************************************************/
        static bool isRightHandedBoundary(const BoundingMeshSegments<dim>& bnd,const Plane<dim>& plane)
        {
            if(bnd.size())
            {
                VectorDim nA(VectorDim::Zero());
                const VectorDim P0(bnd.front()->P0);
                for(const auto& seg : bnd)
                {
                    nA+= 0.5*(seg->P0-P0).cross(seg->P1-seg->P0);
                }
                return nA.dot(plane.unitNormal)>0.0;
            }
            else
            {
                return false;
            }
        }
        
        /**********************************************************************/
        void addMeshIntersections(const BoundingMeshSegments<dim>& bms)
        {
            for(size_t k=0;k<bms.size();++k)
            {// compute 2d points on PeriodicPlane, and coneect them
                //size_t k1(k0!=bms.size()-1? k0+1 : 0);
                std::shared_ptr<PeriodicPlaneNode<dim>> source(periodicPlane->getSharedNode(bms[k]->P0-shift));
                std::shared_ptr<PeriodicPlaneNode<dim>>   sink(periodicPlane->getSharedNode(bms[k]->P1-shift));
                this->emplace_back(new PeriodicPlaneEdge<dim>(this,source,sink,bms[k]));
            }
        }
        
        /**********************************************************************/
        const PeriodicPlaneEdgeContainerType& edges() const
        {
            return *this;
        }
        

    };
    
    
    
    template<int dim>
    struct TypeTraits<PeriodicGlidePlane<dim>>
    {
        typedef Eigen::Matrix<double,dim,1> KeyType;
        typedef PeriodicPlanePatch<dim> ValueType;
        typedef CompareVectorsByComponent<double,dim,float> CompareType;
    };
    

    template<int dim>
    struct PeriodicGlidePlaneBase : private std::map<Eigen::Matrix<double,dim-1,1>,PeriodicPlaneNode<dim>* const,CompareVectorsByComponent<double,dim-1,float>>
    /*                           */,private std::set<const PeriodicPlaneEdge<dim>*>

    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef std::map<Eigen::Matrix<double,dim-1,1>,PeriodicPlaneNode<dim>* const,CompareVectorsByComponent<double,dim-1,float>> NodeCointainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;

        BoundariesContainerType _outerBoundaries;
        BoundariesContainerType _innerBoundaries;

        
        /**********************************************************************/
        static MatrixDim getL2G(const VectorDim& x,
                                const VectorDim& z)
        {
            const double xNorm(x.norm());
            const double zNorm(z.norm());
            assert(xNorm>FLT_EPSILON);
            assert(zNorm>FLT_EPSILON);
            assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
            MatrixDim temp(Eigen::Matrix3d::Identity());
            temp.col(2)=z/zNorm;
            temp.col(0)=x/xNorm;
            temp.col(1)=temp.col(2).cross(temp.col(0));
            return temp;
        }
        
        GlidePlaneFactory<dim>& glidePlaneFactory;
        const std::shared_ptr<GlidePlane<dim>> referencePlane;
        const MatrixDim L2G;
        
        PeriodicGlidePlaneBase(GlidePlaneFactory<dim>& glidePlaneFactory_in,
                           const GlidePlaneKey<dim>& referencePlaneKey,
                           const VectorDim& globalX) :
        /* init */ glidePlaneFactory(glidePlaneFactory_in)
        /* init */,referencePlane(glidePlaneFactory.get(referencePlaneKey))
        /* init */,L2G(getL2G(globalX,referencePlane->unitNormal))
        {
        }
        
//        ~PeriodicGlidePlaneBase()
//        {
//            std::cout<<"DESTROYING PeriodicGlidePlaneBase"<<std::endl;
//        }
        
        



        
        /**********************************************************************/
        GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift)
        {
            return GlidePlaneKey<dim>(referencePlane->grain.grainID,referencePlane->P+shift,referencePlane->n);
        }
        
        std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift)
        {
            return glidePlaneFactory.get(getGlidePlaneKey(shift));
        }
        
        BoundariesContainerType& outerBoundaries()
        {
            return _outerBoundaries;
        }
        
        const BoundariesContainerType& outerBoundaries() const
        {
            return _outerBoundaries;
        }
        
        BoundariesContainerType& innerBoundaries()
        {
            return _innerBoundaries;
        }
        
        const BoundariesContainerType& innerBoundaries() const
        {
            return _innerBoundaries;
        }
        
        const UntwinnedEdgeContainerType& untwinnedEdges() const
        {
            return *this;
        }
        
        UntwinnedEdgeContainerType& untwinnedEdges()
        {
            return *this;
        }
        
        void addUntwinnedEdge(const PeriodicPlaneEdge<dim>* link)
        {
            untwinnedEdges().insert(link);
        }
        
        void removeUntwinnedEdge(const PeriodicPlaneEdge<dim>* link)
        {
            const size_t erased(untwinnedEdges().erase(link));
            assert(erased==1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
        }
        
        /**********************************************************************/
        NodeCointainerType& nodes()
        {
            return *this;
        }
        
        const NodeCointainerType& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        int isInsideOuterBoundary(const VectorLowerDim& test)
        {
            if(isCompact())
            {
                const auto& outerBoundary(outerBoundaries()[0]);
                size_t i, j, c = 0;
                for (i = 0, j = outerBoundary.size()-1; i < outerBoundary.size(); j = i++)
                {
                    const VectorLowerDim& Pi(*outerBoundary[i]->source);
                    const VectorLowerDim& Pj(*outerBoundary[j]->source);

                    if ( ((Pi(1)>test(1)) != (Pj(1)>test(1))) &&
                        (test(0) < (Pj(0)-Pi(0)) * (test(1)-Pi(1)) / (Pj(1)-Pi(1)) + Pi(0)) )
                    {
                        c = !c;
                    }
                }
                return c;
            }
            else
            {
                return 0;
            }
        }
        
        bool isCompact() const
        {
            return outerBoundaries().size()==1 && innerBoundaries().size()==0;
        }
        
        void updateBoundaries()
        {
//            std::cout<<"Updating OuterBoundary"<<std::endl;
            outerBoundaries().clear();
            innerBoundaries().clear();
            UntwinnedEdgeContainerType untwinnedCopy(untwinnedEdges());
            while(untwinnedCopy.size())
            {
                BoundaryContainerType temp;
                //                outerBoundaries().push_back(); // add new vector
                temp.reserve(untwinnedCopy.size());
                temp.push_back(*untwinnedCopy.begin());
                const size_t erased(untwinnedCopy.erase(temp.back()));
                if(erased!=1)
                {
                    std::cout<<"Trying to erase "<<temp.back()->tag()<<std::endl;
                    std::cout<<"untwinnedCopy is"<<std::endl;
                    for(const auto& edgePtr : untwinnedCopy)
                    {
                        std::cout<<"    "<<edgePtr->tag()<<std::endl;
                    }
                    assert(erased==1 && "could not find link in untwinnedEdges 1");
                }
                while(temp.back()->sink->sID!=temp.front()->source->sID)
                {
                    if(temp.back()->sink->outEdge())
                    {// there is a unique outLink from the sink of the current edge, pick that one
                        temp.push_back(temp.back()->sink->outEdge());
                        const size_t erased(untwinnedCopy.erase(temp.back()));
                        if(erased!=1)
                        {
                            std::cout<<"Trying to erase "<<temp.back()->tag()<<std::endl;
                            std::cout<<"untwinnedCopy is"<<std::endl;
                            for(const auto& edgePtr : untwinnedCopy)
                            {
                                std::cout<<"    "<<edgePtr->tag()<<std::endl;
                            }
                            assert(erased==1 && "could not find link in untwinnedEdges 2");
                        }
                    }
                    else
                    {// There could be two or more patches sharing node temp.back()->sink, follow current patch
                        if(temp.back()->next)
                        {
                            if(!temp.back()->next->twin)
                            {// the next link on current patch is untwinned, pick that one
                                temp.push_back(temp.back()->next);
                                const size_t erased(untwinnedCopy.erase(temp.back()));
                                if(erased!=1)
                                {
                                    std::cout<<"Trying to erase "<<temp.back()->tag()<<std::endl;
                                    std::cout<<"untwinnedCopy is"<<std::endl;
                                    for(const auto& edgePtr : untwinnedCopy)
                                    {
                                        std::cout<<"    "<<edgePtr->tag()<<std::endl;
                                    }
                                    assert(erased==1 && "could not find link in untwinnedEdges 2");
                                }

                            }
                            else
                            {
                                std::cout<<"Could not close OuterBoundary 1"<<std::endl;
                                assert(false && "Could not construct OuterBoundary");
                            }
                        }
                        else
                        {
                            std::cout<<"Could not close OuterBoundary 2"<<std::endl;
                            assert(false && "Could not construct OuterBoundary");
                        }
                    }
                }
                if(temp.size())
                {
                    if(isRightHandedBoundary(temp))
                    {
                        _outerBoundaries.push_back(temp);
                    }
                    else
                    {
                        _innerBoundaries.push_back(temp);

                    }
                }
            }

//            std::cout<<"OuterBoundaries #"<<outerBoundaries().size()<<std::endl;
//            for(const auto& bnd : outerBoundaries())
//            {
//                std::cout<<"size="<<bnd.size()<<", area="<<rightHandedArea(bnd)<<std::endl;
//                for (const auto bnd2 : bnd)
//                std::cout<<bnd2->tag()<<std::endl;
//            }
//
//            std::cout<<"InnerBoundaries #"<<innerBoundaries().size()<<std::endl;
//            for(const auto& bnd : innerBoundaries())
//            {
//                std::cout<<"size="<<bnd.size()<<", area="<<rightHandedArea(bnd)<<std::endl;
//                for (const auto bnd2 : bnd)
//                std::cout<<bnd2->tag()<<std::endl;
//                assert(0 && "Code in the right Handed Area");
//            }
            
        }
        
        /**********************************************************************/
        static VectorDim rightHandedNormal(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
        {// TO DO: THIS SHOULD USE THE 2D POINTS
            if(bnd.size()>=3)
            {
                VectorDim nA(VectorDim::Zero());
                const VectorDim P0(bnd.front()->meshIntersection->P0-(bnd.front())->patch->shift);
                for(const auto& seg : bnd)
                {
                	VectorDim temp_shift=seg->patch->shift;
                	VectorDim P0_temp=seg->meshIntersection->P0-temp_shift;
                	VectorDim P1_temp=seg->meshIntersection->P1-temp_shift;
                    nA+= 0.5*(P0_temp-P0).cross(P1_temp-P0_temp);
                }
                return nA;
            }
            else
            {
                return VectorDim::Zero();
            }
        }
        
        /**********************************************************************/
        double rightHandedArea(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
        {
            return rightHandedNormal(bnd).dot(referencePlane->unitNormal);
        }

        /**********************************************************************/
        bool isRightHandedBoundary(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
        {
            return rightHandedArea(bnd)>=0.0;
        }
        

        /**********************************************************************/
        VectorLowerDim getLocalPosition(const VectorDim& point) const
        {
            const VectorDim pointLocal(L2G.transpose()*(point-referencePlane->P));
            assert(fabs(pointLocal(2))<FLT_EPSILON);
            return pointLocal.template segment<2>(0);
        }
        
        
        /**********************************************************************/
        std::shared_ptr<PeriodicPlaneNode<dim> > getSharedNode(const VectorDim& pointDim)
        {
            const VectorLowerDim point(getLocalPosition(pointDim));
            const auto iter(nodes().find(point));
            if(iter==nodes().end())
            {// point does not exist
                std::shared_ptr<PeriodicPlaneNode<dim>> newNode(new PeriodicPlaneNode<dim>(point));
                nodes().emplace(point,newNode.get());
                return newNode;
            }
            else
            {
                assert(iter->second->patchConnectivities().size()>0 && "EXISTING NODE IS ISOLATED");
                const auto patchConnectivity(iter->second->patchConnectivities().begin()->second);
                if(patchConnectivity.outEdge)
                {
                    return patchConnectivity.outEdge->source;
                }
                else
                {
                    if(patchConnectivity.inEdge)
                    {
                        return patchConnectivity.inEdge->sink;
                    }
                    else
                    {
                        assert(false && "EXISTING POSITION IS NOT CONNECTED");
                        return std::shared_ptr<PeriodicPlaneNode<dim> >(nullptr);
                    }
                }
            }
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class PeriodicGlidePlane : public PeriodicGlidePlaneBase<dim>
    /*                      */,public KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>>
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>> PatchContainerType;

    public:


        PeriodicGlidePlane(GlidePlaneFactory<dim>& glidePlaneFactory_in,
                               const GlidePlaneKey<dim>& referencePlaneKey,
                               const VectorDim& globalX) :
        /* init */ PeriodicGlidePlaneBase<dim>(glidePlaneFactory_in,referencePlaneKey,globalX)
        {
            getPatch(VectorDim::Zero());
        }

//        ~PeriodicGlidePlane()
//        {
//            std::cout<<cyanColor<<"Destroying PeriodicGlidePlane"<<defaultColor<<std::endl;
//        }

        /**********************************************************************/
        const PatchContainerType& patches() const
        {
            return *this;
        }

        PatchContainerType& patches()
        {
            return *this;
        }

        /**********************************************************************/
        VectorDim getShift(const PeriodicPlaneEdge<dim>& edge) const
        {

            VectorDim shift(VectorDim::Zero());
            for(const PlanarMeshFace<dim>* const face : edge.meshIntersection->faces)
            {
                const auto parallelFaceID(edge.patch->glidePlane->grain.region.parallelFaces().at(face->sID));
                const auto parallelFace(edge.patch->glidePlane->grain.region.faces().at(parallelFaceID));
                shift+=(parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal();

            }
            return shift;
        }

        void fill_unfilled_patch()
        {
        		assert(0 && "Unfilled Patch detected");
        }

        /**********************************************************************/
        template<typename NodeType>
//        void addPatchesContainingPolygon(const std::vector<std::shared_ptr<NodeType>>& polyPoints)
        void addPatchesContainingPolygon(const std::vector<NodeType>& polyPoints)
        {

            //const double dt(1.0);

            std::cout<<"this->referencePlane->meshIntersections.size()="<<this->referencePlane->meshIntersections.size()<<std::endl;
                // Compute a reference point inside the outerBoundary
                VectorLowerDim insideReferencePoint(VectorLowerDim::Zero());
                for(const auto& seg : this->referencePlane->meshIntersections)
                {
                    insideReferencePoint+=this->getLocalPosition(seg->P0);
                }
                insideReferencePoint/=this->referencePlane->meshIntersections.size();
                assert(this->isInsideOuterBoundary(insideReferencePoint));

                for(const auto& polyPoint : polyPoints)
                {
                    std::cout<<"Finding node "<<polyPoint.sID<<std::endl;
                    const VectorLowerDim P(this->getLocalPosition(polyPoint.P));
                    //                const VectorLowerDim P(getLocalPosition(polyPoint.P+polyPoint.V*dt));
                    while(!this->isInsideOuterBoundary(P))
                    {
                        std::cout<<"point "<<polyPoint.sID<<" outside"<<std::endl;

//                        const std::vector<const PeriodicPlaneEdge<dim>*> currentOuterBnd(this->outerBoundaries()[0]);//Should not this be a multiple loop case
                        for (const auto& currentOuterBnd : this->outerBoundaries())
                        {

                        		for(const auto& bndEdge : currentOuterBnd)
                        		{
                        			std::cout<<"intersecting edge "<<bndEdge->tag()<<std::endl;

                        			SegmentSegmentDistance<dim-1> ssd(P,insideReferencePoint,*bndEdge->source,*bndEdge->sink);
                        			if(ssd.dMin<FLT_EPSILON)
                        			{// intersection with current boundary found
                        				std::cout<<"intersection found "<<std::endl;

                        				getPatch(getShift(*bndEdge)+bndEdge->patch->shift);
                        				break;
                        			}

                        	}
                        		if (this->innerBoundaries().size())
                        		{
//                        			assert(0 && "inside patch detected");
                        			for (const auto& currinBnds :  this->innerBoundaries())
                        			{
                        				for(const auto& inbndEdge : currinBnds)
                        				{
                        					getPatch(getShift(*inbndEdge)+inbndEdge->patch->shift);
                        					break;
                        				}
//                        				assert(0 && "inside patch detected");

                        			}

                        		}


                        }

                    }
                    std::cout<<"point "<<polyPoint.sID<<" inside"<<std::endl;


                }
                std::cout<<"Inserting the patches corresponding to the diagonal positions "<<std::endl;
                for(size_t i=0;i<polyPoints.size();i++)
                {
                	size_t j=i+1;
                	j=(j==polyPoints.size()? 0 : j);
                	std::cout<<"node size "<<polyPoints.size()<<std::endl;
//                	std::cout<<"Nodes connected are "<<polyPoints[i].sID<<"    ----->     "<<polyPoints[j].sID<<"\n";
                	bool seg_inside=false;
        			while(!seg_inside)
        			{
        				const std::vector<const PeriodicPlaneEdge<dim>*> currentOuterBnd(this->outerBoundaries()[0]);
    					seg_inside=true;
        				for(const auto& bndEdge : currentOuterBnd)
        				{
        					std::cout<<"intersecting edge "<<bndEdge->tag()<<std::endl;


        					SegmentSegmentDistance<dim-1> ssd(this->getLocalPosition(polyPoints[i].P),this->getLocalPosition(polyPoints[j].P),*bndEdge->source,*bndEdge->sink);
        					if(ssd.dMin<FLT_EPSILON)
        					{// intersection with current boundary found
        						std::cout<<"intersection found "<<std::endl;
        						seg_inside=false;
        						getPatch(getShift(*bndEdge)+bndEdge->patch->shift);
        						break;
        					}

        				}
                		if (this->innerBoundaries().size())
                		{
//                			assert(0 && "inside patch detected");
                			for (const auto& currinBnds :  this->innerBoundaries())
                			{
                				for(const auto& inbndEdge : currinBnds)
                				{
                					getPatch(getShift(*inbndEdge)+inbndEdge->patch->shift);
                					break;
                				}
//                				assert(0 && "inside patch detected");

                			}

                		}

        			}

                }






        }


        /**********************************************************************/
        void print()
        {
            std::cout<<"PeriodiPlane nodes:"<<std::endl;
            for(const auto& node : this->nodes())
            {
                std::cout<<node.second->sID<<" "<<node.first.transpose()<<std::endl;
            }

            std::cout<<"PeriodiPlane patches:"<<std::endl;
            for(const auto& patch : patches())
            {
                std::cout<<"Patch "<<patch.shift.transpose()<<std::endl;
                for(const auto& edge : patch.second.edges())
                {
                    std::cout<<edge->source->sID<<"->"<<edge->sink->sID<<std::endl;
                }
            }


        }

        /**********************************************************************/
        std::shared_ptr<PeriodicPlanePatch<dim>> getPatch(const VectorDim& shift)
        {
            return patches().get(shift);
        }



    };

    template<int dim>
    struct PeriodicPlanePatchIO
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        size_t glidePlaneID;
        size_t referencePlaneID;
        VectorDim shift;
        
        PeriodicPlanePatchIO() :
        /* init */ glidePlaneID(0)
        /* init */,referencePlaneID(0)
        /* init */,shift(VectorDim::Zero())
        {
            
        }
        
        PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch) :
        /* init */ glidePlaneID(patch.glidePlane->sID)
        /* init */,referencePlaneID(patch.periodicPlane->referencePlane->sID)
        /* init */,shift(patch.shift)
        {
            
        }
        
        /**********************************************************************/
        PeriodicPlanePatchIO(std::stringstream& ss) :
        /* init */ glidePlaneID(0)
        /* init */,referencePlaneID(0)
        /* init */,shift(VectorDim::Zero())
        {
            ss>>glidePlaneID;
            ss>>referencePlaneID;
            for(int d=0;d<dim;++d)
            {
                ss>>shift(d);
            }
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds)
        {
            os  << ds.glidePlaneID<<"\t"
            /**/<< ds.referencePlaneID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.shift.transpose();
            return os;
        }
        
    };
    
    /**********************************************************************/
    template<int dim>
    PeriodicPlaneEdge<dim>::PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                                              const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in) :
    /* init */ patch(patch_in)
    /* init */,source(source_in)
    /* init */,sink(sink_in)
    /* init */,meshIntersection(meshIntersection_in)
    /* init */,next(nullptr)
    /* init */,prev(nullptr)
    /* init */,twin(nullptr)
    {
//        std::cout<<"Creating PeriodicPlaneEdge "<<this->tag()<<std::endl;
        source->addLink(this);
        sink->addLink(this);
        if(twin)
        {
            patch->periodicPlane->removeUntwinnedEdge(twin);
        }
        else
        {
            patch->periodicPlane->addUntwinnedEdge(this);
        }
    }
    
    /**********************************************************************/
    template<int dim>
    PeriodicPlaneEdge<dim>::~PeriodicPlaneEdge()
    {
//        std::cout<<"Destroying PeriodicPlaneEdge "<<this->tag()<<std::endl;

        source->removeLink(this);
        sink->removeLink(this);
        if(next)
        {
            next->prev=nullptr;
        }
        if(prev)
        {
            prev->next=nullptr;
        }
        if(twin)
        {
            twin->twin=nullptr;
            patch->periodicPlane->addUntwinnedEdge(twin);
        }
        else
        {
            patch->periodicPlane->removeUntwinnedEdge(this);
        }
    }
    
    /**********************************************************************/
    template<int dim>
    std::string PeriodicPlaneEdge<dim>::tag() const
    {/*!\returns the string "i->j" where i is source()->sID and j=sink()->sID
      */
        return std::to_string(source->sID) + "->" + std::to_string(sink->sID)+" ("+std::to_string(patch->sID)+")";
    }
}
#endif
