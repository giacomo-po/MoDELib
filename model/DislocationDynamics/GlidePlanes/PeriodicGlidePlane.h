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
#include <list>

#include <GlidePlane.h>
#include <GlidePlaneFactory.h>
#include <TerminalColors.h>




namespace model
{
    

    
    template<int dim>
    class PeriodicGlidePlane;
    
    
    template<int dim>
    struct PeriodicGlidePlaneFactory;
    
    
    
    template<int dim>
    struct TypeTraits<PeriodicGlidePlaneFactory<dim>>
    {
        typedef PeriodicGlidePlane<dim> ValueType;
        typedef GlidePlaneKey<dim> KeyType;
        typedef std::less<std::array<long int,dim+3>> CompareType;
    };
    
//    template<int dim>
//    struct TypeTraits<PeriodicGlidePlaneFactoryTraits<dim>> : public PeriodicGlidePlaneFactoryTraits<dim>
//    {
//
//    };
    
//    template<int dim>
//    struct TypeTraits<PeriodicGlidePlane<dim>> : public PeriodicGlidePlaneFactoryTraits<dim>
//    {
//
//    };
    
    /**************************************************************************/
    template<int dim>
    struct PeriodicGlidePlaneFactory : public KeyConstructableWeakPtrFactory<PeriodicGlidePlaneFactory<dim>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
        typedef KeyConstructableWeakPtrFactory<PeriodicGlidePlaneFactoryType> BaseType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef PeriodicGlidePlane<dim> PeriodicGlidePlaneType;
        typedef typename GlidePlaneType::KeyType GlidePlaneKeyType;
        typedef std::shared_ptr<PeriodicGlidePlaneType> PeriodicGlidePlaneSharedPtrType;
        
    public:
        
        const Polycrystal<dim>& poly;
        GlidePlaneFactory<dim>& glidePlaneFactory;
        const Eigen::Matrix<long int,dim,dim> N;
        
        PeriodicGlidePlaneFactory(const Polycrystal<dim>& poly_in,
                                GlidePlaneFactory<dim>& glidePlaneFactory_in) :
        /* init */ poly(poly_in)
        /* init */,glidePlaneFactory(glidePlaneFactory_in)
        /* init */,N(get_N(poly))
        {
            
        }
        
        GlidePlaneKeyType periodicPlaneKey(GlidePlaneKeyType temp) const
        {
            temp[dim+2]=temp[dim+2]%LatticeGCD<dim>::gcd(temp.r.transpose()*N);
            return temp;
        }
        
        /**********************************************************************/
        PeriodicGlidePlaneSharedPtrType get(const GlidePlaneType& plane)
        {
            return BaseType::get(periodicPlaneKey(plane.key));
        }
        
        /**********************************************************************/
        PeriodicGlidePlaneSharedPtrType get(const GlidePlaneKeyType& temp)
        {
            return BaseType::get(periodicPlaneKey(temp));
        }
    
        /**********************************************************************/
        BaseType& periodicGlidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BaseType& periodicGlidePlanes() const
        {
            return *this;
        }
        
        
        /**********************************************************************/
        static Eigen::Matrix<long int,dim,dim> get_N(const Polycrystal<dim>& poly)
        {
            assert(poly.grains().size()==1 && "Periodic simulations only supported in single crystals");
            const auto& grain(poly.grains().begin()->second);
            const auto& meshRegion(grain.region);
            
            Eigen::Matrix<double,dim,dim> B(Eigen::Matrix<double,dim,dim>::Zero());
            int col=0;
            for(const auto& pair : meshRegion.parallelFaces())
            {
                model::cout<<"Checking if parallel faces "<<pair.first<<"<->"<<pair.second<<" are commensurate"<<std::endl;
                const PlanarMeshFace<dim>& face1(*meshRegion.faces().at(pair.first));
                const PlanarMeshFace<dim>& face2(*meshRegion.faces().at(pair.second));
                const VectorDim cc(face1.center()-face2.center());
                const VectorDim ccc(cc.dot(face1.outNormal())*face1.outNormal());
                if(col<dim)
                {
                    B.col(col)=ccc;
                    col++;
                }
                const LatticeDirection<dim> ld(grain.latticeDirection(face1.outNormal()));
                const double normRatio(ccc.norm()/ld.cartesian().norm());
                if(std::fabs(std::round(normRatio)-normRatio)>FLT_EPSILON)
                {
                    //                            std::cout<<"Face outNormal="<<std::setprecision(15)<<std::scientific<<face1.outNormal().transpose()<<std::endl;
                    std::cout<<"Mesh in direction "<< std::setprecision(15)<<std::scientific<<ld.cartesian().normalized().transpose()<<" is not commensurate for periodicity"<<std::endl;
                    std::cout<<"Mesh size in that direction must be a multiple of "<< std::setprecision(15)<<std::scientific<<ld.cartesian().norm()<<std::endl;
                    std::cout<<"Size detected="<< std::setprecision(15)<<std::scientific<<ccc.norm()<<std::endl;
                    std::cout<<"Closest commensurate size="<< std::setprecision(15)<<std::scientific<<std::round(normRatio)*ld.cartesian().norm()<<std::endl;
                    assert(false && "MESH NOT COMMENSURATE");
                }
            }
            
            Eigen::Matrix<double,dim,dim> Nd(grain.latticeBasis.inverse()*B);
            Eigen::Matrix<double,dim,dim> Ni(Nd.array().round().matrix());
            if((Nd-Ni).norm()>FLT_EPSILON)
            {
                assert(0 && "Mesh not commensurate with the lattice");
            }
            return Ni.template cast<long int>();
        }
        
    };
    
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
    /*                     */,public Eigen::Matrix<double,dim-1,1>
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
        InOutEdgeContainerType inEdges() const
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
            return temp;
        }
        
        /**********************************************************************/
        InOutEdgeContainerType outEdges() const
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
            return temp;
        }
        
        //        /**********************************************************************/
        //        PeriodicPlaneEdgeType* inEdge() const
        //        {
        //            InOutEdgeContainerType temp;
        //            for(const auto& connectivity : neighborConnectivities())
        //            {
        //                if(connectivity.second.inEdge && !connectivity.second.outEdge)
        //                {
        //                    assert(connectivity.second.inEdge->twin==nullptr);
        //                    temp.insert(connectivity.second.inEdge);
        //                }
        //            }
        //            return temp.size()==1? *temp.begin() : nullptr;
        //        }
        //
        //        /**********************************************************************/
        //        PeriodicPlaneEdgeType* outEdge() const
        //        {
        //            InOutEdgeContainerType temp;
        //            for(const auto& connectivity : neighborConnectivities())
        //            {
        //                if(connectivity.second.outEdge && !connectivity.second.inEdge)
        //                {
        //                    assert(connectivity.second.outEdge->twin==nullptr);
        //                    temp.insert(connectivity.second.outEdge);
        //                }
        //            }
        //            return temp.size()==1? *temp.begin() : nullptr;
        //        }
        
        
        
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlanePatch : public StaticID<PeriodicPlanePatch<dim>>
    /*                      */, private std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;

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
        
        /**********************************************************************/
        int contains(const VectorLowerDim& test)
        {
//                const auto& outerBoundary(outerBoundaries()[0]);
                size_t i, j, c = 0;
                for (i = 0, j = edges().size()-1; i < edges().size(); j = i++)
                {
                    const VectorLowerDim& Pi(*edges()[i]->source);
                    const VectorLowerDim& Pj(*edges()[j]->source);
                    
                    if ( ((Pi(1)>test(1)) != (Pj(1)>test(1))) &&
                        (test(0) < (Pj(0)-Pi(0)) * (test(1)-Pi(1)) / (Pj(1)-Pi(1)) + Pi(0)) )
                    {
                        c = !c;
                    }
                }
                return c;
        }
        
        
    };
    
    
    
    
    
    template<int dim>
    struct PeriodicGlidePlaneBase : public StaticID<PeriodicGlidePlaneBase<dim>>
    /*                           */,private std::map<Eigen::Matrix<double,dim-1,1>,const std::weak_ptr<PeriodicPlaneNode<dim>>,CompareVectorsByComponent<double,dim-1,float>>
    /*                           */,private std::set<const PeriodicPlaneEdge<dim>*>
    
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef std::map<Eigen::Matrix<double,dim-1,1>,const std::weak_ptr<PeriodicPlaneNode<dim>>,CompareVectorsByComponent<double,dim-1,float>> NodeCointainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
        
        BoundariesContainerType _outerBoundaries;
        BoundariesContainerType _innerBoundaries;
        
        
        /**********************************************************************/
        static MatrixDim getL2G(VectorDim z)
        {
//            const double xNorm(x.norm());
            const double zNorm(z.norm());
            assert(zNorm>FLT_EPSILON);
            z/=zNorm;
            
            VectorDim x(VectorDim::UnitX().cross(z));
            double xNorm(x.norm());
            if(xNorm>FLT_EPSILON)
            {
                x=x/xNorm;
            }
            else
            {
                x=VectorDim::UnitY().cross(z);
                xNorm=x.norm();
                if(xNorm>FLT_EPSILON)
                {
                    x=x/xNorm;
                }
                else
                {
                    x=VectorDim::UnitZ().cross(z);
                    xNorm=x.norm();
                    if(xNorm>FLT_EPSILON)
                    {
                        x=x/xNorm;
                    }
                    else
                    {
                        assert(false && "CANNOT FIND VECTOR ORTHOGONAL TO z");
                    }
                }
            }

            assert(std::fabs(x.norm()-1.0)<FLT_EPSILON);
            assert(std::fabs(z.norm()-1.0)<FLT_EPSILON);
            assert(fabs(x.dot(z)<FLT_EPSILON));
            MatrixDim temp(Eigen::Matrix3d::Identity());
            temp.col(2)=z;
            temp.col(0)=x;
            temp.col(1)=temp.col(2).cross(temp.col(0));
            return temp;
        }
        
        GlidePlaneFactory<dim>& glidePlaneFactory;
        const std::shared_ptr<GlidePlane<dim>> referencePlane;
        const MatrixDim L2G;
        
        PeriodicGlidePlaneBase(GlidePlaneFactory<dim>& glidePlaneFactory_in,
                               const GlidePlaneKey<dim>& referencePlaneKey) :
        /* init */ glidePlaneFactory(glidePlaneFactory_in)
        /* init */,referencePlane(glidePlaneFactory.get(referencePlaneKey))
        /* init */,L2G(getL2G(referencePlane->unitNormal))
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
            const VectorDim relP(point-referencePlane->P);
            const VectorDim pointLocal(L2G.transpose()*(relP-referencePlane->unitNormal*referencePlane->unitNormal.dot(relP)));
            if(fabs(pointLocal(2))>FLT_EPSILON)
            {
                std::cout<<"L2G=\n"<<L2G<<std::endl;
                std::cout<<"pointLocal"<<pointLocal.transpose()<<std::endl;
                assert(false && "local point has non-zero z-coordinate");
            }
            return pointLocal.template segment<2>(0);
        }
        
        /**********************************************************************/
        VectorDim getGlobalPosition(const VectorLowerDim& point) const
        {// terurns the position on the plane in global goordinates
            return L2G.template block<dim,dim-1>(0,0)*point+referencePlane->P;
        }
        
        
        /**********************************************************************/
        std::shared_ptr<PeriodicPlaneNode<dim> > getSharedNode(const VectorDim& pointDim)
        {
            const VectorLowerDim point(getLocalPosition(pointDim));
            const auto iter(nodes().find(point));
            if(iter==nodes().end())
            {// point does not exist
//                std::shared_ptr<PeriodicPlaneNode<dim>> newNode(new PeriodicPlaneNode<dim>(point));
//                nodes().emplace(point,newNode.get());
                return nodes().emplace(point,std::shared_ptr<PeriodicPlaneNode<dim> >(new PeriodicPlaneNode<dim>(point))).first->second.lock();
//                return newNode;
            }
            else
            {// point exists
                if(iter->second.expired())
                {// node deleted elsewhere
                    nodes().erase(iter);
                    return nodes().emplace(point,std::shared_ptr<PeriodicPlaneNode<dim> >(new PeriodicPlaneNode<dim>(point))).first->second.lock();
                }
                else
                {
                    return iter->second.lock();
                }
                
//                assert(iter->second->patchConnectivities().size()>0 && "EXISTING NODE IS ISOLATED");
//                const auto patchConnectivity(iter->second->patchConnectivities().begin()->second);
//                if(patchConnectivity.outEdge)
//                {
//                    return patchConnectivity.outEdge->source;
//                }
//                else
//                {
//                    if(patchConnectivity.inEdge)
//                    {
//                        return patchConnectivity.inEdge->sink;
//                    }
//                    else
//                    {
//                        assert(false && "EXISTING POSITION IS NOT CONNECTED");
//                        return std::shared_ptr<PeriodicPlaneNode<dim> >(nullptr);
//                    }
//                }
            }
        }
        
    };
    
    template<int dim>
    struct TypeTraits<PeriodicGlidePlane<dim>>
    {
        typedef Eigen::Matrix<double,dim,1> KeyType;
        typedef PeriodicPlanePatch<dim> ValueType;
        typedef CompareVectorsByComponent<double,dim,float> CompareType;
    };

    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class PeriodicGlidePlane : public PeriodicGlidePlaneBase<dim>
    /*                      */,public KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>> // container of patches
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>> PatchContainerType;
        typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
        typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
        typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
        typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef typename GlidePlaneType::KeyType GlidePlaneKeyType;

    public:
        
        PeriodicGlidePlaneFactoryType& periodicGlidePlaneFactory;

        
        /**********************************************************************/
        PeriodicGlidePlane(PeriodicGlidePlaneFactoryType& pgpf,
                           const GlidePlaneKeyType& key_in) :
        /* init */ PeriodicGlidePlaneBase<dim>(pgpf.glidePlaneFactory,key_in)
        /* init */,periodicGlidePlaneFactory(pgpf)
        {
            
        }
        
//        PeriodicGlidePlane(GlidePlaneFactory<dim>& glidePlaneFactory_in,
//                           const GlidePlaneKey<dim>& referencePlaneKey,
//                           const VectorDim& globalX) :
//        /* init */ PeriodicGlidePlaneBase<dim>(glidePlaneFactory_in,referencePlaneKey,globalX)
//        {
//            getPatch(VectorDim::Zero());
//        }
        
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
        void createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy)
        {
//            std::cout<<"createNewBoundary"<<std::endl;
            BoundaryContainerType temp;
            temp.reserve(untwinnedCopy.size());
            while(true)
            {
//                std::cout<<"currentEdge "<<currentEdge->tag()<<std::endl;
                temp.push_back(currentEdge);
                const size_t erased(untwinnedCopy.erase(currentEdge));
                if(erased!=1)
                {
                    
                    std::cout<<"Trying to erase "<<currentEdge->tag()<<std::endl;
                    std::cout<<"untwinnedCopy is"<<std::endl;
                    for(const auto& edgePtr : untwinnedCopy)
                    {
                        std::cout<<"    "<<edgePtr->tag()<<std::endl;
                    }
                    
                    std::ofstream pointsFile("points.txt");
                    std::cout<<"points"<<std::endl;
                    for(const auto& node : this->nodes())
                    {
                        if(!node.second.expired())
                        {
                            pointsFile<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                            std::cout<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                        }
                    }
                    
                    
                    std::ofstream edgesFile("edges.txt");
                    for(const auto& patch : patches())
                    {
                        for(const auto& edge : patch.second->edges())
                        {
                            edgesFile<<edge->tag()<<std::endl;
                        }
                    }
                    
                    
                    assert(erased==1 && "could not find link in untwinnedEdges 2");
                }
                if(temp.back()->sink->sID==temp.front()->source->sID)
                {
                    break;
                }
                currentEdge=currentEdge->next;
                while(currentEdge->twin)
                {
                    currentEdge=currentEdge->twin->next;
                }
                
            }
            
            if(temp.size())
            {// if temp is not empty a loop was closed.
                if(this->isRightHandedBoundary(temp))
                {// assign temp to _outerBoundaries
                    this->_outerBoundaries.push_back(temp);
                }
                else
                {// assign temp to _innerBoundaries
                    this->_innerBoundaries.push_back(temp);
                    
                }
            }
            //            }
        }
        
        /**********************************************************************/
        void updateBoundaries()
        {
            //            std::cout<<"Updating OuterBoundary"<<std::endl;
            this->outerBoundaries().clear();
            this->innerBoundaries().clear();
            UntwinnedEdgeContainerType untwinnedCopy(this->untwinnedEdges());
            while(untwinnedCopy.size())
            {
                for(const auto& edge : untwinnedCopy)
                {
                    if(edge->source->outEdges().size()==1 && edge->source->inEdges().size()==1)
                    {// must start from node with only one edge in and one edge out
                        createNewBoundary(edge,untwinnedCopy);
                        break;
                    }
                }
            }
            
//            std::cout<<"OuterBoundaries #"<<this->outerBoundaries().size()<<std::endl;
            //            for(const auto& bnd : outerBoundaries())
            //            {
            //                std::cout<<"size="<<bnd.size()<<", area="<<rightHandedArea(bnd)<<std::endl;
            //                for (const auto bnd2 : bnd)
            //                std::cout<<bnd2->tag()<<std::endl;
            //            }
            //
//            std::cout<<"InnerBoundaries #"<<this->innerBoundaries().size()<<std::endl;
            //            for(const auto& bnd : innerBoundaries())
            //            {
            //                std::cout<<"size="<<bnd.size()<<", area="<<rightHandedArea(bnd)<<std::endl;
            //                for (const auto bnd2 : bnd)
            //                std::cout<<bnd2->tag()<<std::endl;
            //                assert(0 && "Code in the right Handed Area");
            //            }
            
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
        
        /**********************************************************************/
        void fillHoles()
        {
            while(this->innerBoundaries().size())
            {
                const PeriodicPlaneEdge<dim>* const holeEdge(*this->innerBoundaries().front().begin());
                getPatch(getShift(*holeEdge)+holeEdge->patch->shift);
            }
        }
        
        /**********************************************************************/
//        template<typename NodeType>
        void addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints)
        {
            
            if(polyPoints.size()>=3)
            {
                
                std::set<long int> pointsHeights;
                for(const auto& point : polyPoints)
                {
                    const auto heightPair(LatticePlane::computeHeight(this->referencePlane->n,point));
                    assert(heightPair.first && "Point not on a lattice plane");
                    pointsHeights.insert(heightPair.second);
                }
                assert(pointsHeights.size()==1 && "polyPoints on different planes");
                const GlidePlaneKey<dim> pointsPlaneKey(this->referencePlane->grain.grainID,polyPoints[0],this->referencePlane->n);
                const auto pointsPlane(this->glidePlaneFactory.get(pointsPlaneKey));
                const VectorDim pointsShift(pointsPlane->P-this->referencePlane->P);
                getPatch(pointsShift);
                
                VectorLowerDim insideReferencePoint(VectorLowerDim::Zero());
                for(const auto& seg : pointsPlane->meshIntersections) // problem is here for referencePlane not cutting mesh
                {
                    insideReferencePoint+=this->getLocalPosition(seg->P0);
                }
                insideReferencePoint/=pointsPlane->meshIntersections.size();
                assert(this->isInsideOuterBoundary(insideReferencePoint));

                
//                std::deque<std::shared_ptr<PeriodicPlanePatch<dim>>> tempPatches;
                
                const VectorLowerDim P0(this->getLocalPosition(polyPoints[0]));

                while(!this->isInsideOuterBoundary(P0))
                {
                    std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
                    for(const auto& bndEdge : this->untwinnedEdges())
                    {// loop over outer boundaries and holes
                        SegmentSegmentDistance<dim-1> ssd(P0,insideReferencePoint,*bndEdge->source,*bndEdge->sink);
                        if(ssd.dMin<FLT_EPSILON)
                        {// intersection with current boundary found
                            crossdEdges.insert(bndEdge);
                        }
                    }
                    
                    switch (crossdEdges.size())
                    {
                        case 1:
                        {
                            getPatch(getShift(**crossdEdges.begin())+(*crossdEdges.begin())->patch->shift);
                            break;
                        }
                            
                        case 2:
                        {
                            assert((*crossdEdges.begin())->patch==(*crossdEdges.rbegin())->patch);
                            getPatch(getShift(**crossdEdges.begin())+getShift(**crossdEdges.rbegin())+(*crossdEdges.begin())->patch->shift);
                            break;
                        }
                            
                        default:
                        {
                            std::cout<<"crossdEdges.size()="<<crossdEdges.size()<<std::endl;
                            assert(false && "1 or 2 edges must be crossed");
                            break;
                        }
                    }
                    
                    assert(this->outerBoundaries().size()==1 && "THERE MUST BE ONLY ONE OUTER BOUNDARY");
                }
                
                // erase patches needed to find P0
                std::list<VectorDim> eraseKeys;
                for(const auto& patch : patches())
                {
                    if(!patch.second->contains(P0))
                    {
                        eraseKeys.push_back(patch.first);
                    }
                }
                for(const auto& key : eraseKeys)
                {
                    patches().erase(key);
                }
                
                for(size_t k=0;k<polyPoints.size();++k)
                {
                    const VectorLowerDim startPoint(this->getLocalPosition(polyPoints[k]));
                    const VectorLowerDim endPoint(k==polyPoints.size()-1? this->getLocalPosition(polyPoints[0]) : this->getLocalPosition(polyPoints[k+1]));
                    while(true)
                    {
//                        std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
                        std::map<const PeriodicPlanePatch<dim>*,std::set<const PeriodicPlaneEdge<dim>*>> crossdEdges;
                        for(const auto& bndEdge : this->untwinnedEdges())
                        {// loop over outer boundaries and holes
                            SegmentSegmentDistance<dim-1> ssd(startPoint,endPoint,*bndEdge->source,*bndEdge->sink);
                            if(ssd.dMin<FLT_EPSILON)
                            {// intersection with current boundary found
                                crossdEdges[bndEdge->patch].insert(bndEdge);
                            }
                        }
                        
                        if(crossdEdges.size()==0)
                        {// No untwinned edges have been crossed, break and move to next pair of polyPoints
                            break;
                        }
                        else
                        {
                            for(const auto& pair : crossdEdges)
                            {
                                VectorDim shift(pair.first->shift);
                                for(const auto& edge : pair.second)
                                {
                                    shift+=getShift(*edge);
                                }
                                getPatch(shift); // this will change untwinnedEdges. Remain in while loop

                            }
                        }
                        assert(this->outerBoundaries().size()==1 && "THERE MUST BE ONLY ONE OUTER BOUNDARY");
                    }
                }
                
                fillHoles();
                
                if(true)
                {
                    
                    
                    std::ofstream polyFile("poly.txt");
                    std::ofstream poly3DFile("poly3D.txt");

                    polyFile<<insideReferencePoint.transpose()<<std::endl;
                    for(const auto& node : polyPoints)
                    {
                        polyFile<<"    "<<this->getLocalPosition(node).transpose()<<std::endl;
                        poly3DFile<<node.transpose()<<std::endl;
                    }
                    
                    std::ofstream pointsFile("points.txt");
                    for(const auto& node : this->nodes())
                    {
                        if(!node.second.expired())
                        {
                            pointsFile<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                        }
                    }
                                        
                    std::ofstream edgesFile("edges.txt");
                    for(const auto& patch : patches())
                    {
                        for(const auto& edge : patch.second->edges())
                        {
                            edgesFile<<edge->tag()<<std::endl;
                        }
                    }
                    assert(this->isCompact() && "Plane not compact");
                }
            }
        }
        
        
        /**********************************************************************/
        void print()
        {
            std::cout<<"PeriodiPlane nodes:"<<std::endl;
            for(const auto& node : this->nodes())
            {
                if(!node.second.expired())
                {
                    std::cout<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                }
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
        return std::to_string(source->sID) + " " + std::to_string(sink->sID)+" "+std::to_string(patch->sID);
    }
}
#endif
