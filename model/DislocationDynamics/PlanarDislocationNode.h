/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PLANARDISLOCATIONNODE_H_
#define model_PLANARDISLOCATIONNODE_H_

#include <algorithm>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <tuple>

#include <DislocationNetworkTraits.h>
#include <SplineNode.h>
#include <GramSchmidt.h>
#include <Simplex.h>
#include <LatticeMath.h>
//#include <SimplexBndNormal.h>
#include <Grain.h>
#include <GlidePlane.h>
#include <PlanePlaneIntersection.h>
#include <PlaneLineIntersection.h>
#include <FiniteLineSegment.h>
#include <DislocationNodeIO.h>
//#include <BoundingLineSegments.h>
#include <DefectiveCrystalParameters.h>
#include <ConfinedDislocationObject.h>
#include <DislocationLoopIO.h>


#ifndef NDEBUG
#define VerbosePlanarDislocationNode(N,x) if(verbosePlanarDislocationNode>=N){std::cout<<x;}
#else
#define VerbosePlanarDislocationNode(N,x)
#endif

namespace model
{
    
    template <typename Derived,typename InterpolationType>
    class PlanarDislocationNode : public SplineNode<Derived,TypeTraits<Derived>::dim,TypeTraits<Derived>::corder,InterpolationType>
    /*                         */,public ConfinedDislocationObject<TypeTraits<Derived>::dim>
    {
        
    public:
        
        constexpr static int dim=TypeTraits<Derived>::dim; // make dim available outside class
        constexpr static int corder=TypeTraits<Derived>::corder; // make dim available outside class
        typedef Derived NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
        typedef ConfinedDislocationObject<dim> ConfinedDislocationObjectType;
        typedef typename NodeBaseType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef typename TypeTraits<NodeType>::LoopNetworkType LoopNetworkType;
        constexpr static int NdofXnode=NodeBaseType::NdofXnode;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef typename TypeTraits<NodeType>::MeshLocation MeshLocation;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef typename NodeBaseType::NeighborContainerType NeighborContainerType;
        typedef typename NodeBaseType::LoopLinkContainerType LoopLinkContainerType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        static const double bndTol;
        static int verbosePlanarDislocationNode;
        
    private:
        
        
        
        /**********************************************************************/
        VectorDim snapToBoundingBox(const VectorDim& P) const
        {/*!\param[in] P position to be snapped to the bounding box
          * \returns a point on the bounding box close to P. The returned point
          * is the closest to the bounding box, unless the closest point causes
          * boundarySegments to become interior. In that case the closest boundary
          * vertex is returned.
          */
            
            VerbosePlanarDislocationNode(4,"snapping P="<<P.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            
            
            std::multimap<double,VectorDim> snapMap;
            
            // Collect possible snap points, sorted by distance to P
            //            for(size_t k=0;k<boundingBoxSegments().size();++k)
            for(const auto& seg : this->boundingBoxSegments())
            {
                
                snapMap.emplace((P-seg.P0).squaredNorm(),seg.P0);
                snapMap.emplace((P-seg.P1).squaredNorm(),seg.P1);
                const VectorDim x(seg.snap(P));
                snapMap.emplace((P-x).squaredNorm(),x);
                
                
                //                const auto& vertexPair(boundingBoxSegments()[k]);
                //                const VectorDim segm(vertexPair.P1-vertexPair.P0);
                //                const double segmNorm2(segm.squaredNorm());
                //                if(segmNorm2>FLT_EPSILON)
                //                {
                //                    snapMap.emplace((P-vertexPair.P0).squaredNorm(),vertexPair.P0);
                //                    snapMap.emplace((P-vertexPair.P1).squaredNorm(),vertexPair.P1);
                //
                //                    double u((P-vertexPair.P0).dot(segm)/segmNorm2);
                //                    if(u>0.0 && u<1.0)
                //                    {
                //                        const VectorDim x(vertexPair.P0+u*segm);
                //                        snapMap.emplace((P-x).squaredNorm(),x);
                //                    }
                //                }
                //                else
                //                {
                //                    const VectorDim x(0.5*(vertexPair.P1+vertexPair.P0));
                //                    snapMap.emplace((P-x).squaredNorm(),x);
                //                }
            }
            
            VerbosePlanarDislocationNode(4,"there are "<<snapMap.size()<<" possible snap points."<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            
            
            // Return the first point to which we can snap
            VectorDim   snapPoint(P);
            bool success(false);
            for(const auto& pair : snapMap)
            {
                VerbosePlanarDislocationNode(4,"Checking snap point "<<pair.second.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
                if(isMovableTo(pair.second))
                {
                    VerbosePlanarDislocationNode(4,"Snapping to "<<pair.second.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
                    snapPoint=pair.second;
                    success=true;
                    break;
                    //                    return pair.second;
                }
            }
            
            assert(success && "snapToBoundingBox FAILED.");
            return snapPoint;
        }
        
        
        
    public:
        
        
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        //! The current velocity vector of *this PlanarDislocationNode
        VectorDofType velocity;
        //! The previous velocity vector of *this PlanarDislocationNode
        VectorDofType vOld;
        double velocityReductionCoeff;
        //! The normal unit vector of the boundary on which *this PlanarDislocationNode is moving on
        std::shared_ptr<NodeType> virtualNode;
        
        
        const NodeType* const masterNode;
        const std::set<const PlanarMeshFace<dim>*> masterSegmentFaces;

        std::map<std::set<size_t>,std::shared_ptr<NodeType>> imageNodeContainer;

//        const bool isVirtualBoundaryNode;
        
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            use_velocityFilter=TextFileParser(fileName).readScalar<double>("use_velocityFilter",true);
            velocityReductionFactor=TextFileParser(fileName).readScalar<double>("velocityReductionFactor",true);
            assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
            verbosePlanarDislocationNode=TextFileParser(fileName).readScalar<int>("verbosePlanarDislocationNode",true);
        }
        
        
        /**********************************************************************/
        PlanarDislocationNode(LoopNetworkType* const ln,
                              const VectorDim& Pin,
                              const VectorDofType& Vin,
                              const double& vrc) :
        /* base */ NodeBaseType(ln,Pin)
        /* base */,ConfinedDislocationObjectType(this->get_P())
        /* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
        /* init */,velocity(Vin)
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(vrc)
        /* init */,masterNode(nullptr)
//        /* init */,isVirtualBoundaryNode(false)
        {/*! Constructor from DOF
          */
            VerbosePlanarDislocationNode(1,"Creating PlanarDislocationNode "<<this->sID<<" from position"<<std::endl;);
        }
        
        /**********************************************************************/
        PlanarDislocationNode(const LinkType& pL,
                              const double& u) :
        /* init */ NodeBaseType(pL.loopNetwork,pL.get_r(u))
        /* base */,ConfinedDislocationObjectType(this->get_P())
        /* init */,p_Simplex(get_includingSimplex(this->get_P(),pL.source->includingSimplex()))
        /* init */,velocity((pL.source->velocity+pL.sink->velocity)*0.5) // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(std::min(pL.source->velocityReduction(),pL.sink->velocityReduction()))
        /* init */,masterNode(nullptr)
//        /* init */,isVirtualBoundaryNode(false)
        {/*! Constructor from ExpandingEdge and DOF
          */
            VerbosePlanarDislocationNode(1,"Creating PlanarDislocationNode "<<this->sID<<" from expanding "<<pL.source->sID<<"->"<<pL.sink->sID<<std::endl;);
            //            if(pL.isBoundarySegment())
            //            {
            //                assert(isBoundaryNode());
            //            }
            
        }
        
        /**********************************************************************/
        PlanarDislocationNode(LoopNetworkType* const ln,
                              const VectorDim& Pin,
                              const NodeType* const master) :
        /* base */ NodeBaseType(ln,Pin)
        /* base */,ConfinedDislocationObjectType(this->get_P())
//        /* init */,p_Simplex(this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC? get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL) : NULL)
        /* init */,p_Simplex(NULL)
        /* init */,velocity(VectorDim::Zero())
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(1.0)
        /* init */,masterNode(master)
//        /* init */,isVirtualBoundaryNode(true)
        {/*! Constructor from DOF
          */
            VerbosePlanarDislocationNode(1,"Creating VirtualPlanarDislocationNode "<<this->sID<<" (of master "<<masterNode->sID<<")"<<std::endl;);
        }
        
        /**********************************************************************/
        PlanarDislocationNode(LoopNetworkType* const ln,
                              const NodeType* const master,
                              const LinkType* const masterSegment
//                              const std::set<const PlanarMeshFace<dim>*>& masterSegmentFaces
                              ) :
        /* base */ NodeBaseType(ln,imagePosition(master,masterSegment->meshFaces()))
        /* base */,ConfinedDislocationObjectType(this->get_P())
        //        /* init */,p_Simplex(this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC? get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL) : NULL)
        /* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
        /* init */,velocity(VectorDim::Zero())
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(1.0)
        /* init */,masterNode(master)
        /* init */,masterSegmentFaces(masterSegment->meshFaces())
//        /* init */,isVirtualBoundaryNode(true)
        {/*! Constructor from DOF
          */
            VerbosePlanarDislocationNode(1,"Creating ImagePlanarDislocationNode "<<this->sID<<" (of master "<<masterNode->sID<<")"<<std::endl;);
            
            assert(this->network().mesh.regions().size()==1);
            const auto& region(*this->network().mesh.regions().begin()->second);

            for(const auto& face : masterSegmentFaces)
            {
//                const PlanarMeshFace<dim>& currentFace(**masterSegmentFaces.begin());
//                const PlanarMeshFace<dim>& imageFace(*);
                this->meshFaces().insert(region.faces().at(region.parallelFaces().at(face->sID)).get());
            }
            this->updateConfinement();
        }
        
        /**********************************************************************/
        ~PlanarDislocationNode()
        {
            VerbosePlanarDislocationNode(1,"Destroying PlanarDislocationNode "<<this->sID<<" ("<<this<<")"<<std::endl;);
            //            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<", virtual node count="<<virtualNode.use_count()<<std::endl;);
            //            if(virtualNode)
            //            {
            //                this->network().remove(virtualNode->sID);
            //            }
        }
        
        
        bool isVirtualBoundaryNode() const
        {
            return masterNode && masterSegmentFaces.size()==0;
        }
        
        bool isImageNode() const
        {
            return masterNode && masterSegmentFaces.size()>0;
        }
        
        
        static VectorDim imagePosition(const NodeType* const master,
                                       const std::set<const PlanarMeshFace<dim>*>& masterSegmentFaces)
        {
            assert(master->network().mesh.regions().size()==1);
            const auto& region(*master->network().mesh.regions().begin()->second);

            
            switch(masterSegmentFaces.size())
            {
                    
                case 1:
                {
                    
                    const PlanarMeshFace<dim>& currentFace(**masterSegmentFaces.begin());
                    const PlanarMeshFace<dim>& imageFace(*region.faces().at(region.parallelFaces().at(currentFace.sID)));
                    return imageFace.asPlane().snapToPlane(master->get_P());
                    break;
                }
                    
                case 2:
                {
                    return VectorDim::Zero();

                    break;
                }
                    
                default:
                {
                    assert(false && "meshFaces must be 1 or 2");
                    return VectorDim::Zero();
                    break;
                }
                    
                    
            }
        }
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const VectorDim& X,const Simplex<dim,dim>* const guess) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if (guess==NULL)
            {
                temp=this->network().mesh.search(X);
            }
            else
            {
                const auto grains(this->grains());
                if(grains.size()==1)
                {// node only in one region
                    if((*grains.begin())->grainID!=guess->region->regionID)
                    {
                        temp=this->network().mesh.searchRegion((*grains.begin())->grainID,X);
                    }
                    else
                    {
                        temp=this->network().mesh.searchRegionWithGuess(X,guess);
                    }
                }
                else
                {
                    temp=this->network().mesh.searchWithGuess(X,guess);
                }
            }
            if(!temp.first) // PlanarDislocationNode not found inside mesh
            {// Detect if the PlanarDislocationNode is sligtly outside the boundary
                int faceID;
                const double baryMin(temp.second->pos2bary(X).minCoeff(&faceID));
                const bool isApproxOnBoundary(std::fabs(baryMin)<1.0e3*FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
                if(!isApproxOnBoundary)
                {
                    model::cout<<"PlanarDislocationNode "<<this->sID<<" @ "<<X.transpose()<<std::endl;
                    std::cout<<"Simplex "<<temp.second->xID<<std::endl;
                    model::cout<<"bary "<<temp.second->pos2bary(X)<<std::endl;
                    std::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                    model::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                    model::cout<<"face of barymin is region-boundary Simplex? "<<temp.second->child(faceID).isRegionBoundarySimplex()<<std::endl;
                    assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
                }
            }
            return temp.second;
        }
        
        /**********************************************************************/
        const std::shared_ptr<NodeType>& virtualBoundaryNode() const
        {
            return virtualNode;
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" addLoopLink "<<pL->tag()<<std::endl;);
            NodeBaseType::addLoopLink(pL); // forward to base class
            this->addGlidePlane(pL->loop()->glidePlane.get());
        }
        
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          * This functin overrides LoopNode::removeLoopLink
          */
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
            NodeBaseType::removeLoopLink(pL); // forward to base class
            ConfinedDislocationObjectType::clear();
            for(const auto& loopLink : this->loopLinks())
            {
                this->addGlidePlane(loopLink->loop()->glidePlane.get());
            }
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" finished removeLoopLink "<<pL->tag()<<std::endl;);
        }
        
        /**********************************************************************/
        VectorOfNormalsType constraintNormals() const __attribute__ ((deprecated)) // REMOVE THIS FUNCTION AFTER CHANGING REMESH AND DISLOCATION NETWORK COMPONENT
        {
            VectorOfNormalsType temp;
            for(const auto& plane : this->glidePlanes())
            {
                temp.push_back(plane->unitNormal);
            }
            for(const auto& face : this->meshFaces())
            {
                temp.push_back(face->asPlane().unitNormal);
            }
            
            //            if(isGlissile())
            //            {
            //                for(const auto& plane : this->glidePlanes())
            //                {
            //                    temp.push_back(plane->unitNormal);
            //                }
            //                for(const auto& face : this->meshFaces())
            //                {
            //                    temp.push_back(face->asPlane().unitNormal);
            //                }
            ////                temp.push_back(boundaryNormal);
            //                GramSchmidt::orthoNormalize(temp);
            //                assert(temp.size()>=1 && "GLIDING NODE MUST HAVE AT LEAST ONE CONSTRAINT.");
            //            }
            //            else
            //            {
            //                temp.push_back((VectorDim()<<1.0,0.0,0.0).finished());
            //                temp.push_back((VectorDim()<<0.0,1.0,0.0).finished());
            //                temp.push_back((VectorDim()<<0.0,0.0,1.0).finished());
            //            }
            
            return temp;
        }
        
        /**********************************************************************/
        void projectVelocity()
        {
            
            VectorOfNormalsType temp;
            
            for(const auto& loop : this->loops())
            {
                if(loop->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
                {
                    if(loop->glidePlane)
                    {
                        temp.push_back(loop->glidePlane->unitNormal);
                    }
                }
                else if(loop->loopType==DislocationLoopIO<dim>::SESSILELOOP)
                {
                    velocity.setZero();
                    break;
                }
//                else
//                {
//                    velocity.setZero();
//                    break;
//                }
            }
            
            if(this->glidePlanes().size()>=dim)
            {
                velocity.setZero();
            }
            
            if(velocity.squaredNorm()>FLT_EPSILON)
            {
                
                for(const auto& face : this->meshFaces())
                {
                    temp.push_back(face->asPlane().unitNormal);
                }
                
                //                temp.push_back(boundaryNormal);
                //
                //                for(const auto& gb : grainBoundaries())
                //                {
                //                    temp.push_back(gb->unitNormal);
                //                }
                
                GramSchmidt::orthoNormalize(temp);
                
                for(const auto& vec : temp)
                {
                    velocity-=velocity.dot(vec)*vec;
                }
                
            }
            //            }
        }
        
        
        /**********************************************************************/
        void set_V(const VectorDofType& vNew)
        {
            vOld=velocity; // store current value of velocity before updating
            
            //            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
            velocity=vNew;
            projectVelocity();
            
            if(use_velocityFilter)
            {
                const double filterThreshold=0.05*velocity.norm()*vOld.norm();
                
                if(velocity.dot(vOld)<-filterThreshold)
                {
                    velocityReductionCoeff*=velocityReductionFactor;
                }
                else if(velocity.dot(vOld)>filterThreshold)
                {
                    velocityReductionCoeff/=velocityReductionFactor;
                }
                else
                {
                    // don't change velocityReductionCoeff
                }
                if(velocityReductionCoeff>1.0)
                {
                    velocityReductionCoeff=1.0;
                }
                if(velocityReductionCoeff<0.005)
                {
                    velocityReductionCoeff=0.005;
                }
                velocity*=velocityReductionCoeff;
            }
        }
        
        /**********************************************************************/
        const VectorDofType& get_V() const
        {/*! The nodal velocity vector
          */
            return velocity;
        }
        
        /**********************************************************************/
        const Simplex<dim,dim>* includingSimplex() const
        {/*!\returns A pointer to the const Simplex imcluding *this PlanarDislocationNode
          */
            return p_Simplex;
        }
        
        /**********************************************************************/
        MeshLocation meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            MeshLocation temp = MeshLocation::outsideMesh;
            
            if(!isVirtualBoundaryNode())
            {
                if(isBoundaryNode())
                {
                    temp=MeshLocation::onMeshBoundary;
                }
                else
                {
                    if(isGrainBoundaryNode())
                    {
                        temp=MeshLocation::onRegionBoundary;
                    }
                    else
                    {
                        temp=MeshLocation::insideMesh;
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isBoundaryNode() const
        {
            return this->isOnExternalBoundary();
            //
            //            return this->isOnBoundary() ;
        }
        
        /**********************************************************************/
        bool isGrainBoundaryNode() const
        {
            return this->isOnInternalBoundary();
            //
            //            return this->isOnGrainBoundary();
            //            return grainBoundaries().size();
        }
        
        /**********************************************************************/
        bool isPureBoundaryNode() const
        {
            bool temp(isBoundaryNode());
            if(temp)
            {
                for (const auto& neighborIter : this->neighbors())
                {
                    temp*=std::get<1>(neighborIter.second)->isBoundarySegment();
                    if(!temp)
                    {
                        break;
                    }
                }
            }
            return  temp;
        }
        
        /**********************************************************************/
        bool isConnectedToBoundaryNodes() const
        {
            bool temp(true);
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=(std::get<0>(neighborIter.second)->isBoundaryNode() || std::get<1>(neighborIter.second)->hasZeroBurgers());
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isConnectedToGrainBoundaryNodes() const
        {
            bool temp(true);
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=(std::get<0>(neighborIter.second)->isGrainBoundaryNode() || std::get<1>(neighborIter.second)->hasZeroBurgers());
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleBoundaryNode() const
        {
            VerbosePlanarDislocationNode(4,"PlanarDislocationNode "<<this->sID<<" isSimpleBoundaryNode "<<std::flush;);
            bool temp=false;
            if(this->isOnBoundary())
            {
                temp=true; // true if all non-virtual neighbors are boundary
                std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> chordDeq;
                
                for (const auto& neighborIter : this->neighbors())
                {
                    if (!std::get<1>(neighborIter.second)->isVirtualBoundarySegment()
                        //                            !std::get<1>(neighborIter.second)->hasZeroBurgers()
                        )
                    {
                        temp*=std::get<1>(neighborIter.second)->isBoundarySegment();
                        chordDeq.push_back(std::get<1>(neighborIter.second)->chord());
                    }
                }
                
                if(temp && chordDeq.size())
                {
                    for(const auto& chord : chordDeq)
                    {
                        temp*=(chord.cross(chordDeq[0]).squaredNorm()<FLT_EPSILON);
                    }
                }
                else
                {
                    temp=false;
                }
            }
            VerbosePlanarDislocationNode(4,temp<<std::endl;);
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleGrainBoundaryNode() const
        {
            VerbosePlanarDislocationNode(4,"PlanarDislocationNode "<<this->sID<<" isSimpleGrainBoundaryNode "<<std::flush;);
            bool temp=false;
            if(isGrainBoundaryNode())
            {
                if(this->isSimple())
                {
                    temp=true;
                    std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> chordDeq;
                    
                    for (const auto& neighborIter : this->neighbors())
                    {
                        if (!std::get<1>(neighborIter.second)->hasZeroBurgers())
                        {
                            temp*=std::get<1>(neighborIter.second)->isGrainBoundarySegment();
                            chordDeq.push_back(std::get<1>(neighborIter.second)->chord());
                        }
                    }
                    
                    if(temp && chordDeq.size())
                    {
                        for(const auto& chord : chordDeq)
                        {
                            temp*=(chord.cross(chordDeq[0]).squaredNorm()<(FLT_EPSILON*chord.squaredNorm()*chordDeq[0].squaredNorm()));
                        }
                    }
                    else
                    {
                        temp=false;
                    }
                }
            }
            VerbosePlanarDislocationNode(4,temp<<std::endl;);
            return temp;
        }
        
        /**********************************************************************/
        bool isSessileNode() const
        {
            VerbosePlanarDislocationNode(4,"PlanarDislocationNode "<<this->sID<<" isSessileNode "<<std::flush;);
            bool temp=true;
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=std::get<1>(neighborIter.second)->isSessile();
                if(!temp)
                {
                    break;
                }
            }
            VerbosePlanarDislocationNode(4,temp<<std::endl;);
            return temp;
        }
        
        /**********************************************************************/
        bool isPureSessile() const
        {
            bool temp(true);
            for(const auto& loopLink : this->loopLinks())
            {
                temp*=loopLink->loop()->loopType==DislocationLoopIO<dim>::SESSILELOOP;
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isZeroBurgersNode() const
        {
            
            bool temp=true;
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=std::get<1>(neighborIter.second)->hasZeroBurgers();
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleZeroBurgersNode() const
        {
            VerbosePlanarDislocationNode(4,"PlanarDislocationNode "<<this->sID<<" isSimpleZeroBurgersNode "<<std::flush;);
            bool temp(this->isSimple() && isZeroBurgersNode());
            if(temp)
            {// make sure attached sessile segments are aligned
                const LinkType* firstLink(std::get<1>(this->neighbors().begin()->second));
                const LinkType* secondLink(std::get<1>(this->neighbors().rbegin()->second));
                temp*=(firstLink->chord().normalized().cross(secondLink->chord().normalized()).norm()<FLT_EPSILON);
            }
            VerbosePlanarDislocationNode(4,temp<<std::endl;);
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleSessileNode() const
        {
            VerbosePlanarDislocationNode(4,"PlanarDislocationNode "<<this->sID<<" isSimpleSessileNode "<<std::flush;);
            bool temp(this->isSimple() && isSessileNode());
            if(temp)
            {// make sure attached sessile segments are aligned
                const LinkType* firstLink(std::get<1>(this->neighbors().begin()->second));
                const LinkType* secondLink(std::get<1>(this->neighbors().rbegin()->second));
                temp*=(firstLink->chord().normalized().cross(secondLink->chord().normalized()).norm()<FLT_EPSILON);
            }
            VerbosePlanarDislocationNode(4,temp<<std::endl;);
            return temp;
        }
        
        /**********************************************************************/
        NeighborContainerType nonZeroNeighbors() const
        {
            NeighborContainerType temp;
            for (const auto& neighborIter : this->neighbors())
            {
                if(!std::get<1>(neighborIter.second)->hasZeroBurgers())
                {
                    temp.emplace(neighborIter.first,neighborIter.second);
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isRemovable(const double& Lmin,const double& cosRemove) const
        {
            
            
            
            bool temp(   !isVirtualBoundaryNode()
                      && (   isSimpleBoundaryNode()
                          || isSimpleGrainBoundaryNode()
                          || isSimpleSessileNode()
                          || isSimpleZeroBurgersNode()
                          || isGeometricallyRemovable(Lmin,cosRemove)
                          )
                      );
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" isRemovable "<<temp<<std::endl;);
            
            
            //            if(temp && usePeriodic && virtualNode)
            //            {
            //                    temp*=(virtualNode->isSimpleBoundaryNode() ) ;
            //            }
            
            return temp;
        }
        
        
        /**********************************************************************/
        bool isGeometricallyRemovable(const double& Lmin,const double& cosRemove) const
        {
            bool temp=false;
            if(!this->isBoundaryNode() && !this->isGrainBoundaryNode())
            {
                const auto linksMap=this->linksByLoopID();
                if(linksMap.size()==1)
                {
                    const LoopLinkContainerType& linkSet(linksMap.begin()->second);
                    assert(linkSet.size()==2);
                    const LoopLinkType& link0(**linkSet. begin());
                    const LoopLinkType& link1(**linkSet.rbegin());
                    
                    if(link0.loop()->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
                    {
                        const VectorDim chord0(link0.sink()->get_P()-link0.source()->get_P());
                        const VectorDim chord1(link1.sink()->get_P()-link1.source()->get_P());
                        const double chord0Norm(chord0.norm());
                        const double chord1Norm(chord1.norm());
                        
                        if(chord0Norm<Lmin || chord1Norm<Lmin)
                        {
                            if(chord0.dot(chord1)>cosRemove*chord0Norm*chord1Norm)
                            {
                                const VectorDim dv0(link0.sink()->get_V()-link0.source()->get_V());
                                const VectorDim dv1(link1.sink()->get_V()-link1.source()->get_V());
                                if(chord0.dot(dv0)<0.0 || chord1.dot(dv1)<0.0) // at least one of the two segments is getting shorter
                                {
                                    temp=true;
                                }
                            }
                        }
                        
                        if(chord0Norm+chord1Norm<Lmin)
                        {
                            temp=true;
                        }
                    }
                }
            }
            
            return temp;
        }
        
        
        /**********************************************************************/
        const double& velocityReduction() const
        {
            return velocityReductionCoeff;
        }
        
        /**********************************************************************/
        void resetImageNodes()
        {
            
            if(   this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
               || this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
            {
                if(isBoundaryNode() && !masterNode)
                {
//                    assert(this->bndNormal().squaredNorm()>FLT_EPSILON && "BOUNDARY NODE MUST HAVE NON-ZERO NORMAL");
                    
                    
                    for(auto& pair : this->neighbors())
                    {
                        const LinkType& link(*std::get<1>(pair.second));
                        NodeType& otherNode(*std::get<0>(pair.second));

                        if(link.isBoundarySegment())
                        {
                            std::set<size_t> key;
                            for(const auto& face : link.meshFaces())
                            {
                                key.insert(face->sID);
                            }
                            
                            const auto iter(imageNodeContainer.find(key));
                            if(iter==imageNodeContainer.end())
                            {
                                const auto iterPair(imageNodeContainer.emplace(key,new NodeType(&this->network(),this->p_derived(),&link)));
                                iterPair.first->second->resetVirtualBoundaryNode();
                            }
                            else
                            {
                                const VectorDim imageP(imagePosition(iter->second->masterNode,iter->second->masterSegmentFaces));
                                static_cast<NodeBaseType*>(iter->second.get())->set_P(imageP);
                                iter->second->resetVirtualBoundaryNode();
                            }
                            
                            const auto iterOther(otherNode.imageNodeContainer.find(key));
                            if(iterOther==otherNode.imageNodeContainer.end())
                            {
                                const auto iterPair(otherNode.imageNodeContainer.emplace(key,new NodeType(&this->network(),&otherNode,&link)));
                                iterPair.first->second->resetVirtualBoundaryNode();
                            }
                            else
                            {
                                const VectorDim imageP(imagePosition(iterOther->second->masterNode,iterOther->second->masterSegmentFaces));
                                static_cast<NodeBaseType*>(iterOther->second.get())->set_P(imageP);
                                iterOther->second->resetVirtualBoundaryNode();
                            }
                            
                        }
                    }
                }
            }
        }
        
        /**********************************************************************/
        void resetVirtualBoundaryNode()
        {
            
            if(   this->network().simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM
               || this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
               || this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
            {
                if(isBoundaryNode() && !isVirtualBoundaryNode())
                {
                    assert(this->bndNormal().squaredNorm()>FLT_EPSILON && "BOUNDARY NODE MUST HAVE NON-ZERO NORMAL");
                    
                    if(virtualNode)
                    {
                        static_cast<NodeBaseType*>(virtualNode.get())->set_P(this->get_P()+this->network().simulationParameters.virtualSegmentDistance*this->bndNormal());
                    }
                    else
                    {
                        VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
                        virtualNode.reset(new NodeType(&this->network(),this->get_P()+this->network().simulationParameters.virtualSegmentDistance*this->bndNormal(),this->p_derived()));
                    }
                }
            }
            
            
//
//            if(isBoundaryNode() && !isVirtualBoundaryNode)
//            {
//                assert(this->bndNormal().squaredNorm()>FLT_EPSILON && "BOUNDARY NODE MUST HAVE NON-ZERO NORMAL");
//
//
//
//
//                switch (this->network().simulationParameters.simulationType)
//                {
//                    case DefectiveCrystalParameters::FINITE_FEM:
//                    {
//                        if(virtualNode)
//                        {
//                            static_cast<NodeBaseType*>(virtualNode.get())->set_P(this->get_P()+this->network().simulationParameters.virtualSegmentDistance*this->bndNormal());
//                        }
//                        else
//                        {
//                            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
//                            virtualNode.reset(new NodeType(&this->network(),this->get_P()+this->network().simulationParameters.virtualSegmentDistance*this->bndNormal(),this->p_derived()));
//                        }
//                        break;
//                    }
//                    case DefectiveCrystalParameters::PERIODIC:
//                    {
//                        if(virtualNode)
//                        {
//                            static_cast<NodeBaseType*>(virtualNode.get())->set_P(this->get_P() - ((this->network().mesh.xMax()-this->network().mesh.xMin()).cwiseProduct(this->bndNormal())));
//                        }
//                        else
//                        {
//                            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
//                            virtualNode.reset(new NodeType(&this->network(),this->get_P()-((this->network().mesh.xMax()-this->network().mesh.xMin()).cwiseProduct(this->bndNormal())),this->p_derived()));
//                        }
//                        break;
//                    }
//                    default:
//                        break;
//                }
//            }
        }
        
        /**********************************************************************/
        void setToBoundary(const VectorDim& X)
        {
            VerbosePlanarDislocationNode(5,"PlanarDislocationNode "<<this->sID<<" setToBoundary @"<< X.transpose()<<std::endl;);
            NodeBaseType::set_P(X); // in turn this calls PlanarDislocationSegment::updateGeometry, so the boundaryNormal must be computed before this line
            ConfinedDislocationObjectType::updateGeometry(this->get_P());
            VerbosePlanarDislocationNode(5,"containingSegments "<< this->boundingBoxSegments().containingSegments(this->get_P()).size()<<std::endl;);
            p_Simplex=get_includingSimplex(this->get_P(),p_Simplex);
            VerbosePlanarDislocationNode(5,"boundaryNormal "<< this->bndNormal().transpose()<<std::endl;);
            VerbosePlanarDislocationNode(5,"current boundingBox\n "<< this->boundingBoxSegments()<<std::endl;);
        }
        
        /**********************************************************************/
        bool set_P(const VectorDim& newP)
        {
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P to "<<newP.transpose()<<std::endl;);
            // make sure that node is on glide planes
            //            bool glidePlanesContained=true;
            //            for(const auto& gp : glidePlanes())
            //            {
            //                glidePlanesContained*=gp->contains(newP);
            //            }
            
            //            if(glidePlanesContained)
            //            {
            if(this->isOnBoundary() || this->boundingBoxSegments().contains(newP))
            {// node was on bounding box, it must remain on bounding box
                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P, case A "<<std::endl;);
                const VectorDim X(snapToBoundingBox(newP));
                setToBoundary(X);
            }
            else
            {// internal node
                std::pair<bool,const Simplex<dim,dim>*> temp(this->network().mesh.searchRegionWithGuess(newP,p_Simplex));
                if(temp.first)
                {// internal node, and newP is inside current grain
                    if(   (isConnectedToBoundaryNodes() || isConnectedToGrainBoundaryNodes())
                       && this->boundingBoxSegments().size()==2
                       && this->glidePlaneIntersections())
                    {// force special case to boundary to get rid of small debris
                        if((newP-this->glidePlaneIntersections()->P0).norm()<this->network().surfaceAttractionDistance)
                        {
                            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P, case B "<<std::endl;);
                            setToBoundary(this->glidePlaneIntersections()->P0);
                        }
                        else if((newP-this->glidePlaneIntersections()->P1).norm()<this->network().surfaceAttractionDistance)
                        {
                            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P, case C "<<std::endl;);
                            setToBoundary(this->glidePlaneIntersections()->P1);
                        }
                        else
                        {
                            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P, case D "<<std::endl;);
                            NodeBaseType::set_P(newP);
                            ConfinedDislocationObjectType::updateGeometry(this->get_P());
                        }
                    }
                    else
                    {
                        VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P, case E "<<std::endl;);
                        NodeBaseType::set_P(newP);
                        ConfinedDislocationObjectType::updateGeometry(this->get_P());
                    }
                }
                else
                {// internal node, and newP is outside current grain
                    VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<"::set_P. current P="<< this->get_P().transpose()<<"set_P, case F "<<std::endl;);
                    const VectorDim X(snapToBoundingBox(newP));
                    setToBoundary(X);
                }
            }
            
            p_Simplex=get_includingSimplex(this->get_P(),p_Simplex); // update including simplex
            
            resetImageNodes();
            resetVirtualBoundaryNode();
            
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" this->isOnBoundary()="<<this->isOnBoundary()<<std::endl;);
            const double posDelta((this->get_P()-newP).norm());
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" posDelta="<<posDelta<<std::endl;);
            
            return posDelta<FLT_EPSILON;
        }
        
        /**********************************************************************/
        bool isMovableTo(const VectorDim& X) const
        {
            bool isMovable=true;
            
            VerbosePlanarDislocationNode(4,"checking if PlanarDislocationNode "<<this->sID<< " isMovable:"<<std::endl;);
            
            for(const auto& gp : this->glidePlanes())
            {// X must be contained by all glidePlanes
                isMovable*=gp->contains(X);
            }
            VerbosePlanarDislocationNode(4,"  meshPlanes contains X? "<<isMovable<<std::endl;);
            
            if(isMovable)
            {
                for(const auto& pair : this->neighbors())
                {
                    if(std::get<1>(pair.second)->isSessile())
                    {// sessile segments cannot change direction if this node is moved
                        const double currentNorm((std::get<0>(pair.second)->get_P()-this->get_P()).norm());
                        const double newNorm((std::get<0>(pair.second)->get_P()-X).norm());
                        if(currentNorm>FLT_EPSILON && newNorm>FLT_EPSILON)
                        {
                            const bool sessileNeighborMovable=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*currentNorm*newNorm);
                            VerbosePlanarDislocationNode(4,"  sessileNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<sessileNeighborMovable<<std::endl;);
                            isMovable*=sessileNeighborMovable;
                            if(!isMovable)
                            {
                                break;
                            }
                        }
                    }
                }
            }
            
            if(isMovable && this->isOnBoundary())
            {
                
                
                isMovable*=this->boundingBoxSegments().contains(X);
                
                //                if(this->isOnBoundary())
                //                {// preliminarily check that the current bounding box contains X
                //                    // Check that all bounding lines that contain this->get_P() will also contain X
                //                    //                    std::set<const MeshBoundarySegment<dim>*> containingSegments(boundingBoxSegments().containingSegments(this->get_P()));
                //                    //                    for(const auto& seg : containingSegments)
                //                    //                    {
                //                    //                        isMovable*=seg->contains(X);
                //                    //                    }
                //                    //
                //                }
                
                
                
            }
            
            //            if(isMovable)
            //            {
            //
            //
            //
            //                for(const auto& pair : this->neighbors())
            //                {
            ////                    if(std::get<1>(pair.second)->isBoundarySegment())
            ////                    {// boundary segments other than must remain boundary if this node is moved
            ////
            ////                        //                        const bool bndNeighborMovable=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X));
            ////                        //                        VerbosePlanarDislocationNode(4,"  boundaryNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<bndNeighborMovable<<std::endl;);
            ////                        //                        isMovable*=bndNeighborMovable;
            ////
            ////
            ////                        const auto containingSegments(this->boundingBoxSegments().containingSegments(0.5*(std::get<0>(pair.second)->get_P()+this->get_P()))); // bounding box lines containing center of segment
            ////                        for(const auto& seg : containingSegments)
            ////                        {
            ////                            const bool bndNeighborMovable(seg->contains(0.5*(std::get<0>(pair.second)->get_P()+X)));
            ////                            VerbosePlanarDislocationNode(4,"  boundaryNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<bndNeighborMovable<<std::endl;);
            ////                            isMovable*=bndNeighborMovable;
            ////                        }
            ////
            ////                        if(!isMovable)
            ////                        {
            ////                            break;
            ////                        }
            ////
            ////                        //                        const double currentNorm((std::get<0>(pair.second)->get_P()-this->get_P()).norm());
            ////                        //                        const double newNorm((std::get<0>(pair.second)->get_P()-X).norm());
            ////                        //                        const bool bndNeighborStraight=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*currentNorm*newNorm);
            ////                        //                        VerbosePlanarDislocationNode(4,"  bndNeighborStraight "<<std::get<1>(pair.second)->tag()<< " straight?"<<bndNeighborStraight<<std::endl;);
            ////                        //                        isMovable*=bndNeighborStraight;
            ////                        //                        if(!isMovable)
            ////                        //                        {
            ////                        //                            break;
            ////                        //                        }
            ////                    }
            //
            ////                    if(std::get<1>(pair.second)->isGrainBoundarySegment())
            ////                    {// grain-boundary segments must remain grain-boundary if this node is moved
            ////                        for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
            ////                        {
            ////                            const bool gbNeighborMovable=gb->contains(X);
            ////                            VerbosePlanarDislocationNode(4,"  gbNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<gbNeighborMovable<<std::endl;);
            ////                            isMovable*=gbNeighborMovable;
            ////                            if(!isMovable)
            ////                            {
            ////                                break;
            ////                            }
            ////                        }
            ////                        if(!isMovable)
            ////                        {
            ////                            break;
            ////                        }
            ////                    }
            //
            //                    if(std::get<1>(pair.second)->isSessile())
            //                    {// sessile segments cannot change direction if this node is moved
            //                        const double currentNorm((std::get<0>(pair.second)->get_P()-this->get_P()).norm());
            //                        const double newNorm((std::get<0>(pair.second)->get_P()-X).norm());
            //                        if(currentNorm>FLT_EPSILON && newNorm>FLT_EPSILON)
            //                        {
            //                            const bool sessileNeighborMovable=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*currentNorm*newNorm);
            //                            VerbosePlanarDislocationNode(4,"  sessileNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<sessileNeighborMovable<<std::endl;);
            //                            isMovable*=sessileNeighborMovable;
            //                            if(!isMovable)
            //                            {
            //                                break;
            //                            }
            //                        }
            //                    }
            //                }
            //            }
            
            return isMovable;
        }
        
        /**********************************************************************/
        void moveGlide(const double & dt)
        {
            VerbosePlanarDislocationNode(3,"moving PlanarDislocationNode "<<this->sID<<std::endl;);
            const VectorDim P_old(this->get_P());

            if(!masterNode)
            {
                VectorDim dX=velocity.template segment<dim>(0)*dt;
                VerbosePlanarDislocationNode(3,"moving PlanarDislocationNode "<<this->sID<<", dX="<<dX.transpose()<<std::endl;);
                const VectorDim newP(this->snapToGlidePlanes(this->get_P()+dX));
                set_P(newP);
            }
//            else
//            {
//                const VectorDim imageP(imagePosition(masterNode,masterSegmentFaces));
//                static_cast<NodeBaseType*>(this)->set_P(imageP);
//            }
            
            
            //            if (dX.squaredNorm()>0.0 /*&& this->isGlissile()*/) // move a node only if |v|!=0
            //            {
            //                // Make sure that new position is at intersection of glidePlanes
            //                const VectorDim newP(this->snapToGlidePlanes(this->get_P()+dX));
            //                set_P(newP);
            //            }
            //            else
            //            {
            //                resetVirtualBoundaryNode();
            //                velocity.setZero();
            //            }
            
            // Store actual velocity
            if(dt>0.0)
            {
                velocity=(this->get_P()-P_old)/dt;
            }
        }
        
        
        
    };
    
    // static data
    template <typename Derived,typename InterpolationType>
    bool PlanarDislocationNode<Derived,InterpolationType>::use_velocityFilter=true;
    
    template <typename Derived,typename InterpolationType>
    double PlanarDislocationNode<Derived,InterpolationType>::velocityReductionFactor=0.75;
    
    template <typename Derived,typename InterpolationType>
    const double PlanarDislocationNode<Derived,InterpolationType>::bndTol=1.0e-4;
    
    template <typename Derived,typename InterpolationType>
    int PlanarDislocationNode<Derived,InterpolationType>::verbosePlanarDislocationNode=0;
    
}
#endif
