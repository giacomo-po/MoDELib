/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_PlanarDislocationSegment_H
#define model_PlanarDislocationSegment_H

#include <memory>
#include <set>
#include <algorithm>    // std::set_intersection, std::sort


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <Quadrature.h>
#include <QuadPow.h>
#include <DislocationNetworkTraits.h>
#include <SplineSegment.h>
#include <Material.h>
#include <GlidePlaneObserver.h>
#include <GlidePlane.h>
#include <Coeff2Hermite.h>
#include <DislocationParticle.h>
#include <UniqueOutputFile.h>
#include <GrainBoundary.h>
#include <LineSimplexIntersection.h>
#include <BoundingLineSegments.h>
#include <BoundingLineSegments.h>
#include <MeshPlane.h>
#include <DislocationQuadraturePoint.h>
#include <StraightDislocationSegment.h>


#ifndef NDEBUG
#define VerbosePlanarDislocationSegment(N,x) if(verbosePlanarDislocationSegment>=N){model::cout<<x;}
#else
#define VerbosePlanarDislocationSegment(N,x)
#endif


namespace model
{
    
    
//    template<typename LinkType>
//    class VirtualBoundaryLoopsManager
//    {
//        
//        typedef typename LinkType::LoopType LoopType;
//        typedef typename LinkType::NodeType NodeType;
//        typedef LoopLink<LinkType> LoopLinkType;
//        
//        LinkType& link;
//        //        const std::shared_ptr<LoopType> innermostLoop;
//        //        const std::shared_ptr<LoopType> outermostLoop;
//        
//        //        const std::unique_ptr<LoopType> loop;
//        
//        //        const std::unique_ptr<LoopLinkType> sinkTosource;
//        //        const std::unique_ptr<LoopLinkType> sourceToVirtualSource;
//        //        const std::unique_ptr<LoopLinkType> virtualSourceToSink;
//        //        const std::unique_ptr<LoopLinkType> sinkToVirtualSource;
//        //        const std::unique_ptr<LoopLinkType> virtualSourceToVirtualSink;
//        //        const std::unique_ptr<LoopLinkType> virtualSinkToSink;
//        
//        
//        //        const std::unique_ptr<LoopLink<LinkType>>
//        
//        static std::pair<size_t,size_t> getLoopIDs(LinkType& link)
//        {
//            assert(link.grains().size()==1);
//            assert(link.burgers().norm()>FLT_EPSILON);
//            const Grain<LinkType::dim>& grain(**link.grains().begin());
//            
////            std::shared_ptr<NodeType> virtualSource(link.source->virtualBoundaryNode());
////            std::shared_ptr<NodeType> virtualSink(link.sink->virtualBoundaryNode());
//            
//            std::vector<std::shared_ptr<NodeType>> innerNodes{link.sink,link.source,link.source->virtualBoundaryNode()};
//            std::shared_ptr<LoopType> innerLoop=link.network().insertLoop(innerNodes,link.burgers(),grain.grainID);
//            
//            std::vector<std::shared_ptr<NodeType>> outerNodes{link.sink,link.source->virtualBoundaryNode(),link.sink->virtualBoundaryNode()};
//            std::shared_ptr<LoopType> outerLoop=link.network().insertLoop(outerNodes,innerLoop->burgers(),grain.grainID);
//            
//            innerLoop->update();
//            outerLoop->update();
//            
//            std::cout<<"innerLoop "<<innerLoop.use_count()<<std::endl;
//            std::cout<<"outerLoop "<<outerLoop.use_count()<<std::endl;
//            
//            assert(innerLoop->burgers().norm()>FLT_EPSILON);
//            assert(outerLoop->burgers().norm()>FLT_EPSILON);
//
//            return std::make_pair(innerLoop->sID,outerLoop->sID);
//        }
//        
//        const std::pair<size_t,size_t> loopIDs;
//        
//    public:
//        
//        
//        
//        VirtualBoundaryLoopsManager(LinkType& _link) :
//        /* init */ link(_link)
//        /* init */,loopIDs(getLoopIDs(link))
//        //        /* init */ innermostLoop(new LoopType(&link.network(),link.burgers(),(*link.grains().begin())->grainID))
//        //        /* init */,outermostLoop(new LoopType(&link.network(),link.burgers(),(*link.grains().begin())->grainID))
//        //        /* init */,sinkTosource(new LoopLinkType(link.sink,link.source,innermostLoop))
//        //        /* init */,sourceToVirtualSource(new LoopLinkType(link.source,link.source->virtualBoundaryNode(),innermostLoop))
//        //        /* init */,virtualSourceToSink(new LoopLinkType(link.source->virtualBoundaryNode(),link.sink,innermostLoop))
//        //        /* init */,sinkToVirtualSource(new LoopLinkType(link.sink,link.source->virtualBoundaryNode(),outermostLoop))
//        //        /* init */,virtualSourceToVirtualSink(new LoopLinkType(link.source->virtualBoundaryNode(),link.sink->virtualBoundaryNode(),outermostLoop))
//        //        /* init */,virtualSinkToSink(new LoopLinkType(link.sink->virtualBoundaryNode(),link.sink,outermostLoop))
//        //        /* init */ sinkTosource(new LoopLinkType(link.sink,link.source,std::shared_ptr<LoopType>(new LoopType(&link.network(),link.burgers(),(*link.grains().begin())->grainID))))
//        //        /* init */,sourceToVirtualSource(new LoopLinkType(link.source,link.source->virtualBoundaryNode(),sinkTosource->loop()))
//        //        /* init */,virtualSourceToSink(new LoopLinkType(link.source->virtualBoundaryNode(),link.sink,sinkTosource->loop()))
//        //        /* init */,sinkToVirtualSource(new LoopLinkType(link.sink,link.source->virtualBoundaryNode(),std::shared_ptr<LoopType>(new LoopType(&link.network(),link.burgers(),(*link.grains().begin())->grainID))))
//        //        /* init */,virtualSourceToVirtualSink(new LoopLinkType(link.source->virtualBoundaryNode(),link.sink->virtualBoundaryNode(),sinkToVirtualSource->loop()))
//        //        /* init */,virtualSinkToSink(new LoopLinkType(link.sink->virtualBoundaryNode(),link.sink,sinkToVirtualSource->loop()))
//        {
//            
//            std::cout<<"Created VirtualBoundaryLoops "<<loopIDs.first<<","<<loopIDs.second<<" for link "<<link.tag()<<std::endl;
//            
//            
//        }
//        
////        void update()
////        {
////            link.network().deleteLoop(loopIDs.first);
////            link.network().deleteLoop(loopIDs.second);
////            loopIDs(getLoopIDs(link));
////        }
//        
//        /******************************************************************/
//        ~VirtualBoundaryLoopsManager()
//        {
//            std::cout<<"Destroying VirtualBoundaryLoopsManager "<<link.tag()<<std::endl;
//            
//            link.network().deleteLoop(loopIDs.first);
//            link.network().deleteLoop(loopIDs.second);
////            assert(link.network().loops().find(loopIDs.first )==link.network().loops().end());
////            assert(link.network().loops().find(loopIDs.second)==link.network().loops().end());
//        }
//        
//    };
    
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived>
    class PlanarDislocationSegment : public SplineSegment<Derived,TypeTraits<Derived>::dim,TypeTraits<Derived>::corder>
    /*                      */,private std::set<const MeshPlane<TypeTraits<Derived>::dim>*>
    /*                      */,private std::set<const GrainBoundary<TypeTraits<Derived>::dim>*>
    /*                      */,private std::set<const Grain<TypeTraits<Derived>::dim>*>
    /*                      */,private BoundingLineSegments<TypeTraits<Derived>::dim>
    /*                      */,public DislocationQuadraturePointContainer<TypeTraits<Derived>::dim,TypeTraits<Derived>::corder>
    {
        
        
    public:
        
        static constexpr int dim=TypeTraits<Derived>::dim; // make dim available outside class
        static constexpr int corder=TypeTraits<Derived>::corder; // make dim available outside class
        typedef SplineSegmentBase<dim,corder> SplineSegmentBaseType;
        typedef Derived LinkType;
        typedef StraightDislocationSegment<dim> StraightDislocationSegmentType;
        typedef typename TypeTraits<LinkType>::LoopType LoopType;
        typedef typename TypeTraits<LinkType>::LoopNetworkType NetworkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef typename TypeTraits<LinkType>::NodeType NodeType;
        typedef SplineSegment<LinkType,dim,corder> SplineSegmentType;
        static constexpr int Ncoeff=SplineSegmentBaseType::Ncoeff;
        static constexpr int Ndof=SplineSegmentType::Ndof;
        typedef typename SplineSegmentType::MatrixNdof MatrixNdof;
        typedef typename SplineSegmentType::VectorNdof VectorNdof;
        typedef typename SplineSegmentType::VectorDim VectorDim;
        typedef typename SplineSegmentType::MatrixDim MatrixDim;
        typedef typename SplineSegmentType::MatrixDimNdof MatrixDimNdof;
        typedef typename SplineSegmentType::MatrixNcoeff MatrixNcoeff;
        typedef typename SplineSegmentType::MatrixNcoeffDim MatrixNcoeffDim;
        static constexpr int pOrder=SplineSegmentType::pOrder;
        typedef Eigen::Matrix<double,Ncoeff,1>     VectorNcoeff;
        typedef DislocationQuadraturePoint<dim,corder> DislocationQuadraturePointType;
//        typedef typename DislocationQuadraturePointType::QuadPowDynamicType QuadPowDynamicType;
//        typedef typename DislocationQuadraturePointType::QuadratureDynamicType QuadratureDynamicType;
        typedef DislocationParticle<dim> DislocationParticleType;
//        typedef std::vector<DislocationParticleType*> QuadratureParticleContainerType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef std::set<const GrainBoundary<dim>*> GrainBoundaryContainerType;
        typedef std::set<const Grain<dim>*> GrainContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef std::set<const MeshPlaneType*> MeshPlaneContainerType;
        typedef typename TypeTraits<LinkType>::MeshLocation MeshLocation;
        
        
    private:
        
        std::map<size_t,
        /*    */ std::pair<VectorNcoeff,VectorDim>,
        /*    */ std::less<size_t>,
        /*    */ Eigen::aligned_allocator<std::pair<size_t, std::pair<VectorNcoeff,VectorDim>> >
        /*    */ > h2posMap;
        Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg;
        MatrixNdof Kqq; //! Segment Stiffness Matrix
        VectorNdof Fq; //! Segment Nodal Force Vector
        VectorDim Burgers; //! The Burgers vector
        double BurgersNorm;
        bool _isBoundarySegment;
//        std::unique_ptr<VirtualBoundaryLoopsManager<LinkType>> virtualLinks;
        
    public:
        
        
        static const Eigen::Matrix<double,TypeTraits<Derived>::dim,TypeTraits<Derived>::dim> I;
        static const Eigen::Matrix<double,TypeTraits<Derived>::dim,1> zeroVector;
        static double quadPerLength;
        static double virtualSegmentDistance;
        static int verbosePlanarDislocationSegment;
        
//        QuadratureParticleContainerType quadratureParticleContainer;
        StraightDislocationSegment<dim> straight;
        
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            
            LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
            
            quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",true);
            assert((LinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
            
            verbosePlanarDislocationSegment=TextFileParser("inputFiles/DD.txt").readScalar<int>("verbosePlanarDislocationSegment",true);
            virtualSegmentDistance=TextFileParser(fileName).readScalar<double>("virtualSegmentDistance",true);
            
        }
        
        /******************************************************************/
        PlanarDislocationSegment(const std::shared_ptr<NodeType>& nI,
                                 const std::shared_ptr<NodeType>& nJ) :
        /* init */ SplineSegmentType(nI,nJ)
        /* init */,Burgers(VectorDim::Zero())
        /* init */,BurgersNorm(Burgers.norm())
        /* init */,_isBoundarySegment(this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first)
        /* init */,straight(this->source->get_P(),this->sink->get_P(),Burgers,this->chordLength(),this->unitDirection())
        //        /* init */ qOrder(QuadPowDynamicType::lowerOrder(quadPerLength*this->chord().norm()))
        {/*! Constructor with pointers to source and sink, and flow
          *  @param[in] NodePair_in the pair of source and sink pointers
          *  @param[in] Flow_in the input flow
          */
            
            VerbosePlanarDislocationSegment(1,"Creating PlanarDislocationSegment "<<this->tag()<<std::endl;);
            VerbosePlanarDislocationSegment(2,"_isBoundarySegment "<<_isBoundarySegment<<std::endl;);
            VerbosePlanarDislocationSegment(3,"source->isBoundaryNode() "<<this->source->isBoundaryNode()<<std::endl;);
            VerbosePlanarDislocationSegment(3,"sink->isBoundaryNode() "<<this->sink->isBoundaryNode()<<std::endl;);
            VerbosePlanarDislocationSegment(3,"midpoint is boundary "<<boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first<<std::endl;);
            
            
            
            //            pkGauss.setZero(dim,this->quadraturePoints().size()); // necessary if this is not assembled
            //            _isBoundarySegment=this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first;
            
        }
        
        /******************************************************************/
        ~PlanarDislocationSegment()
        {
            VerbosePlanarDislocationSegment(1,"Destroying PlanarDislocationSegment "<<this->tag()<<std::endl;);
        }
        
        
        /**********************************************************************/
        void updateGeometry()
        {
            SplineSegmentType::updateGeometry();
            straight.updateGeometry();
            
            _isBoundarySegment=this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first;
            
//            if(this->network().useVirtualExternalLoops && _isBoundarySegment && !virtualLinks && !hasZeroBurgers())
//            {
//                VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", creating VirtualBoundaryLoopsManager in updateGeometry"<<std::endl;);
//                virtualLinks.reset(new VirtualBoundaryLoopsManager<LinkType>(this->derived()));
//            }
        }
        
        /**********************************************************************/
        const MeshPlaneContainerType& meshPlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        MeshPlaneContainerType& meshPlanes()
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
        GrainBoundaryContainerType& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundaryContainerType& grainBoundaries() const
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
        bool addMeshPlane(const MeshPlaneType& gp)
        {
            const bool success=meshPlanes().insert(&gp).second;
            if(success)
            {
                const bool sourceContained(gp.contains(this->source->get_P()));
                const bool   sinkContained(gp.contains(this->  sink->get_P()));
                if(!(sourceContained && sinkContained))
                {
                    model::cout<<"PlanarDislocationSegment "<<this->source->sID<<"->"<<this->sink->sID<<std::endl;
                    model::cout<<"sourceContained "<<sourceContained<<std::endl;
                    model::cout<<"  sinkContained "<<sinkContained<<std::endl;
                    assert(false &&  "Glide Plane does not contain source or sink");
                }
                boundingBoxSegments().updateWithMeshPlane(gp); // Update _boundingBoxSegments. This must be called before updateGlidePlaneIntersections
                grains().insert(&this->network().poly.grain(gp.regionIDs.first));    // Insert new grain in grainSet
                grains().insert(&this->network().poly.grain(gp.regionIDs.second));   // Insert new grain in grainSet
            }
            return success;
        }
        
        /**********************************************************************/
        void addLink(LoopLinkType* const pL)
        {
            VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", adding LoopLink "<<pL->tag()<<std::endl;);
            
            SplineSegmentType::addLink(pL); // forward to base class
            
            // Modify Burgers vector
            if(pL->source()->sID==this->source->sID)
            {
                Burgers+=pL->flow().cartesian();
            }
            else
            {
                Burgers-=pL->flow().cartesian();
            }
            BurgersNorm=Burgers.norm();
            
            
            if(!pL->loop()->isVirtualBoundaryLoop())
            {
                
                if(pL->loop()->glidePlane)
                {
                    addMeshPlane(*pL->loop()->glidePlane.get());
                }
                _isBoundarySegment=this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first;
                
                
//                if(this->network().useVirtualExternalLoops && _isBoundarySegment && !hasZeroBurgers())
//                {// If a non-virtual loop is added to the segment, the virtual boundary structure must be reset
//                    
//                    VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", creating VirtualBoundaryLoopsManager in addLink"<<std::endl;);
//                    //virtualLinks->update();
//                    virtualLinks.reset(new VirtualBoundaryLoopsManager<LinkType>(this->derived()));
//                }
            }
        }
        
        /**********************************************************************/
        void removeLink(LoopLinkType* const pL)
        {
            VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", removing LoopLink "<<pL->tag()<<std::endl;);
            
            
            SplineSegmentType::removeLink(pL);  // forward to base class
            
            // Modify Burgers vector
            if(pL->source()->sID==this->source->sID)
            {
                Burgers-=pL->flow().cartesian();
            }
            else
            {
                Burgers+=pL->flow().cartesian();
            }
            BurgersNorm=Burgers.norm();
            
            
            if(!pL->loop()->isVirtualBoundaryLoop())
            {
                
                meshPlanes().clear();
                boundingBoxSegments().clear();
                for(const auto& loopLink : this->loopLinks())
                {
                    if(loopLink->loop()->glidePlane)
                    {
                        addMeshPlane(*loopLink->loop()->glidePlane.get());
                    }
                }
                
                addGrainBoundaryPlanes();
                _isBoundarySegment=this->source->isBoundaryNode() && this->sink->isBoundaryNode() && boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first;
                
                
//                virtualLinks.reset(nullptr); // If a non-virtual loop is added to the segment, the virtual boundary structure must be reset
//                
//                
//                if(this->network().useVirtualExternalLoops && _isBoundarySegment && !hasZeroBurgers())
//                {
//                    VerbosePlanarDislocationSegment(2,"PlanarDislocationSegment "<<this->tag()<<", creating VirtualBoundaryLoopsManager in removeLink"<<std::endl;);
//                    virtualLinks.reset(new VirtualBoundaryLoopsManager<LinkType>(this->derived()));
//                }
                
            }
            
        }
        
        /**********************************************************************/
        size_t addGrainBoundaryPlanes()
        {
            size_t addedGp=0;
            for(const auto& gb : this->source->grainBoundaries())
            {
                if(this->sink->grainBoundaries().find(gb)!=this->sink->grainBoundaries().end())
                {
                    grainBoundaries().insert(gb);
                    //                    const auto& gp(gb->glidePlanes().begin()->second);// HERE BEGIN IS TEMPORARY, UNTIL WE STORE THE GLIDE PLANE OF THE CSL AND DSCL
                    //                    addedGp+=addMeshPlane(*gp.get());
                    addedGp+=addMeshPlane(*gb);
                }
            }
            
            return addedGp;
        }
        
//        /**********************************************************************/
//        void addToStressStraight(std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>>& straightSegmentsDeq) const
//        {
//            if(!hasZeroBurgers())
//            {
//                
//                
//                switch (this->network().simulationParameters.simulationType)
//                {// Initilization based on type of simulation
//                        
//                        
//                    case DefectiveCrystalParameters::FINITE_NO_FEM:
//                    {
//                        if(!isBoundarySegment())
//                        {
//                            straightSegmentsDeq.emplace_back(this->source->get_P(),
//                                                             this->sink->get_P(),
//                                                             burgers());
//                        }
//                        break;
//                    }
//                        
//                    case DefectiveCrystalParameters::FINITE_FEM:
//                    {
//                        assert(0 && "FINISH HERE");
//                        break;
//                    }
//                        
//                    case DefectiveCrystalParameters::PERIODIC:
//                    {
//                        if(!isBoundarySegment())
//                        {
//                            straightSegmentsDeq.emplace_back(this->source->get_P(),
//                                                             this->sink->get_P(),
//                                                             burgers());
//                        }
//                        break;
//                    }
//                        
//                    default:
//                    {
//                        model::cout<<"simulationType MUST BE 0,1, or 2. EXITING."<<std::endl;
//                        exit(EXIT_FAILURE);
//                        break;
//                    }
//                }
//                
//            }
//        }
        

        
        /**********************************************************************/
        void assemble()
        {
            this->updateForcesAndVelocities(*this,quadPerLength);
//            quadratureParticleContainer.clear(); // OLD QUADRATURE PARTICLE METHOD
//            quadratureParticleContainer.reserve(this->quadraturePoints().size());  // OLD QUADRATURE PARTICLE METHOD
            
            Fq= this->quadraturePoints().size()? this->nodalVelocityVector() : VectorNdof::Zero();
            Kqq=this->nodalVelocityMatrix(*this);
            h2posMap=this->hermite2posMap();
            
            Mseg.setZero(Ncoeff*dim,h2posMap.size()*dim);
            size_t c=0;
            for(const auto& pair : h2posMap)
            {
                for(int r=0;r<Ncoeff;++r)
                {
                    Mseg.template block<dim,dim>(r*dim,c*dim)=pair.second.first(r)*MatrixDim::Identity();
                }
                c++;
            }
        }
        
        

        
        /**********************************************************************/
        void addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,
                                 Eigen::VectorXd& FQ) const
        {/*!\param[in] kqqT the stiffness matrix of the network component
          * \param[in] FQ the force vector of the network component
          */
            
            if(!hasZeroBurgers())
            {
                const Eigen::MatrixXd tempKqq(Mseg.transpose()*Kqq*Mseg); // Create the temporaty stiffness matrix and push into triplets
                
                
                size_t localI=0;
                for(const auto& pairI : h2posMap)
                {
                    for(int dI=0;dI<dim;++dI)
                    {
                        const size_t globalI=pairI.first*dim+dI;
                        
                        size_t localJ=0;
                        for(const auto& pairJ : h2posMap)
                        {
                            for(int dJ=0;dJ<dim;++dJ)
                            {
                                const size_t globalJ=pairJ.first*dim+dJ;
                                
                                if (std::fabs(tempKqq(localI,localJ))>FLT_EPSILON)
                                {
                                    //                                    kqqT.push_back(Eigen::Triplet<double>(globalI,globalJ,tempKqq(localI,localJ)));
                                    kqqT.emplace_back(globalI,globalJ,tempKqq(localI,localJ));
                                }
                                
                                localJ++;
                            }
                        }
                        
                        localI++;
                    }
                }
                
                const Eigen::VectorXd tempFq(Mseg.transpose()*Fq); // Create temporary force vector and add to global FQ
                
                localI=0;
                for(const auto& pairI : h2posMap)
                {
                    for(int dI=0;dI<dim;++dI)
                    {
                        const size_t globalI=pairI.first*dim+dI;
                        
                        FQ(globalI)+=tempFq(localI);
                        
                        localI++;
                    }
                }
                
            }
        }
        
        /**********************************************************************/
        const MatrixNdof& get_Kqq() const
        {/*!\returns the stiffness matrix of this segment
          */
            return Kqq;
        }
        
        /**********************************************************************/
        const VectorNdof& get_Fq() const
        {/*!\returns the nodal force vector for this segment
          */
            return Fq;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortionRate() const
        {/*!\returns the plastic strain rate generated by *this segment:
          *	\f[
          *		\beta^P_{ij}=\int_{\mathbb{R}^3}\int d\ell dV=-b_i\int\epsilon_{jkl}w_kd\ell_l
          *	\f]
          */
            
            //\todo this integral should be calculated using shape functions
            
            const VectorDim V((this->source->get_V().template segment<dim>(0)+this->sink->get_V().template segment<dim>(0))*0.5);
            return -Burgers*V.cross(this->chord()).transpose()*(!isBoundarySegment())/this->network().mesh.volume();
        }
        
        /**********************************************************************/
        MatrixDim plasticStrainRate() const
        {
            const VectorDim temp(plasticDistortionRate());
            return (temp+temp.transpose())*0.5;
        }
        
        /**********************************************************************/
        bool isGrainBoundarySegment() const
        {
            return grainBoundaries().size();
        }
        

        
        /**********************************************************************/
        const VectorDim& burgers() const
        {
            return Burgers;
        }
        
        /**********************************************************************/
        const VectorDim& glidePlaneNormal() const
        {
            return meshPlanes().size()==1? (*meshPlanes().begin())->unitNormal : zeroVector;
        }
        
        /**********************************************************************/
        bool isVirtualBoundarySegment() const
        {//!\returns true if all loops of this segment are virtualBoundaryLoops
            bool temp(true);
            for(const auto& loopLink : this->loopLinks())
            {
                temp*=loopLink->loop()->isVirtualBoundaryLoop();
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        std::set<const LoopType*> virtualLoops() const
        {//!\returns a set of pointers to the virtualBoundaryLoops of this segment
            std::set<const LoopType*> temp;
            for(const auto& loopLink : this->loopLinks())
            {
                if(loopLink->loop()->isVirtualBoundaryLoop())
                {
                    temp.insert(loopLink->loop().get());
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSessile() const
        {
            return    !isGlissile()
            /*  */ && !isBoundarySegment()
            /*  */ && !isGrainBoundarySegment()
            /*  */ && !hasZeroBurgers()
            /*  */ && !isVirtualBoundarySegment();
        }
        
        /**********************************************************************/
        bool isGlissile() const
        {/*\returns true if ALL the following conditions are met
          * - the segment is confined by only one plane
          * - its Burgers vector is non-zero
          * - all loops containing this segment are glissile
          */
            bool temp(meshPlanes().size()==1 && !hasZeroBurgers() && !isVirtualBoundarySegment());
            if(temp)
            {
                for(const auto& loopLink : this->loopLinks())
                {
                    temp*=loopLink->loop()->isGlissile;
                }
            }
            return temp;
        }
        
        //        /**********************************************************************/
        //        bool hasZeroBurgers() const
        //        {
        //            return Burgers.squaredNorm()<FLT_EPSILON;
        //        }
        
        /**********************************************************************/
        bool hasZeroBurgers() const
        {
            return BurgersNorm<FLT_EPSILON;
        }
        
        //        /**********************************************************************/
        //        double arcLength() const
        //        {
        //            return SplineSegmentType::template arcLength<16,UniformOpen>();
        //        }
        
        //        /**********************************************************************/
        //        VectorDim velocity(const double& u) const
        //        {
        //            return this->source->get_V().template segment<dim>(0)*(1.0-u)+this->sink->get_V().template segment<dim>(0)*u;
        //        }
        

        
        /**********************************************************************/
        const bool& isBoundarySegment() const // THIS IS CALLED MANY TIMES< CONSIDER STORING
        {/*!\returns true if both nodes are boundary nodes, and the midpoint is
          * on the boundary.
          */
            return _isBoundarySegment;
            //            HERE RETURN REFERENCE
            //
            //            return this->source->isBoundaryNode() &&
            //            /*  */ this->sink->isBoundaryNode() &&
            //            /*  */ boundingBoxSegments().contains(0.5*(this->source->get_P()+this->sink->get_P())).first;
        }
        
        /**********************************************************************/
        MeshLocation meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            
            MeshLocation temp = MeshLocation::outsideMesh;
            
            
            if(isBoundarySegment())
            {
                temp=MeshLocation::onMeshBoundary;
            }
            else
            {
                if(isGrainBoundarySegment())
                {
                    temp=MeshLocation::onRegionBoundary;
                }
                else
                {
                    temp=MeshLocation::insideMesh;
                }
            }
            
            return temp;
        }
        
//        /**********************************************************************/
//        template <class T>
//        friend T& operator << (T& os, const LinkType& ds)
//        {
//            os  << ds.source->sID<<"\t"<< ds.sink->sID<<"\t"
//            /**/<< std::setprecision(15)<<std::scientific<<ds.Burgers.transpose()<<"\t"
//            /**/<< std::setprecision(15)<<std::scientific<<ds.glidePlaneNormal().transpose()<<"\t"
//            /**/<<SplineSegmentBase<dim,corder>::sourceT(ds).transpose()<<"\t"
//            /**/<<SplineSegmentBase<dim,corder>::sinkT(ds).transpose()<<"\t"
//            /**/<<ds.meshLocation();
//            return os;
//        }
        
    };
    
    // Static Data
    template <typename Derived>
    const Eigen::Matrix<double,TypeTraits<Derived>::dim,TypeTraits<Derived>::dim> PlanarDislocationSegment<Derived>::I=Eigen::Matrix<double,TypeTraits<Derived>::dim,TypeTraits<Derived>::dim>::Identity();
    
    template <typename Derived>
    const Eigen::Matrix<double,TypeTraits<Derived>::dim,1> PlanarDislocationSegment<Derived>::zeroVector=Eigen::Matrix<double,TypeTraits<Derived>::dim,1>::Zero();
    
    
    template <typename Derived>
    double PlanarDislocationSegment<Derived>::quadPerLength=0.2;
    
    template <typename Derived>
    double PlanarDislocationSegment<Derived>::virtualSegmentDistance=200.0;
    
    template <typename Derived>
    int PlanarDislocationSegment<Derived>::verbosePlanarDislocationSegment=0;
    
    
}
#endif

