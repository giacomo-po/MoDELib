/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po             <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez       <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Mamdouh Mohamed        <msm07d@fsu.edu>
 * Copyright (C) 2011 by Tamer Crsoby           <tamercrosby@gmail.com>
 * Copyright (C) 2011 by Can Erel               <canerel55@gmail.com>
 * Copyright (C) 2011 by Yinan Cui              <cuiyinan@ucla.edu>
 * Copyright (C) 2017 by Sabyasachi Chatterjee  <sabyasac@andrew.cmu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// valgrind --leak-check=full --show-leak-kinds=all ./DDomp

#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif

// template header cpp
// https://www.codeproject.com/Articles/48575/How-to-Define-a-Template-Class-in-a-h-File-and-Imp

#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <vector>
#include <chrono>
#include <map>
#include <memory>

#include <Eigen/Dense>

#include <TerminalColors.h>

#include <DislocationNetworkTraits.h>
#include <LoopNetwork.h>
//#include <DislocationNetworkComponent.h>
#include <DislocationLoop.h>
#include <DislocationLoopNode.h>
#include <DislocationLoopLink.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
#include <Hermite.h>

#include <DislocationNetworkRemesh.h>
#include <DislocationJunctionFormation.h>
//#include <DislocationCrossSlip.h>
////#include <Material.h>
#include <UniqueOutputFile.h>
#include <DislocationNetworkIO.h>
//#include <DislocationParticle.h>
#include <DislocationStress.h>
////#include <ParticleSystem.h>
//#include <MPIcout.h>
////#include <SingleFieldPoint.h>
//#include <DDtimeIntegrator.h>
//#include <EqualIteratorRange.h>
////#include <BoundingLineSegments.h>
//#include <GrainBoundaryTransmission.h>
////#include <GrainBoundaryDissociation.h>
#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <BVPsolver.h>
#include <ExternalLoadControllerBase.h>
#include <GlidePlaneModule.h>

#include <DislocationNodeContraction.h>
#include <EshelbyInclusion.h>
#include <DDconfigIO.h>
#include <DislocationGlideSolver.h>
//#include <TextFileParser.h>
////#include <DisplacementPoint.h>


////#include <ExternalLoadController.h>
//#include <DislocationInjector.h>
//#include <PeriodicDislocationLoop.h>
////#include <PeriodicLoopObserver.h>
//
////#include <PeriodicDislocationSuperLoop.h>
////#include <PlanarDislocationSuperLoop.h>
//
//#ifdef _MODEL_GREATWHITE_
//#include <MooseSolution.h>
//#endif

#ifndef NDEBUG
#define VerboseDislocationNetwork(N,x) if(verboseDislocationNetwork>=N){std::cout<<x;}
#else
#define VerboseDislocationNetwork(N,x)
#endif

namespace model
{
    
    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationNetwork :public LoopNetwork<DislocationNetwork<dim,corder,InterpolationType> >
    //    /* base                 */ public ParticleSystem<DislocationParticle<_dim> >,
    /*                      */,public std::map<size_t,EshelbyInclusion<dim>>
    {
        
    public:
        
        typedef TypeTraits<DislocationNetwork<dim,corder,InterpolationType>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::VectorLowerDim VectorLowerDim;
        typedef typename TraitsType::MatrixDim MatrixDim;
        typedef DislocationNetworkIO<LoopNetworkType> DislocationNetworkIOType;
        typedef std::map<size_t,EshelbyInclusion<dim>> EshelbyInclusionContainerType;
        typedef BVPsolver<dim,2> BvpSolverType;
        typedef typename BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;

        
        static int verboseDislocationNetwork;


    public:

        const DefectiveCrystalParameters& simulationParameters;
        const SimplicialMesh<dim>& mesh;
        const Polycrystal<dim>& poly;
        GlidePlaneFactory<dim> glidePlaneFactory;
        std::shared_ptr<PeriodicGlidePlaneFactory<dim>> periodicGlidePlaneFactory;
//        const std::unique_ptr<PeriodicDislocationLoopFactory<DislocationNetworkType>> periodicDislocationLoopFactory;
        const std::unique_ptr<BVPsolver<dim,2>>& bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>>& externalLoadController;
        const std::vector<VectorDim>& periodicShifts;
        DislocationNetworkRemesh<LoopNetworkType> networkRemesher;
        DislocationJunctionFormation<LoopNetworkType> junctionsMaker;
        DislocationNodeContraction<LoopNetworkType> nodeContractor;
//        GrainBoundaryTransmission<DislocationNetworkType> gbTransmission;
        //        MatrixDimD _plasticDistortionFromVelocities;
        std::pair<double,MatrixDim> oldPlasticDistortionFromAreas;
        MatrixDim _plasticDistortionRateFromVelocities;
        MatrixDim _plasticDistortionRateFromAreas;
        int ddSolverType;
        bool computeDDinteractions;
        int crossSlipModel;
        int  outputFrequency;
        bool outputBinary;
        bool outputGlidePlanes;
        //        bool outputSpatialCells;
        bool outputElasticEnergy;
        bool outputMeshDisplacement;
        bool outputFEMsolution;
        bool outputDislocationLength;
//        bool outputPlasticDistortion;
        bool outputPlasticDistortionRate;
        bool outputQuadraturePoints;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        const bool computeElasticEnergyPerLength;
//        bool outputPeriodicConfiguration;
        //        unsigned int _userOutputColumn;
        bool use_stochasticForce;
        double surfaceAttractionDistance;
        //        bool computePlasticDistortionRateFromVelocities;
        std::string folderSuffix;
        std::set<const LoopNodeType*> danglingBoundaryLoopNodes;

        /**********************************************************************/
        DislocationNetwork(int& argc, char* argv[],
                           const DefectiveCrystalParameters& _simulationParameters,
                           const SimplicialMesh<dim>& _mesh,
                           const Polycrystal<dim>& _poly,
                           const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                           const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                           const std::vector<VectorDim>& _periodicShifts,
                           long int& runID);
        
        void setConfiguration(const DDconfigIO<dim>&);
        void createEshelbyInclusions();
        const MatrixDim& plasticDistortionRate() const;
        MatrixDim plasticDistortion() const;
        MatrixDim plasticStrainRate() const;
        void updateGeometry();//
        void updatePlasticDistortionRateFromAreas();
        void dummyMove(const int&);
        DislocationNetworkIOType io();
        DislocationNetworkIOType io() const;
        std::tuple<double,double,double,double> networkLength() const;
        const EshelbyInclusionContainerType& eshelbyInclusions() const;
        EshelbyInclusionContainerType& eshelbyInclusions();
        VectorDim displacement(const VectorDim& x) const;
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const;
        MatrixDim stress(const VectorDim& x) const;
        void stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const;
        void assembleAndSolveGlide(const long int& runID);
        void moveGlide(const double & dt_in);
        void singleGlideStepDiscreteEvents(const long int& runID);
        void updateBoundaryNodes();
        bool contract(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);

        
        //
        //    public:
        //
        //
        //        static constexpr int dim=_dim; // make dim available outside class
        //        static constexpr int corder=_corder; // make dim available outside class
        //        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        //        typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
        //        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        //        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        //        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        //        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        //        typedef Eigen::Matrix<double,dim,dim>    MatrixDimD;
        //        typedef Eigen::Matrix<double,dim,1>        VectorDim;
        //        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
        //        //        typedef DislocationParticle<_dim> DislocationParticleType;
        //        //        typedef typename DislocationParticleType::StressField StressField;
        //        //        typedef typename DislocationParticleType::DisplacementField DisplacementField;
        //        //        typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
        //        //        typedef typename ParticleSystemType::SpatialCellType SpatialCellType;
        //        //        typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
        //        typedef BVPsolver<dim,2> BvpSolverType;
        //        typedef typename BvpSolverType::FiniteElementType FiniteElementType;
        //        typedef typename FiniteElementType::ElementType ElementType;
        //        typedef typename LoopNetworkType::IsNodeType IsNodeType;
        //        typedef DislocationNetworkIO<DislocationNetworkType> DislocationNetworkIOType;
        //        typedef Polycrystal<dim> PolycrystalType;
        //        //        typedef ExternalLoadControllerBase<dim> ExternalLoadControllerType;
        //        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        //        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        //        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        //
        //#ifdef DislocationNucleationFile
        //#include DislocationNucleationFile
        //#endif
        //
        //#ifdef _MODEL_GREATWHITE_
        //#include <DislocationNetworkGreatWhite.h>
        //#endif
        //
        //    private:
        //
        //
        //        /**********************************************************************/
        //        void updateVirtualBoundaryLoops()
        //        {
        //
        //
        //            if(   simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM
        //               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
        //            {
        //                //                    if(useVirtualExternalLoops)
        //                //                    {
        //
        //                const auto t0= std::chrono::system_clock::now();
        //                model::cout<<"        Updating virtual boundary loops "<<std::flush;
        //
        //
        //                // First clean up outdated boundary loops
        //                std::set<size_t> removeLoops;
        //                for(const auto& loop : this->loops())
        //                {
        //                    if((loop.second->isVirtualBoundaryLoop() && loop.second->links().size()!=4) || loop.second->isPureVirtualBoundaryLoop())
        //                    {// clean up left over loops from topological operations
        //                        removeLoops.insert(loop.second->sID);
        //                    }
        //                }
        //
        //                for(const size_t& loopID : removeLoops)
        //                {// Remove the virtual loops with ID in removeLoops
        //                    this->deleteLoop(loopID);
        //                }
        //
        //
        //
        //                // Now reconstruct virtual boundary loops
        //                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> virtualLoopVector;
        //                for(const auto& link : this->links())
        //                {
        //                    if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers())
        //                    {
        //                        virtualLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink,link.second->source,link.second->source->virtualBoundaryNode(),link.second->sink->virtualBoundaryNode()},
        //                                                       link.second->burgers(),
        //                                                       (*link.second->grains().begin())->grainID,
        //                                                       DislocationLoopIO<dim>::VIRTUALLOOP);
        //                    }
        //                }
        //
        //                for(const auto& tup : virtualLoopVector)
        //                {// Insert the new virtual loops
        //                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
        //                }
        //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        //            }
        //        }
        //
        //
        //
        
//
//        /**********************************************************************/
//        void setConfiguration(const DDconfigIO<dim>& evl)
//        {
//            this->loopLinks().clear(); // erase base network
//            std::map<size_t,std::shared_ptr<NodeType>> tempNodes;
//            createVertices(evl,tempNodes);
//            createEdges(evl,tempNodes);
//            updatePlasticDistortionRateFromAreas();
//#ifdef _MODEL_MPI_
//            // Avoid that a processor starts writing before other are done reading
//            MPI_Barrier(MPI_COMM_WORLD);
//#endif
//            // Initializing configuration
////            this->io().output(simulationParameters.runID);
//            moveGlide(0.0);    // initial configuration
//        }
//
//        /**********************************************************************/
//        void createVertices(const DDconfigIO<dim>& evl,std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
//        {/*!Creates DislocationNode(s) based on the data read by the DDconfigIO<dim>
//          * object.
//          */
//            size_t kk(1);
//            for (const auto& node : evl.nodes())
//            {
//                const size_t nodeIDinFile(node.sID);
//                NodeType::set_count(nodeIDinFile);
//                if(node.sID==node.masterID)
//                {// a regular node is created
//                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
//
//                    const size_t nodeID(StaticID<NodeType>::nextID());
//                    const auto inserted(tempNodes.emplace(std::piecewise_construct,
//                                                          std::make_tuple(nodeID),
//                                                          std::make_tuple(new NodeType(this,node.P,node.V,node.velocityReduction))));
//                    assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
//                    assert(inserted.first->first == nodeID && "KEY != nodeID");
//                    assert(inserted.first->second->sID == nodeID && "sID != nodeID");
//                    assert(nodeID==nodeIDinFile);
//                }
//                else
//                {
//                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<"), virtual of "<<node.masterID<<std::endl;
//                    const auto isNode(this->node(node.masterID));
//                    assert(isNode.first);
//                    isNode.second->resetVirtualBoundaryNode();
//                }
//                kk++;
//            }
//        }
//
//        std::shared_ptr<NodeType> getSharedNode(const size_t& nodeID,
//                                                const std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
//        {
//            const auto isNode(this->node(nodeID));
//            assert(isNode.first);
//            if(isNode.second->masterNode)
//            {// a virtual node
//                return isNode.second->masterNode->virtualBoundaryNode();
//                //                loopNodes.push_back(isNode.second->masterNode->virtualBoundaryNode());
//            }
//            else
//            {
//                const auto tempNodesFound(tempNodes.find(nodeID));
//                if(tempNodesFound==tempNodes.end())
//                {
//                    model::cout<<"node "<<nodeID<<" not found"<<std::endl;
//                    assert(false && "node shared pointer not found");
//                    return nullptr;
//                }
//                else
//                {
//                    return tempNodesFound->second;
//                }
//                //                loopNodes.push_back(isSharedNode.second);
//            }
//        }
//
//        /**********************************************************************/
//        void createEdges(const DDconfigIO<dim>& evl,const std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
//        {/*!
//          */
//
//
////            std::map<size_t,std::shared_ptr<PeriodicDislocationLoopType> > periodicDislocationLoopMap;
////            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
////               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
////            {
////                for(const auto& pLoop : evl.periodicLoops())
////                {// Collect LoopLinks by loop IDs
////                    StaticID<PeriodicDislocationLoopType>::set_count(pLoop.sID);
////                    periodicDislocationLoopMap.emplace(pLoop.sID,new PeriodicDislocationLoopType(this));
////                }
////            }
//
//            std::map<size_t,std::map<size_t,size_t>> loopMap;
//            for(const auto& looplink : evl.links())
//            {// Collect LoopLinks by loop IDs
//                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
//            }
//
//            assert(loopMap.size()==evl.loops().size());
//
//            size_t loopLumber=1;
//            for(const auto& loop : evl.loops())
//            {// for each loop in the DDconfigIO<dim> object
//
//                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
//                assert(loopFound!=loopMap.end());
//                std::vector<std::shared_ptr<NodeType>> loopNodes;
//                loopNodes.push_back(getSharedNode(loopFound->second.begin()->first,tempNodes));
//                for(size_t k=0;k<loopFound->second.size();++k)
//                {
//                    const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
//                    if(k<loopFound->second.size()-1)
//                    {
//                        loopNodes.push_back(getSharedNode(nodeFound->second,tempNodes));
//                    }
//                    else
//                    {
//                        assert(nodeFound->second==loopNodes[0]->sID);
//                    }
//                }
//
//                LoopType::set_count(loop.sID);
//
//                const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
//
//                model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
//                switch (loop.loopType)
//                {
//                    case DislocationLoopIO<dim>::GLISSILELOOP:
//                    {
////                        LatticePlane loopPlane(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
////                        GlidePlaneKey<dim> loopPlaneKey(loop.grainID,loopPlane);
//                        GlidePlaneKey<dim> loopPlaneKey(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
//                        if(simulationParameters.isPeriodicSimulation())
//                        {
//                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),loop.periodicShift)->sID;
//                            assert(loop.sID==newLoopID);
//                        }
//                        else
//                        {
//                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
//                            assert(loop.sID==newLoopID);
//                        }
////                        if(loop.periodicLoopID>=0)
////                        {// a loop belonging to a periodic loop
////                            const auto pLoopIter(periodicDislocationLoopMap.find(loop.periodicLoopID));
////                            if(pLoopIter!=periodicDislocationLoopMap.end())
////                            {
////                                const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),pLoopIter->second,loop.periodicShift)->sID;
////                                assert(loop.sID==newLoopID);
////                            }
////                            else
////                            {
////                                assert(false && "PeriodicLoop not found in map");
////                            }
////                        }
////                        else
////                        {
////                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
////                            assert(loop.sID==newLoopID);
////                        }
//                        break;
//                    }
//
//                    case DislocationLoopIO<dim>::SESSILELOOP:
//                    {
//                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
//                        assert(loop.sID==newLoopID);
//                        break;
//                    }
//
//                    case DislocationLoopIO<dim>::VIRTUALLOOP:
//                    {
//                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
//                        assert(loop.sID==newLoopID);
//                        break;
//                    }
//
//                    default:
//                        assert(false && "Unknown DislocationLoop type");
//                        break;
//                }
//
//                loopLumber++;
//            }
//
//            model::cout<<std::endl;
//        }
//
        /**********************************************************************/

        
//
//        /**********************************************************************/
//        bool remove(const size_t& nodeID)
//        {
//            const auto isNode=this->node(nodeID);
//            if(isNode.first)
//            {// remove virtual node together with current node
//                if(isNode.second->virtualBoundaryNode())
//                {
//                    LoopNetworkType::remove(isNode.second->virtualBoundaryNode()->sID);
//                }
//            }
//            return LoopNetworkType::remove(nodeID);
//        }
//
        /**********************************************************************/
//        /**********************************************************************/
//        void removeZeroAreaLoops()
//        {
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"        Removing zero-area loops "<<std::flush;
//            std::deque<size_t> loopIDs;
//            for(const auto& loop : this->loops())
//            {
//                if(loop.second->slippedArea()<FLT_EPSILON)
//                {
//                    loopIDs.push_back(loop.second->sID);
//                }
//            }
//
//            for(const auto& loopID : loopIDs)
//            {
//                this->deleteLoop(loopID);
//            }
//            model::cout<<"("<<loopIDs.size()<<" removed)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//        }
//
//        /**********************************************************************/
//        void singleGlideStepDiscreteEvents(const long int& runID)
//        {
//            //! A simulation step consists of the following:
//            //            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
//            //            /*                    */<< ", time="<<totalTime
//            //            /*                    */<< ": nodes="<<this->nodes().size()
//            //            /*                    */<< ", segments="<<this->links().size()
//            //            /*                    */<< ", loopSegments="<<this->loopLinks().size()
//            //            /*                    */<< ", loops="<<this->loops().size()
//            //            /*                    */<< ", components="<<this->components().size()
//            //            /*                    */<< defaultColor<<std::endl;
//
//            //! 1- Check that all nodes are balanced
//            //            checkBalance();
//
//            //! 2 - Update quadrature points
//            //            updateQuadraturePoints();
//            //            updateStressStraightSegments();
//
//            //            for(auto& loop : this->loops()) // TODO: PARALLELIZE THIS LOOP
//            //            {// copmute slipped areas and right-handed normal
//            //                loop.second->updateGeometry();
//            //            }
//            //            updatePlasticDistortionRateFromAreas();
//
//            //! 3- Calculate BVP correction
//            //            updateLoadControllers(runID);
//
//            //#ifdef DislocationNucleationFile
//            //            if(use_bvp && !(runID%use_bvp))
//            //            {
//            //                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
//            //                updateQuadraturePoints();
//            //            }
//            //#endif
//            //            assembleAndSolve(runID,straightSegmentsDeq);
//            //            computeNodaVelocities(runID);
//
//
//            //! 4- Solve the equation of motion
//
//
//            //! 5- Compute time step dt (based on max nodal velocity) and increment totalTime
//            // make_dt();
//
//
//
//            //! 6- Output the current configuration before changing it
//            //            output(runID);
//            //            io().output(runID);
//
//
//            //! 7- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
//            //            move(dt);
//
//            //! 8- Update accumulated quantities (totalTime and plasticDistortion)
//            //            totalTime+=dt;
//            //            updatePlasticDistortionRateFromVelocities();
//
//
//            //! 9- Contract segments of zero-length
//            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
//
//            //            if(runID>0)
//            //            {
//            //                removeZeroAreaLoops();
//            //            }
//
//            //! 10- Cross Slip (needs upated PK force)
//            DislocationCrossSlip<DislocationNetworkType> crossSlip(*this);
//            crossSlip.execute();
//
//
//            //                        gbTransmission.transmit();
//            //gbTransmission.directTransmit();
//
//            //
//            //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
//
//            //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
//
//
//            //! 11- detect loops that shrink to zero and expand as inverted loops
//            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
//
//            //! 12- Form Junctions
//            junctionsMaker.formJunctions(DDtimeIntegrator<0>::dxMax);
//
//            //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
//            //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
//
//
//
//            //! 13- Node redistribution
//            networkRemesher.remesh(runID);
//
//            //            createImageLoops();
//            //            updateImageLoops();
//            updateVirtualBoundaryLoops();
//
//            //            mergeLoopsAtNodes();
//
//            //            DislocationInjector<DislocationNetworkType>(*this).insertRandomStraightDislocation();
//
//            //! 9- Contract segments of zero-length
//            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
//
//            //! 14- If BVP solver is not used, remove DislocationSegment(s) that exited the boundary
//            //            removeBoundarySegments();
//
//            //            removeSmallComponents(3.0*dx,4);
//
//            //            make_bndNormals();
//
//            //! 16 - Increment runID counter
//            //            ++runID;     // increment the runID counter
//        }
//
//        //        const VectorDim periodicVector(const Eigen::Array<int,dim,1>& cellID) const
//        //        {
//        //            return (meshDimensions.array()*cellID.template cast<double>()).matrix();
//        //        }
//
        /**********************************************************************/

//
//
//        /**********************************************************************/
//        bool contract(std::shared_ptr<NodeType> nA,
//                      std::shared_ptr<NodeType> nB)
//        {
//            return nodeContractor.contract(nA,nB);
//        }
//
        /**********************************************************************/

//
//
//
//
//        /**********************************************************************/
//        void moveGlide(const double & dt_in)
//        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
//          */
//
//
//            if(simulationParameters.isPeriodicSimulation())
//            {
//
//
//
//                std::map<Eigen::Matrix<double, DislocationNetworkType::dim,1>, const std::shared_ptr<NodeType>, CompareVectorsByComponent<double,DislocationNetworkType::dim,float>> rveNodesMap;
//                for(const auto& pair : *periodicDislocationLoopFactory)
//                {// output periodic glide planes too
//
//                    if(!pair.second.expired())
//                    {
//                        const auto periodicLoop(pair.second.lock());
////                        periodicLoop->updateRVEloops(*this,dt_in,rveNodesMap);
//                    }
//                }
//
//
//            }
//            else
//            {
//                model::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
//                const auto t0= std::chrono::system_clock::now();
//                for (auto& nodeIter : this->nodes())
//                {
//                    nodeIter.second->moveGlide(dt_in);
//                }
//                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//            }
//
//        }
//
//
//
//        /**********************************************************************/
//        void updatePlasticDistortionRateFromAreas()
//        {
//            //            if(!computePlasticDistortionRateFromVelocities)
//            //            {
//
//            const double dt(simulationParameters.totalTime-_plasticDistortionFromAreas.first);
//            if(dt>=0.0)
//            {
//                const MatrixDimD pd(plasticDistortion());
//                if(dt>0.0)
//                {
//                    _plasticDistortionRateFromAreas=(pd-_plasticDistortionFromAreas.second)/dt;
//                }
//                _plasticDistortionFromAreas=std::make_pair(simulationParameters.totalTime,pd); // update old values
//            }
////
////            const MatrixDimD old(_plasticDistortionFromAreas);
//////            _plasticDistortionFromAreas.setZero();
//////            for(const auto& loop : this->loops())
//////            {
//////                _plasticDistortionFromAreas+= loop.second->plasticDistortion();
//////            }
////            _plasticDistortionFromAreas=plasticDistortion();
////            if(dt>0.0)
////            {
////                _plasticDistortionRateFromAreas=(_plasticDistortionFromAreas-old)/dt;
////            }
//            //            }
//        }
//

//
//        /**********************************************************************/
//        std::tuple<double,double,double,double> networkLength() const
//        {/*!\returns the total line length of the DislocationNetwork. The return
//          * value is a tuple, where the first value is the length of bulk glissile
//          * dislocations, the second value is the length of bulk sessile
//          * dislocations, and the third value is the length accumulated on
//          * the mesh boundary.
//          */
//            double bulkGlissileLength(0.0);
//            double bulkSessileLength(0.0);
//            double boundaryLength(0.0);
//            double grainBoundaryLength(0.0);
//
//            for(auto& loop : this->loops())
//            {
//                for(const auto& loopLink : loop.second->links())
//                {
//                    if(!loopLink.second->pLink->hasZeroBurgers())
//                    {
//                        if(loopLink.second->pLink->isBoundarySegment())
//                        {
//                            boundaryLength+=loopLink.second->pLink->chord().norm();
//                        }
//                        else if(loopLink.second->pLink->isGrainBoundarySegment())
//                        {
//                            grainBoundaryLength+=loopLink.second->pLink->chord().norm();
//                        }
//                        else
//                        {
//                            if(loopLink.second->pLink->isSessile())
//                            {
//                                bulkSessileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
//                            }
//                            else
//                            {
//                                bulkGlissileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
//                            }
//                        }
//                    }
//                }
//            }
//            return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
//        }
//
//
//
        /**********************************************************************/

        
    };
    
}
#endif
