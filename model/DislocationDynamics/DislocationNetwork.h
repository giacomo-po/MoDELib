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


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <vector>
#include <chrono>
#include <map>
#include <memory>

#include <Eigen/Dense>



#include <LoopNetwork.h>
#include <TerminalColors.h>
#include <DislocationNetworkTraits.h>
#include <DislocationNetworkComponent.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
#include <DislocationLoop.h>
#include <DislocationNetworkRemesh.h>
#include <DislocationJunctionFormation.h>
#include <DislocationCrossSlip.h>
//#include <Material.h>
#include <DislocationNetworkIO.h>
#include <DislocationParticle.h>
#include <DislocationStress.h>
//#include <ParticleSystem.h>
#include <MPIcout.h>
//#include <SingleFieldPoint.h>
#include <DDtimeIntegrator.h>
#include <EqualIteratorRange.h>
//#include <BoundingLineSegments.h>
#include <GrainBoundaryTransmission.h>
//#include <GrainBoundaryDissociation.h>
#include <BVPsolver.h>
#include <Polycrystal.h>
#include <DislocationNodeContraction.h>
#include <EshelbyInclusion.h>
#include <TextFileParser.h>
//#include <DisplacementPoint.h>
#include <DefectiveCrystalParameters.h>
#include <ExternalLoadControllerBase.h>
//#include <ExternalLoadController.h>
#include <DislocationInjector.h>
#include <PeriodicDislocationLoop.h>
//#include <PeriodicLoopObserver.h>

//#include <PeriodicDislocationSuperLoop.h>
//#include <PlanarDislocationSuperLoop.h>

#ifdef _MODEL_GREATWHITE_
#include <MooseSolution.h>
#endif

namespace model
{
    
    
    
    template <int _dim, short unsigned int _corder, typename InterpolationType>
    class DislocationNetwork :
    //public PeriodicLoopObserver<PeriodicDislocationLoop<DislocationNetwork<_dim,_corder,InterpolationType>>>
    /*                      */public LoopNetwork<DislocationNetwork<_dim,_corder,InterpolationType> >
    //    /* base                 */ public ParticleSystem<DislocationParticle<_dim> >,
    /*                      */,public std::map<size_t,EshelbyInclusion<_dim>>
    {
        
        
    public:
        
        
        static constexpr int dim=_dim; // make dim available outside class
        static constexpr int corder=_corder; // make dim available outside class
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<double,dim,1>		VectorDim;
        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
        //        typedef DislocationParticle<_dim> DislocationParticleType;
        //        typedef typename DislocationParticleType::StressField StressField;
        //        typedef typename DislocationParticleType::DisplacementField DisplacementField;
        //        typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
        //        typedef typename ParticleSystemType::SpatialCellType SpatialCellType;
        //        typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
        typedef BVPsolver<dim,2> BvpSolverType;
        typedef typename BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename LoopNetworkType::IsNodeType IsNodeType;
        typedef DislocationNetworkIO<DislocationNetworkType> DislocationNetworkIOType;
        typedef Polycrystal<dim> PolycrystalType;
        //        typedef ExternalLoadControllerBase<dim> ExternalLoadControllerType;
        typedef std::map<size_t,EshelbyInclusion<_dim>> EshelbyInclusionContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
#ifdef _MODEL_GREATWHITE_
#include <DislocationNetworkGreatWhite.h>
#endif
        
    private:
        
        
        /**********************************************************************/
        void updateVirtualBoundaryLoops()
        {
            
            
            if(   simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
            {
                //                    if(useVirtualExternalLoops)
                //                    {
                
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"        Updating virtual boundary loops "<<std::flush;
                
                
                // First clean up outdated boundary loops
                std::set<size_t> removeLoops;
                for(const auto& loop : this->loops())
                {
                    if((loop.second->isVirtualBoundaryLoop() && loop.second->links().size()!=4) || loop.second->isPureVirtualBoundaryLoop())
                    {// clean up left over loops from topological operations
                        removeLoops.insert(loop.second->sID);
                    }
                }
                
                for(const size_t& loopID : removeLoops)
                {// Remove the virtual loops with ID in removeLoops
                    this->deleteLoop(loopID);
                }
                
                
                
                // Now reconstruct virtual boundary loops
                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> virtualLoopVector;
                for(const auto& link : this->links())
                {
                    if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers())
                    {
                        virtualLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink,link.second->source,link.second->source->virtualBoundaryNode(),link.second->sink->virtualBoundaryNode()},
                                                       link.second->burgers(),
                                                       (*link.second->grains().begin())->grainID,
                                                       DislocationLoopIO<dim>::VIRTUALLOOP);
                    }
                }
                
                for(const auto& tup : virtualLoopVector)
                {// Insert the new virtual loops
                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
                }
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            }
        }
        
        
        
    public:
        
        const DefectiveCrystalParameters& simulationParameters;
        const SimplicialMesh<dim>& mesh;
        const Polycrystal<dim>& poly;
        GlidePlaneFactory<dim> glidePlaneFactory;
        const std::unique_ptr<PeriodicDislocationLoopFactory<DislocationNetworkType>> periodicDislocationLoopFactory;
        const std::unique_ptr<BVPsolver<dim,2>>& bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>>& externalLoadController;
        const std::vector<VectorDim>& periodicShifts;
        DislocationNetworkRemesh<DislocationNetworkType> networkRemesher;
        DislocationJunctionFormation<DislocationNetworkType> junctionsMaker;
        DislocationNodeContraction<DislocationNetworkType> nodeContractor;
        GrainBoundaryTransmission<DislocationNetworkType> gbTransmission;
        //        MatrixDimD _plasticDistortionFromVelocities;
        MatrixDimD _plasticDistortionFromAreas;
        MatrixDimD _plasticDistortionRateFromVelocities;
        MatrixDimD _plasticDistortionRateFromAreas;
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
        bool outputPlasticDistortion;
        bool outputPlasticDistortionRate;
        bool outputQuadraturePoints;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        bool outputPeriodicConfiguration;
        //        unsigned int _userOutputColumn;
        bool use_stochasticForce;
        double surfaceAttractionDistance;
        //        bool computePlasticDistortionRateFromVelocities;
        std::string folderSuffix;
        
        /**********************************************************************/
        DislocationNetwork(int& argc, char* argv[],
                           const DefectiveCrystalParameters& _simulationParameters,
                           const SimplicialMesh<dim>& _mesh,
                           const Polycrystal<dim>& _poly,
                           const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                           const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                           const std::vector<VectorDim>& _periodicShifts,
                           long int& runID) :
        /* init */ simulationParameters(_simulationParameters)
        /* init */,mesh(_mesh)
        /* init */,poly(_poly)
        /* init */,glidePlaneFactory(poly)
        /* init */,periodicDislocationLoopFactory(simulationParameters.isPeriodicSimulation()? new PeriodicDislocationLoopFactory<DislocationNetworkType>(poly,glidePlaneFactory) : nullptr)
        /* init */,bvpSolver(_bvpSolver)
        /* init */,externalLoadController(_externalLoadController)
        /* init */,periodicShifts(_periodicShifts)
        /* init */,networkRemesher(*this)
        /* init */,junctionsMaker(*this)
        /* init */,nodeContractor(*this)
        /* init */,gbTransmission(*this)
        //        /* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
        ///* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
        //        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true)),
        //        /* init */,totalTime(0.0),
        //        /* init */ dt(0.0),
        //        /* init */ vMax(0.0),
        //        /* init */ Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true)),
        //        /* init */,_plasticDistortionFromVelocities(MatrixDimD::Zero())
        /* init */,_plasticDistortionFromAreas(MatrixDimD::Zero())
        /* init */,_plasticDistortionRateFromVelocities(MatrixDimD::Zero())
        /* init */,_plasticDistortionRateFromAreas(MatrixDimD::Zero())
        /* init */,ddSolverType(TextFileParser("inputFiles/DD.txt").readScalar<int>("ddSolverType",true))
        /* init */,computeDDinteractions(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true))
        /* init */,crossSlipModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("crossSlipModel",true))
        /* init */,outputFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFrequency",true))
        /* init */,outputBinary(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputBinary",true))
        /* init */,outputGlidePlanes(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true))
        /* init */,outputElasticEnergy(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true))
        /* init */,outputMeshDisplacement(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true))
        /* init */,outputFEMsolution(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true))
        /* init */,outputDislocationLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true))
        /* init */,outputPlasticDistortion(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true))
        /* init */,outputPlasticDistortionRate(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true))
        /* init */,outputQuadraturePoints(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true))
        /* init */,outputLinkingNumbers(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true))
        /* init */,outputLoopLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLoopLength",true))
        /* init */,outputSegmentPairDistances(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
        /* init */,outputPeriodicConfiguration(simulationParameters.isPeriodicSimulation()? TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPeriodicConfiguration",true) : false)
        //        /* init */ _userOutputColumn(3)
        /* init */,use_stochasticForce(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true))
        /* init */,surfaceAttractionDistance(TextFileParser("inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true))
        //        /* init */,computePlasticDistortionRateFromVelocities(TextFileParser("inputFiles/DD.txt").readScalar<int>("computePlasticDistortionRateFromVelocities",true))
        /* init */,folderSuffix("")
        {
            
            // Some sanity checks
            //            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
            
            // Initialize static variables
            LinkType::initFromFile("inputFiles/DD.txt");
            NodeType::initFromFile("inputFiles/DD.txt");
            LoopType::initFromFile("inputFiles/DD.txt");
            PeriodicDislocationBase::initFromFile("inputFiles/DD.txt");
            DislocationNetworkComponentType::initFromFile("inputFiles/DD.txt");
            DislocationStressBase<dim>::initFromFile("inputFiles/DD.txt");
            DDtimeIntegrator<0>::initFromFile("inputFiles/DD.txt");
            DislocationCrossSlip<DislocationNetworkType>::initFromFile("inputFiles/DD.txt");
            int stochasticForceSeed=TextFileParser("inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true);
            if(stochasticForceSeed<0)
            {
                StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
            }
            else
            {
                StochasticForceGenerator::init(stochasticForceSeed);
            }
            
            if(argc>1)
            {
                folderSuffix=argv[1];
                //                std::cout<<"folderSuffix="<<folderSuffix<<std::endl;
            }
            //            ParticleSystemType::initMPI(argc,argv);
            
            //
            // Read Vertex and Edge information
            //            DDconfigIO<dim> evl;
            DDconfigIO<dim> evl(folderSuffix);
            evl.read(runID);
            setConfiguration(evl);
            //            std::map<size_t,std::shared_ptr<NodeType>> tempNodes;
            //            createVertices(evl,tempNodes);
            //            createEdges(evl,tempNodes);
            //            updatePlasticDistortionFromAreas(simulationParameters.dt);
            createEshelbyInclusions();
        }
        
        /**********************************************************************/
        void setConfiguration(const DDconfigIO<dim>& evl)
        {
            this->loopLinks().clear(); // erase base network
            std::map<size_t,std::shared_ptr<NodeType>> tempNodes;
            createVertices(evl,tempNodes);
            createEdges(evl,tempNodes);
            updatePlasticDistortionFromAreas(simulationParameters.dt);
#ifdef _MODEL_MPI_
            // Avoid that a processor starts writing before other are done reading
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            // Initializing configuration
//            this->io().output(simulationParameters.runID);
            moveGlide(0.0);    // initial configuration
        }
        
        /**********************************************************************/
        void createVertices(const DDconfigIO<dim>& evl,std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
        {/*!Creates DislocationNode(s) based on the data read by the DDconfigIO<dim>
          * object.
          */
            size_t kk(1);
            for (const auto& node : evl.nodes())
            {
                const size_t nodeIDinFile(node.sID);
                NodeType::set_count(nodeIDinFile);
                if(node.sID==node.masterID)
                {// a regular node is created
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
                    
                    const size_t nodeID(StaticID<NodeType>::nextID());
                    const auto inserted(tempNodes.emplace(std::piecewise_construct,
                                                          std::make_tuple(nodeID),
                                                          std::make_tuple(new NodeType(this,node.P,node.V,node.velocityReduction))));
                    assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
                    assert(inserted.first->first == nodeID && "KEY != nodeID");
                    assert(inserted.first->second->sID == nodeID && "sID != nodeID");
                    assert(nodeID==nodeIDinFile);
                }
                else
                {
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<"), virtual of "<<node.masterID<<std::endl;
                    const auto isNode(this->node(node.masterID));
                    assert(isNode.first);
                    isNode.second->resetVirtualBoundaryNode();
                }
                kk++;
            }
        }
        
        std::shared_ptr<NodeType> getSharedNode(const size_t& nodeID,
                                                const std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
        {
            const auto isNode(this->node(nodeID));
            assert(isNode.first);
            if(isNode.second->masterNode)
            {// a virtual node
                return isNode.second->masterNode->virtualBoundaryNode();
                //                loopNodes.push_back(isNode.second->masterNode->virtualBoundaryNode());
            }
            else
            {
                const auto tempNodesFound(tempNodes.find(nodeID));
                if(tempNodesFound==tempNodes.end())
                {
                    model::cout<<"node "<<nodeID<<" not found"<<std::endl;
                    assert(false && "node shared pointer not found");
                    return nullptr;
                }
                else
                {
                    return tempNodesFound->second;
                }
                //                loopNodes.push_back(isSharedNode.second);
            }
        }
        
        /**********************************************************************/
        void createEdges(const DDconfigIO<dim>& evl,const std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
        {/*!
          */
            
            
//            std::map<size_t,std::shared_ptr<PeriodicDislocationLoopType> > periodicDislocationLoopMap;
//            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
//               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
//            {
//                for(const auto& pLoop : evl.periodicLoops())
//                {// Collect LoopLinks by loop IDs
//                    StaticID<PeriodicDislocationLoopType>::set_count(pLoop.sID);
//                    periodicDislocationLoopMap.emplace(pLoop.sID,new PeriodicDislocationLoopType(this));
//                }
//            }
            
            std::map<size_t,std::map<size_t,size_t>> loopMap;
            for(const auto& looplink : evl.links())
            {// Collect LoopLinks by loop IDs
                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
            }
            
            assert(loopMap.size()==evl.loops().size());
            
            size_t loopLumber=1;
            for(const auto& loop : evl.loops())
            {// for each loop in the DDconfigIO<dim> object
                
                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
                assert(loopFound!=loopMap.end());
                std::vector<std::shared_ptr<NodeType>> loopNodes;
                loopNodes.push_back(getSharedNode(loopFound->second.begin()->first,tempNodes));
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                    if(k<loopFound->second.size()-1)
                    {
                        loopNodes.push_back(getSharedNode(nodeFound->second,tempNodes));
                    }
                    else
                    {
                        assert(nodeFound->second==loopNodes[0]->sID);
                    }
                }
                
                LoopType::set_count(loop.sID);
                
                const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
                
                model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
                switch (loop.loopType)
                {
                    case DislocationLoopIO<dim>::GLISSILELOOP:
                    {
                        LatticePlane loopPlane(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
                        GlidePlaneKey<dim> loopPlaneKey(loop.grainID,loopPlane);
                        if(simulationParameters.isPeriodicSimulation())
                        {
                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),loop.periodicShift)->sID;
                            assert(loop.sID==newLoopID);
                        }
                        else
                        {
                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
                            assert(loop.sID==newLoopID);
                        }
//                        if(loop.periodicLoopID>=0)
//                        {// a loop belonging to a periodic loop
//                            const auto pLoopIter(periodicDislocationLoopMap.find(loop.periodicLoopID));
//                            if(pLoopIter!=periodicDislocationLoopMap.end())
//                            {
//                                const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),pLoopIter->second,loop.periodicShift)->sID;
//                                assert(loop.sID==newLoopID);
//                            }
//                            else
//                            {
//                                assert(false && "PeriodicLoop not found in map");
//                            }
//                        }
//                        else
//                        {
//                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
//                            assert(loop.sID==newLoopID);
//                        }
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::SESSILELOOP:
                    {
                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::VIRTUALLOOP:
                    {
                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    default:
                        assert(false && "Unknown DislocationLoop type");
                        break;
                }
                
                loopLumber++;
            }
            
            model::cout<<std::endl;
        }
        
        /**********************************************************************/
        void createEshelbyInclusions()
        {
            IDreader<'E',1,14,double> inclusionsReader;
            inclusionsReader.read(0,true);
            
            const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
            for(const auto& pair : inclusionsReader)
            {
                
                const size_t& inclusionID(pair.first);
                Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
                
                const VectorDim C(row.template segment<dim>(0));
                const double a(row(dim+0));
                MatrixDimD eT(MatrixDimD::Zero());
                const int typeID(row(13));
                int k=dim+1;
                for(int i=0;i<dim;++i)
                {
                    for(int j=0;j<dim;++j)
                    {
                        eT(i,j)=row(k);
                        k++;
                    }
                }
                
                
                
                EshelbyInclusion<dim>::set_count(inclusionID);
                eshelbyInclusions().emplace(std::piecewise_construct,
                                            std::make_tuple(inclusionID),
                                            std::make_tuple(C,a,eT,poly.nu,poly.mu,inclusionsMobilityReduction[typeID],typeID) );
            }
        }
        
        /**********************************************************************/
        bool remove(const size_t& nodeID)
        {
            const auto isNode=this->node(nodeID);
            if(isNode.first)
            {// remove virtual node together with current node
                if(isNode.second->virtualBoundaryNode())
                {
                    LoopNetworkType::remove(isNode.second->virtualBoundaryNode()->sID);
                }
            }
            return LoopNetworkType::remove(nodeID);
        }
        
        /**********************************************************************/
        void updateGeometry(const double& dt)
        {
            for(auto& loop : this->loops())
            {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
                loop.second->updateGeometry();
            }
            updatePlasticDistortionFromAreas(dt);
        }
        
        //        /**********************************************************************/
        //        double get_dt() const
        //        {
        //            switch (timeIntegrationMethod)
        //            {
        //                case 0:
        //                    return DDtimeIntegrator<0>::get_dt(*this);
        //                    break;
        //
        //                default:
        //                    assert(0 && "time integration method not implemented");
        //                    return 0;
        //                    break;
        //            }
        //        }
        
        /**********************************************************************/
        void removeZeroAreaLoops()
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"        Removing zero-area loops "<<std::flush;
            std::deque<size_t> loopIDs;
            for(const auto& loop : this->loops())
            {
                if(loop.second->slippedArea()<FLT_EPSILON)
                {
                    loopIDs.push_back(loop.second->sID);
                }
            }
            
            for(const auto& loopID : loopIDs)
            {
                this->deleteLoop(loopID);
            }
            model::cout<<"("<<loopIDs.size()<<" removed)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void singleGlideStepDiscreteEvents(const long int& runID)
        {
            //! A simulation step consists of the following:
            //            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
            //            /*                    */<< ", time="<<totalTime
            //            /*                    */<< ": nodes="<<this->nodes().size()
            //            /*                    */<< ", segments="<<this->links().size()
            //            /*                    */<< ", loopSegments="<<this->loopLinks().size()
            //            /*                    */<< ", loops="<<this->loops().size()
            //            /*                    */<< ", components="<<this->components().size()
            //            /*                    */<< defaultColor<<std::endl;
            
            //! 1- Check that all nodes are balanced
            //            checkBalance();
            
            //! 2 - Update quadrature points
            //            updateQuadraturePoints();
            //            updateStressStraightSegments();
            
            //            for(auto& loop : this->loops()) // TODO: PARALLELIZE THIS LOOP
            //            {// copmute slipped areas and right-handed normal
            //                loop.second->updateGeometry();
            //            }
            //            updatePlasticDistortionFromAreas();
            
            //! 3- Calculate BVP correction
            //            updateLoadControllers(runID);
            
            //#ifdef DislocationNucleationFile
            //            if(use_bvp && !(runID%use_bvp))
            //            {
            //                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
            //                updateQuadraturePoints();
            //            }
            //#endif
            //            assembleAndSolve(runID,straightSegmentsDeq);
            //            computeNodaVelocities(runID);
            
            
            //! 4- Solve the equation of motion
            
            
            //! 5- Compute time step dt (based on max nodal velocity) and increment totalTime
            // make_dt();
            
            
            
            //! 6- Output the current configuration before changing it
            //            output(runID);
            //            io().output(runID);
            
            
            //! 7- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
            //            move(dt);
            
            //! 8- Update accumulated quantities (totalTime and plasticDistortion)
            //            totalTime+=dt;
            //            updatePlasticDistortionRateFromVelocities();
            
            
            //! 9- Contract segments of zero-length
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //            if(runID>0)
            //            {
            //                removeZeroAreaLoops();
            //            }
            
            //! 10- Cross Slip (needs upated PK force)
            DislocationCrossSlip<DislocationNetworkType>(*this);
            
            
            //                        gbTransmission.transmit();
            //gbTransmission.directTransmit();
            
            //
            //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
            
            //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
            
            
            //! 11- detect loops that shrink to zero and expand as inverted loops
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);
            
            //! 12- Form Junctions
            junctionsMaker.formJunctions(DDtimeIntegrator<0>::dxMax);
            
            //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
            //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
            
            
            
            //! 13- Node redistribution
            networkRemesher.remesh(runID);
            
            //            createImageLoops();
            //            updateImageLoops();
            updateVirtualBoundaryLoops();
            
            //            mergeLoopsAtNodes();
            
            //            DislocationInjector<DislocationNetworkType>(*this).insertRandomStraightDislocation();
            
            //! 9- Contract segments of zero-length
            //            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();
            
            //! 14- If BVP solver is not used, remove DislocationSegment(s) that exited the boundary
            //            removeBoundarySegments();
            
            //            removeSmallComponents(3.0*dx,4);
            
            //            make_bndNormals();
            
            //! 16 - Increment runID counter
            //            ++runID;     // increment the runID counter
        }
        
        //        const VectorDim periodicVector(const Eigen::Array<int,dim,1>& cellID) const
        //        {
        //            return (meshDimensions.array()*cellID.template cast<double>()).matrix();
        //        }
        
        /**********************************************************************/
        const EshelbyInclusionContainerType& eshelbyInclusions() const
        {
            return *this;
        }
        
        EshelbyInclusionContainerType& eshelbyInclusions()
        {
            return *this;
        }
        
        
        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            return nodeContractor.contract(nA,nB);
        }
        
        
        
        
        
        /**********************************************************************/
        void assembleAndSolveGlide(const long int& runID)
        {/*! Performs the following operatons:
          */
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            //! -1 Compute the interaction StressField between dislocation particles
            
            
            if(corder==0)
            {// For straight segments use analytical expression of stress field
                //                const auto t0= std::chrono::system_clock::now();
                //                model::cout<<"		Collecting StressStraight objects: "<<std::flush;
                //
                //                size_t currentSize=0;
                //                if(computeDDinteractions)
                //                {
                //                    for(const auto& link : this->networkLinks())
                //                    {
                //                        link.second->addToStressStraight(straightSegmentsDeq);
                //                    }
                //
                //                    currentSize=straightSegmentsDeq.size();
                //
                //                    const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());
                //
                //                    for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
                //                    {
                //                        for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
                //                        {
                //                            for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
                //                            {
                //
                //                                const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
                //
                //                                if( cellID.squaredNorm()!=0) //skip current cell
                //                                {
                //                                    for (size_t c=0;c<currentSize;++c)
                //                                    {
                //                                        const VectorDim P0=straightSegmentsDeq[c].P0+(meshSize.array()*cellID.cast<double>().array()).matrix();
                //                                        const VectorDim P1=straightSegmentsDeq[c].P1+(meshSize.array()*cellID.cast<double>().array()).matrix();
                //
                //                                        straightSegmentsDeq.emplace_back(P0,P1,straightSegmentsDeq[c].b);
                //                                    }
                //                                }
                //                            }
                //                        }
                //                    }
                //
                //
                //                }
                //                model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                //
                const auto t1= std::chrono::system_clock::now();
                model::cout<<"		Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
#ifdef _OPENMP
                //                const size_t nThreads = omp_get_max_threads();
                EqualIteratorRange<typename NetworkLinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
                //                SOMETHING WRONG HERE? CPU USE SAYS 100% ONLY?
                
#pragma omp parallel for
                for(size_t thread=0;thread<eir.size();thread++)
                {
                    for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                    {
                        linkIter->second->assembleGlide();
                    }
                }
#else
                for (auto& linkIter : this->links())
                {
                    linkIter.second->assembleGlide();
                }
#endif
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                
                
                
                
                
                
                //#ifdef _OPENMP
                //#pragma omp parallel for
                //#endif
                //                for (unsigned int k=0; k<this->particleSystem().size();++k)
                //                {
                //                    if(this->particleSystem()[k].template fieldPointBase<StressField>().enabled)
                //                    {
                //                        this->particleSystem()[k].template fieldPointBase<StressField>() += StressField::addSourceContribution(this->particleSystem()[k],straightSegmentsDeq);
                //                    }
                //                }
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            }
            else
            {// For curved segments use quandrature integration of stress field
                //                assert(0 && "RE-ENABLE THIS. New Material class is incompatible with old implmentation using static objects");
                
                assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
                
                //
                //                const auto t0= std::chrono::system_clock::now();
                //
                //                if(computeDDinteractions)
                //                {
                //
                //                    if(dislocationImages_x!=0 || dislocationImages_y!=0 || dislocationImages_z!=0)
                //                    {
                //                        assert(0 && "FINISH HERE");
                //                    }
                //
                //                    model::cout<<"		Computing numerical stress field at quadrature points ("<<nThreads<<" threads)..."<<std::flush;
                //                    if (use_extraStraightSegments)
                //                    {
                //                        this->template computeNeighborField<StressField>(ssdeq);
                //                    }
                //                    else
                //                    {
                //                        this->template computeNeighborField<StressField>();
                //                    }
                //                    model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
                //                }
                //
                //                //! -2 Loop over DislocationSegments and assemble stiffness matrix and force vector
                //                const auto t1= std::chrono::system_clock::now();
                //                model::cout<<"		Computing segment stiffness matrices and force vectors ("<<nThreads<<" threads)..."<<std::flush;
                //                typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
                //                assert(0 && "THIS CASE MUST BE REWORKED, since LinkType::assemble DOES NOT EXIST ANYMORE");
                //                //LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
                //                //this->parallelExecute(Lmfp);
                //                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                //
                
            }
            
            //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
            const auto t2= std::chrono::system_clock::now();
            model::cout<<"		Assembling NetworkComponents and solving "<<std::flush;
            
            
            
            switch (ddSolverType)
            {
                case 1: // iterative solver
                {
                    model::cout<<"(MINRES "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // MINRES is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).iterativeSolve();
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).iterativeSolve();
                    }
#endif
                    break;
                }
                    
                case 2: // direct solver
                {
#ifdef _MODEL_PARDISO_SOLVER_
                    model::cout<<"(PardisoLDLT "<<nThreads<<" threads)..."<<std::flush;
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).directSolve();
                    }
#else
                    model::cout<<"(SimplicialLDLT "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // SimplicialLDLT is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).directSolve();
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).directSolve();
                    }
#endif
#endif
                    break;
                }
                    
                default: // lumped solver
                {
                    model::cout<<"(lumpedSolver "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).lumpedSolve(runID);
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).lumpedSolve(runID);
                    }
#endif
                    break;
                }
            }
            
            //            if(DislocationNetworkComponentType::use_directSolver)
            //            {
            //
            //
            //            }
            //            else // iterative solver
            //            {
            //
            //            }
            
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        
        
        
        /**********************************************************************/
        void moveGlide(const double & dt_in)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
            
            model::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            if(simulationParameters.isPeriodicSimulation())
            {
                for (auto& nodeIter : this->nodes())
                {
                    static_cast<typename NodeType::NodeBaseType* const>(nodeIter.second)->set_P(nodeIter.second->get_P()+nodeIter.second->get_V()*dt_in);
                }
                
                for (auto& nodeIter : this->nodes())
                {
                    nodeIter.second->meshFaces().clear();
                }
                
                // remove mesh faces from boundary nodes
                
                for (auto& nodeIter : this->nodes())
                {// trigger regular calls in DislocationNode::moveGlide, but with zero motion
                    nodeIter.second->moveGlide(0.0);
                }
            }
            else
            {
                for (auto& nodeIter : this->nodes())
                {
                    nodeIter.second->moveGlide(dt_in);
                }
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        
        
        /**********************************************************************/
        void updatePlasticDistortionFromAreas(const double& dt)
        {
            //            if(!computePlasticDistortionRateFromVelocities)
            //            {
            const MatrixDimD old(_plasticDistortionFromAreas);
            _plasticDistortionFromAreas.setZero();
            for(const auto& loop : this->loops())
            {
                _plasticDistortionFromAreas+= loop.second->plasticDistortion();
            }
            _plasticDistortionRateFromAreas=(_plasticDistortionFromAreas-old)/dt;
            //            }
        }
        
        //        /**********************************************************************/
        //        const MatrixDimD& plasticDistortionRate() const
        //        {
        //            return computePlasticDistortionRateFromVelocities? _plasticDistortionRateFromVelocities : _plasticDistortionRateFromAreas;
        //        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortionRate() const
        {
            return  _plasticDistortionRateFromAreas;
        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortion() const
        {
            return _plasticDistortionFromAreas;
        }
        
        //        /**********************************************************************/
        //        const MatrixDimD& plasticDistortion() const
        //        {
        //            return computePlasticDistortionRateFromVelocities? _plasticDistortionFromVelocities : _plasticDistortionFromAreas;
        //        }
        
        /**********************************************************************/
        MatrixDimD plasticStrainRate() const
        {/*!\returns the plastic strain rate tensor generated during the last time step.
          */
            //const MatrixDimD temp(plasticDistortionRate());
            return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
        }
        
        /**********************************************************************/
        std::tuple<double,double,double,double> networkLength() const
        {/*!\returns the total line length of the DislocationNetwork. The return
          * value is a tuple, where the first value is the length of bulk glissile
          * dislocations, the second value is the length of bulk sessile
          * dislocations, and the third value is the length accumulated on
          * the mesh boundary.
          */
            double bulkGlissileLength(0.0);
            double bulkSessileLength(0.0);
            double boundaryLength(0.0);
            double grainBoundaryLength(0.0);
            
            for(auto& loop : this->loops())
            {
                for(const auto& loopLink : loop.second->links())
                {
                    if(!loopLink.second->pLink->hasZeroBurgers())
                    {
                        if(loopLink.second->pLink->isBoundarySegment())
                        {
                            boundaryLength+=loopLink.second->pLink->chord().norm();
                        }
                        else if(loopLink.second->pLink->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink.second->pLink->chord().norm();
                        }
                        else
                        {
                            if(loopLink.second->pLink->isSessile())
                            {
                                bulkSessileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
                            }
                        }
                    }
                }
            }
            return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
        }
        
        /**********************************************************************/
        MatrixDimD stress(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            MatrixDimD temp(MatrixDimD::Zero());
            for(const auto& link : this->links())
            {// sum stress field per segment
                if(   !link.second->hasZeroBurgers()
                   && !(link.second->isBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
                   //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
                   )
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp+=link.second->straight.stress(x+shift);
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        void stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<fieldPoints.size();++k)
            {
                fieldPoints[k]=stress(fieldPoints[k].P);
            }
        }
        
        /**********************************************************************/
        VectorDim displacement(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            
            VectorDim temp(VectorDim::Zero());
            
            for(const auto& loop : this->loops())
            {// sum solid angle of each loop
                for(const auto& shift : periodicShifts)
                {
                    temp-=loop.second->solidAngle(x+shift)/4.0/M_PI*loop.second->burgers();
                }
                //                if(!(loop.second->isVirtualBoundaryLoop() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES))
                //                {
                //                    for(const auto& shift : periodicShifts)
                //                    {
                //                        temp-=loop.second->solidAngle(x+shift)/4.0/M_PI*loop.second->burgers();
                //                    }
                //                }
            }
            
            for(const auto& link : this->links())
            {// sum line-integral part of displacement field per segment
                if(   !link.second->hasZeroBurgers()
                   //                   && (!link.second->isBoundarySegment() || simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM)
                   //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
                   )
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp+=link.second->straight.displacement(x+shift);
                    }
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<fieldPoints.size();++k)
            {
                fieldPoints[k]=displacement(fieldPoints[k].P);
            }
        }
        
        
        /**********************************************************************/
        DislocationNetworkIOType io()
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io() const
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
    };
    
}
#endif


//        /**********************************************************************/
//        void updateImageLoops()
//        {
//            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
//               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
//            {
//                const auto t0= std::chrono::system_clock::now();
//                model::cout<<"        Updating periodic  loops "<<std::flush;
//
//                assert(mesh.regions().size()==1);
//                const auto& region(*mesh.regions().begin()->second);
//
//
//                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> imageLoopVector;
//
//
//                for(const auto& link : this->links())
//                {
//                    if(    link.second->isBoundarySegment()
//                       && !link.second->hasZeroBurgers()
//                       && !link.second->hasPeriodicLoop()
//                       //&& !link.second->sink->isVirtualBoundaryNode && !link.second->source->isVirtualBoundaryNode
//                       )
//                    {// First make the check whether the nodes connected to the sinks are not virtualNodes and they lie on the boundary
//
//
//
//                        //model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
//                        //                        const size_t nodeID=this->insertDanglingNode(node.P,node.V,node.velocityReduction).first->first;
//
//
//
//                        //                        const VectorDim imageSourceP(NodeType::imagePosition(link.second->source.get(),link.second->meshFaces()));
//                        //                        const VectorDim imageSinkP(NodeType::imagePosition(link.second->sink.get(),link.second->meshFaces()));
//                        //
//                        //                        std::shared_ptr<NodeType> sourceImage;
//                        //
//                        //                        // search for imageSourceP among images of either source or sink
//                        //
//                        //                        if()
//                        //                        {// position imageSourceP found
//                        //                            // reset sourceImage to that image node
//                        //                        }
//
//                        std::shared_ptr<NodeType> sourceImage(new NodeType(this,link.second->source.get(),link.second));
//                        std::shared_ptr<NodeType>   sinkImage(new NodeType(this,link.second->  sink.get(),link.second));
//
//
//                        link.second->source->imageNodeContainer.insert(sourceImage->sID,sourceImage);
//                        link.second->sink->imageNodeContainer.insert(sinkImage->sID,sinkImage);
//
//                        std::vector<std::shared_ptr<NodeType>> imageNodes;
//                        imageNodes.push_back(sourceImage);
//                        imageNodes.push_back(sinkImage);
//                        imageNodes.push_back(link.second->  sink);
//                        imageNodes.push_back(link.second->source);
//
//                        imageLoopVector.emplace_back(imageNodes,
//                                                     VectorDim::Zero(),
//                                                     region.regionID,
//                                                     DislocationLoopIO<dim>::PERIODICLOOP);
//                    }
//                }
//
//                size_t nLoops(0);
//                for(const auto& tup : imageLoopVector)
//                {// Insert the new virtual loops
//                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
//                    nLoops++;
//                }
//
//                model::cout<<"("<<nLoops<<" periodic loops)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//
//            }
//        }


//        void createImageLoops()
//        {
//            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
//               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
//            {
//
//
//                const auto t0= std::chrono::system_clock::now();
//                model::cout<<"        Creating image  loops "<<std::flush;
//
//
//                for(auto& loop : this->loops())
//                {
//                    if(loop.second->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
//                    {
//                        loop.second->updateBoundaryDecomposition();
//                    }
//                }
//
//
//
//                std::set<const PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>*> loopPairContainer;
//                for(auto& loop : this->loops())
//                {
//                    for(const auto& sequence : loop.second->boundaryLinkSequenceMap)
//                    {
//                        loopPairContainer.insert(sequence.second.loopPair.get());
//                    }
//                }
//
////                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,
////                /*                  */ BoundaryLoopLinkSequence<LoopType>,
////                /*                  */ VectorDim>> imageLoopVector;
//
//                for(const auto& loopPair : loopPairContainer)
//                {
//
//                    PeriodicDislocationSuperLoop<DislocationNetworkType> periodicSuperLoop(*this,*loopPair);
//                    periodicSuperLoop.repair();
//
////                    switch (loopPair->size())
////                    {
////                        case 1:
////                        {// nucleate the new loop
////                            std::cout<<"LoopPair size="<<loopPair->size()<<std::endl;
////
////                            std::vector<std::shared_ptr<NodeType>> imageNodes;
////                            const BoundaryLoopLinkSequence<LoopType>& loopLinkSequence(**loopPair->begin());
////                            const VectorDim bndNormal(region.outNormal(loopLinkSequence.faceIDs));
////                            const VectorDim glideNormal(loopLinkSequence.loop->glidePlane->unitNormal);
////                            const VectorDim glideDir(bndNormal-bndNormal.dot(glideNormal)*glideNormal);
////
////
////                            for(const auto& linkSequence : loopLinkSequence)
////                            {
////
////
////                                const auto startImage(linkSequence.front()->source()->sharedImage(loopLinkSequence.faceIDs));
////                                const auto   endImage(linkSequence.back()   ->sink()->sharedImage(loopLinkSequence.faceIDs));
////
////                                const VectorDim centerNodePos(0.5*(startImage->get_P()+endImage->get_P())+glideDir.normalized()*10.0);
////                                std::shared_ptr<NodeType> centerNode(new NodeType(this,centerNodePos,VectorDim::Zero(),1.0));
////
////
////                                imageNodes.push_back(startImage);
////                                imageNodes.push_back(centerNode);
////                                imageNodes.push_back(endImage);
////
////                            }
////
////                            imageLoopVector.emplace_back(imageNodes,
////                                                         loopLinkSequence,
////                                                         glideNormal);
////
////                            break;
////                        }
////
////                        case 2:
////                        {// update loops
////                            std::cout<<"LoopPair size="<<loopPair->size()<<std::endl;
////
////
////
////                            break;
////                        }
////
////                        default:
////                        {
////                            std::cout<<"LoopPair size="<<loopPair->size()<<std::endl;
////                            assert(false && "LoopPair must gave size 1 or 2");
////                            break;
////                        }
////                    }
//                }
//
//
//
////                for(auto& loop : this->loops())
////                {
////                    if(loop.second->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
////                    {
////                        model::cout<<"        Creating image  of loop "<<loop.second->sID<<std::endl;
////
////                    for(const auto& pair : loop.second->boundaryLinkSequenceMap)
////                    {
////                        const auto mirrorLoop(loop.second->imageLoop(pair.first));
////                        if(mirrorLoop)
////                        {// mirror loops exists
////
////
////
////                        }
////                        else
////                        {// nucleate new loops
////
////                            std::vector<std::shared_ptr<NodeType>> imageNodes;
////                            const VectorDim bndNormal(region.outNormal(pair.first));
////                            const VectorDim glideNormal(loop.second->glidePlane->unitNormal);
////                            const VectorDim glideDir(bndNormal-bndNormal.dot(glideNormal)*glideNormal);
////
////                            for(const auto& linkSequence : pair.second)
////                            {
////                                std::cout<<"face "<<std::flush;
////                                for(const auto& val : pair.first)
////                                {
////                                    std::cout<<val<<" "<<std::endl;
////                                }
////
////                                for(const auto& deq : pair.second)
////                                {
////                                    for(const auto& link : deq)
////                                    {
////                                        std::cout<<link->tag()<<std::endl;
////                                    }
////                                    std::cout<<"---------"<<std::endl;
////
////                                }
////
////                                const auto startImage(linkSequence.front()->source()->sharedImage(pair.first));
////                                const auto   endImage(linkSequence.back()   ->sink()->sharedImage(pair.first));
////
////                                const VectorDim centerNodePos(0.5*(startImage->get_P()+endImage->get_P())+glideDir.normalized()*10.0);
////                                std::shared_ptr<NodeType> centerNode(new NodeType(this,centerNodePos,VectorDim::Zero(),1.0));
////
////
////                                imageNodes.push_back(startImage);
////                                imageNodes.push_back(centerNode);
////                                imageNodes.push_back(endImage);
////
////                            }
////
////                            imageLoopVector.emplace_back(imageNodes,
////                                                         loop.second->flow().cartesian(),
////                                                         glideNormal);
////
////
////
////                        }
////                    }
////                }
////                }
//
//
//
//
////                size_t nLoops(0);
////                for(const auto& tup : imageLoopVector)
////                {// Insert the new virtual loops
////                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<0>(tup)[0]->get_P());
////                    nLoops++;
////                }
//                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
//
//            }
//        }
